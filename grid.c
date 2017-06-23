/*
 * Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
 *
 * This program is free software. It comes without any warranty,
 * to the extent permitted by applicable law. You can redistribute
 * it and/or modify it under the terms of the WTFPL, Version 2, as
 * published by Sam Hocevar. See the accompanying LICENSE file or
 * http://www.wtfpl.net/ for more details.
 */

#include <stdlib.h>
#include "kv.h"
#include "grid.h"
#include "special_function.h"
#include "dynamics.h"
#include <string.h>
#include "NM_measure.h"


//This function returns the normalization constant A for the two-photon initial state used for init_cond=3
void calculate_normalization_const(grid * simulation)
{
   if(simulation->identical_photons) //skip the math 
      simulation->A = 1./sqrt(2.); 
   else
   {
      double k1 = simulation->k1 / simulation->Gamma;
      double k2 = simulation->k2 / simulation->Gamma;
      double alpha1 = simulation->alpha1;
      double alpha2 = simulation->alpha2;

      simulation->A = sqrt( (4.*pow(k1-k2,2) + pow(alpha1+alpha2,2)) / \
	                    (4.*pow(k1-k2,2) + (pow(alpha1,2) + 6.*alpha1*alpha2 + pow(alpha2,2))) );
   }
}


// This function returns the solution psi[j][i] in x<-a subject to two-photon plane wave
double complex plane_wave_BC(int j, int i, grid * simulation)
{
    double x = (i-simulation->origin_index)*simulation->Delta; //check!!!
    double t = j*simulation->Delta;
    double td = simulation->nx*simulation->Delta;
    double w0 = simulation->w0;
    double Gamma = simulation->Gamma;
    double complex K = I*simulation->k;
    double complex W = I*w0 + 0.5*Gamma;
    double complex p = -I*(K - W);

    double complex e_t = I*sqrt(0.5*Gamma)*(cexp(-K*t)-cexp(-W*t))/p;
    double complex sum = 0;

    for(int n=1; n<=(j/simulation->nx); n++)
    {
        double complex temp = ( cexp( n*log(t-n*td) - W*(t-n*td) - lgamma(n+1) ) \
                       - (I*K+w0) * incomplete_gamma_e(n+1, -I*p*(t-n*td), n*clog(I) - (n+1)*clog(p) - K*(t-n*td) ) );

        temp *= cpow(0.5*Gamma, n-0.5);

	// based on my observation, the wavefunction should converge very fast, 
	// so one can just cut the summation off if the precision is reached.
	// this also helps prevent some overflow issue a bit.
        if( cabs(temp) < DBL_EPSILON*cabs(sum) || isnan(cabs(temp)) )
           break;
        else
	   sum += temp;
    }
    e_t -= sum;    
    e_t *= sqrt(2.)*cexp(K*(x-t)); // psi(x,t) = sqrt(2)e^{ik(x-t)}*e(t)
    e_t *= cexp(-0.5*K*td);        //TODO: this phase factor can be eliminated by absorbing into the wavepacket

    if(!isnan(cabs(e_t)))
       return e_t;
    else
    {
       fprintf(stderr, "%s: NaN is produced (at j=%i and i=%i). Abort!\n", __func__, j, i);
       exit(EXIT_FAILURE);
    }
}


// This function returns the solution psi[j][i] in x<-a subject to single-photon exponential wavepacket 
double complex exponential_BC(int j, int i, grid * simulation)
{
   // psi(x,t) = psi(x-t, 0) * e(t)
   return simulation->e1[j] * one_photon_exponential(i-simulation->origin_index-j, simulation->k, simulation->alpha, simulation); 
}


// This function returns the solution psi[j][i] in x<-a subject to two-photon exponential wavepacket 
double complex two_exponential_BC(int j, int i, grid * simulation)
{
   // psi(x,t) = [ \varphi^2(x-t, 0) * e0^1(t) + (1<->2) ]/\sqrt{2} when photon 1 != photon 2
   if(simulation->identical_photons)
      return sqrt(2.) * simulation->A * simulation->e0[j] \
	     * one_photon_exponential(i-simulation->origin_index-j, simulation->k, simulation->alpha, simulation);
   else	 
   {
      static double complex varphi_1, varphi_2; //TODO: make sure both are static
      varphi_1 = one_photon_exponential(i-simulation->origin_index-j, simulation->k1, simulation->alpha1, simulation);
      varphi_2 = one_photon_exponential(i-simulation->origin_index-j, simulation->k2, simulation->alpha2, simulation);
      return (simulation->e0_1[j] * varphi_2 + simulation->e0_2[j] * varphi_1) * simulation->A / sqrt(2.); 
   }
}


//TODO: this should be generalized to acommadate different I.C.
void prepare_qubit_wavefunction(grid * simulation)
{
    //do this only if the exponential wavepacket is used
    if(simulation->init_cond == 2 || simulation->init_cond == 3)
    {
        initialize_e0(simulation);
        initialize_e1(simulation);
    }
}


void initialize_e0(grid * simulation)
{
    if(simulation->identical_photons) //one wavepacket or two identical exponential wavepackets
    {
        simulation->e0 = calloc(simulation->Ny, sizeof(*simulation->e0));
        if(!simulation->e0)
        { 
            fprintf(stderr, "%s: cannot allocate memory. Abort!\n", __func__);
            exit(EXIT_FAILURE);
        }

        int progress = 0;
        for(int j=0; j<simulation->Ny; j++)
        {
            simulation->e0[j] = e0(j, simulation);

            if(j%(simulation->Ny/10)==0)
            {
                printf("%s: %i%% prepared...\r", __func__, progress*10); fflush(stdout);
                progress++;
            }
        }
    }
    else //two different exponential wavepackets
    {
        simulation->e0_1 = calloc(simulation->Ny, sizeof(*simulation->e0_1));
        simulation->e0_2 = calloc(simulation->Ny, sizeof(*simulation->e0_2));
        if(!simulation->e0_1 || !simulation->e0_2)
        { 
            fprintf(stderr, "%s: cannot allocate memory. Abort!\n", __func__);
            exit(EXIT_FAILURE);
        }

        int progress = 0;
        for(int j=0; j<simulation->Ny; j++)
        {
	    simulation->k = simulation->k1; simulation->alpha = simulation->alpha1;
            simulation->e0_1[j] = e0(j, simulation);
	    simulation->k = simulation->k2; simulation->alpha = simulation->alpha2;
            simulation->e0_2[j] = e0(j, simulation);

            if(j%(simulation->Ny/10)==0)
            {
                printf("%s: %i%% prepared...\r", __func__, progress*10); fflush(stdout);
                progress++;
            }
        }
    }
    //TODO: add other I.C. here

    //wash out the status report
    //printf("                                                                           \r"); fflush(stdout);
}


//e1 is the solution of spontaneous emission and is independent of incident wavepackets
void initialize_e1(grid * simulation)
{
    simulation->e1 = calloc(simulation->Ny, sizeof(*simulation->e1));
    if(!simulation->e1)
    { 
        fprintf(stderr, "%s: cannot allocate memory. Abort!\n", __func__);
        exit(EXIT_FAILURE);
    }

    int progress = 0;
    for(int j=0; j<simulation->Ny; j++)
    {
        simulation->e1[j] = e1(j, simulation);

        if(j%(simulation->Ny/10)==0)
        {
            printf("%s: %i%% prepared...\r", __func__, progress*10); fflush(stdout);
            progress++;
        }
    }
}


void initial_condition(grid * simulation)
{// the initial condition is given in-between x/Delta = [-Nx, Nx] for simplicity

    simulation->psit0 = calloc(2*simulation->Nx+1, sizeof(*simulation->psit0));
    if(!simulation->psit0)
    { 
        perror("initial_condition: cannot allocate memory. Abort!\n");
        exit(EXIT_FAILURE);
    }
    simulation->psit0_size = 2*simulation->Nx+1;

    //for nonzero initial conditions
    if(simulation->init_cond == 2) //single-photon exponential wavepacket
    {
       for(int i=0; i<simulation->psit0_size; i++)
           simulation->psit0[i] = one_photon_exponential(i-simulation->Nx, simulation->k, simulation->alpha, simulation);
    }
    //TODO: add other I.C. here
}


void boundary_condition(grid * simulation)
{// the boundary conditions is given for first nx+1 columns with 
 // x/Delta=[-(Nx+nx+1),-(Nx+1)] due to the delay term

    simulation->psix0 = malloc( simulation->Ny*sizeof(*simulation->psix0) );
    if(!simulation->psix0)
    { 
        perror("boundary_condition: cannot allocate memory. Abort!\n");
        exit(EXIT_FAILURE);
    }
    simulation->psix0_x_size = simulation->nx+1;
    simulation->psix0_y_size = 0;

    int progress = 0;
    for(int j=0; j<simulation->Ny; j++)
    {
        simulation->psix0[j] = calloc(simulation->nx+1, sizeof(*simulation->psix0[j])); 
        if(!simulation->psix0[j])
        { 
            fprintf(stderr, "%s: cannot allocate memory at t=%d*Delta. Abort!\n", __func__, j);
            exit(EXIT_FAILURE);
        }
        simulation->psix0_y_size++;
    }

    for(int j=0; j<simulation->psix0_y_size; j++)
    {
        switch(simulation->init_cond)
	{
	   case 1: { //two-photon plane wave
                 for(int i=0; i<simulation->psix0_x_size; i++)
                     simulation->psix0[j][i] = plane_wave_BC(j, i, simulation);
	      }
	      break;
	   case 2: { //single-photon exponential wavepacket
                 for(int i=0; i<simulation->psix0_x_size; i++)
                     simulation->psix0[j][i] = exponential_BC(j, i, simulation);
	      }
	      break;
	   case 3: { //two-photon exponential wavepacket
                 for(int i=0; i<simulation->psix0_x_size; i++)
                     simulation->psix0[j][i] = two_exponential_BC(j, i, simulation);
	      }
	      break;
	   default: { //bad input
              fprintf(stderr, "%s: invalid option. Abort!\n", __func__);
              exit(EXIT_FAILURE);
	      } 
        }

        if(j%(simulation->Ny/10)==0)
        {
            printf("%s: %i%% prepared...\r", __func__, progress*10); fflush(stdout);
            progress++;
        }
    }

    //wash out the status report
    printf("                                                                           \r"); fflush(stdout);
}


void initialize_psi(grid * simulation)
{
    simulation->psi = malloc( simulation->Ny*sizeof(*simulation->psi) );
    if(!simulation->psi)
    { 
        perror("initialize_psi: cannot allocate memory. Abort!\n");
        exit(EXIT_FAILURE);
    }
    simulation->psi_x_size = simulation->Ntotal;
    simulation->psi_y_size = 0;
    for(int j=0; j<simulation->Ny; j++)
    {
        simulation->psi[j] = calloc( simulation->Ntotal, sizeof(*simulation->psi[j]) );
        if(!simulation->psi[j])
        { 
            fprintf(stderr, "%s: cannot allocate memory at t=%d*Delta. Abort!\n", __func__, j);
            exit(EXIT_FAILURE);
        }
        simulation->psi_y_size++;

        // take boundary conditions
        for(int i=0; i<simulation->psix0_x_size; i++)
           simulation->psi[j][i] = simulation->psix0[j][i];
    }
 
    // take the initial condition
    for(int i=0; i<simulation->psit0_size; i++)
       simulation->psi[0][i+simulation->psix0_x_size] = simulation->psit0[i];
}


//a set of checks (poka-yoke) that make sure the input file is sane
void sanity_check(grid * simulation)
{
    //nx must be multiple of 2
    if(simulation->nx % 2) 
    {
        fprintf(stderr, "%s: nx must be an integer multiple of 2. Abort!\n", __func__);
        exit(EXIT_FAILURE);
    }   

    //nx<=2Nx
    if(simulation->nx > 2*simulation->Nx)
    {
        fprintf(stderr, "%s: nx must be smaller than, or at most equal to, twice of Nx (nx<=2Nx). Abort!\n", __func__);
        exit(EXIT_FAILURE);
    } 

    //Nyquist limit
    if(simulation->k >= M_PI/simulation->Delta || simulation->w0 >= M_PI/simulation->Delta)
    {
        fprintf(stderr, "%s: k or w0 must be smaller than pi/Delta in order not to reach the Nyquist limit. Abort!\n", __func__);
        exit(EXIT_FAILURE);
    }

    //it is meaningless if one performs the computation without saving any result
    if(!simulation->save_chi && !simulation->save_psi && !simulation->save_psi_binary && !simulation->measure_NM)
    {
        //fprintf(stderr, "%s: either save_chi or save_psi has to be 1. Abort!\n", __func__);
        fprintf(stderr, "%s: need to specify the output options (available: save_chi, save_psi, save_psi_binary, measure_NM). \
                        Abort!\n", __func__);
        exit(EXIT_FAILURE);
    }

    //if Ny is too small then no result will be written to file
    if(simulation->save_chi && (simulation->Ny <= simulation->Nx + simulation->nx/2))
    {
        fprintf(stderr, "%s: Ny needs to be larger than Nx+nx/2, or \"chi\" will not be stored. Abort!\n", __func__);
        exit(EXIT_FAILURE);
    }

    //check if the initial condition is not correctly given
    //currently the allowed values are:
    //1 (two-photon plane wave)
    //2 (single-photon exponential wavepacket)
    //3 (two-photon exponential wavepackets)
    if(simulation->init_cond < 1 || simulation->init_cond > 3)
    {
        fprintf(stderr, "%s: init_cond has to be 1, 2 or 3. Abort!\n", __func__);
        exit(EXIT_FAILURE);
    }

    if((simulation->init_cond == 1 || simulation->init_cond == 2) && !lookupValue(simulation->parameters_key_value_pair, "k"))
    {
        //k is a mandatory parameter when init_cond is 1 or 2
        fprintf(stderr, "%s: the incident frequency k must be given when init_cond is 1 or 2. Abort!\n", __func__);
        exit(EXIT_FAILURE);
    }

    if(simulation->init_cond == 2) 
    {   
        //want to use exponential wavepacket but forget to set alpha's value
        if(!lookupValue(simulation->parameters_key_value_pair, "alpha"))
	{
           fprintf(stderr, "%s: alpha is not given. Abort!\n", __func__);
           exit(EXIT_FAILURE);
	}

	//always 1 in this case, in case the user prepared it wrong
	simulation->identical_photons = 1;
    }

    if(simulation->init_cond == 3) 
    {
       //make it failed if identical_photons is not explicitly specified for this case
       //so that identical_photons can be safely set to 1 as default
       if(!lookupValue(simulation->parameters_key_value_pair, "identical_photons"))
       {
          fprintf(stderr, "%s: for two-photon wavapacket calculations, identical_photons needs to be specified (either 0 or 1). Abort!\n", 
		 __func__);
          exit(EXIT_FAILURE);
       }

       //two identical photons
       if(simulation->identical_photons && (!lookupValue(simulation->parameters_key_value_pair, "k") ||
	     !lookupValue(simulation->parameters_key_value_pair, "alpha")) )
       {
          fprintf(stderr, "%s: for two identical photons, k and alpha need to be specified. Abort!\n", __func__);
          exit(EXIT_FAILURE);
       }

       //two different photons
       if(!simulation->identical_photons && (!lookupValue(simulation->parameters_key_value_pair, "k1")
             || !lookupValue(simulation->parameters_key_value_pair, "k2")
	     || !lookupValue(simulation->parameters_key_value_pair, "alpha1")
	     || !lookupValue(simulation->parameters_key_value_pair, "alpha2")) )
       {
          fprintf(stderr, "%s: for two different photons, k1, k2, alpha1 and alpha2 need to be specified. Abort!\n", __func__);
          exit(EXIT_FAILURE);
       }

    }

    //calculate_NM_measure only supports well-defined (normalized) wavepackets
    //TODO: allow init_cond=3
    if(simulation->measure_NM && (simulation->init_cond!=2))// || simulation->init_cond!=3))
    {
        fprintf(stderr, "%s: to calculate lambda and mu for NM measures, set init_cond to be 2. Abort!\n", __func__);
        exit(EXIT_FAILURE);
    }
}


void free_initial_boundary_conditions(grid * simulation)
{//free psit0 and psix0 to save memory
    free(simulation->psit0);
    for(int j=0; j<simulation->Ny; j++)
       free(simulation->psix0[j]);
    free(simulation->psix0);
  
    //reset
    simulation->psix0 = NULL;
    simulation->psit0 = NULL;
    simulation->psix0_x_size = 0;
    simulation->psix0_y_size = 0;
    simulation->psit0_size = 0;
}


void free_grid(grid * simulation)
{
    freeKVs(simulation->parameters_key_value_pair);

    //free psit0 and psix0
    //free_initial_boundary_conditions(simulation);

    //free psi
    for(int j=0; j<simulation->Ny; j++)
       free(simulation->psi[j]);
    free(simulation->psi);

    //free e0 & e1 
    //TODO: take care of this part if the code grows!
    if(simulation->init_cond == 2 || simulation->init_cond == 3)
    {
       if(simulation->identical_photons)
          free(simulation->e0);
       else
       {
          free(simulation->e0_1);
          free(simulation->e0_2);
       }
       free(simulation->e1);
    }

    free(simulation);
}


grid * initialize_grid(const char * filename)
{
   grid * FDTDsimulation = malloc(sizeof(*FDTDsimulation));
   if(!FDTDsimulation)
   {
      perror("initialize_grid: could not allocate memory. Abort!\n");
      exit(EXIT_FAILURE);
   }
   
   //read from input
   FDTDsimulation->parameters_key_value_pair = readKVs(filename);

   //initialize from the input parameters
   FDTDsimulation->nx            = atoi(lookupValue(FDTDsimulation->parameters_key_value_pair, "nx"));
   FDTDsimulation->Nx            = atoi(lookupValue(FDTDsimulation->parameters_key_value_pair, "Nx"));
   FDTDsimulation->Ntotal        = 2 * FDTDsimulation->Nx + FDTDsimulation->nx + 2;
   FDTDsimulation->Ny            = atoi(lookupValue(FDTDsimulation->parameters_key_value_pair, "Ny"));
   FDTDsimulation->Delta         = strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "Delta"), NULL);
   //FDTDsimulation->k             = strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "k"), NULL);
   FDTDsimulation->k             = (lookupValue(FDTDsimulation->parameters_key_value_pair, "k") ? \
	                           strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "k"), NULL) : 0); //default: 0
   FDTDsimulation->k1            = (lookupValue(FDTDsimulation->parameters_key_value_pair, "k1") ? \
	                           strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "k1"), NULL) : 0); //default: 0
   FDTDsimulation->k2            = (lookupValue(FDTDsimulation->parameters_key_value_pair, "k2") ? \
	                           strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "k2"), NULL) : 0); //default: 0
   FDTDsimulation->w0            = strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "w0"), NULL);
   FDTDsimulation->Gamma         = strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "gamma"), NULL);
   FDTDsimulation->Lx            = 2 * FDTDsimulation->Nx * FDTDsimulation->Delta;
   FDTDsimulation->Ly            = (FDTDsimulation->Ny-1) * FDTDsimulation->Delta;
   FDTDsimulation->plus_a_index  = FDTDsimulation->Nx + 3*FDTDsimulation->nx/2 + 1;
   FDTDsimulation->minus_a_index = FDTDsimulation->Nx + FDTDsimulation->nx/2 + 1;
   FDTDsimulation->origin_index  = FDTDsimulation->Nx + FDTDsimulation->nx + 1;
   FDTDsimulation->save_chi      = (lookupValue(FDTDsimulation->parameters_key_value_pair, "save_chi") ? \
				   atoi(lookupValue(FDTDsimulation->parameters_key_value_pair, "save_chi")) : 0); //default: off
   FDTDsimulation->save_psi      = (lookupValue(FDTDsimulation->parameters_key_value_pair, "save_psi") ? \
				   atoi(lookupValue(FDTDsimulation->parameters_key_value_pair, "save_psi")) : 0); //default: off
   FDTDsimulation->save_psi_binary = (lookupValue(FDTDsimulation->parameters_key_value_pair, "save_psi_binary") ? \
				   atoi(lookupValue(FDTDsimulation->parameters_key_value_pair, "save_psi_binary")) : 0); //default: off
   FDTDsimulation->init_cond     = (lookupValue(FDTDsimulation->parameters_key_value_pair, "init_cond") ? \
	                           atoi(lookupValue(FDTDsimulation->parameters_key_value_pair, "init_cond")) : 0); //default: 0 (unspecified)
   FDTDsimulation->identical_photons = (lookupValue(FDTDsimulation->parameters_key_value_pair, "identical_photons") ? \
	                           atoi(lookupValue(FDTDsimulation->parameters_key_value_pair, "identical_photons")) : 1); //default: 1 (yes)
   FDTDsimulation->alpha         = (lookupValue(FDTDsimulation->parameters_key_value_pair, "alpha") ? \
	                           strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "alpha"), NULL) : 0); //default: 0
   FDTDsimulation->alpha1        = (lookupValue(FDTDsimulation->parameters_key_value_pair, "alpha1") ? \
	                           strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "alpha1"), NULL) : 0); //default: 0
   FDTDsimulation->alpha2        = (lookupValue(FDTDsimulation->parameters_key_value_pair, "alpha2") ? \
	                           strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "alpha2"), NULL) : 0); //default: 0
   FDTDsimulation->Tstep         = (lookupValue(FDTDsimulation->parameters_key_value_pair, "Tstep") ? \
				   atoi(lookupValue(FDTDsimulation->parameters_key_value_pair, "Tstep")) : 0); //default: 0
   FDTDsimulation->measure_NM    = (lookupValue(FDTDsimulation->parameters_key_value_pair, "measure_NM") ? \
	                           atoi(lookupValue(FDTDsimulation->parameters_key_value_pair, "measure_NM")) : 0); //default: off

   //check the validity of parameters
   sanity_check(FDTDsimulation);

   //calculate the normalization constant
   if(FDTDsimulation->init_cond==3) calculate_normalization_const(FDTDsimulation);

   //initialize arrays
   prepare_qubit_wavefunction(FDTDsimulation);
   initial_condition(FDTDsimulation);
   boundary_condition(FDTDsimulation);
   initialize_psi(FDTDsimulation);

   //save memory
   free_initial_boundary_conditions(FDTDsimulation);

   return FDTDsimulation;
}


void print_initial_condition(grid * simulation)
{
    printf("t=0   ");
    for(int i=0; i<simulation->psit0_size; i++)
        printf("%.2f+%.2fI ", creal(simulation->psit0[i]), cimag(simulation->psit0[i]));
    printf("\n      ");
    for(int i=-simulation->Nx; i<=simulation->Nx; i++)
        printf("x=%.2f ", i*simulation->Delta);
    printf("\n");
}


void print_boundary_condition(grid * simulation)
{
    for(int j=simulation->Ny-1; j>=0; j--)
    {
        printf("          t = %f\n", j*simulation->Delta);
        for(int i=0; i<simulation->psix0_x_size; i++)
        {
            const char * c;
            if(cimag(simulation->psi[j][i])<0)
            {
                c = "";
            }
            else
            {
                c = "+";
            }
            printf("x=%.3f: %.7f%s%.7fI\n", -(simulation->Nx+simulation->nx+1-i)*simulation->Delta, \
                   creal(simulation->psi[j][i]), c, cimag(simulation->psi[j][i]));
        }
        printf("\n");
    }
}


void print_psi(grid * simulation)
{
    for(int j=simulation->Ny-1; j>=0; j--)
    {
        printf("          t = %f\n", j*simulation->Delta);
        for(int i=0; i<simulation->Ntotal; i++)
        {
            const char * c;
            if(cimag(simulation->psi[j][i])<0)
            {
                c = "";
            }
            else
            {
                c = "+";
            }
            //char c[] = (cimag(simulation->psi[j][i])<0?" ":"+");
            printf("x=%.3f: %.7f%s%.7fI\n", -(simulation->Nx+simulation->nx+1-i)*simulation->Delta, \
                   creal(simulation->psi[j][i]), c, cimag(simulation->psi[j][i]));
        }
        printf("\n");
    }
}


//this function stores the computed wavefunction into a file;
//the third argument "part" can be any function converting a 
//double complex to a double, e.g., creal, cimag, cabs, etc. 
void save_psi(grid * simulation, const char * filename, double (*part)(double complex))
{
    char * str = strdup(filename);
    str = realloc(str, (strlen(filename)+10)*sizeof(char) );

    if(part == &creal)
       strcat(str, ".re.out");
    else if(part == &cimag)
       strcat(str, ".im.out");
    else if(part == &cabs)
       strcat(str, ".abs.out");
    else
    {
       fprintf(stderr, "%s: Warning: default filename is used.\n", __func__);
       strcat(str, ".out");
    }

    FILE * f = fopen(str, "w");

    for(int j=0; j<simulation->Ny; j+=(simulation->Tstep+1))
    {
        for(int i=0; i<simulation->Ntotal; i++)
        {
            fprintf( f, "%.5g ", part(simulation->psi[j][i]) );
        }
//        for(int i=0; i<simulation->Ntotal; i++)
//        {
//            fprintf( f, "%.5f ", part(simulation->psi[j][i]) );
//        }
        fprintf( f, "\n");
    }

    fclose(f);
    free(str);
}


//this function stores the computed wavefunction into a binary file
//note that each data point is a complex number which takes 16 bytes!
void save_psi_binary(grid * simulation, const char * filename)
{
    char * str = strdup(filename);
    str = realloc(str, (strlen(filename)+10)*sizeof(char) );
    strcat(str, ".bin");

    FILE * f = fopen(str, "wb");
    if(!f)
    {
        fprintf(stderr, "%s: file cannot be created!", __func__);
        exit(EXIT_FAILURE);
    }

    size_t array_size = simulation->Ntotal - simulation->minus_a_index;

    for(int j=0; j<simulation->Ny; j+=(simulation->Tstep+1))
    {
        fwrite(simulation->psi[j] + simulation->minus_a_index, sizeof(double complex), array_size, f);
    }

    fclose(f);
    free(str);
}


//this function computes the two-photon wavefunction on the fly and
//then writes to a file, so no extra memory is allocated;
//the third argument "part" can be any function converting a
//complex to a double, e.g., creal, cimag, cabs, etc.
void save_chi(grid * simulation, const char * filename, double (*part)(double complex))
{
    char * str = strdup(filename);
    str = realloc(str, (strlen(filename)+15)*sizeof(char) );

    if(part == &creal)
       strcat(str, ".re_chi.out");
    else if(part == &cimag)
       strcat(str, ".im_chi.out");
    else if(part == &cabs)
       strcat(str, ".abs_chi.out");
    else
    {
       fprintf(stderr, "%s: Warning: default filename is used.\n", __func__);
       strcat(str, ".chi.out");
    }

    FILE * f = fopen(str, "w");

    //compute chi(a+Delta, a+Delta+tau, t) with tau=i*Delta and t=j*Delta:
    //to make all terms in chi well-defined requires 0 <= i <= Nx-nx/2.
    //
    //Update: To access transient dynamics for two photons, j now starts from 0 instead of minus_a_index (=Nx+nx/2+1)
    //
    //(In the previous version, j >= simulation->minus_a_index in order to let signal from the 1st qubit reach the boundary;
    //put it differently, one cannot take data before the first light cone intersects with the boundary x=Nx*Delta.)
    for(int j=0; j<=simulation->Ny; j+=(simulation->Tstep+1))
    {
        for(int i=0; i<=simulation->Nx-simulation->nx/2; i++)
        {
            double complex chi = 0;
	    double complex temp = 0;

            if(simulation->init_cond == 1 || simulation->init_cond == 3)
	       chi += two_photon_input(simulation->nx/2+1-j, simulation->nx/2+1+i-j, simulation);

            if( j>=(simulation->nx+i+1) ) 
	       temp += simulation->psi[j-(simulation->nx+i+1)][simulation->minus_a_index-i];

            if( j>=(i+1) ) 
	       temp -= simulation->psi[j-(i+1)][simulation->plus_a_index-i];

	    if( j>=(simulation->nx+1) )
	       temp += simulation->psi[j-(simulation->nx+1)][simulation->minus_a_index+i];

	    if( j>=1 )
	       temp -= simulation->psi[j-1][simulation->plus_a_index+i];

            chi -= sqrt(simulation->Gamma)/2.0 * temp;
            fprintf( f, "%.5g ", part(chi) );
        }
        fprintf( f, "\n");
    }

    fclose(f);
    free(str);
}


void print_grid(grid * simulation)
{
    printf("nx = %d\n", simulation->nx); 
    printf("Nx = %d\n", simulation->Nx); 
    printf("Ntotal = %d\n", simulation->Ntotal); 
    printf("Ny = %d\n", simulation->Ny); 
    printf("Delta = %.3f\n", simulation->Delta); 
    printf("k = %.3f\n", simulation->k); 
    printf("w0 = %.3f\n", simulation->w0); 
    printf("Gamma = %.3f\n", simulation->Gamma); 
    printf("Lx = %.3f\n", simulation->Lx); 
    printf("Ly = %.3f\n", simulation->Ly);
 
    for(int j=simulation->Ny-1; j>=0; j--)
    {
        for(int i=0; i<simulation->Ntotal; i++)
        {
            printf("%.2f+%.2fI ", creal(simulation->psi[j][i]), cimag(simulation->psi[j][i]));
        }
        printf("\n");
    }
}


void save_psi_square_integral(grid * simulation, const char * filename)
{
    char * str = strdup(filename);
    str = realloc(str, (strlen(filename)+18)*sizeof(char) );
    strcat(str, ".psi_square.out");

    int Tmax = (simulation->Ny-1 < simulation->Nx - simulation->nx/2 ? simulation->Ny-1 : simulation->Nx - simulation->nx/2);
    FILE * f = fopen(str, "w");

    for(int j=0; j<Tmax; j++)
       fprintf( f, "%.10g\n", psi_square_integral(j, simulation) );

    fclose(f);
    free(str);
}

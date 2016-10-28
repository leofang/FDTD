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
#include <math.h>
#include "kv.h"
#include "grid.h"
#include "dynamics.h"
#include <string.h>
#include <float.h> //for DBL_EPSILON ~ 2.2E-16


// Compute the incomplete Gamma function gamma(n,x). Note that gamma(n,x) 
// is different from Gamma(n,x) (gamma(n,x)+Gamma(n,x)=1), and can be 
// called in Mathematica by GammaRegularized[n, 0, x].
//
// The following implementation is based on Numerical Recipes Ch.6.2, where 
// the infinite series representation Eq.6.2.5 is used for better precision 
// and stability (the finite sum, Eq.10.70, in Arfken, 5th ed, P.661, which 
// was implemented in the earlier version, is not useful as it subtracts two 
// nearly same numbers). See also ASA032, ASA147 and ASA239.
//
// Note the normalization factor (n-1)!.
complex incomplete_gamma(int n, complex x)
{
   if(n<=0)
   {
      fprintf(stderr, "Error in %s: n must >= 1. Abort!\n", __func__);
      exit(EXIT_FAILURE);
   }

   //gamma(n>=1, x=0)=0, no need to do real computation
   if(x==0)
      return 0;

   //compute exp(-x)*x^n/(n-1)!
   complex prefactor = cexp(n*clog(x)-x-lgamma(n));

   //compute the infinite sum
   complex sum = 1.0/n;
   complex temp = 1.0/n;
   for(int i=1; ; i++)
   {
      temp *= x/(double)(n+i);
      sum += temp;

      if(cabs(temp) < DBL_EPSILON*cabs(sum)) //stop the sum when temp is extremely small
         return prefactor*sum;
   }
}


// This function returns the solution psi[j][i] in x<-a subject to incident plane wave
complex plane_wave_BC(int j, int i, grid * simulation)
{
    double x = (i-simulation->origin_index)*simulation->Delta; //check!!!
    double t = j*simulation->Delta;
    double td = simulation->nx*simulation->Delta;
    double k = simulation->k;
    double w0 = simulation->w0;
    double Gamma = simulation->Gamma;
    complex p = k - w0 + 0.5*I*Gamma;

    complex e_t = I*sqrt(0.5*Gamma)*cexp(-0.5*I*k*td)*(cexp(-I*k*t)-cexp(-I*w0*t-0.5*Gamma*t))/p;
    complex sum = 0;

    for(int n=1; n<=(j/simulation->nx); n++)
    {
        complex temp = cpow(0.5*Gamma, n-0.5) * \
                       ( cpow(t-n*td, n) * cexp(n*(I*w0*td+0.5*Gamma*td)-I*w0*t-0.5*Gamma*t-lgamma(n+1) ) \
                       + cpow(I, n)*(k-w0)*incomplete_gamma(n+1, -I*p*(t-n*td))*cexp(I*n*k*td-I*k*t)/cpow(p, n+1));
	sum += temp;
	// based on my observation, for some chosen parameters the wavefunction converges 
	// very fast, so one can just cut the summation off if the precision is reached.
	// this also helps prevent some overflow issue a bit.
        if(cabs(temp) < DBL_EPSILON*cabs(sum))
           break;
    }
    e_t -= cexp(-0.5*I*k*td)*sum;
    e_t *= sqrt(2.)*cexp(I*k*(x-t)); // psi(x,t) = sqrt(2)e^{ik(x-t)}*e(t)

    return e_t;
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

    //TODO: make a flexible option for other initial conditions
//    for(int i=0; i<simulation->psit0_size; i++)
//    {
//        simulation->psit0[i] = 0; //for two-photon plane wave input
//    }
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

    //TODO: make a flexible option for other boundary conditions
    for(int j=0; j<simulation->psix0_y_size; j++)
    {
        for(int i=0; i<simulation->psix0_x_size; i++)
            simulation->psix0[j][i] = plane_wave_BC(j, i, simulation);
    }
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
    free_initial_boundary_conditions(simulation);

    for(int j=0; j<simulation->Ny; j++)
       free(simulation->psi[j]);
    free(simulation->psi);

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
   FDTDsimulation->k             = strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "k"), NULL);
   FDTDsimulation->w0            = strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "w0"), NULL);
   FDTDsimulation->Gamma         = strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "Gamma"), NULL);
   FDTDsimulation->Lx            = 2 * FDTDsimulation->Nx * FDTDsimulation->Delta;
   FDTDsimulation->Ly            = (FDTDsimulation->Ny-1) * FDTDsimulation->Delta;
   FDTDsimulation->plus_a_index  = FDTDsimulation->Nx + 3*FDTDsimulation->nx/2 + 1;
   FDTDsimulation->minus_a_index = FDTDsimulation->Nx + FDTDsimulation->nx/2 + 1;
   FDTDsimulation->origin_index  = FDTDsimulation->Nx + FDTDsimulation->nx + 1;

   //check the validity of parameters
   sanity_check(FDTDsimulation);

   //initialize arrays
   initial_condition(FDTDsimulation);
   boundary_condition(FDTDsimulation);
   initialize_psi(FDTDsimulation);

//   //save memory
//   free_initial_boundary_conditions(FDTDsimulation);

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
//complex to a double, e.g., creal, cimag, cabs, etc. 
void save_psi(grid * simulation, const char * filename, double (*part)(complex))
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

    for(int j=0; j<simulation->Ny; j++)
    {
        for(int i=0; i<simulation->Ntotal; i++)
        {
            fprintf( f, "%.7f ", part(simulation->psi[j][i]) );
        }
        fprintf( f, "\n");
    }

    fclose(f);
    free(str);
}


//this function computes the two-photon wavefunction on the fly and
//then writes to a file, so no extra memory is allocated;
//the third argument "part" can be any function converting a
//complex to a double, e.g., creal, cimag, cabs, etc.
void save_chi(grid * simulation, const char * filename, double (*part)(complex))
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

    //compute chi(a+Delta, a+Delta+tau, t) with tau=i*Delta and t=j*Delta
    //j must >= simulation->minus_a_index in order to let signal from the 1st qubit reach the boundary
    complex chi = 0;
    for(int j=(simulation->Nx+simulation->nx/2+1); j<=simulation->Ny; j++)
    {
        for(int i=0; i<=simulation->Nx-simulation->nx/2; i++)
        {
	    chi = two_photon_input(simulation->plus_a_index+1, simulation->plus_a_index+1+i, simulation)-sqrt(simulation->Gamma)/2.0* \
		  (  simulation->psi[j-(simulation->nx+i+1)][simulation->minus_a_index-i] \
		   - simulation->psi[j-(i+1)][simulation->plus_a_index-i] \
		   + simulation->psi[j-(simulation->nx+1)][simulation->minus_a_index+i] \
		   - simulation->psi[j-1][simulation->plus_a_index+i] \
		  );
            fprintf( f, "%.4f ", part(chi) );
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

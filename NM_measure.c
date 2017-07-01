/*
 * Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
 *
 * This program is free software. It comes without any warranty,
 * to the extent permitted by applicable law. You can redistribute
 * it and/or modify it under the terms of the WTFPL, Version 2, as
 * published by Sam Hocevar. See the accompanying LICENSE file or
 * http://www.wtfpl.net/ for more details.
 */

#include "NM_measure.h"
#include <string.h>
#include "special_function.h"
#include "dynamics.h"

// This function returns the qubit wavefunction e0(t) in the single-excitation sector
// subject to e0(0)=0 and an exponential wavepacket
double complex e0(int j, grid * simulation)
{
    if(j<=0) return 0;

    double t         = j*simulation->Delta;
    double td        = simulation->nx*simulation->Delta;
    double w0        = simulation->w0;
    double Gamma     = simulation->Gamma;
    double alpha     = simulation->alpha;
    double complex K = I*simulation->k + 0.5*alpha*Gamma;
    double complex W = I*w0 + 0.5*Gamma;
    double complex p = -I*(K - W);

    double complex e_t = sqrt(0.5*alpha)*Gamma*(cexp(-W*t)-cexp(-K*t))/p;
    double complex sum = 0;

    for(int n=1; n<=(j/simulation->nx); n++)
    {
        double complex temp = ( cexp( n*log(t-n*td) - W*(t-n*td) - lgamma(n+1) ) \
                       - (I*K+w0) * incomplete_gamma_e(n+1, -I*p*(t-n*td), n*clog(I) - (n+1)*clog(p) - K*(t-n*td) ) );

        temp *= pow(0.5*Gamma, n-0.5);

	// based on my observation, the wavefunction should converge very fast, 
	// so one can just cut the summation off if the precision is reached.
	// this also helps prevent some overflow issue a bit.
        if( cabs(temp) < DBL_EPSILON*cabs(sum) || isnan(cabs(temp)) )
           break;
        else
	   sum += temp;
    }
    e_t -= I*sqrt(alpha*Gamma)*sum;
    e_t *= cexp(-0.5*I*(simulation->k)*td); //TODO: this phase factor can be eliminated by absorbing into the wavepacket

    if(!isnan(cabs(e_t)))
       return e_t;
    else
    {
       fprintf(stderr, "%s: NaN is produced (at j=%i). Abort!\n", __func__, j);
       exit(EXIT_FAILURE);
    }
}


// This function returns the qubit wavefunction e1(t) in the single-excitation sector
// subject to e1(0)=1 and no incidet wavepacket (spontaneous emission)
double complex e1(int j, grid * simulation)
{
    if(j<0) return 0;
    if(j==0) return 1.;

    double t         = j*simulation->Delta;
    double td        = simulation->nx*simulation->Delta;
    double w0        = simulation->w0;
    double Gamma     = simulation->Gamma;
    double complex W = I*w0 + 0.5*Gamma;

    double complex sum = 0;
    for(int n=1; n<=(j/simulation->nx); n++)
    {
        double complex temp = exp(-lgamma(n+1)) * cpow(0.5*Gamma*cexp(W*td)*(t-n*td), n);

	// based on my observation, the wavefunction should converge very fast, 
	// so one can just cut the summation off if the precision is reached.
	// this also helps prevent some overflow issue a bit.
        if( cabs(temp) < DBL_EPSILON*cabs(sum) || isnan(cabs(temp)) )
           break;
        else
	   sum += temp;
    }
    double complex e_t = cexp(-W*t)*(1.0+sum);

    if(!isnan(cabs(e_t)))
       return e_t;
    else
    {
       fprintf(stderr, "%s: NaN is produced (at j=%i). Abort!\n", __func__, j);
       exit(EXIT_FAILURE);
    }
}


//compute the photon wavefunction phi(x,t) in the single-excitation sector
double complex phi(int j, int i, grid * simulation)
{
   //cannot go beyond "the box"
   if(j<0 || j>=simulation->Ny || i<0 || i>simulation->Ntotal)
   {
      fprintf(stderr, "%s: argument is outside the simulation region (j=%i, i=%i). Abort!\n", __func__, j, i);
      exit(EXIT_FAILURE);
   }

   double complex Phi = one_photon_exponential(i-simulation->origin_index-j, simulation->k, simulation->alpha, simulation);

   //convert the index to real position i=x/Delta
   i = i - simulation->minus_a_index - simulation->nx/2;

   if(j-i-simulation->nx/2 >= 0 && i>-simulation->nx/2) 
       Phi -= sqrt(0.5*simulation->Gamma) * simulation->e0[j-i-simulation->nx/2];

   if(j-i+simulation->nx/2 >= 0 && i>simulation->nx/2) 
       Phi += sqrt(0.5*simulation->Gamma) * simulation->e0[j-i+simulation->nx/2];

   return Phi;
}


//lambda(t) = \int dx |psi(x,t)|^2 - |e0(t)|^2
double lambda(int j, grid * simulation)
{
   if(j<0)
   {
      fprintf(stderr, "%s: argument j is negative (j=%i). Abort!\n", __func__, j);
      exit(EXIT_FAILURE);
   }

   if(j==0) return 1.0;

   return psi_square_integral(j, simulation) - pow(cabs(simulation->e0[j]), 2.0);
}


//mu(t) = \int dx \phi^*(x,t)\psi(x,t)
//      = exp(-alpha Gamma t)e1(t) + \int_{-a}^\infty dx \phi^*(x,t)\psi(x,t)
double complex mu(int j, grid * simulation)
{
   if(j<0)
   {
      fprintf(stderr, "%s: argument j is negative (j=%i). Abort!\n", __func__, j);
      exit(EXIT_FAILURE);
   }

   if(j==0) return 1.0;

   double complex Mu = exp(- simulation->alpha * simulation->Gamma * j * simulation->Delta) * simulation->e1[j];
   double complex sum = 0;

   int xmax = j + simulation->plus_a_index;
   //trapezoidal rule
   sum += 0.5 * conj( phi(j, simulation->minus_a_index, simulation) ) * simulation->psi[j][simulation->minus_a_index];
   for(int i = simulation->minus_a_index+1; i<xmax; i++)
      sum += conj( phi(j, i, simulation) ) * simulation->psi[j][i];
   sum += 0.5 * conj( phi(j, xmax, simulation) ) * simulation->psi[j][xmax];
 
   Mu += simulation->Delta * sum;

   return Mu;
}


void save_e0(grid * simulation, const char * filename, double (*part)(double complex))
{
    char * str = strdup(filename);
    str = realloc(str, (strlen(filename)+14)*sizeof(char) );

    if(part == &creal)
       strcat(str, ".re_e0.out");
    else if(part == &cimag)
       strcat(str, ".im_e0.out");
    else if(part == &cabs)
       strcat(str, ".abs_e0.out");
    else
    {
       fprintf(stderr, "%s: Warning: default filename is used.\n", __func__);
       strcat(str, ".out");
    }

    FILE * f = fopen(str, "w");

    for(int j=0; j<simulation->Ny; j++)
       fprintf( f, "%.10g\n", part(simulation->e0[j]) );

    fclose(f);
    free(str);
}


void save_e1(grid * simulation, const char * filename, double (*part)(double complex))
{
    char * str = strdup(filename);
    str = realloc(str, (strlen(filename)+14)*sizeof(char) );

    if(part == &creal)
       strcat(str, ".re_e1.out");
    else if(part == &cimag)
       strcat(str, ".im_e1.out");
    else if(part == &cabs)
       strcat(str, ".abs_e1.out");
    else
    {
       fprintf(stderr, "%s: Warning: default filename is used.\n", __func__);
       strcat(str, ".out");
    }

    FILE * f = fopen(str, "w");

    for(int j=0; j<simulation->Ny; j++)
       fprintf( f, "%.10g\n", part(simulation->e1[j]) );

    fclose(f);
    free(str);
}


void save_mu(grid * simulation, const char * filename, double (*part)(double complex))
{
    char * str = strdup(filename);
    str = realloc(str, (strlen(filename)+14)*sizeof(char) );

    if(part == &creal)
       strcat(str, ".re_mu.out");
    else if(part == &cimag)
       strcat(str, ".im_mu.out");
    else if(part == &cabs)
       strcat(str, ".abs_mu.out");
    else
    {
       fprintf(stderr, "%s: Warning: default filename is used.\n", __func__);
       strcat(str, ".out");
    }

    int Tmax = (simulation->Ny-1 < simulation->Nx - simulation->nx/2 ? simulation->Ny-1 : simulation->Nx - simulation->nx/2);
    FILE * f = fopen(str, "w");

    for(int j=0; j<Tmax; j++)
       fprintf( f, "%.10g\n", part(simulation->mu[j]) );

    fclose(f);
    free(str);
}


void save_lambda(grid * simulation, const char * filename)
{
    char * str = strdup(filename);
    str = realloc(str, (strlen(filename)+14)*sizeof(char) );
    strcat(str, ".lambda.out");

    int Tmax = (simulation->Ny-1 < simulation->Nx - simulation->nx/2 ? simulation->Ny-1 : simulation->Nx - simulation->nx/2);
    FILE * f = fopen(str, "w");

    for(int j=0; j<Tmax; j++)
       fprintf( f, "%.10g\n", simulation->lambda[j] );

    fclose(f);
    free(str);
}


//calculate and save lambda(t) and mu(t)
void calculate_NM_measure(grid * simulation, const char * filename)
{
    int Tmax = (simulation->Ny-1 < simulation->Nx - simulation->nx/2 ? simulation->Ny-1 : simulation->Nx - simulation->nx/2);

    //their memory is not allocated until this function is called
    simulation->lambda = calloc(Tmax, sizeof(*simulation->lambda));
    simulation->mu     = calloc(Tmax, sizeof(*simulation->mu));
    if(!simulation->lambda || !simulation->mu)
    {
       fprintf(stderr, "%s: memory allocation failed. Abort!\n", __func__);
       exit(EXIT_FAILURE);
    }

    for(int j=0; j<Tmax; j++)
    {
       simulation->lambda[j] = lambda(j, simulation);
       simulation->mu[j]     = mu(j, simulation);
    }

    save_lambda(simulation, filename);
    save_mu(simulation, filename, creal);
    save_mu(simulation, filename, cimag);

    free(simulation->lambda);
    free(simulation->mu);
}

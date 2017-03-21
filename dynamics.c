/*
 * Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
 *
 * This program is free software. It comes without any warranty,
 * to the extent permitted by applicable law. You can redistribute
 * it and/or modify it under the terms of the WTFPL, Version 2, as
 * published by Sam Hocevar. See the accompanying LICENSE file or
 * http://www.wtfpl.net/ for more details.
 */

#include "dynamics.h"

// this function computes the average of 4 points that form a square
// j and i are the array indices psi[j][i] of the upper right corner
double complex square_average(int j, int i, grid * simulation)
{
   if(i<1) //beyond the boundary x=-(Nx+nx+1)*Delta
   {
      fprintf(stderr, "Error in %s: beyond left boundary. Abort!\n", __func__);
      exit(EXIT_FAILURE);
   }

   if(i>simulation->Ntotal-1) //beyond the boundary x=Nx*Delta
   {
      fprintf(stderr, "Error in %s: beyond right boundary. Abort!\n", __func__);
      exit(EXIT_FAILURE);
   }

   if(j<1) //everything is zero below t=0
      return 0;

   double complex square = 0;

   for(int m=j-1; m<=j; m++)
   {
      for(int n=i-1; n<=i; n++)
         square += simulation->psi[m][n];
   }
   square/=4.;

   return square;
}


// this function computes the average of 2 points that form a bar
// j and i are the array indices psi[j][i] of the right end
double complex bar_average(int j, int i, grid * simulation)
{
   if(i<1) //beyond the boundary x=-(Nx+nx+1)*Delta
   {
      fprintf(stderr, "Error in %s: beyond left boundary. Abort!\n", __func__);
      exit(EXIT_FAILURE);
   }

   if(i>simulation->Ntotal-1) //beyond the boundary x=Nx*Delta
   {
      fprintf(stderr, "Error in %s: beyond right boundary. Abort!\n", __func__);
      exit(EXIT_FAILURE);
   }

   if(j<0) //everything is zero below t=0
      return 0;

   return (simulation->psi[j][i]+simulation->psi[j][i-1])/2.;
}


//this function computes chi(x1, x2, 0)
double complex two_photon_input(int i1, int i2, grid * simulation)
{
   //TODO: make a flexible choice for different inputs

   double x1 = i1 * simulation->Delta;
   double x2 = (i2+0.5) * simulation->Delta; //shift +0.5 due to Taylor expansion at the center of square //TODO: find a better way to write it

   //Two-photon plane waves
   return cexp( I * simulation->k * (x1+x2) );
}


// this function computes the exponential wavepacket with a sharp wavefront at x=-a
double complex one_photon_exponential(int i, grid * simulation)
{
   if(i>simulation->minus_a_index)
      return 0;

   double x = (i - simulation->origin_index) * simulation->Delta;
   double a_g = simulation->alpha*simulation->Gamma; 
   
   return I*csqrt(a_g) * cexp(I*simulation->k*x + 0.5*a_g*(x+0.5*simulation->nx*simulation->Delta));
}

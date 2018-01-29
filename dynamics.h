/*
 * Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
 *
 * This program is free software. It comes without any warranty,
 * to the extent permitted by applicable law. You can redistribute
 * it and/or modify it under the terms of the WTFPL, Version 2, as
 * published by Sam Hocevar. See the accompanying LICENSE file or
 * http://www.wtfpl.net/ for more details.
 */

#ifndef __DYNAMICS_H__
#define __DYNAMICS_H__

#include <math.h>
#include "grid.h"

/* Update Jun 21, 2017: 
 * To utilize the C99 inline specifier, it is required to put the function definitions 
 * in the header file which is visible to all source files, and then to provide an 
 * "extern inline" declaration in one and only one source file (translation unit to be 
 * precise). The "extern inline" specifier there can be either replaced by "extern" or 
 * removed completely, but it may NOT be replaced by "inline" alone. For further 
 * discussions, see the references below:
 *   http://www.greenend.org.uk/rjk/tech/inline.html
 *   https://gcc.gnu.org/onlinedocs/gcc/Inline.html
 *   https://stackoverflow.com/questions/6312597/
 *   https://gustedt.wordpress.com/2010/11/29/myth-and-reality-about-inline-in-c99/
 *   http://www.drdobbs.com/the-new-c-inline-functions/184401540
 */

// this function computes the average of 4 points that form a square
// j and i are the array indices psi[j][i] of the upper right corner
inline double complex square_average(int j, int i, grid * simulation)
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
inline double complex bar_average(int j, int i, grid * simulation)
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


// this function computes the exponential wavepacket with a sharp wavefront at x=-a
// update: argument x now refer to the "unit-less" coordinates, so "true x" = x * Delta
inline double complex one_photon_exponential(double x, double k, double alpha, grid * simulation)
{
   if(x>-simulation->nx/2)
      return 0;

   double a_g = alpha * simulation->Gamma; 
   
   return I*csqrt(a_g) * cexp( (I*k*x + 0.5*a_g*(x+0.5*simulation->nx)) * simulation->Delta );
}


//this function computes chi(x1, x2, 0)
//update: arguments x1 & x2 now refer to the "unit-less" coordinates, so "true x1" = x1 * Delta and so on
inline double complex two_photon_input(double x1, double x2, grid * simulation)
{
   double complex chi = 0;
   switch(simulation->init_cond)
   {
      //two-photon plane waves
      case 1: { chi = cexp( I * simulation->k * (x1+x2) * simulation->Delta); } break;

      //two-photon exponentail wavepacket (init_cond=3)
      case 3: {
         if(simulation->identical_photons)
            chi = one_photon_exponential(x1, simulation->k, simulation->alpha, simulation) \
                  * one_photon_exponential(x2, simulation->k, simulation->alpha, simulation);
         else
         {
            chi = simulation->A / sqrt(2.) *
                  ( one_photon_exponential(x1, simulation->k1, simulation->alpha1, simulation) \
                    * one_photon_exponential(x2, simulation->k2, simulation->alpha2, simulation) \
                    + one_photon_exponential(x2, simulation->k1, simulation->alpha1, simulation) \
                    * one_photon_exponential(x1, simulation->k2, simulation->alpha2, simulation) );
         }
      } break;

      //TODO: add other different inputs here
      
      default: { exit(EXIT_FAILURE); }
   }

   return chi;
}


//this function calculates \int dx |\psi(x,t)|^2
inline double psi_square_integral(int j, grid * simulation)
{
   if(j<0)
   {
      fprintf(stderr, "%s: argument j is negative (j=%i). Abort!\n", __func__, j);
      exit(EXIT_FAILURE);
   }

   double sum = 0;
   double Lambda = 0;
   double t = j*simulation->Delta;
   switch(simulation->init_cond)
   {
      case 2: {// = e^(-alpha*Gamma*t)|e1(t)|^2 + \int_{-a}^\infty dx |psi(x,t)|^2
            if(j==0) return 1.0;

            Lambda += exp(-simulation->alpha*simulation->Gamma*t);
            Lambda *= pow(cabs(simulation->e1[j]), 2.0);
	 }
	 break;
      case 3: {
            if(j==0) return 0.0;

            if(simulation->identical_photons)
	    {
	       Lambda += exp(-simulation->alpha*simulation->Gamma*t);
               Lambda *= 2.0 * pow(cabs(simulation->e0[j]), 2.0);
	    }
	    else
	    {
	       Lambda += exp(-simulation->alpha1 * simulation->Gamma * t) * pow(cabs(simulation->e0_2[j]), 2.0);
	       Lambda += exp(-simulation->alpha2 * simulation->Gamma * t) * pow(cabs(simulation->e0_1[j]), 2.0);
	       double complex temp; 
	       temp = I*sqrt(simulation->alpha1*simulation->alpha2)*simulation->Gamma*conj(simulation->e0_2[j])*simulation->e0_1[j]
		      *cexp( (I*(simulation->k1-simulation->k2)-0.5*(simulation->alpha1+simulation->alpha2)*simulation->Gamma) * t)
		      /(simulation->k1-simulation->k2+0.5*I*(simulation->alpha1+simulation->alpha2)*simulation->Gamma);
	       Lambda += 2.0 * creal(temp);
	       Lambda *= simulation->A * simulation->A;
	    }
	 }
	 break;
      default: { //bad input
            fprintf(stderr, "%s: invalid option. Abort!\n", __func__);
            exit(EXIT_FAILURE);
         } 
   }

   int xmax = j + simulation->plus_a_index;
   //trapezoidal rule for x>=-a
   sum += 0.5 * pow(cabs(simulation->psi[j][simulation->minus_a_index]), 2.0);
   for(int i = simulation->minus_a_index+1; i<xmax; i++)
      sum += pow(cabs(simulation->psi[j][i]), 2.0);
   sum += 0.5 * pow(cabs(simulation->psi[j][xmax]), 2.0);
   Lambda += simulation->Delta * sum;

   return Lambda;  
}


inline void solver(int j, int i, grid * simulation)
{
           //points (i) right next to the 1st light cone and in tile B1, 
           //and (ii) right next to the 2nd light cone 
           //should be strictly zero under any circumstances
           if( ((j < simulation->nx) && (i==j+simulation->minus_a_index+1))
               || (i==j+simulation->plus_a_index+1) )
	       return;
               //continue; //do nothing, as psi is already zero initailized

           //free propagation (decay included)
           simulation->psi[j][i] = (1./simulation->Delta-0.25*W)*simulation->psi[j-1][i-1]   \
                                   -0.25*W*(simulation->psi[j-1][i]+simulation->psi[j][i-1]);
           
           //delay term: psi(x-2a, t-2a)theta(t-2a)
           if(j>simulation->nx)
               simulation->psi[j][i] += 0.5*simulation->Gamma*square_average(j-simulation->nx, i-simulation->nx, simulation);
   
           //left light cone No.1: psi(-x-2a, t-x-a)theta(x+a)theta(t-x-a)
           if( (i>simulation->minus_a_index) && (j-i>=-simulation->minus_a_index) )
           { 
               double on_light_cone = (j-i == -simulation->minus_a_index?0.5:1.0);
               simulation->psi[j][i] -= 0.5*simulation->Gamma*bar_average(j-(i-simulation->origin_index)-simulation->nx/2, \
                                        2*simulation->origin_index-i-simulation->nx+1, simulation)*on_light_cone; 
           }
   
           //left light cone No.2: -psi(-x, t-x-a)theta(x+a)theta(t-x-a)
           if( (i>simulation->minus_a_index) && (j-i>=-simulation->minus_a_index) )
           { 
               double on_light_cone = (j-i == -simulation->minus_a_index?0.5:1.0);
               simulation->psi[j][i] += 0.5*simulation->Gamma*bar_average(j-(i-simulation->origin_index)-simulation->nx/2, \
                                        2*simulation->origin_index-i+1, simulation)*on_light_cone; 
           }
    
           //right light cone No.1: psi(2a-x, t-x+a)theta(x-a)theta(t-x+a)
           if( (i>simulation->plus_a_index) && (j-i>=-simulation->plus_a_index) )
           { 
               double on_light_cone = (j-i == -simulation->plus_a_index?0.5:1.0);
               simulation->psi[j][i] -= 0.5*simulation->Gamma*bar_average(j-(i-simulation->origin_index)+simulation->nx/2, \
                                        2*simulation->origin_index-i+simulation->nx+1, simulation)*on_light_cone; 
           }
    
           //right light cone No.2: -psi(-x, t-x+a)theta(x-a)theta(t-x+a)
           if( (i>simulation->plus_a_index) && (j-i>=-simulation->plus_a_index) )
           { 
               double on_light_cone = (j-i == -simulation->plus_a_index?0.5:1.0);
               simulation->psi[j][i] += 0.5*simulation->Gamma*bar_average(j-(i-simulation->origin_index)+simulation->nx/2, \
                                        2*simulation->origin_index-i+1, simulation)*on_light_cone; 
           }
   
           //two-photon input: 2*( chi(x-t,-a-t, 0)-chi(x-t,a-t,0) )
           if( (simulation->init_cond == 1 || simulation->init_cond == 3) \
	       && j-i>=-simulation->minus_a_index ) //it's nonzero only when t-x-a>=0 
           { 
               double on_light_cone = (j-i == -simulation->minus_a_index?0.5:1.0);

               //shift +0.5 due to Taylor expansion at the center of square 
               simulation->psi[j][i] += sqrt(simulation->Gamma) * on_light_cone \
                                        * two_photon_input((i-simulation->origin_index)-j, -simulation->nx/2-j+0.5, simulation);

               if(j>simulation->nx)
               {
                   simulation->psi[j][i] -= sqrt(simulation->Gamma) * on_light_cone \
                                            * two_photon_input((i-simulation->origin_index)-j, simulation->nx/2-j+0.5, simulation);
               }
           }
   
           //prefactor
           simulation->psi[j][i] /= (1./simulation->Delta+0.25*W);


}
#endif

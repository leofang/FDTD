/*
 * Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
 *
 * This program is free software. It comes without any warranty,
 * to the extent permitted by applicable law. You can redistribute
 * it and/or modify it under the terms of the WTFPL, Version 2, as
 * published by Sam Hocevar. See the accompanying LICENSE file or
 * http://www.wtfpl.net/ for more details.
 */

#ifndef __GRID_H__
#define __GRID_H__

#include <stdlib.h>
#include <stdio.h>
#include <complex.h> 
#include "kv.h"

/* 
   Create a grid which stores the wavefunction and other relavant information.
   The layout of the grid should look like this:
  
                                                     i=(Nx+nx+1)
        array index i=0 1 2 3 ... ... nx (nx+1) ... ...   \ ... ...  (2Nx+nx+1)
                      | | | |          \   /               \                 \
   t=(Ny-1)*Delta  ^  % % % % ... % % % % @ @ @ @ @ @ @ @ @ @ ... @ @ @ @ @ @ @ 
   t=(Ny-2)*Delta  |  % % % % ... % % % % @ @ @ @ @ @ @ @ @ @ ... @ @ @ @ @ @ @ 
   t=(Ny-3)*Delta  |  % % % % ... % % % % @ @ @ @ @ @ @ @ @ @ ... @ @ @ @ @ @ @ 
                               .          @ @ @ @ @ @ @ @ @ @ ... @ @ @ @ @ @ @ 
                   j           .          @ @ @ @ @ @ @ @ @ @ ... @ @ @ @ @ @ @ 
                               .          @ @ @ @ @ @ @ @ @ @ ... @ @ @ @ @ @ @ 
   t= 2 * Delta    |  % % % % ... % % % % @ @ @ @ @ @ @ @ @ @ ... @ @ @ @ @ @ @ 
   t= 1 * Delta    |  % % % % ... % % % % @ @ @ @ @ @ @ @ @ @ ... @ @ @ @ @ @ @ 
   t= 0            |  % % % % ... % % % % * * * * * * * * * * ... * * * * * * *
                     /                 /   \                 \               /
     x=-(Nx+nx+1)*Delta   x=-(Nx+1)*Delta   x=-Nx*Delta      x=0     x=Nx*Delta
  
   where the t=0 line (symbol *) is given by the initial condition, and the left 
   area (symbol %) is given by the boundary condition, both of which should be
   given before the calculation starts. The rectangular area (symbol @) of size
   (2Nx+1)*(Ny-1) is the region where the wavefunction is to be solved.
  
   Note that in the layout above, the qubit at x=-a has index i=Nx+nx/2+1, and 
   its mirror image at x=+a has index i=Nx+3nx/2+1. This is why we need nx to 
   be an integer multiple of 2, and also require nx <= 2Nx to reach x>=a and to
   make the "boundary condition" well defined in x<=-a."               
*/


struct _grid 
{
   //spacetime parameters
   int nx;            // nx=2a/Delta
   int Nx;            // defined such that total grid points (to be solved) along x is 2Nx+1
   int Ntotal;        // total grid points along x, which is 2Nx+nx+2 
   int Ny;            // total grid points along t
   double Delta;      // grid size
   double Lx;         // Lx = 2Nx*Delta 
   double Ly;         // Ly = (Ny-1)*Delta
   int plus_a_index;  // array index for x=+a
   int minus_a_index; // array index for x=-a
   int origin_index;  // array index for x=0 

   //physics parameters
   double k;      // incident frequency (in units of 1/Delta)   
   double k1;     // incident frequency for photon #1 (in units of 1/Delta)   
   double k2;     // incident frequency for photon #2 (in units of 1/Delta)   
   double w0;     // qubit frequency  (in units of 1/Delta)
   double Gamma;  // qubit decay rate (in units of 1/Delta)
   double alpha;  // exponential tail (dimensionless)
   double alpha1; // exponential tail for photon #1 (dimensionless)
   double alpha2; // exponential tail for photon #2 (dimensionless)
   double A;      // normalization constant for two-photon exponential wavepacket

   //actual info on dynamics
   double complex * psit0;  //initial condition psi(x,0) (stored as psi0[x])
   double complex ** psi;   //wavefunction psi(x,t) to be computed (stored as psi[t][x])
   double complex ** psix0; //boundary condition psi(-L,0) (stored as psix0[t][x])
   double complex * e0;     //qubit wavefunction for I.C. e(0)=0 and an exponential wavepacket
   double complex * e0_1;   //qubit wavefunction for I.C. e(0)=0 and the exponential wavepacket #1
   double complex * e0_2;   //qubit wavefunction for I.C. e(0)=0 and the exponential wavepacket #2
   double complex * e1;     //qubit wavefunction for I.C. e(0)=1 and no incident wavepacket
   double complex * mu;     //mu(t) for calculating NM measures
   double * lambda;         //lambda(t) for calculating NM measures
   
   //auxiliary parameters
   int psit0_size;   //array size of psit0 
   int psi_x_size;   //array size of psi in x
   int psi_y_size;   //array size of psi in t
   int psix0_x_size; //array size of psix0 in x
   int psix0_y_size; //array size of psix0 in t

   //program options
   int save_chi;          //whether or not to save the two-photon wavefunction to file (default: no)
   int save_psi;          //whether or not to save the wavefunction to file (default: no)
   int save_psi_square_integral; //whether or not to save \int dx |psi(x,t)|^2 to file (default: no)
   int save_psi_binary;   //whether or not to save the wavefunction to binary file (default: no)
   int init_cond;         //the initial condition of the wavefunction (default: unspecified)
   int identical_photons; //whether or not the two photons are identical (default: yes; only effective for init_cond=3)
   size_t Tstep;          //for output of save_psi: save psi for every (Tstep+1) temporal steps
   int measure_NM;        //currently it means whether to save e0 and e1 or not //TODO: extend this part
   size_t Nth;            //number of threads

   //input parameters (stored for convenience)
   kvarray_t * parameters_key_value_pair;
};
typedef struct _grid grid;

double complex plane_wave_BC(int j, int i, grid * simulation);
double complex exponential_BC(int j, int i, grid * simulation);
double complex two_exponential_BC(int j, int i, grid * simulation);
void initial_condition(grid * simulation);
void boundary_condition(grid * simulation);
void initialize_psi(grid * simulation);
void sanity_check (grid * simulation);
void free_initial_boundary_conditions(grid * simulation);
void free_grid(grid * simulation);
grid * initialize_grid(const char * filename);
void print_initial_condition(grid * simulation);
void print_boundary_condition(grid * simulation);
void print_grid(grid * simulation);
void print_psi(grid * simulation);
void save_psi(grid * simulation, const char * filename, double (*part)(double complex));
void save_psi_binary(grid * simulation, const char * filename);
void save_chi(grid * simulation, const char * filename, double (*part)(double complex));
void save_psi_square_integral(grid * simulation, const char * filename);
void prepare_qubit_wavefunction(grid * simulation);
void initialize_e0(grid * simulation);
void initialize_e1(grid * simulation);
void calculate_normalization_const(grid * simulation);

// make W = (i*w0+Gamma/2) as a global variable
double complex W;

#endif

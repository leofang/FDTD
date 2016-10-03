#ifndef _FDTD_H
#define _FDTD_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "kv.h"

// Create a grid which stores the wavefunction and other relavant information.
// The layout of the grid should look like this:
//
//      array index i=0 1 2 3 ... ... nx (nx+1) ... ... ... ... ...  (2Nx+nx+1)
//                    | | | |          \   /                                 \
// t=(Ny-1)*Delta  ^  * * * * ... * * * * @ @ @ @ @ @ @ @ @ @ ... @ @ @ @ @ @ @ 
// t=(Ny-2)*Delta  |  * * * * ... * * * * @ @ @ @ @ @ @ @ @ @ ... @ @ @ @ @ @ @ 
// t=(Ny-3)*Delta  |  * * * * ... * * * * @ @ @ @ @ @ @ @ @ @ ... @ @ @ @ @ @ @ 
//                             .          @ @ @ @ @ @ @ @ @ @ ... @ @ @ @ @ @ @ 
//                 j           .          @ @ @ @ @ @ @ @ @ @ ... @ @ @ @ @ @ @ 
//                             .          @ @ @ @ @ @ @ @ @ @ ... @ @ @ @ @ @ @ 
// t= 2 * Delta    |  * * * * ... * * * * @ @ @ @ @ @ @ @ @ @ ... @ @ @ @ @ @ @ 
// t= 1 * Delta    |  * * * * ... * * * * @ @ @ @ @ @ @ @ @ @ ... @ @ @ @ @ @ @ 
// t= 0            |  * * * * ... * * * * * * * * * * * * * * ... * * * * * * *
//                   /                 /   \                                 /
//   x=-(Nx+nx+1)*Delta   x=-(Nx+1)*Delta   x=-Nx*Delta              x=Nx*Delta
//
// where the t=0 line (symbol *) is given by the initial condition, and the left 
// area (symbol *) is given by the boundary condition, both of which should 
// eventually be the input. The rectangular area (symbol @) of size (2Nx+1)*(Ny-1)
// is the region where the wavefunction is to be solved.



struct _grid 
{
   //spacetime parameters
   int nx;       // nx=2a/Delta
   int Nx;       // defined such that total grid points (to be solved) along x is 2Nx+1
   int Ntotal;   // total grid points along x, which is 2Nx+nx+2 
   int Ny;       // total grid points along t
   double Delta; // grid size
   double Lx;    // Lx = 2Nx*Delta 
   double Ly;    // Ly = (Ny-1)*Delta

   //physics parameters
   double k;     // incident frequency (in units of 1/Delta)   
   double w0;    // qubit frequency    (in units of 1/Delta)
   double Gamma; // qubit decay rate   (in units of 1/Delta) 

   //actual info on dynamics
   double * psit0;  //initial condition psi(x,0) (stored as psi0[x])
   double ** psi;   //wavefunction psi(x,t) to be computed (stored as psi[t][x])
   double ** psix0; //boundary condition psi(-L,0) (stored as psix0[t][x])
   
   //auxiliary parameters
   int psit0_size;   //array size of psit0 
   int psi_x_size;   //array size of psi in x
   int psi_y_size;   //array size of psi in t
   int psix0_x_size; //array size of psix0 in x
   int psix0_y_size; //array size of psix0 in t

   //input parameters (stored for convenience)
   kvarray_t * parameters_key_value_pair;
};
typedef struct _grid grid;


void initial_condition(grid * simulation);
void boundary_condition(grid * simulation);
void initialize_psi(grid * simulation);
void free_initial_boundary_conditions(grid * simulation);
void free_grid(grid * simulation);
grid * initialize_grid(const char * filename);
void print_initial_condition(grid * simulation);
void print_grid(grid * simulation);




#endif

#ifndef _FDTD_H
#define _FDTD_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "kv.h"

// create a grid which stores the wavefunction and other relavant information
//
//

struct _grid 
{
   //spacetime parameters
   int nx;    // nx=2a/Delta
   int Nx;    // defined such that total grid points along x is 2Nx+1
   int Ny;    // total grid points along t
   double Delta; // grid size
   double Lx;    // Lx = 2Nx*Delta 
   double Ly;    // Ly = Ny*Delta

   //physics parameters
   double k;     // incident frequency
   double w0;    // qubit frequency
   double Gamma; // qubit decay rate

   //actual info on dynamics
   double * psit0; //initial condition psi(x,0) (stored as psi0[x])
   double ** psi;  //wavefunction psi(x,t) to be computed (stored as psi[t][x])
   double ** psix0; //boundary condition psi(-L,0) (stored as psix0[t][x])

   //input parameters (stored for convenience)
   kvarray_t * parameters_key_value_pair;
};
typedef struct _grid grid;


void initialCondition(grid * simulation);

void boundaryCondition(grid * simulation);

void initializePsi(grid * simulation);

void freeGrid(grid * simulation);

grid * initializeGrid(const char * filename);

void printInitialCondition(grid * simulation);

void printGrid(grid * simulation);




#endif

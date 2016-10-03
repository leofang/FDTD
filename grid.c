#include <stdlib.h>
#include "kv.h"
#include "grid.h"


void initialCondition(grid * simulation)
{
    simulation->psit0 = calloc(2*simulation->Nx+1, sizeof(double));

    //Finish it!
}


void boundaryCondition(grid * simulation)
{
    simulation->psix0 = malloc( simulation->Ny*sizeof(double) );
    for(int j=0; j<simulation->Ny; j++)
        simulation->psix0[j] = calloc(3, sizeof(double)); //need to know the boundary conditions in first 3 columns due to the delay term

    //Finish it!!
}


void initializePsi(grid * simulation)
{
    simulation->psi = malloc( simulation->Ny*sizeof(double) );
    for(int j=0; j<simulation->Ny; j++)
    {
        simulation->psi[j] = calloc( (2*simulation->Nx+1)+3, sizeof(double) ); //include the boundary conditions
        for(int i=0; i<3; i++)
           simulation->psi[j][i] = simulation->psix0[j][i];
    }
}


void freeGrid(grid * simulation)
{
    freeKVs(simulation->parameters_key_value_pair);
    free(simulation->psit0);
    for(int j=0; j<simulation->Ny; j++)
    {
       free(simulation->psi[j]);
       free(simulation->psix0[j]);
    }
    free(simulation->psi);
    free(simulation->psix0);

    free(simulation);
}


grid * initializeGrid(const char * filename)
{
   grid * FDTDsimulation = malloc(sizeof(*FDTDsimulation));
   if(!FDTDsimulation)
   {
      perror("initializeGrid: could not allocate memory. Abort!\n");
      exit(EXIT_FAILURE);
   }
   
   //read from input
   FDTDsimulation->parameters_key_value_pair = readKVs(filename);

   //initialize from the input parameters
   FDTDsimulation->nx    = atoi(lookupValue(FDTDsimulation->parameters_key_value_pair, "nx"));
   FDTDsimulation->Nx    = atoi(lookupValue(FDTDsimulation->parameters_key_value_pair, "Nx"));
   FDTDsimulation->Ny    = atoi(lookupValue(FDTDsimulation->parameters_key_value_pair, "Ny"));
   FDTDsimulation->Delta = strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "Delta"), NULL);
   FDTDsimulation->k     = strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "k"), NULL);
   FDTDsimulation->w0    = strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "w0"), NULL);
   FDTDsimulation->Gamma = strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "Gamma"), NULL);
   FDTDsimulation->Lx    = 2 * FDTDsimulation->Nx * FDTDsimulation->Delta;
   FDTDsimulation->Ly    = (FDTDsimulation->Ny-1) * FDTDsimulation->Delta;
   
   //initialize arrays
   initialCondition(FDTDsimulation);
   boundaryCondition(FDTDsimulation);
   initializePsi(FDTDsimulation);

   return FDTDsimulation;
}


void printInitialCondition(grid * simulation)
{
    printf("t=0   ");
    for(int i=0; i<2*simulation->Lx+1; i++)
        printf("%.3f ", simulation->psit0[i]);
    printf("\n      ");
    for(int i=-simulation->Nx; i<=simulation->Nx; i++)
        printf("x=%.2f ", i*simulation->Delta);
    printf("\n");
}


void printGrid(grid * simulation)
{
    printf("nx = %d\n", simulation->nx); 
    printf("Nx = %d\n", simulation->Nx); 
    printf("Ny = %d\n", simulation->Ny); 
    printf("Delta = %.3f\n", simulation->Delta); 
    printf("k = %.3f\n", simulation->k); 
    printf("w0 = %.3f\n", simulation->w0); 
    printf("Gamma = %.3f\n", simulation->Gamma); 
    printf("Lx = %.3f\n", simulation->Lx); 
    printf("Ly = %.3f\n", simulation->Ly); 
}

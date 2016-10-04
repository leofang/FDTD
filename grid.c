#include <stdlib.h>
#include "kv.h"
#include "grid.h"


void initial_condition(grid * simulation)
{// the initial condition is given in-between x/Delta = [-Nx, Nx] for simplicity

    simulation->psit0 = calloc(2*simulation->Nx+1, sizeof(*simulation->psit0));
    if(!simulation->psit0)
    { 
        perror("initialCondition: cannot allocate memory. Abort!\n");
        exit(EXIT_FAILURE);
    }
    simulation->psit0_size = 2*simulation->Nx+1;

    for(int i=0; i<simulation->psit0_size; i++)
    {
        simulation->psit0[i] = i;
    }

    //Finish it!
}


void boundary_condition(grid * simulation)
{// the boundary conditions is given for first nx+1 columns with 
 // x/Delta=[-(Nx+nx+1),-(Nx+1)] due to the delay term

    simulation->psix0 = malloc( simulation->Ny*sizeof(*simulation->psix0) );
    if(!simulation->psix0)
    { 
        perror("boundaryCondition: cannot allocate memory. Abort!\n");
        exit(EXIT_FAILURE);
    }
    simulation->psix0_x_size = simulation->nx+1;
    simulation->psix0_y_size = 0;

    for(int j=0; j<simulation->Ny; j++)
    {
        simulation->psix0[j] = calloc(simulation->nx+1, sizeof(*simulation->psix0[j])); 
        if(!simulation->psix0[j])
        { 
            fprintf(stderr, "boundaryCondition: cannot allocate memory at t=%d*Delta. Abort!\n", j);
            exit(EXIT_FAILURE);
        }
        simulation->psix0_y_size++;
    }

    for(int j=0; j<simulation->psix0_y_size; j++)
    {
        for(int i=0; i<simulation->psix0_x_size; i++)
            simulation->psix0[j][i] = i*j;
    }

    //Finish it!!
}


void initialize_psi(grid * simulation)
{
    simulation->psi = malloc( simulation->Ny*sizeof(*simulation->psi) );
    if(!simulation->psi)
    { 
        perror("initializePsi: cannot allocate memory. Abort!\n");
        exit(EXIT_FAILURE);
    }
    simulation->psi_x_size = simulation->Ntotal;
    simulation->psi_y_size = 0;
    for(int j=0; j<simulation->Ny; j++)
    {
        simulation->psi[j] = calloc( simulation->Ntotal, sizeof(*simulation->psi[j]) );
        if(!simulation->psi[j])
        { 
            fprintf(stderr, "initializePsi: cannot allocate memory at t=%d*Delta. Abort!\n", j);
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
      perror("initializeGrid: could not allocate memory. Abort!\n");
      exit(EXIT_FAILURE);
   }
   
   //read from input
   FDTDsimulation->parameters_key_value_pair = readKVs(filename);

   //initialize from the input parameters
   FDTDsimulation->nx    = atoi(lookupValue(FDTDsimulation->parameters_key_value_pair, "nx"));
   FDTDsimulation->Nx    = atoi(lookupValue(FDTDsimulation->parameters_key_value_pair, "Nx"));
   FDTDsimulation->Ntotal= 2 * FDTDsimulation->Nx + FDTDsimulation->nx + 2;
   FDTDsimulation->Ny    = atoi(lookupValue(FDTDsimulation->parameters_key_value_pair, "Ny"));
   FDTDsimulation->Delta = strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "Delta"), NULL);
   FDTDsimulation->k     = strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "k"), NULL);
   FDTDsimulation->w0    = strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "w0"), NULL);
   FDTDsimulation->Gamma = strtod(lookupValue(FDTDsimulation->parameters_key_value_pair, "Gamma"), NULL);
   FDTDsimulation->Lx    = 2 * FDTDsimulation->Nx * FDTDsimulation->Delta;
   FDTDsimulation->Ly    = (FDTDsimulation->Ny-1) * FDTDsimulation->Delta;
   
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

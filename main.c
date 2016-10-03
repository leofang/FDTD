#include "grid.h"
#include "kv.h"


int main(int argc, char **argv)
{
   if(argc != 2)
   {
      fprintf(stderr, "Usage: ./FDTD input_parameters\n");
      exit(EXIT_FAILURE);
   }
 
   grid * simulation = initialize_grid(argv[1]);

   for(int j=0; j<simulation->Ny; j++)
   {
       for(int i=simulation->nx+1; i<simulation->Ntotal; i++)
       {
           
            //delay term psi(x-2a, t-2a)theta(t-2a)
            if(j>simulation->nx)
               simulation->psi[j][i] += square_average(j-simulation->nx, i-simulation->nx, simulation);

       }
   }


//   print_grid(simulation);
//   print_initial_condition(simulation);


   free_grid(simulation);

   return EXIT_SUCCESS;
}

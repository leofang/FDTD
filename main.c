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

   print_grid(simulation);
//   print_initial_condition(simulation);


   free_grid(simulation);

   return EXIT_SUCCESS;
}

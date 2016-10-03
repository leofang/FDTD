#include "grid.h"
#include "kv.h"


int main(int argc, char **argv)
{
   if(argc != 2)
   {
      fprintf(stderr, "Usage: ./FDTD input_parameters\n");
      exit(EXIT_FAILURE);
   }
 
   grid * simulation = initializeGrid(argv[1]);

   printGrid(simulation);
   printInitialCondition(simulation);


   freeGrid(simulation);

   return EXIT_SUCCESS;
}

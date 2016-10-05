#include "grid.h"
#include "kv.h"
#include "dynamics.h"


int main(int argc, char **argv)
{
   if(argc != 2)
   {
      fprintf(stderr, "Usage: ./FDTD input_parameters\n");
      exit(EXIT_FAILURE);
   }
 
   grid * simulation = initialize_grid(argv[1]);

   // W = (i*w0+Gamma/2)
   complex W = simulation->w0*I+0.5*simulation->Gamma;

   for(int j=1; j<simulation->Ny; j++) //start from t=1*Delta
   {
       for(int i=simulation->nx+1; i<simulation->Ntotal; i++) //start from x=-Nx*Delta
       {
         //if i is NOT next to the qubits   
 
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
                                         2*simulation->origin_index-i-simulation->nx, simulation)*on_light_cone; //CHECK!!!!
            }

            //left light cone No.2: -psi(-x, t-x-a)theta(x+a)theta(t-x-a)
 
            //right light cone No.1: psi(2a-x, t-x+a)theta(x-a)theta(t-x+a)
 
            //right light cone No.2: -psi(-x, t-x+a)theta(x-a)theta(t-x+a)

            //two-photon input:

            //prefactor
            simulation->psi[j][i] /= (1./simulation->Delta+0.25*W);

         //if i IS next to qubits (so x=-a+Delta or x=+a+Delta)
 
       }
   }


//   print_grid(simulation);
//   print_initial_condition(simulation);


   free_grid(simulation);

   return EXIT_SUCCESS;
}

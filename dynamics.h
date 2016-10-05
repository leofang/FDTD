#ifndef __DYNAMICS_H__
#define __DYNAMICS_H__

// this function computes the average of 4 points that form a square
// j and i are the array indices psi[j][i] of the upper right corner
complex square_average(int j, int i, grid * simulation)
{
   if(j<1) //everything is zero below t=0
      return 0;

   if(i<1) //beyond the boundary x=-(Nx+nx+1)*Delta
   {
      fprintf(stderr, "Error in square_average: beyond left boundary. Abort!\n");
      exit(EXIT_FAILURE);
   }

   if(i>simulation->Ntotal-1) //beyond the boundary x=Nx*Delta
   {
      fprintf(stderr, "Error in square_average: beyond right boundary. Abort!\n");
      exit(EXIT_FAILURE);
   }

   complex square = 0;

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
complex bar_average(int j, int i, grid * simulation)
{
   if(j<0) //everything is zero below t=0
      return 0;

   if(i<1) //beyond the boundary x=-(Nx+nx+1)*Delta
   {
      fprintf(stderr, "Error in bar_average: beyond left boundary. Abort!\n");
      exit(EXIT_FAILURE);
   }

   if(i>simulation->Ntotal-1) //beyond the boundary x=Nx*Delta
   {
      fprintf(stderr, "Error in bar_average: beyond right boundary. Abort!\n");
      exit(EXIT_FAILURE);
   }

   return (simulation->psi[j][i]+simulation->psi[j][i-1])/2.;
}



#endif

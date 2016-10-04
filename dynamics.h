



// this function computes the average of 4 points that form a square
// j and i are the array indices psi[j][i] of the upper right corner
double square_average(int j, int i, grid * simulation)
{
   double square = 0;

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
double bar_average(int j, int i, grid * simulation)
{
   return (simulation->psi[j][i]+simulation->psi[j][i-1])/2.;
}





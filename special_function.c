/*
 * Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
 *
 * This program is free software. It comes without any warranty,
 * to the extent permitted by applicable law. You can redistribute
 * it and/or modify it under the terms of the WTFPL, Version 2, as
 * published by Sam Hocevar. See the accompanying LICENSE file or
 * http://www.wtfpl.net/ for more details.
 */

#include "special_function.h"

// Compute the normalized incomplete Gamma function P(n,x). which is 
// defined as P(n,x)=gamma(n,x)/(n-1)! and satisfies P(n,x)+Q(n,x)=1. 
// P(n,x) can be called in Mathematica by GammaRegularized[n, 0, x].
//
// The following implementation is based on Numerical Recipes Ch.6.2, 3rd Ed,
// where two different but complementary approaches are required:
//    1. the infinite series representation (Eq.6.2.5)
//    2. the continued fraction representation (Eq.6.2.7)
// The combination hopefully covers the complex plane as much as possible.
// See also ASA032, ASA147 and ASA239.
//
// Apr 25 update:
// To deal with real negative argument (x<0), an additional procedure 
// is added based on Eq.2.5 & Eq.2.28 in arXiv:1608.04152.
//
// Note that
// 1. the normalization factor (n-1)!.
// 2. the finite sum, Eq.10.70, in P.661, Arfken 5th ed, which was
//    implemented in the earlier version, is not useful as it subtracts two
//    nearly same numbers.
double complex incomplete_gamma(int n, double complex x)
{
   if(n<=0)
   {
      fprintf(stderr, "Error in %s: n must >= 1. Abort!\n", __func__);
      exit(EXIT_FAILURE);
   }

   //gamma(n>=1, x=0)=0, no need to do real computation
   if(x==0)
      return 0;

   if(fabs(cimag(x))<DBL_EPSILON) 
   {//real negative argument
      if(creal(x)<-50)
      {
//         printf("Poincare expansion used...\n");
         x = -x; //make x>0
         double complex prefactor = pow(-1, n)*cexp(x-lgamma(n));
         double complex temp = 0.;
         double complex sum = 0.;

         for(int i=0; ; i++)
         {
            temp = Pochhammer(1-n, i) * cpow(x, n-i-1);
            sum += temp;

            if(cabs(temp) < cabs(sum)*DBL_EPSILON) //stop the sum when temp is significantly smaller than sum
            {
               if(!isnan(cabs(prefactor*sum)))
                  return prefactor*sum;
               else
                  break;
            }
         }
      }
      else
      {
//         printf("Series expansion for gamma* used...\n");
         x = -x; //make x>0
         double complex prefactor = pow(-1, n)*cexp(n*clog(x)-lgamma(n));
         double complex temp = 0.;
         double complex sum = 0.;

         for(int i=0; ; i++)
         {
            temp = cexp(i*clog(x)-lgamma(i+1))/(n+i);
            sum += temp;

            if(cabs(temp) < cabs(sum)*DBL_EPSILON) //stop the sum when temp is significantly smaller than sum
            {
               if(!isnan(cabs(prefactor*sum)))
                  return prefactor*sum;
               else
                  break;
            }
         }
      }
   }
   else if(cabs(x)<n+1)
//   else if( cabs(x)>=n+1 || (-50 <= creal(x) && creal(x) < 0) )
   {//compute the infinite sum
//      printf("series expansion used...\n");
      double complex prefactor = cexp(n*clog(x)-x-lgamma(n)); //exp(-x)*x^n/(n-1)!
      double complex sum = 1.0/n;
      double complex temp = 1.0/n;
      for(int i=1; ; i++)
      {
         temp *= x/(double)(n+i);
         sum += temp;

         if(cabs(temp) < cabs(sum)*DBL_EPSILON) //stop the sum when temp is significantly smaller than sum
//         //stop the sum when temp is significantly smaller than sum
//         if( (fabs(creal(temp)) < fabs(creal(sum))*DBL_EPSILON) && (fabs(cimag(temp)) < fabs(cimag(sum))*DBL_EPSILON) )
         {
            if(!isnan(cabs(prefactor*sum)))
               return prefactor*sum;
            else
	       break;
         }
      }
   }
   else
   {//use the continued fraction frac=(a1/b1+)(a2/b2+)(a3/b3+)...
    //the notation follows Ch.5.2 of Numerical Recipes 3rd Ed.
//      printf("continued fraction used...\n");
      double complex prefactor = cexp(n*clog(x)-x-lgamma(n)); //exp(-x)*x^n/(n-1)!
      double complex b = x+1.0-n;        //b1
      double complex c = INFINITY;       //C1
      double complex d = 1 / b;          //D1=a1/b1
      double complex frac = d;           //frac_1
      double complex del = 0;            //Delta_i
      for(int i=1; ; i++)
      {
         double complex a = -i*(i-n);    //calculate a_{i+1}
	 b += 2.0;                //calculate b_{i+1}
	 d = b + a * d;           //calculate 1/D_{i+1}
	 if(cabs(d)<DBL_MIN) 
	    d=DBL_MIN;
         #if (__APPLE__ && __MACH__)
	 //on Mac 1/(inf+0I) gives nan+nanI...
	 c = (i==1?b: b + a / c); //calculate C_{i+1}
         #else
	 c = b + a / c;           //calculate C_{i+1}
         #endif
	 if(cabs(c)<DBL_MIN)
	    c=DBL_MIN;
	 d = 1.0 / d;             //calculate D_{i+1}
	 del = d * c;             
	 frac *= del;             //frac_{i+1} = frac_i*C_{i+1}*D_{i+1}
	 if(cabs(del-1.0) < DBL_EPSILON)
	 {
            if(!isnan(cabs(prefactor*frac)))
               return 1.0-prefactor*frac;
            else
	       break;
	 }
      }
   }

   fprintf(stderr, "%s: NaN is produced (at n=%i and x=%.6f+%.6fI). Abort!\n", __func__, n, creal(x), cimag(x));
   exit(EXIT_FAILURE);
}


//return the Pochhammer's symbol (a)_n = a*(a+1)*(a+2)*...*(a+n-1)
//TODO: add a check to skip calculation if n is large enough for a<0.
double Pochhammer(double a, int n)
{
   if(n < 0)
   {
      fprintf(stderr, "%s: negative argument in n=%i. Abort!\n", __func__, n);
      exit(EXIT_FAILURE);
   }

   if(n == 0)
      return 1.;

   double result = a;
   for(int i=1; i<n; i++)
      result *= (a+i);
   
   return result;
}


//calculate P(n,x)e^y
double complex incomplete_gamma_e(int n, double complex x, double complex y)
{
   if(n<=0)
   {
      fprintf(stderr, "Error in %s: n must >= 1. Abort!\n", __func__);
      exit(EXIT_FAILURE);
   }

   //gamma(n>=1, x=0)=0, no need to do real computation
   if(x==0)
      return 0;

   if(fabs(cimag(x))<DBL_EPSILON) 
   {//real negative argument
      if(creal(x)<-50)
      {
//         printf("Poincare expansion used...\n");
         x = -x; //make x>0
         double complex prefactor = pow(-1, n)*cexp(x+y-lgamma(n));
         double complex temp = 0.;
         double complex sum = 0.;

         for(int i=0; ; i++)
         {
            temp = Pochhammer(1-n, i) * cpow(x, n-i-1);
            sum += temp;

            if(cabs(temp) < cabs(sum)*DBL_EPSILON) //stop the sum when temp is significantly smaller than sum
            {
               if(!isnan(cabs(prefactor*sum)))
                  return prefactor*sum;
               else
                  break;
            }
         }
      }
      else
      {
//         printf("Series expansion for gamma* used...\n");
         x = -x; //make x>0
         double complex prefactor = pow(-1, n)*cexp(y+n*clog(x)-lgamma(n));
         double complex temp = 0.;
         double complex sum = 0.;

         for(int i=0; ; i++)
         {
            temp = cexp(i*clog(x)-lgamma(i+1))/(n+i);
            sum += temp;

            if(cabs(temp) < cabs(sum)*DBL_EPSILON) //stop the sum when temp is significantly smaller than sum
            {
               if(!isnan(cabs(prefactor*sum)))
                  return prefactor*sum;
               else
                  break;
            }
         }
      }
   }
   else if(cabs(x)<n+1)
//   else if( cabs(x)>=n+1 || (-50 <= creal(x) && creal(x) < 0) )
   {//compute the infinite sum
//      printf("series expansion used...\n");
      double complex prefactor = cexp(n*clog(x)-x-lgamma(n)+y); //exp(-x)*x^n/(n-1)!
      double complex sum = 1.0/n;
      double complex temp = 1.0/n;
      for(int i=1; ; i++)
      {
         temp *= x/(double)(n+i);
         sum += temp;

         if(cabs(temp) < cabs(sum)*DBL_EPSILON) //stop the sum when temp is significantly smaller than sum
//         //stop the sum when temp is significantly smaller than sum
//         if( (fabs(creal(temp)) < fabs(creal(sum))*DBL_EPSILON) && (fabs(cimag(temp)) < fabs(cimag(sum))*DBL_EPSILON) )
         {
            if(!isnan(cabs(prefactor*sum)))
               return prefactor*sum;
            else
	       break;
         }
      }
   }
   else
   {//use the continued fraction frac=(a1/b1+)(a2/b2+)(a3/b3+)...
    //the notation follows Ch.5.2 of Numerical Recipes 3rd Ed.
//      printf("continued fraction used...\n");
      double complex prefactor = cexp(n*clog(x)-x-lgamma(n)); //exp(-x)*x^n/(n-1)!
      double complex b = x+1.0-n;        //b1
      double complex c = INFINITY;       //C1
      double complex d = 1 / b;          //D1=a1/b1
      double complex frac = d;           //frac_1
      double complex del = 0;            //Delta_i
      for(int i=1; ; i++)
      {
         double complex a = -i*(i-n);    //calculate a_{i+1}
	 b += 2.0;                //calculate b_{i+1}
	 d = b + a * d;           //calculate 1/D_{i+1}
	 if(cabs(d)<DBL_MIN) 
	    d=DBL_MIN;
         #if (__APPLE__ && __MACH__)
	 //on Mac 1/(inf+0I) gives nan+nanI...
	 c = (i==1?b: b + a / c); //calculate C_{i+1}
         #else
	 c = b + a / c;           //calculate C_{i+1}
         #endif
	 if(cabs(c)<DBL_MIN)
	    c=DBL_MIN;
	 d = 1.0 / d;             //calculate D_{i+1}
	 del = d * c;             
	 frac *= del;             //frac_{i+1} = frac_i*C_{i+1}*D_{i+1}
	 if(cabs(del-1.0) < DBL_EPSILON)
	 {
            if(!isnan(cabs(prefactor*frac)))
               return (1.0-prefactor*frac)*cexp(y);
            else
	       break;
	 }
      }
   }

   fprintf(stderr, "%s: NaN is produced (at n=%i and x=%.6f+%.6fI). Abort!\n", __func__, n, creal(x), cimag(x));
   exit(EXIT_FAILURE);
}

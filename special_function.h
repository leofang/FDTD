/*
 * Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
 *
 * This program is free software. It comes without any warranty,
 * to the extent permitted by applicable law. You can redistribute
 * it and/or modify it under the terms of the WTFPL, Version 2, as
 * published by Sam Hocevar. See the accompanying LICENSE file or
 * http://www.wtfpl.net/ for more details.
 */

#ifndef __SPECIAL_FUNCTION_H__
#define __SPECIAL_FUNCTION_H__


#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
//#include <complex.h>
//#include <math.h>
#include <float.h> //for DBL_EPSILON ~ 2.2E-16


// Compute the incomplete Gamma function gamma(n,x). Note that gamma(n,x)
// is different from Gamma(n,x) (gamma(n,x)+Gamma(n,x)=1), and can be
// called in Mathematica by GammaRegularized[n, 0, x].
//
// The following implementation is based on Numerical Recipes Ch.6.2, 3rd Ed,
// where two different but complementary approaches are required:
//    1. the infinite series representation (Eq.6.2.5)
//    2. the continued fraction representation (Eq.6.2.7)
// The combination hopefully covers the complex plane as much as possible.
// See also ASA032, ASA147 and ASA239.
//
// Note that
// 1. the normalization factor (n-1)!.
// 2. the finite sum, Eq.10.70, in P.661, Arfken 5th ed, which was
//    implemented in the earlier version, is not useful as it subtracts two
//    nearly same numbers.
complex incomplete_gamma(int n, complex x)
{
   if(n<=0)
   {
      fprintf(stderr, "Error in %s: n must >= 1. Abort!\n", __func__);
      exit(EXIT_FAILURE);
   }

   //gamma(n>=1, x=0)=0, no need to do real computation
   if(x==0)
      return 0;

   //compute exp(-x)*x^n/(n-1)!
   complex prefactor = cexp(n*clog(x)-x-lgamma(n));

   if(cabs(x)<n+1)
   {//compute the infinite sum
      complex sum = 1.0/n;
      complex temp = 1.0/n;
      for(int i=1; ; i++)
      {
         temp *= x/(double)(n+i);
         sum += temp;

         if(cabs(temp) < cabs(sum)*DBL_EPSILON) //stop the sum when temp is significantly smaller than sum
         {
            if(!isnan(prefactor*sum))
               return prefactor*sum;
            else
	       break;
         }
      }
   }
   else
   {//use the continued fraction frac=(a1/b1+)(a2/b2+)(a3/b3+)...
    //the notation follows Ch.5.2 of Numerical Recipes 3rd Ed.
      complex b = x+1.0-n;        //b1
      complex c = INFINITY;       //C1
      complex d = 1 / b;          //D1=a1/b1
      complex frac = d;           //frac_1
      complex del = 0;            //Delta_i
      for(int i=1; ; i++)
      {
         complex a = -i*(i-n);    //calculate a_{i+1}
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
            if(!isnan(prefactor*frac))
               return 1.0-prefactor*frac;
            else
	       break;
	 }
      }
   }

   fprintf(stderr, "%s: NaN is produced (at n=%i and x=%.6f+%.6fI). Abort!\n", __func__, n, creal(x), cimag(x));
   exit(EXIT_FAILURE);
}


#endif

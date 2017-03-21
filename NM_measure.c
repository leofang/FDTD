/*
 * Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
 *
 * This program is free software. It comes without any warranty,
 * to the extent permitted by applicable law. You can redistribute
 * it and/or modify it under the terms of the WTFPL, Version 2, as
 * published by Sam Hocevar. See the accompanying LICENSE file or
 * http://www.wtfpl.net/ for more details.
 */

#include "NM_measure.h"
#include <string.h>
#include "special_function.h"

// This function returns the qubit wavefunction e0(t) in the single-excitation sector
// subject to e0(0)=0 and an exponential wavepacket
double complex e0(int j, grid * simulation)
{
    if(j<=0) return 0;

    double t         = j*simulation->Delta;
    double td        = simulation->nx*simulation->Delta;
    double w0        = simulation->w0;
    double Gamma     = simulation->Gamma;
    double alpha     = simulation->alpha;
    double complex K = I*simulation->k + 0.5*alpha*Gamma;
    double complex W = I*w0 + 0.5*Gamma;
    double complex p = -I*(K - W);

    double complex e_t = sqrt(0.5*alpha)*Gamma*(cexp(-W*t)-cexp(-K*t))/p;
    double complex sum = 0;

    for(int n=1; n<=(j/simulation->nx); n++)
    {
        double complex temp = ( cexp( n*log(t-n*td) - W*(t-n*td) - lgamma(n+1) ) \
                       - (I*K+w0) * incomplete_gamma(n+1, -I*p*(t-n*td)) * cexp( n*clog(I) - (n+1)*clog(p) - K*(t-n*td) ) );

        temp *= pow(0.5*Gamma, n-0.5);

	// based on my observation, the wavefunction should converge very fast, 
	// so one can just cut the summation off if the precision is reached.
	// this also helps prevent some overflow issue a bit.
        if( cabs(temp) < DBL_EPSILON*cabs(sum) || isnan(cabs(temp)) )
           break;
        else
	   sum += temp;
    }
    e_t -= I*sqrt(alpha*Gamma)*sum;
    e_t *= cexp(-0.5*I*(simulation->k)*td); //TODO: this phase factor can be eliminated by absorbing into the wavepacket

    if(!isnan(cabs(e_t)))
       return e_t;
    else
    {
       fprintf(stderr, "%s: NaN is produced (at j=%i). Abort!\n", __func__, j);
       exit(EXIT_FAILURE);
    }
}


// This function returns the qubit wavefunction e1(t) in the single-excitation sector
// subject to e1(0)=1 and no incidet wavepacket (spontaneous emission)
double complex e1(int j, grid * simulation)
{
    if(j<0) return 0;
    if(j==0) return 1.;

    double t         = j*simulation->Delta;
    double td        = simulation->nx*simulation->Delta;
    double w0        = simulation->w0;
    double Gamma     = simulation->Gamma;
    double complex W = I*w0 + 0.5*Gamma;

    double complex sum = 0;
    for(int n=1; n<=(j/simulation->nx); n++)
    {
        double complex temp = exp(-lgamma(n+1)) * cpow(0.5*Gamma*cexp(W*td)*(t-n*td), n);

	// based on my observation, the wavefunction should converge very fast, 
	// so one can just cut the summation off if the precision is reached.
	// this also helps prevent some overflow issue a bit.
        if( cabs(temp) < DBL_EPSILON*cabs(sum) || isnan(cabs(temp)) )
           break;
        else
	   sum += temp;
    }
    double complex e_t = cexp(-W*t)*(1.0+sum);

    if(!isnan(cabs(e_t)))
       return e_t;
    else
    {
       fprintf(stderr, "%s: NaN is produced (at j=%i). Abort!\n", __func__, j);
       exit(EXIT_FAILURE);
    }
}


void save_e0(grid * simulation, const char * filename, double (*part)(double complex))
{
    char * str = strdup(filename);
    str = realloc(str, (strlen(filename)+14)*sizeof(char) );

    if(part == &creal)
       strcat(str, ".re_e0.out");
    else if(part == &cimag)
       strcat(str, ".im_e0.out");
    else if(part == &cabs)
       strcat(str, ".abs_e0.out");
    else
    {
       fprintf(stderr, "%s: Warning: default filename is used.\n", __func__);
       strcat(str, ".out");
    }

    FILE * f = fopen(str, "w");

    for(int j=0; j<simulation->Ny; j++)
    {
       fprintf( f, "%.10g\n", part(simulation->e0[j]) );
    }

    fclose(f);
    free(str);
}


void save_e1(grid * simulation, const char * filename, double (*part)(double complex))
{
    char * str = strdup(filename);
    str = realloc(str, (strlen(filename)+14)*sizeof(char) );

    if(part == &creal)
       strcat(str, ".re_e1.out");
    else if(part == &cimag)
       strcat(str, ".im_e1.out");
    else if(part == &cabs)
       strcat(str, ".abs_e1.out");
    else
    {
       fprintf(stderr, "%s: Warning: default filename is used.\n", __func__);
       strcat(str, ".out");
    }

    FILE * f = fopen(str, "w");

    for(int j=0; j<simulation->Ny; j++)
    {
       fprintf( f, "%.10g\n", part(simulation->e1[j]) );
    }

    fclose(f);
    free(str);
}

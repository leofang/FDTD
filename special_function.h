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

double complex incomplete_gamma(int n, double complex x);
double complex incomplete_gamma_e(int n, double complex x, double complex y);
double Pochhammer(double a, int n);

#endif

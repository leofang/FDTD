/*
 * Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
 *
 * This program is free software. It comes without any warranty,
 * to the extent permitted by applicable law. You can redistribute
 * it and/or modify it under the terms of the WTFPL, Version 2, as
 * published by Sam Hocevar. See the accompanying LICENSE file or
 * http://www.wtfpl.net/ for more details.
 */

#ifndef __DYNAMICS_H__
#define __DYNAMICS_H__

#include "grid.h"


//a collection of convenient functions

double complex square_average(int j, int i, grid * simulation);
double complex bar_average(int j, int i, grid * simulation);
double complex two_photon_input(int i1, int i2, grid * simulation);
double complex one_photon_exponential(int i, grid * simulation);

#endif

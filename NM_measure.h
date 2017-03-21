/*
 * Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
 *
 * This program is free software. It comes without any warranty,
 * to the extent permitted by applicable law. You can redistribute
 * it and/or modify it under the terms of the WTFPL, Version 2, as
 * published by Sam Hocevar. See the accompanying LICENSE file or
 * http://www.wtfpl.net/ for more details.
 */

#ifndef __NM_MEASURE_H__
#define __NM_MEASURE_H__

#include "grid.h"

double complex e0(int j, grid * simulation);
double complex e1(int j, grid * simulation);
void save_e0(grid * simulation, const char * filename, double (*part)(double complex));
void save_e1(grid * simulation, const char * filename, double (*part)(double complex));

#endif

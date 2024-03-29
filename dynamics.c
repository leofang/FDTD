/*
 * Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
 *
 * This program is free software. It comes without any warranty,
 * to the extent permitted by applicable law. You can redistribute
 * it and/or modify it under the terms of the WTFPL, Version 2, as
 * published by Sam Hocevar. See the accompanying LICENSE file or
 * http://www.wtfpl.net/ for more details.
 */

#include "dynamics.h"

// Update Jun 21, 2017:
// The purpose of putting function declarations here is to use the C99 inline specifier; 
// see the note in dynamics.h.

//a collection of inline functions 
double complex square_average(int j, int i, grid * simulation);
double complex bar_average(int j, int i, grid * simulation);
double complex two_photon_input(double x1, double x2, grid * simulation);
double complex one_photon_exponential(double x, double k, double alpha, grid * simulation);
double complex chi(int j, int x1, int x2, grid * simulation);
double psi_square_integral(int j, grid * simulation);
double chi_square_double_integral(int j, grid * simulation);
double chi_square_single_integral(int j, int x1, grid * simulation);
double complex chi_square_single_integral_mirrored(int j, int x1, grid * simulation);
void solver(int j, int i, grid * simulation);

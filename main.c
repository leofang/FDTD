/*
 * Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
 *
 * This program is free software. It comes without any warranty,
 * to the extent permitted by applicable law. You can redistribute
 * it and/or modify it under the terms of the WTFPL, Version 2, as
 * published by Sam Hocevar. See the accompanying LICENSE file or
 * http://www.wtfpl.net/ for more details.
 */

#include <math.h>
#include <time.h>
#include "grid.h"
#include "kv.h"
#include "dynamics.h"
#include "NM_measure.h"
#ifdef _OPENMP
  #include <omp.h>
#endif
 

int main(int argc, char **argv)
{
   if(argc != 2)
   {
      fprintf(stderr, "Usage: ./FDTD input_parameters\n");
      exit(EXIT_FAILURE);
   }
   
   printf("FDTD: solving 1+1D delay PDE\n");
   //printf("https://arxiv.org/abs/1707.05943v1, https://github.com/leofang/FDTD\n");
   printf("This code is released under the WTFPL without any warranty.\n");
   printf("See LICENSE or http://www.wtfpl.net/ for more details.\n");
   printf("Copyright (C) 2018 Leo Fang\n\n");

   #ifdef _OPENMP
     printf("FDTD: the executable is compiled with OpenMP, thus it runs parallelly with %i threads...\n", omp_get_max_threads());
   #else
     printf("FDTD: the executable is compiled without OpenMP, thus it runs serially...\n");
   #endif
   
   printf("FDTD: preparing the grid...\n");
   grid * simulation = initialize_grid(argv[1]);
//   printf("\033[F\033[2KFDTD: preparing the grid...Done!\n");
   printf("FDTD: simulation starts...\n");// fflush(stdout);

   // W defined in grid.h
   W = simulation->w0*I+0.5*simulation->Gamma;

   //grid dimension in x
   int xmin = simulation->nx+1;   //start from x=-Nx*Delta
   int xmax = simulation->Ntotal;

   // initialize time measurement 
   clock_t clock_start, clock_end;
   clock_start = clock();
   #ifdef _OPENMP
     double omp_start = omp_get_wtime();
     int xmax_step = xmax + simulation->nx * (omp_get_max_threads()-1);
   #else
     int xmax_step = xmax;
   #endif

   //simulation starts
   #ifdef _OPENMP
     for(int j=1; j<simulation->Ny; j+=omp_get_max_threads())
   #else
     for(int j=1; j<simulation->Ny; j++) //start from t=1*Delta
   #endif
   {
       for(int i=xmin; i<xmax_step; i++) //start from x=-Nx*Delta
       {
           #ifdef _OPENMP
	   //march each (delayed) thread within range one step in x simultaneously; see paper
           #pragma omp parallel for
	   for(int n=0; n<omp_get_max_threads(); n++)
	   {
              int x_temp = i - n * simulation->nx;
              int t_temp = j + n;
	      if(xmin<=x_temp && x_temp<xmax && t_temp<simulation->Ny)
	         solver(t_temp, x_temp, simulation);
	   }
	   #else
	     solver(j, i, simulation);
           #endif
       }
   }

   // stop the timers
   #ifdef _OPENMP
     double omp_end = omp_get_wtime();
     printf("FDTD: simulation ends, OpenMP time elapsd: %f s\n", omp_end - omp_start);
   #endif
   clock_end = clock();
   double cpu_time_used = ((double) (clock_end - clock_start)) / CLOCKS_PER_SEC;
   printf("FDTD: simulation ends, clock time elapsd: %f s\n", cpu_time_used);


//   printf("FDTD: writing results to files...\n");// fflush(stdout);
//   print_initial_condition(simulation);
//   print_boundary_condition(simulation);
//   printf("******************************************\n");
//   print_psi(simulation);
//   print_grid(simulation);
   if(simulation->save_psi)
   {
      printf("FDTD: saving the wavefunction psi...\n");
      save_psi(simulation, argv[1], creal);
      save_psi(simulation, argv[1], cimag);
      //save_psi(simulation, argv[1], cabs);
   }
   if(simulation->save_psi_square_integral) //for testing init_cond=3
   {
      printf("FDTD: saving the psi^2 integral...\n");
      save_psi_square_integral(simulation, argv[1]);
   }
   if(simulation->save_psi_binary)
   {
      printf("FDTD: saving the wavefunction psi as binary...\n");
      save_psi_binary(simulation, argv[1]);
   }
   if(simulation->save_chi)
   {
      printf("FDTD: saving absolute value of the two-photon wavefunction |chi|...\n");
      save_chi(simulation, argv[1], cabs);
   }
   if(simulation->measure_NM)
   {
      printf("FDTD: calculating lambda and mu for NM measures...\n"); fflush(stdout);
      calculate_NM_measure(simulation, argv[1]);
      save_e0(simulation, argv[1], creal);
      save_e0(simulation, argv[1], cimag);
      save_e1(simulation, argv[1], creal);
      save_e1(simulation, argv[1], cimag);
   }

   free_grid(simulation);

   return EXIT_SUCCESS;
}

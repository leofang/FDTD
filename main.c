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

double getRealTime(); 

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

   //4 possibilities
   #ifdef __FDTD_OPENMP_SUPPORT__
     #ifdef __FDTD_PTHREAD_SUPPORT__
       printf("FDTD: the executable is compiled with both OpenMP (for loop parallelizations) and pthreads (for multi-solver).\n");
     #else
       printf("FDTD: the executable is compiled with OpenMP (for loop paralleiizations and multi-solver).\n");
     #endif
   #else
     #ifdef __FDTD_PTHREAD_SUPPORT__
       printf("FDTD: the executable is compiled with pthreads (for multi-solver).\n");
     #else
       printf("FDTD: the executable is compiled as a serial (single-thread) program.\n");
     #endif
   #endif
   
   printf("FDTD: preparing the grid...\n");
   grid * simulation = initialize_grid(argv[1]);
//   printf("\033[F\033[2KFDTD: preparing the grid...Done!\n");

   // W defined in grid.h
   W = simulation->w0*I+0.5*simulation->Gamma;

   //grid dimension in x and t
   int delay = simulation->nx;
   int xmin = delay+1;   //start from x=-Nx*Delta
   int xmax = simulation->Ntotal;
   int tmax = simulation->Ny;

   printf("FDTD: simulation starts...\n");// fflush(stdout);

   // initialize time measurement 
   clock_t clock_start, clock_end;
   double start, end;
   clock_start = clock();
   start = getRealTime();
   #ifdef _OPENMP
     double omp_start, omp_end;
     omp_start = omp_get_wtime();
   #endif

   //simulation starts
   #ifdef _OPENMP
     #pragma omp parallel
     {
        //region 1: x<=-a
        int kmax = (simulation->minus_a_index-xmin+1)+(tmax-1);
        for(int k=0; k<kmax; k++)
        {
            #pragma omp for
            for(int j=1; j<tmax; j++)
            {
               int i = k-j+xmin;
               if(i>=xmin && i<=simulation->minus_a_index)
               {
                  solver(j, i, simulation);
               }
            }
        }
  
        //region 2: -a<x<=+a
        #pragma omp single
        {
           for(int j=1; j<tmax; j++) //start from t=1*Delta
           {
               for(int i=simulation->minus_a_index+1; i<=simulation->plus_a_index; i++)
               {
                  solver(j, i, simulation);
               }
           }
        }

        //region 3: x>+a
        kmax = (xmax-simulation->plus_a_index-1)+(tmax-1);
        for(int k=0; k<kmax; k++)
        {
            #pragma omp for
            for(int j=1; j<tmax; j++)
            {
                int i = k-j+simulation->plus_a_index+1;
                if(i>=simulation->plus_a_index+1 && i<xmax)
                {
                   solver(j, i, simulation);
                }
            }
        }
     }
   #else
     for(int j=1; j<tmax; j++) //start from t=1*Delta
     {
         for(int i=xmin; i<xmax; i++) //start from x=-Nx*Delta
         {
               solver(j, i, simulation);
         }
     }
   #endif

   #ifdef _OPENMP
     omp_end = omp_get_wtime(); // stop the timers
     printf("FDTD: simulation ends, omp_get_wtime time elapsd: %f s\n", omp_end - omp_start);
   #endif
   end = getRealTime();
   printf("FDTD: simulation ends, getRealTime time elapsd: %f s\n", end - start);
   clock_end = clock();
   printf("FDTD: simulation ends, clock time elapsd: %f s\n", ((double) (clock_end - clock_start)) / CLOCKS_PER_SEC);

//   printf("FDTD: writing results to files...\n");// fflush(stdout);
//   print_initial_condition(simulation);
//   print_boundary_condition(simulation);
//   printf("******************************************\n");
//   print_psi(simulation);
//   print_grid(simulation);

   if(simulation->save_psi)
   {
      printf("FDTD: saving the wavefunction psi...\n");
      start = getRealTime();
      save_psi(simulation, argv[1], creal);
      save_psi(simulation, argv[1], cimag);
      //save_psi(simulation, argv[1], cabs);
      end = getRealTime();
      printf("FDTD: psi saved, getRealTime time elapsd: %f s\n", end - start);
   }
   if(simulation->save_psi_square_integral) //for testing init_cond=3
   {
      printf("FDTD: saving the psi^2 integral...\n");
      start = getRealTime();
      save_psi_square_integral(simulation, argv[1]);
      end = getRealTime();
      printf("FDTD: psi^2 integral saved, getRealTime time elapsd: %f s\n", end - start);
   }
   if(simulation->save_psi_binary)
   {
      printf("FDTD: saving the wavefunction psi as binary...\n");
      start = getRealTime();
      save_psi_binary(simulation, argv[1]);
      end = getRealTime();
      printf("FDTD: psi binary saved, getRealTime time elapsd: %f s\n", end - start);
   }
   if(simulation->save_chi)
   {
      printf("FDTD: saving absolute value of the two-photon wavefunction |chi|...\n");
      start = getRealTime();
      save_chi(simulation, argv[1], cabs);
      end = getRealTime();
      printf("FDTD: chi saved, getRealTime time elapsd: %f s\n", end - start);
   }
   if(simulation->measure_NM)
   {
      printf("FDTD: calculating lambda and mu for NM measures...\n"); fflush(stdout);
      start = getRealTime();
      calculate_NM_measure(simulation, argv[1]);
      end = getRealTime();
      printf("FDTD: lambda and mu saved, getRealTime time elapsd: %f s\n", end - start);
      start = getRealTime();
      save_e0(simulation, argv[1], creal);
      save_e0(simulation, argv[1], cimag);
      end = getRealTime();
      printf("FDTD: e0 saved, getRealTime time elapsd: %f s\n", end - start);
      start = getRealTime();
      save_e1(simulation, argv[1], creal);
      save_e1(simulation, argv[1], cimag);
      end = getRealTime();
      printf("FDTD: e1 saved, getRealTime time elapsd: %f s\n", end - start);
   }

   free_grid(simulation);

   return EXIT_SUCCESS;
}

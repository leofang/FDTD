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
#include "pthread_solver.h"

double getRealTime(); //defined in getRealTime.c
double complex W;     //W=I*w0+0.5*Gamma, to be accessed by multiple threads


int main(int argc, char **argv)
{
   if(argc != 2)
   {
      fprintf(stderr, "Usage: ./FDTD input_parameters\n");
      exit(EXIT_FAILURE);
   }
   
   printf("FDTD: solving 1+1D delay PDE\n");
   printf("This code is released under the MIT license, see LICENSE for more details.\n");
   printf("Copyright (C) 2018 Leo Fang\n\n");
   //printf("https://arxiv.org/abs/1707.05943v1, https://github.com/leofang/FDTD\n");
   //printf("This code is released under the WTFPL without any warranty.\n");
   //printf("See LICENSE or http://www.wtfpl.net/ for more details.\n");

   //4 possibilities
   #ifdef __FDTD_OPENMP_SUPPORT__
     #ifdef __FDTD_PTHREAD_SUPPORT__
       printf("FDTD: the executable is compiled with both OpenMP (for loop parallelizations) and pthreads (for multi-solver).\n");
     #else
       printf("FDTD: the executable is compiled with OpenMP (for loop parallelizations and multi-solver).\n");
     #endif
   #else
     #ifdef __FDTD_PTHREAD_SUPPORT__
       printf("FDTD: the executable is compiled with pthreads (for multi-solver).\n");
     #else
       printf("FDTD: the executable is compiled as a serial (single-thread) program.\n");
     #endif
   #endif
   
   printf("FDTD: preparing the grid...\n");
   //TODO: execute initialize_grid in another thread and do OpenMP there,
   //so when it's done the OpenMP thread pool can either be eliminated or handed out
   grid * simulation = initialize_grid(argv[1]);

   #ifdef __FDTD_PTHREAD_SUPPORT__ //for pthread solvers
     int Nth = simulation->Nth;
     pthread_t thread_id[Nth];
     solver_info thread_id_list[Nth];

     //initiallize arrays to record the current positions of the solvers
     //which will be locked 
     int solver_x_positions[Nth];
     int solver_t_positions[Nth];
     pthread_cond_t solver_halt[Nth];
     pthread_mutex_t solver_locks[Nth];

     printf("FDTD: %i threads will be used.\n", Nth); //TODO: re-organize this

     //initialize pthreads
     for(int i=0; i < Nth; i++)
     {
	 //initialize arrays
	 solver_x_positions[i] = simulation->nx+1;
	 solver_t_positions[i] = 1+i;
         if(pthread_cond_init(&solver_halt[i], NULL))
	 {
	    fprintf(stderr, "FDTD: pthread_cond_init fails for %i-th thread, abort.\n", i);
	    exit(EXIT_FAILURE);
	 }
         if(pthread_mutex_init(&solver_locks[i], NULL))
	 {
	    fprintf(stderr, "FDTD: pthread_mutex_init fails for %i-th thread, abort.\n", i);
	    exit(EXIT_FAILURE);
	 }

	 //initialize thread arguments
         thread_id_list[i].id  = i;
         thread_id_list[i].Nth = Nth;
         thread_id_list[i].solver_x_positions = solver_x_positions;
         thread_id_list[i].solver_t_positions = solver_t_positions;
         thread_id_list[i].solver_halt = solver_halt;
         thread_id_list[i].solver_locks = solver_locks;
         thread_id_list[i].simulation = simulation;
     }
   #else //for either openmp or single-thread solver
     //grid dimension in x and t
     int xmin = simulation->nx+1;   //start from x=-Nx*Delta
     int xmax = simulation->Ntotal;
     int tmax = simulation->Ny;
   #endif

   printf("FDTD: simulation starts...\n");// fflush(stdout);

   // initialize time measurement
   getRealTime(); //set up the internal struct
   clock_t clock_start, clock_end;
   double start, end;
   clock_start = clock();
   start = getRealTime();
   #ifdef _OPENMP
     double omp_start, omp_end;
     omp_start = omp_get_wtime();
   #endif

   //simulation starts
   #ifdef __FDTD_PTHREAD_SUPPORT__ //using pthread solvers
      for(int i=0; i < Nth; i++)
      {
          if(pthread_create( &thread_id[i], NULL, solver_wrapper, (void *)&thread_id_list[i]))
	  {
	     fprintf(stderr, "FDTD: pthread_create fails for %i-th thread, abort.\n", i);
	     exit(EXIT_FAILURE);
	  }
      }
      for(int i=0; i < Nth; i++)
      {
          if(pthread_join( thread_id[i], NULL))
	  {
	     fprintf(stderr, "FDTD: pthread_join fails for %i-th thread,\n      attempting to proceed...\n", i);
	  }
      }
   #elif defined _OPENMP //using wavefront approach
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
                  solver(j, i, simulation);
            }
        }
  
        //region 2: -a<x<=+a
        #pragma omp single
        {
           for(int j=1; j<tmax; j++) //start from t=1*Delta
           {
               for(int i=simulation->minus_a_index+1; i<=simulation->plus_a_index; i++)
                  solver(j, i, simulation);
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
                   solver(j, i, simulation);
            }
        }
     }
   #else //single-thread version
     for(int j=1; j<tmax; j++) //start from t=1*Delta
     {
         for(int i=xmin; i<xmax; i++) //start from x=-Nx*Delta
            solver(j, i, simulation);
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

   #ifdef __FDTD_PTHREAD_SUPPORT__ 
     //clean up
     for(int i=0; i < Nth; i++)
     {
	if(pthread_cond_destroy(&solver_halt[i]))
	{
	   fprintf(stderr, "FDTD: pthread_cond_destroy fails for %i-th thread,\n      attempting to proceed...\n", i);
	}
	if(pthread_mutex_destroy(&solver_locks[i]))
	{
	   fprintf(stderr, "FDTD: pthread_mutex_destroy fails for %i-th thread,\n      attempting to proceed...\n", i);
	}
     }
   #endif

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

      //TODO: add a seperate flag for save_BIC
      start = getRealTime();
      save_BIC(simulation, argv[1]);
      end = getRealTime();
      printf("FDTD: BIC saved, getRealTime time elapsd: %f s\n", end - start);
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
   if(simulation->save_chi_map)
   {
      printf("FDTD: saving absolute value of the two-photon wavefunction |chi| as a 2D map...\n");
      start = getRealTime();
      save_chi_map(simulation, argv[1], cabs);
      end = getRealTime();
      printf("FDTD: chi map saved, getRealTime time elapsd: %f s\n", end - start);
      start = getRealTime();
      check_normalization(simulation);
      end = getRealTime();
      printf("FDTD: normalization factor checked, getRealTime time elapsd: %f s\n", end - start);
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

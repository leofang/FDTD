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

   #ifdef __FDTD_PTHREAD_SUPPORT__
      printf("FDTD: the executable is compiled with pthreads, so there could be multiple solvers (depending on Nth).\n");
      #ifdef __FDTD_OPENMP_SUPPORT__
         printf("FDTD: the executable is compiled with OpenMP.\n");
      #endif
   #else
      printf("FDTD: the executable is single-threaded.\n");
   #endif

   printf("FDTD: preparing the grid...\n");
   //TODO: execute initialize_grid in another thread and do OpenMP there,
   //so when it's done the OpenMP thread pool is eliminated
   grid * simulation = initialize_grid(argv[1]);

   #ifdef __FDTD_PTHREAD_SUPPORT__
     int Nth = simulation->Nth;
     pthread_t thread_id[Nth];
     solver_info thread_id_list[Nth];

     //initiallize arrays to record the current positions of the solvers
     //which will be locked 
     int solver_x_positions[Nth];
     int solver_t_positions[Nth];
     pthread_cond_t solver_halt[Nth];
     pthread_mutex_t solver_locks[Nth];

     printf("FDTD: %i threads will be used.\n", Nth);

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
   #endif
   
   printf("FDTD: simulation starts...\n");// fflush(stdout);

   // initialize time measurement
   getRealTime(); //set up the internal struct
   clock_t clock_start, clock_end;
   clock_start = clock();
   double start = getRealTime();

   //simulation starts
   #ifdef __FDTD_PTHREAD_SUPPORT__
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
   #else //single thread
      W = simulation->w0*I+0.5*simulation->Gamma; //W defined in grid.h
      for(int j=1; j<simulation->Ny; j++) //start from t=1*Delta
      {
	  for(int i=simulation->nx+1; i<simulation->Ntotal; i++) //start from x=-Nx*Delta
	  {
	     solver(j, i, simulation);
	  }
      }
   #endif

   double end = getRealTime();
   clock_end = clock();
   double cpu_time_used = ((double) (clock_end - clock_start)) / CLOCKS_PER_SEC;
   printf("FDTD: simulation ends, solver spent: %f s (getRealTime)\n", end - start);
   printf("FDTD: simulation ends, solver spent: %f s (clock)\n", cpu_time_used);
     
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

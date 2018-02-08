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
//#include <time.h>   //for clock_gettime
#include "grid.h"
#include "kv.h"
#include "dynamics.h"
#include "NM_measure.h"
#include "pthread_solver.h"

//#ifndef NTHREADS
//#define NTHREADS 4
//#endif
//  double timing[2]={0};
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

   #ifdef _POSIX_THREADS
      printf("FDTD: the executable is compiled with pthreads, so it will be multi-threaded.\n");
   #else
      printf("FDTD: the executable is single-threaded.\n");
   #endif
   printf("FDTD: preparing the grid...\n");
   grid * simulation = initialize_grid(argv[1]);
//   printf("\033[F\033[2KFDTD: preparing the grid...Done!\n");

   #ifdef _POSIX_THREADS
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
         pthread_cond_init(&solver_halt[i], NULL);
         pthread_mutex_init(&solver_locks[i], NULL);

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
   clock_t clock_start, clock_end;
   clock_start = clock();
   #ifdef _POSIX_THREADS
     double start = getRealTime(); //for timing
   #endif

   //simulation starts
   #ifdef _POSIX_THREADS
      for(int i=0; i < Nth; i++)
          pthread_create( &thread_id[i], NULL, solver_wrapper, (void *)&thread_id_list[i]);
      for(int i=0; i < Nth; i++)
          pthread_join( thread_id[i], NULL);
   #else
      //TODO
      //solver_wrapper(thread_id_list[0]);
      //for(int i=)
   #endif

   #ifdef _POSIX_THREADS
     double end = getRealTime();
     printf("FDTD: simulation ends, getRealTime measures: %lf s\n", end - start);
     
     //clean up
     for(int i=0; i < Nth; i++)
     {
	pthread_cond_destroy(&solver_halt[i]);
	pthread_mutex_destroy(&solver_locks[i]);
     }
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

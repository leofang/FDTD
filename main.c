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
  #include <unistd.h> //for execv
#endif
//  double timing[2]={0};
 

int main(int argc, char **argv)
{
   // #ifdef _OPENMP
   //   if(!omp_get_cancellation())
   //   {
   //      printf("Cancellations were not enabled, enabling cancellation and rerunning program\n\n");
   //      putenv("OMP_CANCELLATION=true");
   //      execv(argv[0], argv);
   //   }
   // #endif

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
     int Nth = omp_get_max_threads(); //get available number of threads
     printf("FDTD: the executable is compiled with OpenMP, so it runs parallelly with %i threads...\n", Nth);
     //initiallize an array to record the current x-positions of the solvers
     //the array will be locked 
     int * solver_x_positions = malloc(Nth * sizeof(*solver_x_positions));
     if(!solver_x_positions)
     {
        fprintf(stderr, "%s: malloc fails, abort.\n", __func__);
	exit(EXIT_FAILURE);
     }
     omp_lock_t position_locks[Nth];
     for(int i=0; i<Nth; i++)
     {
        omp_init_lock(&position_locks[i]);
     }
   #else
     printf("FDTD: the executable is compiled without OpenMP, so it runs serially...\n");
   #endif
   
   printf("FDTD: preparing the grid...\n");
   grid * simulation = initialize_grid(argv[1]);
//   printf("\033[F\033[2KFDTD: preparing the grid...Done!\n");

   // W defined in grid.h
   W = simulation->w0*I+0.5*simulation->Gamma;

   //grid dimension in x and t
   int xmin = simulation->nx+1;   //start from x=-Nx*Delta
   int xmax = simulation->Ntotal;
   int tmax = simulation->Ny;       

   printf("FDTD: simulation starts...\n");// fflush(stdout);

   // initialize time measurement 
   clock_t clock_start, clock_end;
   clock_start = clock();
   #ifdef _OPENMP
     double omp_start = omp_get_wtime();
     int dummy = 0; //avoid compiler loop optimization
   #endif

   //simulation starts
   #ifdef _OPENMP
     #pragma omp parallel firstprivate(dummy)
     {
     int id = omp_get_thread_num();
     for(int j=1+id; j<tmax+id; j+=Nth)
     {
	 //solver_x_positions[id] = xmin-id*simulation->nx; //reset
	 //for(int i=xmin-id*simulation->nx; i<xmax+(Nth-id-1)*simulation->nx; i++) //start from x=-Nx*Delta
	 solver_x_positions[id] = xmin; //reset
	 for(int i=xmin; i<xmax; i++) //start from x=-Nx*Delta
	 {  
   #else
     for(int j=1; j<tmax; j++) //start from t=1*Delta
     {
         for(int i=xmin; i<xmax; i++) //start from x=-Nx*Delta
         {
   #endif
             #ifdef _OPENMP
               #pragma omp cancel parallel if(j==tmax-1 && i==xmax)  //stop if reach the grid boundary
               //march each (delayed) thread within range one step in x simultaneously; see paper
               if(xmin<=i && i<xmax && j<tmax)
	       {
		  if(id!=0) //let the first guy runs unrestrainedly
		  {
		     //others must be dealyed by a certain amount, otherwise busy waiting
		     while(solver_x_positions[id-1]-i<simulation->nx && solver_x_positions[id-1]<xmax)
		     {
			dummy++; //if do nothing in the loop, this loop will be mistakenly optimized by the compiler
		     }
		  }
                  solver(j, i, simulation);
	       }
	       solver_x_positions[id]++;
               #pragma omp flush
             #else
               solver(j, i, simulation);
             #endif
         }
         #pragma omp barrier
     }
   #ifdef _OPENMP
     }
   #endif

   #ifdef _OPENMP
     double omp_end = omp_get_wtime(); // stop the timers
     printf("FDTD: simulation ends, OpenMP time elapsd: %f s\n", omp_end - omp_start);
     for(int i=0; i<Nth; i++)
     {
        omp_destroy_lock(&position_locks[i]);
     }
     free(solver_x_positions);
   #endif
   clock_end = clock();
   double cpu_time_used = ((double) (clock_end - clock_start)) / CLOCKS_PER_SEC;
   printf("FDTD: simulation ends, clock time elapsd: %f s\n", cpu_time_used);
   //  printf("FDTD: (part 1 take in total: %f s)\n", timing[0]);
   //  printf("FDTD: (part 2 take in total: %f s)\n", timing[1]);


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

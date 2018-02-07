/*
 * Copyright (C) 2016 Leo Fang <leofang@phy.duke.edu>
 *
 * This program is free software. It comes without any warranty,
 * to the extent permitted by applicable law. You can redistribute
 * it and/or modify it under the terms of the WTFPL, Version 2, as
 * published by Sam Hocevar. See the accompanying LICENSE file or
 * http://www.wtfpl.net/ for more details.
 */

#ifndef __FDTD_PTHREAD_H__
#define __FDTD_PTHREAD_H__

#include <pthread.h>
#include "dynamics.h"


//TODO: code re-factorize
struct _FDTD_barrier
{
   int Nth;
   int barrier_counter;
   pthread_cond_t barrier_cond;
   pthread_mutex_t barrier_counter_lock;
};
typedef struct _FDTD_barrier FDTD_barrier;


//solver_locks[i] protects solver_x_positions[i]
struct _solver_info
{
   int id;
   int Nth;
   int * barrier_counter;
   pthread_cond_t * barrier_cond;
   pthread_mutex_t * barrier_counter_lock;
   int * solver_x_positions;
   pthread_cond_t * solver_halt;
   pthread_mutex_t * solver_locks;
   grid * simulation;  
};
typedef struct _solver_info solver_info;


//turn on Nth(=NTHREADS by default) threads to run the solver
//the argument arg is of type solver_info* 
inline void * solver_wrapper(void * arg)
{
   //copy the solver info 
   int id  = ((solver_info *)arg)->id; 
   int Nth = ((solver_info *)arg)->Nth; 
   int * barrier_counter = ((solver_info *)arg)->barrier_counter;
   pthread_cond_t * barrier_cond = ((solver_info *)arg)->barrier_cond;
   pthread_mutex_t * barrier_counter_lock = ((solver_info *)arg)->barrier_counter_lock;
   int * solver_x_positions = ((solver_info *)arg)->solver_x_positions; 
   pthread_cond_t * solver_halt = ((solver_info *)arg)->solver_halt;
   pthread_mutex_t * solver_locks = ((solver_info *)arg)->solver_locks;

   //copy the grid info
   grid * simulation = ((solver_info *)arg)->simulation;
   int xmin = simulation->nx+1;   //start from x=-Nx*Delta
   int xmax = simulation->Ntotal;
   int tmin = 1;
   int tmax = simulation->Ny; 
   int nx = simulation->nx;
   W = simulation->w0*I+0.5*simulation->Gamma;  // W defined in grid.h

   //a temporary variable for the previous solver's position
   int temp;

   //initialize solver position
   pthread_mutex_lock(&solver_locks[id]);
   solver_x_positions[id]=xmin-id*nx;
   pthread_mutex_unlock(&solver_locks[id]);

   for(int j=tmin + id; j<tmax+id; j+=Nth)
   {
       //do not waste time on busy-looping; EXPERIMENTAL!!!
       if(j>=tmax)
          break;

       for(int i=xmin-id*nx; i<xmax+(Nth-id-1)*nx; i++) //start from x=-Nx*Delta
       {
          //march each (delayed) thread within range one step in x simultaneously; see paper
          if(xmin<=i && i<xmax && j<tmax)
          {
             if(id!=0)
             {
		pthread_mutex_lock(&solver_locks[id-1]);
   	        temp = solver_x_positions[id-1];
                while(temp-i<nx && temp<xmax)
                {
		   pthread_cond_wait(&solver_halt[id], &solver_locks[id-1]);
                   temp = solver_x_positions[id-1];
                }
		pthread_mutex_unlock(&solver_locks[id-1]);
             }
	     solver(j, i, simulation); //TODO: use mutex to enforce memory barrier??
          }

	  //increment the current position and wake up the next solver if it is waiting
          pthread_mutex_lock(&solver_locks[id]);
          solver_x_positions[id]++;
          pthread_mutex_unlock(&solver_locks[id]);
	  if(id<Nth-1)
	     pthread_cond_signal(&solver_halt[id+1]);
       }

       //barrier begins
       pthread_mutex_lock(barrier_counter_lock);
       (*barrier_counter)++;
       if((*barrier_counter)==Nth) //the last guy resets all solver positions and then wakes up everybody
       {
          //reset solver positions
          for(int th=0; th<Nth; th++)
	  {
             pthread_mutex_lock(&solver_locks[th]);
             solver_x_positions[th]=xmin-th*nx;
             pthread_mutex_unlock(&solver_locks[th]);
	  }

	  //reset the counter 
	  *barrier_counter=0;

          pthread_cond_broadcast(barrier_cond);
       }
       else //others go to sleep
       {
          pthread_cond_wait(barrier_cond, barrier_counter_lock);
       }
       pthread_mutex_unlock(barrier_counter_lock);
       //barrier ends
   }

   //the last barrier in case some threads end earlier while others are still working
   pthread_mutex_lock(barrier_counter_lock);
   (*barrier_counter)++;
   if((*barrier_counter)==Nth)
      pthread_cond_broadcast(barrier_cond);
   pthread_mutex_unlock(barrier_counter_lock);

   return NULL;
}

//#endif //_POSIX_THREADS
#endif //__FDTD_PTHREAD_H__

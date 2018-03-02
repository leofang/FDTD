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
#ifdef __FDTD_PTHREAD_SUPPORT__
  //default: turn OFF pthreads unless __FDTD_PTHREAD_SUPPORT__ is explicitly given during compilation
#include <unistd.h> //for _POSIX_THREADS
#define _MULTI_THREADED
#include <pthread.h>
#include "dynamics.h"


//solver_locks[i] protects solver_x_positions[i] and solver_t_positions[i]
//solver_halt[i] is used to wait for the previous solver (i-1 if i!=0 and i=Nth-1 if i==0)
struct _solver_info
{
   int id;
   int Nth;
   int * solver_x_positions;
   int * solver_t_positions;
   pthread_cond_t * solver_halt;
   pthread_mutex_t * solver_locks;
   grid * simulation;  
};
typedef struct _solver_info solver_info;


//turn on Nth threads to run the solver
//the argument arg is of type solver_info* 
inline void * solver_wrapper(void * arg)
{
   //copy the solver info 
   int id  = ((solver_info *)arg)->id; 
   int Nth = ((solver_info *)arg)->Nth; 
   int * solver_x_positions = ((solver_info *)arg)->solver_x_positions; 
   int * solver_t_positions = ((solver_info *)arg)->solver_t_positions; 
   pthread_cond_t * solver_halt = ((solver_info *)arg)->solver_halt;
   pthread_mutex_t * solver_locks = ((solver_info *)arg)->solver_locks;
   grid * simulation = ((solver_info *)arg)->simulation;

   //copy the grid info
   int xmin = simulation->nx+1;   //start from x=-Nx*Delta
   int xmax = simulation->Ntotal;
   int tmin = 1;
   int tmax = simulation->Ny; 
   int nx = simulation->nx;
   W = simulation->w0*I+0.5*simulation->Gamma;  // W defined in grid.h

   int previous_x, previous_t;               //previous solver's positions
   int previous_id = (id==0 ? Nth-1 : id-1); //previous solver's id
   int next_id = (id==Nth-1 ? 0 : id+1);     //next solver's id

   for(int j=tmin + id; j<tmax+id; j+=Nth)
   {
       //do not waste time on busy-looping; EXPERIMENTAL!!!
       if(j>=tmax)
          break;

       for(int i=xmin; i<xmax; i++) //start from x=-Nx*Delta
       {
          //march each (delayed) thread within range one step in x simultaneously; see paper
          if(j>tmin)
          {
	     pthread_mutex_lock(&solver_locks[previous_id]);
   	     previous_x = solver_x_positions[previous_id];
   	     previous_t = solver_t_positions[previous_id];
             while(j>previous_t && previous_x-i<nx && previous_x<xmax)
             {
	        pthread_cond_wait(&solver_halt[id], &solver_locks[previous_id]);
                previous_x = solver_x_positions[previous_id];
                previous_t = solver_t_positions[previous_id];
             }
	     pthread_mutex_unlock(&solver_locks[previous_id]);
          }
	  solver(j, i, simulation); 

	  //update the current position and wake up the next solver if it is waiting
          pthread_mutex_lock(&solver_locks[id]);
	  if(i==xmin)
	  {
	     solver_x_positions[id]=i;
	     solver_t_positions[id]=j;
	  }
	  else
             solver_x_positions[id]++;
          pthread_mutex_unlock(&solver_locks[id]);
	  pthread_cond_signal(&solver_halt[next_id]);
       }
   }

   //place the solver outside the grid to indicate that it is finished
   pthread_mutex_lock(&solver_locks[id]);
   solver_t_positions[id]=tmax+id;
   pthread_mutex_unlock(&solver_locks[id]);
   pthread_cond_signal(&solver_halt[next_id]);

   return NULL;
}

#endif //__FDTD_PTHREAD_SUPPORT__
#endif //__FDTD_PTHREAD_H__

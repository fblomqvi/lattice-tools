/* an_solver.c
   Copyright (C) 2017 Ferdinand Blomqvist

   This program is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 as published by
   the Free Software Foundation. 

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along with
   this program. If not, see <http://www.gnu.org/licenses/>.  
   
   Written by Ferdinand Blomqvist. */

#include "dbg.h"
#include "an_solver.h"
#include <math.h>

struct s_AN_workspace
{
    size_t* idxs;
    double* delta;
    size_t size;
};

static int cmp(const void* left, const void* right, void* arg);

AN_WS* AN_WS_alloc(size_t size)
{
    AN_WS* ws = malloc(sizeof(AN_WS));
    libcheck_mem(ws);

    ws->idxs = malloc(size * sizeof(size_t));
    llibcheck_mem(ws, error_a);

    ws->delta = malloc(size * sizeof(double));
    llibcheck_mem(ws, error_b);

    ws->size = size;
    return ws;

error_b:
    free(ws->idxs);
error_a:
    free(ws);
error:
    return NULL;
}

void AN_WS_free(AN_WS* ws)
{
    if(ws)
    {
        free(ws->delta);
        free(ws->idxs);
        free(ws);
    }
}

void AN_solve(long* clp, double* recv, AN_WS* ws)
{
    // Project recv onto the hyperplane of dimension size-1 where A_n lives.
    double sum = recv[0];
    for(size_t i = 1; i < ws->size; i++)
        sum += recv[i];

    double s = sum / ws->size;
    debug("s: %f", s);
    for(size_t i = 0; i < ws->size; i++)
    {
        recv[i] -= s;
        debug("recv[%zu]: %f", i, recv[i]);
    }

   AN_solve_np(clp, recv, ws); 
}

void AN_solve_np(long* clp, double* recv, AN_WS* ws)
{
    long Delta = clp[0] = round(recv[0]);
    ws->delta[0] = recv[0] - clp[0];
    debug("round(recv[0]): %ld;\t delta[0]: %f", clp[0], ws->delta[0]);
    for(size_t i = 1; i < ws->size; i++)
    {
        Delta += clp[i] = round(recv[i]);
        ws->delta[i] = recv[i] - clp[i];
        debug("round(recv[%zu]): %ld;\t delta[%zu]: %f", i, clp[i], i, ws->delta[i]);
    }

    debug("Delta: %ld", Delta);
    if(Delta == 0)
        return;
        
    // Initialize the array that will be sorted 
    for(size_t i = 0; i < ws->size; i++)
        ws->idxs[i] = i;
    
    // Sort the array
    qsort_r(ws->idxs, ws->size, sizeof(size_t), cmp, ws->delta);

    if(Delta > 0)
    {
        for(size_t i = 0; i < (size_t) Delta; i++)
            clp[ws->idxs[i]]--;
    }
    else
    {
        for(size_t i = ws->size - 1; i >= ws->size - 1 - Delta; i--)
            clp[ws->idxs[i]]++;
    }
}

static int cmp(const void* left, const void* right, void* arg)
{
    size_t l = *((size_t*) left);
    size_t r = *((size_t*) right);
    double* delta = (double*) arg;
    double diff = delta[l] - delta[r];
    if(diff < 0)
        return -1;
    else if(diff > 0)
        return 1;
    else
        return 0;
}

/* rnd_point.h
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

#ifndef FB_LIBLATTICE_RND_POINT_H
#define FB_LIBLATTICE_RND_POINT_H

#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct s_rnd_point_conf
{
    long* hyperplane;
    size_t dimension;
    size_t num_cwords;
    int A_n;
} RND_PNT_CONF;

RND_PNT_CONF* RND_PNT_CONF_alloc(void);
void RND_PNT_CONF_free(RND_PNT_CONF* conf);

RND_PNT_CONF* RND_PNT_CONF_read(FILE* file);
int RND_PNT_CONF_write(FILE* file, RND_PNT_CONF* conf);

void rnd_point_get(double* p, size_t size, gsl_rng* rng);
void rnd_point_get_A_N(double* p, size_t size, gsl_rng* rng);
void rnd_point_get_SPC(double* p, size_t size, gsl_rng* rng, long* coeffs);

void rnd_point_uniform(double* p, size_t size, gsl_rng* rng, unsigned long n);
void rnd_point_uniform_A_N(double* p, size_t size, gsl_rng* rng, unsigned long n);
void rnd_point_uniform_SPC(double* p, size_t size, gsl_rng* rng, 
                            long* coeffs, unsigned long n);

void rnd_point_min_max(double* p, size_t size, gsl_rng* rng, long min, long max);
void rnd_point_min_max_A_N(double* p, size_t size, gsl_rng* rng, long min, long max);
void rnd_point_min_max_SPC(double* p, size_t size, gsl_rng* rng, 
                            long* coeffs, long min, long max);

void rnd_point_min_max_fast(double* p, size_t size, gsl_rng* rng, long min, long max);
void rnd_point_min_max_A_N_fast(double* p, size_t size, gsl_rng* rng, long min, long max);
void rnd_point_min_max_SPC_fast(double* p, size_t size, gsl_rng* rng, 
                                long* coeffs, long min, long max);

#endif /* FB_LIBLATTICE_RND_POINT_H */

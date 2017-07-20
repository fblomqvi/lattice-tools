/* rnd_point.c
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

#include "rnd_point.h"
#include "dbg.h"
#include <assert.h>

RND_PNT_CONF* RND_PNT_CONF_alloc(void)
{
    RND_PNT_CONF* conf = malloc(sizeof(RND_PNT_CONF));
    libcheck_mem(conf);
    conf->hyperplane = NULL;
    return conf;

error:
    return NULL;
}

void RND_PNT_CONF_free(RND_PNT_CONF* conf)
{
    if(conf)
    {
        free(conf->hyperplane);
        free(conf);
    }
}

RND_PNT_CONF* RND_PNT_CONF_read(FILE* file)
{
    RND_PNT_CONF* conf = malloc(sizeof(RND_PNT_CONF));
    libcheck_mem(conf);

    libcheck(fread(conf, sizeof(RND_PNT_CONF), 1, file) == 1, "fread failed");
    if(conf->hyperplane)
    {
        conf->hyperplane = malloc(conf->dimension * sizeof(long));
        llibcheck_mem(conf->hyperplane, error_a);
        llibcheck(fread(conf->hyperplane, sizeof(long), conf->dimension, 
                    file) == conf->dimension, error_b, "fread failed");
    }

    return conf;

error_b:
    free(conf->hyperplane);
error_a:
    free(conf);
error:
    return NULL;
}

int RND_PNT_CONF_write(FILE* file, RND_PNT_CONF* conf)
{
    libcheck(fwrite(conf, sizeof(RND_PNT_CONF), 1, file) == 1, "fwrite failed");
    if(conf->hyperplane)
    {
        libcheck(fwrite(conf->hyperplane, sizeof(long), conf->dimension, 
                    file) == conf->dimension, "fwrite failed");
    }

    return 0;
error:
    return -1;
}

void rnd_point_get(double* p, size_t size, gsl_rng* rng)
{
    for(size_t i = 0; i < size; i++)
        p[i] = gsl_rng_get(rng) + gsl_rng_uniform(rng);
}

void rnd_point_get_A_N(double* p, size_t size, gsl_rng* rng)
{
    double sum = 0.0;
    for(size_t i = 0; i < size - 1; i++)
    {
        p[i] = gsl_rng_get(rng) + gsl_rng_uniform(rng);
        sum += p[i];
    }
    
    p[size-1] = -sum;
}

void rnd_point_get_SPC(double* p, size_t size, gsl_rng* rng, long* coeffs)
{
    double sum = 0.0;
    for(size_t i = 0; i < size - 1; i++)
    {
        p[i] = gsl_rng_get(rng) + gsl_rng_uniform(rng);
        sum += coeffs[i] * p[i];
    }
    
    p[size-1] = -sum / coeffs[size-1];
}

void rnd_point_uniform(double* p, size_t size, gsl_rng* rng, unsigned long n)
{
    for(size_t i = 0; i < size; i++)
        p[i] = gsl_rng_uniform_int(rng, n) + gsl_rng_uniform(rng);
}

void rnd_point_uniform_A_N(double* p, size_t size, gsl_rng* rng, unsigned long n)
{
    double sum = 0.0;
    for(size_t i = 0; i < size - 1; i++)
    {
        p[i] = gsl_rng_uniform_int(rng, n) + gsl_rng_uniform(rng);
        sum += p[i];
    }
    
    p[size-1] = -sum;
}

void rnd_point_uniform_SPC(double* p, size_t size, gsl_rng* rng, 
                            long* coeffs, unsigned long n)
{
    double sum = 0.0;
    for(size_t i = 0; i < size - 1; i++)
    {
        p[i] = gsl_rng_uniform_int(rng, n) + gsl_rng_uniform(rng);
        sum += coeffs[i] * p[i];
    }
    
    p[size-1] = -sum / coeffs[size-1];
}

void rnd_point_min_max(double* p, size_t size, gsl_rng* rng, long min, long max)
{
    assert(max > min);
    unsigned long n = max - min;
    for(size_t i = 0; i < size; i++)
        p[i] = min + (long) gsl_rng_uniform_int(rng, n) + gsl_rng_uniform(rng);
}

void rnd_point_min_max_A_N(double* p, size_t size, gsl_rng* rng, long min, long max)
{
    assert(max > min);
    unsigned long n = max - min;
    double sum = 0.0;
    for(size_t i = 0; i < size - 1; i++)
    {
        p[i] = min + (long) gsl_rng_uniform_int(rng, n) + gsl_rng_uniform(rng);
        sum += p[i];
    }
    
    p[size-1] = -sum;
}

void rnd_point_min_max_SPC(double* p, size_t size, gsl_rng* rng, 
                            long* coeffs, long min, long max)
{
    assert(max > min);
    unsigned long n = max - min;
    double sum = 0.0;
    for(size_t i = 0; i < size - 1; i++)
    {
        p[i] = min + (long) gsl_rng_uniform_int(rng, n) + gsl_rng_uniform(rng);
        sum += coeffs[i] * p[i];
    }
    
    p[size-1] = -sum / coeffs[size-1];
}

void rnd_point_min_max_fast(double* p, size_t size, gsl_rng* rng, long min, long max)
{
    assert(max > min);
    unsigned long n = max - min;
    for(size_t i = 0; i < size; i++)
        p[i] = min + gsl_rng_uniform(rng) * n;
}

void rnd_point_min_max_A_N_fast(double* p, size_t size, gsl_rng* rng, long min, long max)
{
    assert(max > min);
    unsigned long n = max - min;
    double sum = 0.0;
    for(size_t i = 0; i < size - 1; i++)
    {
        p[i] = min + gsl_rng_uniform(rng) * n;
        sum += p[i];
    }
    
    p[size-1] = -sum;
}

void rnd_point_min_max_SPC_fast(double* p, size_t size, gsl_rng* rng, 
                            long* coeffs, long min, long max)
{
    assert(max > min);
    unsigned long n = max - min;
    double sum = 0.0;
    for(size_t i = 0; i < size - 1; i++)
    {
        p[i] = min + gsl_rng_uniform(rng) * n;
        sum += coeffs[i] * p[i];
    }
    
    p[size-1] = -sum / coeffs[size-1];
}

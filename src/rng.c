/* rng.c
   Copyright (C) 2016 Ferdinand Blomqvist

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
#include "rng.h"
#include <string.h>
#include <time.h>

unsigned long get_random_seed(void)
{
    struct timespec tp;
    int rc = clock_gettime(CLOCK_REALTIME, &tp);
    libcheck(rc == 0, "clock_gettime failed");
    srand(tp.tv_sec + tp.tv_nsec);
    return rand();

error:
    // Lower resolution option if CLOCK_REALTIME is not available
    srand(time(NULL));
    return rand();
}

const gsl_rng_type* get_rng_type(const char* rng_name, const gsl_rng_type** rng_types)
{
    while(*rng_types)
    {
        if(!strcmp(rng_name, (*rng_types)->name))
            return *rng_types;

        rng_types++;
    }

    return NULL;
}

int print_rngs(FILE* file, const gsl_rng_type** rng_types)
{
    int rc = fprintf(file, "Available random number generators are:\n");
    libcheck(rc >= 0, "fprintf failed");

    size_t i = 0;
    for(; *rng_types; rng_types++)
    {
        if(++i % 4)
            rc = fprintf(file, "%-18s", (*rng_types)->name);
        else
            rc = fprintf(file, "%-18s\n", (*rng_types)->name);
        libcheck(rc >= 0, "fprintf failed");
    }
    if(!(i % 4 == 0))
        libcheck(fputc('\n', file) >= 0, "fputc failed");

    return EXIT_SUCCESS;

error:
    return EXIT_FAILURE;
}

gsl_rng* rng_alloc_and_seed(const gsl_rng_type* type, unsigned long seed)
{
    gsl_rng* rng = gsl_rng_alloc(type);
    libcheck_mem(rng);

    gsl_rng_set(rng, seed);
    return rng;

error:
    return NULL;
}

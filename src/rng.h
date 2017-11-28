/* rng.h
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

#ifndef FB_LIBLATTICE_RNG_H
#define FB_LIBLATTICE_RNG_H

#include <gsl/gsl_rng.h>

unsigned long get_random_seed(void);

const gsl_rng_type* get_rng_type(const char* rng_name, const gsl_rng_type** rng_types);

int print_rngs(FILE* file, const gsl_rng_type** rng_types);

gsl_rng* rng_alloc_and_seed(const gsl_rng_type* type, unsigned long seed);

#endif /* FB_LIBLATTICE_RNG_H */

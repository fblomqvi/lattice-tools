/* algorithm.h
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

#ifndef FB_LIBLATTICE_ALGORITHM_H
#define FB_LIBLATTICE_ALGORITHM_H

#include <gsl/gsl_matrix.h>

typedef enum algorithm
{
    ALG_BABAI = 0,
    ALG_DPLANE,
    ALG_DPLANE_ITER,
    ALG_DPLANE_VANILLA,
    ALG_SPHERE_POHST,
    ALG_SPHERE_DP,
    ALG_SPHERE_SE,
    ALG_MAX
} Algorithm;

typedef void (*SOLVE_func)(gsl_vector*, const gsl_vector*, const gsl_matrix*, void*);

int algorithm_parse_name(const char* name);
const char* algorithm_get_name(Algorithm alg);
int algorithm_print_names(FILE* file);

int algorithm_get_fp_init_ws(SOLVE_func* f, void** ws,
                            Algorithm alg, const gsl_matrix* basis);

void algorithm_free_ws(void* ws, Algorithm alg);


#endif /* FB_LIBLATTICE_ALGORITHM_H */

/* lattice_gen.h
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

#ifndef FB_LIBLATTICE_LATTICE_GEN_H
#define FB_LIBLATTICE_LATTICE_GEN_H

#include "matrix_mpz.h"
#include <gsl/gsl_rng.h>

typedef enum l_type
{
    LTYPE_A = 1,
    LTYPE_D,
    LTYPE_D_DUAL,
    LTYPE_RANDOM,
    LTYPE_SPC,
} L_type;

typedef struct s_lgen_params
{
    gsl_rng* rng;
    size_t dimension;
    size_t exponent;
    unsigned long range;
    long offset;
} LGEN_PARAMS;

static inline int lattice_gen_type_needs_rng(L_type type)
{ return (type >= LTYPE_RANDOM); }

t_MAT_MPZ* lattice_gen(L_type type, LGEN_PARAMS* params);

t_MAT_MPZ* lattice_gen_random(L_type type, LGEN_PARAMS* params);
t_MAT_MPZ* lattice_gen_standard(L_type type, size_t dimension, size_t exponent);

t_MAT_MPZ* lattice_gen_Anm(size_t dimension, size_t exponent);
t_MAT_MPZ* lattice_gen_D(size_t dimension, size_t exponent);
t_MAT_MPZ* lattice_gen_D_dual(size_t dimension, size_t exponent);

t_MAT_MPZ* lattice_gen_random_square(LGEN_PARAMS* params);
t_MAT_MPZ* lattice_gen_random_spc(LGEN_PARAMS* params);

t_MAT_MPZ* lattice_gen_random_square_gmp(size_t dimension, size_t exponent, 
                                        unsigned long seed, size_t bits);

#endif /* FB_LIBLATTICE_LATTICE_GEN_H */

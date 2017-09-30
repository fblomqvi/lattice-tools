/* lattice_gen.c
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
#include "lattice_gen.h"
#include "gmp.h"

static t_MAT_MPZ* lgen_An_helper(size_t n);
static t_MAT_MPZ* lgen_Anm_helper(size_t n, size_t m);
static int MAT_MPZ_power(t_MAT_MPZ** M, size_t exponent);

t_MAT_MPZ* lattice_gen(L_type type, LGEN_PARAMS* params)
{
    switch(type)
    {
        case LTYPE_A: return lattice_gen_Anm(params->dimension, params->exponent);
        case LTYPE_D: return lattice_gen_D(params->dimension, params->exponent);
        case LTYPE_D_DUAL: return lattice_gen_D_dual(params->dimension, params->exponent);
        case LTYPE_RANDOM: return lattice_gen_random_square(params);
        case LTYPE_SPC: return lattice_gen_random_spc(params);
        default: return NULL;
    }
}

t_MAT_MPZ* lattice_gen_random(L_type type, LGEN_PARAMS* params)
{
    switch(type)
    {
        case LTYPE_RANDOM: return lattice_gen_random_square(params);
        case LTYPE_SPC: return lattice_gen_random_spc(params);
        default: return NULL;
    }
}

t_MAT_MPZ* lattice_gen_standard(L_type type, size_t dimension, size_t exponent)
{
    switch(type)
    {
        case LTYPE_A: return lattice_gen_Anm(dimension, exponent);
        case LTYPE_D: return lattice_gen_D(dimension, exponent);
        case LTYPE_D_DUAL: return lattice_gen_D_dual(dimension, exponent);
        default: return NULL;
    }
}

t_MAT_MPZ* lattice_gen_Anm(size_t dimension, size_t exponent)
{
    if(exponent > 1)
        return lgen_Anm_helper(dimension, exponent);
    else
        return lgen_An_helper(dimension);
}

t_MAT_MPZ* lattice_gen_D(size_t dimension, size_t exponent)
{
    t_MAT_MPZ* M = MAT_MPZ_alloc(dimension, dimension);
    libcheck_mem(M);

    mpz_set_si(MAT_MPZ_get(M, 0, 0), -1);
    mpz_set_si(MAT_MPZ_get(M, 0, 1), -1);
    for(size_t i = 1; i < dimension; i++)
    {
        mpz_set_si(MAT_MPZ_get(M, i, i-1), 1);
        mpz_set_si(MAT_MPZ_get(M, i, i), -1);
    }

    if(exponent > 1)
        MAT_MPZ_power(&M, exponent);

    return M;

error:
    return NULL;
}

t_MAT_MPZ* lattice_gen_D_dual(size_t dimension, size_t exponent)
{
    t_MAT_MPZ* M = MAT_MPZ_alloc(dimension, dimension);
    libcheck_mem(M);

    for(size_t i = 0; i < dimension-1; i++)
        mpz_set_si(MAT_MPZ_get(M, i, i), 2);

    for(size_t i = 0; i < dimension; i++)
        mpz_set_si(MAT_MPZ_get(M, dimension-1, i), 1);

    if(exponent > 1)
        MAT_MPZ_power(&M, exponent);

    return M;

error:
    return NULL;
}

t_MAT_MPZ* lattice_gen_random_square(LGEN_PARAMS* params)
{
    t_MAT_MPZ* M = MAT_MPZ_alloc(params->dimension, params->dimension);
    libcheck_mem(M);

    for(size_t i = 0; i < params->dimension; i++)
        for(size_t j = 0; j < params->dimension; j++)
        {
            long val = params->offset + gsl_rng_uniform_int(params->rng, params->range);
            mpz_set_si(MAT_MPZ_get(M, i, j), val);
        }

    if(params->exponent > 1)
        MAT_MPZ_power(&M, params->exponent);

    return M;

error:
    return NULL;
}

t_MAT_MPZ* lattice_gen_random_spc(LGEN_PARAMS* params)
{
    t_MAT_MPZ* M = MAT_MPZ_alloc(params->dimension, params->dimension);
    libcheck_mem(M);

    for(size_t i = 0; i < params->dimension - 1; i++)
    {
        mpz_set_ui(MAT_MPZ_get(M, i, i), 1);
        mpz_set_ui(MAT_MPZ_get(M, i, params->dimension - 1),
                1 + gsl_rng_uniform_int(params->rng, params->range - 1));
    }
    mpz_set_ui(MAT_MPZ_get(M, params->dimension - 1, params->dimension - 1), params->range);

    if(params->exponent > 1)
        MAT_MPZ_power(&M, params->exponent);

    return M;

error:
    return NULL;
}

t_MAT_MPZ* lattice_gen_random_square_gmp(size_t dimension, size_t exponent, 
                                        unsigned long seed, size_t bits)
{
    t_MAT_MPZ* M = MAT_MPZ_alloc(dimension, dimension);
    libcheck_mem(M);

    mpz_t offset, range;
    mpz_inits(offset, range, NULL);
    mpz_set_ui(range, 1);
    mpz_mul_2exp(range, range, bits);
    mpz_ui_sub(offset, 1, range);
    mpz_sub(range, range, offset);
    debug_do(gmp_fprintf(stderr, "offset: %Zd\n", offset));
    debug_do(gmp_fprintf(stderr, "range: %Zd\n", range));

    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, seed);

    for(size_t i = 0; i < dimension; i++)
        for(size_t j = 0; j < dimension; j++)
        {
            mpz_urandomm(MAT_MPZ_get(M, i, j), state, range);
            mpz_add(MAT_MPZ_get(M, i, j), MAT_MPZ_get(M, i, j), offset);
        }

    if(exponent > 1)
        MAT_MPZ_power(&M, exponent);

    gmp_randclear(state);
    mpz_clears(offset, range, NULL);
    return M;

error:
    return NULL;
}

static t_MAT_MPZ* lgen_An_helper(size_t n)
{
    t_MAT_MPZ* M = MAT_MPZ_alloc(n, n + 1);
    libcheck_mem(M);

    for(size_t i = 0; i < n; i++)
    {
        mpz_set_si(MAT_MPZ_get(M, i, i), 1);
        mpz_set_si(MAT_MPZ_get(M, i, i+1), -1);
    }

    return M;

error:
    return NULL;
}

static t_MAT_MPZ* lgen_Anm_helper(size_t n, size_t m)
{
    t_MAT_MPZ* M = MAT_MPZ_alloc(n + 1, n + 1);
    libcheck_mem(M);

    for(size_t i = 0; i < n; i++)
    {
        mpz_set_si(MAT_MPZ_get(M, i, i), 1);
        mpz_set_si(MAT_MPZ_get(M, i, i+1), -1);
    }
    mpz_set_si(MAT_MPZ_get(M, n, n), 1);
    mpz_set_si(MAT_MPZ_get(M, n, 0), -1);

    MAT_MPZ_power(&M, m);
    MAT_MPZ_resize(M, n, n + 1);

    return M;

error:
    return NULL;
}

/* M must be a square matrix. The exponentiation is done by brute force... */
static int MAT_MPZ_power(t_MAT_MPZ** M, size_t exponent)
{
    t_MAT_MPZ* R = MAT_MPZ_alloc(MAT_MPZ_rows(*M), MAT_MPZ_cols(*M));
    libcheck_mem(R);
    t_MAT_MPZ* T = MAT_MPZ_alloc(MAT_MPZ_rows(*M), MAT_MPZ_cols(*M));
    llibcheck_mem(T, error_a);

    MAT_MPZ_cpy(T, *M);
    for(size_t i = 0; i < exponent - 1; i++)
    {
        MAT_MPZ_mul(R, *M, T);
        t_MAT_MPZ* T = *M;
        *M = R;
        R = T;
    }
     
    MAT_MPZ_free(T);
    MAT_MPZ_free(R);
    return 0;

error_a:
    MAT_MPZ_free(R);
error:
    return -1;
}

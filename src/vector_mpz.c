/* vector_mpz.c
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
#include "vector_mpz.h"
#include <assert.h>

t_VEC_MPZ* VEC_MPZ_alloc2(size_t max_size, size_t size)
{
    t_VEC_MPZ* vec = malloc(sizeof(t_VEC_MPZ));
    libcheck_mem(vec);

    vec->size = size;
    vec->max_size = max_size;
    vec->data = malloc(max_size * sizeof(mpz_t));
    llibcheck_mem(vec->data, error_a);

    for(size_t i = 0; i < vec->max_size; i++)
        mpz_init(vec->data[i]);

    return vec;

error_a:
    free(vec);
error:
    return NULL;
}

void VEC_MPZ_free(t_VEC_MPZ* vec)
{
    if(!vec)
        return;

    for(size_t i = 0; i < vec->max_size; i++)
        mpz_clear(vec->data[i]);

    free(vec->data);
    free(vec);
}

void VEC_MPZ_cpy(t_VEC_MPZ* dest, t_VEC_MPZ* src)
{
    assert(dest->max_size >= src->size);
    for(size_t i = 0; i < src->size; i++)
        mpz_set(dest->data[i], src->data[i]);
    dest->size = src->size;
}

void VEC_MPZ_do_op1(t_VEC_MPZ* r, const t_VEC_MPZ* a,
                    void (*op)(mpz_t, const mpz_t))
{
    assert(r->size == a->size);
    for(size_t i = 0; i < a->size; i++)
        op(r->data[i], a->data[i]);
}

void VEC_MPZ_do_op2(t_VEC_MPZ* r, const t_VEC_MPZ* a, const t_VEC_MPZ* b,
                    void (*op)(mpz_t, const mpz_t, const mpz_t))
{
    assert(r->size == a->size && a->size == b->size);
    for(size_t i = 0; i < a->size; i++)
        op(r->data[i], a->data[i], b->data[i]);
}

void VEC_MPZ_do_op3(t_VEC_MPZ* r, const t_VEC_MPZ* a, const mpz_t b,
                    void (*op)(mpz_t, const mpz_t, const mpz_t))
{
    assert(r->size == a->size);
    for(size_t i = 0; i < a->size; i++)
        op(r->data[i], a->data[i], b);
}

/*
void VEC_MPZ_dot(mpz_t r, const t_VEC_MPZ* a, const t_VEC_MPZ* b)
{
    assert(a->size == b->size);
    mpz_mul(r, a->data[0], b->data[0]);
    for(size_t i = 1; i < a->size; i++)
        mpz_addmul(r, a->data[i], b->data[i]);
}
*/

int VEC_MPZ_equal(t_VEC_MPZ* a, t_VEC_MPZ* b)
{
    assert(a->size == b->size);
    for(size_t i = 0; i < a->size; i++)
        if(mpz_cmp(a->data[i], b->data[i]))
            return 0;

    return 1;
}

void VEC_MPZ_print(FILE* file, const t_VEC_MPZ* vec)
{
    for(size_t i = 0; i < vec->size; i++)
        gmp_fprintf(file, "%Zd\n", vec->data[i]);
}

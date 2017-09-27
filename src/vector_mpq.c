/* vector_mpq.c
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
#include "vector_mpq.h"
#include <assert.h>

static void mpq_addmul(mpq_t r, const mpq_t a, const mpq_t b, mpq_t ws);
static void mpq_submul(mpq_t r, const mpq_t a, const mpq_t b, mpq_t ws);

static void lcm_denom(mpq_t r, const t_VEC_MPQ* a, size_t size);

t_VEC_MPQ* VEC_MPQ_alloc2(size_t max_size, size_t size)
{
    t_VEC_MPQ* vec = malloc(sizeof(t_VEC_MPQ));
    libcheck_mem(vec);
   
    vec->size = size;
    vec->max_size = max_size;
    vec->data = malloc(max_size * sizeof(mpq_t));
    llibcheck_mem(vec->data, error_a);

    for(size_t i = 0; i < vec->max_size; i++)
        mpq_init(vec->data[i]);

    return vec;

error_a:
    free(vec);
error:
    return NULL;
}

void VEC_MPQ_free(t_VEC_MPQ* vec)
{
    if(!vec)
        return;

    for(size_t i = 0; i < vec->max_size; i++)
        mpq_clear(vec->data[i]);

    free(vec->data);
    free(vec);
}

void VEC_MPQ_cpy(t_VEC_MPQ* dest, t_VEC_MPQ* src)
{
    assert(dest->max_size >= src->size);
    for(size_t i = 0; i < src->size; i++)
        mpq_set(dest->data[i], src->data[i]);
    dest->size = src->size;
}

void VEC_MPQ_do_op1(t_VEC_MPQ* r, const t_VEC_MPQ* a,
                    void (*op)(mpq_t, const mpq_t))
{
    assert(r->size == a->size);
    for(size_t i = 0; i < a->size; i++)
        op(r->data[i], a->data[i]);
}

void VEC_MPQ_do_op2(t_VEC_MPQ* r, const t_VEC_MPQ* a, const t_VEC_MPQ* b,
                    void (*op)(mpq_t, const mpq_t, const mpq_t))
{
    assert(r->size == a->size && a->size == b->size);
    for(size_t i = 0; i < a->size; i++)
        op(r->data[i], a->data[i], b->data[i]);
}

void VEC_MPQ_do_op3(t_VEC_MPQ* r, const t_VEC_MPQ* a, const mpq_t b,
                    void (*op)(mpq_t, const mpq_t, const mpq_t))
{
    assert(r->size == a->size);
    for(size_t i = 0; i < a->size; i++)
        op(r->data[i], a->data[i], b);
}

void VEC_MPQ_dot(mpq_t r, const t_VEC_MPQ* a, const t_VEC_MPQ* b, mpq_t ws)
{
    assert(a->size == b->size);
    mpq_mul(r, a->data[0], b->data[0]);
    for(size_t i = 1; i < a->size; i++)
        mpq_addmul(r, a->data[i], b->data[i], ws);
}

int VEC_MPQ_equal(t_VEC_MPQ* a, t_VEC_MPQ* b)
{
    assert(a->size == b->size);
    for(size_t i = 0; i < a->size; i++)
        if(!mpq_equal(a->data[i], b->data[i]))
            return 0;

    return 1;
}

void VEC_MPQ_addmul(t_VEC_MPQ* r, const t_VEC_MPQ* a, const mpq_t b, mpq_t ws)
{
    assert(r->size == a->size);
    for(size_t i = 0; i < a->size; i++)
        mpq_addmul(r->data[i], a->data[i], b, ws);
}

void VEC_MPQ_submul(t_VEC_MPQ* r, const t_VEC_MPQ* a, const mpq_t b, mpq_t ws)
{
    assert(r->size == a->size);
    for(size_t i = 0; i < a->size; i++)
        mpq_submul(r->data[i], a->data[i], b, ws);
}

void VEC_MPQ_print(FILE* file, const t_VEC_MPQ* vec)
{
    for(size_t i = 0; i < vec->size; i++)
        gmp_fprintf(file, "%Qd ", vec->data[i]);
    fprintf(file, "\n");
}

void VEC_MPQ_clear_denominators(t_VEC_MPQ* r, const t_VEC_MPQ* a, mpq_t ws)
{
    assert(a->size >= 2);
    lcm_denom(ws, a, a->size);
    VEC_MPQ_mul(r, a, ws);
}

void VEC_MPQ_clear_denominators_size(t_VEC_MPQ* r, const t_VEC_MPQ* a, 
                                    mpq_t ws, size_t len)
{
    assert(len >= 2 && len <= a->max_size);
    lcm_denom(ws, a, len);
    VEC_MPQ_mul_size(r, a, ws, len);
}

void VEC_MPQ_mul_size(t_VEC_MPQ* r, const t_VEC_MPQ* a, const mpq_t b, size_t len)
{
    assert(a->max_size >= len && r->max_size >= len);
    for(size_t i = 0; i < len; i++)
        mpq_mul(r->data[i], a->data[i], b);
}

static void mpq_addmul(mpq_t r, const mpq_t a, const mpq_t b, mpq_t ws)
{
    mpq_mul(ws, a, b);
    mpq_add(r, r, ws);
}

static void mpq_submul(mpq_t r, const mpq_t a, const mpq_t b, mpq_t ws)
{
    mpq_mul(ws, a, b);
    mpq_sub(r, r, ws);
}

static void lcm_denom(mpq_t r, const t_VEC_MPQ* a, size_t size)
{
    mpz_set_ui(mpq_denref(r), 1);
    mpz_lcm(mpq_numref(r), mpq_denref(a->data[0]), mpq_denref(a->data[1]));
    for(size_t i = 2; i < size; i++)
        mpz_lcm(mpq_numref(r), mpq_numref(r), mpq_denref(a->data[i]));

}

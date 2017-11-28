/* vector_mpz.h
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

#ifndef FB_VECTOR_MPZ_H
#define FB_VECTOR_MPZ_H

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <assert.h>

typedef struct s_vector_mpz
{
    size_t size;
    size_t max_size;
    mpz_t* data;
} t_VEC_MPZ;

t_VEC_MPZ* VEC_MPZ_alloc2(size_t max_size, size_t size);

static inline t_VEC_MPZ* VEC_MPZ_alloc(size_t size)
{ return VEC_MPZ_alloc2(size, size); }

void VEC_MPZ_free(t_VEC_MPZ* vec);

// Assumes that there is enough space in dest.
void VEC_MPZ_cpy(t_VEC_MPZ* dest, t_VEC_MPZ* src);

// No checks!
static inline void VEC_MPZ_resize(t_VEC_MPZ* vec, size_t size)
{ vec->size = size; }

void VEC_MPZ_do_op1(t_VEC_MPZ* r, const t_VEC_MPZ* a,
                    void (*op)(mpz_t, const mpz_t));
void VEC_MPZ_do_op2(t_VEC_MPZ* r, const t_VEC_MPZ* a, const t_VEC_MPZ* b,
                    void (*op)(mpz_t, const mpz_t, const mpz_t));

void VEC_MPZ_do_op3(t_VEC_MPZ* r, const t_VEC_MPZ* a, const mpz_t b,
                    void (*op)(mpz_t, const mpz_t, const mpz_t));

int VEC_MPZ_equal(t_VEC_MPZ* a, t_VEC_MPZ* b);

// Calculates r = <a,b>
static inline void VEC_MPZ_dot_size(mpz_t r, const t_VEC_MPZ* a,
                                    const t_VEC_MPZ* b, size_t len)
{
    assert(a->size >= len && b->size >= len);
    mpz_mul(r, a->data[0], b->data[0]);
    for(size_t i = 1; i < len; i++)
        mpz_addmul(r, a->data[i], b->data[i]);
}

// Calculates r = <a,b>
static inline void VEC_MPZ_dot(mpz_t r, const t_VEC_MPZ* a, const t_VEC_MPZ* b)
{
    assert(a->size == b->size);
    VEC_MPZ_dot_size(r, a, b, a->size);
}

// Calculates r = a + b
static inline void VEC_MPZ_add(t_VEC_MPZ* r, const t_VEC_MPZ* a, const t_VEC_MPZ* b)
{ VEC_MPZ_do_op2(r, a, b, mpz_add); }

// Calculates r = a + b
static inline void VEC_MPZ_sub(t_VEC_MPZ* r, const t_VEC_MPZ* a, const t_VEC_MPZ* b)
{ VEC_MPZ_do_op2(r, a, b, mpz_sub); }

// Calculates r = b * a
static inline void VEC_MPZ_mul(t_VEC_MPZ* r, const t_VEC_MPZ* a, const mpz_t b)
{ VEC_MPZ_do_op3(r, a, b, mpz_mul); }

// Calculates r += b * a
static inline void VEC_MPZ_addmul(t_VEC_MPZ* r, const t_VEC_MPZ* a, const mpz_t b)
{ VEC_MPZ_do_op3(r, a, b, mpz_addmul); }

// Calculates r -= b * a
static inline void VEC_MPZ_submul(t_VEC_MPZ* r, const t_VEC_MPZ* a, const mpz_t b)
{ VEC_MPZ_do_op3(r, a, b, mpz_submul); }

// Calculates r -= b * a
static inline void VEC_MPZ_neg(t_VEC_MPZ* r, const t_VEC_MPZ* a)
{ VEC_MPZ_do_op1(r, a, mpz_neg); }

static inline void VEC_MPZ_norm_sqr(mpz_t r, const t_VEC_MPZ* a)
{ VEC_MPZ_dot(r, a, a); }


void VEC_MPZ_print(FILE* file, const t_VEC_MPZ* vec);

#define VEC_MPZ_get(v, i) ((v)->data[(i)])
/*
static inline mpz_t* VEC_MPZ_get(const t_VEC_MPZ* vec, size_t i)
{ return vec->data + i; }
*/

static inline size_t VEC_MPZ_size(const t_VEC_MPZ* vec)
{ return vec->size; }

static inline size_t VEC_MPZ_max_size(const t_VEC_MPZ* vec)
{ return vec->max_size; }

static inline void VEC_MPZ_swap_elems(t_VEC_MPZ* vec, size_t i, size_t j)
{ mpz_swap(vec->data[i], vec->data[j]); }

#endif /* FB_VECTOR_MPZ_H */

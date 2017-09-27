/* vehtor_mpq.c
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

#ifndef FB_VECTOR_MPQ_H
#define FB_VECTOR_MPQ_H

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

typedef struct s_vector_mpq
{
    size_t size;
    size_t max_size;
    mpq_t* data;
} t_VEC_MPQ;

t_VEC_MPQ* VEC_MPQ_alloc2(size_t max_size, size_t size);

static inline t_VEC_MPQ* VEC_MPQ_alloc(size_t size)
{ return VEC_MPQ_alloc2(size, size); }

void VEC_MPQ_free(t_VEC_MPQ* vec);

// Assumes that there is enough space in dest.
void VEC_MPQ_cpy(t_VEC_MPQ* dest, t_VEC_MPQ* src);

// No checks!
static inline void VEC_MPQ_resize(t_VEC_MPQ* vec, size_t size)
{ vec->size = size; }

void VEC_MPQ_do_op1(t_VEC_MPQ* r, const t_VEC_MPQ* a,
                    void (*op)(mpq_t, const mpq_t));
void VEC_MPQ_do_op2(t_VEC_MPQ* r, const t_VEC_MPQ* a, const t_VEC_MPQ* b,
                    void (*op)(mpq_t, const mpq_t, const mpq_t));

void VEC_MPQ_do_op3(t_VEC_MPQ* r, const t_VEC_MPQ* a, const mpq_t b,
                    void (*op)(mpq_t, const mpq_t, const mpq_t));

int VEC_MPQ_equal(t_VEC_MPQ* a, t_VEC_MPQ* b);

// Calculates r = <a,b>
void VEC_MPQ_dot(mpq_t r, const t_VEC_MPQ* a, const t_VEC_MPQ* b, mpq_t ws);

// Calculates r = a + b
static inline void VEC_MPQ_add(t_VEC_MPQ* r, const t_VEC_MPQ* a, const t_VEC_MPQ* b)
{ VEC_MPQ_do_op2(r, a, b, mpq_add); }

// Calculates r = a + b
static inline void VEC_MPQ_sub(t_VEC_MPQ* r, const t_VEC_MPQ* a, const t_VEC_MPQ* b)
{ VEC_MPQ_do_op2(r, a, b, mpq_sub); }

// Calculates r = b * a
static inline void VEC_MPQ_mul(t_VEC_MPQ* r, const t_VEC_MPQ* a, const mpq_t b)
{ VEC_MPQ_do_op3(r, a, b, mpq_mul); }

void VEC_MPQ_mul_size(t_VEC_MPQ* r, const t_VEC_MPQ* a, const mpq_t b, size_t len);

// Calculates r = 1/b * a
static inline void VEC_MPQ_div(t_VEC_MPQ* r, const t_VEC_MPQ* a, const mpq_t b)
{ VEC_MPQ_do_op3(r, a, b, mpq_div); }

// Calculates r += b * a
void VEC_MPQ_addmul(t_VEC_MPQ* r, const t_VEC_MPQ* a, const mpq_t b, mpq_t ws);

// Calculates r -= b * a
void VEC_MPQ_submul(t_VEC_MPQ* r, const t_VEC_MPQ* a, const mpq_t b, mpq_t ws);

// Calculates r -= b * a
static inline void VEC_MPQ_neg(t_VEC_MPQ* r, const t_VEC_MPQ* a)
{ VEC_MPQ_do_op1(r, a, mpq_neg); }

static inline void VEC_MPQ_norm_sqr(mpq_t r, const t_VEC_MPQ* a, mpq_t ws)
{ VEC_MPQ_dot(r, a, a, ws); }

// Assumes that a has length at least 2.
void VEC_MPQ_clear_denominators(t_VEC_MPQ* r, const t_VEC_MPQ* a, mpq_t ws);
void VEC_MPQ_clear_denominators_size(t_VEC_MPQ* r, const t_VEC_MPQ* a, 
                                    mpq_t ws, size_t len);

void VEC_MPQ_print(FILE* file, const t_VEC_MPQ* vec);


#define VEC_MPQ_get(v, i) ((v)->data[(i)])
/*
static inline mpq_t* VEC_MPQ_get(const t_VEC_MPQ* vec, size_t i)
{ return vec->data + i; }
*/

static inline size_t VEC_MPQ_size(const t_VEC_MPQ* vec)
{ return vec->size; }

static inline size_t VEC_MPQ_max_size(const t_VEC_MPQ* vec)
{ return vec->max_size; }

static inline void VEC_MPQ_swap_elems(t_VEC_MPQ* vec, size_t i, size_t j)
{ mpq_swap(vec->data[i], vec->data[j]); }

#endif /* FB_VECTOR_MPQ_H */

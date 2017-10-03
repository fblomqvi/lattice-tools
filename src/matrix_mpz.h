/* matrix_mpz.h
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

#ifndef FB_MATRIX_MPZ_H
#define FB_MATRIX_MPZ_H

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "vector_mpz.h"
#include "vector_mpq.h"
#include "print_util.h"

typedef struct s_matrix_mpz
{
    size_t rows;
    size_t cols;
    size_t max_rows;
    size_t max_cols;
    t_VEC_MPZ** columns; 
    int owner;
} t_MAT_MPZ;

t_MAT_MPZ* MAT_MPZ_alloc2(size_t max_rows, size_t max_cols, size_t rows, size_t cols);

static inline t_MAT_MPZ* MAT_MPZ_alloc(size_t rows, size_t cols)
{ return MAT_MPZ_alloc2(rows, cols, rows, cols); }

t_MAT_MPZ* MAT_MPZ_alloc_empty(size_t max_rows, size_t max_cols);
t_MAT_MPZ* MAT_MPZ_alloc_X(size_t size);

void MAT_MPZ_free(t_MAT_MPZ* M);

void MAT_MPZ_resize(t_MAT_MPZ* M, size_t rows, size_t cols);

void MAT_MPZ_cpy(t_MAT_MPZ* dest, const t_MAT_MPZ* src);
void MAT_MPZ_shallow_cpy(t_MAT_MPZ* dest, const t_MAT_MPZ* src);

// Do the elementwise operation op
void MAT_MPZ_do_op1(t_MAT_MPZ* R, const t_MAT_MPZ* A,
                    void (*op)(mpz_t, const mpz_t));

// Do the elementwise operation op
void MAT_MPZ_do_op2(t_MAT_MPZ* R, const t_MAT_MPZ* A, const t_MAT_MPZ* B,
                    void (*op)(mpz_t, const mpz_t, const mpz_t));

// Do the elementwise operation op
void MAT_MPZ_do_op3(t_MAT_MPZ* R, const t_MAT_MPZ* A, const mpz_t beta,
                    void (*op)(mpz_t, const mpz_t, const mpz_t));

// Set R to A^T
void MAT_MPZ_trans(t_MAT_MPZ* R, const t_MAT_MPZ* A);

// Transpose R in-place. Only works with square matrices.
void MAT_MPZ_trans_ip(t_MAT_MPZ* R);

int MAT_MPZ_equal(t_MAT_MPZ* A, t_MAT_MPZ* B);

// Set R to -A
static inline void MAT_MPZ_neg(t_MAT_MPZ* R, const t_MAT_MPZ* A)
{ MAT_MPZ_do_op1(R, A, mpz_neg); }

// Set R to A + B
static inline void MAT_MPZ_add(t_MAT_MPZ* R, const t_MAT_MPZ* A, const t_MAT_MPZ* B)
{ MAT_MPZ_do_op2(R, A, B, mpz_add); }

// Set R to A - B
static inline void MAT_MPZ_sub(t_MAT_MPZ* R, const t_MAT_MPZ* A, const t_MAT_MPZ* B)
{ MAT_MPZ_do_op2(R, A, B, mpz_sub); }

// Set R to beta * A
static inline void MAT_MPZ_scale_mul(t_MAT_MPZ* R, const t_MAT_MPZ* A, const mpz_t beta)
{ MAT_MPZ_do_op3(R, A, beta, mpz_mul); }

// Set r to A * b
void MAT_MPZ_mulvec(t_VEC_MPZ* r, const t_MAT_MPZ* A, const t_VEC_MPZ* b);

// Set r to A^T * b
static inline void MAT_MPZ_mulvec_trans_size(t_VEC_MPZ* r, const t_MAT_MPZ* A, 
                                const t_VEC_MPZ* b, size_t A_rows, size_t A_cols)
{
    for(size_t i = 0; i < A_cols; i++)
        VEC_MPZ_dot_size(VEC_MPZ_get(r, i), A->columns[i], b, A_rows);
}

static inline void MAT_MPZ_mulvec_trans(t_VEC_MPZ* r, const t_MAT_MPZ* A, const t_VEC_MPZ* b)
{ MAT_MPZ_mulvec_trans_size(r, A, b, A->rows, A->cols); }

// Set R to A * B
void MAT_MPZ_mul(t_MAT_MPZ* R, const t_MAT_MPZ* A, const t_MAT_MPZ* B);

// Set R to A^T * B
void MAT_MPZ_mul_trans(t_MAT_MPZ* R, const t_MAT_MPZ* A, const t_MAT_MPZ* B);
void MAT_MPZ_mul_trans_size(t_MAT_MPZ* R, const t_MAT_MPZ* A, const t_MAT_MPZ* B,
                            size_t A_rows, size_t A_cols, size_t B_cols);

// Assumes that B is upper triangular. Does not work if A_rows is smaller than B_cols. 
// Or something like that
static inline void MAT_MPZ_mul_trans_UT_size(t_MAT_MPZ* R, const t_MAT_MPZ* A, 
                                            const t_MAT_MPZ* B, size_t A_rows, 
                                            size_t A_cols, size_t B_cols)
{
    for(size_t i = 0; i < B_cols; i++)
        MAT_MPZ_mulvec_trans_size(R->columns[i], A, B->columns[i], i+1, A_cols);
}

void MAT_MPZ_print_no_trans(FILE* file, const t_MAT_MPZ* M, const PRINTING_FMT* fmt);
void MAT_MPZ_print_trans(FILE* file, const t_MAT_MPZ* M, const PRINTING_FMT* fmt);
void MAT_MPZ_print(FILE* file, const t_MAT_MPZ* M, int transpose, const PRINTING_FMT* fmt);

void MAT_MPZ_solve_trans_size(t_VEC_MPZ* x, const t_MAT_MPZ* A, size_t A_rows, 
                            size_t A_cols, t_VEC_MPQ* ws, mpq_t wsq);

static inline void MAT_MPZ_solve_trans(t_VEC_MPZ* x, const t_MAT_MPZ* A, 
                                        t_VEC_MPQ* ws, mpq_t wsq)
{ MAT_MPZ_solve_trans_size(x, A, A->rows, A->cols, ws, wsq); }

// M must be square.
void MAT_MPZ_diag_prod(mpz_t det, const t_MAT_MPZ* M);

static inline size_t MAT_MPZ_rows(const t_MAT_MPZ* M)
{ return M->rows; }

static inline size_t MAT_MPZ_cols(const t_MAT_MPZ* M)
{ return M->cols; }

static inline size_t MAT_MPZ_max_rows(const t_MAT_MPZ* M)
{ return M->max_rows; }

static inline size_t MAT_MPZ_max_cols(const t_MAT_MPZ* M)
{ return M->max_cols; }

#define MAT_MPZ_get(M, r, c) VEC_MPZ_get((M)->columns[(c)], (r))
/*
static inline mpz_t* MAT_MPZ_get(const t_MAT_MPZ* M, size_t r, size_t c)
{ return VEC_MPZ_get(M->columns[c], r); }
*/

static inline void MAT_MPZ_swap_cols(t_MAT_MPZ* M, size_t i, size_t j)
{
    t_VEC_MPZ* temp = M->columns[i];
    M->columns[i] = M->columns[j];
    M->columns[j] = temp;
}

void MAT_MPZ_swap_rows(t_MAT_MPZ* M, size_t i, size_t j);

static inline t_VEC_MPZ* MAT_MPZ_get_col(const t_MAT_MPZ* M, size_t i)
{ return M->columns[i]; }

/* Simply disards the old column. It is up to the caller to manage the memory
 * of the old column. No bounds checking! */
static inline void MAT_MPZ_set_col(t_MAT_MPZ* M, size_t i, t_VEC_MPZ* v)
{ M->columns[i] = v; }

// Does not check that M has enough space for the new column.
static inline void MAT_MPZ_add_col(t_MAT_MPZ* M, t_VEC_MPZ* v)
{ M->columns[M->cols++] = v; }

// No bounds checking!
static inline t_VEC_MPZ* MAT_MPZ_del_col(t_MAT_MPZ* M, size_t i)
{ 
    t_VEC_MPZ* t = M->columns[i];
    M->columns[i] = M->columns[--M->cols];
    return t;
}

static inline void MAT_MPZ_set_ncol(t_MAT_MPZ* M, size_t num_cols)
{ M->cols = num_cols; }

#endif /* FB_MATRIX_MPZ_H */

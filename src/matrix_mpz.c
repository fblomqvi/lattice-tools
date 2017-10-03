/* matrix_mpz.c
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
#include "matrix_mpz.h"
#include <assert.h>

t_MAT_MPZ* MAT_MPZ_alloc2(size_t max_rows, size_t max_cols, size_t rows, size_t cols)
{
    t_MAT_MPZ* M = malloc(sizeof(t_MAT_MPZ));
    libcheck_mem(M);

    M->rows = rows;
    M->cols = cols;
    M->max_rows = max_rows;
    M->max_cols = max_cols;
    M->owner = 1;

    M->columns = calloc(max_cols, sizeof(t_VEC_MPZ*));
    llibcheck_mem(M->columns, error_a);

    for(size_t i = 0; i < M->max_cols; i++)
        llibcheck_mem(M->columns[i] = VEC_MPZ_alloc2(max_rows, rows), error_b);

    return M;

error_b:
    for(size_t i = 0; i < M->max_cols; i++)
        VEC_MPZ_free(M->columns[i]);
    free(M->columns);
error_a:
    free(M);
error:
    return NULL;
}

t_MAT_MPZ* MAT_MPZ_alloc_empty(size_t max_rows, size_t max_cols)
{
    t_MAT_MPZ* M = malloc(sizeof(t_MAT_MPZ));
    libcheck_mem(M);

    M->rows = max_rows;
    M->cols = 0;
    M->max_rows = max_rows;
    M->max_cols = max_cols;
    M->owner = 0;

    M->columns = malloc(max_cols * sizeof(t_VEC_MPZ*));
    llibcheck_mem(M->columns, error_a);

    return M;

error_a:
    free(M);
error:
    return NULL;
}

t_MAT_MPZ* MAT_MPZ_alloc_X(size_t size)
{
    t_MAT_MPZ* M = MAT_MPZ_alloc2(size, size, 2, 1);
    libcheck_mem(M);

    mpz_set_ui(MAT_MPZ_get(M, 0, 0), 1);
    return M;

error:
    return NULL;
}

void MAT_MPZ_free(t_MAT_MPZ* M)
{
    if(!M)
        return;

    if(M->owner)
    {
        for(size_t i = 0; i < M->max_cols; i++)
            VEC_MPZ_free(M->columns[i]);
    }
    free(M->columns);
    free(M);
}

void MAT_MPZ_resize(t_MAT_MPZ* M, size_t rows, size_t cols)
{
    M->cols = cols;
    M->rows = rows;
    for(size_t i = 0; i < M->cols; i++)
        VEC_MPZ_resize(M->columns[i], rows);
}

void MAT_MPZ_cpy(t_MAT_MPZ* dest, const t_MAT_MPZ* src)
{
    assert(dest->owner);
    assert(dest->max_rows >= src->rows && dest->max_cols >= src->cols);
    dest->rows = src->rows;
    dest->cols = src->cols;
    for(size_t i = 0; i < src->cols; i++)
        VEC_MPZ_cpy(dest->columns[i], src->columns[i]);
}

void MAT_MPZ_shallow_cpy(t_MAT_MPZ* dest, const t_MAT_MPZ* src)
{
    assert(dest->max_rows == src->max_rows && dest->max_cols == src->max_cols);
    dest->rows = src->rows;
    dest->cols = src->cols;
    dest->owner = 0;

    for(size_t i = 0; i < src->cols; i++)
        dest->columns[i] = src->columns[i];

}

void MAT_MPZ_do_op1(t_MAT_MPZ* R, const t_MAT_MPZ* A,
                    void (*op)(mpz_t, const mpz_t))
{
    assert(R->rows == A->rows && R->cols == A->cols);
    for(size_t i = 0; i < A->cols; i++)
        VEC_MPZ_do_op1(R->columns[i], A->columns[i], op);
}

void MAT_MPZ_do_op2(t_MAT_MPZ* R, const t_MAT_MPZ* A, const t_MAT_MPZ* B,
                    void (*op)(mpz_t, const mpz_t, const mpz_t))
{
    assert(R->rows == A->rows && R->cols == A->cols);
    assert(B->rows == A->rows && B->cols == A->cols);
    for(size_t i = 0; i < A->cols; i++)
        VEC_MPZ_do_op2(R->columns[i], A->columns[i], B->columns[i], op);
}

void MAT_MPZ_do_op3(t_MAT_MPZ* R, const t_MAT_MPZ* A, const mpz_t beta,
                    void (*op)(mpz_t, const mpz_t, const mpz_t))
{
    assert(R->rows == A->rows && R->cols == A->cols);
    for(size_t i = 0; i < A->cols; i++)
        VEC_MPZ_do_op3(R->columns[i], A->columns[i], beta, op);
}

void MAT_MPZ_trans(t_MAT_MPZ* R, const t_MAT_MPZ* A)
{
    assert(R->rows == A->cols && R->cols == A->rows);
    for(size_t i = 0; i < R->rows; i++)
        for(size_t j = 0; j < R->cols; j++)
            mpz_set(MAT_MPZ_get(R, i, j), MAT_MPZ_get(A, j, i));
}

void MAT_MPZ_trans_ip(t_MAT_MPZ* R)
{
    assert(R->rows == R->cols);
    for(size_t i = 0; i < R->rows; i++)
        for(size_t j = i + 1; j < R->cols; j++)
            mpz_swap(MAT_MPZ_get(R, i, j), MAT_MPZ_get(R, j, i));
}

int MAT_MPZ_equal(t_MAT_MPZ* A, t_MAT_MPZ* B)
{
    assert(A->cols == B->cols && A->rows == B->rows);
    for(size_t i = 0; i < A->cols; i++)
        if(!VEC_MPZ_equal(A->columns[i], B->columns[i]))
            return 0;

    return 1;
}

void MAT_MPZ_mul(t_MAT_MPZ* R, const t_MAT_MPZ* A, const t_MAT_MPZ* B)
{
    assert(A->cols == B->rows && R->rows == A->rows && R->cols == B->cols);
    for(size_t i = 0; i < B->cols; i++)
        MAT_MPZ_mulvec(R->columns[i], A, B->columns[i]);
}

void MAT_MPZ_mul_trans(t_MAT_MPZ* R, const t_MAT_MPZ* A, const t_MAT_MPZ* B)
{
    assert(A->rows == B->rows && R->rows == A->cols && R->cols == B->cols);
    for(size_t i = 0; i < B->cols; i++)
        MAT_MPZ_mulvec_trans(R->columns[i], A, B->columns[i]);
}

void MAT_MPZ_mul_trans_size(t_MAT_MPZ* R, const t_MAT_MPZ* A, const t_MAT_MPZ* B,
                            size_t A_rows, size_t A_cols, size_t B_cols)
{
    for(size_t i = 0; i < B_cols; i++)
        MAT_MPZ_mulvec_trans_size(R->columns[i], A, B->columns[i], A_rows, A_cols);
}

/*
void MAT_MPZ_mul_trans_UT_size(t_MAT_MPZ* R, const t_MAT_MPZ* A, const t_MAT_MPZ* B,
                            size_t A_rows, size_t A_cols, size_t B_cols)
{
    for(size_t i = 0; i < B_cols; i++)
        MAT_MPZ_mulvec_trans_size(R->columns[i], A, B->columns[i], i+1, A_cols);
}
*/

void MAT_MPZ_mulvec(t_VEC_MPZ* r, const t_MAT_MPZ* A, const t_VEC_MPZ* b)
{
    assert(A->cols == VEC_MPZ_size(b) && VEC_MPZ_size(r) == A->rows);
    VEC_MPZ_mul(r, A->columns[0], VEC_MPZ_get(b, 0));
    for(size_t i = 1; i < A->cols; i++)
        VEC_MPZ_addmul(r, A->columns[i], VEC_MPZ_get(b, i));
}

void MAT_MPZ_solve_trans_size(t_VEC_MPZ* x, const t_MAT_MPZ* A, size_t A_rows, 
                            size_t A_cols, t_VEC_MPQ* ws, mpq_t wsq)
{
    assert(A_cols == A_rows - 1); 
    assert(VEC_MPZ_max_size(x) >= A_rows && VEC_MPQ_max_size(ws) >= A_rows);
    size_t n = A_cols;
    mpq_set_ui(VEC_MPQ_get(ws, n), 1, 1);
    for(size_t k = 1; k <= A_cols; k++)
    {
        t_VEC_MPZ* column = A->columns[A_cols - k];
        mpz_set(mpq_numref(VEC_MPQ_get(ws, n-k)), VEC_MPZ_get(column, n));
        mpz_set(mpq_denref(VEC_MPQ_get(ws, n-k)), VEC_MPZ_get(column, n-k));
        mpq_canonicalize(VEC_MPQ_get(ws, n-k));
        // VEC_MPQ_get(ws, n) == 1 so no need to multiply
        //mpq_mul(wsq, wsq, VEC_MPQ_get(ws, n));
        mpq_neg(VEC_MPQ_get(ws, n-k), VEC_MPQ_get(ws, n-k));
        for(size_t i = 1; i < k; i++)
        {
            mpz_set(mpq_numref(wsq), VEC_MPZ_get(column, n-i));
            mpz_set(mpq_denref(wsq), VEC_MPZ_get(column, n-k));
            mpq_canonicalize(wsq);
            mpq_mul(wsq, wsq, VEC_MPQ_get(ws, n-i));
            mpq_sub(VEC_MPQ_get(ws, n-k), VEC_MPQ_get(ws, n-k), wsq);
        }
    }
    VEC_MPQ_clear_denominators_size(ws, ws, wsq, A_rows);
    for(size_t i = 0; i < A_rows; i++)
        mpz_set(VEC_MPZ_get(x, i), mpq_numref(VEC_MPQ_get(ws, i)));
}

void MAT_MPZ_diag_prod(mpz_t det, const t_MAT_MPZ* M)
{
    assert(M->rows == M->cols && M->rows > 1);
    mpz_mul(det, MAT_MPZ_get(M, 0, 0), MAT_MPZ_get(M, 1, 1));
    for(size_t i = 2; i < MAT_MPZ_cols(M); i++)
        mpz_mul(det, det, MAT_MPZ_get(M, i, i));
}

void MAT_MPZ_print(FILE* file, const t_MAT_MPZ* M, int transpose, const PRINTING_FMT* fmt)
{ 
    if(transpose)
        MAT_MPZ_print_trans(file, M, fmt);
    else
        MAT_MPZ_print_no_trans(file, M, fmt);
}

static void print_row(FILE* file, const t_MAT_MPZ* M, size_t row, const PRINTING_FMT* fmt)
{
    fprintf(file, "%s", fmt->vstart_bracket);
    for(size_t j = 0; j < M->cols-1; j++)
    {
        gmp_fprintf(file, "%Zd", VEC_MPZ_get(M->columns[j], row));
        fprintf(file, "%s", fmt->element_separator);
    }
    gmp_fprintf(file, "%Zd", VEC_MPZ_get(M->columns[M->cols-1], row));
    fprintf(file, "%s", fmt->vend_bracket);
}

void MAT_MPZ_print_no_trans(FILE* file, const t_MAT_MPZ* M, const PRINTING_FMT* fmt)
{
    fprintf(file, "%s", fmt->mstart_bracket);
    for(size_t i = 0; i < M->rows-1; i++)
    {
        print_row(file, M, i, fmt);
        fprintf(file, "%s", fmt->vector_separator);
    }
    print_row(file, M, M->rows-1, fmt);
    fprintf(file, "%s\n", fmt->mend_bracket);
}

static void print_column(FILE* file, const t_MAT_MPZ* M, size_t col, const PRINTING_FMT* fmt)
{
    fprintf(file, "%s", fmt->vstart_bracket);
    for(size_t j = 0; j < M->rows-1; j++)
    {
        gmp_fprintf(file, "%Zd", VEC_MPZ_get(M->columns[col], j));
        fprintf(file, "%s", fmt->element_separator);
    }
    gmp_fprintf(file, "%Zd", VEC_MPZ_get(M->columns[col], M->rows-1));
    fprintf(file, "%s", fmt->vend_bracket);
}

void MAT_MPZ_print_trans(FILE* file, const t_MAT_MPZ* M, const PRINTING_FMT* fmt)
{
    fprintf(file, "%s", fmt->mstart_bracket);
    for(size_t i = 0; i < M->cols-1; i++)
    {
        print_column(file, M, i, fmt);
        fprintf(file, "%s", fmt->vector_separator);
    }
    print_column(file, M, M->cols-1, fmt);
    fprintf(file, "%s\n", fmt->mend_bracket);
}

void MAT_MPZ_swap_rows(t_MAT_MPZ* M, size_t i, size_t j)
{
    for(size_t c = 0; c < M->cols; c++)
        VEC_MPZ_swap_elems(M->columns[c], i, j);
}

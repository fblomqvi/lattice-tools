/* babai.c
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
#include "lt_errno.h"
#include "babai.h"
#include "utility.h"
#include <math.h>
#include <assert.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

struct s_babai_ws
{
    gsl_matrix* R;
    double* data;
    double* s;
    double* y;
    gsl_vector_view v_s;
    gsl_vector_view v_y;
    gsl_matrix_view m_R;
    gsl_matrix_view Q1;
};

int BABAI_WS_alloc_and_init(BABAI_WS** ws_ptr, const gsl_matrix* B)
{
    int lt_errno = LT_SUCCESS;
    size_t n = B->size1;
    size_t m = B->size2;
    llibcheck_se(n >= m, error, lt_errno, LT_ELINDEP);

    BABAI_WS* ws = malloc(sizeof(BABAI_WS));
    llibcheck_se(ws, error, lt_errno, LT_ENOMEM);

    size_t Q_size = n * n;
    size_t R_size = n * m;
    ws->data = malloc((Q_size + R_size + 2 * m) * sizeof(double));
    llibcheck_se(ws->data, error_a, lt_errno, LT_ENOMEM);

    ws->s = ws->data + Q_size + R_size;
    ws->y = ws->s + m;
    ws->v_s = gsl_vector_view_array(ws->s, m);
    ws->v_y = gsl_vector_view_array(ws->y, m);

    ws->m_R = gsl_matrix_view_array(ws->data + Q_size, n, m);
    ws->R = &ws->m_R.matrix;

    gsl_matrix_view m_Q = gsl_matrix_view_array(ws->data, n, n);
    gsl_matrix* Q = &m_Q.matrix;

    lt_errno = utility_compute_QR_decomposition(Q, ws->R, B);
    llibcheck(lt_errno == 0, error_b, "utility_compute_QR_decomposition failed");
    llibcheck_se(utility_Rmm_is_not_singular(ws->R, 10E-10), error_b, lt_errno, LT_ELINDEP);

    ws->Q1 = gsl_matrix_submatrix(Q, 0, 0, n, m);

    *ws_ptr = ws;
    return lt_errno;

error_b:
    free(ws->data);
error_a:
    free(ws);
error:
    *ws_ptr = NULL;
    return lt_errno;
}

void BABAI_WS_free(BABAI_WS* ws)
{
    if(ws)
    {
        free(ws->data);
        free(ws);
    }
}

static double calc_yhat(size_t k, const gsl_matrix* R, 
                        const double* y, const double* s) 
{
    size_t m = R->size2;
    double sum = y[k];
    for (size_t j = k+1; j < m; j++)
        sum -= gsl_matrix_get(R, k, j) * s[j];
    
    return sum;
}

void babai(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, BABAI_WS* ws)
{
    gsl_blas_dgemv(CblasTrans, 1, &ws->Q1.matrix, t, 0, &ws->v_y.vector);
    
    for(int k = B->size2-1; k >= 0; k--)
        ws->s[k] = round(calc_yhat(k, ws->R, ws->y, ws->s) / gsl_matrix_get(ws->R, k, k));

    gsl_blas_dgemv(CblasNoTrans, 1, B, &ws->v_s.vector, 0, clp);
}

void babai_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{ babai(clp, t, B, ws); }

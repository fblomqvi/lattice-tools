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
#include "babai.h"
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

BABAI_WS* BABAI_WS_alloc_and_init(const gsl_matrix* B)
{
    BABAI_WS* ws = malloc(sizeof(BABAI_WS));
    libcheck_mem(ws);

    size_t n = B->size1;
    size_t m = B->size2;
    assert(n >= m);
 
    size_t Q_size = n * n;
    size_t R_size = n * m;
    ws->data = malloc((Q_size + R_size + 2 * m) * sizeof(double));
    llibcheck_mem(ws->data, error_a);

    ws->s = ws->data + Q_size + R_size;
    ws->y = ws->s + m;
    ws->v_s = gsl_vector_view_array(ws->s, m);
    ws->v_y = gsl_vector_view_array(ws->y, m);

    ws->m_R = gsl_matrix_view_array(ws->data + Q_size, n, m);
    ws->R = &ws->m_R.matrix;

    /* Tau is a vector used only for the gsl QR decomposition function.
     * The size of tau must be min of n,m; this is always m for now.*/
    size_t tausize = (n < m) ? n : m;
    /* Compute the QR-decomposition of lattice basis B and store it to matrices
     * Q and R, don't alter B. */
    gsl_matrix* B_copy = gsl_matrix_alloc(n, m);
    llibcheck_mem(B_copy, error_e);

    gsl_vector* tau = gsl_vector_alloc(tausize);
    llibcheck_mem(tau, error_f);

    gsl_matrix_view m_Q = gsl_matrix_view_array(ws->data, n, n);
    gsl_matrix* Q = &m_Q.matrix;

    gsl_matrix_memcpy(B_copy, B);
    gsl_linalg_QR_decomp(B_copy, tau);
    gsl_linalg_QR_unpack(B_copy, tau, Q, ws->R);
    gsl_vector_free(tau);
    gsl_matrix_free(B_copy);

    // Make the diagonal of R positive, as required by the algorithm.
    for(size_t i = 0; i < m; i++) 
    {
        if(gsl_matrix_get(ws->R, i, i) < 0) 
        {
            gsl_vector_view Rrow = gsl_matrix_row(ws->R, i);
            gsl_vector_view Qcol = gsl_matrix_column(Q,i);
            gsl_vector_scale(&Rrow.vector, -1);
            gsl_vector_scale(&Qcol.vector, -1);
        }
    }

    ws->Q1 = gsl_matrix_submatrix(Q, 0, 0, n, m);
    return ws;

error_f:
    gsl_matrix_free(B_copy);
error_e:
    free(ws->data);
error_a:
    free(ws);
error:
    return NULL;
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

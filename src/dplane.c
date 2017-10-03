/* dplane.c
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
#include "utility.h"
#include "lt_errno.h"
#include <math.h>
#include <assert.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "dplane.h"

struct s_dp_ws
{
    gsl_matrix* R;
    double* s;
    double* x;
    double* y;
    double* data;
    gsl_vector_view v_s;
    gsl_vector_view v_y;
    gsl_vector_view v_x;
    gsl_matrix_view m_R;
    gsl_matrix_view Q1;
    double dist_min;
    size_t m;
};

int DP_WS_alloc_and_init(DP_WS** ws_ptr, const gsl_matrix *B)
{
    int lt_errno = LT_SUCCESS;
    size_t n = B->size1;
    size_t m = B->size2;
    libcheck_se(n >= m, lt_errno, LT_EINVAL, 
            "the basis vectors are not linearly independent");

    DP_WS* ws = malloc(sizeof(DP_WS));
    libcheck_se_mem(ws, lt_errno, LT_ESYSTEM);

    ws->m = m;
    size_t Q_size = n * n;
    size_t R_size = n * m;
    ws->data = malloc((Q_size + R_size + 3 * m) * sizeof(double));
    llibcheck_se_mem(ws->data, error_a, lt_errno, LT_ESYSTEM);

    ws->s = ws->data + Q_size + R_size;
    ws->x = ws->s + m;
    ws->y = ws->x + m;

    ws->v_s = gsl_vector_view_array(ws->s, m);
    ws->v_x = gsl_vector_view_array(ws->x, m);
    ws->v_y = gsl_vector_view_array(ws->y, m);

    gsl_matrix_view m_Q = gsl_matrix_view_array(ws->data, n, n);
    gsl_matrix_view m_R = gsl_matrix_view_array(ws->data + Q_size, n, m);
    gsl_matrix* Q = &m_Q.matrix;
    gsl_matrix* R = &m_R.matrix;

    lt_errno = utility_compute_QR_decomposition(Q, R, B);
    lt_llibcheck(lt_errno, error_b, "utility_compute_QR_decomposition failed");
    int rc = utility_Rmm_is_not_singular(R, 10E-10);
    llibcheck_se(rc, error_b, lt_errno, LT_EINVAL, 
            "the basis vectors are not linearly independent");

    ws->Q1 = gsl_matrix_submatrix(Q, 0, 0, n, m);
    ws->m_R = gsl_matrix_submatrix(R, 0, 0, m, m);
    ws->R = &ws->m_R.matrix;

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

void DP_WS_free(DP_WS *ws) 
{
    if (ws) 
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

static void dplane_helper(DP_WS* ws, size_t m, double dist)
{
    if(m == 0)
    {
        double e = calc_yhat(m, ws->R, ws->y, ws->s) / gsl_matrix_get(ws->R, m, m);
        ws->s[m] = round(e);
        double d = (e - ws->s[m]) * gsl_matrix_get(ws->R, m, m);
        double new_dist = dist + d * d;
        if(new_dist < ws->dist_min)
        {
            memcpy(ws->x, ws->s, ws->m * sizeof(double));
            ws->dist_min = new_dist;
        }
    }
    else
    {
        double e = calc_yhat(m, ws->R, ws->y, ws->s) / gsl_matrix_get(ws->R, m, m);
        ws->s[m] = round(e);
        double d = (e - ws->s[m]) * gsl_matrix_get(ws->R, m, m);
        double new_dist = dist + d * d;
        if(new_dist < ws->dist_min)
            dplane_helper(ws, m-1, new_dist);

        if(d >= 0.0)
        {
            ws->s[m] += 1.0;
            d -= gsl_matrix_get(ws->R, m, m);
        }
        else
        {
            ws->s[m] -= 1.0;
            d += gsl_matrix_get(ws->R, m, m);
        }
        new_dist = dist + d * d;
        if(new_dist < ws->dist_min)
            dplane_helper(ws, m-1, new_dist);
    }
}

static void dplane(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, DP_WS *ws)
{
    gsl_blas_dgemv(CblasTrans, 1, &ws->Q1.matrix, t, 0, &ws->v_y.vector);
    ws->dist_min = INFINITY;
    dplane_helper(ws, B->size2-1, 0.0);
    gsl_blas_dgemv(CblasNoTrans, 1, B, &ws->v_x.vector, 0, clp);
}

static void dplane_vanilla_helper(DP_WS* ws, size_t m, double dist)
{
    if(m == 0)
    {
        double e = calc_yhat(m, ws->R, ws->y, ws->s) / gsl_matrix_get(ws->R, m, m);
        ws->s[m] = round(e);
        double d = (e - ws->s[m]) * gsl_matrix_get(ws->R, m, m);
        double new_dist = dist + d * d;
        if(new_dist < ws->dist_min)
        {
            memcpy(ws->x, ws->s, ws->m * sizeof(double));
            ws->dist_min = new_dist;
        }
    }
    else
    {
        // The floor version
        double e = calc_yhat(m, ws->R, ws->y, ws->s) / gsl_matrix_get(ws->R, m, m);
        ws->s[m] = floor(e);
        double d = (e - ws->s[m]) * gsl_matrix_get(ws->R, m, m);
        dplane_vanilla_helper(ws, m-1, dist + d * d);

        // The ceil version
        ws->s[m] += 1.0;
        d -= gsl_matrix_get(ws->R, m, m);
        dplane_vanilla_helper(ws, m-1, dist + d * d);
    }
}

static void dplane_vanilla(gsl_vector* clp, const gsl_vector* t,
                            const gsl_matrix* B, DP_WS *ws)
{
    gsl_blas_dgemv(CblasTrans, 1, &ws->Q1.matrix, t, 0, &ws->v_y.vector);
    ws->dist_min = INFINITY;
    dplane_vanilla_helper(ws, B->size2-1, 0.0);
    gsl_blas_dgemv(CblasNoTrans, 1, B, &ws->v_x.vector, 0, clp);
}

static void dplane_2d(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, DP_WS *ws)
{
    gsl_blas_dgemv(CblasTrans, 1, &ws->Q1.matrix, t, 0, &ws->v_y.vector);

    double e1 = ws->y[1] / gsl_matrix_get(ws->R, 1, 1);
    ws->s[1] = floor(e1);
    double d1 = (e1 - ws->s[1]) * gsl_matrix_get(ws->R, 1, 1);

    double yhat = ws->y[0] - gsl_matrix_get(ws->R, 0, 1) * ws->s[1];
    double e2 = yhat / gsl_matrix_get(ws->R, 0, 0);
    ws->s[0] = round(e2);
    double d2 = (e2 - ws->s[0]) * gsl_matrix_get(ws->R, 0, 0);
    double dist1 = d1 * d1 + d2 * d2;

    ws->x[1] = ws->s[1] + 1.0;
    d1 -= gsl_matrix_get(ws->R, 1, 1);
    yhat -= gsl_matrix_get(ws->R, 0, 1);
    e2 = yhat / gsl_matrix_get(ws->R, 0, 0);
    ws->x[0] = round(e2);
    d2 = (e2 - ws->x[0]) * gsl_matrix_get(ws->R, 0, 0);
    double dist2 = d1 * d1 + d2 * d2;

    if(dist2 < dist1)
        gsl_blas_dgemv(CblasNoTrans, 1, B, &ws->v_x.vector, 0, clp);
    else
        gsl_blas_dgemv(CblasNoTrans, 1, B, &ws->v_s.vector, 0, clp);
}

void dplane_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{
    switch(B->size2)
    {
        case 2:
            dplane_2d(clp, t, B, ws);
            break;
        default:
            dplane(clp, t, B, ws);
            break;
    }
}

void dplane_vanilla_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{
    switch(B->size2)
    {
        case 2:
            dplane_2d(clp, t, B, ws);
            break;
        default:
            dplane_vanilla(clp, t, B, ws);
            break;
    }
}

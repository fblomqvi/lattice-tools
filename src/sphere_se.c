/* sphere_se.h
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
#include "sphere_se.h"
#include "utility.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

struct s_sdse_ws
{
    gsl_matrix* Q;
    gsl_matrix* R;
    double* dist;
    double* e;
    double* s;
    double* x;
    double* y;
    double* step;
    double* data;
    gsl_vector_view v_y;
    gsl_vector_view v_x;
    gsl_matrix_view m_Q;
    gsl_matrix_view m_R;
    gsl_matrix_view Q1;
};

int SDSE_WS_alloc_and_init(SDSE_WS** ws_ptr, const gsl_matrix* B)
{
    int lt_errno = LT_SUCCESS;
    const size_t n = B->size1;
    const size_t m = B->size2;
    libcheck_se(n >= m, lt_errno, LT_EINVAL,
            "the basis vectors are not linearly independent");

    SDSE_WS* ws = malloc(sizeof(SDSE_WS));
    libcheck_se_mem(ws, lt_errno, LT_ESYSTEM);

    const size_t Q_size = n * n;
    const size_t R_size = n * m;
    ws->data = malloc((Q_size + R_size + 6 * m) * sizeof(double));
    llibcheck_se_mem(ws->data, error_a, lt_errno, LT_ESYSTEM);

    ws->dist = ws->data + Q_size + R_size;
    ws->e = ws->dist + m;
    ws->s = ws->e + m;
    ws->x = ws->s + m;
    ws->y = ws->x + m;
    ws->step = ws->y + m;

    ws->v_x = gsl_vector_view_array(ws->x, m);
    ws->v_y = gsl_vector_view_array(ws->y, m);

    ws->m_Q = gsl_matrix_view_array(ws->data, n, n);
    ws->m_R = gsl_matrix_view_array(ws->data + Q_size, n, m);
    ws->Q = &ws->m_Q.matrix;
    ws->R = &ws->m_R.matrix;

    lt_errno = utility_compute_QR_decomposition(ws->Q, ws->R, B);
    lt_llibcheck(lt_errno, error_b, "utility_compute_QR_decomposition failed");
    int rc = utility_Rmm_is_not_singular(ws->R, 10E-10);
    llibcheck_se(rc, error_b, lt_errno, LT_EINVAL,
            "the basis vectors are not linearly independent");

    ws->Q1 = gsl_matrix_submatrix(ws->Q, 0, 0, n, m);

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

void SDSE_WS_free(SDSE_WS *ws)
{
    if(ws)
    {
        free(ws->data);
        free(ws);
    }
}

static double calc_yhat(size_t k, const gsl_matrix* R, const double* y, const double* s)
{
    size_t m = R->size2;
    double sum = y[k];
    for (size_t j = k+1; j < m; j++)
        sum -= gsl_matrix_get(R, k, j) * s[j];

    return sum;
}

static double sgn(double y)
{ return y <= 0.0 ? -1 : 1; }

static void move_down_calculations(size_t k, double* y, SDSE_WS* ws)
{
    ws->e[k] = calc_yhat(k, ws->R, ws->y, ws->s) / gsl_matrix_get(ws->R, k, k);
    ws->s[k] = round(ws->e[k]);
    *y = (ws->e[k] - ws->s[k]) * gsl_matrix_get(ws->R, k, k);
    ws->step[k] = sgn(*y);
}

static void move_up_calculations(size_t k, double* y, SDSE_WS* ws)
{
    ws->s[k] += ws->step[k];
    *y = (ws->e[k] - ws->s[k]) * gsl_matrix_get(ws->R, k, k);
    ws->step[k] = -ws->step[k] - sgn(ws->step[k]);
}

void sphere_se(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, SDSE_WS* ws)
{   
    gsl_blas_dgemv(CblasTrans, 1, &ws->Q1.matrix, t, 0, &ws->v_y.vector);

    size_t m = ws->R->size2;
    double dist_min = INFINITY;
    size_t k = m - 1;
    ws->dist[k] = 0;

    double y;
    move_down_calculations(k, &y, ws);

    while(1)
    {
        double new_dist = ws->dist[k] + y * y;
        if(new_dist < dist_min)
        {
            if(k > 0)
            {
                k--;        // Move down
                ws->dist[k] = new_dist;
                move_down_calculations(k, &y, ws);
            }
            else
            {
                memcpy(ws->x, ws->s, m * sizeof(double));
                dist_min = new_dist;
                k++;        // Move up
                move_up_calculations(k, &y, ws);
            }
        }
        else
        {
            if(k == m - 1)
                break;
            else
            {
                k++;        // Move up
                move_up_calculations(k, &y, ws);
            }
        }
    }

    gsl_blas_dgemv(CblasNoTrans, 1, B, &ws->v_x.vector, 0, clp);
}

void sphere_se_svp(gsl_vector* clp, const gsl_matrix* B, SDSE_WS* ws)
{
    size_t m = ws->R->size2;
    double dist_min = INFINITY;
    size_t k = m - 1;
    ws->dist[k] = 0;
    memset(ws->y, 0, m * sizeof(double));

    ws->e[k] = ws->s[k] = 0.0;
    ws->s[k] = round(ws->e[k]);
    double y = 0.0;
    ws->step[k] = sgn(y);

    while(1)
    {
        double new_dist = ws->dist[k] + y * y;
        if(new_dist < dist_min)
        {
            if(k > 0)
            {
                k--;        // Move down
                ws->dist[k] = new_dist;
                move_down_calculations(k, &y, ws);
            }
            else
            {
                if(new_dist != 0.0)
                {
                    memcpy(ws->x, ws->s, m * sizeof(double));
                    dist_min = new_dist;
                    k++;        // Move up
                }
                move_up_calculations(k, &y, ws);
            }
        }
        else
        {
            if(k == m - 1)
                break;
            else
            {
                k++;        // Move up
                move_up_calculations(k, &y, ws);
            }
        }
    }

    gsl_blas_dgemv(CblasNoTrans, 1, B, &ws->v_x.vector, 0, clp);
}

void sphere_se_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{ sphere_se(clp, t, B, ws); }

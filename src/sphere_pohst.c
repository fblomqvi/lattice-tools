/* sphere_pohst.c
   Copyright (C) 2017 Ferdinand Blomqvist, Juuso Korvuo

   This program is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 as published by
   the Free Software Foundation. 

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along with
   this program. If not, see <http://www.gnu.org/licenses/>.  
   
   Written by Ferdinand Blomqvist and Juuso Korvuo. */

#include "dbg.h"
#include "sphere_pohst.h"
#include "utility.h"
#include "lt_errno.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

struct s_sd_ws
{
    gsl_matrix* Q;
    gsl_matrix* R;
    double* ub;
    double* d2;
    double* yhat;
    double* s;
    double* x;
    double* y;
    double* data;
    gsl_vector_view v_s;
    gsl_vector_view v_y;
    gsl_vector_view v_x;
    gsl_matrix_view m_Q;
    gsl_matrix_view m_R;
    gsl_matrix_view Q1;
    gsl_matrix_view Rsub;
};

static void sd_dp_common(SD_WS* ws, const gsl_vector* t, double* d_sqr, size_t m);

static void set_bounds_dp(SD_WS* ws, size_t k);
static void set_bounds_sd(SD_WS* ws, size_t k);

static void step_5a(SD_WS* ws, double* d_sqr, size_t m);
static void step_5b(SD_WS* ws, size_t* k);

/*NOTE: n is used to denote the amount of rows of the basis matrix B and m for
 *      columns. This notation was chosen in order to be consistent with
 *      the documentation and it is the OPPOSITE of how gsl uses them
 *      (in error messages, for example). */

int SD_WS_alloc_and_init(SD_WS** ws_ptr, const gsl_matrix* B)
{
    int lt_errno = LT_SUCCESS;
    size_t n = B->size1;
    size_t m = B->size2;
    llibcheck_se(n >= m, error, lt_errno, LT_ELINDEP);

    SD_WS* ws = malloc(sizeof(SD_WS));
    llibcheck_se(ws, error, lt_errno, LT_ENOMEM);

    size_t Q_size = n * n;
    size_t R_size = n * m;
    ws->data = malloc((Q_size + R_size + 6 * m) * sizeof(double));
    llibcheck_se(ws->data, error_a, lt_errno, LT_ENOMEM);

    ws->ub = ws->data + Q_size + R_size;
    ws->d2 = ws->ub + m;
    ws->yhat = ws->d2 + m;
    ws->s = ws->yhat + m;
    ws->x = ws->s + m;
    ws->y = ws->x + m;

    ws->v_s = gsl_vector_view_array(ws->s, m);
    ws->v_x = gsl_vector_view_array(ws->x, m);
    ws->v_y = gsl_vector_view_array(ws->y, m);

    ws->m_Q = gsl_matrix_view_array(ws->data, n, n);
    ws->m_R = gsl_matrix_view_array(ws->data + Q_size, n, m);
    ws->Q = &ws->m_Q.matrix;
    ws->R = &ws->m_R.matrix;

    lt_errno = utility_compute_QR_decomposition(ws->Q, ws->R, B);
    llibcheck(lt_errno == 0, error_b, "utility_compute_QR_decomposition failed");
    llibcheck_se(utility_Rmm_is_not_singular(ws->R, 10E-10), error_b, lt_errno, LT_ELINDEP);

    ws->Q1 = gsl_matrix_submatrix(ws->Q, 0, 0, n, m);
    ws->Rsub = gsl_matrix_submatrix(ws->R, 0, 0, m, m);
    
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

void SD_WS_free(SD_WS *ws)
{
    if(ws)
    {
        free(ws->data);
        free(ws);
    }
}

/* Calculate the end return the value of ŷ_k. ŷ_k is y_k from which have been
 * subtracted the kth row values of R multiplied by corresponding values of s 
 * excluding the diagonal. Example with k=2:
 *   [d0  #  #  #  #]      [s0]
 *   [ 0 d1  #  #  #]      [s1]
 * R=[ 0  0 d2  a  b]    s=[s2]  =>  ŷ_2 = y_2 - a*s3 - b*s4
 *   [ 0  0  0 d3  #]      [s3]
 *   [ 0  0  0  0 d4]      [s4]
 */

static double calc_yhat(size_t k, const gsl_matrix* R, const double* y, const double* s)
{
    size_t m = R->size2;
    double sum = y[k];
    for (size_t j = k+1; j < m; j++)
        sum -= gsl_matrix_get(R, k, j) * s[j];

    return sum;
}

/* Calculate and return the value of d_k^2. s_k has to satisfy
 * d_k^2 >= (ŷ_k - r_(k,k)*s_k)^2 for the current point to be inside the
 * hypersphere and d_k^2 is calculated as
 * d_k^2 = d_(k+1)^2 - (ŷ_(k+1) - r(k+1,k+1)*s_(k+1))^2
 */
static double calc_d2(size_t k, SD_WS* ws)
{
    size_t idx = k + 1;
    double res = ws->yhat[idx] - gsl_matrix_get(ws->R, idx, idx) * ws->s[idx];
    return ws->d2[idx] - res * res;
}

/* Finds the closest vector in the lattice with basis B to the target vector t
 * by finding an integer vector s so that ||t-B*s|| is minimized. The algorithm
 * looks for lattice points inside a hypersphere of radius d and picks the values
 * of s which minimize the distance.
 * The best vector s is found by calculating upper and lower bounds for kth
 * element of s (denoted by s_k) and systematically incrementing s_k starting
 * from the lower bound. Whenever B*s is inside the hypersphere, move to k-1:th
 * element of s and repeat the process until k=0 (all coordinates of s have 
 * been determined).
 * If the value of s_k causes the current point not to be inside the sphere,
 * move to previous coordinate s_(k+1) and repeat the process there.
 * The algorithm terminates when s_(m-1) (the last value of s) becomes greater than
 * its upper bound; all possible solutions have been found.
 */
void spheredecode(gsl_vector* clp, const gsl_vector* t, const gsl_matrix *B, SD_WS *ws)
{   
    assert(B->size1 == t->size);
    size_t m = B->size2;
    double d_sqr;
    sd_dp_common(ws, t, &d_sqr, m);

    /*******************************************************************/
    /* ACTUAL ALGORITHM STARTS HERE                                    */
    /*******************************************************************/
    
    // Step 1: set initial values.
    int set_new_bounds = 1;
    size_t solutions = 0;
    size_t k = m - 1;
    ws->d2[k] = d_sqr * 1.000000001;
    ws->yhat[k] = calc_yhat(k, ws->R, ws->y, ws->s);

    while(k < m)
    {
        if(set_new_bounds)
        {
            set_bounds_sd(ws, k);
            set_new_bounds = 0;
        }
        else
            ws->s[k] += 1.0;

        if(ws->s[k] > ws->ub[k])
            k++;
        else
        {
            if(k == 0)
            {
                solutions++;
                step_5a(ws, &d_sqr, m);
            }
            else 
            {
                step_5b(ws, &k);
                set_new_bounds = 1;
            }
        }
    }

    gsl_blas_dgemv(CblasNoTrans, 1, B, &ws->v_x.vector, 0, clp);
    //fprintf(stderr, "SD: Tried %d points\n", solutions);
}

void sd_dp(gsl_vector* clp, const gsl_vector* t, const gsl_matrix *B, SD_WS *ws)
{   
    assert(B->size1 == t->size);
    size_t m = B->size2;
    double d_sqr;
    sd_dp_common(ws, t, &d_sqr, m);

    // Step 1: set initial values.
    int set_new_bounds = 1;
    size_t solutions = 0;
    size_t k = m - 1;
    ws->d2[k] = d_sqr * 1.000000001;
    ws->yhat[k] = calc_yhat(k, ws->R, ws->y, ws->s);

    while(k < m)
    {
        if(set_new_bounds)
        {
            set_bounds_dp(ws, k);
            set_new_bounds = 0;
        }
        else
            ws->s[k] += 1.0;

        if(ws->s[k] > ws->ub[k])
            k++;
        else
        {
            if(k == 0)
            {
                solutions++;
                step_5a(ws, &d_sqr, m);
            }
            else 
            {
                step_5b(ws, &k);
                set_new_bounds = 1;
            }
        }
    }

    gsl_blas_dgemv(CblasNoTrans, 1, B, &ws->v_x.vector, 0, clp);
    //fprintf(stderr, "SD: Tried %d points\n", solutions);
}
void spheredecode_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{ spheredecode(clp, t, B, ws); }

void sd_dp_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{ sd_dp(clp, t, B, ws); }

static void sd_dp_common(SD_WS* ws, const gsl_vector* t, double* d_sqr, size_t m)
{
    // First we do Babai and calculate the distance to the Babai point
    gsl_blas_dgemv(CblasTrans, 1, &ws->Q1.matrix, t, 0, &ws->v_y.vector);
    double dist = 0.0;
    for(int k = m-1; k >= 0; k--)
    {
        double e = calc_yhat(k, ws->R, ws->y, ws->x) / gsl_matrix_get(ws->R, k, k);
        ws->x[k] = round(e);
        double d = (e - ws->x[k]) * gsl_matrix_get(ws->R, k, k);
        dist += d * d;
    }
    *d_sqr = dist;
    
    /*
    // Set Q2Tx_len2 = ||(Q2^T)*x||^2.
    double Q2Tx_len2 = 0;
    if (n > m) {
        gsl_matrix_view Q2 = gsl_matrix_submatrix(ws->Q, 0, m, n, n-m);
        gsl_vector *Q2Tx = gsl_vector_alloc(n-m);
        gsl_blas_dgemv(CblasTrans, 1, &Q2.matrix, t, 0, Q2Tx);
        gsl_blas_ddot(Q2Tx, Q2Tx, &Q2Tx_len2);
        gsl_vector_free(Q2Tx);
    }
    */
}

static inline void set_bounds_sd(SD_WS* ws, size_t k)
{
    double Rkk = gsl_matrix_get(ws->R, k, k);
    double d2_sqrt_div_rkk = sqrt(ws->d2[k]) / Rkk;
    double yhat_div_rkk = ws->yhat[k] / Rkk;
    ws->ub[k] = floor(d2_sqrt_div_rkk + yhat_div_rkk);
    ws->s[k] = ceil(-d2_sqrt_div_rkk + yhat_div_rkk);
}

static inline void set_bounds_dp(SD_WS* ws, size_t k)
{
    double Rkk = gsl_matrix_get(ws->R, k, k);
    double d2_sqrt_div_rkk = sqrt(ws->d2[k]) / Rkk;
    double yhat_div_rkk = ws->yhat[k] / Rkk;
    double max1 = floor(d2_sqrt_div_rkk + yhat_div_rkk);
    double min1 = ceil(-d2_sqrt_div_rkk + yhat_div_rkk);
    double min2 = floor(yhat_div_rkk);
    double max2 = min2 + 1;
    ws->s[k] = min1 < min2 ? min2 : min1;
    ws->ub[k] = max1 > max2 ? max2 : max1;
}

static inline void step_5a(SD_WS* ws, double* d_sqr, size_t m)
{
    // Step 5a: Solution found.
    double res = ws->yhat[0] - gsl_matrix_get(ws->R, 0, 0) * ws->s[0];
    double y_minus_Rs_len = ws->d2[m-1] - ws->d2[0] + res * res;

    // If new s gives a better solution to the CVP, store it to x.
    if(y_minus_Rs_len < *d_sqr)
    {
        //gsl_vector_memcpy(&ws->v_x.vector, &ws->v_s.vector);
        memcpy(ws->x, ws->s, m * sizeof(double));
        *d_sqr = y_minus_Rs_len;
    }
}

static inline void step_5b(SD_WS* ws, size_t* k)
{
    // Step 5b: Decrease k and calculate new ŷ_k and d_k^2
    (*k)--;
    ws->yhat[*k] = calc_yhat(*k, ws->R, ws->y, ws->s);
    ws->d2[*k] = calc_d2(*k, ws);
}

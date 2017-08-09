#include "dbg.h"
#include <math.h>
#include <assert.h>
#include <gsl/gsl_blas.h>
#include "doubleplane.h"

struct s_dp_ws
{
    gsl_matrix* R;
    double* s;
    double* x;
    double* y;
    double* Rs;
    double* data;
    gsl_vector_view v_s;
    gsl_vector_view v_y;
    gsl_vector_view v_x;
    gsl_vector_view v_Rs;
    gsl_matrix_view m_R;
    gsl_matrix_view Q1;
    double dist_min;
    size_t m;
};

DP_WS* DP_WS_alloc_and_init(const gsl_matrix *B)
{
    DP_WS* ws = malloc(sizeof(DP_WS));
    libcheck_mem(ws);

    size_t n = B->size1;
    size_t m = B->size2;
    assert(n >= m);

    ws->m = m;
    size_t Q_size = n * n;
    size_t R_size = n * m;
    ws->data = malloc((Q_size + R_size + 4 * m) * sizeof(double));
    llibcheck_mem(ws->data, error_a);

    ws->s = ws->data + Q_size + R_size;
    ws->x = ws->s + m;
    ws->y = ws->x + m;
    ws->Rs = ws->y + m;

    ws->v_s = gsl_vector_view_array(ws->s, m);
    ws->v_x = gsl_vector_view_array(ws->x, m);
    ws->v_y = gsl_vector_view_array(ws->y, m);
    ws->v_Rs = gsl_vector_view_array(ws->Rs, m);

    gsl_matrix_view m_Q = gsl_matrix_view_array(ws->data, n, n);
    gsl_matrix_view m_R = gsl_matrix_view_array(ws->data + Q_size, n, m);
    gsl_matrix* Q = &m_Q.matrix;
    gsl_matrix* R = &m_R.matrix;

    ws->Q1 = gsl_matrix_submatrix(Q, 0, 0, n, m);
    ws->m_R = gsl_matrix_submatrix(R, 0, 0, m, m);
    ws->R = &ws->m_R.matrix;

    /* Tau is a vector used only for the gsl QR decomposition function.
     * The size of tau must be min of n,m; this is always m for now.*/
    size_t tausize = (n < m) ? n : m;
    /* Compute the QR-decomposition of lattice basis B and store it to matrices
     * Q and R, don't alter B. */
    gsl_matrix* B_copy = gsl_matrix_alloc(n, m);
    llibcheck_mem(B_copy, error_d);

    gsl_vector* tau = gsl_vector_alloc(tausize);
    llibcheck_mem(tau, error_e);

    gsl_matrix_memcpy(B_copy, B);
    gsl_linalg_QR_decomp(B_copy, tau);
    gsl_linalg_QR_unpack(B_copy, tau, Q, R);
    gsl_vector_free(tau);
    gsl_matrix_free(B_copy);

    // Make the diagonal of R positive, as required by the algorithm.
    for(size_t i = 0; i < m; i++)
    {
        if(gsl_matrix_get(R, i, i) < 0)
        {
            gsl_vector_view Rrow = gsl_matrix_row(R, i);
            gsl_vector_view Qcol = gsl_matrix_column(Q,i);
            gsl_vector_scale(&Rrow.vector, -1);
            gsl_vector_scale(&Qcol.vector, -1);
        }
    }

    return ws;

error_e:
    gsl_matrix_free(B_copy);
error_d:
    free(ws->data);
error_a:
    free(ws);
error:
    return NULL;
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

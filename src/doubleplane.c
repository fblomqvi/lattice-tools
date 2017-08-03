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
    int solutions;
    double dist;
    double dist_min;
};

DP_WS* DP_WS_alloc_and_init(const gsl_matrix *B)
{
    DP_WS* ws = malloc(sizeof(DP_WS));
    libcheck_mem(ws);

    size_t n = B->size1;
    size_t m = B->size2;
    assert(n >= m);

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
    double sum = 0;
    for (size_t j = k+1; j < m; j++)
        sum += gsl_matrix_get(R, k, j) * s[j];
    
    return y[k] - sum;
}
void doubleplane_helper(DP_WS* ws, size_t m) 
{
    if(m == 0)
    {
        ws->s[m] = round(calc_yhat(m, ws->R, ws->y, ws->s) / gsl_matrix_get(ws->R, m, m));

        if(ws->solutions)
        {
            // Caclulate ||y-Rs|| and replace ws->x with the new value if needed.
            gsl_vector_memcpy(&ws->v_Rs.vector, &ws->v_s.vector);
            gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, 
                            ws->R, &ws->v_Rs.vector);
            gsl_vector_sub(&ws->v_Rs.vector, &ws->v_y.vector);
            gsl_blas_ddot(&ws->v_Rs.vector, &ws->v_Rs.vector, &ws->dist);
            if(ws->dist < ws->dist_min)
            {
                ws->dist_min = ws->dist;
                gsl_vector_memcpy(&ws->v_x.vector, &ws->v_s.vector);
            }
        }
        else
        {
            gsl_vector_memcpy(&ws->v_Rs.vector, &ws->v_s.vector);
            gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, 
                            ws->R, &ws->v_Rs.vector);
            gsl_vector_sub(&ws->v_Rs.vector, &ws->v_y.vector);
            gsl_blas_ddot(&ws->v_Rs.vector, &ws->v_Rs.vector, &ws->dist_min);
            gsl_vector_memcpy(&ws->v_x.vector, &ws->v_s.vector);
            ws->solutions = 1;
        }
    }
    else
    {
        // The floor version
        ws->s[m] = floor(calc_yhat(m, ws->R, ws->y, ws->s) / gsl_matrix_get(ws->R, m, m));
        doubleplane_helper(ws, m-1);

        // The ceil version
        ws->s[m] += 1.0;
        doubleplane_helper(ws, m-1);
    }
}

void doubleplane(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, DP_WS *ws)
{
    gsl_blas_dgemv(CblasTrans, 1, &ws->Q1.matrix, t, 0, &ws->v_y.vector);
    ws->solutions = 0;
    doubleplane_helper(ws, B->size2-1);
    gsl_blas_dgemv(CblasNoTrans, 1, B, &ws->v_x.vector, 0, clp);
}

static inline void Rx_sub_y(double* Rs, const gsl_matrix* R,
                            const double* x, const double* y)
{
    Rs[0] = x[0] * gsl_matrix_get(R, 0, 0)
                + x[1] * gsl_matrix_get(R, 0, 1) - y[0];
    Rs[1] = x[1] * gsl_matrix_get(R, 1, 1) - y[1];
}

static inline double dist_sqr_2d(const double* v)
{ return v[0] * v[0] + v[1] * v[1]; }

void doubleplane_2d(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, DP_WS *ws)
{
    gsl_blas_dgemv(CblasTrans, 1, &ws->Q1.matrix, t, 0, &ws->v_y.vector);

    ws->s[1] = floor(ws->y[1] / gsl_matrix_get(ws->R, 1, 1));
    double yhat = ws->y[0] - gsl_matrix_get(ws->R, 0, 1) * ws->s[1];
    ws->s[0] = round(yhat / gsl_matrix_get(ws->R, 0, 0));

    ws->x[1] = ws->s[1] + 1.0;
    yhat -= gsl_matrix_get(ws->R, 0, 1);
    ws->x[0] = round(yhat / gsl_matrix_get(ws->R, 0, 0));

    Rx_sub_y(ws->Rs, ws->R, ws->s, ws->y);
    double dist1 = dist_sqr_2d(ws->Rs);

    Rx_sub_y(ws->Rs, ws->R, ws->x, ws->y);
    double dist2 = dist_sqr_2d(ws->Rs);

    if(dist2 < dist1)
        gsl_blas_dgemv(CblasNoTrans, 1, B, &ws->v_x.vector, 0, clp);
    else
        gsl_blas_dgemv(CblasNoTrans, 1, B, &ws->v_s.vector, 0, clp);
}

void doubleplane_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{
    switch(B->size2)
    {
        case 2:
            doubleplane_2d(clp, t, B, ws);
            break;
        default:
            doubleplane(clp, t, B, ws);
            break;
    }
}

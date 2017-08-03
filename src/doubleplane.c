#include "dbg.h"
#include <math.h>
#include <assert.h>
#include <gsl/gsl_blas.h>
#include "doubleplane.h"

struct s_dp_ws
{
    gsl_matrix* Q;
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
    gsl_matrix_view m_Q;
    gsl_matrix_view m_R;
    gsl_matrix_view Q1;
    gsl_matrix_view Rsub;
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

    ws->m_Q = gsl_matrix_view_array(ws->data, n, n);
    ws->m_R = gsl_matrix_view_array(ws->data + Q_size, n, m);
    ws->Q = &ws->m_Q.matrix;
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
    gsl_linalg_QR_unpack(B_copy, tau, ws->Q, ws->R);
    gsl_vector_free(tau);
    gsl_matrix_free(B_copy);

    // Make the diagonal of R positive, as required by the algorithm.
    for(size_t i = 0; i < m; i++)
    {
        if(gsl_matrix_get(ws->R, i, i) < 0)
        {
            gsl_vector_view Rrow = gsl_matrix_row(ws->R, i);
            gsl_vector_view Qcol = gsl_matrix_column(ws->Q,i);
            gsl_vector_scale(&Rrow.vector, -1);
            gsl_vector_scale(&Qcol.vector, -1);
        }
    }

    ws->Q1 = gsl_matrix_submatrix(ws->Q, 0, 0, n, m);
    ws->Rsub = gsl_matrix_submatrix(ws->R, 0, 0, m, m);

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
                            &ws->Rsub.matrix, &ws->v_Rs.vector);
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
                            &ws->Rsub.matrix, &ws->v_Rs.vector);
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

void doubleplane_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{ doubleplane(clp, t, B, ws); }


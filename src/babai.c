#include "dbg.h"
#include "babai.h"
#include <math.h>
#include <assert.h>
#include "utility.h"
#include <gsl/gsl_blas.h>

struct s_babai_ws
{
    gsl_matrix* Q;
    gsl_matrix* R;
    gsl_vector* s;
    gsl_vector* y;
};

BABAI_WS* BABAI_WS_alloc_and_init(const gsl_matrix* B)
{
    BABAI_WS* ws = malloc(sizeof(BABAI_WS));
    libcheck_mem(ws);

    size_t n = B->size1;
    size_t m = B->size2;
    assert(n >= m);
 
    ws->Q = gsl_matrix_alloc(n,n);
    llibcheck_mem(ws->Q, error_a);

    ws->R = gsl_matrix_alloc(n,m);    
    llibcheck_mem(ws->R, error_b);

    ws->y = gsl_vector_alloc(m);
    llibcheck_mem(ws->y, error_c);

    ws->s = gsl_vector_alloc(m);
    llibcheck_mem(ws->s, error_d);

    /* Tau is a vector used only for the gsl QR decomposition function.
     * The size of tau must be min of n,m; this is always m for now.*/
    size_t tausize = (n < m) ? n : m;
    /* Compute the QR-decomposition of lattice basis B and store it to matrices
     * Q and R, don't alter B. */
    gsl_matrix *B_copy = gsl_matrix_alloc(n, m);
    gsl_matrix_memcpy(B_copy, B);
    gsl_vector *tau = gsl_vector_alloc(tausize);
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
    return ws;

error_d:
    gsl_vector_free(ws->y);
error_c:
    gsl_matrix_free(ws->R);
error_b:
    gsl_matrix_free(ws->Q);
error_a:
    free(ws);
error:
    return NULL;
}

void BABAI_WS_free(BABAI_WS* ws)
{
    if(ws)
    {
        gsl_matrix_free(ws->Q);
        gsl_matrix_free(ws->R);
        gsl_vector_free(ws->y);
        gsl_vector_free(ws->s);
        free(ws);
    }
}

static double calc_yhat(size_t k, const gsl_matrix* R, 
                        const gsl_vector* y, const gsl_vector* s) 
{
    size_t m = R->size2;
    double sum = 0;
    for (size_t j = k+1; j < m; j++)
        sum += gsl_matrix_get(R, k, j) * gsl_vector_get(s, j);
    
    return gsl_vector_get(y, k) - sum;
}

void babai(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, BABAI_WS* ws)
{
    gsl_matrix_view Q1 = gsl_matrix_submatrix(ws->Q, 0, 0, B->size1, B->size2);
    gsl_blas_dgemv(CblasTrans, 1, &Q1.matrix, t, 0, ws->y);
    
    for(int k = B->size2-1; k >= 0; k--)
    {
        double s_k = calc_yhat(k, ws->R, ws->y, ws->s) / gsl_matrix_get(ws->R, k, k);
        gsl_vector_set(ws->s, k, round(s_k));
    }

    gsl_blas_dgemv(CblasNoTrans, 1, B, ws->s, 0, clp);
}

void babai_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{ babai(clp, t, B, ws); }

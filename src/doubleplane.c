#include "dbg.h"
#include <math.h>
#include <assert.h>
#include <gsl/gsl_blas.h>
#include "utility.h"
#include "doubleplane.h"
#include "babai.h"

struct s_dp_ws
{
    gsl_matrix* Q;
    gsl_matrix* R;
    gsl_vector* s;
    gsl_vector* y;
    gsl_vector* best_s;
    gsl_vector* Rs;
    gsl_vector* Rbest_s;
    double* ub;
    double* d2;
    double* yhat;
    gsl_matrix_view Q1;
    gsl_matrix_view Rsub;
};

DP_WS* DP_WS_alloc_and_init(const gsl_matrix *B)
{
    DP_WS* ws = malloc(sizeof(DP_WS));
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
    gsl_matrix* B_copy = gsl_matrix_alloc(n, m);
    llibcheck_mem(B_copy, error_e);

    gsl_vector* tau = gsl_vector_alloc(tausize);
    llibcheck_mem(tau, error_f);

    gsl_matrix_memcpy(B_copy, B);
    gsl_linalg_QR_decomp(B_copy, tau);
    gsl_linalg_QR_unpack(B_copy, tau, ws->Q, ws->R);
    gsl_vector_free(tau);
    gsl_matrix_free(B_copy);

    ws->best_s = gsl_vector_alloc(m);
    ws->Rs = gsl_vector_alloc(m);
    ws->Rbest_s = gsl_vector_alloc(m);
    ws->ub = malloc(3 * m * sizeof(double));
    ws->d2 = ws->ub + m;
    ws->yhat = ws->d2 + m;

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

error_f:
    gsl_matrix_free(B_copy);
error_e:
    gsl_vector_free(ws->s);
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

void DP_WS_free(DP_WS *ws) 
{
    if(ws)
    {
        gsl_matrix_free(ws->Q);
        gsl_matrix_free(ws->R);
        gsl_vector_free(ws->y);
        gsl_vector_free(ws->s);
        gsl_vector_free(ws->best_s);
        gsl_vector_free(ws->Rs);
        gsl_vector_free(ws->Rbest_s);
        free(ws->ub);
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

static double calc_d2(size_t k, DP_WS* ws)
{
    size_t idx = k + 1;
    double res = ws->yhat[idx] - gsl_matrix_get(ws->R, idx, idx) * gsl_vector_get(ws->s, idx);
    return -1 * res * res + ws->d2[idx];
}

void doubleplane(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, DP_WS *ws)
{
    assert(B->size1 == t->size);

    // First we do Babai
    gsl_blas_dgemv(CblasTrans, 1, &ws->Q1.matrix, t, 0, ws->y);
    
    for(int k = B->size2-1; k >= 0; k--)
    {
        double s_k = calc_yhat(k, ws->R, ws->y, ws->best_s) / gsl_matrix_get(ws->R, k, k);
        gsl_vector_set(ws->best_s, k, round(s_k));
    }

    // Calculate ||y-Rs||^2
    double d_sqr;
    gsl_blas_dgemv(CblasNoTrans, 1, &ws->Rsub.matrix, ws->best_s, 0, ws->Rbest_s);
    gsl_vector_sub(ws->Rbest_s, ws->y);
    gsl_blas_ddot(ws->Rbest_s, ws->Rbest_s, &d_sqr);

    // Now we can do the combination of doubleplane and spheredecoding
    //size_t n = B->size1;
    size_t m = B->size2;
    double y_minus_Rbest_s_len = d_sqr;
    
    // Step 1: set initial values.
    size_t k = m - 1;
    ws->d2[k] = d_sqr;
    ws->yhat[k] = calc_yhat(k, ws->R, ws->y, ws->s);
    
    int set_new_bounds = 1;
    size_t solutions = 0;
    
    while(k < m) 
    {
        if(set_new_bounds) 
        {
            // Step 2: Calculate upper bound for s_k and set s_k to lower bound.
            double Rkk = gsl_matrix_get(ws->R, k, k);
            double d2_sqrt_div_rkk = sqrt(ws->d2[k]) / Rkk;
            double yhat_div_rkk = ws->yhat[k] / Rkk;
            double max1 = floor(d2_sqrt_div_rkk + yhat_div_rkk);
            double min1 = ceil(-d2_sqrt_div_rkk + yhat_div_rkk);
            double min2 = floor(yhat_div_rkk);
            double max2 = min2 + 1;
            gsl_vector_set(ws->s, k, min1 < min2 ? min2 : min1);
            ws->ub[k] = max1 > max2 ? max2 : max1;
            // Don't set new bounds on next iteration unless required by step 5b.
            set_new_bounds = 0;
        }
        else
        {
            // Step 3: Set s_k = s_k + 1.
            gsl_vector_set(ws->s, k, gsl_vector_get(ws->s, k) + 1);
        }

        if(gsl_vector_get(ws->s, k) > ws->ub[k]) 
        {
            // Step 4: Increase k by 1. If k = m + 1, the algorithm will
            //         terminate (by the while loop).
            k++;
        } 
        else 
        {
            // Step 5
            if(k == 0) 
            {
                // Step 5a: Solution found.
                solutions++;

                gsl_blas_dgemv(CblasNoTrans, 1, &ws->Rsub.matrix, ws->s, 0, ws->Rs);
                gsl_vector_sub(ws->Rs, ws->y);

                // Rs = R*s-y; Rbest_s = R*best_s-y
                double y_minus_Rs_len;
                gsl_blas_ddot(ws->Rs, ws->Rs, &y_minus_Rs_len);

                // If new s gives a better solution to the CVP, store it to best_s.
                if(y_minus_Rs_len < y_minus_Rbest_s_len) 
                {
                    gsl_vector_memcpy(ws->best_s, ws->s);
                    gsl_vector_memcpy(ws->Rbest_s, ws->Rs);
                    y_minus_Rbest_s_len = y_minus_Rs_len;
                }
            } 
            else 
            {
                // Step 5b: Decrease k and calculate new ŷ_k and d_k^2
                k--;
                ws->yhat[k] = calc_yhat(k, ws->R, ws->y, ws->s);
                ws->d2[k] = calc_d2(k, ws);
                // Calculate new bounds for s_k.
                set_new_bounds = 1;
            }
        }
    }

    gsl_blas_dgemv(CblasNoTrans, 1, B, ws->best_s, 0, clp);
    //fprintf(stderr, "DP: Tried %zu points\n", solutions);
}

void doubleplane_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{ doubleplane(clp, t, B, ws); }


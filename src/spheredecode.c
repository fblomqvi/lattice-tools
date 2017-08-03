#include "dbg.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "spheredecode.h"
#include <gsl/gsl_blas.h>

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
    gsl_vector_view Rs;
    gsl_matrix_view m_Q;
    gsl_matrix_view m_R;
    gsl_matrix_view Q1;
    gsl_matrix_view Rsub;
};

/*NOTE: n is used to denote the amount of rows of the basis matrix B and m for
 *      columns. This notation was chosen in order to be consistent with
 *      the documentation and it is the OPPOSITE of how gsl uses them
 *      (in error messages, for example). */

SD_WS* SD_WS_alloc_and_init(const gsl_matrix* B)
{
    SD_WS* ws = malloc(sizeof(SD_WS));
    libcheck_mem(ws);

    size_t n = B->size1;
    size_t m = B->size2;
    assert(n >= m);

    size_t Q_size = n * n;
    size_t R_size = n * m;
    ws->data = malloc((Q_size + R_size + 7 * m) * sizeof(double));
    llibcheck_mem(ws->data, error_a);

    ws->ub = ws->data + Q_size + R_size;
    ws->d2 = ws->ub + m;
    ws->yhat = ws->d2 + m;
    ws->s = ws->yhat + m;
    ws->x = ws->s + m;
    ws->y = ws->x + m;

    ws->v_s = gsl_vector_view_array(ws->s, m);
    ws->v_x = gsl_vector_view_array(ws->x, m);
    ws->v_y = gsl_vector_view_array(ws->y, m);
    ws->Rs = gsl_vector_view_array(ws->y + m, m);

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
    double sum = 0;
    for (size_t j = k+1; j < m; j++)
        sum += gsl_matrix_get(R, k, j) * s[j];

    return y[k] - sum;
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
    return -1 * res * res + ws->d2[idx];
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

    // First we do Babai
    gsl_blas_dgemv(CblasTrans, 1, &ws->Q1.matrix, t, 0, &ws->v_y.vector);
    for(int k = m-1; k >= 0; k--)
        ws->x[k] = round(calc_yhat(k, ws->R, ws->y, ws->x) / gsl_matrix_get(ws->R, k, k));
    
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
        
    // Calculate ||y-Rx||^2
    double d_sqr;
    gsl_blas_dgemv(CblasNoTrans, 1, &ws->Rsub.matrix,
            &ws->v_x.vector, 0, &ws->Rs.vector);
    gsl_vector_sub(&ws->Rs.vector, &ws->v_y.vector);
    gsl_blas_ddot(&ws->Rs.vector, &ws->Rs.vector, &d_sqr);

    /*******************************************************************/
    /* ACTUAL ALGORITHM STARTS HERE                                    */
    /*******************************************************************/
    
    // Step 1: set initial values.
    int set_new_bounds = 1;
    size_t solutions = 0;
    size_t k = m - 1;
    ws->d2[k] = d_sqr;
    ws->yhat[k] = calc_yhat(k, ws->R, ws->y, ws->s);

    while(k < m)
    {
        if(set_new_bounds)
        {
            // Step 2: Calculate upper bound for s_k and set s_k to lower bound.
            double Rkk = gsl_matrix_get(ws->R, k, k);
            double d2_sqrt_div_rkk = sqrt(ws->d2[k]) / Rkk;
            double yhat_div_rkk = ws->yhat[k] / Rkk;
            ws->ub[k] = floor(d2_sqrt_div_rkk + yhat_div_rkk);
            ws->s[k] = ceil(-d2_sqrt_div_rkk + yhat_div_rkk);
            // Don't set new bounds on next iteration unless required by step 5b.
            set_new_bounds = 0;
        }
        else
            ws->s[k] += 1.0;    // Step 3

        if(ws->s[k] > ws->ub[k])
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

                gsl_blas_dgemv(CblasNoTrans, 1, &ws->Rsub.matrix, 
                            &ws->v_s.vector, 0, &ws->Rs.vector);
                gsl_vector_sub(&ws->Rs.vector, &ws->v_y.vector);

                double y_minus_Rs_len;
                gsl_blas_ddot(&ws->Rs.vector, &ws->Rs.vector, &y_minus_Rs_len);

                // If new s gives a better solution to the CVP, store it to best_s.
                if(y_minus_Rs_len < d_sqr)
                {
                    gsl_vector_memcpy(&ws->v_x.vector, &ws->v_s.vector);
                    d_sqr = y_minus_Rs_len;
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

    gsl_blas_dgemv(CblasNoTrans, 1, B, &ws->v_x.vector, 0, clp);
    //fprintf(stderr, "SD: Tried %d points\n", solutions);
}

void spheredecode_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{ spheredecode(clp, t, B, ws); }

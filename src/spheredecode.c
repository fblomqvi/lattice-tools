#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "babai.h"
#include "spheredecode.h"
#include "utility.h"
#include <gsl/gsl_blas.h>

/*NOTE: n is used to denote the amount of rows of the basis matrix B and m for
 *      columns. This notation was chosen in order to be consistent with
 *      the documentation and it is the OPPOSITE of how gsl uses them
 *      (in error messages, for example). */

SD_WS* SD_WS_alloc_and_init(const gsl_matrix* B)
{
    // Set n and m to be the amounts of rows and columns of lattice basis B.
    
    int n = B->size1;
    int m = B->size2;
    assert(n >= m);
    //assert(t->size == n);
    
    SD_WS* ws = malloc(sizeof(SD_WS));

    /* Allocate memory for the spheredecoder function.
     * Matrices Q and R store the QR decomposition of B.
     * y is later set to (Q_1^T)*x.
     * Integer vector s stores possible solutions for min( ||t-B*s|| ), 
     * best_s stores the best solution for the CVP found so far.
     * ub stores the upper bound values used in calculating the interval for s_k.
     * d2 and yhat are used in determining the interval in which s_k must be (see
     * calc_d2 and calc_yhat)
     * Rs and Rbest_s store values used to determine wether a new solution is
     * better than the previous best one.*/
 
    ws->Q = gsl_matrix_alloc(n,n);
    ws->R = gsl_matrix_alloc(n,m);    
    ws->y = gsl_vector_alloc(m);
    ws->best_s = gsl_vector_alloc(m);
    ws->s = gsl_vector_alloc(m);
    ws->ub = malloc(m * sizeof(double));
    ws->d2 = malloc(m * sizeof(double));
    ws->yhat = malloc(m * sizeof(double));
    ws->Rs = gsl_vector_alloc(m);
    ws->Rbest_s = gsl_vector_alloc(m);
    ws->babai_clp = gsl_vector_alloc(n);
    ws->babai_ws = BABAI_WS_alloc_and_init(B);

    /* Tau is a vector used only for the gsl QR decomposition function.
     * The size of tau must be min of n,m; this is always m for now.*/
    int tausize = (n < m) ? n : m;
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
    for (int i = 0; i < m; i++) {
        if (gsl_matrix_get(ws->R, i, i) < 0) {
            gsl_vector_view Rrow = gsl_matrix_row(ws->R, i);
            gsl_vector_view Qcol = gsl_matrix_column(ws->Q,i);
            gsl_vector_scale(&Rrow.vector, -1);
            gsl_vector_scale(&Qcol.vector, -1);
        }
    }
     
    return ws;
}

void SD_WS_free(SD_WS *ws)
{
    if(ws)
    {
        gsl_vector_free(ws->best_s);
        gsl_vector_free(ws->s);
        gsl_vector_free(ws->y);
        gsl_matrix_free(ws->Q);
        gsl_matrix_free(ws->R);
        free(ws->ub);
        free(ws->d2);
        free(ws->yhat);
        gsl_vector_free(ws->Rs);
        gsl_vector_free(ws->Rbest_s);
        gsl_vector_free(ws->babai_clp);
        BABAI_WS_free(ws->babai_ws);
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

double calc_yhat(int k, SD_WS *ws) {
    int m = ws->R->size2;
    double sum = 0;
    for (int j = k+1; j < m; j++) {
        sum += gsl_matrix_get(ws->R, k, j) * gsl_vector_get(ws->s, j);
    }
    return gsl_vector_get(ws->y, k) - sum;
}

/* Calculate and return the value of d_k^2. s_k has to satisfy
 * d_k^2 >= (ŷ_k - r_(k,k)*s_k)^2 for the current point to be inside the
 * hypersphere and d_k^2 is calculated as
 * d_k^2 = d_(k+1)^2 - (ŷ_(k+1) - r(k+1,k+1)*s_(k+1))^2
 */
double calc_d2(int k, SD_WS *ws) {
    double res = ws->yhat[k+1] - gsl_matrix_get(ws->R, k+1, k+1) * gsl_vector_get(ws->s, k+1);
    res = -1*res*res + ws->d2[k+1];
    return res;
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
void spheredecode(gsl_vector* clp, const gsl_vector* t, const gsl_matrix *B, double d, SD_WS *ws)
{   
    int n = B->size1;
    int m = B->size2;

    /* Split Q to Q1 and Q2 for the following vector operations. (If Q2 doesn't
     * exist, then Q1 = Q) */
    gsl_matrix_view Q1 = gsl_matrix_submatrix(ws->Q, 0, 0, n, m);
    
    // Set y = (Q1^T)*x.
    gsl_blas_dgemv(CblasTrans, 1, &Q1.matrix, t, 0, ws->y);
    
    // Set Q2Tx_len2 = ||(Q2^T)*x||^2.
    double Q2Tx_len2 = 0;
    if (n > m) {
        gsl_matrix_view Q2 = gsl_matrix_submatrix(ws->Q, 0, m, n, n-m);
        gsl_vector *Q2Tx = gsl_vector_alloc(n-m);
        gsl_blas_dgemv(CblasTrans, 1, &Q2.matrix, t, 0, Q2Tx);
        gsl_blas_ddot(Q2Tx, Q2Tx, &Q2Tx_len2);
        gsl_vector_free(Q2Tx);
    }
        
    /*******************************************************************/
    /* ACTUAL ALGORITHM STARTS HERE                                    */
    /*******************************************************************/
    
    // Step 1: set initial values.
    int k = m-1;
    ws->d2[k] = d*d - Q2Tx_len2;
    ws->yhat[k] = calc_yhat(k, ws);
    
    int set_new_bounds = 1;
    int solutions = 0;
    
    if (ws->d2[k] < 0) {
        // If B is not a square matrix, it is possible for the initial distance
        // requirement to be negative. In this case, terminate.
        fprintf(stderr, "Initial radius too small. Terminating. (d = %f, d2_k = %f)\n", d, ws->d2[k]);
        //res = NULL;
    } else {
        while (k < m) {
            if (set_new_bounds) {
                // Step 2: Calculate upper bound for s_k and set s_k to lower bound.
                double Rkk = gsl_matrix_get(ws->R, k, k);
                ws->ub[k] = floor( (sqrt(ws->d2[k]) + ws->yhat[k])/Rkk );
                gsl_vector_set(ws->s, k, ceil( (-sqrt(ws->d2[k]) + ws->yhat[k])/Rkk ) - 1);
                // Don't set new bounds on next iteration unless required by step 5b.
                set_new_bounds = 0;
            }

            // Step 3: Set s_k = s_k + 1.
            gsl_vector_set(ws->s, k, gsl_vector_get(ws->s, k) + 1);

            if (gsl_vector_get(ws->s, k) > ws->ub[k]) {
                // Step 4: Increase k by 1. If k = m + 1, the algorithm will
                //         terminate (by the while loop).
                k++;
            } else {
                // Step 5
                if (k == 0) {
                    // Step 5a: Solution found.
                    solutions++;

                    // Rs = non-zero part of R times s; Rbest_s = non-zero part of R times best_s. 
                    gsl_matrix_view Rsub = gsl_matrix_submatrix(ws->R, 0, 0, m, m);

                    if (solutions == 1) {
                        // First found solution is stored to best_s and Rbest_s is set to R*best_s-y
                        gsl_vector_memcpy(ws->best_s, ws->s);
                        gsl_blas_dgemv(CblasNoTrans, 1, &Rsub.matrix, ws->best_s, 0, ws->Rbest_s);
                        gsl_vector_sub(ws->Rbest_s, ws->y);
                    } else {
                        gsl_blas_dgemv(CblasNoTrans, 1, &Rsub.matrix, ws->s, 0, ws->Rs);
                        gsl_vector_sub(ws->Rs, ws->y);

                        // Rs = R*s-y; Rbest_s = R*best_s-y
                        double y_minus_Rs_len;
                        double y_minus_Rbest_s_len;
                        // y_minus_Rs_len = ||Rs||^2; y_minus_Rbest_s_len = ||Rbest_s||^2
                        gsl_blas_ddot(ws->Rs, ws->Rs, &y_minus_Rs_len);
                        gsl_blas_ddot(ws->Rbest_s, ws->Rbest_s, &y_minus_Rbest_s_len);

                        // If new s gives a better solution to the CVP, store it to best_s.
                        if (y_minus_Rs_len < y_minus_Rbest_s_len) {
                            gsl_vector_memcpy(ws->best_s, ws->s);
                            gsl_vector_memcpy(ws->Rbest_s, ws->Rs);
                        }
                    }
                } else {
                    // Step 5b: Decrease k and calculate new ŷ_k and d_k^2
                    k--;
                    ws->yhat[k] = calc_yhat(k, ws);
                    ws->d2[k] = calc_d2(k, ws);
                    // Calculate new bounds for s_k.
                    set_new_bounds = 1;
                }
            }
        }

        /* Calculate the closest lattice vector to x, i.e., B*s using the solution
         * of min(||t-B*s||) and store it to res. 
         */
        gsl_blas_dgemv(CblasNoTrans, 1, B, ws->best_s, 0, clp);
    }
}

void spheredecode_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{ 
    SD_WS *sd_ws = (SD_WS*) ws;
    babai_g(sd_ws->babai_clp, t, B, sd_ws->babai_ws);
    gsl_vector_sub(sd_ws->babai_clp, t);
    double d = 0;
    gsl_blas_ddot(sd_ws->babai_clp, sd_ws->babai_clp, &d);
    d = sqrt(d);
    if (d == 0) {
        d = 0.001;
    } else {
        d = 1.001*d;
    }
    spheredecode(clp, t, B, d, sd_ws);
}

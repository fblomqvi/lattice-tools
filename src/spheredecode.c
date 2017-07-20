#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "spheredecode.h"
#include "utility.h"
#include <gsl/gsl_blas.h>

/*NOTE: n is used to denote the amount of rows of the basis matrix B and m for
 *      columns. This notation was chosen in order to be consistent with
 *      the documentation and it is the OPPOSITE of how gsl uses them
 *      (in error messages, for example).
 */

/* Calculate the end return the value of ŷ_k. ŷ_k is y_k from which have been
 * subtracted the kth row values of R multiplied by corresponding values of s 
 * excluding the diagonal. Example with k=2:
 *   [d0  #  #  #  #]      [s0]
 *   [ 0 d1  #  #  #]      [s1]
 * R=[ 0  0 d2  a  b]    s=[s2]  =>  ŷ_2 = y_2 - a*s3 - b*s4
 *   [ 0  0  0 d3  #]      [s3]
 *   [ 0  0  0  0 d4]      [s4]
 */
double calc_yhat(int k, gsl_vector *y, gsl_matrix *R, gsl_vector *s) {
    int m = R->size2;
    double sum = 0;
    for (int j = k+1; j < m; j++) {
        sum += gsl_matrix_get(R, k, j) * gsl_vector_get(s, j);
    }
    return gsl_vector_get(y, k) - sum;
}

/* Calculate and return the value of d_k^2. s_k has to satisfy
 * d_k^2 >= (ŷ_k - r_(k,k)*s_k)^2 for the current point to be inside the
 * hypersphere and d_k^2 is calculated as
 * d_k^2 = d_(k+1)^2 - (ŷ_(k+1) - r(k+1,k+1)*s_(k+1))^2
 */
double calc_d2(int k, gsl_vector *y, gsl_matrix *R, gsl_vector *s, double *d2, double *yhat) {
    double res = yhat[k+1] - gsl_matrix_get(R, k+1, k+1) * gsl_vector_get(s, k+1);
    res = -1*res*res + d2[k+1];
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
void spheredecode(gsl_vector* clp, const gsl_vector* t, const gsl_matrix *B, double d)
{
//gsl_vector *spheredecode(gsl_matrix *B, gsl_vector *t, double d) {
    /* The first part of the function does the necessary pre-work for the
     * algorithm. Look below for the actual start of the algorithm.
     */
    
    // Set n and m to be the amounts of rows and columns of lattice basis B.
    int n = B->size1;
    int m = B->size2;
    
    assert(n >= m);
    assert(t->size == n);
    // Tau is a vector used only for the gsl QR decomposition function.
    // The size of tau must be min of n,m; this is always m for now.
    int tausize = (n < m) ? n : m;
    // Compute the QR-decomposition of lattice basis B and store it to matrices
    // Q and R, don't alter B.
    gsl_matrix *B_copy = gsl_matrix_alloc(n, m);
    gsl_matrix_memcpy(B_copy, B);
    gsl_vector *tau = gsl_vector_alloc(tausize);
    gsl_linalg_QR_decomp(B_copy, tau);
    gsl_matrix *Q = gsl_matrix_alloc(n, n);
    gsl_matrix *R = gsl_matrix_alloc(n, m);
    gsl_linalg_QR_unpack(B_copy, tau, Q, R);
    gsl_vector_free(tau);
    gsl_matrix_free(B_copy);
    
    // Make the diagonal of R positive, as required by the algorithm.
    for (int i = 0; i < m; i++) {
        if (gsl_matrix_get(R, i, i) < 0) {
            gsl_vector_view Rrow = gsl_matrix_row(R, i);
            gsl_vector_view Qcol = gsl_matrix_column(Q,i);
            gsl_vector_scale(&Rrow.vector, -1);
            gsl_vector_scale(&Qcol.vector, -1);
        }
    }
    
    /* Split Q to Q1 and Q2 for the following vector operations. (If Q2 doesn't
     * exist, then Q1 = Q)
     */
    gsl_matrix_view Q1 = gsl_matrix_submatrix(Q, 0, 0, n, m);
    
    // Set y = (Q1^T)*x.
    gsl_vector *y = gsl_vector_alloc(Q1.matrix.size2);
    gsl_blas_dgemv(CblasTrans, 1, &Q1.matrix, t, 0, y);
    
    // Set Q2tt_len2 = ||(Q2^T)*x||^2.
    double Q2tt_len2 = 0;
    if (n > m) {
        gsl_matrix_view Q2 = gsl_matrix_submatrix(Q, 0, m, n, n-m);
        gsl_vector *Q2tt = gsl_vector_alloc(Q2.matrix.size2);
        gsl_blas_dgemv(CblasTrans, 1, &Q2.matrix, t, 0, Q2tt);
        gsl_blas_ddot(Q2tt, Q2tt, &Q2tt_len2);
        gsl_vector_free(Q2tt);
    }
    
    /* Uncomment to calculate B = Q*R to check that the QR-decomposition function
     * and altering of R-matrix's diagonal signs worked as intended; that the
     * resulting matrix B is the original matrix B. Used for debugging.
     */
    /*
    printf("Original B:\n");
    print_matrix(B);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, Q, R, 0, B);
    printf("B = QR:\n");
    print_matrix(B);
    */
    
    /* Allocate memory for integer vectors s and best_s and the result vector.
     * s stores possible solutions for min( ||t-B*s|| ), where t is the
     * target vector and B the lattice basis. The values of s change constantly
     * as the algorithm progresses. best_s stores the best solution for the CVP
     * found so far. Eventually the value res = B*best_s is returned as a solution
     * to the CVP.
     */
    gsl_vector *s = gsl_vector_calloc(m);
    gsl_vector *best_s = gsl_vector_calloc(m);
    //gsl_vector *res = gsl_vector_alloc(t->size);
    
    /* Allocate memory for the lists that store the upper bounds, maximum
     * distances and ŷ-values.
     */
    double *ub = calloc(m, sizeof(double));
    double *d2 = calloc(m, sizeof(double));
    double *yhat = calloc(m, sizeof(double));
    
    /*******************************************************************/
    /* ACTUAL ALGORITHM STARTS HERE                                    */
    /*******************************************************************/
    
    // Step 1: set initial values.
    int k = m-1;
    d2[k] = d*d - Q2tt_len2;
    yhat[k] = calc_yhat(k, y, R, s);
    
    int set_new_bounds = 1;
    int solutions = 0;
    
    if (d2[k] < 0) {
        // If B is not a square matrix, it is possible for the initial distance
        // requirement to be negative. In this case, terminate.
        fprintf(stderr, "Initial radius too small. Terminating.\n");
        //res = NULL;
    } else {
        while (k < m) {
            if (set_new_bounds) {
                // Step 2: Calculate upper bound for s_k and set s_k to lower bound.
                double Rkk = gsl_matrix_get(R, k, k);
                ub[k] = floor( (sqrt(d2[k]) + yhat[k])/Rkk );
                gsl_vector_set(s, k, ceil( (-sqrt(d2[k]) + yhat[k])/Rkk ) - 1);
                // Don't set new bounds on next iteration unless required by step 5b.
                set_new_bounds = 0;
            }

            // Step 3: Set s_k = s_k + 1.
            gsl_vector_set(s, k, gsl_vector_get(s, k) + 1);

            if (gsl_vector_get(s, k) > ub[k]) {
                // Step 4: Increase k by 1. If k = m + 1, the algorithm will
                //         terminate (by the while loop).
                k++;
            } else {
                // Step 5
                if (k == 0) {
                    // Step 5a: Solution found.
                    solutions++;

                    // Calculate R*s and compare ||y-R*s|| to previous best solution.
                    gsl_vector *Rs = gsl_vector_alloc(s->size);
                    gsl_vector *Rbest_s = gsl_vector_alloc(s->size);

                    // Rs = non-zero part of R times s; Rbest_s = non-zero part of R times best_s.
                    gsl_matrix_view Rsub = gsl_matrix_submatrix(R, 0, 0, m, m);
                    gsl_blas_dgemv(CblasNoTrans, 1, &Rsub.matrix, s, 0, Rs);
                    gsl_blas_dgemv(CblasNoTrans, 1, &Rsub.matrix, best_s, 0, Rbest_s);

                    // Rs = R*s-y; Rbest_s = R*best_s-y
                    gsl_vector_sub(Rs, y);
                    gsl_vector_sub(Rbest_s, y);
                    double y_minus_Rs_len;
                    double y_minus_Rbest_s_len;
                    // y_minus_Rs_len = ||Rs||^2; y_minus_Rbest_s_len = ||Rbest_s||^2
                    gsl_blas_ddot(Rs, Rs, &y_minus_Rs_len);
                    gsl_blas_ddot(Rbest_s, Rbest_s, &y_minus_Rbest_s_len);

                    // If new s gives a better solution to the CVP, store it to best_s.
                    if (y_minus_Rs_len < y_minus_Rbest_s_len) {
                        gsl_vector_memcpy(best_s, s);
                    }

                    gsl_vector_free(Rs);
                    gsl_vector_free(Rbest_s);
                } else {
                    // Step 5b: Decrease k and calculate new ŷ_k and d_k^2
                    k--;
                    yhat[k] = calc_yhat(k, y, R, s);
                    d2[k] = calc_d2(k, y, R, s, d2, yhat);
                    // Calculate new bounds for s_k.
                    set_new_bounds = 1;
                }
            }
        }

        /* Calculate the closest lattice vector to x, i.e., B*s using the solution
         * of min(||t-B*s||) and store it to res. 
         */
        gsl_blas_dgemv(CblasNoTrans, 1, B, best_s, 0, clp);
        /*
        printf("Sphere decoder found %d solution", solutions);
        if (solutions != 1) printf("s");
        printf(" within d=%f.\n", d);
        */
    }
    
    /* Free all memory allocated by the spheredecode-function, except for
     * the result vector res.
     */
    gsl_vector_free(best_s);
    gsl_vector_free(s);
    gsl_vector_free(y);
    gsl_matrix_free(Q);
    gsl_matrix_free(R);
    free(d2);
    free(ub);
    free(yhat);
    
    //return res;
}

void spheredecode_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{ 
    double d = 5.0;
    spheredecode(clp, t, B, d);
}

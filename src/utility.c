#include "dbg.h"
#include "lt_errno.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_blas.h>
#include "utility.h"

/* Allocates memory for a new vector c, copies the contents of the argument
 * vector v into c and returns a pointer to c.
 */
/*
gsl_vector *clone_vector(const gsl_vector *v) {
    gsl_vector *c = gsl_vector_alloc(v->size);
    gsl_vector_memcpy(c, v);
    return c;
}
*/

int utility_compute_QR_decomposition(gsl_matrix* Q, gsl_matrix* R, const gsl_matrix* B)
{
    int lt_errno = LT_SUCCESS;
    size_t n = B->size1;
    size_t m = B->size2;
    size_t tausize = (n < m) ? n : m;
    gsl_matrix* B_copy = gsl_matrix_alloc(n, m);
    libcheck_se_mem(B_copy, lt_errno, GSL_ENOMEM);

    gsl_vector* tau = gsl_vector_alloc(tausize);
    llibcheck_se_mem(B_copy, error_a, lt_errno, GSL_ENOMEM);

    gsl_matrix_memcpy(B_copy, B);
    lt_errno = gsl_linalg_QR_decomp(B_copy, tau);
    lt_llibcheck(lt_errno, error_b, "gsl_linalg_QR_decomp failed");
    lt_errno = gsl_linalg_QR_unpack(B_copy, tau, Q, R);
    lt_llibcheck(lt_errno, error_b, "gsl_linalg_QR_unpack failed");

    // Make the diagonal of R positive, as required by the algorithms.
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

error_b:
    gsl_vector_free(tau);
error_a:
    gsl_matrix_free(B_copy);
error:
    return lt_errno;
}

int utility_Rmm_is_not_singular(const gsl_matrix* R, double epsilon)
{
    for(size_t i = 0; i < R->size2; i++)
        if(fabs(gsl_matrix_get(R, i, i)) < epsilon)
            return 0;

    return 1;
}

void print_matrix(gsl_matrix *A) {
    for (size_t i = 0; i < A->size1; i++) {
        gsl_vector_view w = gsl_matrix_row(A, i);
        for (size_t j = 0; j < w.vector.size; j++) {
            printf("%10.10g ", gsl_vector_get(&w.vector, j));
        }
        printf("\n");
    }
}

void fprintf_matrix(FILE *f, gsl_matrix *A) {
    for (size_t i = 0; i < A->size1; i++) {
        gsl_vector_view w = gsl_matrix_row(A, i);
        for (size_t j = 0; j < w.vector.size; j++) {
            fprintf(f, "%f ", gsl_vector_get(&w.vector, j));
        }
        fprintf(f, "\n");
    }
}

void fprintf_result(FILE *res_file, gsl_matrix *B, gsl_vector *t, gsl_vector *res) {
    fprintf(res_file, "B=\n");
    fprintf_matrix(res_file, B);
    fprintf(res_file, "t=\n");
    gsl_vector_fprintf(res_file, t, "%f");
    fprintf(res_file, "res=\n");
    gsl_vector_fprintf(res_file, res, "%f");
    fprintf(res_file, "\n------------------------------------\n\n");
}

/* Create and return a new gsl_matrix from file f and return a pointer to it.
 * Reads the file until the ending brackets ]]. If the ]] are followed
 * by T, the returned matrix is the transpose of the one in file f:
 * [[1 2] returns 1 2 and [[1 2]  returns 1 3
 * [3 4]]         3 4     [3 4]]T         2 4
 */
/*
gsl_matrix *read_matrix(FILE *f, int transpose) {
    char previous = fgetc(f);
    char current = fgetc(f);
    int entries = 0;
    int rows = 0;
    int read = 2;

    while(!(current == ']' && previous == ']')) {
        if (isdigit(previous) && !(isdigit(current) || current == '.')) entries++;
        if (current == ']') rows++;
        previous = current;
        current = fgetc(f);
        read++;
        if (feof(f)) {
            fprintf(stderr, "Matrix file ended before closing ]] was found.\n");
            return NULL;
        }
    }

    fseek(f, -read, SEEK_CUR);

    int cols = entries/rows;
    if (entries % rows != 0) {
        fprintf(stderr, "Amount of entries does not match matrix dimensions (expected %d, found %d).\n", (cols+1)*rows, entries);
        return NULL;
    }

    gsl_matrix* M;
    if(transpose)
    {
        M = gsl_matrix_alloc(cols, rows);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                double current;
                while(!fscanf(f, "%lf", &current)) {
                    fseek(f, 1, SEEK_CUR);
                }
                gsl_matrix_set(M, j, i, current);
            }
        }
    }
    else
    {
        M = gsl_matrix_alloc(rows, cols);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                double current;
                while(!fscanf(f, "%lf", &current)) {
                    fseek(f, 1, SEEK_CUR);
                }
                gsl_matrix_set(M, i, j, current);
            }
        }
    }

    return M;
}
*/

/* Create and return a new gsl_vector from file f and return a pointer to it.
 * Reads the file until the first ]-character.
 */
/*
gsl_vector *read_vector(FILE *f) {
    char previous = '\0';
    char current = fgetc(f);
    int entries = 0;
    int read = 1;

    while (current != ']') {
        if (isdigit(previous) && !(isdigit(current) || current == '.')) entries++;
        previous = current;
        current = fgetc(f);
        read++;
        if (feof(f)) {
            fprintf(stderr, "Vector file ended before closing ] was found.\n");
            return NULL;
        }
    }
    entries++;

    fseek(f, -read, SEEK_CUR);
    gsl_vector *v = gsl_vector_alloc(entries);

    for (int i = 0; i < entries; i++) {
        double current;
        while(!fscanf(f, "%lf", &current)) {
            fseek(f, 1, SEEK_CUR);
        }
        gsl_vector_set(v, i, current);
    }
    return v;
}
*/

/* Generates a q-ary lattice from the parity-check matrix H. H has to be of the
 * form H = [A I] where I is the identity matrix.
 * Outputs B = [ I  0]
 *             [-A qI]
 */
/*
gsl_matrix *generate_basis(gsl_matrix *H, int q) {

    // Set the dimensions of the corresponding generator matrix G.
    int n = H->size2;
    int k = n - H->size1;

    // Check that the parity-check matrix H is of the correct form.
    gsl_matrix_view I_view = gsl_matrix_submatrix(H, 0, k, n-k, n-k);
    gsl_matrix *I = gsl_matrix_alloc(H->size1, H->size1);
    gsl_matrix_set_identity(I);
    assert(gsl_matrix_equal(&I_view.matrix, I));
    gsl_matrix_free(I);

    // Clone the non-identity part of H to A.
    gsl_matrix_view A_view = gsl_matrix_submatrix(H, 0, 0, n-k, k);
    gsl_matrix *A = gsl_matrix_alloc(A_view.matrix.size1, A_view.matrix.size2);
    gsl_matrix_memcpy(A, &A_view.matrix);

    // Set G_sub to be the additive inverse of the non-identity part of H.
    for (size_t i = 0; i < A->size1; i++) {
        for (size_t j = 0; j < A->size2; j++) {
            int aij = (int) round(gsl_matrix_get(A, i, j));
            if (aij % q == 0) {
                gsl_matrix_set(A, i, j, 0);
            } else {
                gsl_matrix_set(A, i, j, q - (aij % q));
            }
        }
    }

    // Allocate memory for the basis matrix (the result) and...
    gsl_matrix *B = gsl_matrix_alloc(n, n);
    // ...set diagonal to identity.
    gsl_matrix_set_identity(B);

    // ...set lower left part to be the additive inverse of the non-identity part of H.
    for (int i = k; i < n; i++) {
        for (int j = 0; j < k; j++) {
            gsl_matrix_set(B, i, j, gsl_matrix_get(A, i-k, j));
        }
    }

    // ...set the lower right part to be q times identity.
    for (int i = k; i < n; i++) {
            gsl_matrix_set(B, i, i, q);
    }

    gsl_matrix_free(A);

    return B;
}
*/

/* Calculates the value of µ = <v,u>/||u||^2 used in calculating the projection
 * of v onto u.
 */
/*
double calc_mu (const gsl_vector *v, const gsl_vector *u) {
    double inner_product, u_length_squared, mu;
    gsl_blas_ddot(v, u, &inner_product);    // Sets inner_product = <v,u>.
    gsl_blas_ddot(u, u, &u_length_squared); // Sets u_length_squared = <u,u> = ||u||^2.
    mu = inner_product/u_length_squared;    // Sets µ = <v,u>/||u||^2.
    return mu;
}
*/

/* Allocates memory for and returns the vector that is the projection of v onto u:
 * projection(gsl_vector *v, gsl_vector *u) = µ*u = <v,u>/||u||^2 * u
 */
/*
gsl_vector *projection(const gsl_vector *v, const gsl_vector *u) {
    assert(v->size == u->size);
    double mu = calc_mu(v, u);
    gsl_vector *res = clone_vector(u);
    gsl_vector_scale(res, mu);
    return res;
}
*/

/* Calculates the orthogonal non-normalized basis of the input B using the
 * modified Gram-Schmidt process. Memory is allocated for the result matrix
 * and the original input matrix is left untouched. Is used as a subroutine
 * in the Babai-algorithm.
 */
/*
gsl_matrix *gram_schmidt(const gsl_matrix *B) {
    // Allocate memory for the result matrix.
    gsl_matrix *res = gsl_matrix_alloc(B->size1, B->size2);

    for (size_t i = 0; i < B->size2; i++) {
        //Initialize new vector v and copy i:th column of input
        // matrix B to vector v.
        gsl_vector_const_view b = gsl_matrix_const_column(B, i);
        gsl_vector *v = clone_vector(&b.vector);

        for (size_t j = 0; j < i; j++) {
            // u is the j:th already calculated orthogonal basis vector.
            gsl_vector_const_view u = gsl_matrix_const_column(res, j);
            if (!gsl_vector_isnull(&u.vector)) {
                // Calculate the projection of v onto u.
                gsl_vector *proj_v_to_u = projection(v, &u.vector);
                // Subtract the projection vector from v.
                gsl_vector_sub(v, proj_v_to_u);
                // Free memory allocated by the projection function.
                gsl_vector_free(proj_v_to_u);
            }
        }
        // Set the i:th column of the result matrix to be v.
        gsl_matrix_set_col(res, i, v);
        // Free memory allocated to v.
        gsl_vector_free(v);
    }

    return res;
}
*/

/*
int is_singular(gsl_matrix *B) {
    int n = B->size1;
    int m = B->size2;
    gsl_matrix *U = gsl_matrix_alloc(n, m);
    gsl_matrix_memcpy(U, B);
    gsl_vector *s = gsl_vector_alloc(m);
    gsl_matrix *V = gsl_matrix_alloc(m, m);
    gsl_vector *work = gsl_vector_alloc(m);
    gsl_linalg_SV_decomp(U, V, s, work);

    gsl_matrix_free(U);
    gsl_matrix_free(V);
    gsl_vector_free(work);

    for (int i = 0; i < m; i++) {
        if (gsl_vector_get(s, i) < 0.01) {
            gsl_vector_free(s);
            return 1;
        }
    }
    gsl_vector_free(s);
    return 0;
}

gsl_matrix *random_basis(int n, int m, gsl_rng *rng, int maxval) {
    int rngmax = maxval+1;
    gsl_matrix *B = gsl_matrix_alloc(n,m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            gsl_matrix_set(B, i, j, gsl_rng_uniform_int(rng, rngmax));
        }
    }

    if (is_singular(B)) {
        gsl_matrix_free(B);
        return random_basis(n, m, rng, maxval);
    } else {
        return B;
    }
}

gsl_matrix *random_basis_s(int m, gsl_rng *rng, int maxval) {
    int rngmax = maxval+1;
    gsl_matrix *B = gsl_matrix_alloc(m+1,m);
    for (int i = 0; i < m; i++) {
        gsl_matrix_set(B, 0, i, gsl_rng_uniform_int(rng, rngmax));
    }
    for (int j = 1; j < m+1; j++) {
        gsl_matrix_set(B, j, j-1, 1);
    }
    return B;
}

gsl_vector *random_target(int n, gsl_rng *rng, int maxval) {
    int rngmax = maxval+1;
    gsl_vector *t = gsl_vector_alloc(n);
    for (int i = 0; i < n; i++) {
        gsl_vector_set(t, i, gsl_rng_uniform_int(rng, rngmax));
    }
    return t;
}
*/

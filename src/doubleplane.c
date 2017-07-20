#include <math.h>
#include <assert.h>
#include "utility.h"
#include <gsl/gsl_blas.h>

gsl_vector *doubleplane_helper(gsl_matrix *B, gsl_matrix *W, gsl_vector *t) {    
    int n = B->size1;
    int m = B->size2;
    
    // Get the last column vectors of B and W into vectors b and w.
    gsl_vector *b = gsl_vector_alloc(n);
    gsl_vector *w = gsl_vector_alloc(n);
    gsl_vector *res = gsl_vector_alloc(n);
    gsl_matrix_get_col(b, B, m-1);
    gsl_matrix_get_col(w, W, m-1);
    // Calculate mu = <t,w>/<w,w>
    double mu = calc_mu(t, w);
    
    if (m == 1) {
        // Only one basis vector left, closest lattice point is round(µ)*b_1
        gsl_vector_scale(b, round(mu));
        gsl_vector_memcpy(res, b);
    } else {
        double f = floor(mu);
        double c = ceil(mu);

        // Set s1 = f*b and s2 = c*b
        gsl_vector *s1 = clone_vector(b);
        gsl_vector *s2 = clone_vector(b);
        gsl_vector_scale(s1, f);
        gsl_vector_scale(s2, c);

        // Set new target vectors newt1 = t-s1 and newt2 = t-s2
        gsl_vector *newt1 = clone_vector(t);
        gsl_vector *newt2 = clone_vector(t);
        gsl_vector_sub(newt1, s1);
        gsl_vector_sub(newt2, s2);

        // Create a submatrix view of B and W. The views contain all but the last
        // column of the original matrix.
        gsl_matrix_view Bsub = gsl_matrix_submatrix(B, 0, 0, n, m-1);
        gsl_matrix_view Wsub = gsl_matrix_submatrix(W, 0, 0, n, m-1);

        // Set x1 and x2 to be the solutions to m-1 -dimensional CVP with targets
        // being newt1 and newt2 (t-floor(µ)*b and t-ceil(µ)*b).
        gsl_vector *x1 = doubleplane_helper(&Bsub.matrix, &Wsub.matrix, newt1);
        gsl_vector *x2 = doubleplane_helper(&Bsub.matrix, &Wsub.matrix, newt2);

        // Set x1=x1+s1, x2=x2+s2 and newt1=t-x1, newt2=t-x2.
        gsl_vector_memcpy(newt1, t);
        gsl_vector_memcpy(newt2, t);
        gsl_vector_add(x1, s1);
        gsl_vector_add(x2, s2);
        gsl_vector_sub(newt1, x1);
        gsl_vector_sub(newt2, x2);

        // Calculate the norms of t-x1 and t-x2.
        double newt1normsquared = 0;
        double newt2normsquared = 0;
        gsl_blas_ddot(newt1, newt1, &newt1normsquared);
        gsl_blas_ddot(newt2, newt2, &newt2normsquared);

        // Set result to be x1 or x2, whichever is closer to t.
        if (newt1normsquared <= newt2normsquared) {
            gsl_vector_memcpy(res, x1);
        } else {
            gsl_vector_memcpy(res, x2);
        }
        
        // Free memory.
        gsl_vector_free(s1);
        gsl_vector_free(s2);
        gsl_vector_free(x1);
        gsl_vector_free(x2);
        gsl_vector_free(newt1);
        gsl_vector_free(newt2);
    }
    gsl_vector_free(b);
    gsl_vector_free(w);

    return res;
}

void doubleplane(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B)
{
    assert(B->size1 == t->size);
    
    // Calculate orthogonal basis for span(B) using the modified Gram-Schmidt process.
    gsl_matrix *W = gram_schmidt(B);
    
    // Call the helper function, which does all the work.
    gsl_vector *res = doubleplane_helper(B, W, t);

    // Copy the result to clp. This is a bad workaround done in order to obtain
    // the correct function signature.
    gsl_vector_memcpy(clp, res);
    
    // Free Gram-Schmidt basis matrix.
    gsl_matrix_free(W);
    gsl_vector_free(res);
}

void doubleplane_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{ doubleplane(clp, t, B); }


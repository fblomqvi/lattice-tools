#include "dbg.h"
#include "babai.h"
#include <math.h>
#include <assert.h>
#include "utility.h"
#include <gsl/gsl_blas.h>

struct s_babai_ws
{
    gsl_matrix* W;
    gsl_vector* t;
};

BABAI_WS* BABAI_WS_alloc_and_init(const gsl_matrix* B)
{
    BABAI_WS* ws = malloc(sizeof(BABAI_WS));
    libcheck_mem(ws);

    // Calculate orthogonal basis for span(B) using the modified Gram-Schmidt process.
    ws->W = gram_schmidt(B);
    llibcheck_mem(ws->W, error_a);

    ws->t = gsl_vector_alloc(B->size1);
    llibcheck_mem(ws->t, error_b);
    return ws;

error_b:
    gsl_matrix_free(ws->W);
error_a:
    free(ws);
error:
    return NULL;
}

void BABAI_WS_free(BABAI_WS* ws)
{
    if(ws)
    {
        gsl_matrix_free(ws->W);
        gsl_vector_free(ws->t);
        free(ws);
    }
}

/* Calculates an approximation for the nearest lattice point in the lattice
 * with basis B to target vector orig_t by repeatedly projecting the target
 * vector to the nearest hyperplane defined by the Gram-Schmidt orthogonalization
 * of B. Allocates memory for the result vector, the original matrix B and
 * vector orig_t are left untouched.
 */
void babai(gsl_vector* clp, const gsl_vector* orig_t, const gsl_matrix* B, BABAI_WS* ws)
{
    assert(B->size1 == orig_t->size);
    
    gsl_vector_set_zero(clp);
    gsl_vector_memcpy(ws->t, orig_t);
    
    for(int i = B->size2-1; i >= 0; i--)
    {
        /* Set b to be the i:th column of the lattice basis matrix B and
         * w the i:th column of the Gram-Schmidt orthogonalized basis matrix W.
         */
        gsl_vector_const_view v_b = gsl_matrix_const_column(B, i);
        gsl_vector_const_view v_w = gsl_matrix_const_column(ws->W, i);
        
        // µ = <t,w>/||w||^2 and c is µ rounded to the closest integer.
        double c = round(calc_mu(ws->t, &v_w.vector));
        
        /* Calculate the new target t = t - c*b and add
         * c*b to result: res = res + c*b.
         */
        gsl_blas_daxpy(-c, &v_b.vector, ws->t); // t = t - c * b
        gsl_blas_daxpy(c, &v_b.vector, clp);    // clp = clp - c * b
    }
}

void babai_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{ babai(clp, t, B, ws); }

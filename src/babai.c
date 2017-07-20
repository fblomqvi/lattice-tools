#include "dbg.h"
#include "babai.h"
#include <math.h>
#include <assert.h>
#include "utility.h"

struct s_babai_ws
{
    gsl_matrix* W;
    gsl_vector* b;
    gsl_vector* t;
    gsl_vector* w;
};

BABAI_WS* BABAI_WS_alloc_and_init(const gsl_matrix* B)
{
    BABAI_WS* ws = malloc(sizeof(BABAI_WS));
    libcheck_mem(ws);

    // Calculate orthogonal basis for span(B) using the modified Gram-Schmidt process.
    ws->W = gram_schmidt(B);
    llibcheck_mem(ws->W, error_a);

    ws->b = gsl_vector_alloc(B->size1);
    llibcheck_mem(ws->b, error_b);

    ws->t = gsl_vector_alloc(B->size1);
    llibcheck_mem(ws->t, error_c);

    ws->w = gsl_vector_alloc(B->size1);
    llibcheck_mem(ws->w, error_d);

    return ws;

error_d:
    gsl_vector_free(ws->t);
error_c:
    gsl_vector_free(ws->b);
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
        gsl_vector_free(ws->b);
        gsl_vector_free(ws->t);
        gsl_vector_free(ws->w);
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
    
    for (int i = B->size2-1; i >= 0; i--) {
        /* Set b to be the i:th column of the lattice basis matrix B and
         * w the i:th column of the Gram-Schmidt orthogonalized basis matrix W.
         */
        gsl_matrix_get_col(ws->b, B, i);
        gsl_matrix_get_col(ws->w, ws->W, i);
        
        // µ = <t,w>/||w||^2 and c is µ rounded to the closest integer.
        double mu = calc_mu(ws->t,ws->w);
        double c = round(mu);
        
        /* Calculate the new target t = t - (µ-c)*w - c*b and add
         * c*b to result: res = res + c*b.
         */
        gsl_vector_scale(ws->w, mu-c);  //w = (µ-c)*w
        gsl_vector_scale(ws->b, c);     //b = c*b
        
        gsl_vector_sub(ws->t, ws->w);   //t = t - w
        gsl_vector_sub(ws->t, ws->b);   //t = t - b
        gsl_vector_add(clp, ws->b); //res = res + b
    }
}

void babai_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{ babai(clp, t, B, ws); }

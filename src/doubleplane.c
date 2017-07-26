#include <math.h>
#include <assert.h>
#include <gsl/gsl_blas.h>
#include "utility.h"
#include "doubleplane.h"


DP_WS* DP_WS_alloc_and_init(const gsl_matrix *B)
{
    DP_WS *ws = malloc(sizeof(DP_WS));

    ws->W = gram_schmidt(B);
    ws->w = gsl_vector_alloc(B->size1);
    ws->b = gsl_vector_alloc(B->size1);
    ws->tx1 = gsl_vector_alloc(B->size1);
    ws->tx2 = gsl_vector_alloc(B->size1);
    ws->s1 = malloc(B->size2 * sizeof(gsl_vector *));
    ws->s2 = malloc(B->size2 * sizeof(gsl_vector *));
    ws->x1 = malloc(B->size2 * sizeof(gsl_vector *));
    ws->x2 = malloc(B->size2 * sizeof(gsl_vector *));
    ws->targs = malloc(B->size2 * sizeof(gsl_vector *));

    for (size_t i = 0; i < B->size2; i++) {
        ws->s1[i] = gsl_vector_alloc(B->size1);
        ws->s2[i] = gsl_vector_alloc(B->size1);
        ws->x1[i] = gsl_vector_alloc(B->size1);
        ws->x2[i] = gsl_vector_alloc(B->size1);
        ws->targs[i] = gsl_vector_alloc(B->size1);
    }
    return ws;
}

void DP_WS_free(DP_WS *ws) {
    if (ws) {
        size_t m = ws->W->size2;
        gsl_matrix_free(ws->W);
        gsl_vector_free(ws->w);
        gsl_vector_free(ws->b);
        gsl_vector_free(ws->tx1);
        gsl_vector_free(ws->tx2);
        for (size_t i = 0; i < m; i++) {
            gsl_vector_free(ws->s1[i]);
            gsl_vector_free(ws->s2[i]);
            gsl_vector_free(ws->x1[i]);
            gsl_vector_free(ws->x2[i]);
            gsl_vector_free(ws->targs[i]);
        }
        free(ws->s1);
        free(ws->s2);
        free(ws->x1);
        free(ws->x2);
        free(ws->targs);
        free(ws);
    }
}

/* Doubleplane helper function. The i argument (which is either 0 or 1) decides whether the result is
 * stored to x1 or x2. If the call is the top level call (i.e. m == B->size2), the result vector
 * is stored to clp. */
void doubleplane_helper(gsl_vector *clp, const gsl_matrix *B, const gsl_vector *t, 
                        DP_WS *ws, size_t m, int i) 
{
    // Get the last column vectors of B and W into vectors b and w.
    gsl_matrix_get_col(ws->b, B, m-1);
    gsl_matrix_get_col(ws->w, ws->W, m-1);
    gsl_vector_memcpy(ws->targs[m-1], t);
    // Calculate mu = <t,w>/<w,w>
    double mu = calc_mu(t, ws->w);
    
    if (m == 1) {
        // Only one basis vector left, closest lattice point is round(µ)*b_1
        gsl_vector_scale(ws->b, round(mu));
        if (i == 0) {
            gsl_vector_memcpy(ws->x1[1], ws->b);
        } else {
            gsl_vector_memcpy(ws->x2[1], ws->b);
        }
    } else {
        double f = floor(mu);
        double c = ceil(mu);

        // Set s1 = f*b and s2 = c*b
        gsl_vector_memcpy(ws->s1[m-1], ws->b);
        gsl_vector_memcpy(ws->s2[m-1], ws->b);
        gsl_vector_scale(ws->s1[m-1], f);
        gsl_vector_scale(ws->s2[m-1], c);

        /* Set new target vector tx1 = t-s1 and call doubleplane
         * recursively using tx1 as the new target. The result vector
         * of the recursive call will be stored in x1[m-1] */
        gsl_vector_memcpy(ws->tx1, ws->targs[m-1]);
        gsl_vector_sub(ws->tx1, ws->s1[m-1]);
        doubleplane_helper(clp, B, ws->tx1, ws, m-1, 0);
        // Do the same with tx2
        gsl_vector_memcpy(ws->tx2, ws->targs[m-1]);
        gsl_vector_sub(ws->tx2, ws->s2[m-1]);
        doubleplane_helper(clp, B, ws->tx2, ws, m-1, 1);

        /* Reset tx1 and tx2 to be the original target vectors of 
         * the current step. Set x1=x1+s1, x2=x2+s2 and tx1=t-x1, tx2=t-x2.*/
        gsl_vector_memcpy(ws->tx1, ws->targs[m-1]);
        gsl_vector_memcpy(ws->tx2, ws->targs[m-1]);
        gsl_vector_add(ws->x1[m-1], ws->s1[m-1]);
        gsl_vector_add(ws->x2[m-1], ws->s2[m-1]);
        gsl_vector_sub(ws->tx1, ws->x1[m-1]);
        gsl_vector_sub(ws->tx2, ws->x2[m-1]);

        // Calculate the norms of t-x1 and t-x2.
        double tx1normsquared = 0;
        double tx2normsquared = 0;
        gsl_blas_ddot(ws->tx1, ws->tx1, &tx1normsquared);
        gsl_blas_ddot(ws->tx2, ws->tx2, &tx2normsquared);

        /* If the current call is the top level call, store the result
         * to clp. Otherwise store it to the previous level's x-vector
         * (x1[m] if i == 0 and x2[m] if i == 1) */
        if (m == B->size2) {
            if (tx1normsquared <= tx2normsquared) {
                gsl_vector_memcpy(clp, ws->x1[m-1]);
            } else {
                gsl_vector_memcpy(clp, ws->x2[m-1]);
            }
        } else if (i == 0) {
            if (tx1normsquared <= tx2normsquared) {
                gsl_vector_memcpy(ws->x1[m], ws->x1[m-1]);
            } else {
                gsl_vector_memcpy(ws->x1[m], ws->x2[m-1]);
            }
        } else {
            if (tx1normsquared <= tx2normsquared) {
                gsl_vector_memcpy(ws->x2[m], ws->x1[m-1]);
            } else {
                gsl_vector_memcpy(ws->x2[m], ws->x2[m-1]);
            }
        }
    }
}

void doubleplane(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, DP_WS *ws)
{
    assert(B->size1 == t->size);
    // Call the helper function, which does all the work.
    doubleplane_helper(clp, B, t, ws, B->size2, 0);
}

void doubleplane_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{ doubleplane(clp, t, B, ws); }


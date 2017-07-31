#include <math.h>
#include <assert.h>
#include <gsl/gsl_blas.h>
#include "utility.h"
#include "doubleplane.h"

struct s_dp_ws
{
    gsl_matrix* W;
    gsl_vector* t;
    gsl_vector* clp;
    gsl_vector* x;
    gsl_vector_view* v_B;
    gsl_vector_view* v_W;
    double dist;
    double dist_min;
    int solutions;
};

DP_WS* DP_WS_alloc_and_init(const gsl_matrix *B)
{
    DP_WS *ws = malloc(sizeof(DP_WS));

    ws->W = gram_schmidt(B);
    ws->t = gsl_vector_alloc(B->size1);
    ws->x = gsl_vector_alloc(B->size1);
    ws->clp = gsl_vector_alloc(B->size1);
    ws->v_B = malloc(B->size2 * sizeof(gsl_vector_view));
    ws->v_W = malloc(B->size2 * sizeof(gsl_vector_view));

    for(size_t i = 0; i< B->size2; i++)
    {
        ws->v_B[i] = gsl_matrix_column(B, i);
        ws->v_W[i] = gsl_matrix_column(ws->W, i);
    }
    return ws;
}

void DP_WS_free(DP_WS *ws) 
{
    if (ws) 
    {
        gsl_matrix_free(ws->W);
        free(ws->v_B);
        free(ws->v_W);
        free(ws->clp);
        free(ws->x);
        free(ws->t);
        free(ws);
    }
}

static void vector_sub(gsl_vector* r, const gsl_vector* a, const gsl_vector* b)
{
    assert(r->size == a->size && r->size == b->size);
    for(size_t i = 0; i < a->size; i++)
        gsl_vector_set(r, i, gsl_vector_get(a, i) - gsl_vector_get(b, i));
}

void doubleplane_helper(gsl_vector* clp, const gsl_vector* t, DP_WS* ws, size_t m) 
{
    if(m == 0)
    {
        double c = round(calc_mu(ws->t, &ws->v_W[m].vector));
        gsl_blas_daxpy(c, &ws->v_B[m].vector, ws->clp);    // clp = clp + c * b

        if(ws->solutions)
        {
            // Caclulate ||t-x|| and replace clp with the new value if needed.
            vector_sub(ws->x, ws->clp, t);
            gsl_blas_ddot(ws->x, ws->x, &ws->dist);
            if(ws->dist < ws->dist_min)
            {
                ws->dist_min = ws->dist;
                gsl_vector_memcpy(clp, ws->clp);
            }
        }
        else
        {
            vector_sub(ws->x, ws->clp, t);
            gsl_blas_ddot(ws->x, ws->x, &ws->dist_min);
            gsl_vector_memcpy(clp, ws->clp);
            ws->solutions = 1;
        }

        // Reset ws->clp
        gsl_blas_daxpy(-c, &ws->v_B[m].vector, ws->clp);    // clp = clp + c * b
    }
    else
    {
        // The floor version
        double c = floor(calc_mu(ws->t, &ws->v_W[m].vector));
        gsl_blas_daxpy(-c, &ws->v_B[m].vector, ws->t); // t = t - c * b
        gsl_blas_daxpy(c, &ws->v_B[m].vector, ws->clp);    // clp = clp - c * b

        doubleplane_helper(clp, t, ws, m-1);

        gsl_blas_daxpy(-1.0, &ws->v_B[m].vector, ws->t); // t = t - b
        gsl_blas_daxpy(1.0, &ws->v_B[m].vector, ws->clp);    // clp = clp + b
        doubleplane_helper(clp, t, ws, m-1);

        // Reset ws->clp and ws->t
        gsl_blas_daxpy(c + 1, &ws->v_B[m].vector, ws->t);
        gsl_blas_daxpy(-(c + 1), &ws->v_B[m].vector, ws->clp);
    }
}

void doubleplane(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, DP_WS *ws)
{
    assert(B->size1 == t->size);

    gsl_vector_set_zero(ws->clp);
    gsl_vector_memcpy(ws->t, t);
    ws->solutions = 0;
    doubleplane_helper(clp, t, ws, B->size2-1);
}

void doubleplane_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws)
{ doubleplane(clp, t, B, ws); }


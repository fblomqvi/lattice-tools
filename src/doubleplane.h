#ifndef DOUBLEPLANE_H
#define DOUBLEPLANE_H

#include <gsl/gsl_linalg.h>

typedef struct s_dp_ws
{
    gsl_matrix *W;
    gsl_vector *w;
    gsl_vector *b;
    gsl_vector **targs;
    gsl_vector **s1;
    gsl_vector **s2;
    gsl_vector **x1;
    gsl_vector **x2;
    gsl_vector *tx1;
    gsl_vector *tx2;
} DP_WS;

//gsl_vector *doubleplane_helper(gsl_matrix *B, gsl_matrix *W, gsl_vector *t);
//gsl_vector *doubleplane(gsl_matrix *B, gsl_vector *t);
DP_WS *DP_WS_alloc_and_init(gsl_matrix *B);
void DP_WS_free(DP_WS *ws);
void doubleplane(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, DP_WS *ws);

void doubleplane_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws);

#endif /* DOUBLEPLANE_H */

#ifndef DOUBLEPLANE_H
#define DOUBLEPLANE_H

#include <gsl/gsl_linalg.h>

//gsl_vector *doubleplane_helper(gsl_matrix *B, gsl_matrix *W, gsl_vector *t);
//gsl_vector *doubleplane(gsl_matrix *B, gsl_vector *t);
void doubleplane(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B);

void doubleplane_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws);

#endif /* DOUBLEPLANE_H */

#ifndef SPHEREDECODE_H
#define SPHEREDECODE_H

#include <gsl/gsl_linalg.h>

typedef struct s_sd_ws
{
    double a;
} SD_WS;

typedef struct s_sd_args
{
    gsl_matrix* B;
    double d;
} SD_ARGS;

//double calc_yhat(int k, gsl_vector *y, gsl_matrix *R, gsl_vector *s);
//double calc_d2(int k, gsl_vector *y, gsl_matrix *R, gsl_vector *s, double *d2, double *yhat);
//gsl_vector *spheredecode(gsl_matrix *B, gsl_vector *t, double d);
void spheredecode(gsl_vector* clp, const gsl_vector* t, const gsl_matrix *B, double d);

void spheredecode_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws);

#endif /* SPHEREDECODE_H */


#ifndef SPHEREDECODE_H
#define SPHEREDECODE_H

#include <gsl/gsl_linalg.h>

typedef struct s_sd_ws
{
    gsl_matrix *Q;
    gsl_matrix *R; 
    gsl_vector *y;
    gsl_vector *best_s;
    gsl_vector *s;
    double *ub;
    double *d2;
    double *yhat;
    gsl_vector *Rs;
    gsl_vector *Rbest_s;
    gsl_vector *babai_clp;
    BABAI_WS *babai_ws;
} SD_WS;

typedef struct s_sd_args
{
    gsl_matrix* B;
    double d;
} SD_ARGS;

SD_WS *SD_WS_alloc_and_init(const gsl_matrix* B);
void SD_WS_free(SD_WS *ws);

void spheredecode(gsl_vector* clp, const gsl_vector* t, const gsl_matrix *B, double d, SD_WS *ws);

void spheredecode_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws);

#endif /* SPHEREDECODE_H */


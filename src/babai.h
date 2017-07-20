#ifndef BABAI_H
#define BABAI_H

#include <gsl/gsl_linalg.h>

typedef struct s_babai_ws BABAI_WS;

BABAI_WS* BABAI_WS_alloc_and_init(const gsl_matrix* B);
void BABAI_WS_free(BABAI_WS* ws);

void babai(gsl_vector* clp, const gsl_vector* orig_t, const gsl_matrix* B, BABAI_WS* ws);
void babai_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws);

#endif /* BABAI_H */

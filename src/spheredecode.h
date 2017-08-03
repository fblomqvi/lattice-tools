#ifndef SPHEREDECODE_H
#define SPHEREDECODE_H

#include <gsl/gsl_linalg.h>

typedef struct s_sd_ws SD_WS;

SD_WS *SD_WS_alloc_and_init(const gsl_matrix* B);
void SD_WS_free(SD_WS *ws);

void spheredecode(gsl_vector* clp, const gsl_vector* t, const gsl_matrix *B, SD_WS *ws);

void spheredecode_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws);

#endif /* SPHEREDECODE_H */


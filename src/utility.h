#ifndef UTILITY_H
#define UTILITY_H

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>

int utility_compute_QR_decomposition(gsl_matrix* Q, gsl_matrix* R, const gsl_matrix* B);
int utility_Rmm_is_not_singular(const gsl_matrix* R, double epsilon);

void print_matrix(gsl_matrix *A);
void fprintf_matrix(FILE *f, gsl_matrix *A);
void fprintf_result(FILE *res_file, gsl_matrix *B, gsl_vector *t, gsl_vector *res);

/*
gsl_matrix *read_matrix(FILE *f, int transpose);
gsl_vector *read_vector(FILE *f);

gsl_vector *clone_vector(const gsl_vector *v);

double calc_mu (const gsl_vector *v, const gsl_vector *u);
gsl_vector *projection(const gsl_vector *v, const gsl_vector *u);
gsl_matrix *gram_schmidt(const gsl_matrix *B);

gsl_matrix *generate_basis(gsl_matrix *H, int q);

int is_singular(gsl_matrix *B);
gsl_matrix *random_basis(int n, int m, gsl_rng *rng, int maxval);
gsl_matrix *random_basis_s(int m, gsl_rng *rng, int maxval);
gsl_vector *random_target(int n, gsl_rng *rng, int maxval);
*/

#endif /* UTILITY_H */


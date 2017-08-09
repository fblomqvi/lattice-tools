/* sphere_pohst.h
   Copyright (C) 2017 Ferdinand Blomqvist, Juuso Korvuo

   This program is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 as published by
   the Free Software Foundation.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along with
   this program. If not, see <http://www.gnu.org/licenses/>.  
   
   Written by Ferdinand Blomqvist and Juuso Korvuo. */

#ifndef FB_LIBLATTICE_SPHERE_POHST_H
#define FB_LIBLATTICE_SPHERE_POHST_H

#include <gsl/gsl_matrix.h>

typedef struct s_sd_ws SD_WS;

SD_WS *SD_WS_alloc_and_init(const gsl_matrix* B);
void SD_WS_free(SD_WS *ws);

void spheredecode(gsl_vector* clp, const gsl_vector* t, const gsl_matrix *B, SD_WS *ws);

void spheredecode_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws);

void sd_dp_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws);

#endif /* FB_LIBLATTICE_SPHERE_POHST_H */


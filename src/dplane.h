/* dplane.h
   Copyright (C) 2017 Ferdinand Blomqvist

   This program is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License version 2 as published by
   the Free Software Foundation. 

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along with
   this program. If not, see <http://www.gnu.org/licenses/>.  
   
   Written by Ferdinand Blomqvist. */

#ifndef FB_LIBLATTICE_DPLANE_H
#define FB_LIBLATTICE_DPLANE_H

#include <gsl/gsl_matrix.h>

typedef struct s_dp_ws DP_WS;

int DP_WS_alloc_and_init(DP_WS** ws_ptr, const gsl_matrix *B);
void DP_WS_free(DP_WS *ws);

void dplane_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws);
void dplane_vanilla_g(gsl_vector* clp, const gsl_vector* t, const gsl_matrix* B, void* ws);

#endif /* FB_LIBLATTICE_DPLANE_H */

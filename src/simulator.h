/* simulator.h
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

#ifndef FB_LDPC_TOOLS_SIMULATOR_H
#define FB_LDPC_TOOLS_SIMULATOR_H

#include "algorithm.h"
#include <gsl/gsl_rng.h>

typedef struct s_simulator SIMULATOR;

typedef struct s_simulator_options
{
    char* infile;
    char* outfile;
    size_t min_err;
    double vnr_begin;
    double vnr_step;
    double vnr_end;
    double bit_err_cutoff;
    double frame_err_cutoff;
    unsigned long seed;
    int zero_cwords;
    const gsl_rng_type* rng_type;
} SIM_OPTIONS;

/*
typedef struct s_simulation_status
{
    size_t dimension;
    size_t frames;
    size_t frame_errs;
    size_t bit_errs;
    double sigma;
    double vnr;
} SIM_STATUS;
*/

/* Takes ownership of the basis matrix. */
int SIMULATOR_from_basis(SIMULATOR** sim_ptr, gsl_matrix* basis, Algorithm alg);

void SIMULATOR_free(SIMULATOR* sim);

int SIMULATOR_run(SIMULATOR* sim, SIM_OPTIONS* opt);

#endif /* FB_LDPC_TOOLS_SIMULATOR_H */

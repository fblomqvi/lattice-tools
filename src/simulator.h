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

enum
{
    SIM_CONTINUE = 0,
    SIM_TERMINATION_VNR = 1,
    SIM_TERMINATION_BER,
    SIM_TERMINATION_FER
};

typedef struct s_simulator SIMULATOR;

typedef struct s_simulator_options
{
    size_t min_err;
    double vnr_begin;
    double vnr_step;
    double vnr_end;
    double ser_cutoff;
    double fer_cutoff;
    unsigned long seed;
    int zero_cwords;
    const gsl_rng_type* rng_type;
} SIM_OPTIONS;

typedef struct s_simulation_status
{
    size_t n;
    size_t m;
    size_t frames;
    size_t frame_errs;
    size_t bit_errs;
    size_t total;
    double rate;
    double sigma;
    double vnr;
    double vol;
    int termination_reason;
} SIM_STATUS;

typedef int (*SimCallback)(const SIM_STATUS*, void*);

/* Takes ownership of the basis matrix. */
int SIMULATOR_from_basis(SIMULATOR** sim_ptr, gsl_matrix* basis, Algorithm alg);

void SIMULATOR_set_callbacks(SIMULATOR* sim,
                            SimCallback vnr_callback,
                            SimCallback start_callback,
                            SimCallback end_callback,
                            void* args);

void SIMULATOR_free(SIMULATOR* sim);

size_t SIMULATOR_get_dimension(const SIMULATOR* sim);
size_t SIMULATOR_get_rank(const SIMULATOR* sim);
double SIMULATOR_get_rate(const SIMULATOR* sim);

int SIMULATOR_run(SIMULATOR* sim, SIM_OPTIONS* opt);

#endif /* FB_LDPC_TOOLS_SIMULATOR_H */

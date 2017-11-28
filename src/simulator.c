/* simulator.c
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

#include "dbg.h"
#include "defs.h"
#include "simulator.h"
#include "channel.h"
#include "rng.h"
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

typedef struct s_decoder_workspace
{
    double* decoded;
    double* received;
    double* transmitted;
    gsl_vector_view v_d;
    gsl_vector_view v_r;
} SIM_WS;

struct s_simulator
{
    gsl_matrix* basis;
    void* alg_ws;
    FILE* outfile;
    CHANNEL* channel;
    void (*get_received)(CHANNEL*);
    SOLVE_func decode;
    size_t (*is_decoding_error)(SIM_WS*, size_t len, size_t*);
    SimCallback vnr_callback;
    SimCallback start_callback;
    SimCallback end_callback;
    void* callback_args;
    size_t dt_size;
    SIM_STATUS status;
    Algorithm alg;
};

static void init_function_pointers(SIMULATOR* sim, SIM_OPTIONS* opt);

static size_t is_not_zero_codeword(SIM_WS* ws, size_t len, size_t* bit_errs);
static size_t is_frame_error(SIM_WS* ws, size_t len, size_t* bit_errs);

static void process_cword(SIMULATOR* sim, SIM_WS* ws, SIM_STATUS* status);

static int read_cword_binary(void* cword, size_t cword_len, size_t dt_size);

static int simulation_termination_status(SIM_STATUS* status, SIM_OPTIONS* opt);
static void simulation_status_init(SIM_STATUS* status, SIMULATOR* sim);

static int simulation_zero_cword_run_vnr(SIMULATOR* sim, SIM_OPTIONS* opt, SIM_WS* ws);
static int simulation_read_cword_run_vnr(SIMULATOR* sim, SIM_OPTIONS* opt, SIM_WS* ws);

static SIM_WS* SIM_WS_alloc(size_t n, size_t size, int zero_cwords);
static void SIM_WS_free(SIM_WS* ws);

static int compute_volume(double* vol, const gsl_matrix* B);

int SIMULATOR_from_basis(SIMULATOR** sim_ptr, gsl_matrix* basis, Algorithm alg)
{
    int lt_errno = LT_SUCCESS;
    SIMULATOR* sim = calloc(1, sizeof(SIMULATOR));
    libcheck_se_mem(sim, lt_errno, LT_ESYSTEM);

    lt_errno = algorithm_get_fp_init_ws(&sim->decode, &sim->alg_ws, alg, basis);
    lt_llibcheck(lt_errno, error_a, "algorithm_get_fp_init_ws failed");

    sim->basis = basis;
    sim->alg = alg;
    sim->status.n =  basis->size1;
    sim->status.m = basis->size2;
    sim->status.rate = (double) sim->status.m / sim->status.n;

    lt_errno = compute_volume(&sim->status.vol, basis);
    lt_llibcheck(lt_errno, error_b, "compute_volume failed");

    debug("rate: %f", sim->status.rate);
    debug("volume: %f", sim->status.vol);

    // Temporary
    sim->dt_size = sizeof(double);
    *sim_ptr = sim;
    return lt_errno;

error_b:
    algorithm_free_ws(sim->alg_ws, sim->alg);
error_a:
    free(sim);
error:
    *sim_ptr = NULL;
    return lt_errno;
}

void SIMULATOR_set_callbacks(SIMULATOR* sim,
                            SimCallback vnr_callback,
                            SimCallback start_callback,
                            SimCallback end_callback,
                            void* args)
{
    sim->vnr_callback = vnr_callback;
    sim->start_callback = start_callback;
    sim->end_callback = end_callback;
    sim->callback_args = args;
}

void SIMULATOR_free(SIMULATOR* sim)
{
    if(sim)
    {
        algorithm_free_ws(sim->alg_ws, sim->alg);
        gsl_matrix_free(sim->basis);
        free(sim);
    }
}

size_t SIMULATOR_get_dimension(const SIMULATOR* sim)
{ return sim->status.n; }

size_t SIMULATOR_get_rank(const SIMULATOR* sim)
{ return sim->status.m; }

double SIMULATOR_get_rate(const SIMULATOR* sim)
{ return sim->status.rate; }

int SIMULATOR_run(SIMULATOR* sim, SIM_OPTIONS* opt)
{
    int ret = LT_FAILURE;
    SIM_WS* ws = SIM_WS_alloc(sim->status.n, sim->dt_size, opt->zero_cwords);
    check_mem(ws);

    gsl_rng* rng = rng_alloc_and_seed(opt->rng_type, opt->seed);
    lcheck_mem(rng, error_a);

    sim->channel = CHANNEL_alloc(ws->received, ws->transmitted, sim->status.n, 0, rng);
    lcheck_mem(sim->channel, error_b);

    init_function_pointers(sim, opt);

    int (*simulation_run_vnr)(SIMULATOR*, SIM_OPTIONS*, SIM_WS*)
        = opt->zero_cwords
        ? simulation_zero_cword_run_vnr : simulation_read_cword_run_vnr;

    sim->status.vnr = opt->vnr_begin;
    sim->status.total = 0;

    int rc;
    if(sim->start_callback)
    {
        rc = sim->start_callback(&sim->status, sim->callback_args);
        llibcheck(rc == 0, error_c, "simulation callback failed");
    }

    while((sim->status.termination_reason =
                simulation_termination_status(&sim->status, opt))
            == SIM_CONTINUE)
    {
        rc = simulation_run_vnr(sim, opt, ws);
        llibcheck(rc == 0, error_c, "simulation_run_vnr failed");
        rc = sim->vnr_callback(&sim->status, sim->callback_args);
        llibcheck(rc == 0, error_c, "simulation callback failed");
        sim->status.vnr += opt->vnr_step;
    }

    if(sim->end_callback)
    {
        rc = sim->end_callback(&sim->status, sim->callback_args);
        llibcheck(rc == 0, error_c, "simulation callback failed");
    }

    ret = LT_SUCCESS;

error_c:
    CHANNEL_free(sim->channel);
error_b:
    gsl_rng_free(rng);
error_a:
    SIM_WS_free(ws);
error:
    return ret;
}

static void init_function_pointers(SIMULATOR* sim, SIM_OPTIONS* opt)
{
    sim->get_received = opt->zero_cwords ?
        channel_get_zero_cword_AWGN : channel_get_cword_AWGN;
    sim->is_decoding_error = opt->zero_cwords ?
            is_not_zero_codeword : is_frame_error;
}


static size_t is_not_zero_codeword(SIM_WS* ws, size_t len, size_t* bit_errs)
{
    double* cword = ws->decoded;
    size_t errors = 0;
    for(size_t i = 0; i < len; i++)
        if(fabs(cword[i]) > EPSILON_EQUAL)
            errors++;

    *bit_errs += errors;
    return (errors ? 1 : 0);
}

static size_t is_frame_error(SIM_WS* ws, size_t len, size_t* bit_errs)
{
    double* cword = ws->decoded;
    double* transmitted = ws->transmitted;
    size_t errors = 0;
    for(size_t i = 0; i < len; i++)
        if(fabs(cword[i] - transmitted[i]) > EPSILON_EQUAL)
            errors++;

    *bit_errs += errors;
    return (errors ? 1 : 0);
}

static void process_cword(SIMULATOR* sim, SIM_WS* ws, SIM_STATUS* status)
{
    sim->get_received(sim->channel);
    sim->decode(&ws->v_d.vector, &ws->v_r.vector, sim->basis, sim->alg_ws);
    status->frame_errs += sim->is_decoding_error(ws, sim->status.n, &status->bit_errs);
    status->frames++;
}

static int read_cword_binary(void* cword, size_t cword_len, size_t dt_size)
{
    size_t rc = fread(cword, dt_size, cword_len, stdin);
    check(rc == cword_len, "fread failed");
    return 0;

error:
    if(feof(stdin))
    {
        log_err_ne("Could not read the next codeword; end of file reached");
        clearerr(stdin);
    }
    return -1;
}

static void simulation_status_init(SIM_STATUS* status, SIMULATOR* sim)
{
    status->frames = 0;
    status->frame_errs = 0;
    status->bit_errs = 0;
    status->sigma = sqrt(sim->status.vol / (status->rate * pow(10.0, status->vnr / 10.0)));
    CHANNEL_set_sigma(sim->channel, status->sigma);
    debug("sigma: %f", status->sigma);
}

static int simulation_termination_status(SIM_STATUS* status, SIM_OPTIONS* opt)
{
    if(status->vnr > (opt->vnr_end + EPSILON_EQUAL))
        return SIM_TERMINATION_VNR;
    if((double) status->frame_errs / status->frames <= opt->fer_cutoff)
        return SIM_TERMINATION_FER;

    const size_t bits = status->frames * status->n;
    const double ber = (double) status->bit_errs / bits;
    if(ber <= opt->ser_cutoff)
        return SIM_TERMINATION_BER;

    return SIM_CONTINUE;
}

static int simulation_zero_cword_run_vnr(SIMULATOR* sim, SIM_OPTIONS* opt, SIM_WS* ws)
{
        SIM_STATUS* status = &sim->status;
        simulation_status_init(status, sim);

        while(sim->status.frame_errs < opt->min_err)
            process_cword(sim, ws, status);

        status->total += status->frames;
        return 0;
}

static int simulation_read_cword_run_vnr(SIMULATOR* sim, SIM_OPTIONS* opt, SIM_WS* ws)
{
        SIM_STATUS* status = &sim->status;
        simulation_status_init(status, sim);

        while(sim->status.frame_errs < opt->min_err)
        {
            int rc = read_cword_binary(ws->transmitted, sim->status.n, sim->dt_size);
            libcheck(rc == 0, "read_cword_binary failed");

            process_cword(sim, ws, status);
        }

        sim->status.total += sim->status.frames;
        return 0;

error:
    return -1;
}

static SIM_WS* SIM_WS_alloc(size_t n, size_t dt_size, int zero_cwords)
{
    SIM_WS* ws = malloc(sizeof(SIM_WS));
    libcheck_mem(ws);

    const size_t t = zero_cwords ? 2 : 3;
    ws->decoded = malloc(t * n * dt_size);
    llibcheck_mem(ws->decoded, error_a);

    ws->received = ws->decoded + n;
    ws->transmitted = zero_cwords ? NULL : ws->received + n;

    ws->v_d = gsl_vector_view_array(ws->decoded, n);
    ws->v_r = gsl_vector_view_array(ws->received, n);
    return ws;

error_a:
    free(ws);
error:
    return NULL;
}

static void SIM_WS_free(SIM_WS* ws)
{
    free(ws->decoded);
    free(ws);
}

static int compute_volume(double* vol, const gsl_matrix* B)
{
    int lt_errno = LT_SUCCESS;
    const size_t m = B->size2;

    gsl_matrix* G = gsl_matrix_alloc(m, m);
    libcheck_se_mem(G, lt_errno, GSL_ENOMEM);

    gsl_permutation* perm = gsl_permutation_alloc(m);
    llibcheck_se_mem(perm, error_a, lt_errno, GSL_ENOMEM);

    lt_errno = gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, B, B, 0.0, G);
    lt_llibcheck(lt_errno, error_b, "gsl_blas_dgemm failed");

    int signum;
    lt_errno = gsl_linalg_LU_decomp(G, perm, &signum);
    lt_llibcheck(lt_errno, error_b, "gsl_linalg_LU_decomp failed");

    *vol = pow(gsl_linalg_LU_det(G, signum), 2.0 / m);

error_b:
    gsl_permutation_free(perm);
error_a:
    gsl_matrix_free(G);
error:
    return lt_errno;
}

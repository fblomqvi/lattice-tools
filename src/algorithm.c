/* algorithm.c
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
#include "lt_errno.h"
#include "algorithm.h"
#include "babai.h"
#include "dplane.h"
#include "sphere_pohst.h"
#include "sphere_se.h"

#define ALG_NAME_BABAI "babai"
#define ALG_NAME_DPLANE "dplane"
#define ALG_NAME_DPLANE_VANILLA "dplane-vanilla"
#define ALG_NAME_SPHERE_DP "dplane-pohst"
#define ALG_NAME_SPHERE_SE "sphere"
#define ALG_NAME_SPHERE_POHST "sphere-pohst"

static const char* alg_names[] = {
    ALG_NAME_BABAI,
    ALG_NAME_DPLANE,
    ALG_NAME_DPLANE_VANILLA,
    ALG_NAME_SPHERE_POHST,
    ALG_NAME_SPHERE_DP,
    ALG_NAME_SPHERE_SE
};

int algorithm_parse_name(const char* name)
{
    for(int i = ALG_BABAI; i < ALG_MAX; i++)
        if(!strcmp(name, alg_names[i])) return i;

    return LT_FAILURE;
}

const char* algorithm_get_name(Algorithm alg)
{ return alg_names[alg]; }

int algorithm_print_names(FILE* file)
{
    libcheck(fprintf(file, "Available algorithms are:\n") > 0, "printing error");
    for(int i = ALG_BABAI; i < ALG_MAX; i++)
        libcheck(fprintf(file, "%s\n", alg_names[i]) > 0, "printing error");

    return LT_SUCCESS;

error:
    return LT_FAILURE;
}

int algorithm_get_fp_init_ws(SOLVE_func* f, void** ws, 
                            Algorithm alg, const gsl_matrix* basis)
{
    int lt_errno;
    switch(alg)
    {
        case ALG_BABAI:
            *f = babai_g;
            lt_errno = BABAI_WS_alloc_and_init((BABAI_WS**) ws, basis);
            break;
        case ALG_DPLANE:
            *f = dplane_g;
            lt_errno = DP_WS_alloc_and_init((DP_WS**) ws, basis);
            break;
        case ALG_DPLANE_VANILLA:
            *f = dplane_vanilla_g;
            lt_errno = DP_WS_alloc_and_init((DP_WS**) ws, basis);
            break;
        case ALG_SPHERE_POHST:
            *f = spheredecode_g;
            lt_errno = SD_WS_alloc_and_init((SD_WS**) ws, basis);
            break;
        case ALG_SPHERE_DP:
            *f = sd_dp_g;
            lt_errno = SD_WS_alloc_and_init((SD_WS**) ws, basis);
            break;
        case ALG_SPHERE_SE:
            *f = sphere_se_g;
            lt_errno = SDSE_WS_alloc_and_init((SDSE_WS**) ws, basis);
            break;
        default:
            return LT_FAILURE;
    }

    return lt_errno;
}

void algorithm_free_ws(void* ws, Algorithm alg)
{
    switch(alg)
    {
        case ALG_BABAI:
            BABAI_WS_free(ws);
            break;
        case ALG_DPLANE:
            DP_WS_free(ws);
            break;
        case ALG_SPHERE_POHST:
        case ALG_SPHERE_DP:
            SD_WS_free(ws);
            break;
        case ALG_SPHERE_SE:
            SDSE_WS_free(ws);
            break;
        default:
            return;
    }
}


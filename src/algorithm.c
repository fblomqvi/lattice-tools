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
#include "algorithm.h"
#include "babai.h"
#include "doubleplane.h"
#include "spheredecode.h"

int algorithm_parse_name(const char* alg)
{
    if(!strcmp(alg, ALG_NAME_BABAI)) return ALG_BABAI;
    if(!strcmp(alg, ALG_NAME_DPLANE)) return ALG_DPLANE;
    if(!strcmp(alg, ALG_NAME_SPHERE)) return ALG_SPHERE;
    if(!strcmp(alg, ALG_NAME_SD_DP)) return ALG_SD_DP;
    return 0;
}

int algorithm_get_fp_init_ws(SOLVE_func* f, void** ws, 
                            Algorithm alg, const gsl_matrix* basis)
{
    switch(alg)
    {
        case ALG_BABAI:
            *f = babai_g;
            *ws = BABAI_WS_alloc_and_init(basis);
            break;
        case ALG_DPLANE:
            *f = doubleplane_g;
            *ws = DP_WS_alloc_and_init(basis);
            break;
        case ALG_SPHERE:
            *f = spheredecode_g;
            *ws = SD_WS_alloc_and_init(basis);
            break;
        case ALG_SD_DP:
            *f = sd_dp_g;
            *ws = SD_WS_alloc_and_init(basis);
            break;
        default:
            return -1;
    }

    libcheck(*ws, "*_alloc_and_init failed");
    return 0;

error:
    return -1;
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
        case ALG_SPHERE:
        case ALG_SD_DP:
            SD_WS_free(ws);
            break;
        default:
            return;
    }
}

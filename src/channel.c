/* channel.c
   Copyright (C) 2016 Ferdinand Blomqvist

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
#include "channel.h"
#include <gsl/gsl_randist.h>
#include <math.h>

struct s_channel
{
    double* received;
    void* transmitted;
    size_t n;
    size_t q;
    gsl_rng* rng;
    double sigma;
};

CHANNEL* CHANNEL_alloc(double* received, void* transmitted,
                            size_t n, size_t q, gsl_rng* rng)
{
    CHANNEL* c = malloc(sizeof(CHANNEL));
    libcheck_mem(c);

    c->received = received;
    c->transmitted = transmitted;
    c->n = n;
    c->q = q;
    c->rng = rng;
    return c;

error:
    return NULL;
}

void CHANNEL_free(CHANNEL* c)
{ free(c); }

void CHANNEL_set_sigma(CHANNEL* c, double sigma)
{ c->sigma = sigma; }


/*
void channel_get_zero_cword_BPSK_AWGN(CHANNEL* c)
{
    for(size_t i = 0; i < c->n; i++)
        c->received[i] = 1.0 + gsl_ran_gaussian(c->rng, c->sigma);
}

void channel_get_cword_BPSK_AWGN(CHANNEL* c)
{
    b_type* transmitted = (b_type*) c->transmitted;
    for(size_t i = 0; i < c->n; i++)
        c->received[i] = (transmitted[i] ? -1 : 1) + gsl_ran_gaussian(c->rng, c->sigma);
}

void channel_get_cword_BPSK_AWGN_Q(CHANNEL* c)
{
    q_type* transmitted = (q_type*) c->transmitted;
    for(size_t i = 0; i < c->n; i++)
        c->received[i] = (transmitted[i] ? -1 : 1) + gsl_ran_gaussian(c->rng, c->sigma);
}

void channel_get_zero_cword_QARY_AWGN(CHANNEL* c)
{
    for(size_t i = 0; i < c->n; i++)
    {
        fp_type temp = MOD(gsl_ran_gaussian(c->rng, c->sigma), c->q);
        if(temp < 0.0)
            temp += c->q;
        c->received[i] = temp;
    }
}

void channel_get_cword_QARY_AWGN(CHANNEL* c)
{
    q_type* transmitted = (q_type*) c->transmitted;
    for(size_t i = 0; i < c->n; i++)
    {
        fp_type temp = MOD(transmitted[i] + gsl_ran_gaussian(c->rng, c->sigma), c->q);
        if(temp < 0.0)
            temp += c->q;
        c->received[i] = temp;
    }
}
*/

void channel_get_zero_cword_AWGN(CHANNEL* c)
{
    for(size_t i = 0; i < c->n; i++)
        c->received[i] = gsl_ran_gaussian(c->rng, c->sigma);
}

void channel_get_cword_AWGN(CHANNEL* c)
{
    double* transmitted = (double*) c->transmitted;
    for(size_t i = 0; i < c->n; i++)
        c->received[i] = transmitted[i] + gsl_ran_gaussian(c->rng, c->sigma);
}

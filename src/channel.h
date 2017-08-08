/* channel.h
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

#ifndef FB_LDPC_TOOLS_CHANNEL_H
#define FB_LDPC_TOOLS_CHANNEL_H

//#include "decoder_defs.h"
#include <gsl/gsl_rng.h>

typedef struct s_channel CHANNEL;

CHANNEL* CHANNEL_alloc(double* received, void* transmitted,
                            size_t n, size_t q, gsl_rng* rng);
void CHANNEL_free(CHANNEL* ch);
void CHANNEL_set_sigma(CHANNEL* ch, double sigma);

/*
void channel_get_zero_cword_BPSK_AWGN(CHANNEL* c);
void channel_get_cword_BPSK_AWGN(CHANNEL* c);
void channel_get_cword_BPSK_AWGN_Q(CHANNEL* c);

void channel_get_zero_cword_QARY_AWGN(CHANNEL* c);
void channel_get_cword_QARY_AWGN(CHANNEL* c);
*/

void channel_get_zero_cword_AWGN(CHANNEL* c);
void channel_get_cword_AWGN(CHANNEL* c);

#endif /* FB_LDPC_TOOLS_CHANNEL_H */

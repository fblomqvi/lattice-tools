/* lat_sim.h
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

#ifndef FB_LATTICE_TOOLS_LAT_SIM_H
#define FB_LATTICE_TOOLS_LAT_SIM_H

#include "simulator.h"

typedef struct s_lat_sim_callback_args
{
    FILE* file;
} LSC_ARGS;

int lat_sim_end_callback(const SIM_STATUS* status, void* args);

int lat_sim_start_callback_std(const SIM_STATUS* status, void* args);
int lat_sim_start_callback_quiet(const SIM_STATUS* status, void* args);

int lat_sim_vnr_callback_std(const SIM_STATUS* status, void* args);
int lat_sim_vnr_callback_quiet(const SIM_STATUS* status, void* args);

#endif /* FB_LATTICE_TOOLS_LAT_SIM_H */

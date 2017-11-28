/* an_solver.h
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

#ifndef FB_LIBLATTICE_SOLVER_H
#define FB_LIBLATTICE_SOLVER_H

#include <stdlib.h>

typedef struct s_AN_workspace AN_WS;

AN_WS* AN_WS_alloc(size_t size);
void AN_WS_free(AN_WS* ws);

void AN_solve(long* clp, double* recv, AN_WS* ws);
void AN_solve_np(long* clp, double* recv, AN_WS* ws);

#endif /* FB_LIBLATTICE_SOLVER_H */

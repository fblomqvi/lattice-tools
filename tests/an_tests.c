/* an_tests.c
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

#include "minunit.h"
#include "an_solver.h"
#include <stdio.h>

static char* test_alloc()
{
    AN_WS* ws = AN_WS_alloc(5);
    mu_assert(ws, "AN_WS_alloc failed!");

    double recv[5] = {0,0,0,0,0};
    long clp[5];
    AN_solve_np(clp, recv, ws);

    for(size_t i = 0; i < 5; i++)
        printf("%ld\n", clp[i]);

    AN_WS_free(ws);

    return NULL;
}

char* all_tests()
{
    /* Running tests */
    mu_suite_start();
    mu_run_test(test_alloc);

    return NULL;
}

RUN_TESTS(all_tests);

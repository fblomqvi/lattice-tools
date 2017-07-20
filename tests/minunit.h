/* Macros for unit testing.
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

/* These unit testing macros are essentially the same as the unit testing
   macros found in chapter 30 of 'Learn C the hard way' by Zed A. Shaw
   (http://c.learncodethehardway.org/book/ex30.html). I have tweaked and
   expanded them to my liking. Furthermore I have fixed some possible
   'swallowing the semicolon' issues by adding 'do {...} * while(0)' statements
   to the right places. */
#undef NDEBUG
#ifndef FB_UTILITY_MINUNIT_H
#define FB_UTILITY_MINUNIT_H

#include <stdio.h>
#include <stdlib.h>

/* We set the PROGRAM_NAME so that debug() outputs a useful program name. */
//#define PROGRAM_NAME program_name
#include "dbg.h"

#define mu_suite_start() char* message = NULL

#define mu_assert(test, message) \
    do {\
        if(!(test))\
        {\
            log_err(message);\
            return message;\
        }\
    } while(0)

#define mu_run_test(test) \
    do {\
        debug("----- running '%s'", #test);\
        message = test();\
        tests_run++;\
        if(message)\
            return message;\
    } while(0)

#define RUN_TESTS(name) \
    int main(int argc, char* argv[])\
    {\
        tests_run = 0;\
        program_name = argv[0];\
        printf("----\nRUNNING: %s\n", program_name);\
        char* result = name();\
        if(result != 0)\
        {\
            printf("FAILED: %s\n", result);\
        }\
        else\
        {\
            printf("ALL TESTS PASSED\n");\
        }\
        printf("Tests run: %d\n", tests_run);\
        exit(result != 0);\
    }

static int tests_run;
static const char* program_name;

#endif /* FB_UTILITY_MINUNIT_H */

/* A program that solves the SVP for the given lattice.
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
#include "version.h"
#include "parse.h"
#include "sphere_se.h"
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <stdarg.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

typedef struct s_options 
{
    char* basis_file;
    char* output;
    int rows_as_basis;
} OPT;

static const OPT OPT_default = {
    .basis_file = NULL, .output = NULL, 
    .rows_as_basis = 0
};

static int print_help(FILE* file)
{
    static const char* formatstr = 
"Usage: %s [OPTION]... INPUT\n"
"  or:  %s [OPTION]... INPUT OUTPUT\n\n%s\n";

    static const char* helpstr = 

"Solves the shortest vector problem for the lattice read from INPUT. Outputs\n"
"to stdout if no output file is given.\n\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -t, --transpose              Transpose the basis read from INPUT.\n"
"      --help                   Display this help and exit.\n"
"      --version                Output version information and exit.";
    
    return (fprintf(file, formatstr, PROGRAM_NAME, PROGRAM_NAME, helpstr) < 0) 
                ? EXIT_FAILURE : EXIT_SUCCESS;
}

static void parse_cmdline(int argc, char* const argv[], OPT* opt)
{
    static const char* optstring = "t";
    static struct option longopt[] = {
        {"transpose", no_argument, NULL, 't'},
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'V'},
        {0, 0, 0, 0}
    };


    // Parsing the command line
    int ch;
    //char* endptr;
    while((ch = getopt_long(argc, argv, optstring, longopt, NULL)) != -1)
    {
        switch(ch)
        {
            case 't':
                opt->rows_as_basis = 1;
                break;
            case 'h':
                exit(print_help(stdout));
            case 'V':
                exit(print_version(stdout));
            default:
                goto error;
        }
    }

    lcheck_pf(optind < argc, log_plain, error, "No basis file given");
    if(strcmp(argv[optind], "-"))
        opt->basis_file = argv[optind];
    optind++;
    
    if(optind < argc)
        opt->output = argv[optind];

    return;

error:
    fprintf(stderr, "Try '%s --help' for more information.\n", PROGRAM_NAME);
    exit(EXIT_FAILURE);
}

static int print_lattice_point(FILE* file, const double* cword, size_t cword_len)
{
    libcheck(fprintf(file, "[") > 0, "printing failed");
    for(size_t i = 0; i < cword_len-1; i++)
        libcheck(fprintf(file, " %f", cword[i]) > 0, "printing failed");
    libcheck(fprintf(file, " %f ]\n", cword[cword_len-1]) > 0, "printing failed");
    return 0;
    
error:
    return -1;
}

static int solve_svp(FILE* outfile, OPT* opt)
{
    int ret = EXIT_FAILURE;
    debug("basis file: %s", opt->basis_file);
    FILE* infile = stdin;
    if(opt->basis_file)
    {
        infile = fopen(opt->basis_file, "r");
        check(infile, "Could not open '%s' for reading", opt->basis_file);
    }

    gsl_matrix* basis = parse_fpLLL_matrix(infile, opt->rows_as_basis);
    if(opt->basis_file)
        fclose(infile);
    lcheck(basis, error_a, "error when processing input file '%s'", opt->basis_file);

    SDSE_WS* ws = SDSE_WS_alloc_and_init(basis);
    lcheck_mem(ws, error_a);

    double* slp = malloc(basis->size1 * sizeof(double));
    lcheck_mem(slp, error_b);

    gsl_vector_view v_slp = gsl_vector_view_array(slp, basis->size1);
    sphere_se_svp(&v_slp.vector, basis, ws);

    int rc = print_lattice_point(outfile, slp, basis->size1);
    lcheck(rc == 0, error_c, "print_lattice_point failed");
    ret = EXIT_SUCCESS;

error_c:
    free(slp);
error_b:
    SDSE_WS_free(ws);
error_a:
    gsl_matrix_free(basis);
error:
    return ret;
}

int main(int argc, char* argv[])
{
    OPT opt = OPT_default;
    int ret = EXIT_FAILURE;
    argv[0] = PROGRAM_NAME = "lat-svp";
    parse_cmdline(argc, argv, &opt);

    FILE* outfile = stdout;
    if(opt.output)
    {
        outfile = fopen(opt.output, "w");
        check(outfile, "Could not open '%s' for writing", opt.output);
    }

    int rc = solve_svp(outfile, &opt);
    llibcheck(rc == 0, error_a, "solve_svp failed");
    ret = EXIT_SUCCESS;

error_a:
    if(opt.output)
        fclose(outfile);
error:
    return ret;
}

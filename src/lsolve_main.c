/* A program that solves the CVP on the given lattice.
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
#include "rng.h"
#include "rnd_point.h"
#include "algorithm.h"
#include "utility.h"
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <stdarg.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 10E-10

typedef enum enum_mode
{
    MODE_STANDARD = 1,
    MODE_COMPARE
} Mode;

typedef struct s_options 
{
    char* basis_file;
    char* output;
    size_t cword_num;
    Algorithm alg;
    Algorithm alg_cmp;
    Mode mode;
    int no_config;
    int binary_out;
    int rows_as_basis;
    gsl_matrix* basis;
    SOLVE_func solve;
    SOLVE_func solve_cmp;
    void* ws;
    void* ws_cmp;
} OPT;

static const OPT OPT_default = {
    .basis_file = NULL, .output = NULL, 
    .cword_num = 0, .rows_as_basis = 0,
    .no_config = 0, .binary_out = 1,
    .alg = ALG_SPHERE_SE, 
    .mode = MODE_STANDARD, .basis = NULL,
    .ws = NULL, .ws_cmp = NULL };

static int print_help(FILE* file)
{
    static const char* formatstr = 
"Usage: %s [OPTION]... INPUT\n"
"  or:  %s [OPTION]... INPUT OUTPUT\n\n%s\n";

    static const char* helpstr = 

"Solves the closest vector problem on (R^n, L), where L is the lattice defined\n"
"by the basis read from INPUT. Reads the points to decode from stdin in the\n"
"binary format produced by rnd-point. Outputs to stdout if no output file is given.\n\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -a, --algorithm=ALG1         Select the decoding algorithm. To see a list of all\n"
"                                 available algorithms give 'list' as argument.\n"
"                                 The default algorithm is '" ALG_NAME_SPHERE "'.\n"
"  -c, --compare=ALG2           Compare ALG2 to ALG1.\n"
"  -n, --num-points=NUM         The number of codewords to decode. Zero (0) makes the\n"
"                                 solver run until it runs out of input. NOTE: this\n"
"                                 option is not yet implemented, it always behaves as\n"
"                                 zero had been given.\n"
"  -C, --no-config              Do not try to read the configuration from stdin.\n"
"  -R, --readable-output        Produce readable output instead of binary output.\n"
"  -t, --transpose              Transpose the basis read from INPUT.\n"
"      --help                   Display this help and exit.\n"
"      --version                Output version information and exit.";
    
    return (fprintf(file, formatstr, PROGRAM_NAME, PROGRAM_NAME, helpstr) < 0) 
                ? EXIT_FAILURE : EXIT_SUCCESS;
}

static void parse_cmdline(int argc, char* const argv[], OPT* opt)
{
    static const char* optstring = "a:c:n:tRC";
    static struct option longopt[] = {
        {"algorithm", required_argument, NULL, 'a'},
        {"compare", required_argument, NULL, 'c'},
        {"num-points", required_argument, NULL, 'n'},
        {"readable-output", no_argument, NULL, 'R'},
        {"transpose", no_argument, NULL, 't'},
        {"no-config", no_argument, NULL, 'C'},
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'V'},
        {0, 0, 0, 0}
    };


    // Parsing the command line
    int ch;
    char* endptr;
    while((ch = getopt_long(argc, argv, optstring, longopt, NULL)) != -1)
    {
        switch(ch)
        {
            case 'a':
                if(!strcmp(optarg, "list"))
                    exit(algorithm_print_names(stdout));
                else
                {
                    opt->alg = algorithm_parse_name(optarg);
                    check(opt->alg > 0, "invalid argument to option '%c': '%s'",
                            ch, optarg);
                }
                break;
            case 'c':
                opt->alg_cmp = algorithm_parse_name(optarg);
                check(opt->alg_cmp > 0, "invalid argument to option '%c': '%s'", ch, optarg);
                opt->mode = MODE_COMPARE;
                break;
            case 'n':
                opt->cword_num = strtoul(optarg, &endptr, 10);
                check(*endptr == '\0' && !(errno == ERANGE && opt->cword_num == ULONG_MAX),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 'C':
                opt->no_config = 1;
                break;
            case 'R':
                opt->binary_out = 0;
                break;
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
    opt->basis_file = argv[optind++];
    
    if(optind < argc)
        opt->output = argv[optind];

    return;

error:
    fprintf(stderr, "Try '%s --help' for more information.\n", PROGRAM_NAME);
    exit(EXIT_FAILURE);
}

static int print_codeword(FILE* file, const double* cword, size_t cword_len)
{
    libcheck(fprintf(file, "(") > 0, "printing failed");
    for(size_t i = 0; i < cword_len-1; i++)
        libcheck(fprintf(file, "%f, ", cword[i]) > 0, "printing failed");
    libcheck(fprintf(file, "%f) ---> ", cword[cword_len-1]) > 0, "printing failed");
    return 0;
    
error:
    return -1;
}

static int print_lattice_point(FILE* file, const double* cword, size_t cword_len)
{
    libcheck(fprintf(file, "(") > 0, "printing failed");
    for(size_t i = 0; i < cword_len-1; i++)
        libcheck(fprintf(file, "%f, ", cword[i]) > 0, "printing failed");
    libcheck(fprintf(file, "%f)\n", cword[cword_len-1]) > 0, "printing failed");
    return 0;
    
error:
    return -1;
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

static RND_PNT_CONF* get_config_binary(FILE* infile, const OPT* opt)
{
    RND_PNT_CONF* conf = opt->no_config ? RND_PNT_CONF_alloc() : RND_PNT_CONF_read(infile);
    libcheck(conf, "RND_PNT_CONF_* failed");

    if(opt->no_config)
    {
        conf->dimension = opt->basis->size1;
        conf->num_cwords = opt->cword_num;
    }
    else
    {
        if(opt->cword_num != OPT_default.cword_num)
            conf->num_cwords = opt->cword_num;
    }
    debug("dimension: %zu", conf->dimension);
    debug("num_cwords: %zu", conf->num_cwords);

    return conf;

error:
    return NULL;
}

static int parse_basis_init_ws(OPT* opt)
{
    debug("basis file: %s", opt->basis_file);
    FILE* infile = fopen(opt->basis_file, "r");
    check(infile, "Could not open '%s' for reading", opt->basis_file);

    // TODO: Add checks
    opt->basis = read_matrix(infile, opt->rows_as_basis);
    fclose(infile);

    algorithm_get_fp_init_ws(&opt->solve, &opt->ws, opt->alg, opt->basis);
    if(opt->mode == MODE_COMPARE)
        algorithm_get_fp_init_ws(&opt->solve_cmp, &opt->ws_cmp, opt->alg_cmp, opt->basis);

    //print_matrix(opt->basis);

    return 0;

error:
    return -1;
}

static void free_basis_and_ws(OPT* opt)
{
    gsl_matrix_free(opt->basis);
    algorithm_free_ws(opt->ws, opt->alg);
    if(opt->mode == MODE_COMPARE)
        algorithm_free_ws(opt->ws_cmp, opt->alg_cmp);
}

static int solve(FILE* outfile, OPT* opt)
{
    int ret = EXIT_FAILURE;
    int rc = parse_basis_init_ws(opt);
    libcheck(rc == 0, "parse_basis_init_ws failed");

    RND_PNT_CONF* conf = get_config_binary(stdin, opt);
    llibcheck(conf, error_a, "get_config_binary failed");

    double* cword = malloc(2 * conf->dimension * sizeof(double));
    llibcheck_mem(cword, error_b);

    double* clp = cword + conf->dimension;
    gsl_vector_view v_cword = gsl_vector_view_array(cword, conf->dimension);
    gsl_vector_view v_clp = gsl_vector_view_array(clp, conf->dimension);

    if(opt->binary_out)
    {
        while(1)
        {
            int rc = read_cword_binary(cword, conf->dimension, sizeof(double));
            llibcheck(rc == 0, error_c, "read_cword_binary failed");
            
            opt->solve(&v_clp.vector, &v_cword.vector, opt->basis, opt->ws);
            llibcheck(fwrite(clp, sizeof(long), conf->dimension, outfile) 
                    == conf->dimension, error_c, "fwrite failed");
        }
    }
    else
    {
        while(1)
        {
            int rc = read_cword_binary(cword, conf->dimension, sizeof(double));
            llibcheck(rc == 0, error_c, "read_cword_binary failed");
            
            rc = print_codeword(outfile, cword, conf->dimension);
            llibcheck(rc == 0, error_c, "print_codeword failed");

            opt->solve(&v_clp.vector, &v_cword.vector, opt->basis, opt->ws);
            rc = print_lattice_point(outfile, clp, conf->dimension);
            llibcheck(rc == 0, error_c, "print_lattice_point failed");
        }
    }

    ret = EXIT_SUCCESS;

error_c:
    free_basis_and_ws(opt);
error_b:
    free(cword);
error_a:
    RND_PNT_CONF_free(conf);
error:
    return ret;
}

static int solutions_not_equal(double* a, double* b, size_t len)
{
    for(size_t i = 0; i < len; i++)
        if(fabs(a[i] - b[i]) > EPSILON)
            return 1;

    return 0;
}

static int compare(FILE* outfile, OPT* opt)
{
    int ret = EXIT_FAILURE;
    int rc = parse_basis_init_ws(opt);
    libcheck(rc == 0, "parse_basis_init_ws failed");

    RND_PNT_CONF* conf = get_config_binary(stdin, opt);
    llibcheck(conf, error_a, "get_config_binary failed");

    double* cword = malloc(3 * conf->dimension * sizeof(double));
    llibcheck_mem(cword, error_b);

    double* clp1 = cword + conf->dimension;
    double* clp2 = clp1 + conf->dimension;
    gsl_vector_view v_cword = gsl_vector_view_array(cword, conf->dimension);
    gsl_vector_view v_clp1 = gsl_vector_view_array(clp1, conf->dimension);
    gsl_vector_view v_clp2 = gsl_vector_view_array(clp2, conf->dimension);

    size_t num_checked = 0;
    size_t num_different = 0;
    while(1)
    {
        int rc = read_cword_binary(cword, conf->dimension, sizeof(double));
        llibcheck(rc == 0, error_c, "read_cword_binary failed");

        opt->solve(&v_clp1.vector, &v_cword.vector, opt->basis, opt->ws);
        opt->solve_cmp(&v_clp2.vector, &v_cword.vector, opt->basis, opt->ws_cmp);
        if(solutions_not_equal(clp1, clp2, conf->dimension))
            num_different++;
        num_checked++;
    }

    ret = EXIT_SUCCESS;

error_c:
    fprintf(outfile, "Compared %zu solutions and found %zu differences\n",
            num_checked, num_different);
    free_basis_and_ws(opt);
error_b:
    free(cword);
error_a:
    RND_PNT_CONF_free(conf);
error:
    return ret;
}

int main(int argc, char* argv[])
{
    OPT opt = OPT_default;
    int ret = EXIT_FAILURE;
    argv[0] = PROGRAM_NAME = "lsolve";
    parse_cmdline(argc, argv, &opt);

    FILE* outfile = stdout;
    if(opt.output)
    {
        outfile = fopen(opt.output, "w");
        check(outfile, "Could not open '%s' for writing", opt.output);
    }

    switch(opt.mode)
    {
        case MODE_STANDARD:
        {
            int rc = solve(outfile, &opt);
            llibcheck(rc == 0, error_a, "solve failed");
            break;
        }
        case MODE_COMPARE:
        {
            int rc = compare(outfile, &opt);
            llibcheck(rc == 0, error_a, "compare failed");
            break;
        }
    }
    ret = EXIT_SUCCESS;

error_a:
    if(opt.output)
        fclose(outfile);
error:
    return ret;
}

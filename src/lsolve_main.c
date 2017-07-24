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
#include "babai.h"
#include "spheredecode.h"
#include "doubleplane.h"
#include "utility.h"
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <stdarg.h>
#include <time.h>
#include <stdlib.h>

#define ALG_NAME_BABAI "babai"
#define ALG_NAME_DPLANE "dplane"
#define ALG_NAME_SPHERE "sphere"

typedef enum algorithm
{
    ALG_BABAI = 1,
    ALG_DPLANE,
    ALG_SPHERE,
} Algorithm;

typedef struct s_options 
{
    char* basis_file;
    char* output;
    size_t cword_num;
    Algorithm alg;
    int no_config;
    int binary_out;
    void (*solve)(gsl_vector*, const gsl_vector*, const gsl_matrix*, void*);
    gsl_matrix* basis;
    void* ws;
} OPT;

static const OPT OPT_default = {
    .basis_file = NULL, .output = NULL, 
    .cword_num = 0, 
    .no_config = 0, .binary_out = 1, 
    .alg = ALG_SPHERE,
    .basis = NULL, .ws = NULL };

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
"  -a, --algorithm              Select the decoding algorithm. Valid values are\n"
"                                 '" ALG_NAME_BABAI "', '" ALG_NAME_DPLANE "' and '" 
                                ALG_NAME_SPHERE "'. The default is \n"
"                                 '" ALG_NAME_SPHERE "'.\n"
"  -n, --num-points=NUM         The number of codewords to decode. Zero (0) makes the\n"
"                                 solver run until it runs out of input. NOTE: this\n"
"                                 option is not yet implemented, it always behaves as\n"
"                                 zero had been given.\n"
"  -C, --no-config              Do not try to read the configuration from stdin.\n"
"  -R, --readable-output        Produce readable output instead of binary output.\n"
"      --help                   Display this help and exit.\n"
"      --version                Output version information and exit.";
    
    return (fprintf(file, formatstr, PROGRAM_NAME, PROGRAM_NAME, helpstr) < 0) 
                ? EXIT_FAILURE : EXIT_SUCCESS;
}

static int parse_alg(const char* alg)
{
    if(!strcmp(alg, ALG_NAME_BABAI)) return ALG_BABAI;
    if(!strcmp(alg, ALG_NAME_DPLANE)) return ALG_DPLANE;
    if(!strcmp(alg, ALG_NAME_SPHERE)) return ALG_SPHERE;
    return 0;
}

static void parse_cmdline(int argc, char* const argv[], OPT* opt)
{
    static const char* optstring = "a:n:RC";
    static struct option longopt[] = {
        {"algorithm", required_argument, NULL, 'a'},
        {"num-points", required_argument, NULL, 'n'},
        {"readable-output", no_argument, NULL, 'R'},
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
                opt->alg = parse_alg(optarg); 
                check(opt->alg > 0, "invalid argument to option '%c': '%s'", ch, optarg);
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
    opt->basis = read_matrix(infile);
    fclose(infile);

    //print_matrix(basis);
    switch(opt->alg)
    {
        case ALG_BABAI:
            opt->solve = babai_g;
            opt->ws = BABAI_WS_alloc_and_init(opt->basis);
            break;
        case ALG_DPLANE:
            opt->solve = doubleplane_g;
            opt->ws = DP_WS_alloc_and_init(opt->basis);
            break;
        case ALG_SPHERE:
            opt->solve = spheredecode_g;
            opt->ws = SD_WS_alloc_and_init(opt->basis);
            break;
        default:
            return -1;
    }

    return 0;

error:
    return -1;
}

static void free_basis_and_ws(OPT* opt)
{
    gsl_matrix_free(opt->basis);
    switch(opt->alg)
    {
        case ALG_BABAI:
            BABAI_WS_free(opt->ws);
            break;
        case ALG_DPLANE:
            DP_WS_free(opt->ws);
            break;
        case ALG_SPHERE:
            SD_WS_free(opt->ws);
            break;
        default:
            return;
    }
}

static int solve(FILE* outfile, OPT* opt)
{
    int ret = EXIT_FAILURE;
    int rc = parse_basis_init_ws(opt);
    libcheck(rc == 0, "parse_basis_init_ws failed");

    RND_PNT_CONF* conf = get_config_binary(stdin, opt);
    llibcheck(conf, error_a, "get_config_binary failed");

    gsl_vector* cword = gsl_vector_alloc(conf->dimension);
    llibcheck_mem(cword, error_b);

    gsl_vector* clp = gsl_vector_alloc(conf->dimension);
    llibcheck_mem(clp, error_c);

    if(opt->binary_out)
    {
        while(1)
        {
            // A dirty hack to comply with existing code
            int rc = read_cword_binary(cword->data, conf->dimension, sizeof(double));
            llibcheck(rc == 0, error_d, "read_cword_binary failed");
            
            opt->solve(clp, cword, opt->basis, opt->ws);
            // A dirty hack to comply with existing code
            llibcheck(fwrite(clp->data, sizeof(long), conf->dimension, outfile) 
                    == conf->dimension, error_d, "fwrite failed");
        }
    }
    else
    {
        while(1)
        {
            // A dirty hack to comply with existing code
            int rc = read_cword_binary(cword->data, conf->dimension, sizeof(double));
            llibcheck(rc == 0, error_d, "read_cword_binary failed");
            
            // A dirty hack to comply with existing code
            rc = print_codeword(outfile, cword->data, conf->dimension);
            llibcheck(rc == 0, error_d, "print_codeword failed");

            opt->solve(clp, cword, opt->basis, opt->ws);
            // A dirty hack to comply with existing code
            rc = print_lattice_point(outfile, clp->data, conf->dimension);
            llibcheck(rc == 0, error_d, "print_lattice_point failed");
        }
    }

    ret = EXIT_SUCCESS;

error_d:
    free_basis_and_ws(opt);
error_c:
    gsl_vector_free(clp);
error_b:
    gsl_vector_free(cword);
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

    int rc = solve(outfile, &opt);
    llibcheck(rc == 0, error_a, "generate_* failed");
    ret = EXIT_SUCCESS;

error_a:
    if(opt.output)
        fclose(outfile);
error:
    return ret;
}

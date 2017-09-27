/* A program that generates lattices.
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
#include "lattice_gen.h"
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <stdarg.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#define LNAME_TYPE_A "A"
#define LNAME_TYPE_D "D"
#define LNAME_TYPE_D_DUAL "D-dual"
#define LNAME_TYPE_RANDOM "R"
#define LNAME_TYPE_SPC "SPC"

typedef struct s_options 
{
    char* output;
    const gsl_rng_type* rng_type;
    size_t dimension;
    size_t exponent;
    size_t q;
    L_type type;
    unsigned long seed;
    int cols_as_basis;
} OPT;

static const OPT OPT_default = {
    .output = NULL, .exponent = 1,
    .dimension = 0, .cols_as_basis = 0,
    .seed = 0, .type = LTYPE_RANDOM, .q = 0
};

static int gen_and_print(FILE* file, OPT* opt);
static void parse_cmdline(int argc, char* const argv[], OPT* opt);
static int print_help(FILE* file);

static int type_parse_name(const char* type);
static int type_print_names(FILE* file);


int main(int argc, char* argv[])
{
    OPT opt = OPT_default;
    int ret = EXIT_FAILURE;
    argv[0] = PROGRAM_NAME = "lgen";
    parse_cmdline(argc, argv, &opt);

    FILE* outfile = stdout;
    if(opt.output)
    {
        outfile = fopen(opt.output, "w");
        check(outfile, "Could not open '%s' for writing", opt.output);
    }

    int rc = gen_and_print(outfile, &opt);
    llibcheck(rc == 0, error_a, "gen_and_print failed");
    ret = EXIT_SUCCESS;

error_a:
    if(opt.output)
        fclose(outfile);
error:
    return ret;
}

static int gen_and_print(FILE* file, OPT* opt)
{
    LGEN_PARAMS params = { 
        .dimension = opt->dimension, .exponent = opt->exponent,
        .max = 0, .rng = NULL
    };

    if(lattice_gen_type_needs_rng(opt->type))
    {
        opt->seed = opt->seed ? opt->seed : (unsigned long) time(NULL) + clock();
        params.rng = rng_alloc_and_seed(opt->rng_type, opt->seed);
        libcheck_mem(params.rng);
    }

    t_MAT_MPZ* M = lattice_gen(opt->type, &params);
    lcheck(M, error_a, "lattice_gen failed");
    MAT_MPZ_print_fpLLL(file, M, opt->cols_as_basis);
    MAT_MPZ_free(M);
    gsl_rng_free(params.rng);

    return 0;

error_a:
    gsl_rng_free(params.rng);
error:
    return -1;
}

static void parse_cmdline(int argc, char* const argv[], OPT* opt)
{
    static const char* optstring = "d:e:r:S:tT:";
    static struct option longopt[] = {
        {"dimension", required_argument, NULL, 'd'},
        {"exponent", required_argument, NULL, 'e'},
        {"rng", required_argument, NULL, 'r'},
        {"seed", required_argument, NULL, 'S'},
        {"transpose", no_argument, NULL, 't'},
        {"type", required_argument, NULL, 'T'},
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'V'},
        {0, 0, 0, 0}
    };

    opt->rng_type = gsl_rng_default;

    // Parsing the command line
    int ch;
    char* endptr;
    while((ch = getopt_long(argc, argv, optstring, longopt, NULL)) != -1)
    {
        switch(ch)
        {
            case 'd':
                opt->dimension = strtoul(optarg, &endptr, 10);
                check(*endptr == '\0' && !(errno == ERANGE && opt->dimension == ULONG_MAX),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 'e':
                opt->exponent = strtoul(optarg, &endptr, 10);
                check(*endptr == '\0' && !(errno == ERANGE && opt->exponent == ULONG_MAX),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 'r':
            {
                const gsl_rng_type** rng_types = gsl_rng_types_setup();
                if(!strcmp(optarg, "list"))
                    exit(print_rngs(stdout, rng_types));
                else
                {
                    opt->rng_type = get_rng_type(optarg, rng_types);
                    check(opt->rng_type, "invalid random number generator");
                }
                break;
            }
            case 'S':
            {
                opt->seed = strtoul(optarg, &endptr, 10);
                check(*endptr == '\0' && !(errno == ERANGE && opt->seed == ULONG_MAX), 
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            }
            case 't':
                opt->cols_as_basis = 1;
                break;
            case 'T':
                if(!strcmp(optarg, "list"))
                    exit(type_print_names(stdout));
                else
                {
                    opt->type = type_parse_name(optarg);
                    check(opt->type > 0, "invalid argument to option '%c': '%s'",
                            ch, optarg);
                }
                break;
            case 'h':
                exit(print_help(stdout));
            case 'V':
                exit(print_version(stdout));
            default:
                goto error;
        }
    }

    if(optind < argc)
        opt->output = argv[optind];

    return;

error:
    fprintf(stderr, "Try '%s --help' for more information.\n", PROGRAM_NAME);
    exit(EXIT_FAILURE);
}

static int print_help(FILE* file)
{
    static const char* formatstr = 
"Usage: %s [OPTION]...\n"
"  or:  %s [OPTION]... OUTPUT\n\n%s\n";

    static const char* helpstr = 
"Generates lattices of the given type. Outputs to stdout if no output file is given.\n\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -d, --dimension=DIM          The dimension of the lattice to generate.\n"
"  -e, --exponent=VAL           The value of the exponent when creating a lattice\n"
"                                 of type A_d^e.\n"
"  -r, --rng=RNG                The random number generator to use. To see a list of all\n"
"                                 available generators give 'list' as argument.\n"
"  -S, --seed=SEED              The seed for the random number generator.\n"
"  -t, --transpose              Transpose the output.\n"
"  -T, --type=TYPE              The type of lattice to generate.\n"
"      --help                   Display this help and exit.\n"
"      --version                Output version information and exit.";
    
    return (fprintf(file, formatstr, PROGRAM_NAME, PROGRAM_NAME, helpstr) < 0) 
                ? EXIT_FAILURE : EXIT_SUCCESS;
}

static int type_parse_name(const char* type)
{
    if(!strcmp(type, LNAME_TYPE_A)) return LTYPE_A;
    if(!strcmp(type, LNAME_TYPE_D)) return LTYPE_D;
    if(!strcmp(type, LNAME_TYPE_D_DUAL)) return LTYPE_D_DUAL;
    if(!strcmp(type, LNAME_TYPE_RANDOM)) return LTYPE_RANDOM;
    if(!strcmp(type, LNAME_TYPE_SPC)) return LTYPE_SPC;
    return 0;
}

static int type_print_names(FILE* file)
{
    libcheck(fprintf(file, "Available lattice types are:\n") > 0, "printing error");
    libcheck(fprintf(file, "%s\n", LNAME_TYPE_A) > 0, "printing error");
    libcheck(fprintf(file, "%s\n", LNAME_TYPE_D) > 0, "printing error");
    libcheck(fprintf(file, "%s\n", LNAME_TYPE_D_DUAL) > 0, "printing error");
    libcheck(fprintf(file, "%s\n", LNAME_TYPE_RANDOM) > 0, "printing error");
    libcheck(fprintf(file, "%s\n", LNAME_TYPE_SPC) > 0, "printing error");
    return 0;

error:
    return -1;
}

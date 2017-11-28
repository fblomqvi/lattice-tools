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
#include "configuration.h"
#include "rng.h"
#include "lattice_gen.h"
#include "print_util.h"
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <stdarg.h>
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
    LGEN_PARAMS par;
    L_type type;
    unsigned long seed;
    size_t bits;
    long min;
    long max;
    int min_set;
    int max_set;
    int cols_as_basis;
    int no_config;
    int use_gmp_rand;
    enum PrintingFmt format;
} OPT;

static const OPT OPT_default = {
    .output = NULL, .cols_as_basis = 0,
    .seed = 0, .type = LTYPE_RANDOM,
    .no_config = 0, .bits = 0,
    .use_gmp_rand = 0,
    .min = 0, .max = 0,
    .format = PRINTING_FMT_DEFAULT,
    .par = {
        .dimension = 0, .exponent = 1,
        .offset = 0, .range = 0, .rng = NULL
    }
};

static int gen_and_print(FILE* file, OPT* opt);
static void parse_cmdline(int argc, char* const argv[], OPT* opt);
static int print_help(FILE* file);
static int print_config(FILE* file, const OPT* opt, gsl_rng* rng);
static int set_rng_range_and_offset(OPT* opt);
static const char* type_get_name(L_type type);
static int type_parse_name(const char* type);
static int type_print_names(FILE* file);


int main(int argc, char* argv[])
{
    OPT opt = OPT_default;
    int ret = EXIT_FAILURE;
    argv[0] = PROGRAM_NAME = "lat-gen";
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
    if(lattice_gen_type_needs_rng(opt->type))
    {
        // Check if we should use gmp for randomness.
        opt->seed = opt->seed ? opt->seed : get_random_seed();

        if(!opt->use_gmp_rand)
        {
            opt->par.rng = rng_alloc_and_seed(opt->rng_type, opt->seed);
            libcheck_mem(opt->par.rng);

            int rc = set_rng_range_and_offset(opt);
            libcheck(rc == 0, "set_rng_range_and_offset failed");
        }
    }
    else
    {
        // Ignore the min and max options when printing the config.
        opt->bits = opt->min_set = opt->max_set = 0;
    }

    if(!opt->no_config)
    {
        int rc = print_config(file, opt, opt->par.rng);
        check(rc == 0, "print_config failed");
    }


    t_MAT_MPZ* M;
    if(opt->use_gmp_rand)
        M = lattice_gen_random_square_gmp(opt->par.dimension, opt->par.exponent,
                                        opt->seed, opt->bits);
    else
        M = lattice_gen(opt->type, &opt->par);
    check(M, "lattice_gen failed");

    MAT_MPZ_print(file, M, opt->cols_as_basis, printing_fmt_get(opt->format));
    MAT_MPZ_free(M);
    gsl_rng_free(opt->par.rng);

    return 0;

error:
    gsl_rng_free(opt->par.rng);
    return -1;
}

static void parse_cmdline(int argc, char* const argv[], OPT* opt)
{
    static const char* optstring = "b:Cd:e:f:Gm:M:r:S:tT:";
    static struct option longopt[] = {
        {"bits", required_argument, NULL, 'b'},
        {"output-format", required_argument, NULL, 'f'},
        {"no-print-config", no_argument, NULL, 'C'},
        {"dimension", required_argument, NULL, 'd'},
        {"exponent", required_argument, NULL, 'e'},
        {"gmp-for-rand", no_argument, NULL, 'G'},
        {"min", required_argument, NULL, 'm'},
        {"max", required_argument, NULL, 'M'},
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
            case 'b':
                opt->bits = strtoul(optarg, &endptr, 10);
                check(*endptr == '\0' && !(errno == ERANGE && opt->bits == ULONG_MAX)
                        && opt->bits >= 1,
                    "invalid argument to option '%c': '%s'", ch, optarg);
                opt->min_set = opt->max_set = 0;
                break;
            case 'C':
                opt->no_config = 1;
                break;
            case 'd':
                opt->par.dimension = strtoul(optarg, &endptr, 10);
                check(*endptr == '\0' &&
                        !(errno == ERANGE && opt->par.dimension == ULONG_MAX),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 'e':
                opt->par.exponent = strtoul(optarg, &endptr, 10);
                check(*endptr == '\0' &&
                        !(errno == ERANGE && opt->par.exponent == ULONG_MAX),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 'f':
                if(!strcmp(optarg, "list"))
                    exit(printing_fmt_print_names(stdout));
                else
                {
                    int format_id = printing_fmt_parse_name(optarg);
                    check(format_id >= 0, "invalid argument to option '%c': '%s'",
                            ch, optarg);
                    opt->format = format_id;
                }
                break;
            case 'G':
                opt->use_gmp_rand = 1;
                break;
            case 'm':
                opt->min = strtol(optarg, &endptr, 10);
                check(*endptr == '\0' &&
                    !(errno == ERANGE && (opt->min == LONG_MAX || opt->min == LONG_MIN)),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                opt->min_set = 1;
                opt->bits = 0;
                break;
            case 'M':
                opt->max = strtol(optarg, &endptr, 10);
                check(*endptr == '\0' &&
                    !(errno == ERANGE && (opt->max == LONG_MAX || opt->max == LONG_MIN)),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                opt->max_set = 1;
                opt->bits = 0;
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

    // Check for mandatory arguments.
    check(opt->par.dimension > 0, "missing mandatory option -- '%c'", 'd');

    check(opt->type != LTYPE_SPC || opt->max_set,
            "option '%c' is mandatory for type %s", 'M', type_get_name(LTYPE_SPC));

    opt->use_gmp_rand = opt->type == LTYPE_RANDOM && opt->bits
                        && (opt->use_gmp_rand || opt->bits >= 32);

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
"Generates integer lattices of the given type. Outputs a basis matrix (with the rows\n"
"as basis vectors) for the generated lattice. The outputs is written to stdout if no\n"
"output file is given.\n\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -b, --bits=BITS              Generates matrix entries with absolute value smaller\n"
"                                 than 2^BITS.\n"
"  -C, --no-print-config        Do not output the configuration the program was run with.\n"
"  -d, --dimension=DIM          The dimension of the lattice to generate. This option\n"
"                                 is mandatory.\n"
"  -e, --exponent=VAL           Output the VAL:th power of the basis matrix.\n"
"  -f, --output-format=FMT      The output format that will be used. Give 'list' as\n"
"                                 argument to get a list of all formats.\n"
"  -G, --gmp-for-rand           Use the gmp library for randomness. This is the default\n"
"                                 if BITS > 31. Therefore the purpose of this option is\n"
"                                 to force the use of gmp as the source of randomness\n"
"                                 when BITS <= 31. This option is ignored unless the\n"
"                                 lattice type is 'R' and the range of the matrix entries\n"
"                                 is specified with the 'b' option. However, whenever\n"
"                                 this option is eligible, it takes precedence over the\n"
"                                 'r' option.\n"
"  -m, --min=NUM                The minimum value for the matrix entries (inclusive).\n"
"                                 Defaults to the min-value of the chosen rng.\n"
"  -M, --max=NUM                The maximum value for the matrix entries (exclusive).\n"
"                                 Defaults to the max-value of the chosen rng.\n"
"  -r, --rng=RNG                The random number generator to use. To see a list of all\n"
"                                 available generators give 'list' as argument.\n"
"  -S, --seed=SEED              The seed for the random number generator.\n"
"  -t, --transpose              Transpose the output.\n"
"  -T, --type=TYPE              The type of lattice to generate. Give 'list' as argument\n"
"                                 to see a list of all available lattice types.\n"
"      --help                   Display this help and exit.\n"
"      --version                Output version information and exit.\n\n"
"The options 'b' and 'm' combined with 'M' are alternative ways of specifying\n"
"the range of the matrix entries. If the user gives options from both categories,\n"
"then the one specified last takes precedence.\n\n"
"See the man page for further details.";

    return (fprintf(file, formatstr, PROGRAM_NAME, PROGRAM_NAME, helpstr) < 0)
                ? EXIT_FAILURE : EXIT_SUCCESS;
}

static int print_config(FILE* file, const OPT* opt, gsl_rng* rng)
{
    struct config conf[] = {
        {"bits", (union value) opt->bits, type_size, opt->bits},
        {"dimension", (union value) opt->par.dimension, type_size, 1},
        {"exponent", (union value) opt->par.exponent, type_size, opt->par.exponent > 1},
        {"output-format", (union value) printing_fmt_get_name(opt->format),
                        type_str, opt->format != PRINTING_FMT_DEFAULT},
        {"gmp-for-rand", (union value) opt->use_gmp_rand, type_bool, opt->bits < 32},
        {"min", (union value) opt->min, type_long, opt->min_set},
        {"max", (union value) opt->max, type_long, opt->max_set},
        {"rng", (union value) (rng ? gsl_rng_name(rng) : ""), type_str,
                        rng != NULL && rng->type != gsl_rng_default},
        {"seed", (union value) opt->seed, type_ulong, opt->seed != 0},
        {"transpose", (union value) opt->cols_as_basis, type_bool, 1},
        {"type", (union value) type_get_name(opt->type), type_str, 1},
        {0, (union value) 0, type_bool, 0}
    };

    return util_print_genby_and_config(file, PROGRAM_NAME, NULL, NULL, conf);
}

static int set_rng_range_and_offset(OPT* opt)
{
    if(opt->type == LTYPE_SPC)
    {
        // Ignore the min option.
        opt->par.offset = 0;
        opt->min_set = 0;
        check(opt->max > 0, "the maximum value (-M) must be positive");
        opt->par.range = opt->max;
    }
    else
    {
        if(opt->bits)
        {
            debug("opt->bits is set");
            unsigned long val= 1L << opt->bits;
            opt->par.offset = -val + 1;
            opt->par.range = val - opt->par.offset;
        }
        else
        {
            if(opt->min_set || opt->max_set)
            {
                check(opt->min < opt->max, "The minimum value (-m) must be "
                        "smaller than the maximu value (-M)");
            }
            opt->par.offset = opt->min_set ?
                opt->min : (long) gsl_rng_min(opt->par.rng);
            opt->par.range = (opt->max_set ? opt->max
                    : (long) gsl_rng_max(opt->par.rng)) - opt->par.offset;
        }
    }
    debug("min: %ld", opt->par.offset);
    debug("range: %lu", opt->par.range);
    debug("max: %ld", opt->par.offset + opt->par.range);
    debug("rng_min: %lu", gsl_rng_min(opt->par.rng));
    debug("rng_max: %lu", gsl_rng_max(opt->par.rng));
    check(opt->par.range <= gsl_rng_max(opt->par.rng),
            "The selected range (min: %ld, max: %lu) is to large for the selected rng",
            opt->par.offset, opt->par.offset + opt->par.range);
    return 0;

error:
    return -1;
}

static const char* type_get_name(L_type type)
{
    switch(type)
    {
        case LTYPE_A: return LNAME_TYPE_A;
        case LTYPE_D: return LNAME_TYPE_D;
        case LTYPE_D_DUAL: return LNAME_TYPE_D_DUAL;
        case LTYPE_RANDOM: return LNAME_TYPE_RANDOM;
        case LTYPE_SPC: return LNAME_TYPE_SPC;
        default: return NULL;
    }
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
    libcheck(fprintf(file, "%s (default)\n", LNAME_TYPE_RANDOM) > 0, "printing error");
    libcheck(fprintf(file, "%s\n", LNAME_TYPE_SPC) > 0, "printing error");
    return 0;

error:
    return -1;
}

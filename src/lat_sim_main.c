/* A program that runs decoding simulations with lattices.
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
#include "parse.h"
#include "simulator.h"
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <stdarg.h>
#include <math.h>

typedef struct s_opt
{
    int transpose;
    Algorithm alg;
    SIM_OPTIONS sim;
} OPT;


static int print_help(FILE* file)
{
    static const char* formatstr = 
"Usage: %s [OPTION]... INPUT OUTPUT\n\n%s\n";

    static const char* helpstr = 
"Decoding simulations for lattices. The basis of the lattice is read from INPUT.\n\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -a, --algorithm=ALG          Select the decoding algorithm. To see a list of all\n"
"                                 available algorithms give 'list' as argument.\n"
"                                 The default algorithm is '" ALG_NAME_SPHERE "'.\n"
"  -E, --min-errors=ERRORS      The minimum number of frame errors per SNR. The default\n"
"                                 is 50.\n"
"  -r, --rng=RNG                The random number generator to use. To see a list of all\n"
"                                 available generators give 'list' as argument.\n"
"  -S, --seed=SEED              The seed for the random number generator.\n"
"  -e, --snr-end=END            The final SNR. The default is 3.0.\n"
"  -b, --snr-begin=BEGIN        The initial SNR. The default is 1.0.\n"
"  -s, --snr-step=STEP          The step size used when incrementing the SNR. The\n"
"                                 default step is 0.2.\n"
"  -t, --transpose              Transpose the basis read from INPUT.\n"
"      --help                   Display this help and exit.\n"
"      --version                Output version information and exit.";
    
    return (fprintf(file, formatstr, PROGRAM_NAME, helpstr) < 0) 
                ? EXIT_FAILURE : EXIT_SUCCESS;
}

static void parse_cmdline(int argc, char* const argv[], OPT* opt)
{
    static const char* optstring = "a:E:r:S:e:b:s:tc";
    static struct option longopt[] = {
        {"algorithm", required_argument, NULL, 'a'},
        {"min-error", required_argument, NULL, 'E'},
        {"rng", required_argument, NULL, 'r'},
        {"seed", required_argument, NULL, 'S'},
        {"snr-end", required_argument, NULL, 'e'},
        {"snr-begin", required_argument, NULL, 'b'},
        {"snr-step", required_argument, NULL, 's'},
        {"transpose", no_argument, NULL, 't'},
        {"cwords-from-stdin", no_argument, NULL, 'c'},
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'V'},
        {0, 0, 0, 0}
    };
    
    // Setting default options
    *opt = (OPT) { .alg = ALG_SPHERE_SE, .transpose = 0 };
    opt->sim = (SIM_OPTIONS) {
        .min_err = 50, .snr_begin = 1.0,
        .snr_step = 0.2, .snr_end = 3.0, .seed = 0,
        .zero_cwords = 1, .infile = NULL, .outfile = NULL,
        .rng_type = gsl_rng_default };

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
                    check(opt->alg > 0, "invalid argument to option '%c': '%s'", ch, optarg);
                }
                break;
            case 'E':
                opt->sim.min_err = strtoul(optarg, &endptr, 10);
                check(*endptr == '\0' && !(errno == ERANGE && opt->sim.min_err == ULONG_MAX),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 'r':
            {
                const gsl_rng_type** rng_types = gsl_rng_types_setup();
                if(!strcmp(optarg, "list"))
                    exit(print_rngs(stdout, rng_types));
                else
                {
                    opt->sim.rng_type = get_rng_type(optarg, rng_types);
                    check(opt->sim.rng_type, "invalid random number generator");
                }
                break;
            }
            case 'S':
                opt->sim.seed = strtoul(optarg, &endptr, 10);
                check(*endptr == '\0' && !(errno == ERANGE && opt->sim.seed == ULONG_MAX),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 'e':
                opt->sim.snr_end = strtod(optarg, &endptr);
                check(*endptr == '\0' && opt->sim.snr_end >= 0 
                        && !(errno == ERANGE && opt->sim.snr_end == HUGE_VAL),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 'b':
                opt->sim.snr_begin = strtod(optarg, &endptr);
                check(*endptr == '\0' && opt->sim.snr_begin >= 0 
                        && !(errno == ERANGE && opt->sim.snr_begin == HUGE_VAL),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 's':
                opt->sim.snr_step = strtod(optarg, &endptr);
                check(*endptr == '\0' && opt->sim.snr_step >= 0 
                        && !(errno == ERANGE && opt->sim.snr_step == HUGE_VAL),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 't':
                opt->transpose = 1;
                break;
            case 'c':
                opt->sim.zero_cwords = 0;
                break;
            case 'h':
                exit(print_help(stdout));
            case 'V':
                exit(print_version(stdout));
            default:
                goto error;
        }
    }

    // Check for mandatory arguments, i.e., input and output files.
    lcheck_pf(optind < argc, log_plain, error, "No input file given");
    opt->sim.infile = argv[optind++];

    lcheck_pf(optind < argc, log_plain, error, "No output file given");
    opt->sim.outfile = argv[optind];

    return;

error:
    fprintf(stderr, "Try '%s --help' for more information.\n", PROGRAM_NAME);
    exit(EXIT_FAILURE);
}

static int parse_and_simulate(FILE* file, OPT* opt)
{
    gsl_matrix* basis = parse_fpLLL_matrix(file, opt->transpose);
    check(basis, "error when processing input file '%s'", opt->sim.infile);

    SIMULATOR* sim = SIMULATOR_from_basis(basis, opt->alg);
    lcheck_mem(sim, error_a);

    int rc = SIMULATOR_run(sim, &opt->sim);
    SIMULATOR_free(sim);
    check(rc == 0, "simulation failed");

    return 0;

error_a:
    gsl_matrix_free(basis);
error:
    return -1;
}

int main(int argc, char* argv[])
{
    OPT opt;
    int ret = EXIT_FAILURE;
    argv[0] = PROGRAM_NAME = "lat-sim";
    parse_cmdline(argc, argv, &opt);

    FILE* infile = fopen(opt.sim.infile, "r");
    check(infile, "Could not open '%s'", opt.sim.infile);

    int rc = parse_and_simulate(infile, &opt);
    llibcheck(rc == 0, error_a, "parse_and_simulate failed");
    ret = EXIT_SUCCESS;

error_a:
    fclose(infile);
error:
    return ret;
}

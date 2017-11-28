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
#include "configuration.h"
#include "lat_sim.h"
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <stdarg.h>
#include <math.h>

typedef struct s_opt
{
    char* infile;
    char* outfile;
    int transpose;
    int alg;
    int quiet;
    SIM_OPTIONS sim;
} OPT;


static const OPT OPT_default = {
    .infile = NULL, .outfile = NULL,
    .transpose = 1, .alg = ALG_SPHERE_SE,
    .quiet = 0,
    .sim = {
        .min_err = 50, .vnr_begin = 5.0,
        .vnr_step = 1.0, .vnr_end = 20.0,
        .seed = 0, .zero_cwords = 1,
        .ser_cutoff = 1E-10,
        .fer_cutoff = 1E-10
    }
};

static int print_help(FILE* file)
{
    static const char* formatstr =
"Usage: %s [OPTION]... INPUT OUTPUT\n\n"
"Decoding simulations for lattices. The basis of the lattice is read from INPUT.\n\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -a, --algorithm=ALG          Select the decoding algorithm. To see a list of all\n"
"                                 available algorithms give 'list' as argument.\n"
"                                 The default algorithm is '%s'.\n"
"  -F, --fer-cutoff=FER         Stop the simulation when the frame error rate drops\n"
"                                 below FER. The default is %.2e.\n"
"  -E, --min-errors=ERRORS      The minimum number of frame errors per VNR. The default\n"
"                                 is %zu.\n"
"  -q, --quiet                  Supress output.\n"
"  -r, --rng=RNG                The random number generator to use. To see a list of all\n"
"                                 available generators give 'list' as argument.\n"
"  -S, --seed=SEED              The seed for the random number generator.\n"
"  -B, --ser-cutoff=SER         Stop the simulation when the symbol error rate drops\n"
"                                 below SER. The default is %.2e.\n"
"  -t, --transpose              Transpose the basis read from INPUT.\n"
"  -e, --vnr-end=END            The final VNR. The default is %.1f.\n"
"  -b, --vnr-begin=BEGIN        The initial VNR. The default is %.1f.\n"
"  -s, --vnr-step=STEP          The step size used when incrementing the VNR. The\n"
"                                 default step is %.1f.\n"
"      --help                   Display this help and exit.\n"
"      --version                Output version information and exit.\n";

    return (fprintf(file, formatstr, PROGRAM_NAME,
                algorithm_get_name(ALG_SPHERE_SE),
                OPT_default.sim.fer_cutoff,
                OPT_default.sim.min_err,
                OPT_default.sim.ser_cutoff,
                OPT_default.sim.vnr_end,
                OPT_default.sim.vnr_begin,
                OPT_default.sim.vnr_step) < 0)
                ? EXIT_FAILURE : EXIT_SUCCESS;
}

static void parse_cmdline(int argc, char* const argv[], OPT* opt)
{
    static const char* optstring = "a:B:E:F:qr:S:e:b:s:tc";
    static struct option longopt[] = {
        {"algorithm", required_argument, NULL, 'a'},
        {"fer-cutoff", required_argument, NULL, 'F'},
        {"min-error", required_argument, NULL, 'E'},
        {"quiet", no_argument, NULL, 'q'},
        {"rng", required_argument, NULL, 'r'},
        {"seed", required_argument, NULL, 'S'},
        {"ser-cutoff", required_argument, NULL, 'B'},
        {"transpose", no_argument, NULL, 't'},
        {"vnr-end", required_argument, NULL, 'e'},
        {"vnr-begin", required_argument, NULL, 'b'},
        {"vnr-step", required_argument, NULL, 's'},
        {"cwords-from-stdin", no_argument, NULL, 'c'},
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'V'},
        {0, 0, 0, 0}
    };

    opt->sim.rng_type = gsl_rng_default;

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
                    check(opt->alg >= 0,
                            "invalid argument to option '%c': '%s'", ch, optarg);
                }
                break;
            case 'B':
                opt->sim.ser_cutoff = strtod(optarg, &endptr);
                check(*endptr == '\0' && opt->sim.ser_cutoff >= 0
                        && !(errno == ERANGE && opt->sim.ser_cutoff == HUGE_VAL),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 'E':
                opt->sim.min_err = strtoul(optarg, &endptr, 10);
                check(*endptr == '\0' && !(errno == ERANGE && opt->sim.min_err == ULONG_MAX),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 'F':
                opt->sim.fer_cutoff = strtod(optarg, &endptr);
                check(*endptr == '\0' && opt->sim.fer_cutoff >= 0
                        && !(errno == ERANGE && opt->sim.fer_cutoff == HUGE_VAL),
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
                opt->sim.vnr_end = strtod(optarg, &endptr);
                check(*endptr == '\0' && opt->sim.vnr_end >= 0
                        && !(errno == ERANGE && opt->sim.vnr_end == HUGE_VAL),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 'b':
                opt->sim.vnr_begin = strtod(optarg, &endptr);
                check(*endptr == '\0' && opt->sim.vnr_begin >= 0
                        && !(errno == ERANGE && opt->sim.vnr_begin == HUGE_VAL),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 's':
                opt->sim.vnr_step = strtod(optarg, &endptr);
                check(*endptr == '\0' && opt->sim.vnr_step >= 0
                        && !(errno == ERANGE && opt->sim.vnr_step == HUGE_VAL),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 'q':
                opt->quiet = 1;
                break;
            case 't':
                opt->transpose = 0;
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
    opt->infile = argv[optind++];

    lcheck_pf(optind < argc, log_plain, error, "No output file given");
    opt->outfile = argv[optind];

    return;

error:
    fprintf(stderr, "Try '%s --help' for more information.\n", PROGRAM_NAME);
    exit(EXIT_FAILURE);
}

static int print_config(FILE* file, const OPT* opt)
{
    struct config conf[] = {
        {"algorithm", (union value) algorithm_get_name(opt->alg),
            type_str, opt->alg != ALG_SPHERE_SE},
        {"fer-cutoff", (union value) opt->sim.fer_cutoff,
            type_dbl, opt->sim.fer_cutoff != OPT_default.sim.fer_cutoff},
        {"min-error", (union value) opt->sim.min_err,
            type_size, opt->sim.min_err != OPT_default.sim.min_err},
        {"rng", (union value) opt->sim.rng_type->name,
            type_str, opt->sim.rng_type != gsl_rng_default},
        {"seed", (union value) opt->sim.seed, type_ulong, opt->sim.seed != 0},
        {"ser-cutoff", (union value) opt->sim.ser_cutoff,
            type_dbl, opt->sim.ser_cutoff != OPT_default.sim.ser_cutoff},
        {"vnr-end", (union value) opt->sim.vnr_end,
            type_dbl, opt->sim.vnr_end != OPT_default.sim.vnr_end},
        {"vnr-begin", (union value) opt->sim.vnr_begin,
            type_dbl, opt->sim.vnr_begin != OPT_default.sim.vnr_begin},
        {"vnr-step", (union value) opt->sim.vnr_step,
            type_dbl, opt->sim.vnr_step != OPT_default.sim.vnr_step},
        {"transpose", (union value) !opt->transpose, type_bool, 1},
        {0, (union value) 0, type_bool, 0}
    };

    return util_print_genby_and_config(file, PROGRAM_NAME, NULL, NULL, conf);
}

static int print_info(FILE* file, const SIMULATOR* sim, const OPT* opt)
{
    const char* formatstr =
        "\nDecoding simulations with lattices\n"
        "VNR (dB) from %.2f to %.2f in increments of %.2f\n"
        "%zu frame errors awaited per VNR\n"
        "Symbol error rate cutoff: %.2e\n"
        "Frame error rate cutoff: %.2e\n"
        "Lattice parameters:\n"
        "  dimension: %zu\n"
        "  rank: %zu\n"
        "  rate: %.4f\n"
        "Decoding algorithm: %s\n"
        "RNG: %s\n"
        "Seed: %lu\n\n";

    int rc = fprintf(file, formatstr,
        opt->sim.vnr_begin, opt->sim.vnr_end, opt->sim.vnr_step,
        opt->sim.min_err, opt->sim.ser_cutoff,
        opt->sim.fer_cutoff,
        SIMULATOR_get_dimension(sim), SIMULATOR_get_rank(sim),
        SIMULATOR_get_rate(sim),
        algorithm_get_name(opt->alg),
        opt->sim.rng_type->name, opt->sim.seed);

    return rc > 0 ? LT_SUCCESS : LT_ESYSTEM;
}

static int parse_and_simulate(FILE* infile, FILE* outfile, OPT* opt)
{
    int lt_errno = LT_FAILURE;
    LSC_ARGS callback_args = { .ui = stdout, .file = outfile };

    gsl_matrix* basis = parse_fpLLL_matrix(infile, opt->transpose);
    check(basis, "error when processing input file '%s'", opt->infile);

    SIMULATOR* sim;
    lt_errno = SIMULATOR_from_basis(&sim, basis, opt->alg);
    lt_lcheck(lt_errno, error_a, "SIMULATOR_from_basis failed");

    check_se(callback_args.file, lt_errno, LT_ESYSTEM,
            "could not open '%s' for writing", opt->outfile);

    opt->sim.seed = opt->sim.seed ? opt->sim.seed : get_random_seed();

    if(opt->quiet)
    {
        SIMULATOR_set_callbacks(sim, lat_sim_vnr_callback_quiet,
                                lat_sim_start_callback_quiet,
                                NULL,
                                &callback_args);
    }
    else
    {
        SIMULATOR_set_callbacks(sim, lat_sim_vnr_callback_std,
                                lat_sim_start_callback_std,
                                lat_sim_end_callback,
                                &callback_args);
        lt_errno = print_info(callback_args.ui, sim, opt);
        lt_check(lt_errno, "print_info failed");
    }

    lt_errno = print_config(callback_args.file, opt);
    lt_check(lt_errno, "print_config failed");

    lt_errno = SIMULATOR_run(sim, &opt->sim);
    SIMULATOR_free(sim);
    lt_check(lt_errno, "SIMULATOR_run failed");

    return lt_errno;

error_a:
    gsl_matrix_free(basis);
error:
    return lt_errno;
}

int main(int argc, char* argv[])
{
    OPT opt = OPT_default;
    int ret = EXIT_FAILURE;
    argv[0] = PROGRAM_NAME = "lat-sim";

    gsl_set_error_handler_off();
    parse_cmdline(argc, argv, &opt);

    FILE* infile = fopen(opt.infile, "r");
    check(infile, "Could not open '%s'", opt.infile);

    FILE* outfile = fopen(opt.outfile, "w");
    lcheck(outfile, error_a, "Could not open '%s'", opt.outfile);

    int rc = parse_and_simulate(infile, outfile, &opt);
    llibcheck(rc == 0, error_b, "parse_and_simulate failed");
    ret = EXIT_SUCCESS;

error_b:
    fclose(outfile);
error_a:
    fclose(infile);
error:
    return ret;
}

/* A program that generates random points in Euclidian space.
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
#include "parse.h"
#include "print_util.h"
#include "utility.h"
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <stdarg.h>
#include <time.h>
#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <assert.h>

#define HYP_MODE_NO 0
#define HYP_MODE_A 1
#define HYP_MODE_H 2
#define HYP_MODE_B 3

typedef enum enum_mode
{
    MODE_STANDARD = 1,
    MODE_ZERO
} Mode;

typedef struct s_options
{
    char* input;
    char* output;
    long* hyperplane;
    gsl_matrix* basis;
    const gsl_rng_type* rng_type;
    int (*print_point)(FILE*, const double*, size_t, const PRINTING_FMT*);
    size_t cword_len;
    size_t cword_num;
    unsigned long seed;
    long min;
    long max;
    Mode mode;
    enum PrintingFmt format;
    int output_conf;
    int hyp_mode;
    int get_point_switch;
    int sloppy;
    int min_set;
    int max_set;
    int rows_as_basis;
} RP_OPT;

static int print_help(FILE* file)
{
    static const char* formatstr =
"Usage: %s [OPTION]...\n"
"  or:  %s [OPTION]... OUTPUT\n\n%s\n";

    static const char* helpstr =
"Generates random points in Euclidean space. Outputs to stdout\n"
"if no output file is given.\n\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -A, --A-n                    Generates points on the hyperplane with normal vector\n"
"                                 (1, ..., 1).\n"
"  -C, --no-config              Do not output configuration.\n"
"  -d, --dimension=DIM          The dimesion of the Euclidean space from which to\n"
"                                 sample points.\n"
"  -H, --hyperplane=FILE        Generates points on the hyperplane defined by the\n"
"                                 coefficients read from FILE. Give '-' as argument \n"
"                                 to read from stdin.\n"
"  -B, --basis=FILE             Generates points in the space spanned by the given basis.\n"
"                                 Give '-' as argument to read from stdin.\n"
"  -f, --output-format=FMT      The output format that will be used if 'R' is given.\n"
"                                 Give 'list' as argument to get a list of all formats.\n"
"  -n, --num-points=NUM         The number of codewords to generate. Zero (0) makes the\n"
"                                 encoder run until it is killed.\n"
"  -m, --min=NUM                The smallest possible value for a coordinate (inclusive).\n"
"                                 Defaults to the min-value of the chosen rng.\n"
"  -M, --max=NUM                The largest possible value for a coordinate (exclusive).\n"
"                                 Defaults to the max-value of the chosen rng.\n"
"  -r, --rng=RNG                The random number generator to use. To see a list of all\n"
"                                 available generators give 'list' as argument.\n"
"  -R, --readable-output        Produce readable output instead of binary output.\n"
"  -s, --sloppy                 Lower quality output (especially when the range selected\n"
"                                 by -m and -M is large). However, this option speeds up\n"
"                                 the generator considerably.\n"
"  -S, --seed=SEED              The seed for the random number generator.\n"
"  -t, --transpose              Transpose the basis read from FILE.\n"
"  -Z, --force-zero-binary      Makes the encoder output all zero codewords in binary\n"
"                                 format. Note that this is a developer option meant\n"
"                                 for testing.\n"
"      --help                   Display this help and exit.\n"
"      --version                Output version information and exit.\n\n"
"If more that one of 'A', 'B', and 'H' are given, then the one specified last\n"
"takes precedence.\n";

    return (fprintf(file, formatstr, PROGRAM_NAME, PROGRAM_NAME, helpstr) < 0)
                ? EXIT_FAILURE : EXIT_SUCCESS;
}

static void parse_cmdline(int argc, char* const argv[], RP_OPT* opt)
{
    static const char* optstring = "AB:Cd:f:H:n:m:M:r:RsS:tZ";
    static struct option longopt[] = {
        {"A-n", no_argument, NULL, 'A'},
        {"no-config", no_argument, NULL, 'C'},
        {"dimension", required_argument, NULL, 'd'},
        {"basis", required_argument, NULL, 'B'},
        {"hyperplane", required_argument, NULL, 'H'},
        {"output-format", required_argument, NULL, 'f'},
        {"num-points", required_argument, NULL, 'n'},
        {"min", required_argument, NULL, 'm'},
        {"max", required_argument, NULL, 'M'},
        {"rng", required_argument, NULL, 'r'},
        {"readable-output", no_argument, NULL, 'R'},
        {"sloppy", no_argument, NULL, 's'},
        {"seed", required_argument, NULL, 'S'},
        {"transpose", no_argument, NULL, 't'},
        {"force-zero-binary", no_argument, NULL, 'Z'},
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'V'},
        {0, 0, 0, 0}
    };

    // Setting default options
    *opt = (RP_OPT) {
        .input = NULL, .output = NULL, .cword_len = 0,
        .cword_num = 0, .seed = 0,
        .min = 0, .max = 0, .sloppy = 0,
        .output_conf = 1, .hyp_mode = HYP_MODE_NO,
        .min_set = 0, .max_set = 0,
        .print_point = print_dvector_binary, .hyperplane = NULL,
        .basis = NULL, .rows_as_basis = 1,
        .mode = MODE_STANDARD, .rng_type = gsl_rng_default,
        .format = PRINTING_FMT_DEFAULT
    };

    // Parsing the command line
    int ch;
    char* endptr;
    while((ch = getopt_long(argc, argv, optstring, longopt, NULL)) != -1)
    {
        switch(ch)
        {
            case 'A':
                opt->hyp_mode = HYP_MODE_A;
                break;
            case 'C':
                opt->output_conf = 0;
                break;
            case 'd':
                opt->cword_len = strtoul(optarg, &endptr, 10);
                check(*endptr == '\0' && !(errno == ERANGE && opt->cword_len == ULONG_MAX),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 'B':
                opt->hyp_mode = HYP_MODE_B;
                opt->input = !strcmp(optarg, "-") ? NULL : optarg;
                break;
            case 'H':
                opt->hyp_mode = HYP_MODE_H;
                opt->input = !strcmp(optarg, "-") ? NULL : optarg;
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
            case 'n':
                opt->cword_num = strtoul(optarg, &endptr, 10);
                check(*endptr == '\0' && !(errno == ERANGE && opt->cword_num == ULONG_MAX),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 'm':
                opt->min = strtol(optarg, &endptr, 10);
                check(*endptr == '\0' &&
                    !(errno == ERANGE && (opt->min == LONG_MAX || opt->min == LONG_MIN)),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                opt->min_set = 1;
                break;
            case 'M':
                opt->max = strtol(optarg, &endptr, 10);
                check(*endptr == '\0' &&
                    !(errno == ERANGE && (opt->max == LONG_MAX || opt->max == LONG_MIN)),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                opt->max_set = 1;
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
            case 'R':
                opt->print_point = print_dvector_std;
                break;
            case 's':
                opt->sloppy = 1;
                break;
            case 'S':
            {
                opt->seed = strtoul(optarg, &endptr, 10);
                check(*endptr == '\0' && !(errno == ERANGE && opt->seed == ULONG_MAX),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            }
            case 't':
                opt->rows_as_basis = 0;
                break;
            case 'Z':
                opt->mode = MODE_ZERO;
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
    if(opt->hyp_mode != HYP_MODE_B)
        check(opt->cword_len > 0, "missing mandatory option -- '%c'", 'd');

    return;

error:
    fprintf(stderr, "Try '%s --help' for more information.\n", PROGRAM_NAME);
    exit(EXIT_FAILURE);
}

static int generate_zero_helper(FILE* file, void* cword, size_t size,
                                size_t n, size_t num_cwords)
{
    if(num_cwords)
    {
        size_t i = 0;
        do {
            check(fwrite(cword, size, n, file) == n, "fwrite failed");
        } while(++i < num_cwords);
    }
    else
    {
        while(1)
            check(fwrite(cword, size, n, file) == n, "fwrite failed");
    }
    return 0;

error:
    return -1;
}

static int generate_zero_pnt(FILE* file, size_t size, size_t n, size_t num_cwords)
{
    void* cword = calloc(n, size);
    check_mem(cword);

    int rc = generate_zero_helper(file, cword, size, n, num_cwords);
    free(cword);
    libcheck(rc == 0, "generate_zero_helper failed");
    return 0;

error:
    return -1;
}

static void get_point(double* p, gsl_rng* rng, RP_OPT* opt)
{
    switch(opt->get_point_switch)
    {
        case 0:
            rnd_point_get(p, opt->cword_len, rng);
            break;
        case 1:
            rnd_point_min_max(p, opt->cword_len, rng, opt->min, opt->max);
            break;
        case 2:
            rnd_point_get_A_N(p, opt->cword_len, rng);
            break;
        case 3:
            rnd_point_min_max_A_N(p, opt->cword_len, rng, opt->min, opt->max);
            break;
        case 4:
            rnd_point_get_SPC(p, opt->cword_len, rng, opt->hyperplane);
            break;
        case 5:
            rnd_point_min_max_SPC(p, opt->cword_len, rng, opt->hyperplane,
                                    opt->min, opt->max);
            break;
        case 8:
        case 9:
            rnd_point_min_max_fast(p, opt->cword_len, rng, opt->min, opt->max);
            break;
        case 10:
        case 11:
            rnd_point_min_max_A_N_fast(p, opt->cword_len, rng, opt->min, opt->max);
            break;
        case 12:
        case 13:
            rnd_point_min_max_SPC_fast(p, opt->cword_len, rng, opt->hyperplane,
                                        opt->min, opt->max);
            break;
        default:
            break;
    }
}

static int set_range(int custom_range, gsl_rng* rng, RP_OPT* opt)
{
    if(custom_range)
    {
        debug("Custom range!");
        if(!opt->min_set)
            opt->min = gsl_rng_min(rng);
        if(!opt->max_set)
            opt->max = gsl_rng_max(rng);

        check(opt->min < opt->max, "min (%ld) must be smaller than max (%ld)",
                opt->min, opt->max);
        long max_range = gsl_rng_max(rng) - gsl_rng_min(rng);
        check(opt->max - opt->min <= max_range,
                "Selected range too large for the selected rng;");
    }
    else if(opt->sloppy)
    {
        debug("Sloppy!");
        opt->min = gsl_rng_min(rng);
        opt->max = gsl_rng_max(rng);
    }

    return LT_SUCCESS;

error:
    return LT_EINVAL;
}

static int generate_standard(FILE* outfile, RP_OPT* opt)
{
    int ret = -1;
    double* cword = malloc(opt->cword_len * sizeof(double));
    check_mem(cword);

    opt->seed = opt->seed ? opt->seed : get_random_seed();
    gsl_rng* rng = rng_alloc_and_seed(opt->rng_type, opt->seed);
    lcheck_mem(rng, error_a);

    int custom_range = (opt->min_set || opt->max_set);
    opt->get_point_switch = custom_range + 2 * opt->hyp_mode + 8 * opt->sloppy;
    llibcheck(set_range(custom_range, rng, opt) == 0, error_b, "set_range failed");

    const PRINTING_FMT* fmt = printing_fmt_get(opt->format);
    size_t i = 0;
    do {
        get_point(cword, rng, opt);
        lcheck(opt->print_point(outfile, cword, opt->cword_len, fmt) == 0,
                error_b, "print_point failed");
    } while(!opt->cword_num || ++i < opt->cword_num);
    ret = 0;

error_b:
    gsl_rng_free(rng);
error_a:
    free(cword);
error:
    return ret;
}

static int generate_with_basis(FILE* outfile, RP_OPT* opt)
{
    int lt_errno = LT_SUCCESS;
    const gsl_matrix* B = opt->basis;
    const size_t n = B->size1;
    const size_t m = B->size2;
    check_se(n >= m, lt_errno, LT_EINVAL,
            "the basis vectors are not linearly independent");

    const size_t Q_size = n * n;
    const size_t R_size = n * m;
    double* data = malloc((Q_size + R_size + n + m) * sizeof(double));
    check_se_mem(data, lt_errno, LT_ESYSTEM);

    double* r = data + Q_size + R_size;
    double* x = r + n;

    gsl_vector_view v_r = gsl_vector_view_array(r, n);
    gsl_vector_view v_x = gsl_vector_view_array(x, m);

    gsl_matrix_view m_Q = gsl_matrix_view_array(data, n, n);
    gsl_matrix_view m_R = gsl_matrix_view_array(data + Q_size, n, m);

    lt_errno = utility_compute_QR_decomposition(&m_Q.matrix, &m_R.matrix, B);
    lt_lcheck(lt_errno, error_a, "utility_compute_QR_decomposition failed");
    int rc = utility_Rmm_is_not_singular(&m_R.matrix, 10E-10);
    lcheck_se(rc, error_a, lt_errno, LT_EINVAL,
            "the basis vectors are not linearly independent");
    gsl_matrix_view m_Q1 = gsl_matrix_submatrix(&m_Q.matrix, 0, 0, n, m);

    opt->seed = opt->seed ? opt->seed : get_random_seed();
    gsl_rng* rng = rng_alloc_and_seed(opt->rng_type, opt->seed);
    lcheck_se_mem(data, error_a, lt_errno, GSL_ENOMEM);

    int custom_range = (opt->min_set || opt->max_set);
    opt->get_point_switch = custom_range + 8 * opt->sloppy;
    lt_errno = set_range(custom_range, rng, opt);
    lt_llibcheck(lt_errno, error_b, "set_range failed");

    const PRINTING_FMT* fmt = printing_fmt_get(opt->format);
    opt->cword_len = m;
    size_t i = 0;
    do {
        get_point(x, rng, opt);
        gsl_blas_dgemv(CblasNoTrans, 1, &m_Q1.matrix, &v_x.vector, 0, &v_r.vector);
        lt_errno = opt->print_point(outfile, r, n, fmt);
        lt_lcheck(lt_errno, error_b, "print_point failed");
    } while(!opt->cword_num || ++i < opt->cword_num);

error_b:
    gsl_rng_free(rng);
error_a:
    free(data);
error:
    return lt_errno;
}

static int print_conf_binary(FILE* file, RP_OPT* opt)
{
    RND_PNT_CONF conf = {
        .hyperplane = opt->hyperplane, .dimension = opt->cword_len,
        .num_cwords = opt->cword_num, .A_n = opt->hyp_mode };

    return RND_PNT_CONF_write(file, &conf);
}

static int print_conf_plain(FILE* file, RP_OPT* opt)
{
    // TODO
    return 0;
}

static int read_hyperplane_coeffs(RP_OPT* opt)
{
    int ret = -1;
    FILE* infile = stdin;
    if(opt->input)
    {
        infile = fopen(opt->input, "r");
        check(infile, "Could not open '%s' for reading", opt->input);
    }

    opt->hyperplane = malloc(opt->cword_len * sizeof(long));
    lcheck_mem(opt->hyperplane, error_a);

    int rc = parse_hyperplane_coeffs(infile, opt->hyperplane, opt->cword_len);
    llibcheck(rc == 0, error_a, "parse_hyperplane_coeffs failed");
    ret = 0;

error_a:
    if(opt->input)
        fclose(infile);
error:
    return ret;
}

static int read_basis(RP_OPT* opt)
{
    int ret = -1;
    FILE* infile = stdin;
    if(opt->input)
    {
        infile = fopen(opt->input, "r");
        check(infile, "Could not open '%s' for reading", opt->input);
    }

    opt->basis = parse_fpLLL_matrix(infile, opt->rows_as_basis);
    llibcheck(opt->basis, error_a, "parse_fpLLL_matrix failed");
    opt->cword_len = opt->basis->size1;
    ret = 0;

error_a:
    if(opt->input)
        fclose(infile);
error:
    return ret;
}

int main(int argc, char* argv[])
{
    RP_OPT opt;
    int ret = EXIT_FAILURE;
    argv[0] = PROGRAM_NAME = "rnd-point";

    gsl_set_error_handler_off();
    parse_cmdline(argc, argv, &opt);

    FILE* outfile = stdout;
    if(opt.output)
    {
        outfile = fopen(opt.output, "w");
        check(outfile, "Could not open '%s' for writing", opt.output);
    }

    int rc;
    if(opt.hyp_mode == HYP_MODE_H)
    {
        rc = read_hyperplane_coeffs(&opt);
        llibcheck(rc == 0, error_a, "read_hyperplane_coeffs failed");
    }
    else if(opt.hyp_mode == HYP_MODE_B)
    {
        rc = read_basis(&opt);
        llibcheck(rc == 0, error_a, "read_basis failed");
    }

    if(opt.output_conf)
    {
        rc = (opt.print_point == print_dvector_binary ? print_conf_binary(outfile, &opt)
                : print_conf_plain(outfile, &opt));
        llibcheck(rc == 0, error_a, "print_conf_* failed");
    }

    switch(opt.mode)
    {
        default:
        case MODE_STANDARD:
            if(opt.hyp_mode == HYP_MODE_B)
                rc = generate_with_basis(outfile, &opt);
            else
                rc = generate_standard(outfile, &opt);
            break;
        case MODE_ZERO:
            rc = generate_zero_pnt(outfile, sizeof(double), opt.cword_len, opt.cword_num);
            break;
    }
    llibcheck(rc == 0, error_a, "generate_* failed");
    ret = EXIT_SUCCESS;

error_a:
    gsl_matrix_free(opt.basis);
    free(opt.hyperplane);
    if(opt.output)
        fclose(outfile);
error:
    return ret;
}

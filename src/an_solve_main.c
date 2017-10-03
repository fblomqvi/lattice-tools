/* A program that solves the CVP on (R^n, A_n)
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
#include "an_solver.h"
#include "rnd_point.h"
#include "print_util.h"
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <stdarg.h>
#include <time.h>
#include <stdlib.h>

typedef struct s_options 
{
    char* input;
    char* output;
    size_t cword_len;
    size_t cword_num;
    enum PrintingFmt format;
    int no_config;
    int binary_out;
    void (*an_solve)(long*, double*, AN_WS*);
    int (*plp)(FILE*, const long*, size_t, const PRINTING_FMT*);
} OPT;

static const OPT OPT_default = {
    .input = NULL, .output = NULL, 
    .cword_len = 0, .cword_num = 0, 
    .no_config = 0, .binary_out = 1, 
    .an_solve = AN_solve,
    .plp = print_lvector_std,
    .format = PRINTING_FMT_DEFAULT

};

static int print_help(FILE* file)
{
    static const char* formatstr = 
"Usage: %s [OPTION]...\n"
"  or:  %s [OPTION]... OUTPUT\n\n%s\n";

    static const char* helpstr = 

"Solves the closest vector problem on (R^d, A_{d-1}). Reads the points\n"
"to decode from stdin in the binary format produced by rnd-point.\n"
"Outputs to stdout if no output file is given.\n\n"
"Mandatory arguments to long options are mandatory for short options too.\n"
"  -d, --dimension=NUM          Assumes that the input points lie in R^NUM.\n"
"  -f, --output-format=FMT      The output format that will be used if 'R' is given.\n"
"                                 Give 'list' as argument to get a list of all formats.\n"
"  -n, --num-points=NUM         The number of codewords to decode. Zero (0) makes the\n"
"                                 solver run until it runs out of input.\n"
"  -A, --A-n                    Assumes that the input points lie on the hyperplane\n"
"                                 with normal vector (1, ..., 1).\n"
"  -C, --no-config              Do not try to read the configuration from stdin.\n"
"                                 Use of this option makes '-d' mandatory.\n"
"  -R, --readable-output        Produce readable output instead of binary output.\n"
"  -F, --output-float           Output floating point numbers instead of integers.\n"
"      --help                   Display this help and exit.\n"
"      --version                Output version information and exit.";
    
    return (fprintf(file, formatstr, PROGRAM_NAME, PROGRAM_NAME, helpstr) < 0) 
                ? EXIT_FAILURE : EXIT_SUCCESS;
}

static void parse_cmdline(int argc, char* const argv[], OPT* opt)
{
    static const char* optstring = "d:f:n:AFRC";
    static struct option longopt[] = {
        {"dimension", required_argument, NULL, 'd'},
        {"output-format", required_argument, NULL, 'f'},
        {"num-points", required_argument, NULL, 'n'},
        {"A-n", no_argument, NULL, 'A'},
        {"readable-output", no_argument, NULL, 'R'},
        {"output-float", no_argument, NULL, 'F'},
        {"no-config", no_argument, NULL, 'C'},
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'V'},
        {0, 0, 0, 0}
    };

    int float_out = 0;

    // Parsing the command line
    int ch;
    char* endptr;
    while((ch = getopt_long(argc, argv, optstring, longopt, NULL)) != -1)
    {
        switch(ch)
        {
            case 'd':
                opt->cword_len = strtoul(optarg, &endptr, 10);
                check(*endptr == '\0' && !(errno == ERANGE && opt->cword_len == ULONG_MAX),
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
            case 'n':
                opt->cword_num = strtoul(optarg, &endptr, 10);
                check(*endptr == '\0' && !(errno == ERANGE && opt->cword_num == ULONG_MAX),
                    "invalid argument to option '%c': '%s'", ch, optarg);
                break;
            case 'A':
                opt->an_solve = AN_solve_np;
                break;
            case 'C':
                opt->no_config = 1;
                break;
            case 'F':
                float_out = 1;
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

    //lcheck_pf(optind < argc, log_plain, error, "No input file given");
    //opt->input = argv[optind++];
    
    check(!opt->no_config || opt->cword_len, "'-C' makes '-d' mandatory");

    if(optind < argc)
        opt->output = argv[optind];

    if(float_out)
    {
        opt->plp = opt->binary_out ? 
            print_lvector_binary_float : print_lvector_as_double_std;
    }
    else
    {
        opt->plp = opt->binary_out ? 
            print_lvector_binary : print_lvector_std;
    }

    return;

error:
    fprintf(stderr, "Try '%s --help' for more information.\n", PROGRAM_NAME);
    exit(EXIT_FAILURE);
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
        conf->dimension = opt->cword_len;
        conf->num_cwords = opt->cword_num;
        conf->A_n = (opt->an_solve != OPT_default.an_solve);
    }
    else
    {
        if(opt->cword_len != OPT_default.cword_len)
            conf->dimension = opt->cword_len;
        if(opt->cword_num != OPT_default.cword_num)
            conf->num_cwords = opt->cword_num;
        if(opt->an_solve != OPT_default.an_solve)
            conf->A_n = 1;
    }
    debug("dimension: %zu", conf->dimension);
    debug("num_cwords: %zu", conf->num_cwords);
    debug("A-n: %d", conf->A_n);

    return conf;

error:
    return NULL;
}

static int solve(FILE* outfile, const OPT* opt)
{
    int ret = EXIT_FAILURE;
    RND_PNT_CONF* conf = get_config_binary(stdin, opt);
    libcheck(conf, "get_config_binary failed");

    double* cword = malloc(conf->dimension * sizeof(double));
    llibcheck_mem(cword, error_a);

    long* clp = malloc(conf->dimension * sizeof(long));
    llibcheck_mem(clp, error_b);

    AN_WS* ws = AN_WS_alloc(conf->dimension);
    llibcheck_mem(ws, error_c);

    if(opt->binary_out)
    {
        size_t i = 0;
        while(!conf->num_cwords || i < conf->num_cwords)
        {
            int rc = read_cword_binary(cword, conf->dimension, sizeof(double));
            llibcheck(rc == 0, error_d, "read_cword_binary failed");
            
            opt->an_solve(clp, cword, ws);
            rc = opt->plp(outfile, clp, conf->dimension, NULL);
            llibcheck(rc == 0, error_d, "print_lvector* failed");
            i++;
        }
    }
    else
    {
        const PRINTING_FMT* fmt = printing_fmt_get(opt->format);
        size_t i = 0;
        while(!conf->num_cwords || i < conf->num_cwords)
        {
            int rc = read_cword_binary(cword, conf->dimension, sizeof(double));
            llibcheck(rc == 0, error_d, "read_cword_binary failed");

            rc = print_dvector_cword(outfile, cword, conf->dimension, fmt);
            llibcheck(rc == 0, error_d, "print_codeword failed");
            
            opt->an_solve(clp, cword, ws);
            rc = opt->plp(outfile, clp, conf->dimension, fmt);
            llibcheck(rc == 0, error_d, "print_lvector* failed");
            i++;
        }
    }

    ret = EXIT_SUCCESS;

error_d:
    AN_WS_free(ws);
error_c:
    free(clp);
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
    argv[0] = PROGRAM_NAME = "an-solve";
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

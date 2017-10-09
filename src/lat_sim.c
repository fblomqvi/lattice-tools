/* lat_sim.c
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
#include "lat_sim.h"

#define TABLE_VLINE "+----------+--------------+---------------+"\
                    "------------+------------------+-------------------+"

static const char* table_headers[] = {
    "VNR (dB)",
    "Frame-errors",
    "Symbol-errors",
    "Frames",
    "Frame-error rate",
    "Symbol-error rate"
};

static const char* get_termination_reason_str(int reason);

static int print_table_footer(FILE* file);
static int print_table_headers(FILE* file);

static int print_file_headers(FILE* file);

static int print_vnr_result_file(FILE* file, double vnr, size_t frame_errs, 
                        size_t bit_errs, size_t frames, size_t bits);
static int print_vnr_result_user(FILE* file, double vnr, size_t frame_errs, 
                        size_t bit_errs, size_t frames, size_t bits);

int lat_sim_end_callback(const SIM_STATUS* status, void* args)
{
    LSC_ARGS* a = (LSC_ARGS*) args;
    const char* formatstr =
        "\nSimulation finished\n"
        "Processed %zu lattice points (frames)\n"
        "Simulation termination reason: %s\n";

    int lt_errno = print_table_footer(a->ui);
    lt_libcheck(lt_errno, "print_table_footer failed");

    lt_errno = fprintf(a->ui, formatstr, status->total,
                get_termination_reason_str(status->termination_reason)) > 0 ?
                LT_SUCCESS : LT_ESYSTEM;
error:
    return lt_errno;
}

int lat_sim_start_callback_std(const SIM_STATUS* status __attribute__((__unused__)),
                                void* args)
{
    LSC_ARGS* a = (LSC_ARGS*) args;
    int lt_errno = print_table_headers(a->ui);
    lt_libcheck(lt_errno, "print_table_headers failed");
    lt_errno = print_file_headers(a->file);
    lt_libcheck(lt_errno, "print_table_headers failed");

error:
    return lt_errno;
}

int lat_sim_start_callback_quiet(const SIM_STATUS* status __attribute__((__unused__)),
                                void* args)
{ 
    LSC_ARGS* a = (LSC_ARGS*) args;
    return print_file_headers(a->file); 
}

int lat_sim_vnr_callback_std(const SIM_STATUS* status, void* args)
{
    LSC_ARGS* a = (LSC_ARGS*) args;
    int lt_errno = print_vnr_result_user(a->ui, status->vnr, 
                        status->frame_errs, status->bit_errs, 
                        status->frames, status->frames * status->n);
    lt_libcheck(lt_errno, "print_vnr_result_user failed");
    lt_errno = print_vnr_result_file(a->file, status->vnr, 
                        status->frame_errs, status->bit_errs, 
                        status->frames, status->frames * status->n);
    lt_libcheck(lt_errno, "print_vnr_result_file failed");
error:
    return lt_errno;
}

int lat_sim_vnr_callback_quiet(const SIM_STATUS* status, void* args)
{
    LSC_ARGS* a = (LSC_ARGS*) args;
    return print_vnr_result_file(a->file, status->vnr, 
                        status->frame_errs, status->bit_errs, 
                        status->frames, status->frames * status->n);
}

static const char* get_termination_reason_str(int reason)
{
    switch(reason)
    {
        case SIM_TERMINATION_VNR:
            return "the final VNR was reached";
        case SIM_TERMINATION_BER:
            return "the symbol error rate cutoff was reached";
        case SIM_TERMINATION_FER:
            return "the frame error rate cutoff was reached";
        default:
            return "unknown termination reason";
    }
}

static int print_table_footer(FILE* file)
{ return fprintf(file, TABLE_VLINE "\n") > 0 ? LT_SUCCESS : LT_ESYSTEM; }

static int print_table_headers(FILE* file)
{
    return fprintf(file, TABLE_VLINE "\n| %-9s| %-13s| %-14s| %-11s| %-17s| %-18s|\n" 
            TABLE_VLINE "\n", 
            table_headers[0], table_headers[1], table_headers[2],
            table_headers[3], table_headers[4], table_headers[5]) > 0 ?
        LT_SUCCESS : LT_ESYSTEM;
}

static int print_file_headers(FILE* file)
{
    return fprintf(file, "%s\t%s\t%s\t%s\t%s\t%s\n",
            table_headers[0], table_headers[1], table_headers[2],
            table_headers[3], table_headers[4], table_headers[5]) > 0 ?
        LT_SUCCESS : LT_ESYSTEM;
}
static int print_vnr_result_user(FILE* file, double vnr, size_t frame_errs, 
                        size_t bit_errs, size_t frames, size_t bits)
{
    const double fer = (double) frame_errs / frames;
    const double ber = (double) bit_errs / bits;
    int rc =  fprintf(file, "|%9.4f |%13zu |%14zu |%11zu |%17.8e |%18.8e |\n", 
                    vnr, frame_errs, bit_errs, frames, fer, ber); 
    libcheck(rc >= 0, "printing error");
    return LT_SUCCESS;

error:
    return LT_ESYSTEM;
}

static int print_vnr_result_file(FILE* file, double vnr, size_t frame_errs, 
                        size_t bit_errs, size_t frames, size_t bits)
{
    const double fer = (double) frame_errs / frames;
    const double ber = (double) bit_errs / bits;
    int rc =  fprintf(file, "%f\t%zu\t%zu\t%zu\t%.8e\t%.8e\n", 
                    vnr, frame_errs, bit_errs, frames, fer, ber);
    libcheck(rc >= 0, "printing error");
    libcheck(fflush(file) == 0, "fflush failed");
    return LT_SUCCESS;

error:
    return LT_ESYSTEM;
}

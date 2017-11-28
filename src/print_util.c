/* print_util.c
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
#include "print_util.h"
#include "lt_errno.h"

#define FMT_NAME_DEFAULT "default"
#define FMT_NAME_FPLLL "fplll"
#define FMT_NAME_MATLAB "matlab"
#define FMT_NAME_MATHEMATICA "mathematica"
#define FMT_NAME_PLAIN "plain"

static const PRINTING_FMT standard_formats[] = {
    {"[", "]", "[ ", " ]", " ", "\n", "%f", "%ld"},
    {"[", "]", "[", "]", " ", "\n", "%f", "%ld"},
    {"[", "]", "[", "]", ",", ";\n", "%f", "%ld"},
    {"{", "}", "{", "}", ",", ",\n", "%f", "%ld"},
    {"", "", "", "", " ", "\n", "%f", "%ld"}
};

static const char* format_names[] = {
    FMT_NAME_DEFAULT,
    FMT_NAME_FPLLL,
    FMT_NAME_MATLAB,
    FMT_NAME_MATHEMATICA,
    FMT_NAME_PLAIN
};

const char* printing_fmt_get_name(enum PrintingFmt fmt)
{ return format_names[fmt]; }

const PRINTING_FMT* printing_fmt_get(enum PrintingFmt fmt)
{ return standard_formats + fmt; }

int printing_fmt_parse_name(const char* name)
{
    for(int i = PRINTING_FMT_DEFAULT; i < PRINTING_FMT_MAX; i++)
        if(!strcmp(name, format_names[i])) return i;

    return LT_FAILURE;
}

int printing_fmt_print_names(FILE* file)
{
    libcheck(fprintf(file, "Available printing formats are:\n") > 0, "printing error");
    for(int i = PRINTING_FMT_DEFAULT; i < PRINTING_FMT_MAX; i++)
        libcheck(fprintf(file, "%s\n", format_names[i]) > 0, "printing error");

    return LT_SUCCESS;

error:
    return LT_FAILURE;
}

int print_dvector(FILE* file, const double* vec,
                size_t vec_len, const PRINTING_FMT* fmt,
                const char* prefix, const char* postfix)
{
    libcheck(fprintf(file, "%s%s", prefix, fmt->vstart_bracket) >= 0, "printing failed");
    for(size_t i = 0; i < vec_len-1; i++)
    {
        libcheck(fprintf(file, fmt->double_format, vec[i]) >= 0, "printing failed");
        libcheck(fprintf(file, "%s", fmt->element_separator) >= 0, "printing failed");
    }
    libcheck(fprintf(file, fmt->double_format, vec[vec_len-1]) >= 0, "printing failed");
    libcheck(fprintf(file, "%s%s", fmt->vend_bracket, postfix) >= 0, "printing failed");
    return LT_SUCCESS;

error:
    return LT_FAILURE;
}

int print_dvector_binary(FILE* file, const double* vec, size_t vec_len,
                        const PRINTING_FMT* fmt __attribute__((__unused__)))
{ return fwrite(vec, sizeof(double), vec_len, file) == vec_len ? 0 : -1; }

int print_lvector(FILE* file, const long* vec,
                size_t vec_len, const PRINTING_FMT* fmt,
                const char* prefix, const char* postfix)
{
    libcheck(fprintf(file, "%s%s", prefix, fmt->vstart_bracket) >= 0, "printing failed");
    for(size_t i = 0; i < vec_len-1; i++)
    {
        libcheck(fprintf(file, fmt->integer_format, vec[i]) >= 0, "printing failed");
        libcheck(fprintf(file, "%s", fmt->element_separator) >= 0, "printing failed");
    }
    libcheck(fprintf(file, fmt->integer_format, vec[vec_len-1]) >= 0, "printing failed");
    libcheck(fprintf(file, "%s%s", fmt->vend_bracket, postfix) >= 0, "printing failed");
    return LT_SUCCESS;

error:
    return LT_FAILURE;
}

int print_lvector_as_double(FILE* file, const long* vec,
                            size_t vec_len, const PRINTING_FMT* fmt,
                            const char* prefix, const char* postfix)
{
    libcheck(fprintf(file, "%s%s", prefix, fmt->vstart_bracket) >= 0, "printing failed");
    for(size_t i = 0; i < vec_len-1; i++)
    {
        libcheck(fprintf(file, fmt->double_format, (double) vec[i]) >= 0, "printing failed");
        libcheck(fprintf(file, "%s", fmt->element_separator) >= 0, "printing failed");
    }
    libcheck(fprintf(file, fmt->double_format, (double) vec[vec_len-1]) >= 0, "printing failed");
    libcheck(fprintf(file, "%s%s", fmt->vend_bracket, postfix) >= 0, "printing failed");
    return LT_SUCCESS;

error:
    return LT_FAILURE;
}

int print_lvector_binary(FILE* file, const long* vec, size_t vec_len,
                        const PRINTING_FMT* fmt __attribute__((__unused__)))
{ return fwrite(vec, sizeof(long), vec_len, file) == vec_len ? 0 : -1; }

int print_lvector_binary_float(FILE* file, const long* vec, size_t vec_len,
                        const PRINTING_FMT* fmt __attribute__((__unused__)))
{
    size_t buf_size = vec_len < 1024 ? vec_len : 1024;
    size_t left = vec_len;
    double buffer[buf_size];

    while(left)
    {
        size_t write_size = left < buf_size ? left : buf_size;
        for(size_t i = 0; i < write_size; i++)
            buffer[i] = (double) vec[i];

        libcheck(fwrite(buffer, sizeof(double), write_size, file)
                == write_size, "fwrite failed");
        left -= write_size;
    }
    return LT_SUCCESS;

error:
    return LT_FAILURE;
}

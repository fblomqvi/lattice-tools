/* print_util.h
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

#ifndef FB_LATTICE_TOOLS_PRINT_UTIL_H
#define FB_LATTICE_TOOLS_PRINT_UTIL_H

#include <stdio.h>

enum PrintingFmt
{
    PRINTING_FMT_DEFAULT = 0,
    PRINTING_FMT_FPLLL,
    PRINTING_FMT_MATLAB,
    PRINTING_FMT_MATHEMATICA,
    PRINTING_FMT_MAX
};

typedef struct s_printing_format
{
    const char* mstart_bracket;
    const char* mend_bracket;
    const char* vstart_bracket;
    const char* vend_bracket;
    const char* element_separator;
    const char* vector_separator;
    const char* double_format;
    const char* integer_format;
} PRINTING_FMT;


const char* printing_fmt_get_name(enum PrintingFmt fmt);
const PRINTING_FMT* printing_fmt_get(enum PrintingFmt fmt);
int printing_fmt_parse_name(const char* name);
int printing_fmt_print_names(FILE* file);

int print_dvector(FILE* file, const double* vec, 
                size_t vec_len, const PRINTING_FMT* fmt, 
                const char* prefix, const char* postfix);

int print_dvector_binary(FILE* file, const double* vec, size_t vec_len, 
                        const PRINTING_FMT* fmt __attribute__((__unused__)));

static inline int print_dvector_cword(FILE* file, const double* vec, 
                                    size_t vec_len, const PRINTING_FMT* fmt)
{ return print_dvector(file, vec, vec_len, fmt, "", " --> "); }

static inline int print_dvector_std(FILE* file, const double* vec, 
                                    size_t vec_len, const PRINTING_FMT* fmt)
{ return print_dvector(file, vec, vec_len, fmt, "", "\n"); }


int print_lvector(FILE* file, const long* vec,
                size_t vec_len, const PRINTING_FMT* fmt,
                const char* prefix, const char* postfix);

int print_lvector_binary(FILE* file, const long* vec, size_t vec_len, 
                        const PRINTING_FMT* fmt __attribute__((__unused__)));

int print_lvector_binary_float(FILE* file, const long* vec, size_t vec_len,
                        const PRINTING_FMT* fmt __attribute__((__unused__)));

static inline int print_lvector_std(FILE* file, const long* vec, 
                                    size_t vec_len, const PRINTING_FMT* fmt)
{ return print_lvector(file, vec, vec_len, fmt, "", "\n"); }

int print_lvector_as_double(FILE* file, const long* vec,
                            size_t vec_len, const PRINTING_FMT* fmt,
                            const char* prefix, const char* postfix);

static inline int print_lvector_as_double_std(FILE* file, const long* vec, 
                                            size_t vec_len, const PRINTING_FMT* fmt)
{ return print_lvector_as_double(file, vec, vec_len, fmt, "", "\n"); }

#endif /* FB_LATTICE_TOOLS_PRINT_UTIL_H */

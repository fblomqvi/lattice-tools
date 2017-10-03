/* Macros for debugging and error handling.
   Copyright (C) 2016 Ferdinand Blomqvist

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

/* These debug macros are inspired by the debug macros found in chapter 20 of
   'Learn C the hard way' by Zed A. Shaw
   (http://c.learncodethehardway.org/book/ex20.html). I have tweaked and
   expanded them to my liking. Furthermore I have fixed some possible
   'swallowing the semicolon' issues by adding 'do {...} * while(0)' statements
   to the right places. */

#ifndef FB_UTILITY_DBG_H
#define FB_UTILITY_DBG_H

#include "prog_name.h"
#include <stdio.h>
#include <errno.h>
#include <string.h>

/* A convenience function that adds ': <strerror>\n' to the output specified by
 * the format string and the additional arguments. */
void fprintf_we(FILE* file, const char* format, ...);
void fprintf_we_lt(FILE* file, const int lt_errno, const char* format, ...);

#ifdef NDEBUG
#define debug(msg, ...)
#define debug_err(msg, ...)
#define debugf(file, msg, ...)
#define debug_do(statement)
#else /* NDEBUG */
#ifdef DBG_NO_PROG_NAME
#define debug(msg, ...) \
    fprintf(stderr, "[DEBUG] (%s:%d) " msg "\n",\
             __FILE__, __LINE__, ##__VA_ARGS__)

#define debug_err(msg, ...) \
    fprintf_we(stderr, "[DEBUG] (%s:%d) " msg,\
             __FILE__, __LINE__, ##__VA_ARGS__)

#define debugf(file, msg, ...) \
    fprintf(file, "[DEBUG] (%s:%d) " msg "\n",\
             __FILE__, __LINE__, ##__VA_ARGS__)

#else /* DBG_NO_PROG_NAME not defined */
#define debug(msg, ...) \
    fprintf(stderr, "%s: [DEBUG] (%s:%d) " msg "\n",\
             PROGRAM_NAME, __FILE__, __LINE__, ##__VA_ARGS__)

#define debug_err(msg, ...) \
    fprintf_we(stderr, "%s: [DEBUG] (%s:%d) " msg,\
             PROGRAM_NAME, __FILE__, __LINE__, ##__VA_ARGS__)

#define debugf(file, msg, ...) \
    fprintf(file, "%s: [DEBUG] (%s:%d) " msg "\n",\
             PROGRAM_NAME, __FILE__, __LINE__, ##__VA_ARGS__)

#endif /* DBG_NO_PROG_NAME */
#define debug_do(statement) statement
#endif /* NDEBUG */


#ifdef DBG_NO_PROG_NAME
#define log_err(msg, ...) \
    fprintf_we(stderr, "[ERROR] (%s:%d) " msg, \
            __FILE__, __LINE__, ##__VA_ARGS__)

#define log_err_ne(msg, ...) \
    fprintf(stderr, "[ERROR] (%s:%d) " msg "\n", \
            __FILE__, __LINE__, ##__VA_ARGS__)

#define log_warn(msg, ...) \
    fprintf_we(stderr, "[WARN] (%s:%d) " msg, \
            __FILE__, __LINE__, ##__VA_ARGS__)

#define log_info(msg, ...) \
    fprintf(stderr, "[INFO] (%s:%d) " msg "\n", \
            __FILE__, __LINE__, ##__VA_ARGS__)

#define log_plain(msg, ...) \
    fprintf(stderr, msg "\n", ##__VA_ARGS__)

#define log_err_lt(lt_errno, msg, ...) \
    fprintf_we_lt(stderr, lt_errno, "[ERROR] (%s:%d) " msg, \
            __FILE__, __LINE__, ##__VA_ARGS__)
#else

/* The log_<name> macros need the PROGRAM_NAME macro to be defined so we set it
 * to a dummy program name if it is unset. */
#ifndef PROGRAM_NAME 
#define PROGRAM_NAME "dummy-prog-name"
#endif /* PROGRAM_NAME */

#define log_err(msg, ...) \
    fprintf_we(stderr, "%s: [ERROR] (%s:%d) " msg, \
            PROGRAM_NAME, __FILE__, __LINE__, ##__VA_ARGS__)

#define log_err_ne(msg, ...) \
    fprintf(stderr, "%s: [ERROR] (%s:%d) " msg "\n", \
            PROGRAM_NAME, __FILE__, __LINE__, ##__VA_ARGS__)

#define log_warn(msg, ...) \
    fprintf_we(stderr, "%s: [WARN] (%s:%d) " msg, \
            PROGRAM_NAME, __FILE__, __LINE__, ##__VA_ARGS__)

#define log_info(msg, ...) \
    fprintf(stderr, "%s: [INFO] (%s:%d) " msg "\n", \
            PROGRAM_NAME, __FILE__, __LINE__, ##__VA_ARGS__)

#define log_plain(msg, ...) \
    fprintf(stderr, "%s: " msg "\n", PROGRAM_NAME, ##__VA_ARGS__)

#define log_err_lt(lt_errno, msg, ...) \
    fprintf_we_lt(stderr, lt_errno, "%s: [ERROR] (%s:%d) " msg, \
            PROGRAM_NAME, __FILE__, __LINE__, ##__VA_ARGS__)
#endif /* DBG_NO_PROG_NAME */


/* Checks if 'c' is true or false. If 'c' is false, then the error is logged by
 * 'print_func', errno is cleared, and the execution jumps to 'label'. */
#define lcheck_pf(c, print_func, label, msg, ...) \
    do {\
        if(!(c)) \
        {\
            print_func(msg, ##__VA_ARGS__);\
            errno = 0;\
            goto label;\
        }\
    } while(0)

/* Library version of lcheck_pf. Exactly like lcheck_pf(), except that errno is
 * not set to zero. */
#define llibcheck_pf(c, print_func, label, msg, ...) \
    do {\
        if(!(c)) \
        {\
            print_func(msg, ##__VA_ARGS__);\
            goto label;\
        }\
    } while(0)

/* Checks if 'c' is true or false. If 'c' is false, then the error is logged by
 * 'print_func', errno is cleared, and the execution jumps to 'label'. */
#define lcheck_lt_pf(c, print_func, label, lt_errno, msg, ...) \
    do {\
        if(!(c)) \
        {\
            print_func(lt_errno, msg, ##__VA_ARGS__);\
            goto label;\
        }\
    } while(0)

#define llibcheck_se_pf(c, print_func, label, lt_errno, error_code) \
    do {\
        if(!(c)) \
        {\
            lt_errno = error_code;\
            print_func("%s", lt_strerror(error_code));\
            goto label;\
        }\
    } while(0)

/* Convenience macros for different printing functions. */
#define lcheck(c, label, msg, ...) lcheck_pf(c, log_err, label, msg, ##__VA_ARGS__)
#define lcheck_dbg(c, label, msg, ...) lcheck_pf(c, debug_err, label, msg, ##__VA_ARGS__)
#define llibcheck(c, label, msg, ...) llibcheck_pf(c, debug_err, label, msg, ##__VA_ARGS__)

/* Convenience macros that uses the standardized 'error' label. */
#define check(c, msg, ...) lcheck(c, error, msg, ##__VA_ARGS__)
#define check_dbg(c, msg, ...) lcheck_dbg(c, error, msg, ##__VA_ARGS__)
#define libcheck(c, msg, ...) llibcheck(c, error, msg, ##__VA_ARGS__)

/* Convenience macros for memory allocation checking. */
#define lcheck_pf_mem(c, print_func, label) \
    lcheck_pf(c, print_func, label, "Memory allocation error")

#define lcheck_mem(c, label) lcheck_pf_mem(c, log_err, label)
#define lcheck_dbg_mem(c, label) lcheck_pf_mem(c, debug_err, label)
#define llibcheck_mem(c, label) \
    llibcheck_pf(c, debug_err, label, "Memory allocation error")

#define check_mem(c) lcheck_mem(c, error)
#define check_dbg_mem(c) lcheck_dbg_mem(c, error)
#define libcheck_mem(c) llibcheck_mem(c, error)

/* More conveniece macros */
#define lcheck_lt(c, label, lt_errno, msg, ...) \
    lcheck_lt_pf(c, log_err_lt, label, lt_errno, msg, ##__VA_ARGS__)
#define check_lt(c, lt_errno, msg, ...) \
    lcheck_lt_pf(c, log_err_lt, error, lt_errno, msg, ##__VA_ARGS__)

#define llibcheck_se(c, label, lt_errno, error_code) \
    llibcheck_se_pf(c, debug, label, lt_errno, error_code)

//#define check_fatal(c) check(c, "Fatal error")

#define lcheck_parse(c, label, line_num, line, msg, ...) \
    lcheck_pf(c, log_err_ne, label, \
        "Parsing error on line %zu: '%s'; " msg, \
        line_num, line, ##__VA_ARGS__)

#define check_parse(c, line_num, line, msg, ...) \
    lcheck_parse(c, error, line_num, line, msg, ##__VA_ARGS__)

#endif /* FB_UTILITY_DBG_H */

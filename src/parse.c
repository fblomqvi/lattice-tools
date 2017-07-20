/* Various parsing functions.
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
#include "parse.h"
#include <limits.h>
#include <ctype.h>
//#include <sys/stat.h>
#include <assert.h>
//#include <math.h>

#define printable_input(str) (is_printable(str) ? str : "<non-printable input>")

/*
typedef struct s_parser_state
{
    size_t next_edge;
    size_t next_site;
    size_t line_num;
    int eof;
    size_t target_site_degree;
    size_t target_check_degree;
} STATE;
*/

/*
static void print_error_msg(int eof, size_t line_number)
{
    if(eof)
        log_err("End of input file reached at line %zu", line_number);
    else
        log_err("Error in input file on line %zu", line_number);
}
*/

static int is_printable(const char* str)
{ 
    while(*str != '\0')
    {
        if(!isprint(*str))
            return 0;
        str++;
    }
    return 1;
}

/*
static int is_comment(char* ptr, char comment_start)
{ return *ptr == comment_start; }

static char* skip_line(char* ptr, STATE* state)
{
    while(*ptr != '\0' && *ptr++ != '\n');
    if(*ptr != '\0')
        state->line_num++;
    return ptr;
}

static void find_eol(char* line, char** next_line, int* eof)
{
    while(*line != '\0')
    {
        if(*line == '\n')
        {
            *line = '\0';
            *next_line = line + 1;
            return;
        }
        line++;
    }
    *eof = 1;    
}
*/

/*
static int read_size_t_line(char* ptr, char** endptr, size_t* n, STATE* state)
{
    char* line = ptr;
    find_eol(ptr, endptr, &state->eof);

    libcheck(isdigit(*ptr), "syntax error, expected digit");

    *n = strtoul(ptr, &ptr, 10);
    libcheck(*ptr == '\0' && !(errno == ERANGE && *n == ULONG_MAX), 
                "syntax error or input is not an integer");

    state->line_num++;
    return 0;

error:
    log_err_ne("Parsing error on line %zu: '%s'; expected '%%zu'", 
            state->line_num, printable_input(line));
    return -1;
}

static int read_size_t_lines(char* ptr, char** endptr, size_t* n, size_t num, STATE* state)
{
    for(size_t i = 0; i < num; i++)
    {
        int rc = read_size_t_line(ptr, &ptr, n + i, state);
        libcheck(rc == 0, "read_size_t_line failed");
        //debug("n: %zu", n);
    }
    *endptr = ptr;
    return 0;

error:
    return -1;
}
*/

/* Reads at most buf_len - 1 non-whitespace characters from the stream pointed
 * to by file. Note that buf_len must be non-zero */
static int read_until_space(FILE* file, char* buffer, size_t buf_len)
{
    assert(buf_len > 0);

    int ch;
    size_t read_len = 0;
    
    // We skip initial whitespace.
    do {
        ch = fgetc(file);
        check(ch != EOF, "Error parsing: %s", ferror(file) ? 
                "stream error" : "end of file reached");
    } while(isspace(ch));

    do {
        check(read_len < buf_len - 1, "Error parsing: "
                "syntax error, input to long for buffer");

        buffer[read_len++] = ch;
        ch = fgetc(file);
        check(ch != EOF, "Error parsing: %s", ferror(file) ? 
                "stream error" : "end of file reached");
    } while(!isspace(ch));

    buffer[read_len] = '\0';
    return 0;

error:
    return EOF;
}

/* Reads at most buf_len - 1 non-whitespace characters from the stream pointed
 * to by file. Note that buf_len must be non-zero */
/*
static int read_until_newline(FILE* file, char* buffer, size_t buf_len)
{
    assert(buf_len > 0);

    int ch;
    size_t read_len = 0;
    
    // We skip initial whitespace.
    do {
        ch = fgetc(file);
        if(ch == EOF)
            return 1;
    } while(isspace(ch));

    int eof = 0;
    do {
        check(read_len < buf_len - 1, "Error parsing: "
                "syntax error, input to long for buffer");

        buffer[read_len++] = ch;
        ch = fgetc(file);
        if(ch == EOF)
        {
            eof = 1;
            break;
        }
    } while(ch != '\n');

    buffer[read_len] = '\0';
    return eof ? 1 : 0;

error:
    return EOF;
}
*/

static int read_degree(FILE* file, long* degree)
{
    char* endptr;
    char buffer[24];
    int rc = read_until_space(file, buffer, sizeof(buffer));
    libcheck(rc == 0, "read_until_space failed");
    *degree = strtol(buffer, &endptr, 10);
    check(*endptr == '\0' && 
            !(errno == ERANGE && (*degree == LONG_MAX || *degree == LONG_MIN)), 
            "Error parsing '%s'; expected '%%zu'", printable_input(buffer));
    return 0;

error:
    return -1;
}
/* returns zero on success, and EOF in case of an error or EOF */
static int skip_comment_lines(FILE* file, char comment_start)
{
    while(1)
    {
        int ch = fgetc(file);
        check(ch != EOF, "Error parsing: %s", ferror(file) ? 
                "stream error" : "end of file reached");

        if(ch == comment_start)
        {
            do {
                ch = fgetc(file);
                check(ch != EOF, "Error parsing: %s", ferror(file) ? 
                        "stream error" : "end of file reached");
            } while(ch != '\n');
        }
        else
        {
            int rc = ungetc(ch, file);
            check(rc == ch, "Error parsing: ungetc failed");
            break;
        }
    }
    return 0;

error:
    return EOF;
}

int parse_hyperplane_coeffs(FILE* file, long* coeffs, size_t num_coeffs)
{
    int rc = skip_comment_lines(file, '#');
    libcheck(rc == 0, "skip_comment_lines failed");

    for(size_t i = 0; i < num_coeffs; i++)
    {
        rc = read_degree(file, coeffs + i);
        libcheck(rc == 0, "failed to read degree sequence: read_degree failed");
    }
    return 0;

error:
    return -1;
}

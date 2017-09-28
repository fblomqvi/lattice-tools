/* configuration.h
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

#ifndef FB_UTILITY_CONFIGURATION_H
#define FB_UTILITY_CONFIGURATION_H

#include <stdlib.h>
#include <stdio.h>
#include "version.h"

union value
{
    int i;
    long l;
    unsigned long u;
    size_t s;
    const char* str;
    int b;
};

enum value_type
{
    type_int,
    type_long,
    type_ulong,
    type_size,
    type_str,
    type_bool
};

struct config
{
    const char* name;
    union value val;
    enum value_type type;
    int print_flag;
};

void util_get_current_time(char* timestr, size_t max, const char* format);

int util_print_genby_and_config(FILE* file, const char* prog_name,
                                const char* tstr, struct version* ver,
                                struct config* conf);

#endif /* FB_UTILITY_CONFIGURATION_H */

/* version.c
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
#include <stdio.h>
#include <stdlib.h>

const struct version g_current_version = { VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH };

int print_version(FILE* file)
{
    static const char* formatstr =
        "%s (%s) version %u.%u.%u (released %u-%.2u-%.2u, compiled %s %s)\n"
        "Copyright (C) 2017 Ferdinand Blomqvist\n"
        "License GPLv2: GNU GPL version 2 <http://gnu.org/licenses/gpl.html>.\n"
        "This is free software: you are free to change and redistribute it.\n"
        "There is NO WARRANTY, to the extent permitted by law.\n\n"
        "Written by Ferdinand Blomqvist.\n";
    return (fprintf(file, formatstr, PROGRAM_NAME, PACKAGE_NAME, VERSION_MAJOR,
            VERSION_MINOR, VERSION_PATCH, RELEASE_YEAR, RELEASE_MONTH,
            RELEASE_DAY, __DATE__, __TIME__) < 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}

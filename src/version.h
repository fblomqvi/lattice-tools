/* version.h
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

#ifndef FB_LIBLATTICE_VERSION_H
#define FB_LIBLATTICE_VERSION_H

#define VERSION_MAJOR 0
#define VERSION_MINOR 1
#define VERSION_PATCH 0

#define RELEASE_YEAR 2017
#define RELEASE_MONTH 7
#define RELEASE_DAY 18

#define PACKAGE_NAME "lattice-tools"

struct version
{
    unsigned major;
    unsigned minor;
    unsigned patch;
};

extern const struct version g_current_version;

int print_version(FILE* file);

#endif /* FB_LIBLATTICE_VERSION_H */

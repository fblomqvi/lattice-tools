/* lt_errno.h
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

#ifndef FB_LATTICE_TOOLS_LT_ERRNO_H
#define FB_LATTICE_TOOLS_LT_ERRNO_H

#include <errno.h>

enum 
{ 
    LT_SUCCESS = 0, 
    LT_FAILURE = -1,
    LT_ELINDEP = 1024,   /* Basis vectors are linearly dependent */
    LT_ENOMEM,          /* malloc failed */
};

const char* lt_strerror(const int lt_errno);

#endif /* FB_LATTICE_TOOLS_LT_ERRNO_H */

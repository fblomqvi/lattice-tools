/* lt_errno.c
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

#include "lt_errno.h"
#include <string.h>
#include <errno.h>
#include <gsl/gsl_errno.h>

const char* lt_strerror(const int lt_errno)
{
    switch(lt_errno)
    {
        case LT_SUCCESS:
            return "success";
        case LT_FAILURE:
            return "failure";
            /*
        case LT_ELINDEP:
            return "the basis vectors are linearly dependent";
        case LT_ENOMEM:
            return "malloc failed";
            */
        case LT_EINVAL:
            return "invalid argument supplied by user";
        case LT_ESYSTEM:
            return strerror(errno);
        default:
            return gsl_strerror(lt_errno);
    }
}

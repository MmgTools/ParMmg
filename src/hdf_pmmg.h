/* =============================================================================
**  This file is part of the parmmg software package for parallel tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux, 2017-
**
**  parmmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  parmmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with parmmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the parmmg distribution only if you accept them.
** =============================================================================
*/

#ifndef HDF_PMMG_H
#define HDF_PMMG_H
/**
 * \file hdf_pmmg.h
 * \brief HDF5 tools that are used in different part of parMmg
 * \author Gabriel Suau (Bx INP/Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 */
#include "hdf5.h"

#define HDF_CHECK(func_call,on_failure) do {                            \
    int hdf_ret_val;                                                    \
                                                                        \
    hdf_ret_val = func_call;                                            \
                                                                        \
    if ( hdf_ret_val < 0 ) {                                            \
      fprintf(stderr," ## Error: %s:%d: HDF5 error\n.",                 \
              __func__,__LINE__ );                                      \
      on_failure;                                                       \
    }                                                                   \
  } while(0)

#endif

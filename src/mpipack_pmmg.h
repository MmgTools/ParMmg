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

#ifndef MPIPACK_PMMG_H
#define MPIPACK_PMMG_H
/**
 * \file mpipack_pmmg.h
 * \brief header of packing functions (pack into a char buffer)
 * \author Luca Cirrottola (Inria)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 */
#include "libmmgtypes.h"

int PMMG_mpisizeof_grp ( PMMG_pGrp grp );
int PMMG_mpisizeof_grp4finalmerge ( PMMG_pGrp grp );
int PMMG_mpipack_grp ( PMMG_pGrp grp,char **buffer );
int PMMG_mpipack_grp4finalmerge ( PMMG_pGrp grp,char **buffer );

#endif

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

#ifndef MPIUNPACK_PMMG_H
#define MPIUNPACK_PMMG_H
/**
 * \file mpiunpack_pmmg.h
 * \brief header of unpacking functions (unpack from a char buffer)
 * \author Luca Cirrottola (Inria)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 */
#include "libmmgtypes.h"

int PMMG_mpiunpack_grp ( PMMG_pParMesh,PMMG_pGrp,int,char **buffer);

int PMMG_mpiunpack_parmesh ( PMMG_pParMesh,PMMG_pGrp,int,PMMG_pInt_comm,int*,
                             PMMG_pExt_comm*,char** );

#endif

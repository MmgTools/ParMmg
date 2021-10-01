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

/**
 * \file coorcell_pmmg.h
 * \brief coorcell_pmmg.c header file
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#ifndef COORCELL_PMMG_H

#define COORCELL_PMMG_H

#include "parmmg.h"

/* numerical accuracy for coordinate checks */
#define PMMG_EPSCOOR  1.e-10
#define PMMG_EPSCOOR2 3.e-20

/**
 * \struct PMMG_coorCell
 *
 * \brief Cell containing a point coordinate and an index
 *
 */
typedef struct {
  double  c[3]; /*!< point coordinates */
  int     idx;  /*!< index associated to the point */
  int     grp;  /*!< a group to which belong the point */
} PMMG_coorCell;


int PMMG_compare_coorCell (const void * a, const void * b);
int PMMG_find_coorCellListBoundingBox(PMMG_coorCell*,int,double*,double*,double*);
int PMMG_scale_coorCellList (PMMG_coorCell*,int,double*,double*,double*);
int PMMG_unscale_coorCellList (PMMG_coorCell*,int,double*,double*,double);

#endif

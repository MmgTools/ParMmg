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
 * \file interpmesh_pmmg.h
 * \brief Interpolate data from a background mesh to the current mesh.
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \author Luca Cirrottola (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#ifndef INTERPMESH_PMMG_H

#define INTERPMESH_PMMG_H

#include "locate_pmmg.h"

int PMMG_interp4bar_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTetra pt,int,PMMG_barycoord*);
int PMMG_interp4bar_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTetra pt,int,PMMG_barycoord*);
int PMMG_interp3bar_iso(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTria ptr,int,PMMG_barycoord*);
int PMMG_interp3bar_ani(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTria ptr,int,PMMG_barycoord*);
int PMMG_interp2bar_iso( MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTria ptr,int ip,int l,PMMG_barycoord *barycoord );
int PMMG_interp2bar_ani( MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTria ptr,int ip,int l,PMMG_barycoord *barycoord );

extern int (*PMMG_interp4bar)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTetra pt,int,PMMG_barycoord*);
extern int (*PMMG_interp3bar)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTria ptr,int,PMMG_barycoord*);
extern int (*PMMG_interp2bar)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTria ptr,int ip,int l,PMMG_barycoord *barycoord);

#endif

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
 * \file locate_pmmg.h
 * \brief Point localization for interpolation on a new mesh.
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \author Luca Cirrottola (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#ifndef LOCATE_PMMG_H

#define LOCATE_PMMG_H

/** \struct PMMG_baryCoord
 *
 * \brief Struct containing the index and value of a barycentric coordinate
 *
 */
typedef struct {
  int    idx; /*!< direction */
  double val; /*!< coordinate value */
} PMMG_baryCoord;

double PMMG_quickarea(double *a,double *b,double *c,double *n);
void PMMG_get_baryCoord( double *val,PMMG_baryCoord *phi,int ndim );
int PMMG_compare_baryCoord( const void *a,const void *b );
int PMMG_compute_baryCoord2d( MMG5_pMesh mesh,MMG5_pTria ptr,int k,double *coord,
                              double *normal,PMMG_baryCoord *barycoord );
int PMMG_compute_baryCoord3d( MMG5_pMesh mesh,MMG5_pTetra pt,double *coord,
                              double *faceAreas, PMMG_baryCoord *barycoord );
int PMMG_locatePointInClosestTria( MMG5_pMesh mesh,int k,MMG5_pPoint ppt,
                                   PMMG_baryCoord *barycoord );
int PMMG_locatePointInTria( MMG5_pMesh mesh,MMG5_pTria ptr,int k,MMG5_pPoint ppt,
                            double *triaNormal,PMMG_baryCoord *barycoord,
                            double *closestDist,int *closestTria );
int PMMG_locatePointInTetra( MMG5_pMesh mesh,MMG5_pTetra pt,MMG5_pPoint ppt,
                             double *faceAreas,PMMG_baryCoord *barycoord );
int PMMG_locatePointBdy( MMG5_pMesh mesh,MMG5_pPoint ppt,int init,
                         double *triaNormals,PMMG_baryCoord *barycoord,
                         int ip,int igrp );
int PMMG_locatePointVol( MMG5_pMesh mesh,MMG5_pPoint ppt,int init,
                         double *faceAreas,PMMG_baryCoord *barycoord,
                         int ip,int igrp );
void PMMG_locate_setStart( MMG5_pMesh mesh,MMG5_pMesh meshOld );

#endif

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

/** \struct PMMG_barycoord
 *
 * \brief Struct containing the index and value of a barycentric coordinate
 *
 */
typedef struct {
  int    idx; /*!< direction */
  double val; /*!< coordinate value */
} PMMG_barycoord;

/** \struct PMMG_locateStats
 *
 * \brief Struct containing the statistics of localization searches
 *
 */
typedef struct {
  int    nexhaust; /*!< number of exhaustive searches */
  int    stepmax;  /*!< maximum number of steps on the search paths */
  int    stepmin;  /*!< minimum number of steps on the search paths */
  double stepav;   /*!< average number of steps on the search paths */
} PMMG_locateStats;


double PMMG_quickarea(double *a,double *b,double *c,double *n);
void PMMG_barycoord_get( double *val,PMMG_barycoord *phi,int ndim );
int  PMMG_barycoord_compare( const void *a,const void *b );
int  PMMG_barycoord2d_compute( MMG5_pMesh mesh,MMG5_pTria ptr,int k,double *coord,
                               double *normal,PMMG_barycoord *barycoord );
int  PMMG_barycoord3d_compute( MMG5_pMesh mesh,MMG5_pTetra pt,double *coord,
                               double *faceAreas, PMMG_barycoord *barycoord );
int PMMG_locatePointInTria( MMG5_pMesh mesh,MMG5_pTria ptr,int k,MMG5_pPoint ppt,
                            double *triaNormal,PMMG_barycoord *barycoord,
                            double *closestDist,int *closestTria );
int PMMG_locatePointInTetra( MMG5_pMesh mesh,MMG5_pTetra pt,MMG5_pPoint ppt,
                             double *faceAreas,PMMG_barycoord *barycoord );
int PMMG_locatePointBdy( MMG5_pMesh mesh,MMG5_pPoint ppt,int init,
                         double *triaNormals,PMMG_barycoord *barycoord,
                         int ip );
int PMMG_locatePointVol( MMG5_pMesh mesh,MMG5_pPoint ppt,int init,
                         double *faceAreas,PMMG_barycoord *barycoord,
                         int ip );
int PMMG_locatePoint_errorCheck( MMG5_pMesh mesh,int ip,int kfound,int myrank,int igrp );
void PMMG_locate_setStart( MMG5_pMesh mesh,MMG5_pMesh meshOld );
void PMMG_locate_postprocessing( MMG5_pMesh mesh,MMG5_pMesh meshOld,PMMG_locateStats *locStats );
void PMMG_locate_print( PMMG_locateStats *locStats,int ngrp,int myrank );

#endif

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
 * \file barycoord_pmmg.h
 * \brief Barycentric coordinates for point localization in a mesh.
 * \author Luca Cirrottola (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#ifndef BARYCOORD_PMMG_H

#define BARYCOORD_PMMG_H

/** \struct PMMG_barycoord
 *
 * \brief Struct containing the index and value of a barycentric coordinate
 *
 */
typedef struct {
  int    idx; /*!< direction */
  double val; /*!< coordinate value */
} PMMG_barycoord;


double PMMG_quickarea(double *a,double *b,double *c,double *n);
void PMMG_barycoord_get( double *val,PMMG_barycoord *phi,int ndim );
int  PMMG_barycoord_compare( const void *a,const void *b );
int  PMMG_barycoord2d_compute( MMG5_pMesh mesh,MMG5_pTria ptr,int k,double *coord,
                               double *normal,PMMG_barycoord *barycoord );
double PMMG_barycoord2d_compute1( MMG5_pMesh mesh,MMG5_pTria ptr,int k,int ia,
                                  double *proj,double *normal );
int PMMG_barycoord2d_project( MMG5_pMesh mesh,MMG5_pTria ptr,double *coord,
                                 double *proj,double dist,double *normal );
int  PMMG_barycoord3d_compute( MMG5_pMesh mesh,MMG5_pTetra pt,double *coord,
                               double *faceAreas, PMMG_barycoord *barycoord );
int  PMMG_barycoord2d_evaluate( MMG5_pMesh mesh,MMG5_pTria ptr,int k,
                                double *coord,double *triaNormal,
                                PMMG_barycoord *barycoord );
int  PMMG_barycoord3d_evaluate( MMG5_pMesh mesh,MMG5_pTetra pt,
                                double *coord,double *faceAreas,
                                PMMG_barycoord *barycoord );

#endif

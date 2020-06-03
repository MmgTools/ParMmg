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
#include "parmmg.h"
#include "barycoord_pmmg.h"

/**
 * \param a coordinates of the first point of face.
 * \param b coordinates of the second point of face.
 * \param c coordinates of the third point of face.
 * \param n pointer to store the computed normal.
 * \return 1
 *
 * Compute triangle area given three points on the surface.
 *
 */
double PMMG_quickarea(double *a,double *b,double *c,double *n) {
  double        area[3],abx,aby,abz,acx,acy,acz;

  /* area */
  abx = b[0] - a[0];
  aby = b[1] - a[1];
  abz = b[2] - a[2];

  acx = c[0] - a[0];
  acy = c[1] - a[1];
  acz = c[2] - a[2];

  area[0] = aby*acz - abz*acy;
  area[1] = abz*acx - abx*acz;
  area[2] = abx*acy - aby*acx;

  return area[0]*n[0]+area[1]*n[1]+area[2]*n[2];
}

/**
 * \param coord pointer to a double array
 * \param barycoord pointer to the point barycentric coordinates in the current
 * tetra
 * \param ndim space dimension of the manifold
 *
 *  Get barycentric coordinates.
 *
 */
void PMMG_barycoord_get( double *val,PMMG_barycoord *phi,int ndim ) {
  int i;

  for( i = 0; i < ndim; i++ )
    val[phi[i].idx] = phi[i].val;

}

/**
 * \param a pointer to point barycentric coordinates
 * \param b pointer to point barycentric coordinates
 *
 * \return -1 if (a < b), +1 if (a > b), 0 if equal
 *
 *  Compare the barycentric coordinates of a given point in a given tetrahedron.
 *
 */
int PMMG_barycoord_compare( const void *a,const void *b ) {
  PMMG_barycoord *coord_a;
  PMMG_barycoord *coord_b;

  coord_a = (PMMG_barycoord *)a;
  coord_b = (PMMG_barycoord *)b;

  if( coord_a->val > coord_b-> val ) return 1;
  if( coord_a->val < coord_b-> val ) return -1;

  return 0;
}

int PMMG_barycoord_isInside( PMMG_barycoord *phi ) {
  if( phi[0].val > -MMG5_EPS )
    return 1;
  else
    return 0;
}

int PMMG_barycoord_isBorder( PMMG_barycoord *phi,int *ifoundEdge,int *ifoundVertex ) {
  if( phi[0].val < MMG5_EPS )
    if( phi[1].val < MMG5_EPS )
      *ifoundVertex = phi[2].idx;
    else
      *ifoundEdge = phi[0].idx;

  return 1;
}

/**
 * \param mesh pointer to the mesh structure
 * \param ptr pointer to the current triangle
 * \param coord coordinates of the point to project
 * \param proj tangentially projected point coordinates
 * \param dist orthogonal distance of the point
 * \param normal unit normal of the current triangle
 *
 * \return 0 if fail, 1 if success
 *
 *  Compute tangent and normal projection on triangle..
 *
 */
int PMMG_barycoord2d_project( MMG5_pMesh mesh,MMG5_pTria ptr,double *coord,
                                 double *proj,double dist,double *normal ) {
  double *c0;
  int    d;

  /* Project point on the triangle plane */
  c0 = mesh->point[ptr->v[0]].c;

  dist = 0.0;
  for( d = 0; d < 3; d++ )
    dist += (coord[d]-c0[d])*normal[d];

  for( d = 0; d < 3; d++ )
    proj[d] = coord[d] - dist*normal[d];

  return 1;
}

/**
 * \param mesh pointer to the mesh structure
 * \param ptr pointer to the current triangle
 * \param k index of the triangle
 * \param ia edge index
 * \param normal unit normal of the current triangle
 *
 * \return 0 if fail, 1 if success
 *
 *  Compute the barycentric coordinates of a given point in a given triangle
 *  only with respect to a given edge.
 *
 */
double PMMG_barycoord2d_compute1( MMG5_pMesh mesh,MMG5_pTria ptr,int k,int ia,
                                  double *proj,double *normal ) {
  double *c1,*c2;

  /* Retrieve face areas and compute barycentric coordinate */
  c1 = mesh->point[ptr->v[MMG5_inxt2[ia]]].c;
  c2 = mesh->point[ptr->v[MMG5_inxt2[ia+1]]].c;

  return PMMG_quickarea( proj, c1, c2, normal )/ptr->qual;
}

/**
 * \param mesh pointer to the mesh structure
 * \param ptr pointer to the current triangle
 * \param k index of the triangle
 * \param coord pointer to the point coordinates
 * \param normal unit normal of the current triangle
 * \param barycoord pointer to the point barycentric coordinates in the current
 * tetra
 *
 * \return 0 if fail, 1 if success
 *
 *  Compute the barycentric coordinates of a given point in a given triangle.
 *
 */
int PMMG_barycoord2d_compute( MMG5_pMesh mesh,MMG5_pTria ptr,int k,double *coord,
                              double *normal,PMMG_barycoord *barycoord ) {
  double dist,proj[3],*c1,*c2,vol;
  int    ia,i;

  /* Project point on the triangle plane */
  c1 = mesh->point[ptr->v[0]].c;

  dist = 0.0;
  for( i = 0; i < 3; i++ )
    dist += (coord[i]-c1[i])*normal[i];

  for( i = 0; i < 3; i++ )
    proj[i] = coord[i] - dist*normal[i];

  /* Retrieve tria area */
  vol = ptr->qual;

  /* Retrieve face areas and compute barycentric coordinates */
  for( ia = 0; ia < 3; ia++ ) {
    c1 = mesh->point[ptr->v[MMG5_inxt2[ia]]].c;
    c2 = mesh->point[ptr->v[MMG5_inxt2[ia+1]]].c;

    barycoord[ia].val = PMMG_quickarea( proj, c1, c2, normal )/vol;
    barycoord[ia].idx = ia;
  }

  /* Store normal distance in the third coordinate */
  barycoord[3].val = dist;
  barycoord[3].idx = 3;

  return 1;
}

/**
 * \param mesh pointer to the mesh structure
 * \param pt pointer to the current tetra
 * \param coord pointer to the point coordinates
 * \param faceAreas oriented face areas of the current tetrahedron
 * \param barycoord pointer to the point barycentric coordinates in the current
 * tetra
 *
 * \return 0 if fail, 1 if success
 *
 *  Compute the barycentric coordinates of a given point in a given tetrahedron.
 *
 */
int PMMG_barycoord3d_compute( MMG5_pMesh mesh,MMG5_pTetra pt,double *coord,
                              double *faceAreas,PMMG_barycoord *barycoord ) {
  double *c0,*normal,vol;
  int    ifac;

  /* Retrieve tetra volume */
  vol = pt->qual;

  /* Retrieve face areas and compute barycentric coordinates */
  for( ifac = 0; ifac < 4; ifac++ ) {
    normal = &faceAreas[3*ifac];
    c0 = mesh->point[pt->v[MMG5_idir[ifac][0]]].c;
    barycoord[ifac].val = -( (coord[0]-c0[0])*normal[0] +
                             (coord[1]-c0[1])*normal[1] +
                             (coord[2]-c0[2])*normal[2] )/vol;
    barycoord[ifac].idx = ifac;
  }

  return 1;
}

/**
 * \param mesh pointer to the mesh structure
 * \param ptr pointer to the current triangle
 * \param k index of the triangle
 * \param coord pointer to the point coordinates
 * \param normal unit normal of the current triangle
 * \param barycoord pointer to the point barycentric coordinates in the current
 * tetra
 *
 * \return 0 if the point is outside the element, 1 if inside.
 *
 *  Compute the barycentric coordinates of a given point in a given triangle,
 *  sort them and evaluate if the point is inside.
 *
 */
int PMMG_barycoord2d_evaluate( MMG5_pMesh mesh,MMG5_pTria ptr,int k,
                               double *coord,double *triaNormal,
                               PMMG_barycoord *barycoord ) {

  /* Get barycentric coordinates and sort them in ascending order */
  PMMG_barycoord2d_compute(mesh, ptr, k, coord, triaNormal, barycoord);
  qsort(barycoord,3,sizeof(PMMG_barycoord),PMMG_barycoord_compare);

  /* Return inside/outside status */
  return PMMG_barycoord_isInside( barycoord );
}

/**
 * \param mesh pointer to the mesh structure
 * \param pt pointer to the current tetra
 * \param coord pointer to the point coordinates
 * \param faceAreas oriented face areas of the current tetrahedron
 * \param barycoord pointer to the point barycentric coordinates in the current
 * tetra
 *
 * \return 0 if the point is outside the element, 1 if inside.
 *
 *  Compute the barycentric coordinates of a given point in a given tetrahedron,
 *  sort them and evaluate if the point is inside.
 *
 */
int PMMG_barycoord3d_evaluate( MMG5_pMesh mesh,MMG5_pTetra pt,
                               double *coord,double *faceAreas,
                               PMMG_barycoord *barycoord ) {

  /* Get barycentric coordinates and sort them in ascending order */
  PMMG_barycoord3d_compute(mesh, pt, coord, faceAreas, barycoord);
  qsort(barycoord,4,sizeof(PMMG_barycoord),PMMG_barycoord_compare);

  /* Return inside/outside status */
  return PMMG_barycoord_isInside( barycoord );
}

/**
 * \param mesh pointer to the background mesh structure
 * \param k index of the triangle to analyze
 * \param ppt pointer to the point to locate
 * \param barycoord barycentric coordinates of the point to be located
 *
 * \return 1 if found; 0 if not found
 *
 *  Locate a surface point in its closest background triangles, and find its
 *  closest point.
 *
 */
int PMMG_barycoord2d_getClosest( MMG5_pMesh mesh,int k,MMG5_pPoint ppt,
                                 PMMG_barycoord *barycoord ) {
  MMG5_pTria ptr;
  double *c,dist[3],norm,min;
  int i,d,itarget;

  ptr = &mesh->tria[k];

  c = mesh->point[ptr->v[0]].c;
  for( d = 0; d < 3; d++ )
    dist[d] = ppt->c[d] - c[d];
  norm = sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
  min = norm;
  itarget = 0;

  for( i = 1; i < 3; i++ ) {
    c = mesh->point[ptr->v[i]].c;
    for( d = 0; d < 3; d++ )
      dist[d] = ppt->c[d] - c[d];
    norm = sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]);
    if( norm < min ) {
      min = norm;
      itarget = i;
    }
  }

  for( i = 0; i < 3; i++ ) {
    barycoord[i].val = 0.0;
    barycoord[i].idx = i;
  }
  barycoord[itarget].val = 1.0;

  return 1;
}

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
 * \author Luca Cirrottola (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */
#include "parmmg.h"
#include "locate_pmmg.h"

/**
 * \param minNew lower bounds of the new box in each space direction
 * \param maxNew upper bounds of the new box in each space direction
 * \param minOld lower bounds of the old box in each space direction
 * \param maxOld upper bounds of the old box in each space direction
 *
 * \return 0 if not intersectiong, 1 otherwise
 *
 *  Check that the intersection of two bounding boxes is not empty.
 *
 */
int PMMG_intersect_boundingBox( double *minNew, double *maxNew,
                                double *minOld, double *maxOld ) {
  int idim;

  for( idim = 0; idim < 3; idim++ ) {
    if( maxNew[idim] <= minOld[idim] )
      return 0;
    if( maxOld[idim] <= minNew[idim] )
      return 0;
  }

  return 1;
}

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
void PMMG_get_baryCoord( double *val,PMMG_baryCoord *phi,int ndim ) {
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
int PMMG_compare_baryCoord( const void *a,const void *b ) {
  PMMG_baryCoord *coord_a;
  PMMG_baryCoord *coord_b;

  coord_a = (PMMG_baryCoord *)a;
  coord_b = (PMMG_baryCoord *)b;

  if( coord_a->val > coord_b-> val ) return 1;
  if( coord_a->val < coord_b-> val ) return -1;

  return 0;
}

int PMMG_allPositive_baryCoord( PMMG_baryCoord *phi ) {
  if( phi[0].val > -MMG5_EPS )
    return 1;
  else
    return 0;
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
int PMMG_compute_baryCoord2d( MMG5_pMesh mesh,MMG5_pTria ptr,int k,double *coord,
                              double *normal,PMMG_baryCoord *barycoord ) {
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
int PMMG_compute_baryCoord3d( MMG5_pMesh mesh,MMG5_pTetra pt,double *coord,
                              double *faceAreas,PMMG_baryCoord *barycoord ) {
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
int PMMG_locatePointInClosestTria( MMG5_pMesh mesh,int k,MMG5_pPoint ppt,
                                   PMMG_baryCoord *barycoord ) {
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

/**
 * \param mesh pointer to the background mesh structure
 * \param ptr pointer to the triangle to analyze
 * \param k index of the triangle
 * \param ppt pointer to the point to locate
 * \param triaNormal unit normal of the current triangle
 * \param barycoord barycentric coordinates of the point to be located
 *
 * \return 1 if found; 0 if not found
 *
 *  Locate a point in a background triangles, and provide its barycentric
 *  coordinates.
 *
 */
int PMMG_locatePointInTria( MMG5_pMesh mesh,MMG5_pTria ptr,int k,MMG5_pPoint ppt,
                            double *triaNormal,PMMG_baryCoord *barycoord,
                            double *closestDist,int *closestTria ) {
  MMG5_pPoint    ppt0,ppt1;
  double         h,hmax;
  int            j,d,found = 0;

  /** Mark tria */
  ptr->flag = mesh->base;

  /** Get barycentric coordinates and sort them in ascending order */
  PMMG_compute_baryCoord2d(mesh, ptr, k, ppt->c, triaNormal, barycoord);
  qsort(barycoord,3,sizeof(PMMG_baryCoord),PMMG_compare_baryCoord);

  /* Rough check on the distance from the surface */
  hmax = 0.0;
  for( j = 0; j < 3; j++ ) {
    ppt0 = &mesh->point[ptr->v[j]];
    ppt1 = &mesh->point[ptr->v[MMG5_inxt2[j]]];
    h = 0.0;
    for( d = 0; d < 3; d++ )
      h += (ppt0->c[d]-ppt1->c[d])*(ppt0->c[d]-ppt1->c[d]);
    h = sqrt(h);
    if( h > hmax ) hmax = h;
  }
  if( fabs(barycoord[3].val) > hmax ) return 0;

  /** Save element index (with negative sign) if it is the closest one */
  if( fabs(barycoord[3].val) < *closestDist ) {
    *closestDist = fabs(barycoord[3].val);
    *closestTria = -k;
  }

  /** Exit if inside the element */
  if( PMMG_allPositive_baryCoord( barycoord ) ) found = 1;

  return found;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param pt pointer to the tetra to analyze
 * \param ppt pointer to the point to locate
 * \param faceAreas oriented face areas of the current tetrahedra
 * \param barycoord barycentric coordinates of the point to be located
 *
 * \return 1 if found; 0 if not found
 *
 *  Locate a point in a background tetrahedron, and provide its barycentric
 *  coordinates.
 *
 */
int PMMG_locatePointInTetra( MMG5_pMesh mesh,MMG5_pTetra pt,MMG5_pPoint ppt,
                             double *faceAreas,PMMG_baryCoord *barycoord ) {
  int found = 0;

  /** Mark tetra */
  pt->flag = mesh->base;

  /** Get barycentric coordinates and sort them in ascending order */
  PMMG_compute_baryCoord3d(mesh, pt, ppt->c, faceAreas, barycoord);
  qsort(barycoord,4,sizeof(PMMG_baryCoord),PMMG_compare_baryCoord);

  /** Exit if inside the element */
  if( PMMG_allPositive_baryCoord( barycoord ) ) found = 1;

  return found;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ppt pointer to the point to locate
 * \param init index of the starting element
 * \param triaNormals unit normals of the all triangles in the mesh
 * \param barycoord barycentric coordinates of the point to be located
 * \param ip current mesh point
 * \param igrp current mesh group
 *
 * \return ie if positive, index of the target element; if negative, index of
 * the closest element; 0 if not found
 *
 *  Locate a point in a background mesh surface by traveling the triangles
 *  adjacency.
 *
 */
int PMMG_locatePointBdy( MMG5_pMesh mesh,MMG5_pPoint ppt,int init,
                         double *triaNormals,PMMG_baryCoord *barycoord,
                         int ip,int igrp ) {
  MMG5_pTria     ptr,ptr1;
  int            *adjt,iel,i,idxTria,step,closestTria;
  double         vol,eps,closestDist;
  static int     mmgWarn0=0,mmgWarn1=0;

  if(!init)
    idxTria = 1;
  else
    idxTria = init;

  assert( idxTria <= mesh->nt );

  step = 0;
  ++mesh->base;

  closestTria = 0;
  closestDist = 1.0e10;

  while(step <= mesh->nt) {
    step++;

    /** Get tria */
    ptr = &mesh->tria[idxTria];
    if ( !MG_EOK(ptr) ) continue;

    /** Exit the loop if you find the element */
    if( PMMG_locatePointInTria( mesh, ptr, idxTria, ppt, &triaNormals[3*idxTria],
                                barycoord, &closestDist, &closestTria ) ) break;

    /** Compute new direction */
    adjt = &mesh->adjt[3*(idxTria-1)+1];
    for( i=0; i<3; i++ ) {
      iel = adjt[barycoord[i].idx]/3;

      /* Skip if on boundary */
      if (!iel) continue;

      /* Skip if already marked */
      ptr1 = &mesh->tria[iel];
      if(ptr1->flag == mesh->base) continue;

      /* Get next otherwise */
      idxTria = iel;
      break;
    }

    /** Stuck: Start exhaustive research */
    if (i == 3) step = mesh->nt+1;

  }

  /** Boundary hit or cyclic path: Perform exhaustive research */
  if( step == (mesh->nt+1) ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      if ( mesh->info.imprim > PMMG_VERB_DETQUAL ) {
        fprintf(stderr,"\n  ## Warning %s: Cannot locate point,"
                " performing exhaustive research.\n",__func__);
      }
    }

    for( idxTria=1; idxTria<mesh->nt+1; idxTria++ ) {

      /** Get tetra */
      ptr = &mesh->tria[idxTria];
      if ( !MG_EOK(ptr) ) continue;

      /*¨Skip already analized tetras */
      if( ptr->flag == mesh->base ) continue;

      /** Exit the loop if you find the element */
      if( PMMG_locatePointInTria( mesh, ptr, idxTria, ppt,
                                  &triaNormals[3*idxTria], barycoord,
                                  &closestDist, &closestTria ) ) break;

    }

    /** Element not found: Return the closest one with negative sign (if found) */
    if ( idxTria == mesh->nt+1 ) {
      if ( !mmgWarn1 ) {
        mmgWarn1 = 1;
        if ( mesh->info.imprim > PMMG_VERB_VERSION ) {
          fprintf(stderr,"\n  ## Warning %s: Point not located, smallest external area %e.",
                  __func__,closestDist);
        }
      }
      return closestTria;
    }

  }
  return idxTria;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ppt pointer to the point to locate
 * \param init index of the starting element
 * \param faceAreas oriented face areas of the all tetrahedra in the mesh
 * \param barycoord barycentric coordinates of the point to be located
 * \param ip current mesh point
 * \param igrp current mesh group
 *
 * \return ie if positive, index of the target element; if negative, index of
 * the closest element; 0 if not found
 *
 *  Locate a point in a background mesh by traveling the elements adjacency.
 *
 */
int PMMG_locatePointVol( MMG5_pMesh mesh,MMG5_pPoint ppt,int init,
                         double *faceAreas,PMMG_baryCoord *barycoord,
                         int ip,int igrp ) {
  MMG5_pTetra    pt,pt1;
  int            *adja,iel,i,idxTet,step,closestTet;
  double         vol,eps,closestDist;
  static int     mmgWarn0=0,mmgWarn1=0;

  if(!init)
    idxTet = 1;
  else
    idxTet = init;

  assert( idxTet <= mesh->ne );

  step = 0;
  ++mesh->base;
  while(step <= mesh->ne) {
    step++;

    /** Get tetra */
    pt = &mesh->tetra[idxTet];
    if ( !MG_EOK(pt) ) continue;

    /** Exit the loop if you find the element */
    if( PMMG_locatePointInTetra( mesh, pt, ppt,&faceAreas[12*idxTet],
                                 barycoord ) ) break;

    /** Compute new direction */
    adja = &mesh->adja[4*(idxTet-1)+1];
    for( i=0; i<4; i++ ) {
      iel = adja[barycoord[i].idx]/4;

      /* Skip if on boundary */
      if (!iel) continue;

      /* Skip if already marked */
      pt1 = &mesh->tetra[iel];
      if(pt1->flag == mesh->base) continue;

      /* Get next otherwise */
      idxTet = iel;
      break;
    }

    /** Stuck: Start exhaustive research */
    if (i == 4) step = mesh->ne+1;

  }

  /** Boundary hit or cyclic path: Perform exhaustive research */
  if( step == (mesh->ne+1) ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      if ( mesh->info.imprim > PMMG_VERB_DETQUAL ) {
        fprintf(stderr,"\n  ## Warning %s: Cannot locate point,"
                " performing exhaustive research.\n",__func__);
      }
    }

    closestTet = 0;
    closestDist = 1.0e10;
    for( idxTet=1; idxTet<mesh->ne+1; idxTet++ ) {

      /** Get tetra */
      pt = &mesh->tetra[idxTet];
      if ( !MG_EOK(pt) ) continue;

      /*¨Skip already analized tetras */
      if( pt->flag == mesh->base ) continue;

      /** Exit the loop if you find the element */
      if( PMMG_locatePointInTetra( mesh, pt, ppt,&faceAreas[12*idxTet],
                                   barycoord ) ) break;

      /** Save element index (with negative sign) if it is the closest one */
      vol = pt->qual;
      if( fabs(barycoord[0].val)*vol < closestDist ) {
        closestDist = fabs(barycoord[0].val)*vol;
        closestTet = -idxTet;
      }

    }

    /** Element not found: Return the closest one with negative sign (if found) */
    if ( idxTet == mesh->ne+1 ) {
      if ( !mmgWarn1 ) {
        mmgWarn1 = 1;
        if ( mesh->info.imprim > PMMG_VERB_VERSION ) {
          fprintf(stderr,"\n  ## Warning %s: Point not located, smallest external volume %e.",
                  __func__,closestDist);
        }
      }
      return closestTet;
    }

  }
  return idxTet;
}

/**
 * \param mesh pointer to the current mesh structure
 * \param meshOld pointer to the background mesh structure
 *
 *  For each point in the background mesh, store the index of a neighbouring
 *  triangle or tetrahedron.
 *
 */
void PMMG_locate_setStart( MMG5_pMesh mesh,MMG5_pMesh meshOld ) {
  MMG5_pPoint ppt;
  MMG5_pTria  ptr;
  MMG5_pTetra pt;
  int         ie,iloc,ip;

  /* Reset indices */
  for( ip = 1; ip <= meshOld->np; ip++ )
    meshOld->point[ip].s = 0;
  for( ip = 1; ip <= mesh->np; ip++ )
    mesh->point[ip].s = 0;


//  /* Store triangle index */
//  for( ie = 1; ie <= meshOld->nt; ie++ ) {
//    ptr = &meshOld->tria[ie];
//    for( iloc = 0; iloc < 3; iloc++ ) {
//      ip = ptr->v[iloc];
//      ppt = &meshOld->point[ip];
//      assert( ppt->tag & MG_BDY );
//      if( ppt->s ) continue;
//      ppt->s = -ie;
//    }
//  }
//
//  /* Fetch triangle index */
//  for( ip = 1; ip <= mesh->np; ip++ ) {
//    ppt = &mesh->point[ip];
//    if( !(ppt->tag & MG_BDY) ) continue;
//    ppt->s = -meshOld->point[ppt->src].s;
//  }


  /* Store tetra index */
  for( ie = 1; ie <= meshOld->ne; ie++ ) {
    pt = &meshOld->tetra[ie];
    for( iloc = 0; iloc < 4; iloc++ ) {
      ip = pt->v[iloc];
      ppt = &meshOld->point[ip];
      if( ppt->s > 0 ) continue;
      ppt->s = ie;
    }
  }

  /* Fetch tetra index */
  for( ip = 1; ip <= mesh->np; ip++ ) {
    ppt = &mesh->point[ip];
    if( ppt->tag & MG_BDY ) continue;
    ppt->s = meshOld->point[ppt->src].s;
  }


}

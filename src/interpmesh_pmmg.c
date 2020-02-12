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
 * \file interpmesh_pmmg.c
 * \brief Interpolate data from a background mesh to the current mesh.
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \author Nikos Pattakos (Inria)
 * \author Luca Cirrottola (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */
#include "parmmg.h"
#include "locate_pmmg.h"

/**
 * \param coord pointer to a double array
 * \param barycoord pointer to the point barycentric coordinates in the current
 * tetra
 *
 *  Get barycentric coordinates.
 *
 */
void PMMG_get_baryCoord( double *val, PMMG_baryCoord *phi ) {
  int i;

  for( i = 0; i < 4; i++ )
    val[phi[i].idx] = phi[i].val;

}

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
 * \param mesh pointer to the mesh structure
 * \param ptr pointer to the current triangle
 * \param k index of the triangle
 * \param coord pointer to the point coordinates
 * \param faceAreas oriented face areas of the current tetrahedron
 * \param barycoord pointer to the point barycentric coordinates in the current
 * tetra
 *
 * \return 0 if fail, 1 if success
 *
 *  Compute the barycentric coordinates of a given point in a given triangle.
 *
 */
int PMMG_compute_baryCoord2d( MMG5_pMesh mesh, MMG5_pTria ptr, int k,
                              double *coord, double *triaNormals,
                              PMMG_baryCoord *barycoord ) {
  double *c1,*c2,*normal,vol;
  int    ia;

  /* Retrieve tria area */
  vol = ptr->qual;

  /* Retrieve tria unit normal */
  normal = &triaNormals[3*k];
 
  /* Retrieve face areas and compute barycentric coordinates */
  for( ia = 0; ia < 3; ia++ ) {
    c1 = mesh->point[ptr->v[MMG5_inxt3[ia]]].c;
    c2 = mesh->point[ptr->v[MMG5_inxt3[ia+1]]].c;

    barycoord[ia].val = MMG2D_quickarea( coord, c1, c2 )/vol;
    barycoord[ia].idx = ia;
  }

  /* Store normal distance in the third coordinate */
  barycoord[3].val = (coord[0]-c1[0])*normal[0] +
                     (coord[1]-c1[1])*normal[1] +
                     (coord[2]-c1[2])*normal[2];
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
int PMMG_compute_baryCoord3d( MMG5_pMesh mesh, MMG5_pTetra pt,
                    double *coord, double *faceAreas, PMMG_baryCoord *barycoord ) {
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

/**
 * \param mesh pointer to the background mesh structure
 * \param ptr pointer to the triangle to analyze
 * \param k index of the triangle
 * \param ppt pointer to the point to locate
 * \param triaNormalsAreas unit normals of the current triangle
 * \param barycoord barycentric coordinates of the point to be located
 *
 * \return 1 if found; 0 if not found
 *
 *  Locate a point in a background triangles, and provide its barycentric
 *  coordinates.
 *
 */
int PMMG_locatePointInTria( MMG5_pMesh mesh, MMG5_pTria ptr, int k,
                            MMG5_pPoint ppt, double *triaNormals,
                            PMMG_baryCoord *barycoord ) {
  double         eps;
  int            found = 0;

  eps = MMG5_EPS;

  /** Mark tria */
  ptr->flag = mesh->base;

  /** Get barycentric coordinates and sort them in ascending order */
  PMMG_compute_baryCoord2d(mesh, ptr, k, ppt->c, triaNormals, barycoord);
  qsort(barycoord,3,sizeof(PMMG_baryCoord),PMMG_compare_baryCoord);

  /** Exit if inside the element */
  if( barycoord[0].val > -eps ) found = 1;

  return found;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ptr pointer to the tetra to analyze
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
int PMMG_locatePointInTetra( MMG5_pMesh mesh, MMG5_pTetra ptr, MMG5_pPoint ppt,
                             double *faceAreas, PMMG_baryCoord *barycoord ) {
  double         eps;
  int            found = 0;

  eps = MMG5_EPS;

  /** Mark tetra */
  ptr->flag = mesh->base;

  /** Get barycentric coordinates and sort them in ascending order */
  PMMG_compute_baryCoord3d(mesh, ptr, ppt->c, faceAreas, barycoord);
  qsort(barycoord,4,sizeof(PMMG_baryCoord),PMMG_compare_baryCoord);

  /** Exit if inside the element */
  if( barycoord[0].val > -eps ) found = 1;

  return found;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ppt pointer to the point to locate
 * \param init index of the starting element
 * \param triaNormals unit normals of the all triangles in the mesh
 * \param barycoord barycentric coordinates of the point to be located
 *
 * \return ie if positive, index of the target element; if negative, index of
 * the closest element; 0 if not found
 *
 *  Locate a point in a background mesh surface by traveling the triangles
 *  adjacency.
 *
 */
int PMMG_locatePointBdy( MMG5_pMesh mesh, MMG5_pPoint ppt, int init,
                         double *triaNormals, PMMG_baryCoord *barycoord ) {
  MMG5_pTria     ptr,ptr1;
  int            *adjt,iel,i,idxTria,step,closestTria;
  double         vol,eps,closestDist;
  static int     mmgWarn0=0,mmgWarn1=0;

  if(!init)
    idxTria = 1;
  else
    idxTria = init;

  step = 0;
  ++mesh->base;
  while(step <= mesh->nt) {
    step++;

    /** Get tria */
    ptr = &mesh->tria[idxTria];
    if ( !MG_EOK(ptr) ) continue;

    /** Exit the loop if you find the element */
    if( PMMG_locatePointInTria( mesh, ptr, idxTria, ppt,
                                &triaNormals[3*idxTria], barycoord ) ) break;

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

    closestTria = 0;
    closestDist = 1.0e10;
    for( idxTria=1; idxTria<mesh->nt+1; idxTria++ ) {

      /** Get tetra */
      ptr = &mesh->tria[idxTria];
      if ( !MG_EOK(ptr) ) continue;

      /*¨Skip already analized tetras */
      if( ptr->flag == mesh->base ) continue;

      /** Exit the loop if you find the element */
      if( PMMG_locatePointInTria( mesh, ptr, idxTria, ppt,
                                  &triaNormals[3*idxTria], barycoord ) ) break;

      /** Save element index (with negative sign) if it is the closest one */
      vol = ptr->qual;
      if( fabs(barycoord[0].val)*vol < closestDist ) {
        closestDist = fabs(barycoord[0].val)*vol;
        closestTria = -idxTria;
      }

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
 *
 * \return ie if positive, index of the target element; if negative, index of
 * the closest element; 0 if not found
 *
 *  Locate a point in a background mesh by traveling the elements adjacency.
 *
 */
int PMMG_locatePoint( MMG5_pMesh mesh, MMG5_pPoint ppt, int init,
                      double *faceAreas, PMMG_baryCoord *barycoord ) {
  MMG5_pTetra    ptr,pt1;
  int            *adja,iel,i,idxTet,step,closestTet;
  double         vol,eps,closestDist;
  static int     mmgWarn0=0,mmgWarn1=0;

  if(!init)
    idxTet = 1;
  else
    idxTet = init;

  step = 0;
  ++mesh->base;
  while(step <= mesh->ne) {
    step++;

    /** Get tetra */
    ptr = &mesh->tetra[idxTet];
    if ( !MG_EOK(ptr) ) continue;

    /** Exit the loop if you find the element */
    if( PMMG_locatePointInTetra( mesh, ptr, ppt,&faceAreas[12*idxTet],
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
      ptr = &mesh->tetra[idxTet];
      if ( !MG_EOK(ptr) ) continue;

      /*¨Skip already analized tetras */
      if( ptr->flag == mesh->base ) continue;

      /** Exit the loop if you find the element */
      if( PMMG_locatePointInTetra( mesh, ptr, ppt,&faceAreas[12*idxTet],
                                   barycoord ) ) break;

      /** Save element index (with negative sign) if it is the closest one */
      vol = ptr->qual;
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
 * \param m pointer to the 2x2 symmetric matrix
 * \param im pointer to the inverse matrix
 *
 * \return 0 if fail, 1 if success
 *
 *  Invert 2x2 symmetric matrix.
 */
int PMMG_invmat22( double *m, double *im ) {
  double det;

  det = m[0]*m[2] - m[1]*m[1];
  if ( fabs(det) < MMG5_EPS*MMG5_EPS ) {
    fprintf(stderr,"\n  ## Error: %s: null metric det : %E \n",
            __func__,det);
    return 0;
  }
  det = 1.0 / det;

  im[0] =  det*m[3];
  im[1] = -det*m[1];
  im[2] =  det*m[0];

  return 1;
}

/**
 * \param mesh pointer to the current mesh
 * \param met pointer to the current metrics
 * \param oldMet pointer to the background metrics
 * \param ptr pointer to the target background triangle
 * \param ip index of the current point
 * \param phi barycentric coordinates of the point to be interpolated
 *
 * \return 0 if fail, 1 if success
 *
 *  Linearly interpolate point metrics on a target background triangle.
 *  This function is analogous to the MMG5_interp4bar_iso() function.
 */
int PMMG_interp3bar_iso( MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,
                         MMG5_pTria ptr,int ip,double *phi ) {
  int iloc,i,ier;

  assert( met->size == 1 );

  /** Linear interpolation of the squared size */
  met->m[ip] = phi[0]*oldMet->m[ptr->v[0]] +
               phi[1]*oldMet->m[ptr->v[1]] + 
               phi[2]*oldMet->m[ptr->v[2]];

  return 1;
}

/**
 * \param mesh pointer to the current mesh
 * \param met pointer to the current metrics
 * \param oldMet pointer to the background metrics
 * \param ptr pointer to the target background triangle
 * \param ip index of the current point
 * \param phi barycentric coordinates of the point to be interpolated
 *
 * \return 0 if fail, 1 if success
 *
 *  Linearly interpolate point the metrics inverse on a target background
 *  triangle.
 *  This function is analogous to the MMG5_interp4barintern() function.
 *
 */
int PMMG_interp3bar_ani( MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,
                         MMG5_pTria ptr,int ip,double *phi ) {
  double dm[3][3],mi[3][3],m[3];
  int    iloc,i,isize,nsize,ier;

  assert( met->size == 3 );
  nsize  = met->size;

  for( i=0; i<3; i++ ) {
    for(isize = 0; isize < nsize; isize++ )
      dm[i][isize] = oldMet->m[nsize*ptr->v[i]+isize];
    if( !PMMG_invmat22(dm[i],mi[i]) ) return 0;
  }

  /** Linear interpolation of the metrics */
  for( isize = 0; isize < nsize; isize++ ) {
    m[isize] = phi[0]*mi[0][isize]+
               phi[1]*mi[1][isize]+
               phi[2]*mi[2][isize];
  }

  if( !PMMG_invmat22(m,mi[0]) ) return 0;
  for( isize = 0; isize < nsize; isize++ )
    met->m[nsize*ip+isize] = mi[0][isize];

  return 1;
}

/**
 * \param mesh pointer to the current mesh
 * \param met pointer to the current metrics
 * \param oldMet pointer to the background metrics
 * \param pt pointer to the target background tetrahedron
 * \param ip index of the current point
 * \param phi barycentric coordinates of the point to be interpolated
 *
 * \return 0 if fail, 1 if success
 *
 *  Linearly interpolate point metrics on a target background tetrahedron.
 *  This function is analogous to the MMG5_interp4bar_iso() function.
 */
int PMMG_interp4bar_iso( MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,
                                  MMG5_pTetra pt,int ip,double *phi ) {
  int iloc,i,ier;

  assert( met->size == 1 );

  /** Linear interpolation of the squared size */
  met->m[ip] += phi[0]*oldMet->m[pt->v[0]]+
                phi[1]*oldMet->m[pt->v[1]]+
                phi[2]*oldMet->m[pt->v[2]]+
                phi[3]*oldMet->m[pt->v[3]];

  return 1;
}

/**
 * \param mesh pointer to the current mesh
 * \param met pointer to the current metrics
 * \param oldMet pointer to the background metrics
 * \param pt pointer to the target background tetrahedron
 * \param ip index of the current point
 * \param phi barycentric coordinates of the point to be interpolated
 *
 * \return 0 if fail, 1 if success
 *
 *  Linearly interpolate point the metrics inverse on a target background
 *  tetrahedron.
 *  This function is analogous to the MMG5_interp4barintern() function.
 *
 */
int PMMG_interp4bar_ani( MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,
                                  MMG5_pTetra pt,int ip,double *phi ) {
  double dm[4][6],mi[4][6],m[6];
  int    iloc,i,isize,nsize,ier;

  assert( met->size == 6 );
  nsize  = met->size;

  for( i=0; i<4; i++ ) {
    for(isize = 0; isize < nsize; isize++ )
      dm[i][isize] = oldMet->m[nsize*pt->v[i]+isize];
    if( !MMG5_invmat(dm[i],mi[i]) ) return 0;
  }

  /** Linear interpolation of the metrics */
  for( isize = 0; isize < nsize; isize++ ) {
    m[isize] = phi[0]*mi[0][isize] + phi[1]*mi[1][isize] +
               phi[2]*mi[2][isize] + phi[3]*mi[3][isize];
  }

  if( !MMG5_invmat(m,mi[0]) ) return 0;
  for( isize = 0; isize < nsize; isize++ )
    met->m[nsize*ip+isize] = mi[0][isize];

  return 1;
}

/**
 * \param mesh pointer to the current mesh.
 * \param oldMesh pointer to the background mesh.
 * \param met pointer to the current metrics.
 * \param oldMet pointer to the background metrics.
 * \param permNodGlob permutation array for nodes.
 *
 * \return 0 if fail, 1 if success
 *
 * Copy the metric of a freezed interface point.
 *
 */
int PMMG_copyMetrics_point( MMG5_pMesh mesh,MMG5_pMesh oldMesh,
                            MMG5_pSol met,MMG5_pSol oldMet,int* permNodGlob) {
  MMG5_pPoint    ppt;
  int            isize,nsize,ip;

  if ( !mesh->info.inputMet || mesh->info.hsiz > 0.0 ) return 1;

  nsize   = met->size;

#warning Luca: when surface adapt will be ready, distinguish BDY from PARBDY

  /** Freezed points: Copy the  metrics  */
  if ( (!oldMesh->info.renum) || !permNodGlob ) {
    /* No permutation array: simple copy */
    for( ip = 1; ip <= oldMesh->np; ++ip ) {
      ppt = &oldMesh->point[ip];
      if( !MG_VOK(ppt) ) continue;

      ppt->flag = mesh->base;
      if( ppt->tag & MG_REQ ) {
        for( isize = 0; isize<nsize; isize++ ) {
          met->m[nsize*ip+isize] = oldMet->m[nsize*ip+isize];
        }
      }
    }
  }
  else {
    /* Due to the scotch renumbering we must copy from the old mesh to the
     * new one */
    for( ip = 1; ip <= oldMesh->np; ++ip ) {
      ppt = &oldMesh->point[ip];
      if( !MG_VOK(ppt) ) continue;

      ppt->flag = mesh->base;
      if( ppt->tag & MG_REQ ) {
        for( isize = 0; isize<nsize; isize++ ) {
          met->m[nsize*permNodGlob[ip]+isize] = oldMet->m[nsize*ip+isize];
        }
      }
    }
  }

  return 1;
}


/**
 * \param mesh pointer to the current mesh structure.
 * \param oldMesh pointer to the background mesh structure.
 * \param met pointer to the current metrics structure.
 * \param oldMet pointer to the background metrics structure.
 * \param faceAreas pointer to the array of oriented face areas.
 * \param permNodGlob permutation array of nodes.
 *
 * \return 0 if fail, 1 if success
 *
 *  Interpolate metrics for all groups from background to current meshes.
 *  Do nothing if no metrics is provided (info.inputMet == 0), otherwise:
 *  - if the metrics is constant, recompute it;
 *  - else, interpolate the non-constant metrics.
 *
 *  Oriented face areas are pre-computed in this function before proceeding
 *  with the localization.
 *
 */
int PMMG_interpMetrics_mesh( MMG5_pMesh mesh,MMG5_pMesh oldMesh,
                             MMG5_pSol met,MMG5_pSol oldMet,
                             double *faceAreas,int *permNodGlob ) {
  MMG5_pTetra pt;
  MMG5_pPoint ppt;
  PMMG_baryCoord barycoord[4];
  double      coord[4],*normal;
  int         ip,istart,ie,ifac,ia,ib,ic,iloc;
  int         ier;
  static int  mmgWarn=0;

  if( mesh->info.inputMet != 1 ) {

    /* Nothing to do */
    return 1;

  } else {

    if( mesh->info.hsiz > 0.0 ) {

      /* Compute constant metrics */
      if ( !MMG3D_Set_constantSize(mesh,met) ) return 0;

    } else {

      /** Pre-compute oriented face areas */
      for( ie = 1; ie <= oldMesh->ne; ie++ ) {
        pt = &oldMesh->tetra[ie];
        /* Store tetra volume in the qual field */
        pt->qual = MMG5_orvol( oldMesh->point, pt->v );
        /* Store oriented face normals */
        for( ifac = 0; ifac < 4; ifac++ ) {
          normal = &faceAreas[12*ie+3*ifac];
          ia = pt->v[MMG5_idir[ifac][0]];
          ib = pt->v[MMG5_idir[ifac][1]];
          ic = pt->v[MMG5_idir[ifac][2]];
          ier = MMG5_nonUnitNorPts( oldMesh,ia,ib,ic,normal );
        }
      }


      /** Interpolate metrics */
      oldMesh->base = 0;
      for ( ie = 1; ie < oldMesh->ne+1; ie++ ) {
        pt = &oldMesh->tetra[ie];
        if ( !MG_EOK(pt) ) continue;
        pt->flag = oldMesh->base;
      }

      mesh->base++;
      istart = 1;
      for( ie = 1; ie <= mesh->ne; ie++ ) {
        pt = &mesh->tetra[ie];
        if( !MG_EOK(pt) ) continue;
        for( iloc = 0; iloc < 4; iloc++ ) {
          ip = pt->v[iloc];
          ppt = &mesh->point[ip];
          if( !MG_VOK(ppt) ) continue;

          /* Skip already interpolated points */
          if( ppt->flag == mesh->base ) continue;

          if( ppt->tag & MG_REQ ) {
            /* Flag point as interpolated */
            ppt->flag = mesh->base;
            continue; // treated by copyMetric_points
          } else {

            /** Locate point in the old mesh */
            istart = PMMG_locatePoint( oldMesh, ppt, istart,
                                       faceAreas, barycoord );
            if( !istart ) {
              fprintf(stderr,"\n  ## Error: %s:"
                      " point %d not found, coords %e %e %e\n",__func__,
                      ip, mesh->point[ip].c[0],
                      mesh->point[ip].c[1],mesh->point[ip].c[2]);
              return 0;
            } else if( istart < 0 ) {
              if ( !mmgWarn ) {
                mmgWarn = 1;
                if ( mesh->info.imprim > PMMG_VERB_VERSION ) {
                  fprintf(stderr,"\n  ## Warning: %s: point %d not"
                          " found, coords %e %e %e\n",__func__,
                          ip, mesh->point[ip].c[0],mesh->point[ip].c[1],
                          mesh->point[ip].c[2]);
                }
              }
              istart = -istart;
            }

            /** Interpolate point metrics */
            PMMG_get_baryCoord( coord, barycoord );
            ier = PMMG_interp4bar(mesh,met,oldMet,
                                           &oldMesh->tetra[istart],
                                           ip,coord);

            /* Flag point as interpolated */
            ppt->flag = mesh->base;
          }
        }
      }
    }
  }

  return 1;
}

/**
 * \param parmesh pointer to the parmesh structure.
 * \param permNodGlob permutation array of nodes.
 *
 * \return 0 if fail, 1 if success
 *
 *  Interpolate metrics for all groups from background to current meshes.
 *  Do nothing if no metrics is provided (info.inputMet == 0), otherwise:
 *  - if the metrics is constant, recompute it;
 *  - else, interpolate the non-constant metrics.
 *
 */
int PMMG_interpMetrics( PMMG_pParMesh parmesh,int *permNodGlob ) {
  PMMG_pGrp   grp,oldGrp;
  MMG5_pMesh  mesh,oldMesh;
  MMG5_pSol   met,oldMet;
  MMG5_Hash   hash;
  size_t      memAv,oldMemMax;
  double      *faceAreas;
  int         igrp,ier;

  PMMG_TRANSFER_AVMEM_TO_PARMESH(parmesh,memAv,oldMemMax);

  /** Loop on current groups */
  ier = 1;
  for( igrp = 0; igrp < parmesh->ngrp; igrp++ ) {

    grp  = &parmesh->listgrp[igrp];
    mesh = grp->mesh;
    met  = grp->met;

    oldGrp  = &parmesh->old_listgrp[igrp];
    oldMesh = oldGrp->mesh;
    oldMet  = oldGrp->met;

    if( !oldMesh->adjt ) {
      PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh,memAv,oldMemMax);
      /* create surface adjacency */
      memset ( &hash, 0x0, sizeof(MMG5_Hash));
      if ( !MMG3D_hashTria(mesh,&hash) ) {
        MMG5_DEL_MEM(mesh,hash.item);
        fprintf(stderr,"\n  ## Hashing problem (2). Exit program.\n");
        return 0;
      }
      PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,mesh,memAv,oldMemMax);
    }

    /** Pre-allocate oriented face areas */
    if( ( mesh->info.inputMet == 1 ) && ( mesh->info.hsiz <= 0.0 ) ) {
      ier = 1;
      PMMG_MALLOC( parmesh,faceAreas,12*(oldMesh->ne+1),double,"faceAreas",ier=0 );
      if( !ier ) {
        PMMG_DEL_MEM(parmesh,faceAreas,double,"faceAreas");
        return 0;
      }
    }

    if( !PMMG_interpMetrics_mesh( mesh, oldMesh, met, oldMet,
                                  faceAreas, permNodGlob ) )
      ier = 0;

    /** Deallocate oriented face areas */
    if( ( mesh->info.inputMet == 1 ) && ( mesh->info.hsiz <= 0.0 ) ) {
      PMMG_DEL_MEM(parmesh,faceAreas,double,"faceAreas");
    }

    if( oldMesh->adjt ) {
      PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh,memAv,oldMemMax);
      MMG5_DEL_MEM(mesh,mesh->adjt);
      if( hash.item )
        MMG5_DEL_MEM(mesh,hash.item);
      PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,mesh,memAv,oldMemMax);
    }

  }

  return ier;
}

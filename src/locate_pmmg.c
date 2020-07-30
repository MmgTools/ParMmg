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
 * \param mesh pointer to the current mesh structure.
 * \param triaNormals pointer to the array of non-normalized triangle normals.
 *
 * \return 1.
 *
 *  Precompute non-normalized triangle normals.
 *
 */
int PMMG_precompute_triaNormals( MMG5_pMesh mesh,double *triaNormals ) {
  MMG5_pTria  ptr;
  double      *normal,dd;
  int         k,ia,ib,ic;
  int         ier;

  for( k = 1; k <= mesh->nt; k++ ) {
    ptr = &mesh->tria[k];
    ia = ptr->v[0];
    ib = ptr->v[1];
    ic = ptr->v[2];
    normal = &triaNormals[3*k];
    /* Store triangle unit normal and volume */
    ier = MMG5_nonUnitNorPts( mesh,ia,ib,ic,normal );
    ptr->qual = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
    dd = 1.0/ptr->qual;
    normal[0] *= dd;
    normal[1] *= dd;
    normal[2] *= dd;
  }

  return 1;
}

/**
 * \param mesh pointer to the current mesh structure.
 * \param faceAreas pointer to the array of oriented face areas.
 *
 * \return 1.
 *
 *  Precompute oriented face areas on tetrahedra.
 *
 */
int PMMG_precompute_faceAreas( MMG5_pMesh mesh,double *faceAreas ) {
  MMG5_pTetra pt;
  double      *normal;
  int         ie,ifac,ia,ib,ic;
  int         ier;

  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    /* Store tetra volume in the qual field */
    pt->qual = MMG5_orvol( mesh->point, pt->v );
    /* Store oriented face normals */
    for( ifac = 0; ifac < 4; ifac++ ) {
      normal = &faceAreas[12*ie+3*ifac];
      ia = pt->v[MMG5_idir[ifac][0]];
      ib = pt->v[MMG5_idir[ifac][1]];
      ic = pt->v[MMG5_idir[ifac][2]];
      ier = MMG5_nonUnitNorPts( mesh,ia,ib,ic,normal );
    }
  }

  return 1;
}

/**
 * \param parmesh pointer to the parmesh structure.
 * \param mesh pointer to the current mesh structure.
 * \param nodeTrias double pointer to the node triangles graph.
 *
 * \return 1.
 *
 *  Precompute node triangles graph on the surface.
 *
 */
int PMMG_precompute_nodeTrias( PMMG_pParMesh parmesh,MMG5_pMesh mesh,int **nodeTrias ) {
  MMG5_pTria  ptr;
  MMG5_pPoint ppt;
  int         ip,k,iloc;
  int         nt;

  for( ip = 1; ip <= mesh->np; ip++ ) {
    ppt = &mesh->point[ip];
    ppt->tmp = 0;
    ppt->flag = 0;
  }

  for( k = 1; k <= mesh->nt; k++ ) {
    ptr = &mesh->tria[k];
    for( iloc = 0; iloc < 3; iloc++ ) {
      mesh->point[ptr->v[iloc]].tmp++;
    }
  }

  /* Sum */
  for( ip = 2; ip <= mesh->np; ip++ )
    mesh->point[ip].tmp += mesh->point[ip-1].tmp;
  nt = mesh->point[mesh->np].tmp;

  /* Shift */
  for( ip = mesh->np; ip > 1; ip-- )
    mesh->point[ip].tmp = mesh->point[ip-1].tmp;
  mesh->point[1].tmp = 0;

  /* Allocate and fill */
  PMMG_MALLOC( parmesh,*nodeTrias,nt,int,"nodeTrias",return 0 );

  for( k = 1; k <= mesh->nt; k++ ) {
    ptr = &mesh->tria[k];
    for( iloc = 0; iloc < 3; iloc++ ) {
      ip = ptr->v[iloc];
      ppt = &mesh->point[ip];
      (*nodeTrias)[ppt->tmp+ppt->flag++] = k;
    }
  }

  return 1;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param iel index of the background triangle
 * \param iloc local index of the cone point in the background triangle
 * \param ppt pointer to the point to locate
 * \param list empty list of point neighbours
 *
 * \return 1 if found; 0 if not found
 *
 *  Locate a point in the shadow cone of a background point.
 *
 */
int PMMG_locatePointInCone( MMG5_pMesh mesh,int *nodeTrias,int iel,int iloc,
                            MMG5_pPoint ppt ) {
  MMG5_pTria     ptr;
  MMG5_pPoint    ppt0,ppt1;
  double         p[3],a[3],dist,alpha;
  int            *myTrias,nt,ip,jp,k,jloc,d;

  ppt0 = &mesh->point[mesh->tria[iel].v[iloc]];

  /* Mark point */
  ppt0->flag = mesh->base;

  /* Target point vector */
  for( d = 0; d < 3; d++ ) p[d] = ppt->c[d]-ppt0->c[d];
  dist = 0.0;
  for( d = 0; d < 3; d++ ) dist += p[d]*p[d];
  dist = sqrt(dist);

  /* Scan the neighbours */
  ip = mesh->tria[iel].v[iloc];
  nt = mesh->point[ip+1].tmp;
  myTrias = &nodeTrias[mesh->point[ip].tmp];
  for( k = 0; k < nt; k++ ){
    ptr = &mesh->tria[myTrias[k]];
    for( jloc = 0; jloc < 3; jloc++ ) {
      jp = ptr->v[jloc];
      ppt1 = &mesh->point[jp];
      if( ppt1->flag == ip ) continue;
      ppt1->flag = ip;
      /* Edge vector */
      for( d = 0; d < 3; d++ ) a[d] = ppt1->c[d]-ppt0->c[d];
      /* Rough check on maximum distance */
      if( dist > mesh->info.hausd ) {
        return 0;
      }
      /* Scalar product of the target vector with the edge vector */
      alpha = 0.0;
      for( d = 0; d < 3; d++ ) alpha += a[d]*p[d];
      if( alpha > 0.0 ) {
        return 0;
      }
    }
  }

  /* Found if all the neighbours have been scanned */
  return 1;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ptr pointer to the background triangle
 * \param k index of the background triangle
 * \param l local index of the edge
 * \param ppt pointer to the point to locate
 *
 * \return -1 if the point is too distant; 4 if found; the local index of the
 * node to check in triangle k otherwise.
 *
 *  Locate a point in the shadow wedge of a background edge. The adjacent
 *  triangles have already been tested.
 *
 */
int PMMG_locatePointInWedge( MMG5_pMesh mesh,MMG5_pTria ptr,int k,int l,MMG5_pPoint ppt,PMMG_barycoord *barycoord ) {
  MMG5_pTria  ptr1;
  MMG5_pPoint ppt0,ppt1;
  double      a[3],p[3],norm2,dist,alpha;
  int         i0,i1,d;

  /* Check nodal sides */
  i0 = MMG5_inxt2[l];
  i1 = MMG5_iprv2[l];
  ppt0 = &mesh->point[ptr->v[i0]];
  ppt1 = &mesh->point[ptr->v[i1]];

  /* Target point vector */
  for( d = 0; d < 3; d++ ) p[d] = ppt->c[d]-ppt0->c[d];

  /* Edge vector and norm */
  for( d = 0; d < 3; d++ ) a[d] = ppt1->c[d]-ppt0->c[d];
  norm2 = 0.0;
  for( d = 0; d < 3; d++ ) norm2 += a[d]*a[d];

  /* Scalar product of the target vector with the edge vector and orthogonal
   * projection */
  alpha = 0.0;
  for( d = 0; d < 3; d++ ) alpha += a[d]*p[d];
  for( d = 0; d < 3; d++ ) p[d] -= (alpha/norm2)*a[d];

  /* Rough check on maximum distance */
  dist = 0.0;
  for( d = 0; d < 3; d++ ) dist += p[d]*p[d];
  dist = sqrt(dist);
  if( dist > mesh->info.hausd ) return PMMG_UNSET;

  /* Check scalar product */
  if( alpha < 0.0 ) {
    ppt1->flag = mesh->base;
    return i0;
  } else if( alpha > norm2 ) {
    ppt0->flag = mesh->base;
    return i1;
  }

  /* Store inner product in the barycentric coordinate of the edge (forget
   * reordering) */
  for( d = 0; d < 3; d++ ) barycoord[d].idx = d;
  barycoord[l].val  = 0.0;
  barycoord[i0].val = 1.0-alpha/norm2;
  barycoord[i1].val = alpha/norm2;

  return 4;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ptr pointer to the triangle to analyze
 * \param ppt pointer to the point to locate
 * \param triaNormal unit normal of the current triangle
 *
 * \return 0 if too far, 1 otherwise.
 *
 *  Check point orthogonal distance from triangle.
 *
 */
int PMMG_locateChkDistTria( MMG5_pMesh mesh,MMG5_pTria ptr,MMG5_pPoint ppt,
                            double *triaNormal ) {
  MMG5_pPoint ppt0;
  double      norm,dist[3];
  int         d;

  /* Orthogonal distance */
  ppt0 = &mesh->point[ptr->v[0]];
  for( d = 0; d < 3; d++ )
    dist[d] = ppt->c[d]-ppt0->c[d];

  norm = 0.0;
  for( d = 0; d < 3; d++ )
    norm += dist[d]*triaNormal[d];
  norm = fabs(norm);

  if( norm > mesh->info.hausd ) return 0;

  return 1;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ptr pointer to the triangle to analyze
 * \param k index of the triangle
 * \param ppt pointer to the point to locate
 * \param triaNormal unit normal of the current triangle
 * \param barycoord barycentric coordinates of the point to be located
 * \param h distance from the triangle
 * \param closestDist distance from the closest triangle
 * \param closestTria index of the closest triangle (with negative sign)
 *
 * \return 1 if found; 0 if not found
 *
 *  Locate a point in a background triangles, and provide its barycentric
 *  coordinates.
 *
 */
int PMMG_locatePointInTria( MMG5_pMesh mesh,MMG5_pTria ptr,int k,MMG5_pPoint ppt,
                            double *triaNormal,PMMG_barycoord *barycoord,
                            double *h,double *closestDist,int *closestTria ) {
  MMG5_pPoint    ppt0,ppt1;
  double         hmax,a;
  double         norm,dist[3];
  int            j,d,found;

  /* Mark tria */
  ptr->flag = mesh->base;

  /* Evaluate point in tetra through barycentric coordinates */
  found = PMMG_barycoord2d_evaluate( mesh,ptr,k,ppt->c,triaNormal,barycoord );
  if( !found) return 0;

  /* Distance from center of mass */
  for( d = 0; d < 3; d++ )
    dist[d] = ppt->c[d];
  for( j = 0; j < 3; j++ ) {
    ppt0 = &mesh->point[ptr->v[j]];
    for( d = 0; d < 3; d++ )
      dist[d] -= ppt0->c[d]/3.0;
  }
  norm = 0;
  for( d = 0; d < 3; d++ )
    norm += dist[d]*dist[d];
  norm = sqrt(norm);

  /* Save element index if it is the closest one */
  if( norm < *closestDist ) {
    *closestDist = norm;
    *closestTria = k;
  }
  assert(*closestTria);

  /* Rough check on the distance from the surface */
  if( !PMMG_locateChkDistTria( mesh,ptr,ppt,triaNormal ) ) return 0;

  return found;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param pt pointer to the tetra to analyze
 * \param k index of the tetra
 * \param ppt pointer to the point to locate
 * \param faceAreas oriented face areas of the current tetrahedron
 * \param barycoord barycentric coordinates of the point to be located
 * \param closestDist pointer to the distance from the closest tetrahedron
 * \param closestTet pointer to the index of the closest tetrahedron
 *
 * \return 1 if found; 0 if not found
 *
 *  Locate a point in a background tetrahedron, and provide its barycentric
 *  coordinates.
 *
 */
int PMMG_locatePointInTetra( MMG5_pMesh mesh,MMG5_pTetra pt,int k,MMG5_pPoint ppt,
                             double *faceAreas,PMMG_barycoord *barycoord,
                             double *closestDist,int *closestTet) {
  double vol;
  int    found;

  /* Mark tetra */
  pt->flag = mesh->base;

  /* Evaluate point in tetra through barycentric coordinates */
  found = PMMG_barycoord3d_evaluate( mesh,pt,ppt->c,faceAreas,barycoord );

  /** Save element index if it is the closest one */
  vol = pt->qual;
  if( fabs(barycoord[0].val)*vol < *closestDist ) {
    *closestDist = fabs(barycoord[0].val)*vol;
    *closestTet = k;
  }

  return found;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ppt pointer to the point to locate
 * \param triaNormals non-normalized triangle normals of all mesh triangles
 * \param iTria pointer to the index of the found triangle
 * \param closestTria pointer to the index of the closest triangle
 * \param closestDist pointer to the closest distance
 *
 * \return 1 if found, 0 otherwise.
 *
 *  Exhaustive point search on the background triangles. If the point is not
 *  found, the triangle pointers points to the closest triangle.
 *
 */
int PMMG_locatePoint_exhaustTria( MMG5_pMesh mesh,MMG5_pPoint ppt,
                                  double *triaNormals,PMMG_barycoord *barycoord,
                                  int *iTria,int *closestTria,double *closestDist ) {
  MMG5_pTria     ptr;
  double         h;

  for( *iTria = 1; *iTria <= mesh->nt; (*iTria)++ ) {

    /* Increase step counter */
    ppt->s--;

    /** Get tetra */
    ptr = &mesh->tria[*iTria];
    if ( !MG_EOK(ptr) ) continue;

    /*¨Skip already analized tetras */
    if( ptr->flag == mesh->base ) continue;

    /** Exit the loop if you find the element */
    if( PMMG_locatePointInTria( mesh, ptr, *iTria, ppt,
                                &triaNormals[3*(*iTria)], barycoord,
                                &h, closestDist, closestTria ) ) break;

  }

  if( *iTria <= mesh->nt ) {
    return 1;
  } else {
    *iTria = *closestTria;
    /* Recompute barycentric coordinates to the closest point */
    PMMG_barycoord2d_getClosest( mesh,*iTria,ppt,barycoord );
    return 0;
  }
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ppt pointer to the point to locate
 * \param kfound pointer to the index of the starting element
 * \param triaNormals unit normals of the all triangles in the mesh
 * \param baryfound barycentric coordinates of the point in the found triangle
 *
 * \return 1 if the triangle has been updated, 0 otherwise.
 *
 *  Once the point is located in a candidate triangle, check whether its
 *  neighbours also see him (convex configurations) and update triangle index
 *  and barycentric coordinates if this is the case.
 *
 */
int PMMG_locatePoint_foundConvex( MMG5_pMesh mesh,MMG5_pPoint ppt,int *kfound,
                                  double *triaNormals,PMMG_barycoord *baryfound,
                                  double *h,double *closestDist,int *closestTria ) {
  MMG5_pTria ptr;
  PMMG_barycoord barycoord[4];
  int    *adjt,l,i,k,kmin,updated;
  double hmin;

#warning Luca: check distance computation
  adjt = &mesh->adjt[3*(*kfound-1)+1];
  hmin = *h;

  updated = 0;
  for( l = 0; l < 3; l++ ) {
    k = adjt[l]/3;
    if( ! k ) continue;

    ptr = &mesh->tria[k];
    if( !MG_EOK(ptr) ) continue;
    /* Visited triangles don't see the point or have already been listed here */
    if( ptr->flag == mesh->base ) continue;

    /** Exit the loop if you find the element */
    if( PMMG_locatePointInTria( mesh, ptr, k, ppt, &triaNormals[3*k],
                                barycoord, h, closestDist, closestTria ) ) {
      if( *h < hmin ) {
        updated = 1;
        hmin = *h;
        *kfound = k;
        for( i = 0; i < 4; i++ ) {
          baryfound[i].val = barycoord[i].val;
          baryfound[i].idx = barycoord[i].idx;
        }
      }
    }

  }

  return updated;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ppt pointer to the point to locate
 * \param triaNormals unit normals of the all triangles in the mesh
 * \param barycoord barycentric coordinates of the point to be located
 * \param iTria pointer to the index of the triangle
 * \param ifoundEdge pointer to the index of the local edge
 * \param ifoundVertex pointer to the index of the local vertex
 *
 * \return 0 if not found (closest), 1 if found, -1 if found through exhaustive
 * search.
 *
 *  Locate a point in a background mesh surface by traveling the triangles
 *  adjacency.
 *
 */
int PMMG_locatePointBdy( MMG5_pMesh mesh,MMG5_pPoint ppt,
                         double *triaNormals,int *nodeTrias,PMMG_barycoord *barycoord,
                         int *iTria,int *ifoundEdge,int *ifoundVertex ) {
  MMG5_pTria     ptr,ptr1;
  int            *adjt,i,k,k1,kprev,iprev,step,closestTria,stuck,backward;
  int            iloc;
  double         vol,eps,h,closestDist;
  static int     mmgWarn0=0,mmgWarn1=0;
  int            ier;

  if(!(*iTria))
    k = 1;
  else
    k = *iTria;

  assert( k <= mesh->nt );

  kprev = 0;
  stuck = 0;
  step = 0;
  ++mesh->base;

  closestTria = 0;
  closestDist = 1.0e10;

  *ifoundEdge   = PMMG_UNSET;
  *ifoundVertex = PMMG_UNSET;

  while( (step <= mesh->nt) && (!stuck) ) {
    step++;

    assert(kprev != k ) ;
    /** Get tria */
    ptr = &mesh->tria[k];
    if ( !MG_EOK(ptr) ) continue;

    /** Exit the loop if you find the element */
    if( PMMG_locatePointInTria( mesh, ptr, k, ppt, &triaNormals[3*k],
                                barycoord, &h, &closestDist, &closestTria ) ) {
      PMMG_barycoord_isBorder( barycoord, ifoundEdge, ifoundVertex );
      break;
    }

    /** Compute new direction */
    adjt = &mesh->adjt[3*(k-1)+1];
    kprev = k;
    for( i=0; i<3; i++ ) {
      k1 = adjt[barycoord[i].idx]/3;

      /* Skip if on boundary */
      if( !k1 ) continue;

      /* Test shadow regions if the tria has already been visited */
      ptr1 = &mesh->tria[k1];
      if(ptr1->flag == mesh->base) {
        iloc = PMMG_locatePointInWedge( mesh,ptr,kprev,iprev,ppt,barycoord );
        if( iloc == PMMG_UNSET ) continue;
        if( iloc == 4 ) {
          *ifoundEdge = iprev;
          ppt->s = step;
          *iTria = k;
          return 1;
        } else {
          ier = PMMG_locatePointInCone( mesh,nodeTrias,kprev,iloc,ppt );
          if( ier ) {
            *ifoundVertex = iloc;
            ppt->s = step;
            *iTria = k;
            return 1;
          }
        }
        continue;
      }

      /* Get next otherwise */
      k = k1;
      break;
    }
    iprev = barycoord[i].idx;

    /** Stuck: Start exhaustive research */
    if (i == 3) stuck = 1;

  }

  /* Store number of steps in the path for postprocessing */
  if( stuck )
    ppt->s = -step;
  else
    ppt->s = step;

  if( step > mesh->nt ) {
    /* Recompute barycentric coordinates to the closest point */
    *iTria = closestTria;
    PMMG_barycoord2d_getClosest( mesh,*iTria,ppt,barycoord );
    return 0;
  }


  /* If a candidate triangle has been found, check convex configurations */
  if( !stuck )
    PMMG_locatePoint_foundConvex( mesh,ppt,&k,triaNormals,barycoord,
                                  &h,&closestDist,&closestTria);

  /* Return the index of the tria */
  *iTria = k;

  /** Boundary hit or cyclic path: Perform exhaustive research */
  if( stuck ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      if ( mesh->info.imprim > PMMG_VERB_DETQUAL ) {
        fprintf(stderr,"\n  ## Warning %s: Cannot locate point,"
                " performing exhaustive research.\n",__func__);
      }
    }

    ier = PMMG_locatePoint_exhaustTria( mesh, ppt,triaNormals,barycoord,
                                           iTria,&closestTria,&closestDist );
    if( ier ) {
      return -1;
    } else {
    /** Element not found: Return the closest one */
      if ( !mmgWarn1 ) {
        mmgWarn1 = 1;
        if ( mesh->info.imprim > PMMG_VERB_VERSION ) {
          fprintf(stderr,"\n  ## Warning %s: Point not located, smallest external area %e.",
                  __func__,closestDist);
        }
      }
      return 0;
    }

  }

  return 1;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ppt pointer to the point to locate
 * \param faceAreas oriented face areas of the all tetrahedra in the mesh
 * \param idxTet pointer to the index of the found tetrahedron
 * \param closestTet pointer to the index of the closest tetrahedron
 * \param closestDist pointer to the distance from the closest tetrahedron
 *
 * \return 1 if found, 0 otherwise.
 *  Exhaustive point search on the background tetrahedra.
 *
 */
int PMMG_locatePoint_exhaustTetra( MMG5_pMesh mesh,MMG5_pPoint ppt,
                                   double *faceAreas,PMMG_barycoord *barycoord,
                                   int *idxTet,int *closestTet,double *closestDist ) {
  MMG5_pTetra    pt;
  double         vol;

  for( *idxTet = 1; *idxTet <= mesh->ne; (*idxTet)++ ) {

    /* Increase step counter */
    ppt->s--;

    /** Get tetra */
    pt = &mesh->tetra[*idxTet];
    if ( !MG_EOK(pt) ) continue;

    /*¨Skip already analized tetras */
    if( pt->flag == mesh->base ) continue;

    /** Exit the loop if you find the element */
    if( PMMG_locatePointInTetra( mesh, pt, *idxTet, ppt,&faceAreas[12*(*idxTet)],
                                 barycoord, closestDist, closestTet ) ) break;

  }

  if( *idxTet <= mesh->ne ) {
    return 1;
  } else {
    *idxTet = *closestTet;
    /* Recompute barycentric coordinates to the closest point */
    PMMG_barycoord3d_getClosest( mesh,*idxTet,ppt,barycoord );
    return 0;
  }
  return 1;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ppt pointer to the point to locate
 * \param init index of the starting element
 * \param faceAreas oriented face areas of the all tetrahedra in the mesh
 * \param barycoord barycentric coordinates of the point to be located
 * \param idxTet pointer to the index of the found tetrahedron.
 *
 * \return 0 if not found (closest), 1 if found, -1 if found through exhaustive
 * search.
 *
 *  Locate a point in a background mesh by traveling the elements adjacency.
 *
 */
int PMMG_locatePointVol( MMG5_pMesh mesh,MMG5_pPoint ppt,
                         double *faceAreas,PMMG_barycoord *barycoord,
                         int *idxTet ) {
  MMG5_pTetra    pt,pt1;
  int            *adja,iel,i,step,closestTet,stuck;
  double         vol,eps,closestDist;
  static int     mmgWarn0=0,mmgWarn1=0;
  int            ier;

  if(!(*idxTet))
    *idxTet = 1;

  assert( *idxTet <= mesh->ne );

  closestTet = 0;
  closestDist = 1.0e10;

  stuck = 0;
  step = 0;
  ++mesh->base;
  while( (step <= mesh->ne) && (!stuck) ) {
    step++;

    /** Get tetra */
    pt = &mesh->tetra[*idxTet];
    if ( !MG_EOK(pt) ) continue;

    /** Exit the loop if you find the element */
    if( PMMG_locatePointInTetra( mesh, pt, *idxTet,ppt,&faceAreas[12*(*idxTet)],
                                 barycoord,&closestDist,&closestTet ) ) break;

    /** Compute new direction (barycentric coordinates are sorted in increasing
     *  order) */
    adja = &mesh->adja[4*(*idxTet-1)+1];
    for( i=0; i<4; i++ ) {
      iel = adja[barycoord[i].idx]/4;

      /* Skip if on boundary */
      if (!iel) continue;

      /* Skip if already marked */
      pt1 = &mesh->tetra[iel];
      if(pt1->flag == mesh->base) continue;

      /* Get next otherwise */
      *idxTet = iel;
      break;
    }

    /** Stuck: Start exhaustive research */
    if (i == 4) stuck = 1;

  }

  /* Store number of steps in the path for postprocessing */
  if( stuck )
    ppt->s = -step;
  else
    ppt->s = step;

  if( step > mesh->ne ) {
    /* Recompute barycentric coordinates to the closest point */
    *idxTet = closestTet;
    PMMG_barycoord3d_getClosest( mesh,*idxTet,ppt,barycoord );
    return 0;
  }


  /** Boundary hit or cyclic path: Perform exhaustive research */
  if( stuck ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      if ( mesh->info.imprim > PMMG_VERB_DETQUAL ) {
        fprintf(stderr,"\n  ## Warning %s: Cannot locate point,"
                " performing exhaustive research.\n",__func__);
      }
    }

    ier = PMMG_locatePoint_exhaustTetra( mesh,ppt,faceAreas,barycoord,
                                         idxTet,&closestTet,&closestDist );

    if( ier ) {
      return -1;
    } else {
      /** Element not found: Return the closest one */
      if ( !mmgWarn1 ) {
        mmgWarn1 = 1;
        if ( mesh->info.imprim > PMMG_VERB_VERSION ) {
          fprintf(stderr,"\n  ## Warning %s: Point not located, smallest external volume %e.",
                  __func__,closestDist);
        }
      }
      return 0;
    }

  }
  return 1;
}

/**
 * \param mesh pointer to the current mesh structure
 * \param ip point index
 * \param ier error code
 * \param myrank process rank
 * \param igrp mesh group index
 *
 * Analise the found element and display warnings for exhaustive searches and
 * points not found.
 *
 */
void PMMG_locatePoint_errorCheck( MMG5_pMesh mesh,int ip,int ier,
                                 int myrank,int igrp ) {
  MMG5_pPoint ppt;

  ppt = &mesh->point[ip];

  if( !ier ) {
    fprintf(stderr,"\n  ## Warning: %s (rank %d, grp %d): closest element for"
            " point %d, coords %e %e %e\n",__func__,myrank,igrp,
            ip,ppt->tag,ppt->c[0],ppt->c[1],ppt->c[2]);
  } else if ( ier < 0 ) {
    fprintf(stderr,"\n  ## Warning: %s (rank %d, grp %d): exhaustive search for"
            " point %d, tag %d, coords %e %e %e\n",__func__,myrank,igrp,
            ip,ppt->tag,ppt->c[0],ppt->c[1],ppt->c[2]);
  }
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
    if( !MG_VOK(ppt) ) continue;
    if( ppt->tag & MG_BDY ) continue;
    ppt->s = meshOld->point[ppt->src].s;
    assert(ppt->s);
  }

}

/**
 * \param mesh pointer to the current mesh structure
 * \param meshOld pointer to the background mesh structure
 * \param locStats localization statistics structure
 *
 *  Compute localization statistics.
 *
 */
void PMMG_locate_postprocessing( MMG5_pMesh mesh,MMG5_pMesh meshOld,PMMG_locateStats *locStats ) {
  MMG5_pPoint ppt;
  int         ip,np;

  locStats->stepmin  = meshOld->ne;
  locStats->stepmax  = 0;
  locStats->stepav   = 0;
  locStats->nexhaust = 0;
  np = 0;

  /* Get the number of steps from the ppt->s field */
  for( ip = 1; ip <= mesh->np; ip++ ) {
    ppt = &mesh->point[ip];
    if( !MG_VOK(ppt) ) continue;
    if( ppt->tag & MG_BDY ) continue;
    /* Exhaustive searches are identified by a negative number of steps
     * (counting both the number of steps before and after the exhaustive
     * search) */
    if( ppt->s < 0 ) {
      locStats->nexhaust++;
      ppt->s *= -1;
    }
    np++;
    assert( ppt->s );
    if( ppt->s < locStats->stepmin ) locStats->stepmin = ppt->s;
    if( ppt->s > locStats->stepmax ) locStats->stepmax = ppt->s;
    locStats->stepav += ppt->s;
  }
  locStats->stepav *= 1.0/np;

}

/**
 * \param locStats localization statistics structure
 * \param ngrp number of meshes
 * \param myrank process rank
 *
 *  Print localization statistics.
 *
 */
void PMMG_locate_print( PMMG_locateStats *locStats,int ngrp,int myrank ) {
  int igrp;

  for( igrp = 0; igrp < ngrp; igrp++ )
    printf("         Localization report (rank %d grp %d): nexhaust %d max step %d, min step %d, av %f\n",myrank,igrp,locStats[igrp].nexhaust,locStats[igrp].stepmax,locStats[igrp].stepmin,locStats[igrp].stepav);
}

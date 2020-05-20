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
int PMMG_locatePointInCone( MMG5_pMesh mesh,int iel,int iloc,MMG5_pPoint ppt,
                            int *list ) {
  MMG5_pPoint    ppt0,ppt1;
  double         p[3],a[3],dist2,norm2,alpha;
  int            ilist,lon,jp,d,found;

  ppt0 = &mesh->point[mesh->tria[iel].v[iloc]];

  /* Mark point */
  ppt0->flag = mesh->base;

  /* Target point vector */
  for( d = 0; d < 3; d++ ) p[d] = ppt->c[d]-ppt0->c[d];
  dist2 = 0.0;
  for( d = 0; d < 3; d++ ) dist2 += p[d]*p[d];

  /* Get neighbour nodes */
  lon = MMG5_boulep(mesh,iel,iloc,mesh->adjt,list);

  /* Initialize as if the target is found in the cone */
  found = 1;

  /* Scan the list of neighbours */
  for( ilist = 1; ilist <= lon; ilist++ ){
    jp = list[ilist];
    ppt1 = &mesh->point[jp];
    /* Edge vector */
    for( d = 0; d < 3; d++ ) a[d] = ppt1->c[d]-ppt0->c[d];
    norm2 = 0.0;
    for( d = 0; d < 3; d++ ) norm2 += a[d]*a[d];
    /* Rough check on maximum distance */
    if( dist2 >= norm2 ) {
      found = 0;
      break;
    }
    /* Scalar product of the target vector with the edge vector */
    alpha = 0.0;
    for( d = 0; d < 3; d++ ) alpha += a[d]*p[d];
    if( alpha > 0.0 ) {
      found = 0;
      break;
    }
  }

  return found;
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
  double      a[3],p[3],norm2,dist2,alpha;
  int         i0,i1,d;

  /* Check nodal sides */
  i0 = MMG5_inxt2[l];
  i1 = MMG5_iprv2[l];
  ppt0 = &mesh->point[ptr->v[i0]];
  ppt1 = &mesh->point[ptr->v[i1]];

  /* Target point vector */
  for( d = 0; d < 3; d++ ) p[d] = ppt->c[d]-ppt0->c[d];
  dist2 = 0.0;
  for( d = 0; d < 3; d++ ) dist2 += p[d]*p[d];

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
  dist2 = 0.0;
  for( d = 0; d < 3; d++ ) dist2 += p[d]*p[d];
  if( dist2 >= norm2 ) return -1;

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
  barycoord->val = alpha/norm2;

  return 4;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ptr pointer to the triangle to analyze
 * \param k index of the triangle
 * \param ppt pointer to the point to locate
 * \param triaNormal unit normal of the current triangle
 * \param barycoord barycentric coordinates of the point to be located
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
                            double *closestDist,int *closestTria ) {
  MMG5_pPoint    ppt0,ppt1;
  double         h,hmax;
  double         proj,dist[3];
  int            j,d,found;

  /* Mark tria */
  ptr->flag = mesh->base;

  /* Evaluate point in tetra through barycentric coordinates */
  found = PMMG_barycoord2d_evaluate( mesh,ptr,k,ppt->c,triaNormal,barycoord );


  /* Rough check on the distance from the surface
   * (avoid being on the opposite pole) */
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

  /* Distance from the first vertex */
  ppt0 = &mesh->point[ptr->v[0]];
  for( d = 0; d < 3; d++ )
    dist[d] = ppt->c[d]-ppt0->c[d];

  /* Project */
  proj = 0;
  for( d = 0; d < 3; d++ )
    proj += dist[d]*triaNormal[d];

  /* Orthogonalize to get distance in the plane */
  for( d = 0; d < 3; d++ )
    dist[d] -= proj*triaNormal[d];

  h = 0.0;
  for( d = 0; d < 3; d++ )
    h += dist[d]*dist[d];
  h = sqrt(h);

  /* Save element index (with negative sign) if it is the closest one */
  if( h < *closestDist ) {
    *closestDist = h;
    *closestTria = -k;
  }

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
                             double *faceAreas,PMMG_barycoord *barycoord ) {

  /* Mark tetra */
  pt->flag = mesh->base;

  /* Evaluate point in tetra through barycentric coordinates */
  return PMMG_barycoord3d_evaluate( mesh,pt,ppt->c,faceAreas,barycoord );
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ppt pointer to the point to locate
 * \param triaNormals non-normalized triangle normals of all mesh triangles
 * \param closestTria index of the closest triangle
 *
 * \return k index of the target element; if higher than the number of
 *           triangles, the point has not been found.
 *
 *  Exhaustive point search on the background triangles.
 *
 */
int PMMG_locatePoint_exhaustTria( MMG5_pMesh mesh,MMG5_pPoint ppt,
                                  double *triaNormals,PMMG_barycoord *barycoord,
                                  double *closestDist,int *closestTria ) {
  MMG5_pTria     ptr;
  int            k;

  for( k = 1; k <= mesh->nt; k++ ) {

    /* Increase step counter */
    ppt->s--;

    /** Get tetra */
    ptr = &mesh->tria[k];
    if ( !MG_EOK(ptr) ) continue;

    /*¨Skip already analized tetras */
    if( ptr->flag == mesh->base ) continue;

    /** Exit the loop if you find the element */
    if( PMMG_locatePointInTria( mesh, ptr, k, ppt,
                                &triaNormals[3*k], barycoord,
                                closestDist, closestTria ) ) break;

  }

  return k;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ppt pointer to the point to locate
 * \param init index of the starting element
 * \param triaNormals unit normals of the all triangles in the mesh
 * \param barycoord barycentric coordinates of the point to be located
 * \param ip current mesh point
 *
 * \return ie if positive, index of the target element; if negative, index of
 * the closest element; 0 if not found
 *
 *  Locate a point in a background mesh surface by traveling the triangles
 *  adjacency.
 *
 */
int PMMG_locatePointBdy( MMG5_pMesh mesh,MMG5_pPoint ppt,int init,
                         double *triaNormals,PMMG_barycoord *barycoord,int ip,
                         int *foundWedge,int *foundCone ) {
  MMG5_pTria     ptr,ptr1;
  int            *adjt,i,k,k1,kprev,iprev,step,closestTria,stuck,backward;
  int            list[MMG5_LMAX],iloc;
  double         vol,eps,closestDist;
  static int     mmgWarn0=0,mmgWarn1=0;

  if(!init)
    k = 1;
  else
    k = init;

  assert( k <= mesh->nt );

  kprev = 0;
  stuck = 0;
  step = 0;
  ++mesh->base;

  closestTria = 0;
  closestDist = 1.0e10;

  *foundWedge = -1;
  *foundCone  = -1;

  while( (step <= mesh->nt) && (!stuck) ) {
    step++;

    /** Get tria */
    ptr = &mesh->tria[k];
    if ( !MG_EOK(ptr) ) continue;

    /** Exit the loop if you find the element */
    if( PMMG_locatePointInTria( mesh, ptr, k, ppt, &triaNormals[3*k],
                                barycoord, &closestDist, &closestTria ) ) {
      break;
    }

    /** Compute new direction */
    adjt = &mesh->adjt[3*(k-1)+1];
    kprev = k;
    for( i=0; i<3; i++ ) {
      k1 = adjt[barycoord[i].idx]/3;

      /* Skip if on boundary */
      if( !k1 ) continue;

      if( (k1 == kprev) && !i && barycoord[i+1].val > -MMG5_EPS  ) {
        /* You should check it only if the tria check tells you to ONLY go back! */
        iloc = PMMG_locatePointInWedge( mesh,ptr,k,barycoord[i].idx,ppt,barycoord );
        if( iloc == 4 ) {
          *foundWedge = barycoord[i].idx;
          break;
        } else if( iloc > -1 ) {
          if( PMMG_locatePointInCone( mesh,kprev,iloc,ppt,list ) ) {
            /* Return triangle index and local node index */
            *foundCone = iloc;
            break;
          }
        }
      }

      /* Skip if already marked */
      ptr1 = &mesh->tria[k1];
      if(ptr1->flag == mesh->base) continue;

      /* Get next otherwise */
      k = k1;
      break;
    }
    iprev = barycoord[i].idx;

    /** Stuck: Start exhaustive research */
    if (i == 3) stuck = 1;

  }

  /* Store number of steps in the path for postprocessing */
  if( stuck ) step *= -1;
  ppt->s = step;

  /** Boundary hit or cyclic path: Perform exhaustive research */
  if( stuck ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      if ( mesh->info.imprim > PMMG_VERB_DETQUAL ) {
        fprintf(stderr,"\n  ## Warning %s: Cannot locate point,"
                " performing exhaustive research.\n",__func__);
      }
    }

    k = PMMG_locatePoint_exhaustTria( mesh, ppt,triaNormals,barycoord,
                                      &closestDist,&closestTria );
    /** Element not found: Return the closest one with negative sign (if found) */
    if ( k == mesh->nt+1 ) {
      if ( !mmgWarn1 ) {
        mmgWarn1 = 1;
        if ( mesh->info.imprim > PMMG_VERB_VERSION ) {
          fprintf(stderr,"\n  ## Warning %s: Point not located, smallest external area %e.",
                  __func__,closestDist);
        }
      }
      /* Recompute barycentric coordinates to the closest point */
      PMMG_barycoord2d_getClosest( mesh,-closestTria,ppt,barycoord );
      return closestTria;
    }

  }
  return k;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ppt pointer to the point to locate
 * \param faceAreas oriented face areas of the all tetrahedra in the mesh
 * \param closestTet index of the closest tetrahedron
 *
 * \return ie index of the target element; if higher than the number of
 *            tetrahedra, the point has not been found.
 *
 *  Exhaustive point search on the background tetrahedra.
 *
 */
int PMMG_locatePoint_exhaustTetra( MMG5_pMesh mesh,MMG5_pPoint ppt,
                                   double *faceAreas,PMMG_barycoord *barycoord,
                                   int *closestTet ) {
  MMG5_pTetra    pt;
  int            ie;
  double         vol,closestDist;

  *closestTet = 0;
  closestDist = 1.0e10;
  for( ie = 1; ie <= mesh->ne; ie++ ) {

    /* Increase step counter */
    ppt->s--;

    /** Get tetra */
    pt = &mesh->tetra[ie];
    if ( !MG_EOK(pt) ) continue;

    /*¨Skip already analized tetras */
    if( pt->flag == mesh->base ) continue;

    /** Exit the loop if you find the element */
    if( PMMG_locatePointInTetra( mesh, pt, ppt,&faceAreas[12*ie],
                                 barycoord ) ) break;

    /** Save element index (with negative sign) if it is the closest one */
    vol = pt->qual;
    if( fabs(barycoord[0].val)*vol < closestDist ) {
      closestDist = fabs(barycoord[0].val)*vol;
      *closestTet = -ie;
    }
  }

  return ie;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ppt pointer to the point to locate
 * \param init index of the starting element
 * \param faceAreas oriented face areas of the all tetrahedra in the mesh
 * \param barycoord barycentric coordinates of the point to be located
 * \param ip current mesh point
 *
 * \return ie if positive, index of the target element; if negative, index of
 * the closest element; 0 if not found
 *
 *  Locate a point in a background mesh by traveling the elements adjacency.
 *
 */
int PMMG_locatePointVol( MMG5_pMesh mesh,MMG5_pPoint ppt,int init,
                         double *faceAreas,PMMG_barycoord *barycoord,
                         int ip ) {
  MMG5_pTetra    pt,pt1;
  int            *adja,iel,i,idxTet,step,closestTet,stuck;
  double         vol,eps,closestDist;
  static int     mmgWarn0=0,mmgWarn1=0;

  if(!init)
    idxTet = 1;
  else
    idxTet = init;

  assert( idxTet <= mesh->ne );

  stuck = 0;
  step = 0;
  ++mesh->base;
  while( (step <= mesh->ne) && (!stuck) ) {
    step++;

    /** Get tetra */
    pt = &mesh->tetra[idxTet];
    if ( !MG_EOK(pt) ) continue;

    /** Exit the loop if you find the element */
    if( PMMG_locatePointInTetra( mesh, pt, ppt,&faceAreas[12*idxTet],
                                 barycoord ) ) break;

    /** Compute new direction (barycentric coordinates are sorted in increasing
     *  order) */
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
    if (i == 4) stuck = 1;

  }

  /* Store number of steps in the path for postprocessing */
  if( stuck ) step *= -1;
  ppt->s = step;

  /** Boundary hit or cyclic path: Perform exhaustive research */
  if( stuck ) {
    if ( !mmgWarn0 ) {
      mmgWarn0 = 1;
      if ( mesh->info.imprim > PMMG_VERB_DETQUAL ) {
        fprintf(stderr,"\n  ## Warning %s: Cannot locate point,"
                " performing exhaustive research.\n",__func__);
      }
    }

    idxTet = PMMG_locatePoint_exhaustTetra( mesh,ppt,faceAreas,barycoord,
                                            &closestTet );

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
 * \param ip point index
 * \param kfound index of the found element
 * \param myrank process rank
 * \param igrp mesh group index
 *
 * \return kfound for a valid element, 0 otherwise
 *
 * Analise the found element and display warnings for exhaustive searches.
 *
 */
int PMMG_locatePoint_errorCheck( MMG5_pMesh mesh,int ip,int kfound,
                                 int myrank,int igrp ) {
  MMG5_pPoint ppt;

  ppt = &mesh->point[ip];

  if( !kfound ) {
    fprintf(stderr,"\n  ## Error: %s (rank %d, grp %d): point %d not found, tag %d,"
            " coords %e %e %e\n",__func__,myrank,igrp,
            ip,ppt->tag,ppt->c[0],ppt->c[1],ppt->c[2]);
  } else if( kfound < 0 ) {
    if ( mesh->info.imprim > PMMG_VERB_VERSION ) {
      fprintf(stderr,"\n  ## Warning: %s (rank %d, grp %d): closest element for"
              " point %d, coords %e %e %e\n",__func__,myrank,igrp,
              ip,ppt->tag,ppt->c[0],ppt->c[1],ppt->c[2]);
    }
    return -kfound;
  }

  return kfound;
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

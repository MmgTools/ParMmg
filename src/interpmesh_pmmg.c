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

void PMMG_storeScalingParam( PMMG_pParMesh parmesh,int igrp ) {
  MMG5_pMesh mesh, meshOld;
  int i;

  mesh    = parmesh->listgrp[igrp].mesh;
  meshOld = parmesh->old_listgrp[igrp].mesh;

  meshOld->info.delta = mesh->info.delta;
  for( i = 0; i < mesh->dim; i++ )
    meshOld->info.min[i] = mesh->info.min[i];

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
int PMMG_compute_baryCoord( MMG5_pMesh mesh, MMG5_pTetra pt,
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
 * \param ptr pointer to the tetra to analyze
 * \param ppt pointer to the point to locate
 * \param faceAreas oriented face areas of the current tetrahedra
 * \param barycoord barycentric coordinates of the point to be located
 *
 * \return 1 if found; 0 if not found
 *
 *  Locate a point in a background tetrahedron, and provide its barycentric
 *  coordinates..
 *
 */
int PMMG_locatePointInTetra( MMG5_pMesh mesh, MMG5_pTetra ptr, MMG5_pPoint ppt,
                             double *faceAreas, PMMG_baryCoord *barycoord ) {
  double         vol,eps;
  int            found = 0;

  vol = ptr->qual;
  eps = MMG5_EPS;

  /** Mark tetra */
  ptr->flag = mesh->base;

  /** Get barycentric coordinates and sort them in ascending order */
  PMMG_compute_baryCoord(mesh, ptr, ppt->c, faceAreas, barycoord);
  qsort(barycoord,4,sizeof(PMMG_baryCoord),PMMG_compare_baryCoord);

  /** Exit if inside the element */
  if( barycoord[0].val > -eps ) found = 1;

  return found;
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
int PMMG_locatePoint( MMG5_pMesh mesh, MMG3D_pPROctree q, int useOctree,
                      MMG5_pPoint ppt, int init,
                      double *faceAreas, PMMG_baryCoord *barycoord ) {
  MMG5_pTetra    ptr,pt1;
  int            *adja,iel,i,idxTet,step,closestTet;
  double         vol,eps,closestDist,coord[3];
  static int     mmgWarn0=0,mmgWarn1=0;
  int64_t        zcoord;
  int            ileaf,nleaves,*leaves;

  if( useOctree ) {
    /* Compute z-ordered coordinate */
    for( i = 0; i < mesh->dim; i++ )
      coord[i] = ( 1.0 / mesh->info.delta )*( ppt->c[i] - mesh->info.min[i] );
    zcoord = MMG3D_getPROctreeCoordinate( q, coord, mesh->dim );
    nleaves = PMMG_getPROctree_leaves( q, zcoord, &leaves );
    if( !nleaves ) return 0;
  }

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
    if( useOctree ) {

      for( ileaf=0; ileaf<nleaves; ileaf++ ) {

        /** Get tetra */
        idxTet = leaves[ileaf]/4;
        ptr = &mesh->tetra[idxTet];
        if ( !MG_EOK(ptr) ) continue;
  
        /** Exit the loop if you find the element */
        if( PMMG_locatePointInTetra( mesh, ptr, ppt,&faceAreas[12*idxTet],
                                     barycoord ) ) break;
  
  
        /** Save element index (with negative sign) if it is the closest one */
        if( fabs(barycoord[0].val)*vol < closestDist ) {
          closestDist = fabs(barycoord[0].val)*vol;
          closestTet = -idxTet;
        }

      }

      /** Element not found: Return the closest one with negative sign (if found) */
      if ( ileaf == nleaves ) {
        if ( !mmgWarn1 ) {
          mmgWarn1 = 1;
          if ( mesh->info.imprim > PMMG_VERB_VERSION ) {
            fprintf(stderr,"\n  ## Warning %s: Point not located, smallest external volume %e.",
                    __func__,closestDist);
          }
        }
        return closestTet;
      }
    } else {
      for( idxTet=1; idxTet<=mesh->ne; idxTet++ ) {

        /** Get tetra */
        ptr = &mesh->tetra[idxTet];
        if ( !MG_EOK(ptr) ) continue;
  
        /** Exit the loop if you find the element */
        if( PMMG_locatePointInTetra( mesh, ptr, ppt,&faceAreas[12*idxTet],
                                     barycoord ) ) break;


        /** Save element index (with negative sign) if it is the closest one */
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

  }
  
  return idxTet;
}

/**
 * \param grp pointer to the current group structure
 * \param oldGrp pointer to the bqckground group structure
 * \param pt pointer to the target background tetrahedron
 * \param ip index of the current point
 * \param phi barycentric coordinates of the point to be interpolated
 *
 * \return 0 if fail, 1 if success
 *
 *  Linearly interpolate point metrics on a target background tetrahedron..
 *
 */
int PMMG_interpMetrics_point( PMMG_pGrp grp,PMMG_pGrp oldGrp,MMG5_pTetra pt,
                              int ip,PMMG_baryCoord *phi ) {
  MMG5_pMesh     mesh;
  MMG5_pSol      met,oldMet;
  MMG5_pPoint    ppt;
  int            iloc,i,isize,nsize,ier;

  met    = grp->met;
  oldMet = oldGrp->met;
  mesh   = grp->mesh;
  ppt    = &mesh->point[ip];
  nsize  = met->size;

  /** Linear interpolation of the metrics */
  for( isize = 0; isize<nsize; isize++ ) {
    met->m[nsize*ip+isize] = 0.0;
    /* Barycentric coordinates could be permuted */
    for( i=0; i<4; i++ ) {
      iloc = phi[i].idx;
      met->m[nsize*ip+isize] += phi[i].val*oldMet->m[nsize*pt->v[iloc]+isize];
    }
  }

  return 1;
}

int PMMG_interpMetrics_grp( PMMG_pParMesh parmesh,MMG3D_pPROctree *qgrps,int igrp ) {
  PMMG_pGrp   grp,oldGrp;
  MMG5_pMesh  mesh,oldMesh;
  MMG3D_pPROctree q;
  MMG5_pTetra pt;
  MMG5_pPoint ppt;
  PMMG_baryCoord barycoord[4];
  double      *faceAreas,*normal;
  int         ip,istart,ie,ifac,ia,ib,ic,iloc;
  int         ier;
  static int  mmgWarn=0;


  grp = &parmesh->listgrp[igrp];
  mesh = grp->mesh;
  oldGrp = &parmesh->old_listgrp[igrp];
  oldMesh = oldGrp->mesh;
  q = qgrps[igrp];

  if( mesh->info.inputMet != 1 ) {

    /* Nothing to do */
    return 1;

  } else {

    if( mesh->info.hsiz > 0.0 ) {

      /* Compute constant metrics */
      if ( !MMG3D_Set_constantSize(mesh,grp->met) ) return 0;

    } else {

      /** Pre-compute oriented face areas */   
      ier = 1;
      PMMG_MALLOC( parmesh,faceAreas,12*(oldMesh->ne+1),double,"faceAreas",ier=0 );
      if( !ier ) {
        PMMG_DEL_MEM(parmesh,faceAreas,double,"faceAreas");
        return 0;
      }
  
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

          /** Locate point in the old mesh */
          istart = PMMG_locatePoint( oldMesh, q, parmesh->info.PROctree_mode,
                                     ppt, istart,
                                     faceAreas, barycoord );
          if( !istart ) {
            fprintf(stderr,"\n  ## Error: %s: proc %d (grp %d),"
                    " point %d not found, coords %e %e %e\n",__func__,
                    parmesh->myrank,igrp,ip, mesh->point[ip].c[0],
                    mesh->point[ip].c[1],mesh->point[ip].c[2]);
            return 0;
          } else if( istart < 0 ) {
            if ( !mmgWarn ) {
              mmgWarn = 1;
              if ( mesh->info.imprim > PMMG_VERB_VERSION ) {
                fprintf(stderr,"\n  ## Warning: %s: proc %d (grp %d), point %d not"
                        " found, coords %e %e %e\n",__func__,parmesh->myrank,
                        igrp,ip, mesh->point[ip].c[0],mesh->point[ip].c[1],
                        mesh->point[ip].c[2]);
              }
            }
            istart = -istart;
          }

          /** Interpolate point metrics */
          ier = PMMG_interpMetrics_point(grp,oldGrp,&oldMesh->tetra[istart],
                                         ip,barycoord);

          /* Flag point as interpolated */
          ppt->flag = mesh->base;
        }
      }
    }
  }

  PMMG_DEL_MEM( parmesh,faceAreas,double,"faceAreas");
  return 1;
}

/**
 * \param parmesh pointer to the parmesh structure
 *
 * \return 0 if fail, 1 if success
 *
 *  Interpolate metrics for all groups from background to current meshes.
 *  Do nothing if no metrics is provided (info.inputMet == 0), otherwise:
 *  - if the metrics is constant, recompute it;
 *  - else, interpolate the non-constant metrics.
 *
 */
int PMMG_interpMetrics( PMMG_pParMesh parmesh, MMG3D_pPROctree *q ) {
  int         igrp,ier;

  /** Loop on current groups */
  for( igrp = 0; igrp < parmesh->ngrp; igrp++ )
    if( !PMMG_interpMetrics_grp( parmesh, q, igrp ) )
      return 0;

  return 1;
}

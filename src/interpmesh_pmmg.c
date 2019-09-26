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
  int            *adja,iel,ip,idxTet,step,closestTet;
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

    adja = &mesh->adja[4*(idxTet-1)+1];
    vol = ptr->qual;
    eps = MMG5_EPS;

    /** Mark tetra */
    ptr->flag = mesh->base;

    /** Get barycentric coordinates and sort them in ascending order */
    PMMG_compute_baryCoord(mesh, ptr, ppt->c, &faceAreas[12*idxTet], barycoord);
    qsort(barycoord,4,sizeof(PMMG_baryCoord),PMMG_compare_baryCoord);

    /** Exit if inside the element */
    if( barycoord[0].val > -eps ) break;

    /** Compute new direction */
    for( ip=0; ip<4; ip++ ) {
      iel = adja[barycoord[ip].idx]/4;

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
    if (ip == 4) step = mesh->ne+1;

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

      adja = &mesh->adja[4*(idxTet-1)+1];
      vol = ptr->qual;
      eps = MMG5_EPS;

      /** Mark tetra */
      ptr->flag = mesh->base;

      /** Get barycentric coordinates and sort them in ascending order */
      PMMG_compute_baryCoord(mesh, ptr, ppt->c, &faceAreas[12*idxTet], barycoord);
      qsort(barycoord,4,sizeof(PMMG_baryCoord),PMMG_compare_baryCoord);

      /** Exit if inside the element */
      if( barycoord[0].val > -eps ) break;

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

/**
 * \param grp pointer to the current group structure.
 * \param oldGrp pointer to the bqckground group structure.
 * \param permNodGlob permutation array for nodes.
 *
 * \return 0 if fail, 1 if success
 *
 * Copy the metric of a freezed interface point.
 *
 */
int PMMG_copyMetrics_point( PMMG_pGrp grp,PMMG_pGrp oldGrp, int* permNodGlob) {
  MMG5_pMesh     mesh,oldMesh;
  MMG5_pSol      met,oldMet;
  MMG5_pPoint    ppt;
  int            isize,nsize,ip;

  mesh    = grp->mesh;
  met     = grp->met;
  oldMesh = oldGrp->mesh;
  oldMet  = oldGrp->met;
  nsize   = met->size;

#warning Luca: when surface adapt will be ready, distinguish BDY from PARBDY

  /** Freezed points: Copy the  metrics  */
  if ( (!oldGrp->mesh->info.renum) || !permNodGlob ) {
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
int PMMG_interpMetrics_grps( PMMG_pParMesh parmesh,int *permNodGlob ) {
  PMMG_pGrp   grp,oldGrp;
  MMG5_pMesh  mesh,oldMesh;
  MMG5_pTetra pt;
  MMG5_pPoint ppt;
  PMMG_baryCoord barycoord[4];
  double      **faceAreas,*normal;
  int         igrp,ip,istart,ie,ifac,ia,ib,ic,iloc;
  int         ier;
  static int  mmgWarn=0;

  /** Pre-compute oriented face areas */
  PMMG_MALLOC( parmesh,faceAreas,parmesh->nold_grp,double*,"faceAreas pointer",return 0);

  for( igrp = 0; igrp < parmesh->nold_grp; igrp++ ) {
    grp = &parmesh->old_listgrp[igrp];
    mesh = grp->mesh;

    ier = 1;
    PMMG_MALLOC( parmesh,faceAreas[igrp],12*(mesh->ne+1),double,"faceAreas",ier=0 );
    if( !ier ) {
      int igrp1;
      for( igrp1 = 0; igrp1 < igrp; igrp1++ )
        PMMG_DEL_MEM(parmesh,faceAreas[igrp1],double,"faceAreas");
      PMMG_DEL_MEM(parmesh,faceAreas,double*,"faceAreas pointer");
      return 0;
    }

    for( ie = 1; ie <= mesh->ne; ie++ ) {
      pt = &mesh->tetra[ie];
      /* Store tetra volume in the qual field */
      pt->qual = MMG5_orvol( mesh->point, pt->v );
      /* Store oriented face normals */
      for( ifac = 0; ifac < 4; ifac++ ) {
        normal = &faceAreas[igrp][12*ie+3*ifac];
        ia = pt->v[MMG5_idir[ifac][0]];
        ib = pt->v[MMG5_idir[ifac][1]];
        ic = pt->v[MMG5_idir[ifac][2]];
        ier = MMG5_nonUnitNorPts( mesh,ia,ib,ic,normal );
      }
    }
  }


  /** Loop on current groups */
  for( igrp = 0; igrp < parmesh->ngrp; igrp++ ) {
    grp = &parmesh->listgrp[igrp];
    mesh = grp->mesh;

    if( mesh->info.inputMet != 1 ) {

      /* Nothing to do */
      continue;

    } else {

      if( mesh->info.hsiz > 0.0 ) {

        /* Compute constant metrics */
        if ( !MMG3D_Set_constantSize(mesh,grp->met) ) return 0;

      } else {

        /* Interpolate metrics */
        oldGrp = &parmesh->old_listgrp[igrp];
        oldMesh = oldGrp->mesh;

        oldMesh->base = 0;
        for ( ie = 1; ie < oldMesh->ne+1; ie++ ) {
          pt = &oldMesh->tetra[ie];
          if ( !MG_EOK(pt) ) continue;
          pt->flag = oldMesh->base;
        }

        mesh->base++;
        istart = 1;

        /* Copy the metric of interface points */
        ier = PMMG_copyMetrics_point( grp,oldGrp,permNodGlob );

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
              continue; // treated by copyMetric_points
            } else {

              /** Locate point in the old mesh */
              istart = PMMG_locatePoint( oldMesh, ppt, istart,
                                         faceAreas[igrp], barycoord );
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
            }

            /* Flag point as interpolated */
            ppt->flag = mesh->base;
          }
        }
      }
    }

  }
  for( igrp = 0; igrp < parmesh->nold_grp; igrp++)
    PMMG_DEL_MEM( parmesh,faceAreas[igrp],double,"faceAreas");
  PMMG_DEL_MEM( parmesh,faceAreas,double*,"faceAreas pointer");
  return 1;
}

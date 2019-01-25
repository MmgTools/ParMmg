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
 * \param barycoord pointer to the point barycentric coordinates in the current
 * tetra
 *
 * \return 0 if fail, 1 if success
 *
 *  Compute the barycentric coordinates of a given point in a given tetrahedron.
 *
 */
int PMMG_compute_baryCoord( MMG5_pMesh mesh, MMG5_pTetra pt,
                    double *coord, PMMG_baryCoord *barycoord ) {
  double *c0,*c1,*c2,*c3,vol;

  vol = MMG5_orvol( mesh->point, pt->v );

  c0 = mesh->point[pt->v[0]].c;
  c1 = mesh->point[pt->v[1]].c;
  c2 = mesh->point[pt->v[2]].c;
  c3 = mesh->point[pt->v[3]].c;

  barycoord[0].val = MMG5_det4pt( coord, c1,    c2,    c3    )/vol;
  barycoord[1].val = MMG5_det4pt( c0,    coord, c2,    c3    )/vol;
  barycoord[2].val = MMG5_det4pt( c0,    c1,    coord, c3    )/vol;
  barycoord[3].val = MMG5_det4pt( c0,    c1,    c2,    coord )/vol;

  barycoord[0].idx = 0;
  barycoord[1].idx = 1;
  barycoord[2].idx = 2;
  barycoord[3].idx = 3;

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
 *
 * \return ie if positive, index of the target element; if negative, index of
 * the closest element; 0 if not found
 *
 *  Locate a point in a background mesh by traveling the elements adjacency.
 *
 */
int PMMG_locatePoint( MMG5_pMesh mesh, MMG5_pPoint ppt, int init ) {
  MMG5_pTetra    ptr,pt1;
  PMMG_baryCoord barycoord[4];
  int            *adja,iel,ip,idxTet,step,closestTet;
  double         vol,eps,closestDist;

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
    adja = &mesh->adja[4*(idxTet-1)+1];
    vol = MMG5_orvol( mesh->point, ptr->v );
    eps = MMG5_EPS;

    /** Mark tetra */
    ptr->flag = mesh->base;
 
    /** Get barycentric coordinates and sort them in ascending order */
    PMMG_compute_baryCoord(mesh, ptr, ppt->c, barycoord);
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
    fprintf(stderr,"\n ## Warning %s: Cannot locate point, performing exhaustive research.",__func__);

    closestTet = 0;
    closestDist = 1.0e10;
    for( idxTet=1; idxTet<mesh->ne+1; idxTet++ ) {
 
      /** Get tetra */
      ptr = &mesh->tetra[idxTet];
      adja = &mesh->adja[4*(idxTet-1)+1];
      vol = MMG5_orvol( mesh->point, ptr->v );
      eps = MMG5_EPS;

      /** Mark tetra */
      ptr->flag = mesh->base;
 
      /** Get barycentric coordinates and sort them in ascending order */
      PMMG_compute_baryCoord(mesh, ptr, ppt->c, barycoord);
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
      fprintf(stderr,"\n ## Warning %s: Point not located, smallest external volume %e.",closestDist,__func__);
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
 *
 * \return 0 if fail, 1 if success
 *
 *  Linearly interpolate point metrics on a target background tetrahedron..
 *
 */
int PMMG_interpMetrics_point( PMMG_pGrp grp,PMMG_pGrp oldGrp,MMG5_pTetra pt,int ip ) {
  MMG5_pMesh     mesh;
  MMG5_pSol      met,oldMet;
  MMG5_pPoint    ppt;
  PMMG_baryCoord phi[4];
  int            iloc,isize,nsize,ier;

  met    = grp->met;
  oldMet = oldGrp->met;
  mesh   = grp->mesh;
  ppt    = &mesh->point[ip];
  nsize  = met->size;

  /** Get barycentric coordinates **/
  ier = PMMG_compute_baryCoord( oldGrp->mesh, pt, ppt->c, phi );

  /** Linear interpolation of the metrics */
  for( isize = 0; isize<nsize; isize++ ) {
    met->m[nsize*ip+isize] = 0.0;
    for( iloc=0; iloc<4; iloc++ ) {
      met->m[nsize*ip+isize] += phi[iloc].val*oldMet->m[nsize*pt->v[iloc]+isize];
    }
  }

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
int PMMG_interpMetrics_grps( PMMG_pParMesh parmesh ) {
  PMMG_pGrp   grp,oldGrp;
  MMG5_pMesh  mesh,oldMesh;
  MMG5_pTetra pt;
  double      barycoord[4],vol,eps;
  int         *adja;
  int         igrp,ip,ie,ier;


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
        for ( ie = 1; ie < oldMesh->ne+1; ie++ )
          oldMesh->tetra[ie].flag = oldMesh->base;


        ie = 1;
        for( ip=1; ip<mesh->np+1; ip++ ) {
          if( !MG_VOK(&mesh->point[ip]) ) continue;
    
          /** Locate point in the old mesh */
          ie = PMMG_locatePoint( oldMesh, &mesh->point[ip], ie );
          if( !ie ) {
            fprintf(stderr,"\n  ## Error: %s: proc %d (grp %d), point %d not found, coords %e %e %e\n",__func__,parmesh->myrank,igrp,ip, mesh->point[ip].c[0],mesh->point[ip].c[1],mesh->point[ip].c[2]);
            return 0;
          } else if( ie < 0 ) {
            fprintf(stderr,"\n  ## Warning: %s: proc %d (grp %d), point %d not found, coords %e %e %e\n",__func__,parmesh->myrank,igrp,ip, mesh->point[ip].c[0],mesh->point[ip].c[1],mesh->point[ip].c[2]);
            ie = -ie;
          }
    
          pt = &oldMesh->tetra[ie];
    
          /** Interpolate point metrics */
          ier = PMMG_interpMetrics_point(grp,oldGrp,pt,ip);
        }

      }
    }
  
  }
  return 1;
}

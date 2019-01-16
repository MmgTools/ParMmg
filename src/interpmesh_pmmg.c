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

/**
 * \param nelem number of elements in the initial group
 * \param target_mesh_size wanted number of elements per group
 *
 * \return the needed number of groups
 *
 *  Compute the needed number of groups to create groups of \a target_mesh_size
 *  elements from a group of nelem elements.
 *
 */
int PMMG_baryCoord( MMG5_pMesh mesh, MMG5_pTetra pt,
                    double *coord, double *barycoord ) {
  double *c0,*c1,*c2,*c3,vol;

  vol = MMG5_orvol( mesh->point, pt->v );

  c0 = mesh->point[pt->v[0]].c;
  c1 = mesh->point[pt->v[1]].c;
  c2 = mesh->point[pt->v[2]].c;
  c3 = mesh->point[pt->v[3]].c;

  barycoord[0] = MMG5_det4pt( coord, c1,    c2,    c3    )/vol;
  barycoord[1] = MMG5_det4pt( c0,    coord, c2,    c3    )/vol;
  barycoord[2] = MMG5_det4pt( c0,    c1,    coord, c3    )/vol;
  barycoord[3] = MMG5_det4pt( c0,    c1,    c2,    coord )/vol;

  return 1;
}

int PMMG_interpMetrics_point( PMMG_pGrp grp,PMMG_pGrp oldGrp,MMG5_pTetra pt,int ip ) {
  MMG5_pMesh  mesh;
  MMG5_pSol   met,oldMet;
  MMG5_pPoint ppt;
  double      phi[4];
  int         iloc,isize,nsize,ier;

  met    = grp->met;
  oldMet = oldGrp->met;
  mesh   = grp->mesh;
  ppt    = &mesh->point[ip];
  nsize  = met->size;

  /** Get barycentric coordinates **/
  ier = PMMG_baryCoord( oldGrp->mesh, pt, ppt->c, phi );

  /** Linear interpolation of the metrics */
  for( isize = 0; isize<nsize; isize++ ) {
    met->m[nsize*ip+isize] = 0.0;
    for( iloc=0; iloc<4; iloc++ ) {
      met->m[nsize*ip+isize] += phi[iloc]*oldMet->m[nsize*pt->v[iloc]+isize];
    }
  }

  return 1;
}

int PMMG_interpMetrics_grps( PMMG_pParMesh parmesh ) {
  PMMG_pGrp   grp,oldGrp;
  MMG5_pMesh  mesh,oldMesh;
  MMG5_pTetra pt;
  double      barycoord[4],vol,eps;
  int         *adja;
  int         igrp,ip,ie,ier;

  oldGrp  = &parmesh->old_listgrp[0];
  oldMesh = oldGrp->mesh;

  for( igrp=0; igrp<parmesh->ngrp; igrp++ ) {
    grp = &parmesh->listgrp[igrp];
    mesh = grp->mesh;

    if( !grp->met->m ) continue;

    ie = 1;
    for( ip=1; ip<mesh->np+1; ip++ ) {
      if( !MG_VOK(&mesh->point[ip]) ) continue;
      /** Locate point in the old mesh */
//      ie = PMMG_locatePoint( oldMesh, &mesh->point[ip], ie );
#warning Luca: Brute force method (only for testing), to remove when oldMesh->adja and localization from FMG will be ok.
      for( ie=1; ie<oldMesh->ne+1; ie++ ) {
        if( !MG_EOK(&oldMesh->tetra[ie]) ) continue;

        PMMG_baryCoord(oldMesh, &oldMesh->tetra[ie], mesh->point[ip].c, barycoord);
        vol = MMG5_orvol( oldMesh->point, oldMesh->tetra[ie].v );
        if( barycoord[0]>-MMG5_EPS*vol &&
            barycoord[1]>-MMG5_EPS*vol &&
            barycoord[2]>-MMG5_EPS*vol &&
            barycoord[3]>-MMG5_EPS*vol ) break;
      }
      if( ie==oldMesh->ne+1 ) {
        ie = 0;
      }
      if( !ie ) {
        fprintf(stderr,"\n  ## Error: %s: proc %d (grp %d), point %d not found, coords %e %e %e\n",__func__,parmesh->myrank,igrp,ip, mesh->point[ip].c[0],mesh->point[ip].c[1],mesh->point[ip].c[2]);
        return 0;
      }

      pt = &oldMesh->tetra[ie];

      /** Interpolate point metrics */
      ier = PMMG_interpMetrics_point(grp,oldGrp,pt,ip);
    }
  }

  return 1;
}

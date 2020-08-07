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
 * \file analys_pmmg.c
 * \brief Functions for parallel mesh analysis.
 * \author Luca Cirrottola (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Parallel mesh analysis.
 *
 */

#include "parmmg.h"

/**
 * \param parmesh pointer toward the parmesh structure
 * \param mesh pointer toward the mesh structure
 *
 * Check for singularities.
 * \remark Modeled after the MMG5_singul function.
 */
int PMMG_singul(PMMG_pParMesh parmesh,MMG5_pMesh mesh) {
  MMG5_pTria     pt;
  MMG5_pPoint    ppt,p1,p2;
  double         ux,uy,uz,vx,vy,vz,dd;
  int            list[MMG3D_LMAX+2],listref[MMG3D_LMAX+2],k,nc,xp,nr,ns,nre;
  char           i;

  nre = nc = 0;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( !MG_VOK(ppt) || ( ppt->tag & MG_CRN ) || ( ppt->tag & MG_NOM ) )
        continue;
      else if ( MG_EDG(ppt->tag) ) {
        /* Store the number of ridges passing through the point (xp) and the
         * number of ref edges (nr) */
#warning Luca: need to count in parallel and to merge the two half-ball lists, communicate {ng,nr,vx,vy,vz}
        ns = MMG5_bouler(mesh,mesh->adjt,k,i,list,listref,&xp,&nr,MMG3D_LMAX);

        if ( !ns )  continue;
        if ( (xp+nr) > 2 ) {
          ppt->tag |= MG_CRN + MG_REQ;
          ppt->tag &= ~MG_NOSURF;
          nre++;
          nc++;
        }
        else if ( (xp == 1) && (nr == 1) ) {
          ppt->tag |= MG_REQ;
          ppt->tag &= ~MG_NOSURF;
          nre++;
        }
        else if ( xp == 1 && !nr ){
          ppt->tag |= MG_CRN + MG_REQ;
          ppt->tag &= ~MG_NOSURF;
          nre++;
          nc++;
        }
        else if ( nr == 1 && !xp ){
          ppt->tag |= MG_CRN + MG_REQ;
          ppt->tag &= ~MG_NOSURF;
          nre++;
          nc++;
        }
        /* check ridge angle */
        else {
          p1 = &mesh->point[list[1]];
          p2 = &mesh->point[list[2]];
          ux = p1->c[0] - ppt->c[0];
          uy = p1->c[1] - ppt->c[1];
          uz = p1->c[2] - ppt->c[2];
          vx = p2->c[0] - ppt->c[0];
          vy = p2->c[1] - ppt->c[1];
          vz = p2->c[2] - ppt->c[2];
          dd = (ux*ux + uy*uy + uz*uz) * (vx*vx + vy*vy + vz*vz);
          if ( fabs(dd) > MMG5_EPSD ) {
            dd = (ux*vx + uy*vy + uz*vz) / sqrt(dd);
            if ( dd > -mesh->info.dhd ) {
              ppt->tag |= MG_CRN;
              nc++;
            }
          }
        }
      }
    }
  }

  if ( abs(mesh->info.imprim) > 3 && nre > 0 )
    fprintf(stdout,"     %d corners, %d singular points detected\n",nc,nre);
  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param mesh pointer toward the mesh structure
 * \return 0 if fail, 1 if success.
 *
 * Check all boundary triangles.
 */
int PMMG_analys_tria(PMMG_pParMesh parmesh,MMG5_pMesh mesh) {

  /**--- stage 1: data structures for surface */
  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"\n  ** SURFACE ANALYSIS\n");

  /* create tetra adjacency */
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Hashing problem (1). Exit program.\n");
    return 0;
  }

  /* create prism adjacency */
  if ( !MMG3D_hashPrism(mesh) ) {
    fprintf(stderr,"\n  ## Prism hashing problem. Exit program.\n");
    return 0;
  }
  /* compatibility triangle orientation w/r tetras */
  if ( !MMG5_bdryPerm(mesh) ) {
    fprintf(stderr,"\n  ## Boundary orientation problem. Exit program.\n");
    return 0;
  }

  /* identify surface mesh */
  if ( !MMG5_chkBdryTria(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    return 0;
  }
  MMG5_freeXTets(mesh);
  MMG5_freeXPrisms(mesh);

  /* set bdry entities to tetra: create xtetra and set references */
  if ( !MMG5_bdrySet(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    return 0;
  }

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param mesh pointer toward the mesh structure
 *
 * \remark Modeled after the MMG3D_analys function, it doesn't deallocate the
 * tria structure in order to be able to build communicators.
 */
int PMMG_analys(PMMG_pParMesh parmesh,MMG5_pMesh mesh) {
  MMG5_Hash      hash;
  MMG5_HGeom     hpar;

  /* Set surface triangles to required in nosurf mode or for parallel boundaries */
  MMG3D_set_reqBoundaries(mesh);


  /* create surface adjacency */
  if ( !MMG3D_hashTria(mesh,&hash) ) {
    MMG5_DEL_MEM(mesh,hash.item);
    fprintf(stderr,"\n  ## Hashing problem (2). Exit program.\n");
    return 0;
  }

  /* build hash table for geometric edges */
  if ( !MMG5_hGeom(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem (0). Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    MMG5_DEL_MEM(mesh,mesh->htab.geom);
    return 0;
  }

  /**--- stage 2: surface analysis */
  if ( abs(mesh->info.imprim) > 5  || mesh->info.ddebug )
    fprintf(stdout,"  ** SETTING TOPOLOGY\n");

  /* identify connexity */
  if ( !MMG5_setadj(mesh) ) {
    fprintf(stderr,"\n  ## Topology problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /* check for ridges: check dihedral angle using adjacent triangle normals */
#warning Luca: here a parallel edge communicator would be useful
  if ( mesh->info.dhd > MMG5_ANGLIM && !MMG5_setdhd(mesh) ) {
    fprintf(stderr,"\n  ## Geometry problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /* identify singularities */
  if ( !PMMG_singul(parmesh,mesh) ) {
    fprintf(stderr,"\n  ## PMMG_singul problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug )
    fprintf(stdout,"  ** DEFINING GEOMETRY\n");

#warning Luca: check why this is not done before
  /* Pass MG_PARBDY tag from faces to edges */
  {
    assert( parmesh->ngrp == 1);
    PMMG_pGrp    grp = &parmesh->listgrp[0];
    MMG5_pTetra  pt;
    MMG5_pxTetra pxt;
    MMG5_pPoint  ppt;
    int i,iel,ifac,j,ia;
    for ( i=0; i<grp->nitem_int_face_comm; i++ ) {
      iel  =   grp->face2int_face_comm_index1[i] / 12;
      ifac = ( grp->face2int_face_comm_index1[i] % 12 ) / 3;
      pt = &mesh->tetra[iel];
      assert( pt->xt );
      pxt = &mesh->xtetra[pt->xt];
//      /* If already boundary, make it recognizable as a "true" boundary */
//      if( pxt->ftag[ifac] & MG_BDY ) pxt->ftag[ifac] |= MG_PARBDYBDY;
//      /* Tag face */
//      PMMG_tag_par_face(pxt,ifac);
      /* Tag face edges */
      for ( j=0; j<3; j++ ) {
        ia = MMG5_iarf[ifac][j];
        pxt->tag[j] |= MG_PARBDY;
//        PMMG_tag_par_edge(pxt,ia);
      }
//      /* Tag face nodes */
//      for ( j=0 ; j<3 ; j++) {
//        ppt = &mesh->point[pt->v[MMG5_idir[ifac][j]]];
//        PMMG_tag_par_node(ppt);
//      }
    }
  }

  /* Hash parallel edges */
  if( PMMG_hashPar( mesh,&hpar ) != PMMG_SUCCESS ) return 0;

  /* Build edge communicator */
  if( !PMMG_build_edgeComm( parmesh,mesh,&hpar ) ) return 0;

  /* define (and regularize) normals: create xpoints */
#warning Luca: it uses boulen (twice) and boulec
  if ( !MMG5_norver(mesh) ) {
    fprintf(stderr,"\n  ## Normal problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /* set non-manifold edges sharing non-intersecting multidomains as required */
  if ( abs(mesh->info.imprim) > 5  || mesh->info.ddebug )
    fprintf(stdout,"  ** UPDATING TOPOLOGY AT NON-MANIFOLD POINTS\n");

  if ( !MMG5_setNmTag(mesh,&hash) ) {
    fprintf(stderr,"\n  ## Non-manifold topology problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    return 0;
  }

  /* check subdomains connected by a vertex and mark these vertex as corner and required */
#warning Luca: check that parbdy are skipped
  MMG5_chkVertexConnectedDomains(mesh);

  /* build hash table for geometric edges */
  if ( !mesh->na && !MMG5_hGeom(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem (0). Exit program.\n");
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    MMG5_DEL_MEM(mesh,mesh->htab.geom);
    return 0;
  }

  /* Update edges tags and references for xtetras */
  if ( !MMG5_bdryUpdate(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    return 0;
  }

  /* define geometry for non manifold points */
  if ( !MMG3D_nmgeom(mesh) ) return 0;

  /* release memory */
  MMG5_DEL_MEM(mesh,hpar.geom);
  MMG5_DEL_MEM(mesh,mesh->htab.geom);
  MMG5_DEL_MEM(mesh,mesh->adjt);

  if ( mesh->nprism ) MMG5_DEL_MEM(mesh,mesh->adjapr);

  return 1;
}

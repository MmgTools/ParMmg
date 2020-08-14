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
 * \param parmesh pointer toward the parmesh structure.
 * \param mesh pointer toward the mesh structure.
 * \param adjt pointer toward the table of triangle adjacency.
 * \param start index of triangle where we start to work.
 * \param ip index of vertex on which we work.
 * \param list pointer toward the computed list of GEO vertices incident to \a ip.
 * \param listref pointer toward the corresponding edge references
 * \param ng pointer toward the number of ridges.
 * \param nr pointer toward the number of reference edges.
 * \param lmax maxmum size for the ball of the point \a ip.
 * \return The number of edges incident to the vertex \a ip.
 *
 * Store edges and count the number of ridges and reference edges incident to
 * the vertex \a ip.
 * \remark Sams as MMG5_bouler(), but skip edges whose extremity is flagged with
 * a rank lower than myrank.
 *
 */
int PMMG_bouler(PMMG_pParMesh parmesh,MMG5_pMesh mesh,int *adjt,int start,int ip,
                 int *list,int *listref,int *ng,int *nr,int lmax) {
  MMG5_pTria    pt;
  int           *adja,k,ns;
  char          i,i1,i2;

  pt  = &mesh->tria[start];
  if ( !MG_EOK(pt) )  return 0;

  /* check other triangle vertices */
  k  = start;
  i  = ip;
  *ng = *nr = ns = 0;
  do {
    i1 = MMG5_inxt2[i];
    if ( MG_EDG(pt->tag[i1]) && (parmesh->myrank < mesh->point[pt->v[i2]].flag) ) {
      i2 = MMG5_iprv2[i];
      if ( pt->tag[i1] & MG_GEO )
        *ng = *ng + 1;
      else if ( pt->tag[i1] & MG_REF )
        *nr = *nr + 1;
      ns++;
      list[ns] = pt->v[i2];
      listref[ns] = pt->edg[i1];
      if ( ns > lmax-2 )  return -ns;
    }
    adja = &adjt[3*(k-1)+1];
    k  = adja[i1] / 3;
    i  = adja[i1] % 3;
    i  = MMG5_inxt2[i];
    pt = &mesh->tria[k];
  }
  while ( k && k != start );

  /* reverse loop */
  if ( k != start ) {
    k = start;
    i = ip;
    do {
      pt = &mesh->tria[k];
      i2 = MMG5_iprv2[i];
      if ( MG_EDG(pt->tag[i2]) && (parmesh->myrank < mesh->point[pt->v[i1]].flag) ) {
        i1 = MMG5_inxt2[i];
        if ( pt->tag[i2] & MG_GEO )
          *ng = *ng + 1;
        else if ( pt->tag[i2] & MG_REF )
          *nr = *nr + 1;
        ns++;
        list[ns] = pt->v[i1];
        listref[ns] = pt->edg[i2];
        if ( ns > lmax-2 )  return -ns;
      }
      adja = &adjt[3*(k-1)+1];
      k = adja[i2] / 3;
      i = adja[i2] % 3;
      i = MMG5_iprv2[i];
    }
    while ( k && k != start );
  }

  return ns;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param mesh pointer toward the mesh structure
 *
 * Check for singularities.
 * \remark Modeled after the MMG5_singul function.
 */
int PMMG_singul(PMMG_pParMesh parmesh,MMG5_pMesh mesh) {
  PMMG_pGrp      grp;
  PMMG_pInt_comm int_node_comm;
  PMMG_pExt_comm ext_node_comm;
  MPI_Comm       comm;
  MPI_Status     status;
  MMG5_pTria     pt;
  MMG5_pPoint    ppt,p1,p2;
  double         ux,uy,uz,vx,vy,vz,dd;
  int            list[MMG3D_LMAX+2],listref[MMG3D_LMAX+2],k,nc,xp,nr,ns,nre;
  int            ip,idx,iproc;
  int            nitem,color;
  int            *intvalues,*itosend,*itorecv,*iproc2comm;
  double         *doublevalues,*rtosend,*rtorecv;
  char           i,j;

  comm   = parmesh->comm;
  assert( parmesh->ngrp == 1 );
  grp = &parmesh->listgrp[0];
  int_node_comm = parmesh->int_node_comm;

  /* Allocate intvalues to accomodate xp and nr for each point, doublevalues to
   * accomodate two edge vectors at most. */
  PMMG_CALLOC(parmesh,int_node_comm->intvalues,2*int_node_comm->nitem,int,"intvalues",return 0);
  PMMG_CALLOC(parmesh,int_node_comm->doublevalues,6*int_node_comm->nitem,double,"doublevalues",return 0);
  intvalues    = int_node_comm->intvalues;
  doublevalues = int_node_comm->doublevalues;

  /* Store source triangle for every boundary point */
  for( ip = 1; ip <= mesh->np; ip++ ) {
    ppt = &mesh->point[ip];
    if( !MG_VOK(ppt) || !(ppt->tag & MG_PARBDY ) ) continue;
    ppt->flag = 0;
  }
  for( k = 1; k <= mesh->nt; k++ ) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;
    for( i = 0; i < 3; i++ ) {
      ppt = &mesh->point[pt->v[i]];
      if( ppt->flag ) continue;
      ppt->s = 3*k+i;
      ppt->flag++;
    }
  }

  /* Array to reorder communicators */
  PMMG_MALLOC(parmesh,iproc2comm,parmesh->nprocs,int,"iproc2comm",return 0);

  for( iproc = 0; iproc < parmesh->nprocs; iproc++ )
    iproc2comm[iproc] = PMMG_UNSET;

  for( k = 0; k < parmesh->next_node_comm; k++ ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    iproc = ext_node_comm->color_out;
    iproc2comm[iproc] = k;
  }
  for( iproc = parmesh->nprocs-1; iproc >= 0; iproc-- ) {
    k = iproc2comm[iproc];
    if( k == PMMG_UNSET ) continue;
    ext_node_comm = &parmesh->ext_node_comm[k];
    for( i = 0; i < ext_node_comm->nitem; i++ ) {
      idx = ext_node_comm->int_comm_index[i];
      intvalues[idx] = iproc;
    }
  }
  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    ip  = grp->node2int_node_comm_index1[i];
    idx = grp->node2int_node_comm_index2[i];
    ppt = &mesh->point[ip];
    if( !MG_VOK(ppt) ) continue;
    ppt->flag = intvalues[idx];
  }

  /* Local singularity analysis */
  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    ip  = grp->node2int_node_comm_index1[i];
    idx = grp->node2int_node_comm_index2[i];
    ppt = &mesh->point[ip];
    if ( !MG_VOK(ppt) || ( ppt->tag & MG_CRN ) || ( ppt->tag & MG_NOM ) )
      continue;
    else if ( MG_EDG(ppt->tag) ) {
      /* Store the number of ridges passing through the point (xp) and the
       * number of ref edges (nr) */
#warning Luca: need to count in parallel and to merge the two half-ball lists, communicate {ng,nr,vx,vy,vz}
#warning Luca: need to count edges on communicators only once
      ns = PMMG_bouler(parmesh,mesh,mesh->adjt,ppt->s/3,ppt->s%3,list,listref,&xp,&nr,MMG3D_LMAX);
      /* Add nb of geo/ridges to intvalue */
      intvalues[2*idx]   = xp;
      intvalues[2*idx+1] = nr;
      if ( (xp+nr) > 2 ) continue;
      /* Add edge vectors to doublevalues */
      for( j = 0; j < xp+nr; j++ ) {
        p1 = &mesh->point[list[j]];
        doublevalues[6*idx+3*j]   = p1->c[0]-ppt->c[0];
        doublevalues[6*idx+3*j+1] = p1->c[1]-ppt->c[1];
        doublevalues[6*idx+3*j+2] = p1->c[2]-ppt->c[2];
      }
    }
  }

  /* Exchange values on the interfaces among procs */
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    nitem         = ext_node_comm->nitem;
    color         = ext_node_comm->color_out;

    PMMG_CALLOC(parmesh,ext_node_comm->itosend,2*nitem,int,"itosend array",
                return 0);
    PMMG_CALLOC(parmesh,ext_node_comm->itorecv,2*nitem,int,"itorecv array",
                return 0);
    PMMG_CALLOC(parmesh,ext_node_comm->rtosend,6*nitem,double,"rtosend array",
                return 0);
    PMMG_CALLOC(parmesh,ext_node_comm->rtorecv,6*nitem,double,"rtorecv array",
                return 0);
    itosend = ext_node_comm->itosend;
    itorecv = ext_node_comm->itorecv;
    rtosend = ext_node_comm->rtosend;
    rtorecv = ext_node_comm->rtorecv;

    for ( i=0; i<nitem; ++i ) {
      idx  = ext_node_comm->int_comm_index[i];
      for( j = 0; j < 2; j++ ) itosend[2*i+j] = intvalues[2*idx+j];
      for( j = 0; j < 6; j++ ) rtosend[6*i+j] = doublevalues[6*idx+j];
    }

#warning Luca: change this tags
    MPI_CHECK(
      MPI_Sendrecv(itosend,nitem,MPI_INT,color,MPI_PARMESHGRPS2PARMETIS_TAG+1,
                   itorecv,nitem,MPI_INT,color,MPI_PARMESHGRPS2PARMETIS_TAG+1,
                   comm,&status),return 0 );

    MPI_CHECK(
      MPI_Sendrecv(rtosend,nitem,MPI_DOUBLE,color,MPI_PARMESHGRPS2PARMETIS_TAG+2,
                   rtorecv,nitem,MPI_DOUBLE,color,MPI_PARMESHGRPS2PARMETIS_TAG+2,
                   comm,&status),return 0 );
  }

  /* First pass: Sum nb. of singularities */
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    itorecv = ext_node_comm->itorecv;

    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx  = ext_node_comm->int_comm_index[i];

      /* Current nb. of singularities == current nb. of stored edge vectors */
      ns = intvalues[2*idx]+intvalues[2*idx+1];

      /* Increment nb. of singularities */
      intvalues[2*idx]   += itorecv[2*i];
      intvalues[2*idx+1] += itorecv[2*i+1];

      /* Store new edge vectors if not already full */
      if( (ns < 2) && (intvalues[2*idx]+intvalues[2*idx+1] <= 2) ) {
        for( j = 0; j < itorecv[2*i]+itorecv[2*i+1]; j++ ) {
          doublevalues[6*idx+3*(ns+j)]   = rtorecv[6*i+3*j];
          doublevalues[6*idx+3*(ns+j)+1] = rtorecv[6*i+3*j+1];
          doublevalues[6*idx+3*(ns+j)+2] = rtorecv[6*i+3*j+2];
        }
      }
    }
  }

  /* Second pass: Analysis */
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    itosend = ext_node_comm->itosend;
    itorecv = ext_node_comm->itorecv;
    rtosend = ext_node_comm->rtosend;
    rtorecv = ext_node_comm->rtorecv;

    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx  = ext_node_comm->int_comm_index[i];

      xp = intvalues[2*idx];
      nr = intvalues[2*idx+1];
      if ( (xp == PMMG_UNSET) || (nr == PMMG_UNSET) )  continue;

      if ( !(xp+nr) )  continue;
      if ( (xp+nr) > 2 ) {
        intvalues[2*idx]   = PMMG_UNSET;
        intvalues[2*idx+1] = PMMG_UNSET;
        nre++;
        nc++;
      }
      else if ( (xp == 1) && (nr == 1) ) {
        intvalues[2*idx]   = PMMG_UNSET;
        nre++;
      }
      else if ( xp == 1 && !nr ){
        intvalues[2*idx]   = PMMG_UNSET;
        intvalues[2*idx+1] = PMMG_UNSET;
        nre++;
        nc++;
      }
      else if ( nr == 1 && !xp ){
        intvalues[2*idx]   = PMMG_UNSET;
        intvalues[2*idx+1] = PMMG_UNSET;
        nre++;
        nc++;
      }
      /* check ridge angle */
      else {
        ux = doublevalues[6*idx];
        uy = doublevalues[6*idx+1];
        uz = doublevalues[6*idx+2];
        vx = doublevalues[6*idx+3];
        vy = doublevalues[6*idx+4];
        vz = doublevalues[6*idx+5];
        dd = (ux*ux + uy*uy + uz*uz) * (vx*vx + vy*vy + vz*vz);
        if ( fabs(dd) > MMG5_EPSD ) {
          dd = (ux*vx + uy*vy + uz*vz) / sqrt(dd);
          if ( dd > -mesh->info.dhd ) {
            intvalues[2*idx+1] = PMMG_UNSET;
            nc++;
          }
        }
      }
    }
  }

  /* Third pass: Tag points */
  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    ip  = grp->node2int_node_comm_index1[i];
    idx = grp->node2int_node_comm_index2[i];
    ppt = &mesh->point[ip];
    if( intvalues[2*idx]   == PMMG_UNSET ) { /* is required */
      ppt->tag |= MG_REQ;
      ppt->tag &= ~MG_NOSURF;
    }
    if( intvalues[2*idx+1] == PMMG_UNSET ) { /* is corner */
      ppt->tag |= MG_CRN;
    }
  }

  if ( abs(mesh->info.imprim) > 3 && nre > 0 )
    fprintf(stdout,"     %d corners, %d singular points detected\n",nc,nre);


  PMMG_DEL_MEM(parmesh,iproc2comm,int,"iproc2comm");
  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,int,"intvalues");
  PMMG_DEL_MEM(parmesh,int_node_comm->doublevalues,double,"doublevalues");
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    PMMG_DEL_MEM(parmesh,ext_node_comm->itosend,int,"itosend array");
    PMMG_DEL_MEM(parmesh,ext_node_comm->itorecv,int,"itorecv array");
    PMMG_DEL_MEM(parmesh,ext_node_comm->rtosend,double,"rtosend array");
    PMMG_DEL_MEM(parmesh,ext_node_comm->rtorecv,double,"rtorecv array");
  }

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

  /* Hash parallel edges */
  if( PMMG_hashPar( mesh,&hpar ) != PMMG_SUCCESS ) return 0;

  /* Build edge communicator */
  if( !PMMG_build_edgeComm( parmesh,mesh,&hpar ) ) return 0;

  /* check for ridges: check dihedral angle using adjacent triangle normals */
#warning Luca: here a parallel edge communicator would be useful
  if ( mesh->info.dhd > MMG5_ANGLIM && !MMG5_setdhd(mesh) ) {
    fprintf(stderr,"\n  ## Geometry problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /* identify singularities on interior points */
  if ( !MMG5_singul(mesh) ) {
    fprintf(stderr,"\n  ## MMG5_singul problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /* identify singularities on parallel points */
  if ( !PMMG_singul(parmesh,mesh) ) {
    fprintf(stderr,"\n  ## PMMG_singul problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug )
    fprintf(stdout,"  ** DEFINING GEOMETRY\n");

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
  PMMG_edge_comm_free( parmesh );
  MMG5_DEL_MEM(mesh,hpar.geom);
  MMG5_DEL_MEM(mesh,mesh->htab.geom);
  MMG5_DEL_MEM(mesh,mesh->adjt);

  if ( mesh->nprism ) MMG5_DEL_MEM(mesh,mesh->adjapr);

  return 1;
}

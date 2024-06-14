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
 * \file isovalue_pmmg.c
 * \brief Create implicit surface in distribuited mesh.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (InriaSoft)
 * \author Laetitia Mottet (UBordeaux)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Main library functions (parallel remeshing starting from centralized or
 * distributed data.
 *
 */

#include "parmmg.h"
#include "mmgexterns_private.h"

/**
 * \param parmesh pointer toward a parmesh structure
 * \param mesh pointer toward the mesh structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * \todo Fill the funtion
 *
 * Snap values of the level set function very close to 0 to exactly 0,
 * and prevent nonmanifold patterns from being generated.
 *
 */
int PMMG_create_overlap(PMMG_pParMesh parmesh, MPI_Comm comm) {

  PMMG_pInt_comm int_node_comm;
  PMMG_pExt_comm ext_node_comm;
  PMMG_pGrp      grp,grp_init,grp_overlap;
  MMG5_pPoint    p0;
  MMG5_pMesh     mesh,meshOld;
  MMG5_pTetra    pt;
  MPI_Status     status;
  MPI_Request    *request;
  PMMG_pOverlap   overlap;

  int ier = 1; /* initialize error */

  MMG5_int ip,ngrp,ne_init,n2inc_max,f2ifc_max,base_init;
  MMG5_int nt_overlap, np_overlap;

  int i, k, ireq;

  int inode;
  int i_commn;

  int nitem_int_node;
  int nitem_ext_node;
  int next_node_comm;

  int *ip_overlap[2];
  double c_inode[3];
  MMG5_int np_new;

  int      *tetraVertices_ToSend, *tetraVertices_ToRecv;
  int      *tetraIdx_ToSend, *tetraIdx_ToRecv;
  double   *pointCoord_ToSend, *pointCoord_ToRecv;
  int      *pointIdx_ToSend, *pointIdx_ToRecv;
  uint16_t *pointTag_ToSend, *pointTag_ToRecv;

  int color_in_node,color_out_node;
  int nitem_to_share_old;

  int idx_node_ext,idx_node_int,idx_node_mesh;

  double *rtosend,*rtorecv,*doublevalues;
  int    *itosend,*itorecv,*intvalues;

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
      fprintf(stdout,"\n      ## TODO:: PMMG_create_overlap.\n");

  /* For now - creation of overlap works only on packed mesh ,
      i.e. on one group */
  /* Ensure only one group on each proc */
  assert(parmesh->ngrp == 1);

  /* Initialization */
  grp  = &parmesh->listgrp[0];
  mesh = parmesh->listgrp[0].mesh;
  int_node_comm = parmesh->int_node_comm;
  overlap = &parmesh->overlap;
  ireq = 0;
  request     = NULL;
  nitem_to_share_old = 1;

  /* Allocate and reset memory */
  /* allocate doublevalues */
  PMMG_CALLOC(parmesh,int_node_comm->doublevalues,3*10*int_node_comm->nitem,double,"doublevalues",ier = 0);
  memset(int_node_comm->doublevalues,0x00,3*10*int_node_comm->nitem*sizeof(double));

  PMMG_CALLOC(parmesh,ip_overlap[0],1,int,"ip_overlap",ier = 0);
  PMMG_CALLOC(parmesh,ip_overlap[1],1,int,"ip_overlap",ier = 0);
  memset(ip_overlap[0],0x00,1*sizeof(int));
  memset(ip_overlap[1],0x00,1*sizeof(int));

  PMMG_CALLOC(parmesh,tetraVertices_ToSend,4*mesh->ne,int,"tetraVertices_ToSend",ier = 0);
  memset(tetraVertices_ToSend,0x00,4*mesh->ne*sizeof(int));

  PMMG_CALLOC(parmesh,tetraIdx_ToSend,mesh->ne,int,"tetraIdx_ToSend",ier = 0);
  memset(tetraIdx_ToSend,0x00,mesh->ne*sizeof(int));

  doublevalues = int_node_comm->doublevalues;
  intvalues    = int_node_comm->intvalues;

  next_node_comm = parmesh->next_node_comm;  // Number of communicator for nodes

  PMMG_CALLOC(parmesh,overlap->np_overlap,next_node_comm,int,"parmesh->overlap.np_overlap ",return 0);
  PMMG_CALLOC(parmesh,overlap->nt_overlap,next_node_comm,int,"parmesh->overlap.nt_overlap ",return 0);

  /* STEP 1 - Identify nodes and tetra to communicate and MPI_Sendrecv them */

  /* Loop over the number of node communicator */
  for (i_commn=0; i_commn<next_node_comm; i_commn++) {

    /* Get external edge communicator information */
    ext_node_comm  = &parmesh->ext_node_comm[i_commn]; // External node communicator
    color_in_node  = ext_node_comm->color_in;          // Color of the hosting proc - this proc
    color_out_node = ext_node_comm->color_out;         // Color of the remote  proc - the proc to exchange with
    nitem_ext_node = ext_node_comm->nitem;             // Number of nodes in common between these 2 procs

    rtosend = ext_node_comm->rtosend;
    rtorecv = ext_node_comm->rtorecv;

    itosend = ext_node_comm->itosend;
    itorecv = ext_node_comm->itorecv;

    /* Allocate memory */
    /* rtosend   - at max the whole mesh with np points of 3 coordinates is sent */
    /* intvalues - at max the whole mesh with np points i sent*/
    /* itosend   - the number of point to send    np_overlap_send - 1 int */
    /* itorecv   - the number of point to receive np_overlap_recv - 1 int */
    PMMG_CALLOC(parmesh,rtosend,3*mesh->np,double,"rtosend",ier = 0);
    // PMMG_CALLOC(parmesh,intvalues,mesh->np,int,"intvalues",ier = 0);
    PMMG_CALLOC(parmesh,itosend,1,int,"itosend",ier = 0);
    PMMG_CALLOC(parmesh,itorecv,1,int,"itorecv",ier = 0);
    PMMG_CALLOC(parmesh,pointTag_ToSend,mesh->np,uint16_t,"pointTag_ToSend",ier = 0);

    PMMG_CALLOC(parmesh,pointCoord_ToSend,3*mesh->np,double,"pointCoord_ToSend",ier = 0);
    PMMG_CALLOC(parmesh,pointIdx_ToSend,mesh->np,int,"pointIdx_ToSend",ier = 0);


    /* Initialize arrays to send */
    memset(rtosend,0x00,3*mesh->np*sizeof(double));
    // memset(intvalues,0x00,mesh->np*sizeof(int));
    memset(itosend,0x00,sizeof(int));
    memset(itorecv,0x00,sizeof(int));
    // memset(pointTag_ToSend,0x00,mesh->np*sizeof(uint16_t));
    // memset(pointIdx_ToSend,0x00,mesh->np*sizeof(int));

    /* STEP 1.1 - First loop to estimate the number of nodes and tetra to communicate
                  for the CALLOC of itosend, itorecv, rtosend and rtorecv */
    nt_overlap = np_overlap = 0;

    /* Loop over the nodes in the external edge communicator */
    for (inode=0; inode < nitem_ext_node; inode++) {

      /* Get the indices of the nodes in internal communicators */
      idx_node_ext  = ext_node_comm->int_comm_index[inode];
      idx_node_int  = grp->node2int_node_comm_index2[idx_node_ext];
      idx_node_mesh = grp->node2int_node_comm_index1[idx_node_int];

      /* Add the flag 1 to these nodes*/
      p0 = &mesh->point[idx_node_mesh];
      p0->flag=-1;
    }

    /* Loop over number of tetra and assign MG_OVERLAP */
    for (k=1; k<=mesh->ne; k++) {
      /* Get the tetra k */
      pt  = &mesh->tetra[k];
      if (pt->tag & MG_OVERLAP) {
        pt->tag =~ MG_OVERLAP;
      }

      /* Loop over vertex of tetra and assign MG_OVERLAP to tetras */
      for (i=0; i<4; i++) {
        ip = pt->v[i];
        p0 = &mesh->point[ip];
        if ( p0->flag < 0 ) {
          pt->tag |= MG_OVERLAP;
          nt_overlap++;
          break;
        }
      }

      /* Store info of tetra to send */
      tetraVertices_ToSend[4*nt_overlap]  = pt->v[0];
      tetraVertices_ToSend[4*nt_overlap+1]= pt->v[1];
      tetraVertices_ToSend[4*nt_overlap+2]= pt->v[2];
      tetraVertices_ToSend[4*nt_overlap+3]= pt->v[3];
      tetraIdx_ToSend[nt_overlap] = k;

      /* If tetra is now MG_OVERLAP, then assign MG_OVERLAP to nodes */
      if (pt->tag & MG_OVERLAP) {
        for (i=0; i<4; i++) {
          ip = pt->v[i];
          p0 = &mesh->point[ip];
          if ( (p0->flag>=0)  && (p0->flag!=color_out_node+1) ) {
            p0->flag = color_out_node+1;
            pointCoord_ToSend[3*np_overlap]   = p0->c[0];
            pointCoord_ToSend[3*np_overlap+1] = p0->c[1];
            pointCoord_ToSend[3*np_overlap+2] = p0->c[2];
            pointIdx_ToSend[np_overlap] = ip;
            pointTag_ToSend[np_overlap] = p0->tag+MG_OVERLAP;
            np_overlap++;
          }
        }
      }
      pt->tag =~ MG_OVERLAP;
    }

    /* Reinitialise flags of points */
    for (inode=0; inode < nitem_ext_node; inode++) {
      /* Get the indices of the nodes in internal communicators */
      idx_node_ext  = ext_node_comm->int_comm_index[inode];
      idx_node_int  = grp->node2int_node_comm_index2[idx_node_ext];
      idx_node_mesh = grp->node2int_node_comm_index1[idx_node_int];

      /* Add the tag MG_OVERLAP to these nodes*/
      p0 = &mesh->point[idx_node_mesh];
      p0->flag = 0;
    }

    /* Number of item to share */
    itosend[0] = np_overlap;

    /* First send and receive the number of points to exchange/share */
    /* TODO :: Use MPI_Probe instead? */
    MPI_CHECK(
      MPI_Sendrecv(itosend,1,MPI_INT,color_out_node,MPI_OVERLAP_TAG,
                   itorecv,1,MPI_INT,color_out_node,MPI_OVERLAP_TAG,
                   comm,&status),return 0 );

    ext_node_comm->nitem_to_share = itorecv[0];

    // fprintf(stdout, "OVERLAP Proc=%d->Proc=%d :: npoint=%d, ntetra=%d; nptosend=%d, nptorecv=%d \n",color_in_node,color_out_node,np_overlap,nt_overlap,itosend[0],itorecv[0]);

    /* Send and receive the index of points to exchange/share */
    PMMG_CALLOC(parmesh,pointIdx_ToRecv,ext_node_comm->nitem_to_share,int,"pointIdx_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointIdx_ToSend,np_overlap,MPI_INT,color_out_node,MPI_OVERLAP_TAG+1,
                   pointIdx_ToRecv,ext_node_comm->nitem_to_share,MPI_INT,color_out_node,MPI_OVERLAP_TAG+1,
                   comm,&status),return 0 );

    /* Send and receive the tags of points to exchange/share */
    PMMG_CALLOC(parmesh,pointTag_ToRecv,ext_node_comm->nitem_to_share,uint16_t,"pointTag_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointTag_ToSend,np_overlap,MPI_UINT16_T,color_out_node,MPI_OVERLAP_TAG+2,
                   pointTag_ToRecv,ext_node_comm->nitem_to_share,MPI_UINT16_T,color_out_node,MPI_OVERLAP_TAG+2,
                   comm,&status),return 0 );

    /* For rtorecv, we receive nitem_to_share number of points */
    PMMG_CALLOC(parmesh,pointCoord_ToRecv,3*ext_node_comm->nitem_to_share,double,"pointCoord_ToRecv",ier = 0);
    memset(pointCoord_ToRecv,0x00,3*ext_node_comm->nitem_to_share*sizeof(double));

    /* Now send and receive the coordinates of the points */
    MPI_CHECK(
      MPI_Sendrecv(pointCoord_ToSend,3*np_overlap,MPI_DOUBLE,color_out_node,MPI_OVERLAP_TAG+3,
                   pointCoord_ToRecv,3*ext_node_comm->nitem_to_share,MPI_DOUBLE,color_out_node,MPI_OVERLAP_TAG+3,
                   comm,&status),return 0 );

    /* Add these point to the local mesh->point */
    PMMG_REALLOC(parmesh,ip_overlap[0],ext_node_comm->nitem_to_share+nitem_to_share_old,nitem_to_share_old,
                 int,"Re-allocate ip_overlap",return 0);
    PMMG_REALLOC(parmesh,ip_overlap[1],ext_node_comm->nitem_to_share+nitem_to_share_old,nitem_to_share_old,
                 int,"Re-allocate ip_overlap",return 0);

    uint16_t tag_inode;

    mesh->xpmax  = MG_MAX( (long long)(1.5*mesh->xp),mesh->npmax);

    for (inode=0; inode < ext_node_comm->nitem_to_share; inode++) {
      c_inode[0]=pointCoord_ToRecv[3*inode];
      c_inode[1]=pointCoord_ToRecv[3*inode+1];
      c_inode[2]=pointCoord_ToRecv[3*inode+2];
      tag_inode =pointTag_ToRecv[inode];
      if (tag_inode & MG_PARBDY) {
        fprintf(stdout, "OVERLAP Proc=%d->Proc=%d :: inode=%d, tag_inode=%d, ip=%d \n",parmesh->myrank,color_out_node,inode,tag_inode,pointIdx_ToRecv[inode]);
      }
      np_new = MMG3D_newPt(mesh,c_inode,tag_inode,0);

      // For debugging - this array is not really needed for now...
      ip_overlap[0][inode+nitem_to_share_old-1] = pointIdx_ToRecv[inode];
      ip_overlap[1][inode+nitem_to_share_old-1] = np_new;
    }
    nitem_to_share_old += ext_node_comm->nitem_to_share;

    /* Deallocate memory*/
    PMMG_DEL_MEM(parmesh,itosend,int,"itosend");
    PMMG_DEL_MEM(parmesh,itorecv,int,"itorecv");

    PMMG_DEL_MEM(parmesh,pointCoord_ToRecv,double,"pointCoord_ToRecv");
    PMMG_DEL_MEM(parmesh,pointCoord_ToSend,double,"pointCoord_ToSend");
    PMMG_DEL_MEM(parmesh,pointIdx_ToRecv,double,"pointIdx_ToRecv");
    PMMG_DEL_MEM(parmesh,pointIdx_ToSend,double,"pointIdx_ToSend");


  }

  return 1;
}
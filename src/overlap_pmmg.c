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

  int np_in2out, np_out2in;
  int nt_in2out, nt_out2in;
  int np_out,np_in;

  int nitem_int_node;
  int nitem_ext_node;
  int next_node_comm;

  int LOCAL_MPI_TAG;

  int *ip_overlap[2];
  double c_inode[3];
  MMG5_int np_new;

  int      *pointIdxInterface_ToSend, *pointIdxInterface_ToRecv;
  int      *tetraVertices_ToSend, *tetraVertices_ToRecv;
  int      *tetraVertices_ToRecv_outIdx, *tetraVertices_ToRecv_inIdx;
  int      *tetraFlag_ToSend, *tetraFlag_ToRecv;
  int      *hash_overlap_ToSend, *hash_overlap_ToRecv;
  int      *tetraIdx_ToSend, *tetraIdx_ToRecv;
  double   *pointCoord_ToSend, *pointCoord_ToRecv;
  int      *pointIdx_ToSend, *pointIdx_ToRecv;
  uint16_t *pointTag_ToSend, *pointTag_ToRecv;
  int *hash_overlap_new2old, *hash_overlap_old2new;
  int *hash_in2out, *hash_out2in;

  int *tetraVertices_ToAdd;
  int ip1,ip2,icoord;
  int ip_in, ip_out;
  uint16_t tag_inode;

  int color_in_node,color_out_node;

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
  ireq = 0;
  request     = NULL;
  next_node_comm = parmesh->next_node_comm;  // Number of communicator for nodes
  LOCAL_MPI_TAG = 1;

  /* Allocate and reset memory */
  PMMG_CALLOC(parmesh,parmesh->overlap,next_node_comm,PMMG_Overlap,"allocate PMMG_Overlap ",ier = 0);

  // intvalues    = int_node_comm->intvalues;

  /* STEP 1 - Identify nodes and tetra to communicate and MPI_Sendrecv them */

  /* Loop over the number of node communicator */
  for (i_commn=0; i_commn<next_node_comm; i_commn++) {

    /* Get external edge communicator information */
    ext_node_comm  = &parmesh->ext_node_comm[i_commn]; // External node communicator
    color_in_node  = ext_node_comm->color_in;          // Color of the hosting proc - this proc
    color_out_node = ext_node_comm->color_out;         // Color of the remote  proc - the proc to exchange with
    nitem_ext_node = ext_node_comm->nitem;             // Number of nodes in common between these 2 procs

    itosend = ext_node_comm->itosend;
    itorecv = ext_node_comm->itorecv;

    /* Fill overlap values */
    overlap = &parmesh->overlap[i_commn];
    overlap->color_in  = color_in_node;
    overlap->color_out = color_out_node;

    /* Allocate memory */
    PMMG_CALLOC(parmesh,tetraVertices_ToSend,4*mesh->ne,int,"tetraVertices_ToSend",ier = 0);
    PMMG_CALLOC(parmesh,tetraFlag_ToSend,4*mesh->ne,int,"tetraFlag_ToSend",ier = 0);
    PMMG_CALLOC(parmesh,tetraIdx_ToSend,mesh->ne,int,"tetraIdx_ToSend",ier = 0);
    PMMG_CALLOC(parmesh,pointTag_ToSend,mesh->np,uint16_t,"pointTag_ToSend",ier = 0);
    PMMG_CALLOC(parmesh,pointCoord_ToSend,3*mesh->np,double,"pointCoord_ToSend",ier = 0);
    PMMG_CALLOC(parmesh,pointIdxInterface_ToSend,nitem_ext_node,int,"pointIdxInterface_ToSend",ier = 0);

    /* Initialize arrays to send */
    // memset(tetraFlag_ToSend,0x00,4*mesh->ne*sizeof(int));

    /* STEP 1.1 - First loop to estimate the number of nodes and tetra to communicate
                  for the CALLOC of itosend, itorecv, rtosend and rtorecv */
    nt_in2out = np_in2out = 0;
    // nt_out2in = np_out2in = 0;

    /* Loop over the nodes in the external edge communicator */
    for (inode=0; inode < nitem_ext_node; inode++) {

      /* Get the indices of the nodes in internal communicators */
      idx_node_ext  = ext_node_comm->int_comm_index[inode];
      idx_node_int  = grp->node2int_node_comm_index2[idx_node_ext];
      idx_node_mesh = grp->node2int_node_comm_index1[idx_node_int];

      /* Add the flag 1 to these nodes*/
      p0 = &mesh->point[idx_node_mesh];
      p0->flag=-1;
      pointIdxInterface_ToSend[inode] = idx_node_mesh;
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
          nt_in2out++;
          break;
        }
      }

      /* If tetra is now MG_OVERLAP, then assign MG_OVERLAP to nodes */
      if (pt->tag & MG_OVERLAP) {

        for (i=0; i<4; i++) {
          ip = pt->v[i];
          p0 = &mesh->point[ip];
          tetraVertices_ToSend[4*(nt_in2out-1)+i]  = pt->v[i];
          if ( (p0->flag>=0)  && (p0->flag!=color_out_node+1) ) {
            tetraFlag_ToSend[4*(nt_in2out-1)+i]  = 1;
            p0->flag = color_out_node+1;
            pointCoord_ToSend[3*np_in2out]   = p0->c[0];
            pointCoord_ToSend[3*np_in2out+1] = p0->c[1];
            pointCoord_ToSend[3*np_in2out+2] = p0->c[2];
            pointTag_ToSend[np_in2out] = p0->tag+MG_OVERLAP;
            np_in2out++;
          }
          else if (p0->flag==-1) {
            tetraFlag_ToSend[4*(nt_in2out-1)+i]  = -1;
          }
          else {
            tetraFlag_ToSend[4*(nt_in2out-1)+i]  = 0;
          }
          tetraIdx_ToSend[nt_in2out] = k;
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

    /* First send and receive the number of points to exchange/share */
    MPI_CHECK(
      MPI_Sendrecv(&np_in2out,1,MPI_INT,color_out_node,MPI_OVERLAP_TAG+1,
                   &np_out2in,1,MPI_INT,color_out_node,MPI_OVERLAP_TAG+1,
                   comm,&status),return 0 );
    overlap->np_in2out = np_in2out;
    overlap->np_out2in = np_out2in;

    // fprintf(stdout, "OVERLAP Proc=%d-Proc=%d :: np_in2out=%d, np_out2in=%d\n",parmesh->myrank,color_out_node,np_in2out,np_out2in);

    MPI_CHECK(
      MPI_Sendrecv(&nt_in2out,1,MPI_INT,color_out_node,MPI_OVERLAP_TAG+2,
                   &nt_out2in,1,MPI_INT,color_out_node,MPI_OVERLAP_TAG+2,
                   comm,&status),return 0 );
    overlap->nt_in2out = nt_in2out;
    overlap->nt_out2in = nt_out2in;

    // fprintf(stdout, "OVERLAP Proc=%d-Proc=%d :: nt_in2out=%d, nt_out2in=%d\n",
    //                   parmesh->myrank,color_out_node,nt_in2out,nt_out2in);

    np_in = mesh->np;
    MPI_CHECK(
      MPI_Sendrecv(&np_in, 1,MPI_INT,color_out_node,MPI_OVERLAP_TAG+3,
                   &np_out,1,MPI_INT,color_out_node,MPI_OVERLAP_TAG+3,
                   comm,&status),return 0 );

    // fprintf(stdout, "OVERLAP Proc=%d-Proc=%d :: np_in=%d, np_out=%d, np_in2out=%d, np_out2in=%d\n",
    //                   parmesh->myrank,color_out_node,np_in,np_out,np_in2out,np_out2in);

    /* Alloc hash table for point correspondance */
    hash_in2out = overlap->hash_in2out;
    hash_out2in = overlap->hash_out2in;
    PMMG_CALLOC(parmesh,hash_in2out,np_in +np_out2in+1,int,"hash_in2out",ier = 0); // TO check if size OK
    PMMG_CALLOC(parmesh,hash_out2in,np_out          +1,int,"hash_out2in",ier = 0); // TO check if size OK

    // PMMG_CALLOC(parmesh,parmesh->overlap[i_commn].hash_in2out,np_in +np_out2in+1,int,"hash_in2out",ier = 0); // TO check if size OK
    // PMMG_CALLOC(parmesh,parmesh->overlap[i_commn].hash_out2in,np_out          +1,int,"hash_out2in",ier = 0); // TO check if size OK


    /* Send and receive the local indexes of nodes on interface */
    PMMG_CALLOC(parmesh,pointIdxInterface_ToRecv,nitem_ext_node,int,"pointIdxInterface_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointIdxInterface_ToSend,nitem_ext_node,MPI_INT,color_out_node,MPI_OVERLAP_TAG+4,
                   pointIdxInterface_ToRecv,nitem_ext_node,MPI_INT,color_out_node,MPI_OVERLAP_TAG+4,
                   comm,&status),return 0 );

    /* Fill the hash table of overlap with interface points first */
    for (inode=0; inode < nitem_ext_node; inode++) {
      ip_out = pointIdxInterface_ToRecv[inode]; // Index of nodes on color_out (the other partition)
      ip_in  = pointIdxInterface_ToSend[inode]; // Index of nodes on color_in  (this partition)
      hash_out2in[ip_out] = ip_in;  // From index on color_out, I found index on color_in
      hash_in2out[ip_in]  = ip_out; // From index on color_in,  I found index on color_out
    }

    // /* Send and receive the index of points to exchange/share */
    PMMG_CALLOC(parmesh,tetraVertices_ToRecv_outIdx,4*nt_out2in,int,"tetraVertices_ToRecv_outIdx",ier = 0);
    PMMG_CALLOC(parmesh,tetraVertices_ToRecv_inIdx, 4*nt_out2in,int,"tetraVertices_ToRecv_inIdx",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(tetraVertices_ToSend,       4*nt_in2out,MPI_INT,color_out_node,MPI_OVERLAP_TAG+5,
                   tetraVertices_ToRecv_outIdx,4*nt_out2in,MPI_INT,color_out_node,MPI_OVERLAP_TAG+5,
                   comm,&status),return 0 );

    /* Send and receive the index of points to exchange/share */
    PMMG_CALLOC(parmesh,tetraFlag_ToRecv,4*nt_out2in,int,"tetraFlag_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(tetraFlag_ToSend,4*nt_in2out,MPI_INT,color_out_node,MPI_OVERLAP_TAG+6,
                   tetraFlag_ToRecv,4*nt_out2in,MPI_INT,color_out_node,MPI_OVERLAP_TAG+6,
                   comm,&status),return 0 );

    /* Send and receive the tags of points to exchange/share */
    PMMG_CALLOC(parmesh,pointTag_ToRecv,np_out2in,uint16_t,"pointTag_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointTag_ToSend,np_in2out,MPI_UINT16_T,color_out_node,MPI_OVERLAP_TAG+7,
                   pointTag_ToRecv,np_out2in,MPI_UINT16_T,color_out_node,MPI_OVERLAP_TAG+7,
                   comm,&status),return 0 );

    /* Now send and receive the coordinates of the points */
    PMMG_CALLOC(parmesh,pointCoord_ToRecv,3*np_out2in,double,"pointCoord_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointCoord_ToSend,3*np_in2out,MPI_DOUBLE,color_out_node,MPI_OVERLAP_TAG+8,
                   pointCoord_ToRecv,3*np_out2in,MPI_DOUBLE,color_out_node,MPI_OVERLAP_TAG+8,
                   comm,&status),return 0 );

    // /* Fill the hash table of overlap with the other points */
    icoord=0;
    mesh->xpmax  = MG_MAX( (long long)(1.5*mesh->xp),mesh->npmax);

    for (inode=0; inode < 4*nt_out2in; inode++) {
      ip_out=tetraVertices_ToRecv_outIdx[inode];
      if (tetraFlag_ToRecv[inode]==1) {
        c_inode[0]=pointCoord_ToRecv[3*icoord];
        c_inode[1]=pointCoord_ToRecv[3*icoord+1];
        c_inode[2]=pointCoord_ToRecv[3*icoord+2];
        tag_inode =pointTag_ToRecv[icoord];
        icoord += 1;

        ip_in = MMG3D_newPt(mesh,c_inode,tag_inode,0);

        // if (parmesh->myrank==2 && color_out_node==1) {
          // fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: NEW           - inode=%d, ip_out=%d, ip_in=%d\n",
          //                  parmesh->myrank,color_out_node,inode,ip_out,ip_in);
        // }
        hash_out2in[ip_out] = ip_in;  // From index on color_out, I found index on color_in
        hash_in2out[ip_in]  = ip_out; // From index on color_in,  I found index on color_out
        tetraVertices_ToRecv_inIdx[inode]=ip_in;
      }
      else if (tetraFlag_ToRecv[inode]==-1) {
        tetraVertices_ToRecv_inIdx[inode]=hash_out2in[ip_out];
        // if (parmesh->myrank==2 && color_out_node==1) {
        //   fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: INTERFACE     - inode=%d, ip_out=%d, ip_in=%d\n",
        //                    parmesh->myrank,color_out_node,inode,ip_out,hash_out2in[ip_out]);
        // }
      }
      else{
        tetraVertices_ToRecv_inIdx[inode]=hash_out2in[ip_out];
        // if (parmesh->myrank==2 && color_out_node==1) {
        //   fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: ALREADY ADDED - inode=%d, ip_out=%d, ip_in=%d\n",
        //                    parmesh->myrank,color_out_node,inode,ip_out,hash_out2in[ip_out]);
        // }
      }
    }

    if (parmesh->myrank==0 && color_out_node==1) {
      for (inode=0; inode < 4*nt_out2in; inode+=4) {
        fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: inIdx=[%d-%d-%d-%d], outIdx=[%d-%d-%d-%d]\n",
                          parmesh->myrank,color_out_node,
                          tetraVertices_ToRecv_inIdx[inode],tetraVertices_ToRecv_inIdx[inode+1],tetraVertices_ToRecv_inIdx[inode+2],tetraVertices_ToRecv_inIdx[inode+3],
                          tetraVertices_ToRecv_outIdx[inode],tetraVertices_ToRecv_outIdx[inode+1],tetraVertices_ToRecv_outIdx[inode+2],tetraVertices_ToRecv_outIdx[inode+3]);
      }

      fprintf(stdout, "\n\n");

      for (inode=0; inode < mesh->np; inode++) {
        fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: hash_in2out - ip_in=%d, ip_out=%d\n",
                          parmesh->myrank,color_out_node,inode,hash_in2out[inode]);
      }

      fprintf(stdout, "\n\n");

      for (inode=0; inode < np_out; inode++) {
        fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: hash_out2in ip_out=%d, ip_in=%d\n",
                            parmesh->myrank,color_out_node,inode,hash_out2in[inode]);
      }
    }

    /* Deallocate memory*/
      // overlap->hash_in2out = hash_in2out;
      // overlap->hash_out2in = hash_out2in;
    // PMMG_DEL_MEM(parmesh,hash_in2out,int,"hash_in2out");
    // PMMG_DEL_MEM(parmesh,hash_out2in,int,"hash_out2in");

    PMMG_DEL_MEM(parmesh,pointCoord_ToSend,double,"pointCoord_ToSend");
    PMMG_DEL_MEM(parmesh,pointCoord_ToRecv,double,"pointCoord_ToRecv");

    PMMG_DEL_MEM(parmesh,pointTag_ToSend,uint16_t,"pointTag_ToSend");
    PMMG_DEL_MEM(parmesh,pointTag_ToRecv,uint16_t,"pointTag_ToRecv");

    PMMG_DEL_MEM(parmesh,pointIdxInterface_ToSend,int,"pointIdxInterface_ToSend");
    PMMG_DEL_MEM(parmesh,pointIdxInterface_ToRecv,int,"pointIdxInterface_ToRecv");

    PMMG_DEL_MEM(parmesh,tetraVertices_ToSend,int,"tetraVertices_ToSend");
    PMMG_DEL_MEM(parmesh,tetraVertices_ToRecv_inIdx, int,"tetraVertices_ToRecv_inIdx");
    PMMG_DEL_MEM(parmesh,tetraVertices_ToRecv_outIdx,int,"tetraVertices_ToRecv_outIdx");

    PMMG_DEL_MEM(parmesh,tetraFlag_ToSend,int,"tetraFlag_ToSend");
    PMMG_DEL_MEM(parmesh,tetraFlag_ToRecv,int,"tetraFlag_ToRecv");

    PMMG_DEL_MEM(parmesh,tetraIdx_ToSend,int,"tetraIdx_ToSend");

  }

  fprintf(stdout, "\n\n-------> END of OVERLAP \n");


  return 1;
}
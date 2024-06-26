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
#include "inlined_functions_3d_private.h"

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

  PMMG_pInt_comm int_comm;
  PMMG_pExt_comm ext_comm,ext_comm_ter;
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_pSol      met,ls;
  PMMG_pOverlap  overlap;
  MPI_Status     status;
  MMG5_pPoint    p0;
  MMG5_pTetra    pt;

  int ier = 1; /* initialize error */
  int LOCAL_MPI_TAG;
  int print_msg;

  int i,j,k,r;
  int inode,icoord;
  int icomm,icomm_ter;
  int iPBDY,rPBDY;

  MMG5_int ip;
  int ip_in,ip_out,ip_ter;
  uint16_t tag_inode;

  int np_out,np_in;
  int npTot_in2out,  npTot_out2in;
  int npInt_in2out,  npInt_out2in;
  int npPBDY_in2out, npPBDY_out2in;
  int ntTot_in2out, ntTot_out2in;
  int np_PBDY,np_PBDY_in,np_PBDY_out;
  int ndataPBDY_in2out,ndataPBDY_out2in;

  int idx_ext,idx_int;

  int nitem_int;
  int nitem_ext,nitem_ext_ter;
  int next_comm;

  int color_in,color_out,color_ter;

  double c_inode[3];

  int      *pointIdxInterface_ToSend, *pointIdxInterface_ToRecv;
  int      *tetraVertices_ToSend, *tetraVertices_ToRecv;
  int      *tetraVertices_ToRecv_outIdx, *tetraVertices_ToRecv_inIdx;
  int      *tetraFlag_ToSend, *tetraFlag_ToRecv;
  double   *pointCoord_ToSend, *pointCoord_ToRecv;
  double   *pointCoordPBDY_ToSend, *pointCoordPBDY_ToRecv;
  uint16_t *pointTag_ToSend, *pointTag_ToRecv;
  uint16_t *pointTagPBDY_ToSend, *pointTagPBDY_ToRecv;
  int      *pointIdxPBDY_ToSend;
  int      *dataPBDY_ToSend, *dataPBDY_ToRecv;
  int      *pointPBDY_ToRecv_inIdx, *point_PBDY_ToRecv_outIdx;
  int      *pointPBDY_added;
  int      *hash_in2out, *hash_out2in;

  int *n_ToSend;


  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
      fprintf(stdout,"\n      ## TODO:: PMMG_create_overlap.\n");

  /* Creation of overlap works only on packed mesh, i.e. on one group */
  /* Ensure only one group on each proc */
  assert(parmesh->ngrp == 1);

  /* Global initialization */
  grp  = &parmesh->listgrp[0];
  mesh =  parmesh->listgrp[0].mesh;
  met  =  parmesh->listgrp[0].met;

  int_comm  = parmesh->int_node_comm;
  next_comm = parmesh->next_node_comm;  // Number of communicator for nodes

  r = 0;

  /* Global allocation memory */
  PMMG_CALLOC(parmesh,parmesh->overlap,next_comm,PMMG_Overlap,"overlap",ier = 0);
  PMMG_CALLOC(parmesh,pointPBDY_added,mesh->np,int,"pointPBDY_added",ier = 0);

  /* Loop over the number of node communicator */
  for (icomm=0; icomm<next_comm; icomm++) {

    /* Local initialization */
    LOCAL_MPI_TAG = 1;
    npPBDY_in2out = 0;
    npInt_in2out  = 0;
    npTot_in2out  = 0;
    ntTot_in2out  = 0;
    np_PBDY    = 0;
    np_PBDY_in = 0;
    ndataPBDY_in2out = 0;

    print_msg = (parmesh->myrank==0); // && color_out==1);

    /* Get external node communicator information */
    ext_comm  = &parmesh->ext_node_comm[icomm]; // External node communicator
    color_in  = ext_comm->color_in;             // Color of this partition Pcolor_in
    color_out = ext_comm->color_out;            // Color of the remote partition Pcolor_out
    nitem_ext = ext_comm->nitem;                // Nbr of nodes in common between Pcolor_in and Pcolor_out

    /* Overlap variables */
    overlap = &parmesh->overlap[icomm];
    overlap->color_in  = color_in;
    overlap->color_out = color_out;

    /* Local allocation memory */
    PMMG_CALLOC(parmesh,pointCoord_ToSend,       3*mesh->np,double,  "pointCoord_ToSend",       ier = 0);
    PMMG_CALLOC(parmesh,pointCoordPBDY_ToSend,   3*mesh->np,double,  "pointCoordPBDY_ToSend",   ier = 0);
    PMMG_CALLOC(parmesh,pointTag_ToSend,         mesh->np,  uint16_t,"pointTag_ToSend",         ier = 0);
    PMMG_CALLOC(parmesh,pointTagPBDY_ToSend,     mesh->np,  uint16_t,"pointTagPBDY_ToSend",     ier = 0);
    PMMG_CALLOC(parmesh,pointIdxPBDY_ToSend,     mesh->np,  int,     "pointIdxPBDY_ToSend",     ier = 0);
    PMMG_CALLOC(parmesh,pointIdxInterface_ToSend,nitem_ext, int,     "pointIdxInterface_ToSend",ier = 0);
    PMMG_CALLOC(parmesh,tetraVertices_ToSend,    4*mesh->ne,int,     "tetraVertices_ToSend",    ier = 0);
    PMMG_CALLOC(parmesh,tetraFlag_ToSend,        4*mesh->ne,int,     "tetraFlag_ToSend",        ier = 0);
    PMMG_CALLOC(parmesh,n_ToSend,6,int,"n_ToSend",ier = 0);

    /** STEP 1 - Identification of the nodes between Pcolor_in and Pcolor_out.
              a) Assign flag -1 to these nodes
              b) Store the nodes index to be send to the Pcolor_out */
    /* Loop over the nodes in the external node communicator **/
    for (i=0; i < nitem_ext; i++) {
      /* Get the indices of the nodes in internal communicators */
      idx_ext = ext_comm->int_comm_index[i];
      idx_int = grp->node2int_node_comm_index2[idx_ext];
      ip      = grp->node2int_node_comm_index1[idx_int];

      /* Add the flag -1 to these nodes */
      p0 = &mesh->point[ip];
      p0->flag=-1;
      pointIdxInterface_ToSend[i] = ip;
    }

    /** STEP 2 - Identification of the points and the tetra to send to Pcolor_out.
              a) ntTot_in2out: nbr of tetra  from Pcolor_in to send to Pcolor_out.
              b) npInt_in2out: nbr of points from Pcolor_in to send to Pcolor_out.
              b) np_PBDY: nbr of points tag MG_PARBDY in Pcolor_in to send to Pcolor_out.
              c) tetraVertices_ToSend: list of tetra vertices to send to Pcolor_out.
              d) tetraFlag_ToSend: flag the tetra vertices using following rules
                    +1 means the vertex is seen for the first time
                      0 means the vertex has been already seen (from another tetra)
                    -1 means the vertex is on the interface between Pcolor_in and Pcolor_out
              e) pointCoord_ToSend: nodes coordinates to send to Pcolor_out
              f) pointTag_ToSend: nodes tags to send to Pcolor_out
              g) pointIdxPBDY_ToSend: nodes indexes w/ MG_PARBDY tag (except
                  those on interface between Pcolor_in & Pcolor_out) to send to Pcolor_out **/
    /* Loop over the tetra on this partition, i.e. Pcolor_in */
    for (k=1; k<=mesh->ne; k++) {
      pt  = &mesh->tetra[k];

      /* Loop over the vertices. If one vertex if flag -1, assign MG_OVERLAP to this tetra */
      for (i=0; i<4; i++) {
        ip = pt->v[i];
        p0 = &mesh->point[ip];
        if ( p0->flag < 0 ) {
          pt->tag |= MG_OVERLAP;
          break;
        }
      }

      /* If tetra is now MG_OVERLAP, then loop over the vertices */
      if (pt->tag & MG_OVERLAP) {
        for (i=0; i<4; i++) {
          ip = pt->v[i];
          p0 = &mesh->point[ip];
          tetraVertices_ToSend[4*ntTot_in2out+i] = ip;

          /* If this node has never been seen before
             and is not on interface between Pcolor_in and Pcolor_out,
             assign flag 1 to this vertex in tetraFlag_ToSend,
             store the coordinates and the tag of the node */
          if ( (p0->flag>=0)  && (p0->flag!=color_out+1) ) {

            /* Update flag of this vertex to identify it has already been seen */
            p0->flag = color_out+1;

            /* If this vertex is tag MG_PARBDY, store it in pointIdxPBDY_ToSend
               pointIdxPBDY_ToSend stores the points w/ MG_PARBDY tag except the ones
               at interface between Pcolor_in and Pcolor_out */
            if (p0->tag & MG_PARBDY) {
              tetraFlag_ToSend[4*ntTot_in2out+i] = 2;
              pointCoordPBDY_ToSend[3*npPBDY_in2out]   = p0->c[0];
              pointCoordPBDY_ToSend[3*npPBDY_in2out+1] = p0->c[1];
              pointCoordPBDY_ToSend[3*npPBDY_in2out+2] = p0->c[2];
              pointTagPBDY_ToSend[npPBDY_in2out]       = p0->tag+MG_OVERLAP;
              pointIdxPBDY_ToSend[npPBDY_in2out]       = ip;
              npPBDY_in2out++;
            }
            else {
              /* Store useful information */
              tetraFlag_ToSend[4*ntTot_in2out+i]  = 1;
              pointCoord_ToSend[3*npInt_in2out]   = p0->c[0];
              pointCoord_ToSend[3*npInt_in2out+1] = p0->c[1];
              pointCoord_ToSend[3*npInt_in2out+2] = p0->c[2];
              pointTag_ToSend[npInt_in2out]       = p0->tag+MG_OVERLAP;
              npInt_in2out++;
            }
          }
          /* If this node is on the interface between Pcolor_in and Pcolor_out,
             this node will be treated separately,
             so assign the flag -1 to this vertex in tetraFlag_ToSend */
          else if (p0->flag==-1) {
            tetraFlag_ToSend[4*ntTot_in2out+i]  = -1;
          }
          /* If this node has been seen before in another tetra,
             assigne the flag 0 to this vertex in tetraFlag_ToSend */
          else {
            tetraFlag_ToSend[4*ntTot_in2out+i]  = 0;
          }
        }
        ntTot_in2out++;
      }

      /* Remove the tag MG_OVERLAP - not needed here: here it is used as a way
         to identify the tetra to be send. Later, tetras having the tag MG_OVERLAP
         will actually be tetras received from Pcolor_out */
      pt->tag &= ~MG_OVERLAP;
    }

    npTot_in2out = npInt_in2out + npPBDY_in2out;

    /** STEP 3 - Reinitialise flags of points to 0 **/
    for (i=0; i < nitem_ext; i++) {
      idx_ext = ext_comm->int_comm_index[i];
      idx_int = grp->node2int_node_comm_index2[idx_ext];
      ip      = grp->node2int_node_comm_index1[idx_int];
      p0 = &mesh->point[ip];
      p0->flag = 0;
    }

    /** STEP 4 - Special treatment for nodes with tag MG_PARBDY.
                 Identification of the points at the interface between Pcolor_in
                 and another Pcolor_ter (!=Pcolor_out) to be send to Pcolor_out
                 dataPBDY_ToSend: for points with MG_PARBDY, store
                  - Pcolor_in:  this partition color
                  - Pcolor_ter: the other partition w/ which this node is shared
                        (Pcolor_ter!=Pcolor_out)
                  - The position of the point in the external communicator
                  - The local index in the mesh
                NB: This part of the algo might be changed and optimised **/

    /* Allocate the variables dataPBDY_ToSend */
    PMMG_CALLOC(parmesh,dataPBDY_ToSend,4*npPBDY_in2out*parmesh->nprocs,int,"dataPBDY_ToSend",ier = 0);

    /* Loop over the nodes with MG_PARBDY tags */
    for (i=0; i<npPBDY_in2out; i++) {
      ip = pointIdxPBDY_ToSend[i];
      p0 = &mesh->point[ip];

      /* Loop over the partitions having nodes in common w/ Pcolor_in */
      for (icomm_ter=0; icomm_ter<next_comm; icomm_ter++) {

        /* Search on the other partition than Pcolor_out */
        if (icomm_ter != icomm) {

          /* Get external node communicator information */
          ext_comm_ter  = &parmesh->ext_node_comm[icomm_ter]; // External node communicator
          color_ter     = ext_comm_ter->color_out;            // Color of the remote partition Pcolor_ter
          nitem_ext_ter = ext_comm_ter->nitem;                // Nbr of nodes in common between Pcolor_in and Pcolor_ter

          /* Loop over the nodes in the external node communicator Pcolor_ter */
          for (j=0; j < nitem_ext_ter; j++) {
            /* Get the indices of the nodes in internal communicators */
            idx_ext = ext_comm_ter->int_comm_index[j];
            idx_int = grp->node2int_node_comm_index2[idx_ext];
            ip_ter  = grp->node2int_node_comm_index1[idx_int];

            /* Each time the node ip is found being shared w/ another
               partition store Pcolor_in, Pcolor_ter, s and the
               position j of this node in the external comm of Pcolor_ter */
            if (ip==ip_ter) {
              dataPBDY_ToSend[4*ndataPBDY_in2out]   = color_in;
              dataPBDY_ToSend[4*ndataPBDY_in2out+1] = color_ter;
              dataPBDY_ToSend[4*ndataPBDY_in2out+2] = j;  // position in external communicator
              dataPBDY_ToSend[4*ndataPBDY_in2out+3] = ip; // index of point on this partition
              fprintf(stdout, "OVERLAP Proc=%d :: dataPBDY_ToSend=[%d-%d-%d-%d] \n", color_in,
                      dataPBDY_ToSend[4*ndataPBDY_in2out],dataPBDY_ToSend[4*ndataPBDY_in2out+1],
                      dataPBDY_ToSend[4*ndataPBDY_in2out+2],dataPBDY_ToSend[4*ndataPBDY_in2out+3]);
              ndataPBDY_in2out++;
              break;
            }
          }
        }
      }
    }
    fprintf(stdout, "-------------------------------------------- \n\n");

    n_ToSend[0]=npInt_in2out;     // Nbr of interior  point from Pcolor_in to send to Pcolor_out
    n_ToSend[1]=npPBDY_in2out;    // Nbr of MG_PARBDY point from Pcolor_in to send to Pcolor_out
    n_ToSend[2]=npTot_in2out;     // Total nbr of points from Pcolor_in to send to Pcolor_out
    n_ToSend[3]=ndataPBDY_in2out; // Nbr of data for MG_PARBDY points from Pcolor_in to send to Pcolor_out
    n_ToSend[4]=np_in;            // Total nbr of points on mesh Pcolor_in
    n_ToSend[5]=ntTot_in2out;     // Total nbr of tetras from Pcolor_in to send to Pcolor_out

    /** STEP 5 - Send and Receive all the data from the other partitions **/
    /* First send and receive the number of points to exchange/share */
    // TODO :: Do one comme instead of several
    MPI_CHECK(
      MPI_Sendrecv(&npTot_in2out,1,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   &npTot_out2in,1,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;
    overlap->np_in2out = npTot_in2out;
    overlap->np_out2in = npTot_out2in;

    MPI_CHECK(
      MPI_Sendrecv(&ntTot_in2out,1,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   &ntTot_out2in,1,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;
    overlap->nt_in2out = ntTot_in2out;
    overlap->nt_out2in = ntTot_out2in;

    np_in = mesh->np;
    MPI_CHECK(
      MPI_Sendrecv(&np_in, 1,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   &np_out,1,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    MPI_CHECK(
      MPI_Sendrecv(&ndataPBDY_in2out, 1,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   &ndataPBDY_out2in,1,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    // TODO :: Double check all size of variables
    /* Alloc hash table for point correspondance */
    hash_in2out = overlap->hash_in2out;
    hash_out2in = overlap->hash_out2in;
    PMMG_CALLOC(parmesh,hash_in2out,np_in +npTot_out2in+1,int,"hash_in2out",ier = 0); // TpointIdxInterface_ToRecvO check if size OK
    PMMG_CALLOC(parmesh,hash_out2in,np_out          +1,int,"hash_out2in",ier = 0); // TO check if size OK

    /* Send and receive the local indexes of nodes on interface */
    PMMG_CALLOC(parmesh,dataPBDY_ToRecv,4*np_PBDY_out,int,"dataPBDY_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(dataPBDY_ToSend,4*ndataPBDY_in2out,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   dataPBDY_ToRecv,4*ndataPBDY_out2in,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /* Send and receive the local indexes of nodes on interface */
    PMMG_CALLOC(parmesh,pointIdxInterface_ToRecv,nitem_ext,int,"pointIdxInterface_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointIdxInterface_ToSend,nitem_ext,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   pointIdxInterface_ToRecv,nitem_ext,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    // /* Send and receive the index of points to exchange/share */
    PMMG_CALLOC(parmesh,tetraVertices_ToRecv_outIdx,4*ntTot_out2in,int,"tetraVertices_ToRecv_outIdx",ier = 0);
    PMMG_CALLOC(parmesh,tetraVertices_ToRecv_inIdx, 4*ntTot_out2in,int,"tetraVertices_ToRecv_inIdx",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(tetraVertices_ToSend,       4*ntTot_in2out,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   tetraVertices_ToRecv_outIdx,4*ntTot_out2in,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /* Send and receive the index of points to exchange/share */
    PMMG_CALLOC(parmesh,tetraFlag_ToRecv,4*ntTot_out2in,int,"tetraFlag_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(tetraFlag_ToSend,4*ntTot_in2out,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   tetraFlag_ToRecv,4*ntTot_out2in,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /* Send and receive the tags of points to exchange/share */
    PMMG_CALLOC(parmesh,pointTag_ToRecv,npInt_out2in,uint16_t,"pointTag_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointTag_ToSend,npInt_in2out,MPI_UINT16_T,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   pointTag_ToRecv,npInt_out2in,MPI_UINT16_T,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /* Now send and receive the coordinates of the points */
    PMMG_CALLOC(parmesh,pointCoord_ToRecv,3*npInt_out2in,double,"pointCoord_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointCoord_ToSend,3*npInt_in2out,MPI_DOUBLE,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   pointCoord_ToRecv,3*npInt_out2in,MPI_DOUBLE,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    PMMG_CALLOC(parmesh,pointCoordPBDY_ToRecv,3*npPBDY_out2in,double,"pointCoordPBDY_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointCoordPBDY_ToSend,3*npPBDY_in2out,MPI_DOUBLE,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   pointCoordPBDY_ToRecv,3*npPBDY_out2in,MPI_DOUBLE,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /** STEP 6 - BLABLALBA **/
    /* Fill the hash table of overlap with interface points first */
    for (inode=0; inode < nitem_ext; inode++) {
      ip_out = pointIdxInterface_ToRecv[inode]; // Index of nodes on color_out (the other partition)
      ip_in  = pointIdxInterface_ToSend[inode]; // Index of nodes on color_in  (this partition)
      hash_out2in[ip_out] = ip_in;  // From index on color_out, I found index on color_in
      hash_in2out[ip_in]  = ip_out; // From index on color_in,  I found index on color_out
    }

    // TODO :: Take care of PBDY nodes first and remove it from the `big` loop

    /* Fill the hash table of overlap with the other points */
    icoord=0;
    mesh->xpmax  = MG_MAX( (long long)(1.5*mesh->xp),mesh->npmax);

    for (inode=0; inode < 4*ntTot_out2in; inode++) {
      ip_out=tetraVertices_ToRecv_outIdx[inode];
      if (tetraFlag_ToRecv[inode]==1) {
        c_inode[0]=pointCoord_ToRecv[3*icoord];
        c_inode[1]=pointCoord_ToRecv[3*icoord+1];
        c_inode[2]=pointCoord_ToRecv[3*icoord+2];
        tag_inode =pointTag_ToRecv[icoord];
        icoord += 1;

        if (tag_inode & MG_PARBDY) {
          for (iPBDY=0; iPBDY < np_PBDY_out; iPBDY++) {
            int ip_test        = dataPBDY_ToRecv[4*iPBDY+3];
            color_ter          = dataPBDY_ToRecv[4*iPBDY+1];
            int loc_test       = dataPBDY_ToRecv[4*iPBDY+2];
            int min_color_test = MG_MIN(color_out,color_ter);
            int max_color_test = MG_MAX(color_out,color_ter);

            /* Search if this point has already been added */
            int duplicated_point = 0;
            if (ip_test == ip_out) {
              for (rPBDY=0; rPBDY < r; rPBDY++) {
                if ( (color_ter == pointPBDY_added[4*rPBDY]) || (color_ter == pointPBDY_added[4*rPBDY+1]) ) {
                  if (loc_test == pointPBDY_added[4*rPBDY+2]) {
                    ip_in = pointPBDY_added[4*rPBDY+3];
                    hash_out2in[ip_out] = ip_in;  // From index on color_out, I found index on color_in
                    hash_in2out[ip_in]  = ip_out; // From index on color_in,  I found index on color_out
                    tetraVertices_ToRecv_inIdx[inode]=ip_in;
                    duplicated_point =1;
                    break;
                  }
                }
              }

              /* If the point has not been added, yet, then add it */
              if (duplicated_point==0) {
                if ( (ip_test == ip_out) & (color_ter != color_in) ) {
                  ip_in = MMG3D_newPt(mesh,c_inode,tag_inode,0);
                  hash_out2in[ip_out] = ip_in;  // From index on color_out, I found index on color_in
                  hash_in2out[ip_in]  = ip_out; // From index on color_in,  I found index on color_out
                  tetraVertices_ToRecv_inIdx[inode]=ip_in;
                  pointPBDY_added[4*r]  =min_color_test;
                  pointPBDY_added[4*r+1]=max_color_test;
                  pointPBDY_added[4*r+2]=loc_test;
                  pointPBDY_added[4*r+3]=ip_in;
                  r++;
                }
                break;
              }
            }
          }
        }
        else {
          ip_in = MMG3D_newPt(mesh,c_inode,tag_inode,0);
          hash_out2in[ip_out] = ip_in;  // From index on color_out, I found index on color_in
          hash_in2out[ip_in]  = ip_out; // From index on color_in,  I found index on color_out
          tetraVertices_ToRecv_inIdx[inode]=ip_in;
        }
      }
      else if (tetraFlag_ToRecv[inode]==-1) {
        tetraVertices_ToRecv_inIdx[inode]=hash_out2in[ip_out];
      //   fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: INTERFACE     - inode=%d, ip_out=%d, ip_in=%d\n",
      //                    color_in,color_out,inode,ip_out,hash_out2in[ip_out]);
      }
      else{
        tetraVertices_ToRecv_inIdx[inode]=hash_out2in[ip_out];
      //   fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: ALREADY ADDED - inode=%d, ip_out=%d, ip_in=%d\n",
      //                    color_in,color_out,inode,ip_out,hash_out2in[ip_out]);
      }
    }

    if (print_msg) {
      for (inode=0; inode < 4*ntTot_out2in; inode+=4) {
        fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: inIdx=[%d-%d-%d-%d], outIdx=[%d-%d-%d-%d]\n",
                          color_in,color_out,
                          tetraVertices_ToRecv_inIdx[inode],tetraVertices_ToRecv_inIdx[inode+1],tetraVertices_ToRecv_inIdx[inode+2],tetraVertices_ToRecv_inIdx[inode+3],
                          tetraVertices_ToRecv_outIdx[inode],tetraVertices_ToRecv_outIdx[inode+1],tetraVertices_ToRecv_outIdx[inode+2],tetraVertices_ToRecv_outIdx[inode+3]);
      }

      fprintf(stdout, "\n\n");

      // for (inode=0; inode < mesh->np; inode++) {
      //   fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: hash_in2out - ip_in=%d, ip_out=%d\n",
      //                     color_in,color_out,inode,hash_in2out[inode]);
      // }

      // fprintf(stdout, "\n\n");

      // for (inode=0; inode < np_out; inode++) {
      //   fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: hash_out2in ip_out=%d, ip_in=%d\n",
      //                       color_in,color_out,inode,hash_out2in[inode]);
      // }
    }

    /* Add the tetra to the mesh */
    // int v0,v1,v2,v3;
    // MMG5_pTetra         ptnew;
    // // MMG5_xTetra         xt,xt1;

    // MMG5_int iel;
    // // pt0=&mesh->tetra[1];

    // for (int itetra=0; itetra < nt_out2in; itetra++) {
    //   v0 = tetraVertices_ToRecv_inIdx[4*itetra];
    //   v1 = tetraVertices_ToRecv_inIdx[4*itetra+1];
    //   v2 = tetraVertices_ToRecv_inIdx[4*itetra+2];
    //   v3 = tetraVertices_ToRecv_inIdx[4*itetra+3];

    //   iel = MMG3D_newElt(mesh);
    //   if ( !iel ) {
    //     MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
    //                         fprintf(stderr,"\n  ## Error: %s: unable to allocate"
    //                                 " a new element.\n",__func__);
    //                         MMG5_INCREASE_MEM_MESSAGE();
    //                         fprintf(stderr,"  Exit program.\n");
    //                         return 0);
    //   }
    //   ptnew = &mesh->tetra[iel];
    //   // memcpy(ptnew,pt0,sizeof(MMG5_Tetra));
    //   ptnew->v[0] = v0;
    //   ptnew->v[1] = v1;
    //   ptnew->v[2] = v2;
    //   ptnew->v[3] = v3;
    //   ptnew->qual   = MMG5_caltet(mesh,met,ptnew) ;//MMG5_orcal(mesh,met,iel);
    // }

    /* Deallocate memory*/
    PMMG_DEL_MEM(parmesh,pointCoord_ToSend,double,"pointCoord_ToSend");
    PMMG_DEL_MEM(parmesh,pointCoord_ToRecv,double,"pointCoord_ToRecv");

    PMMG_DEL_MEM(parmesh,pointTag_ToSend,uint16_t,"pointTag_ToSend");
    PMMG_DEL_MEM(parmesh,pointTag_ToRecv,uint16_t,"pointTag_ToRecv");

    PMMG_DEL_MEM(parmesh,pointIdxPBDY_ToSend,int,"pointIdxPBDY_ToSend");

    PMMG_DEL_MEM(parmesh,pointIdxInterface_ToSend,int,"pointIdxInterface_ToSend");
    PMMG_DEL_MEM(parmesh,pointIdxInterface_ToRecv,int,"pointIdxInterface_ToRecv");

    PMMG_DEL_MEM(parmesh,tetraVertices_ToSend,int,"tetraVertices_ToSend");
    PMMG_DEL_MEM(parmesh,tetraVertices_ToRecv_inIdx, int,"tetraVertices_ToRecv_inIdx");
    PMMG_DEL_MEM(parmesh,tetraVertices_ToRecv_outIdx,int,"tetraVertices_ToRecv_outIdx");

    PMMG_DEL_MEM(parmesh,tetraFlag_ToSend,int,"tetraFlag_ToSend");
    PMMG_DEL_MEM(parmesh,tetraFlag_ToRecv,int,"tetraFlag_ToRecv");

    PMMG_DEL_MEM(parmesh,dataPBDY_ToSend,int,"dataPBDY_ToSend");
    PMMG_DEL_MEM(parmesh,dataPBDY_ToRecv,int,"dataPBDY_ToRecv");

  }

  fprintf(stdout, "\n\n-------> END of OVERLAP \n");


  return 1;
}
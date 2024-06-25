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

  PMMG_pInt_comm int_node_comm;
  PMMG_pExt_comm ext_node_comm;
  PMMG_pGrp      grp,grp_init,grp_overlap;
  MMG5_pPoint    p0;
  MMG5_pMesh     mesh,meshOld;
  MMG5_pTetra    pt;
  MPI_Status     status;
  MPI_Request    *request;
  PMMG_pOverlap   overlap;
  MMG5_pSol  met,ls;


  int ier = 1; /* initialize error */

  MMG5_int ip,ngrp,ne_init,n2inc_max,f2ifc_max,base_init;
  MMG5_int nt_overlap, np_overlap;

  int i, k, ireq;

  int inode;
  int i_commn;

  int np_in2out, np_out2in;
  int nt_in2out, nt_out2in;
  int np_out,np_in;
  int print_msg;

  int nitem_int_node;
  int nitem_ext_node;
  int next_node_comm;

  int LOCAL_MPI_TAG;
  int np_PBDY_in,np_PBDY_out;

  int *ip_overlap[2];
  double c_inode[3];
  MMG5_int np_new;

  int      *pointIdxInterface_ToSend, *pointIdxInterface_ToRecv;
  int      *tetraVertices_ToSend, *tetraVertices_ToRecv;
  int      *tetraVertices_ToRecv_outIdx, *tetraVertices_ToRecv_inIdx;
  int      *tetraFlag_ToSend, *tetraFlag_ToRecv;
  int      *hash_overlap_ToSend, *hash_overlap_ToRecv;
  double   *pointCoord_ToSend, *pointCoord_ToRecv;
  // int      *pointIdx_ToSend, *pointIdx_ToRecv;
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
  int r;

  int i_commn2,np_PBDY;
  int color_in_node2, color_out_node2, nitem_ext_node2;
  PMMG_pExt_comm ext_node_comm2;
  int *pointIdxPBDY_ToSend;
  int *pointPBDY_ToRecv;
  int *pointPBDY_ToSend, *pointPBDY_ToRecv_inIdx, *point_PBDY_ToRecv,outIdx;
  int i_PBDY;
  int *pointPBDY_added;

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
      fprintf(stdout,"\n      ## TODO:: PMMG_create_overlap.\n");

  /* For now - creation of overlap works only on packed mesh, i.e. on one group */
  /* Ensure only one group on each proc */
  assert(parmesh->ngrp == 1);

  /* Initialization */
  grp  = &parmesh->listgrp[0];
  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  int_node_comm  = parmesh->int_node_comm;
  next_node_comm = parmesh->next_node_comm;  // Number of communicator for nodes

  request = NULL;
  ireq = 0;
  r = 0;

  /* Allocate and reset memory */
  PMMG_CALLOC(parmesh,parmesh->overlap,next_node_comm,PMMG_Overlap,"allocate PMMG_Overlap ",ier = 0);
  PMMG_CALLOC(parmesh,pointPBDY_added,mesh->np,int,"pointIdxInterface_ToRecv",ier = 0);

  /* Loop over the number of node communicator */
  for (i_commn=0; i_commn<next_node_comm; i_commn++) {

    /* Local Initialization */
    LOCAL_MPI_TAG = 1;
    print_msg = (parmesh->myrank==0); // && color_out_node==1);
    nt_in2out  = 0;
    np_in2out  = 0;
    np_PBDY    = 0;
    np_PBDY_in = 0;

    /* Get external edge communicator information */
    ext_node_comm  = &parmesh->ext_node_comm[i_commn]; // External node communicator
    color_in_node  = ext_node_comm->color_in;          // Color of the hosting proc - this proc
    color_out_node = ext_node_comm->color_out;         // Color of the remote  proc - the proc to exchange with
    nitem_ext_node = ext_node_comm->nitem;             // Number of nodes in common between these 2 procs

    /* Overlap variables */
    overlap = &parmesh->overlap[i_commn];
    overlap->color_in  = color_in_node;
    overlap->color_out = color_out_node;

    /* Allocate memory */
    PMMG_CALLOC(parmesh,tetraVertices_ToSend,4*mesh->ne,int,"tetraVertices_ToSend",ier = 0);
    PMMG_CALLOC(parmesh,tetraFlag_ToSend,4*mesh->ne,int,"tetraFlag_ToSend",ier = 0);
    PMMG_CALLOC(parmesh,pointTag_ToSend,mesh->np,uint16_t,"pointTag_ToSend",ier = 0);
    PMMG_CALLOC(parmesh,pointCoord_ToSend,3*mesh->np,double,"pointCoord_ToSend",ier = 0);
    PMMG_CALLOC(parmesh,pointIdxInterface_ToSend,nitem_ext_node,int,"pointIdxInterface_ToSend",ier = 0);
    // PMMG_CALLOC(parmesh,pointIdx_ToSend,mesh->np,int,"pointIdx_ToSend",ier = 0);
    PMMG_CALLOC(parmesh,pointIdxPBDY_ToSend,mesh->np,int,"pointIdxPBDY_ToSend",ier = 0);

    /** STEP 1 - Identification of the nodes between Proc_color_in and Proc_color_out.
                a) Assign flag -1 to these nodes
                b) Store the nodes index to be send to the Proc_color_out */
    /* Loop over the nodes in the external edge communicator **/
    for (inode=0; inode < nitem_ext_node; inode++) {

      /* Get the indices of the nodes in internal communicators */
      idx_node_ext  = ext_node_comm->int_comm_index[inode];
      idx_node_int  = grp->node2int_node_comm_index2[idx_node_ext];
      idx_node_mesh = grp->node2int_node_comm_index1[idx_node_int];

      /* Add the flag -1 to these nodes */
      p0 = &mesh->point[idx_node_mesh];
      p0->flag=-1;
      pointIdxInterface_ToSend[inode] = idx_node_mesh;
    }

    /** STEP 2 - Identification of the points and the tetra to send to Proc_color_out.
                a) nt_in2out: nbr of tetra  from Proc_color_in to send to Proc_color_out.
                b) np_in2out: nbr of points from Proc_color_in to send to Proc_color_out.
                b) np_PBDY: nbr of points tag MG_PARBDY in Proc_color_in to send to Proc_color_out.
                c) tetraVertices_ToSend: list of tetra vertices to send to Proc_color_out.
                d) tetraFlag_ToSend: flag the tetra vertices using following rules
                      +1 means the vertex is seen for the first time
                       0 means the vertex has been already seen (from another tetra)
                      -1 means the vertex is on the interface between Proc_color_in and Proc_color_out
                e) pointCoord_ToSend: nodes coordinates to send to Proc_color_out
                f) pointTag_ToSend: nodes tags to send to Proc_color_out
                g) pointIdxPBDY_ToSend: nodes indexes w/ MG_PARBDY tag to send to Proc_color_out **/
    /* Loop over the tetra on this partition, i.e. Proc_color_in */
    for (k=1; k<=mesh->ne; k++) {
      pt  = &mesh->tetra[k];

      /* Loop over tetra vertices. If one vertex if flag -1, assign MG_OVERLAP to this tetra */
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
          tetraVertices_ToSend[4*nt_in2out+i] = ip;

          /* If this node has never been seen before
             and is not on interface between Proc_color_in and Proc_color_out,
             assign flag 1 to this vertex in tetraFlag_ToSend,
             store the coordinates and the tag of the node */
          if ( (p0->flag>=0)  && (p0->flag!=color_out_node+1) ) {

            /* Update flag of this vertex to identify it has already been seen */
            p0->flag = color_out_node+1;

            /* Store useful information */
            tetraFlag_ToSend[4*nt_in2out+i]  = 1;
            pointCoord_ToSend[3*np_in2out]   = p0->c[0];
            pointCoord_ToSend[3*np_in2out+1] = p0->c[1];
            pointCoord_ToSend[3*np_in2out+2] = p0->c[2];
            pointTag_ToSend[np_in2out] = p0->tag+MG_OVERLAP;
            np_in2out++;

            /* If this vertex is tag MG_PARBDY, then ... */
            if (p0->tag & MG_PARBDY) {
              pointIdxPBDY_ToSend[np_PBDY] = ip;
              np_PBDY++;
            }
          }
          /* If this node is on the interface between Proc_color_in and Proc_color_out,
             this node will be treated separately,
             so assign the flag -1 to this vertex in tetraFlag_ToSend */
          else if (p0->flag==-1) {
            tetraFlag_ToSend[4*nt_in2out+i]  = -1;
          }
          /* If this node has been seen before in another tetra,
             assigne the flag 0 to this vertex in tetraFlag_ToSend */
          else {
            tetraFlag_ToSend[4*nt_in2out+i]  = 0;
          }
        }
        nt_in2out++;
      }

      /* Remove the tag MG_OVERLAP - not needed here: here it is used as a way
         to identify the tetra to be send. Later, tetras having the tag MG_OVERLAP
         will actually be tetras received from Proc_color_out */
      pt->tag &= ~MG_OVERLAP;
    }

    /** STEP 3 - Reinitialise flags of points to 0 **/
    for (inode=0; inode < nitem_ext_node; inode++) {
      idx_node_ext  = ext_node_comm->int_comm_index[inode];
      idx_node_int  = grp->node2int_node_comm_index2[idx_node_ext];
      idx_node_mesh = grp->node2int_node_comm_index1[idx_node_int];
      p0 = &mesh->point[idx_node_mesh];
      p0->flag = 0;
    }

    /** STEP 4 - Special treatment for nodes with tag MG_PARBDY.
                Identification of the points at the interface between Proc_color_in
                and other Proc_color_ter that need to be send to Proc_color_out
                pointPBDY_ToSend: store for each point w/ MG_PARBDY
                  - Proc_color_in
                  - Proc_color_ter
                NB: This part of the algo might be changed and optimised **/
    PMMG_CALLOC(parmesh,pointPBDY_ToSend,4*np_PBDY*parmesh->nprocs,int,"pointPBDY_ToSend",ier = 0);

    /* Loop over the number of */
    for (i_PBDY=0; i_PBDY<np_PBDY; i_PBDY++) {
      ip = pointIdxPBDY_ToSend[i_PBDY];
      p0 = &mesh->point[ip];

      /* Loop over the other partition having nodes in common with this partition */
      for (i_commn2=0; i_commn2<next_node_comm; i_commn2++) {
        if (i_commn2 != i_commn) {

          /* Get external edge communicator information */
          ext_node_comm2  = &parmesh->ext_node_comm[i_commn2]; // External node communicator
          color_out_node2 = ext_node_comm2->color_out;         // Color of the remote  proc - the proc to exchange with
          nitem_ext_node2 = ext_node_comm2->nitem;             // Number of nodes in common between these 2 procs

          /* Loop over the nodes in the external edge communicator */
          for (inode=0; inode < nitem_ext_node2; inode++) {
            /* Get the   indices of the nodes in internal communicators */
            idx_node_ext  = ext_node_comm2->int_comm_index[inode];
            idx_node_int  = grp->node2int_node_comm_index2[idx_node_ext];
            idx_node_mesh = grp->node2int_node_comm_index1[idx_node_int];

            if (ip==idx_node_mesh) {
              pointPBDY_ToSend[4*np_PBDY_in]   = color_in_node ; // min(color_in,color_out)
              pointPBDY_ToSend[4*np_PBDY_in+1] = color_out_node2;  // max(color_in,color_out)
              pointPBDY_ToSend[4*np_PBDY_in+2] = inode; // index in external communicator
              pointPBDY_ToSend[4*np_PBDY_in+3] = ip; // index of point on this partition
              fprintf(stdout, "OVERLAP Proc=%d :: pointPBDY_ToSend=[%d-%d-%d-%d] \n", color_in_node,
                          pointPBDY_ToSend[4*np_PBDY_in],pointPBDY_ToSend[4*np_PBDY_in+1],pointPBDY_ToSend[4*np_PBDY_in+2],pointPBDY_ToSend[4*np_PBDY_in+3]);
              np_PBDY_in++;
              break;
            }
          }
        }
      }
    }

    /** STEP 5 - Send and Receive all the data from the other partitions **/
    /* First send and receive the number of points to exchange/share */
    MPI_CHECK(
      MPI_Sendrecv(&np_in2out,1,MPI_INT,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   &np_out2in,1,MPI_INT,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;
    overlap->np_in2out = np_in2out;
    overlap->np_out2in = np_out2in;

    MPI_CHECK(
      MPI_Sendrecv(&nt_in2out,1,MPI_INT,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   &nt_out2in,1,MPI_INT,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;
    overlap->nt_in2out = nt_in2out;
    overlap->nt_out2in = nt_out2in;

    np_in = mesh->np;
    MPI_CHECK(
      MPI_Sendrecv(&np_in, 1,MPI_INT,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   &np_out,1,MPI_INT,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    MPI_CHECK(
      MPI_Sendrecv(&np_PBDY_in, 1,MPI_INT,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   &np_PBDY_out,1,MPI_INT,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /* Alloc hash table for point correspondance */
    hash_in2out = overlap->hash_in2out;
    hash_out2in = overlap->hash_out2in;
    PMMG_CALLOC(parmesh,hash_in2out,np_in +np_out2in+1,int,"hash_in2out",ier = 0); // TpointIdxInterface_ToRecvO check if size OK
    PMMG_CALLOC(parmesh,hash_out2in,np_out          +1,int,"hash_out2in",ier = 0); // TO check if size OK

    /* Send and receive the local indexes of nodes on interface */
    PMMG_CALLOC(parmesh,pointPBDY_ToRecv,4*np_PBDY_out,int,"pointIdxInterface_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointPBDY_ToSend,4*np_PBDY_in,MPI_INT,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   pointPBDY_ToRecv,4*np_PBDY_out,MPI_INT,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /* Send and receive the local indexes of nodes on interface */
    PMMG_CALLOC(parmesh,pointIdxInterface_ToRecv,nitem_ext_node,int,"pointIdxInterface_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointIdxInterface_ToSend,nitem_ext_node,MPI_INT,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   pointIdxInterface_ToRecv,nitem_ext_node,MPI_INT,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    // /* Send and receive the index of points to exchange/share */
    PMMG_CALLOC(parmesh,tetraVertices_ToRecv_outIdx,4*nt_out2in,int,"tetraVertices_ToRecv_outIdx",ier = 0);
    PMMG_CALLOC(parmesh,tetraVertices_ToRecv_inIdx, 4*nt_out2in,int,"tetraVertices_ToRecv_inIdx",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(tetraVertices_ToSend,       4*nt_in2out,MPI_INT,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   tetraVertices_ToRecv_outIdx,4*nt_out2in,MPI_INT,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /* Send and receive the index of points to exchange/share */
    PMMG_CALLOC(parmesh,tetraFlag_ToRecv,4*nt_out2in,int,"tetraFlag_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(tetraFlag_ToSend,4*nt_in2out,MPI_INT,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   tetraFlag_ToRecv,4*nt_out2in,MPI_INT,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /* Send and receive the tags of points to exchange/share */
    PMMG_CALLOC(parmesh,pointTag_ToRecv,np_out2in,uint16_t,"pointTag_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointTag_ToSend,np_in2out,MPI_UINT16_T,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   pointTag_ToRecv,np_out2in,MPI_UINT16_T,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /* Now send and receive the coordinates of the points */
    PMMG_CALLOC(parmesh,pointCoord_ToRecv,3*np_out2in,double,"pointCoord_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointCoord_ToSend,3*np_in2out,MPI_DOUBLE,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   pointCoord_ToRecv,3*np_out2in,MPI_DOUBLE,color_out_node,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /** STEP 6 - BLABLALBA **/
    /* Fill the hash table of overlap with interface points first */
    for (inode=0; inode < nitem_ext_node; inode++) {
      ip_out = pointIdxInterface_ToRecv[inode]; // Index of nodes on color_out (the other partition)
      ip_in  = pointIdxInterface_ToSend[inode]; // Index of nodes on color_in  (this partition)
      hash_out2in[ip_out] = ip_in;  // From index on color_out, I found index on color_in
      hash_in2out[ip_in]  = ip_out; // From index on color_in,  I found index on color_out
    }

    /* Fill the hash table of overlap with the other points */
    icoord=0;
    mesh->xpmax  = MG_MAX( (long long)(1.5*mesh->xp),mesh->npmax);

    int *pointIdxPBDY;
    PMMG_CALLOC(parmesh,pointIdxPBDY,np_PBDY*parmesh->nprocs,int,"pointPBDY_ToSend",ier = 0);
    for (inode=0; inode < 4*nt_out2in; inode++) {
      ip_out=tetraVertices_ToRecv_outIdx[inode];
      if (tetraFlag_ToRecv[inode]==1) {
        c_inode[0]=pointCoord_ToRecv[3*icoord];
        c_inode[1]=pointCoord_ToRecv[3*icoord+1];
        c_inode[2]=pointCoord_ToRecv[3*icoord+2];
        tag_inode =pointTag_ToRecv[icoord];
        icoord += 1;

        if (tag_inode & MG_PARBDY) {
          // if (parmesh->myrank==0) {
          //   fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: MG_PARBDY   - inode=%d, ip_out=%d, np_PBDY_out=%d \n",
          //                     parmesh->myrank,color_out_node,inode,ip_out,np_PBDY_out);
          // }
          for (i_PBDY=0; i_PBDY < np_PBDY_out; i_PBDY++) {
            int ip_test        = pointPBDY_ToRecv[4*i_PBDY+3];
            int color_ter_test = pointPBDY_ToRecv[4*i_PBDY+1];
            int loc_test       = pointPBDY_ToRecv[4*i_PBDY+2];
            int min_color_test = MG_MIN(color_out_node,color_ter_test);
            int max_color_test = MG_MAX(color_out_node,color_ter_test);
            // if (parmesh->myrank==0) {
            //   fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: MG_PARBDY   - i_PBDY=%d - ip_test=%d, color_ter_test=%d, loc_test=%d, min_color_test=%d, max_color_test=%d \n",
            //                   parmesh->myrank,color_out_node,i_PBDY,ip_test,color_ter_test,loc_test,min_color_test,max_color_test);
            // }

            /* Search if this point has already been added */
            int duplicated_point = 0;
            if (ip_test == ip_out) {
              for (int r_PBDY=0; r_PBDY < r; r_PBDY++) {
                if ( (color_ter_test == pointPBDY_added[4*r_PBDY]) || (color_ter_test == pointPBDY_added[4*r_PBDY+1]) ) {
                  if (loc_test == pointPBDY_added[4*r_PBDY+2]) {
                    ip_in = pointPBDY_added[4*r_PBDY+3];
                    hash_out2in[ip_out] = ip_in;  // From index on color_out, I found index on color_in
                    hash_in2out[ip_in]  = ip_out; // From index on color_in,  I found index on color_out
                    tetraVertices_ToRecv_inIdx[inode]=ip_in;
                    // if (parmesh->myrank==0) {
                    //   fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: MG_PARBDY   - This point exists in pointPBDY_added point pointPBDY_added=[%d-%d-%d-%d], ip_in=%d - ip_out=%d, hash_out2in[ip_out]=%d, hash_in2out[ip_in]=%d\n\n\n",
                    //                    parmesh->myrank,color_out_node,pointPBDY_added[4*r_PBDY],pointPBDY_added[4*r_PBDY+1],pointPBDY_added[4*r_PBDY+2],pointPBDY_added[4*r_PBDY+3],
                    //                    ip_in,ip_out,hash_out2in[ip_out],hash_in2out[ip_in]);
                    // }

                    duplicated_point =1;
                    break;
                  }
                }
              }

              /* If the point has not been added, yet, then add it */
              if (duplicated_point==0) {
                // if (parmesh->myrank==0) {
                //   fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: MG_PARBDY   - Add this point \n", parmesh->myrank,color_out_node);
                // }
                if ( (ip_test == ip_out) & (color_ter_test != color_in_node) ) {
                  ip_in = MMG3D_newPt(mesh,c_inode,tag_inode,0);
                  hash_out2in[ip_out] = ip_in;  // From index on color_out, I found index on color_in
                  hash_in2out[ip_in]  = ip_out; // From index on color_in,  I found index on color_out
                  tetraVertices_ToRecv_inIdx[inode]=ip_in;
                  pointPBDY_added[4*r]  =min_color_test;
                  pointPBDY_added[4*r+1]=max_color_test;
                  pointPBDY_added[4*r+2]=loc_test;
                  pointPBDY_added[4*r+3]=ip_in;
                  // if (parmesh->myrank==0) {
                  //   fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: MG_PARBDY   - Add this point pointPBDY_added=[%d-%d-%d-%d] \n\n\n",
                  //                   parmesh->myrank,color_out_node,pointPBDY_added[4*r],pointPBDY_added[4*r+1],pointPBDY_added[4*r+2],pointPBDY_added[4*r+3]);
                  // }
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
        // if (parmesh->myrank==2 && color_out_node==1) {
          // fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: NEW           - inode=%d, ip_out=%d, ip_in=%d\n",
          //                  parmesh->myrank,color_out_node,inode,ip_out,ip_in);
        // }
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

    if (print_msg) {
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
    // overlap->hash_in2out = hash_in2out;
    // overlap->hash_out2in = hash_out2in;
    // PMMG_DEL_MEM(parmesh,hash_in2out,int,"hash_in2out");
    // PMMG_DEL_MEM(parmesh,hash_out2in,int,"hash_out2in");

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

    PMMG_DEL_MEM(parmesh,pointPBDY_ToSend,int,"pointPBDY_ToSend");

  }

  fprintf(stdout, "\n\n-------> END of OVERLAP \n");


  return 1;
}
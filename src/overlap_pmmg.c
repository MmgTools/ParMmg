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
 * \param comm pointer toward ...
 *
 * \return 1 if success, 0 if fail.
 *
 * \todo Fill the funtion
 *
 * TODO
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
  int ip_in,ip_out,ip_ter,ip;
  int icoord,icomm,icomm_ter;
  int src,ref;

  int ne_in,np_out,np_in;
  int npTot_in2out, npTot_out2in;
  int npInt_in2out, npInt_out2in;
  int npPBDY_in2out, npPBDY_out2in;
  int ntTot_in2out, ntTot_out2in;
  int ndataPBDY_in2out, ndataPBDY_out2in;
  int ndataPBDY_added;

  int loc,duplicated_point;
  int min_color, max_color;

  int idx_ext,idx_int;
  int nitem_ext,nitem_ext_ter,next_comm;
  int color_in,color_out,color_ter;

  double coord[3];

  uint16_t tag;

  int      *n_ToSend, *n_ToRecv;
  int      *pointIdxInterface_ToSend, *pointIdxInterface_ToRecv;
  int      *tetraVertices_ToSend, *tetraVertices_ToRecv;
  int      *tetraVertices_ToRecv_outIdx, *tetraVertices_ToRecv_inIdx;
  int      *tetraFlag_ToSend, *tetraFlag_ToRecv;
  int      *pointIdxPBDY_ToSend, *pointIdxPBDY_ToRecv;
  int      *pointRefPBDY_ToSend, *pointRefPBDY_ToRecv;
  int      *pointRefInt_ToSend, *pointRefInt_ToRecv;
  int      *dataPBDY_ToSend, *dataPBDY_ToRecv;
  int      *pointPBDY_ToRecv_inIdx, *point_PBDY_ToRecv_outIdx;
  int      *dataPBDY_AlreadyAdded;
  int      *hash_in2out, *hash_out2in;
  double   *pointCoordInt_ToSend, *pointCoordInt_ToRecv;
  double   *pointCoordPBDY_ToSend, *pointCoordPBDY_ToRecv;
  uint16_t *pointTagInt_ToSend,  *pointTagInt_ToRecv;
  uint16_t *pointTagPBDY_ToSend, *pointTagPBDY_ToRecv;

  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
      fprintf(stdout,"\n      ## Create Overlap.\n");
      fprintf(stdout,"\n      ## TODO:: PMMG_create_overlap.\n");
  }

  /* Creation of overlap works only on packed mesh, i.e. on one group */
  /* Ensure only one group on each proc */
  assert(parmesh->ngrp == 1);

  /* Global initialization */
  grp  = &parmesh->listgrp[0];
  mesh =  parmesh->listgrp[0].mesh;
  met  =  parmesh->listgrp[0].met;
  ls   =  parmesh->listgrp[0].ls;

  int_comm  = parmesh->int_node_comm;
  next_comm = parmesh->next_node_comm;  // Number of nodes communicator

  ndataPBDY_added = 0;
  ne_in = mesh->ne;

  /* Global allocation memory */
  PMMG_CALLOC(parmesh,parmesh->overlap,next_comm,PMMG_Overlap,"overlap",ier = 0);
  PMMG_CALLOC(parmesh,dataPBDY_AlreadyAdded,5*mesh->np,int,"dataPBDY_AlreadyAdded",ier = 0);

  /* Loop over the number of node communicator */
  for (icomm=0; icomm<next_comm; icomm++) {

    /* Local initialization */
    LOCAL_MPI_TAG = 1;
    npPBDY_in2out = 0;
    npInt_in2out  = 0;
    npTot_in2out  = 0;
    ntTot_in2out  = 0;
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
    PMMG_CALLOC(parmesh,pointCoordInt_ToSend,    3*mesh->np,double,  "pointCoordInt_ToSend",    ier = 0);
    PMMG_CALLOC(parmesh,pointCoordPBDY_ToSend,   3*mesh->np,double,  "pointCoordPBDY_ToSend",   ier = 0);
    PMMG_CALLOC(parmesh,pointTagInt_ToSend,      mesh->np,  uint16_t,"pointTagInt_ToSend",      ier = 0);
    PMMG_CALLOC(parmesh,pointRefInt_ToSend,      mesh->np,  int,     "pointRefInt_ToSend",     ier = 0);
    PMMG_CALLOC(parmesh,pointTagPBDY_ToSend,     mesh->np,  uint16_t,"pointTagPBDY_ToSend",     ier = 0);
    PMMG_CALLOC(parmesh,pointRefPBDY_ToSend,     mesh->np,  int,     "pointRefPBDY_ToSend",     ier = 0);
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
              b) npInt_in2out: nbr of points (except those tag MG_PARBDY) from Pcolor_in to send to Pcolor_out.
              b) npPBDY_in2out: nbr of points tag MG_PARBDY from Pcolor_in to send to Pcolor_out.
              c) tetraVertices_ToSend: list of tetra vertices to send to Pcolor_out.
              d) tetraFlag_ToSend: flag the tetra vertices using following rules
                    +1 means the vertex is seen for the first time
                      0 means the vertex has been already seen (from another tetra)
                        OR is on the interface between Pcolor_in and Pcolor_out
                        OR is tag MG_PARBDY (special treatments for those)
              e) pointCoordInt_ToSend:  nodes coordinates to send to Pcolor_out
              f) pointCoordPBDY_ToSend: nodes (tagged MG_PARBDY) coordinates to send to Pcolor_out
              g) pointTagInt_ToSend: nodes tags to send to Pcolor_out
              h) pointTagPBDY_ToSend: nodes (tagged MG_PARBDY) tags to send to Pcolor_out
              i) pointIdxPBDY_ToSend: nodes indexes w/ MG_PARBDY tag (except
                  those on interface between Pcolor_in & Pcolor_out) to send to Pcolor_out **/
    /* Loop over the tetra on this partition, i.e. Pcolor_in */
    for (k=1; k<=ne_in; k++) {
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
        // TODO : add ref of tetra in this table tetraVertices_ToSend
        // ref = pt->ref;
        // tetraVertices_ToSend[5*ntTot_in2out+5] = ref;

        for (i=0; i<4; i++) {
          ip  = pt->v[i];
          p0  = &mesh->point[ip];
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
              tetraFlag_ToSend[4*ntTot_in2out+i]       = 0;
              pointCoordPBDY_ToSend[3*npPBDY_in2out]   = p0->c[0];
              pointCoordPBDY_ToSend[3*npPBDY_in2out+1] = p0->c[1];
              pointCoordPBDY_ToSend[3*npPBDY_in2out+2] = p0->c[2];
              pointTagPBDY_ToSend[npPBDY_in2out]       = p0->tag;
              pointRefPBDY_ToSend[npPBDY_in2out]       = p0->ref;
              pointIdxPBDY_ToSend[npPBDY_in2out]       = ip;
              npPBDY_in2out++;
            }
            else {
              /* Store useful information */
              tetraFlag_ToSend[4*ntTot_in2out+i]     = 1;
              pointCoordInt_ToSend[3*npInt_in2out]   = p0->c[0];
              pointCoordInt_ToSend[3*npInt_in2out+1] = p0->c[1];
              pointCoordInt_ToSend[3*npInt_in2out+2] = p0->c[2];
              pointTagInt_ToSend[npInt_in2out]       = p0->tag;
              pointRefInt_ToSend[npInt_in2out]       = p0->ref;
              npInt_in2out++;
            }
          }
          /* If this node is on the interface between Pcolor_in and Pcolor_out,
             this node will be treated separately, assign 0 in tetraFlag_ToSend
             OR
             If this node has been seen before in another tetra,
             assign 0 in tetraFlag_ToSend */
          else {
            tetraFlag_ToSend[4*ntTot_in2out+i]  = 0;
          }
        }
        ntTot_in2out++;
        /* Remove the tag MG_OVERLAP - not needed here: here it is used as a way
          to identify the tetra to be send. Later, tetras having the tag MG_OVERLAP
          will actually be tetras received from Pcolor_out */
        pt->tag &= ~MG_OVERLAP;
      }
    }

    /* The total number of point to send to Pcolor_out npTot_in2out is
       the nbr of point not tagged MG_PARBDY npInt_in2out
       plus the ones tagged MG_PARBDY npPBDY_in2out */
    npTot_in2out = npInt_in2out + npPBDY_in2out;

    /** STEP 3 - Reinitialise flags of points to 0 **/
    for (i=0; i < nitem_ext; i++) {
      idx_ext = ext_comm->int_comm_index[i];
      idx_int = grp->node2int_node_comm_index2[idx_ext];
      ip      = grp->node2int_node_comm_index1[idx_int];
      p0 = &mesh->point[ip];
      p0->flag = 0;
      if (p0->tag & MG_OVERLAP) p0->tag =~ MG_OVERLAP;
    }

    /** STEP 4 - Special treatment for nodes with tag MG_PARBDY.
                 Identification of the points at the interface between Pcolor_in
                 and another Pcolor_ter (!=Pcolor_out) to be send to Pcolor_out
                 dataPBDY_ToSend: for points with MG_PARBDY, store
                  - Pcolor_ter: the other partition w/ which this node is shared
                        (Pcolor_ter!=Pcolor_out)
                  - The position of the point in the external communicator
                  - The local index in the mesh
                NB: This part of the algo might be changed and optimised **/

    /* Allocate the variables dataPBDY_ToSend */
    PMMG_CALLOC(parmesh,dataPBDY_ToSend,3*npPBDY_in2out*parmesh->nprocs,int,"dataPBDY_ToSend",ier = 0);

    /* Loop over the nodes with MG_PARBDY tags */
    for (i=0; i<npPBDY_in2out+1; i++) {
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
              dataPBDY_ToSend[3*ndataPBDY_in2out]   = color_ter;
              dataPBDY_ToSend[3*ndataPBDY_in2out+1] = j;  // position in external communicator
              dataPBDY_ToSend[3*ndataPBDY_in2out+2] = ip; // index of point on this partition
              // fprintf(stdout, "OVERLAP Proc=%d :: color_ter=%d, dataPBDY_ToSend=[%d-%d-%d] \n", color_in,color_ter,
              //         dataPBDY_ToSend[3*ndataPBDY_in2out],dataPBDY_ToSend[3*ndataPBDY_in2out+1],dataPBDY_ToSend[3*ndataPBDY_in2out+2]);
              ndataPBDY_in2out++;
              break;
            }
          }
        }
      }
    }

    np_in = mesh->np;

    n_ToSend[0]=npInt_in2out;     // Nbr of interior  point from Pcolor_in to send to Pcolor_out
    n_ToSend[1]=npPBDY_in2out;    // Nbr of MG_PARBDY point from Pcolor_in to send to Pcolor_out
    n_ToSend[2]=npTot_in2out;     // Total nbr of points from Pcolor_in to send to Pcolor_out
    n_ToSend[3]=ndataPBDY_in2out; // Nbr of data for MG_PARBDY points from Pcolor_in to send to Pcolor_out
    n_ToSend[4]=np_in;            // Total nbr of points on mesh Pcolor_in
    n_ToSend[5]=ntTot_in2out;     // Total nbr of tetras from Pcolor_in to send to Pcolor_out

    /** STEP 5 - Send and Receive all the data from the other partitions **/
    // TODO :: Improve number of communications - 10 comm for now -> might be reduced to 7
    /* STEP 5.1 - First send and receive the different sizes to exchange/share */
    PMMG_CALLOC(parmesh,n_ToRecv,6,int,"n_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(n_ToSend,6,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   n_ToRecv,6,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    npInt_out2in     = n_ToRecv[0];
    npPBDY_out2in    = n_ToRecv[1];
    npTot_out2in     = n_ToRecv[2];
    ndataPBDY_out2in = n_ToRecv[3];
    np_out           = n_ToRecv[4];
    ntTot_out2in     = n_ToRecv[5];

    // fprintf(stdout, "OVERLAP Proc %d-Proc %d :: npInt_in2out=%d, npPBDY_in2out=%d, npTot_in2out=%d -- npInt_out2in=%d, npPBDY_out2in=%d, npTot_out2in=%d  \n", color_in,color_out,
    //         npInt_in2out,npPBDY_in2out,npTot_in2out,npInt_out2in,npPBDY_out2in,npTot_out2in);

    // fprintf(stdout, "OVERLAP Proc %d-Proc %d :: ntTot_in2out=%d, ntTot_out2in=%d  \n", color_in,color_out,
    //         ntTot_in2out,ntTot_out2in);

    overlap->np_in2out = npTot_in2out;
    overlap->np_out2in = npTot_out2in;
    overlap->nt_in2out = ntTot_in2out;
    overlap->nt_out2in = ntTot_out2in;

    /* STEP 5.2 - Alloc hash table for nodes index link */
    hash_in2out = overlap->hash_in2out;
    hash_out2in = overlap->hash_out2in;
    PMMG_CALLOC(parmesh,hash_in2out,np_in +npTot_out2in+1,int,"hash_in2out",ier = 0);
    PMMG_CALLOC(parmesh,hash_out2in,np_out             +1,int,"hash_out2in",ier = 0);

    /* STEP 5.3 - Send and receive all the other data */
    /* Send/receive data needed to identify MG_PARBDY nodes */
    PMMG_CALLOC(parmesh,dataPBDY_ToRecv,3*ndataPBDY_out2in,int,"dataPBDY_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(dataPBDY_ToSend,3*ndataPBDY_in2out,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   dataPBDY_ToRecv,3*ndataPBDY_out2in,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    // if (color_in == 0) {
    //   for (i=0; i<ndataPBDY_out2in; i++) {
    //     fprintf(stdout, "OVERLAP Proc=%d :: color_out=%d, dataPBDY_ToRecv=[%d-%d-%d] \n", color_in,color_out,
    //             dataPBDY_ToRecv[3*i],dataPBDY_ToRecv[3*i+1],dataPBDY_ToRecv[3*i+2]);
    //   }
    // }

    /* Send/receive nodes index located on interface between Pcolor_in and Pcolor_out */
    PMMG_CALLOC(parmesh,pointIdxInterface_ToRecv,nitem_ext,int,"pointIdxInterface_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointIdxInterface_ToSend,nitem_ext,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   pointIdxInterface_ToRecv,nitem_ext,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /* Send/receive nodes index located on interface between Pcolor_in and Pcolor_out */
    PMMG_CALLOC(parmesh,pointIdxPBDY_ToRecv,npPBDY_out2in,int,"pointIdxPBDY_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointIdxPBDY_ToSend,npPBDY_in2out,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   pointIdxPBDY_ToRecv,npPBDY_out2in,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /* Send/receive nodes index of tetras vertices */
    // TODO :: Maybe add flag into tetraVertices_ToSend to avoid 2 communications
    PMMG_CALLOC(parmesh,tetraVertices_ToRecv_outIdx,4*ntTot_out2in,int,"tetraVertices_ToRecv_outIdx",ier = 0);
    PMMG_CALLOC(parmesh,tetraVertices_ToRecv_inIdx, 4*ntTot_out2in,int,"tetraVertices_ToRecv_inIdx",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(tetraVertices_ToSend,       4*ntTot_in2out,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   tetraVertices_ToRecv_outIdx,4*ntTot_out2in,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /* Send/receive nodes flag of tetras vertices */
    PMMG_CALLOC(parmesh,tetraFlag_ToRecv,4*ntTot_out2in,int,"tetraFlag_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(tetraFlag_ToSend,4*ntTot_in2out,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   tetraFlag_ToRecv,4*ntTot_out2in,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /* Send/receive nodes tag */
    // TODO :: group the next 2 communications
    PMMG_CALLOC(parmesh,pointTagInt_ToRecv,npInt_out2in,uint16_t,"pointTagInt_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointTagInt_ToSend,npInt_in2out,MPI_UINT16_T,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   pointTagInt_ToSend,npInt_out2in,MPI_UINT16_T,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    PMMG_CALLOC(parmesh,pointTagPBDY_ToRecv,npPBDY_out2in,uint16_t,"pointTagPBDY_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointTagPBDY_ToSend,npPBDY_in2out,MPI_UINT16_T,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   pointTagPBDY_ToRecv,npPBDY_out2in,MPI_UINT16_T,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /* Send/receive nodes references */
    // TODO :: group the next 2 communications
    PMMG_CALLOC(parmesh,pointRefInt_ToRecv,npInt_out2in,int,"pointRefInt_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointRefInt_ToSend,npInt_in2out,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   pointRefInt_ToRecv,npInt_out2in,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    PMMG_CALLOC(parmesh,pointRefPBDY_ToRecv,npPBDY_out2in,int,"pointRefPBDY_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointRefPBDY_ToSend,npPBDY_in2out,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   pointRefPBDY_ToRecv,npPBDY_out2in,MPI_INT,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /* Send/receive nodes coordinates */
    // TODO :: group the next 2 communications
    PMMG_CALLOC(parmesh,pointCoordInt_ToRecv,3*npInt_out2in,double,"pointCoordInt_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointCoordInt_ToSend,3*npInt_in2out,MPI_DOUBLE,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   pointCoordInt_ToRecv,3*npInt_out2in,MPI_DOUBLE,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    PMMG_CALLOC(parmesh,pointCoordPBDY_ToRecv,3*npPBDY_out2in,double,"pointCoordPBDY_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointCoordPBDY_ToSend,3*npPBDY_in2out,MPI_DOUBLE,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   pointCoordPBDY_ToRecv,3*npPBDY_out2in,MPI_DOUBLE,color_out,MPI_OVERLAP_TAG+LOCAL_MPI_TAG,
                   comm,&status),return 0 );
    LOCAL_MPI_TAG++;

    /** STEP 6 - Fill the overlap hash table with interface points **/
    for (i=0; i < nitem_ext; i++) {
      ip_out = pointIdxInterface_ToRecv[i]; // Node index on Pcolor_out (the other partition)
      ip_in  = pointIdxInterface_ToSend[i]; // Node index on Pcolor_in  (this partition)
      hash_out2in[ip_out] = ip_in;  // From index on Pcolor_out, I found index on Pcolor_in
      hash_in2out[ip_in]  = ip_out; // From index on Pcolor_in,  I found index on Pcolor_out
    }

    /** STEP 7 - Add nodes having MG_PARBDY tag to the overlap hash tables **/
    mesh->xpmax  = MG_MAX( (long long)(1.5*mesh->xp),mesh->npmax);

    /* Loop over the points having MG_PARBDY tag */
    for (i=0; i < npPBDY_out2in+1; i++) {
      ip_out=pointIdxPBDY_ToRecv[i];
      coord[0]=pointCoordPBDY_ToRecv[3*i];
      coord[1]=pointCoordPBDY_ToRecv[3*i+1];
      coord[2]=pointCoordPBDY_ToRecv[3*i+2];
      tag     =pointTagPBDY_ToRecv[i]+MG_OVERLAP;
      ref     =pointRefPBDY_ToRecv[i];

      /* Loop over the data allowing to identify these points */
      for (j=0; j < ndataPBDY_out2in+1; j++) {
        color_ter = dataPBDY_ToRecv[3*j];
        loc       = dataPBDY_ToRecv[3*j+1];
        ip        = dataPBDY_ToRecv[3*j+2];
        min_color = MG_MIN(color_out,color_ter);
        max_color = MG_MAX(color_out,color_ter);

        /* Search if this point has already been added to this mesh */
        duplicated_point = 0;
        if (ip == ip_out) {
          /* Search into dataPBDY_AlreadyAdded to seen if it exists already */
          for (r=0; r < ndataPBDY_added+1; r++) {
            if ( (color_out == dataPBDY_AlreadyAdded[5*r+1]) && (color_ter == dataPBDY_AlreadyAdded[5*r])) {
              if ((loc == dataPBDY_AlreadyAdded[5*r+2]) ) {
                ip_in = dataPBDY_AlreadyAdded[5*r+3];
                duplicated_point =1;
                break;
              }
            }
            if ( (color_out == dataPBDY_AlreadyAdded[5*r]) && (ip_out == dataPBDY_AlreadyAdded[5*r+4])) {
                ip_in = dataPBDY_AlreadyAdded[5*r+3];
                duplicated_point =1;
                break;
            }
          }

          /* If the point has not been added, yet, then add it */
          if (duplicated_point==0) {
            if ( color_ter != color_in ) {
#ifdef USE_POINTMAP
              src = 1; //p0->src;
#else
              src = 1;
#endif
              ip_in = MMG3D_newPt(mesh,coord,tag,src);
              // TODO :: add ref here
            }
          }
          hash_out2in[ip_out] = ip_in;  // From index on color_out, I found index on color_in
          hash_in2out[ip_in]  = ip_out; // From index on color_in,  I found index on color_out
          dataPBDY_AlreadyAdded[5*ndataPBDY_added]  = color_out;
          dataPBDY_AlreadyAdded[5*ndataPBDY_added+1]= color_ter;
          dataPBDY_AlreadyAdded[5*ndataPBDY_added+2]= loc;
          dataPBDY_AlreadyAdded[5*ndataPBDY_added+3]= ip_in;
          dataPBDY_AlreadyAdded[5*ndataPBDY_added+4]= ip_out;
          ndataPBDY_added++;
        }
      }
    }

    if (color_in == 0) {
      fprintf(stdout, "OVERLAP Proc %d :: ndataPBDY_added=%d  \n", color_in,ndataPBDY_added);
      for (i=0; i < ndataPBDY_added; i++) {
        fprintf(stdout, "OVERLAP Proc=%d :: dataPBDY_AlreadyAdded=[%d-%d-%d-%d-%d]\n",
                          color_in,
                          dataPBDY_AlreadyAdded[5*i],dataPBDY_AlreadyAdded[5*i+1],dataPBDY_AlreadyAdded[5*i+2],
                          dataPBDY_AlreadyAdded[5*i+3],dataPBDY_AlreadyAdded[5*i+4]);
      }
    }

    /** STEP 8 - Add all the other nodes to the overlap hash tables **/
    icoord=0;

    for (i=0; i < 4*ntTot_out2in; i++) {
      ip_out=tetraVertices_ToRecv_outIdx[i];

      if (tetraFlag_ToRecv[i]==1) {
        coord[0]=pointCoordInt_ToRecv[3*icoord];
        coord[1]=pointCoordInt_ToRecv[3*icoord+1];
        coord[2]=pointCoordInt_ToRecv[3*icoord+2];
        ref     =pointRefInt_ToRecv[icoord];
        tag     =pointTagInt_ToRecv[icoord]+MG_OVERLAP;
        icoord += 1;

#ifdef USE_POINTMAP
        src = 1; //p0->src;
#else
        src = 1;
#endif

        ip_in = MMG3D_newPt(mesh,coord,tag,src);
        mesh->point[ip_in].ref = ref;
        hash_out2in[ip_out] = ip_in;  // From index on color_out, I found index on color_in
        hash_in2out[ip_in]  = ip_out; // From index on color_in,  I found index on color_out
        tetraVertices_ToRecv_inIdx[i]=ip_in;
      }
      else{
        tetraVertices_ToRecv_inIdx[i]=hash_out2in[ip_out];
      }
    }


    if (print_msg) {
      for (i=0; i < 4*ntTot_out2in; i+=4) {
        fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: inIdx=[%d-%d-%d-%d], outIdx=[%d-%d-%d-%d]\n",
                          color_in,color_out,
                          tetraVertices_ToRecv_inIdx[i],tetraVertices_ToRecv_inIdx[i+1],tetraVertices_ToRecv_inIdx[i+2],tetraVertices_ToRecv_inIdx[i+3],
                          tetraVertices_ToRecv_outIdx[i],tetraVertices_ToRecv_outIdx[i+1],tetraVertices_ToRecv_outIdx[i+2],tetraVertices_ToRecv_outIdx[i+3]);
      }

      fprintf(stdout, "\n\n");

      for (i=0; i < np_in+npTot_out2in+1; i++) {
        if (hash_in2out[i] != 0) {
          fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: hash_in2out - ip_in=%d, ip_out=%d\n",
                            color_in,color_out,i,hash_in2out[i]);
        }
      }

      fprintf(stdout, "\n\n");

      for (i=0; i < np_out+1; i++) {
        if (hash_out2in[i] != 0) {
          fprintf(stdout, "OVERLAP Proc=%d<-Proc=%d :: hash_out2in ip_out=%d, ip_in=%d\n",
                              color_in,color_out,i,hash_out2in[i]);
        }
      }
    }

    /* Add the tetra to the mesh */
    int v0,v1,v2,v3;
    MMG5_pTetra         ptnew;
    MMG5_int iel;
    // MMG5_xTetra         xt,xt1;

    // pt0=&mesh->tetra[1];

    for (i=0; i < ntTot_out2in; i++) {
      v0 = tetraVertices_ToRecv_inIdx[4*i];
      v1 = tetraVertices_ToRecv_inIdx[4*i+1];
      v2 = tetraVertices_ToRecv_inIdx[4*i+2];
      v3 = tetraVertices_ToRecv_inIdx[4*i+3];

      iel = MMG3D_newElt(mesh);
      if ( !iel ) {
        MMG3D_TETRA_REALLOC(mesh,iel,mesh->gap,
                            fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                    " a new element.\n",__func__);
                            MMG5_INCREASE_MEM_MESSAGE();
                            fprintf(stderr,"  Exit program.\n");
                            return 0);
      }
      ptnew = &mesh->tetra[iel];
      // memcpy(ptnew,pt0,sizeof(MMG5_Tetra));
      ptnew->v[0] = v0;
      ptnew->v[1] = v1;
      ptnew->v[2] = v2;
      ptnew->v[3] = v3;
      ptnew->qual = MMG5_caltet(mesh,met,ptnew) ;//MMG5_orcal(mesh,met,iel);
      ptnew->tag |= MG_OVERLAP;
    }

    /* Deallocate memory*/
    PMMG_DEL_MEM(parmesh,pointCoordInt_ToSend,double,"pointCoordInt_ToSend");
    PMMG_DEL_MEM(parmesh,pointCoordInt_ToRecv,double,"pointCoordInt_ToRecv");

    PMMG_DEL_MEM(parmesh,pointCoordPBDY_ToSend,double,"pointCoordPBDY_ToSend");
    PMMG_DEL_MEM(parmesh,pointCoordPBDY_ToRecv,double,"pointCoordPBDY_ToRecv");

    PMMG_DEL_MEM(parmesh,pointTagInt_ToSend,uint16_t,"pointTagInt_ToSend");
    PMMG_DEL_MEM(parmesh,pointTagInt_ToRecv,uint16_t,"pointTagInt_ToRecv");

    PMMG_DEL_MEM(parmesh,pointTagPBDY_ToSend,uint16_t,"pointTagPBDY_ToSend");
    PMMG_DEL_MEM(parmesh,pointTagPBDY_ToRecv,uint16_t,"pointTagPBDY_ToRecv");

    PMMG_DEL_MEM(parmesh,pointIdxInterface_ToSend,int,"pointIdxInterface_ToSend");
    PMMG_DEL_MEM(parmesh,pointIdxInterface_ToRecv,int,"pointIdxInterface_ToRecv");

    PMMG_DEL_MEM(parmesh,tetraVertices_ToSend,int,"tetraVertices_ToSend");
    PMMG_DEL_MEM(parmesh,tetraVertices_ToRecv_inIdx, int,"tetraVertices_ToRecv_inIdx");
    PMMG_DEL_MEM(parmesh,tetraVertices_ToRecv_outIdx,int,"tetraVertices_ToRecv_outIdx");

    PMMG_DEL_MEM(parmesh,tetraFlag_ToSend,int,"tetraFlag_ToSend");
    PMMG_DEL_MEM(parmesh,tetraFlag_ToRecv,int,"tetraFlag_ToRecv");

    PMMG_DEL_MEM(parmesh,dataPBDY_ToSend,int,"dataPBDY_ToSend");
    PMMG_DEL_MEM(parmesh,dataPBDY_ToRecv,int,"dataPBDY_ToRecv");

    PMMG_DEL_MEM(parmesh,pointIdxPBDY_ToSend,int,"pointIdxPBDY_ToSend");

  }

  /* Deallocate memory*/
  PMMG_DEL_MEM(parmesh,dataPBDY_AlreadyAdded,int,"dataPBDY_AlreadyAdded");

  /* Realloc np and ne */
  // mesh->ne = mesh->nei;
  // mesh->np = mesh->npi;

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout, "\n\n-------> END of OVERLAP \n");

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param comm pointer toward ...
 *
 * \return 1 if success, 0 if fail.
 *
 * \todo Fill the funtion
 *
 * TODO
 *
 */
int PMMG_delete_overlap(PMMG_pParMesh parmesh, MPI_Comm comm) {

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
      fprintf(stdout,"\n      ## TODO:: PMMG_delete_overlap.\n");

return 1;
}
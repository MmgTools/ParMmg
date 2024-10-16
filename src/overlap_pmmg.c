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
 * \file overlap_pmmg.c
 * \brief Create and delete overlap.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (InriaSoft)
 * \author Laetitia Mottet (UBordeaux)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Functions to create and delete the overlap
 *
 */

#include "parmmg.h"
#include "mmgexterns_private.h"
#include "inlined_functions_3d_private.h"

/**
 * \param parmesh pointer toward a parmesh structure
 * \param comm MPI communicator for ParMmg
 *
 * \return 1 if success, 0 if fail.
 *
 * Create the overlap. The overlap consists in sending to and receiving from
 * neighbour partitions one extra layer of point and associated tetra.
 *
 * \remark Data transferred between partitions:
 *  - mesh->point.c, mesh->point.tag and mesh->point.ref
 *  - mesh->tetra.v and mesh->tetra.ref
 *  - mesh->ls
 * Data NOT transferred between partitions:
 *  - Other mesh->point and mesh->tetra fields
 *  - mesh->xtetra fields
 *
 */
int PMMG_create_overlap(PMMG_pParMesh parmesh, MPI_Comm comm) {

  /* Local variables
     Remark: *Interior_*: pts not tagged MG_PARBDY
             *PBDY_*    : pts     tagged MG_PARBDY
             *in2out: data from color_in sends to color_out
             *out2in: data on color_in receives from color_out */
  PMMG_pInt_comm int_comm;
  PMMG_pExt_comm ext_comm,ext_comm_ter;
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_pSol      ls;
  PMMG_pOverlap  overlap;
  MPI_Status     status;
  MMG5_pPoint    p0;
  MMG5_pTetra    pt;

  int ier = 1;       // Initialize error
  int local_mpi_tag; // Tag used for MPI comm
  int i,j,k,r;
  int idx, ip, ip_in,ip_out,ip_ter,pos;
  int icomm,icomm_ter;
  int ref;
  int duplicated_point;
  int nitem_ext,nitem_ext_ter,next_comm;
  int color_in,color_out,color_ter;
  double coord[3],ls_val;
  uint16_t tag;

  int *n_ToSend, *n_ToRecv;       // Tables to send/receive nbr of points and tetras
  int *hash_in2out, *hash_out2in; // Hash table needed in structure PMMG_pOverlap
  double *lsInterior_ToSend, *lsPBDY_ToSend, *lsInterior_ToRecv, *lsPBDY_ToRecv; // LS values to send/receive

  /* Number of tetras */
  int nt_initial;   // Initial tetras nbr on   Pcolor_in
  int ntTot_in2out; // Total   tetras nbr from Pcolor_in sends    to   Pcolor_out
  int ntTot_out2in; // Total   tetras nbr on   Pcolor_in receives from Pcolor_out

  /* Tetras vertices index and ref to send/receive */
  int *tetraVertices_ToSend;        // Indices of tetra vertices from Pcolor_in sends to Pcolor_out
  int *tetraVertices_ToRecv_inIdx;  // Indices of tetra vertices from Pcolor_out on Pcolor_in
  int *tetraVertices_ToRecv_outIdx; // Indices of tetra vertices from Pcolor_out on Pcolor_out
  int *tetraRef_ToSend, *tetraRef_ToRecv;
  int *tetraVerticesSeen_ToSend, *tetraVerticesSeen_ToRecv; // Flag tetra vertices

  /* Number of points */
  int np_in, np_out;     // Nbr of pts on Pcolor_in or Pcolor_out
  int npInterior_in2out; // Nbr of pts not tagged MG_PARBDY from Pcolor_in  sends    to Pcolor_out
  int npPBDY_in2out;     // Nbr of pts     tagged MG_PARBDY from Pcolor_in  sends    to Pcolor_out
  int npInterior_out2in; // Nbr of pts not tagged MG_PARBDY from Pcolor_out receives on Pcolor_in
  int npPBDY_out2in;     // Nbr of pts     tagged MG_PARBDY from Pcolor_out receives on Pcolor_in
  int npTot_in2out;      // Total nbr of pts from Pcolor_in  sends    to Pcolor_out npTot=npInterior+npPBDY
  int npTot_out2in;      // Total nbr of pts from Pcolor_out receives on Pcolor_in  npTot=npInterior+npPBDY

  /* Points coordinates, tag, index and ref to send/receive */
  double   *pointCoordInterior_ToSend, *pointCoordPBDY_ToSend;
  double   *pointCoordInterior_ToRecv, *pointCoordPBDY_ToRecv;
  uint16_t *pointTagInterior_ToSend, *pointTagPBDY_ToSend;
  uint16_t *pointTagInterior_ToRecv, *pointTagPBDY_ToRecv;
  int      *pointRefInterior_ToSend, *pointRefPBDY_ToSend;
  int      *pointRefInterior_ToRecv, *pointRefPBDY_ToRecv;
  int      *pointIdxPBDY_ToSend, *pointIdxInterface_ToSend;
  int      *pointIdxPBDY_ToRecv, *pointIdxInterface_ToRecv;

  /* Data needed to identify MG_PARBDY pts located on Pcolor_ter
     and ensure they are added only once in Pcolor_in */
  int ndataPBDY_in2out, ndataPBDY_out2in, ndataPBDY_added;        // Nbr of MG_PARBDY points
  int *dataPBDY_ToSend, *dataPBDY_ToRecv, *dataPBDY_AlreadyAdded; // Data to identify MG_PARBDY points

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
      fprintf(stdout,"\n      ## Create Overlap.\n");

  /* Creation of overlap works only when there is one group */
  /* Ensure only one group on each proc */
  assert (  (parmesh->ngrp == 1 || parmesh->ngrp == 0)
            && "Overlap not implemented for more than 1 group per rank");

  /* Global initialization */
  grp  = &parmesh->listgrp[0];
  mesh =  parmesh->listgrp[0].mesh;
  ls   =  parmesh->listgrp[0].ls;
  next_comm = parmesh->next_node_comm; // Nbr of external node communicators
  int_comm  = parmesh->int_node_comm;  // Internal node communicator
  nt_initial      = mesh->ne;
  ndataPBDY_added = 0;

  /* Reset flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* Global allocation memory */
  PMMG_CALLOC(parmesh,int_comm->intvalues,int_comm->nitem,int,"intvalues",return 0);
  PMMG_CALLOC(parmesh,parmesh->overlap,next_comm,PMMG_Overlap,"overlap",ier = 0);
  PMMG_CALLOC(parmesh,dataPBDY_AlreadyAdded,5*mesh->np,int,"dataPBDY_AlreadyAdded",ier = 0);

  PMMG_CALLOC(parmesh,n_ToSend,6,int,"n_ToSend",ier = 0);
  PMMG_CALLOC(parmesh,pointCoordInterior_ToSend, 3*mesh->np, double,  "pointCoordInterior_ToSend",ier = 0);
  PMMG_CALLOC(parmesh,pointCoordPBDY_ToSend,     3*mesh->np, double,  "pointCoordPBDY_ToSend",    ier = 0);
  PMMG_CALLOC(parmesh,pointTagInterior_ToSend,     mesh->np, uint16_t,"pointTagInterior_ToSend",  ier = 0);
  PMMG_CALLOC(parmesh,pointTagPBDY_ToSend,         mesh->np, uint16_t,"pointTagPBDY_ToSend",      ier = 0);
  PMMG_CALLOC(parmesh,pointRefInterior_ToSend,     mesh->np, int,     "pointRefInterior_ToSend",  ier = 0);
  PMMG_CALLOC(parmesh,pointRefPBDY_ToSend,         mesh->np, int,     "pointRefPBDY_ToSend",      ier = 0);
  PMMG_CALLOC(parmesh,pointIdxPBDY_ToSend,         mesh->np, int,     "pointIdxPBDY_ToSend",      ier = 0);
  PMMG_CALLOC(parmesh,pointIdxInterface_ToSend,    mesh->np, int,     "pointIdxInterface_ToSend",ier = 0);
  PMMG_CALLOC(parmesh,tetraVertices_ToSend,      4*mesh->ne, int,     "tetraVertices_ToSend",     ier = 0);
  PMMG_CALLOC(parmesh,tetraVerticesSeen_ToSend,  4*mesh->ne, int,     "tetraVerticesSeen_ToSend", ier = 0);
  PMMG_CALLOC(parmesh,tetraRef_ToSend,             mesh->ne, int,     "tetraRef_ToSend",          ier = 0);
  PMMG_CALLOC(parmesh,lsInterior_ToSend,           mesh->np, double,  "lsInterior_ToSend",        ier = 0);
  PMMG_CALLOC(parmesh,lsPBDY_ToSend,               mesh->np, double,  "lsPBDY_ToSend",            ier = 0);

  /* Store point index in internal communicator intvalues.
     Like we have only one group (see the assert at the beginning of the function),
     we can store directly the id of the node in the internal comm int_comm->intvalue
     as we know it won't be overwrite by another group. Indeed, entities shared
     by the groups have a unique position in this buffer. In short, the next steps
     work because we have only one group when we are in this function. */
  for( i = 0; i < grp->nitem_int_node_comm; i++ ){
    ip  = grp->node2int_node_comm_index1[i];
    pos = grp->node2int_node_comm_index2[i];
    int_comm->intvalues[pos] = ip;
  }

  /* Loop over the number of external node communicator */
  for (icomm=0; icomm<next_comm; icomm++) {

    /* Local initialization */
    local_mpi_tag     = 1;
    npPBDY_in2out     = 0;
    npInterior_in2out = 0;
    npTot_in2out      = 0;
    ntTot_in2out      = 0;
    ndataPBDY_in2out  = 0;
    np_in             = mesh->np;

    /* Get external node communicator information */
    ext_comm  = &parmesh->ext_node_comm[icomm]; // External node communicator
    color_in  = ext_comm->color_in;             // Color of this partition Pcolor_in
    color_out = ext_comm->color_out;            // Color of the remote partition Pcolor_out
    nitem_ext = ext_comm->nitem;                // Nbr of nodes in common between Pcolor_in and Pcolor_out

    /* Overlap variables */
    overlap = &parmesh->overlap[icomm];
    overlap->color_in  = color_in;
    overlap->color_out = color_out;

    /** STEP 1 - Identify nodes at the interface between Pcolor_in and Pcolor_out
              a) Assign flag -1 to these nodes
              b) Store the local point indices in pointIdxInterface_ToSend
                 to be able to fill hash_out2in
          Note: The following works because we have only one group in this function
          (see the assert at the beginning of the function) and we have been able
          to store the local index of node **/
    /* Loop over the nodes in the external node communicator */
    for (i=0; i < nitem_ext; i++) {
      /* Get the index of the node stored in internal communicator */
      pos = ext_comm->int_comm_index[i];
      ip  = int_comm->intvalues[pos];

      /* Add the flag -1 to this point */
      p0 = &mesh->point[ip];
      p0->flag = -1;
      pointIdxInterface_ToSend[i] = ip;
    }

    /** STEP 2 - Identification of points and tetras to send to Pcolor_out
         tetraVerticesSeen_ToSend: flag the tetra vertices using following rules
           +1 means the vertex is seen for the first time
            0 means the vertex has been already seen (from another tetra)
              OR is on the interface between Pcolor_in and Pcolor_out
              OR is tagged MG_PARBDY (special treatments for those) **/
    /* Loop over the tetra on this partition, i.e. Pcolor_in */
    for (k=1; k<=nt_initial; k++) {
      pt  = &mesh->tetra[k];
      if (!MG_EOK(pt)) continue;

      /* If one vertex if flag -1, assign MG_OVERLAP to this tetra */
      for (i=0; i<4; i++) {
        ip = pt->v[i];
        p0 = &mesh->point[ip];
        if ( p0->flag < 0 ) {
          pt->tag |= MG_OVERLAP;
          break;
        }
      }

      /* If tetra has not been identified as MG_OVERLAP, then ignore it */
      if ( !(pt->tag & MG_OVERLAP) ) continue;

      tetraRef_ToSend[ntTot_in2out] = pt->ref;

      /* Loop over the vertices of this tetra, all the nodes belong to the overlap
         that we want to send to Pcolor_out */
      for (i=0; i<4; i++) {
        ip  = pt->v[i];
        p0  = &mesh->point[ip];
        tetraVertices_ToSend[4*ntTot_in2out+i] = ip;

        /* If this node has never been seen before for the partition on Pcolor_out
           (i.e. p0->flag!=color_out+1) and is not on interface between Pcolor_in
           and Pcolor_out (p0->flag>=0):
            a) in VerticesSeen_ToSend : assign 1 if this node is not tagged MG_PARBDY;
               assign 0 otherwise (special treatment for those).
            b) store its coordinates in pointCoord*_ToSend, its tag in
               pointTag*_ToSend; its ref in pointRef*_ToSend; its level-set value
               in ls*_ToSend. `*` being either `Interior` for pts not tagged
               MG_PARBDY or `PBDY` for pts tagged MG_PARBDY */
        if ( (p0->flag>=0)  && (p0->flag!=color_out+1) ) {

          /* Update flag of this vertex to identify it has already been seen */
          p0->flag = color_out+1;

          /* If this vertex is not tagged MG_PARBDY, store needed info in *Interior* to be sent to Pcolor_out */
          if (!(p0->tag & MG_PARBDY)) {
            tetraVerticesSeen_ToSend[4*ntTot_in2out+i] = 1; // Vertex seen for the first time
            pointCoordInterior_ToSend[3*npInterior_in2out]   = p0->c[0];
            pointCoordInterior_ToSend[3*npInterior_in2out+1] = p0->c[1];
            pointCoordInterior_ToSend[3*npInterior_in2out+2] = p0->c[2];
            pointTagInterior_ToSend[npInterior_in2out]       = p0->tag;
            pointRefInterior_ToSend[npInterior_in2out]       = p0->ref;
            lsInterior_ToSend[npInterior_in2out]             = ls->m[ip];
            npInterior_in2out++;
          }
          /* If this vertex is tagged MG_PARBDY, store needed info in *PBDY* to
             be sent to Pcolor_out. pointIdxPBDY_ToSend stores MG_PARBDY points
             except the ones at interface between Pcolor_in and Pcolor_out */
          else {
            tetraVerticesSeen_ToSend[4*ntTot_in2out+i] = 0; // Vertex tagged MG_PARBDY (special treatments for those)
            pointCoordPBDY_ToSend[3*npPBDY_in2out]   = p0->c[0];
            pointCoordPBDY_ToSend[3*npPBDY_in2out+1] = p0->c[1];
            pointCoordPBDY_ToSend[3*npPBDY_in2out+2] = p0->c[2];
            pointTagPBDY_ToSend[npPBDY_in2out]       = p0->tag;
            pointRefPBDY_ToSend[npPBDY_in2out]       = p0->ref;
            pointIdxPBDY_ToSend[npPBDY_in2out]       = ip;
            lsPBDY_ToSend[npPBDY_in2out]             = ls->m[ip];
            npPBDY_in2out++;
          }
        }
        /* If this node is on the interface between Pcolor_in and Pcolor_out (p0->flag=-1),
            this node will be treated separately, assign 0 in tetraVerticesSeen_ToSend
            OR
            If this node has been seen before in another tetra (p0->flag = color_out+1),
            assign 0 in tetraVerticesSeen_ToSend */
        else {
          tetraVerticesSeen_ToSend[4*ntTot_in2out+i]  = 0;
        }
      }
      ntTot_in2out++;
      /* Remove the tag MG_OVERLAP : here it is used as a way to identify the
         tetra to be sent. Later, tetras having the tag MG_OVERLAP will actually
         be tetras received from Pcolor_out */
      pt->tag &= ~MG_OVERLAP;
    }

    /* The total number of point to send to Pcolor_out npTot_in2out is
       the nbr of point not tagged MG_PARBDY npInterior_in2out
       plus the ones tagged MG_PARBDY npPBDY_in2out */
    npTot_in2out = npInterior_in2out + npPBDY_in2out;

    /** STEP 3 - Reinitialize points flags and tag **/
    for (i=0; i < nitem_ext; i++) {
      pos = ext_comm->int_comm_index[i];
      ip  = int_comm->intvalues[pos];
      p0  = &mesh->point[ip];
      p0->flag = 0;
      if (p0->tag & MG_OVERLAP) p0->tag &=~ MG_OVERLAP;
    }

    /** STEP 4 - Special treatment for nodes with tag MG_PARBDY.
                 Identification of the points at the interface between Pcolor_in
                 and another Pcolor_ter (!=Pcolor_out) to be sent to Pcolor_out.
                 Array `dataPBDY_ToSend` (for points with MG_PARBDY) stores
                  - Pcolor_ter: the other partition w/ which this node is shared
                        (Pcolor_ter!=Pcolor_out)
                  - The position of the point in the external communicator
                  - The local index in the mesh
                NB: This part of the algo might be optimised **/

    /* Allocate the variables dataPBDY_ToSend */
    PMMG_CALLOC(parmesh,dataPBDY_ToSend,3*npPBDY_in2out*parmesh->nprocs,int,"dataPBDY_ToSend",ier = 0);

    /* Loop over the nodes with MG_PARBDY tags to store useful data  to identify
       them in dataPBDY_ToSend */
    for (i=0; i<npPBDY_in2out; i++) {
      ip = pointIdxPBDY_ToSend[i];
      p0 = &mesh->point[ip];

      /* Loop over the partitions having nodes in common w/ Pcolor_in */
      for (icomm_ter=0; icomm_ter<next_comm; icomm_ter++) {

        /* Search on the other partition than Pcolor_out */
        if (icomm_ter == icomm) continue;

        /* Get external node communicator information */
        ext_comm_ter  = &parmesh->ext_node_comm[icomm_ter]; // External node communicator
        color_ter     = ext_comm_ter->color_out;            // Color of the remote partition Pcolor_ter
        nitem_ext_ter = ext_comm_ter->nitem;                // Nbr of nodes in common between Pcolor_in and Pcolor_ter

        assert( (color_ter != color_out) && "Unexpected case: duplicated communicator?" );

        /* Loop over the nodes in the external node communicator Pcolor_ter */
        for (j=0; j < nitem_ext_ter; j++) {
          /* Get the indices of the nodes in internal communicators */
          pos    = ext_comm_ter->int_comm_index[j];
          ip_ter = int_comm->intvalues[pos];

          if ( !(ip==ip_ter) ) continue;

          /* Each time the node ip is found being shared w/ another
              partition store Pcolor_ter and the position j of this node
              in the external comm of Pcolor_ter and the local index of the point */
          dataPBDY_ToSend[3*ndataPBDY_in2out]   = color_ter;
          dataPBDY_ToSend[3*ndataPBDY_in2out+1] = j;  // position in external communicator
          dataPBDY_ToSend[3*ndataPBDY_in2out+2] = ip; // index of point on this partition
          ndataPBDY_in2out++;
          break;
        }
      }
    }

    /** STEP 5 - Store all the different sizes to exchange **/
    n_ToSend[0] = npInterior_in2out;// Nbr of interior  point from Pcolor_in to send to Pcolor_out
    n_ToSend[1] = npPBDY_in2out;    // Nbr of MG_PARBDY point from Pcolor_in to send to Pcolor_out
    n_ToSend[2] = npTot_in2out;     // Total nbr of points from Pcolor_in to send to Pcolor_out
    n_ToSend[3] = ndataPBDY_in2out; // Nbr of data for MG_PARBDY points from Pcolor_in to send to Pcolor_out
    n_ToSend[4] = np_in;            // Total nbr of points on mesh Pcolor_in
    n_ToSend[5] = ntTot_in2out;     // Total nbr of tetras from Pcolor_in to send to Pcolor_out

    /** STEP 6 - Send and Receive all the data from the other partitions **/
    // TODO :: Improve number of communications
    /* STEP 6.1 - First send/receive the different sizes to exchange */
    PMMG_CALLOC(parmesh,n_ToRecv,6,int,"n_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(n_ToSend,6,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   n_ToRecv,6,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   comm,&status),return 0 );
    local_mpi_tag++;

    npInterior_out2in = n_ToRecv[0];
    npPBDY_out2in     = n_ToRecv[1];
    npTot_out2in      = n_ToRecv[2];
    ndataPBDY_out2in  = n_ToRecv[3];
    np_out            = n_ToRecv[4];
    ntTot_out2in      = n_ToRecv[5];

    /* Fill overlap variables*/
    overlap->np_in2out = npTot_in2out; // Total pts   nbr sends    from Pcolor_in to   Pcolor_out
    overlap->np_out2in = npTot_out2in; // Total pts   nbr receives on   Pcolor_in from Pcolor_out
    overlap->nt_in2out = ntTot_in2out; // Total tetra nbr sends    from Pcolor_in to   Pcolor_out
    overlap->nt_out2in = ntTot_out2in; // Total tetra nbr receives on   Pcolor_in from Pcolor_out

    /* STEP 6.2 - Alloc hash table */
    hash_in2out = overlap->hash_in2out;
    hash_out2in = overlap->hash_out2in;
    PMMG_CALLOC(parmesh,hash_in2out,np_in +npTot_out2in+1,int,"hash_in2out",ier = 0);
    PMMG_CALLOC(parmesh,hash_out2in,np_out             +1,int,"hash_out2in",ier = 0);

    /* STEP 6.3 - Send and receive all the other data */
    /* Send/receive indices of tetras vertices */
    PMMG_CALLOC(parmesh,tetraVertices_ToRecv_outIdx,4*ntTot_out2in,int,"tetraVertices_ToRecv_outIdx",ier = 0);
    PMMG_CALLOC(parmesh,tetraVertices_ToRecv_inIdx, 4*ntTot_out2in,int,"tetraVertices_ToRecv_inIdx",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(tetraVertices_ToSend,       4*ntTot_in2out,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   tetraVertices_ToRecv_outIdx,4*ntTot_out2in,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   comm,&status),return 0 );
    local_mpi_tag++;

    /* Send/receive flag to know if tetra vertex has already been seen */
    PMMG_CALLOC(parmesh,tetraVerticesSeen_ToRecv,4*ntTot_out2in,int,"tetraVerticesSeen_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(tetraVerticesSeen_ToSend,4*ntTot_in2out,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   tetraVerticesSeen_ToRecv,4*ntTot_out2in,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   comm,&status),return 0 );
    local_mpi_tag++;

    /* Send/receive tetras refs */
    PMMG_CALLOC(parmesh,tetraRef_ToRecv,ntTot_out2in,int,"tetraRef_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(tetraRef_ToSend,ntTot_in2out,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   tetraRef_ToRecv,ntTot_out2in,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   comm,&status),return 0 );
    local_mpi_tag++;

    /* Send/receive points indices */
    PMMG_CALLOC(parmesh,pointIdxInterface_ToRecv,nitem_ext,int,"pointIdxInterface_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointIdxInterface_ToSend,nitem_ext,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   pointIdxInterface_ToRecv,nitem_ext,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   comm,&status),return 0 );
    local_mpi_tag++;

    PMMG_CALLOC(parmesh,pointIdxPBDY_ToRecv,npPBDY_out2in,int,"pointIdxPBDY_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointIdxPBDY_ToSend,npPBDY_in2out,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   pointIdxPBDY_ToRecv,npPBDY_out2in,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   comm,&status),return 0 );
    local_mpi_tag++;

    /* Send/receive points tag */
    PMMG_CALLOC(parmesh,pointTagInterior_ToRecv,npInterior_out2in,uint16_t,"pointTagInterior_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointTagInterior_ToSend,npInterior_in2out,MPI_UINT16_T,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   pointTagInterior_ToRecv,npInterior_out2in,MPI_UINT16_T,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   comm,&status),return 0 );
    local_mpi_tag++;

    PMMG_CALLOC(parmesh,pointTagPBDY_ToRecv,npPBDY_out2in,uint16_t,"pointTagPBDY_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointTagPBDY_ToSend,npPBDY_in2out,MPI_UINT16_T,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   pointTagPBDY_ToRecv,npPBDY_out2in,MPI_UINT16_T,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   comm,&status),return 0 );
    local_mpi_tag++;

    /* Send/receive points references */
    PMMG_CALLOC(parmesh,pointRefInterior_ToRecv,npInterior_out2in,int,"pointRefInterior_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointRefInterior_ToSend,npInterior_in2out,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   pointRefInterior_ToRecv,npInterior_out2in,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   comm,&status),return 0 );
    local_mpi_tag++;

    PMMG_CALLOC(parmesh,pointRefPBDY_ToRecv,npPBDY_out2in,int,"pointRefPBDY_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointRefPBDY_ToSend,npPBDY_in2out,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   pointRefPBDY_ToRecv,npPBDY_out2in,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   comm,&status),return 0 );
    local_mpi_tag++;

    /* Send/receive points coordinates */
    PMMG_CALLOC(parmesh,pointCoordInterior_ToRecv,3*npInterior_out2in,double,"pointCoordInterior_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointCoordInterior_ToSend,3*npInterior_in2out,MPI_DOUBLE,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   pointCoordInterior_ToRecv,3*npInterior_out2in,MPI_DOUBLE,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   comm,&status),return 0 );
    local_mpi_tag++;

    PMMG_CALLOC(parmesh,pointCoordPBDY_ToRecv,3*npPBDY_out2in,double,"pointCoordPBDY_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(pointCoordPBDY_ToSend,3*npPBDY_in2out,MPI_DOUBLE,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   pointCoordPBDY_ToRecv,3*npPBDY_out2in,MPI_DOUBLE,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   comm,&status),return 0 );
    local_mpi_tag++;

    /* Send/receive data needed to identify MG_PARBDY nodes */
    PMMG_CALLOC(parmesh,dataPBDY_ToRecv,3*ndataPBDY_out2in,int,"dataPBDY_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(dataPBDY_ToSend,3*ndataPBDY_in2out,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   dataPBDY_ToRecv,3*ndataPBDY_out2in,MPI_INT,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   comm,&status),return 0 );
    local_mpi_tag++;

    /* Send/receive LS values */
    PMMG_CALLOC(parmesh,lsInterior_ToRecv,npInterior_out2in,double,"lsInterior_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(lsInterior_ToSend,npInterior_in2out,MPI_DOUBLE,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   lsInterior_ToRecv,npInterior_out2in,MPI_DOUBLE,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   comm,&status),return 0 );
    local_mpi_tag++;

    PMMG_CALLOC(parmesh,lsPBDY_ToRecv,npPBDY_out2in,double,"lsPBDY_ToRecv",ier = 0);
    MPI_CHECK(
      MPI_Sendrecv(lsPBDY_ToSend,npPBDY_in2out,MPI_DOUBLE,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   lsPBDY_ToRecv,npPBDY_out2in,MPI_DOUBLE,color_out,MPI_OVERLAP_TAG+local_mpi_tag,
                   comm,&status),return 0 );
    local_mpi_tag++;

    /** STEP 7 - Fill the overlap hash tables with interface points **/
    for (i=0; i < nitem_ext; i++) {
      ip_out = pointIdxInterface_ToRecv[i]; // Point index on Pcolor_out (the other partition)
      ip_in  = pointIdxInterface_ToSend[i]; // Point index on Pcolor_in  (this partition)
      hash_out2in[ip_out] = ip_in;  // From index on Pcolor_out, I found index on Pcolor_in
      hash_in2out[ip_in]  = ip_out; // From index on Pcolor_in,  I found index on Pcolor_out
    }

    /** STEP 8 - Add nodes from Pcolor_out having MG_PARBDY tag to:
     *     a) the local mesh on Pcolor_in
     *     b) the overlap hash tables (hash_out2in and hash_in2out)
     *   Create variable dataPBDY_AlreadyAdded storing:
     *     (1) color_out: neighbouring partition
     *     (2) color_ter: tierce partition w/ which MG_PARBDY nodes is shared
     *     (3) k: position of node ip in external communicator Pcolor_out-Pcolor_ter
     *     (4) ip_in: local index of node on Pcolor_in
     *     (5) ip_out: local index of node on Pcolor_out
     *  Steps are as follows:
     *  8.1. Identify if duplicated_point = 0 or = 1 based on data in dataPBDY_AlreadyAdded
     *  8.2. If duplicated_point = 0, then we add this new point to the mesh,
     *       otherwise do nothing
     *  8.3. For both duplicated_point = 0 and = 1, the hash tables (to find
     *       correspondence between indexes) is updated with the information of this point
     *  8.4. For both duplicated_point = 0 and = 1, dataPBDY_AlreadyAdded is
     *       updated with the information of this point.
     **/
    mesh->xpmax = MG_MAX( (MMG5_int)(1.5*mesh->xp),mesh->npmax);

    /* Loop over the points having MG_PARBDY tag to add them to mesh on Pcolor_in */
    for (i=0; i < npPBDY_out2in; i++) {
      ip_out   = pointIdxPBDY_ToRecv[i]; // Index of point on Pcolor_out
      coord[0] = pointCoordPBDY_ToRecv[3*i];
      coord[1] = pointCoordPBDY_ToRecv[3*i+1];
      coord[2] = pointCoordPBDY_ToRecv[3*i+2];
      tag      = pointTagPBDY_ToRecv[i] | MG_OVERLAP; // Add the tag MG_OVERLAP to this point
      ref      = pointRefPBDY_ToRecv[i];
      ls_val   = lsPBDY_ToRecv[i];

      /* Loop over the array dataPBDY_ToRecv allowing to find the point ip_out */
      for (j=0; j < ndataPBDY_out2in; j++) {

        color_ter = dataPBDY_ToRecv[3*j];
        k         = dataPBDY_ToRecv[3*j+1]; // Position in external communicator Pcolor_out-Pcolor_ter
        ip        = dataPBDY_ToRecv[3*j+2]; // Index of point on Pcolor_out

        if ( !(ip == ip_out) ) continue;

        /* STEP 8.1 - Search into dataPBDY_AlreadyAdded if this point has already
            been added to the mesh
            - if this point exists already, then duplicated_point = 1;
            - otherwise duplicated_point = 0 */
        duplicated_point = 0;
        for (r=0; r < ndataPBDY_added+1; r++) {
          /* If the tuple (Pcolor_out-Pcolor_ter) is the right one */
          if ( (color_ter == dataPBDY_AlreadyAdded[5*r]) && (color_out == dataPBDY_AlreadyAdded[5*r+1]) ) {
            /* And the point is at the same position in external comm */
            if ( k == dataPBDY_AlreadyAdded[5*r+2] ) {
              /* then this point has alreadyh been added */
              ip_in = dataPBDY_AlreadyAdded[5*r+3];
              duplicated_point = 1;
              break;
            }
          }
          /* If this node ip_out from Pcolor_out has already been added:
             get the index ip_in it has on this partition Pcolor_in */
          if ( (color_out == dataPBDY_AlreadyAdded[5*r]) && (ip_out == dataPBDY_AlreadyAdded[5*r+4])) {
              ip_in = dataPBDY_AlreadyAdded[5*r+3];
              duplicated_point = 1;
              break;
          }
        }

        /* STEP 8.2 - If this point has not been added yet duplicated_point = 0,
            then add it to the mesh; otherwise do nothing */
        if (duplicated_point==0) {
          ip_in = MMG3D_newPt(mesh,coord,tag,1);
          mesh->point[ip_in].ref = ref; // Add ref
          mesh->point[ip_in].xp  = 0;   // Assign 0 to xp
          ls->m[ip_in]           = ls_val;
        }

        /* STEP 8.3 - Update the hash tables to know correspondence of
            local index on color_in and color_out */
        hash_out2in[ip_out] = ip_in;  // From index on color_out, I found index on color_in
        hash_in2out[ip_in]  = ip_out; // From index on color_in,  I found index on color_out

        /* STEP 8.4 - Update dataPBDY_AlreadyAdded data necessary to identify which
            MG_PARBDY points have already been treated */
        dataPBDY_AlreadyAdded[5*ndataPBDY_added]  = color_out;
        dataPBDY_AlreadyAdded[5*ndataPBDY_added+1]= color_ter;
        dataPBDY_AlreadyAdded[5*ndataPBDY_added+2]= k;
        dataPBDY_AlreadyAdded[5*ndataPBDY_added+3]= ip_in;
        dataPBDY_AlreadyAdded[5*ndataPBDY_added+4]= ip_out;
        ndataPBDY_added++;
      }
    }

    /** STEP 9 - Add all the other nodes from Pcolor_out to:
     *     a) the local mesh on Pcolor_in
     *     b) the overlap hash tables (hash_out2in and hash_in2out) **/
    j=0;

    for (i=0; i < 4*ntTot_out2in; i++) {
      ip_out=tetraVertices_ToRecv_outIdx[i]; // This point has the index ip_out on Pcolor_out

      if (tetraVerticesSeen_ToRecv[i]==1) {
        coord[0] = pointCoordInterior_ToRecv[3*j];
        coord[1] = pointCoordInterior_ToRecv[3*j+1];
        coord[2] = pointCoordInterior_ToRecv[3*j+2];
        ref      = pointRefInterior_ToRecv[j];
        tag      = pointTagInterior_ToRecv[j] | MG_OVERLAP; // Add the tag MG_OVERLAP to this point
        ls_val   = lsInterior_ToRecv[j];
        j += 1;

        /* New overlapping node is created. src is set equal to 1 by default,
           but the value does not really matter because this point will then be
           deleted by PMMG_delete_overlap once the overlap is not needed anymore */
        ip_in = MMG3D_newPt(mesh,coord,tag,1);
        mesh->point[ip_in].ref = ref; // Assign ref
        mesh->point[ip_in].xp  = 0;   // Assign 0 to xp
        ls->m[ip_in]           = ls_val;
        hash_out2in[ip_out] = ip_in;  // From index on color_out, I found index on color_in
        hash_in2out[ip_in]  = ip_out; // From index on color_in,  I found index on color_out
        tetraVertices_ToRecv_inIdx[i]=ip_in; // This point has the index ip_in on Pcolor_in
      }
      else{
        tetraVertices_ToRecv_inIdx[i]=hash_out2in[ip_out]; // Find the local index of this point from hash table
      }
    }

    /** STEP 10 - Add the tetra to the mesh */
    for (i=0; i < ntTot_out2in; i++) {
      /* Create a new tetra*/
      k = MMG3D_newElt(mesh);
      if ( !k ) {
        MMG3D_TETRA_REALLOC(mesh,k,mesh->gap,
                            fprintf(stderr,"\n  ## Error: %s: unable to allocate"
                                    " a new element.\n",__func__);
                            MMG5_INCREASE_MEM_MESSAGE();
                            fprintf(stderr,"  Exit program.\n");
                            return 0);
      }
      pt = &mesh->tetra[k];

      /* Add vertices index, ref and tag to this new tetra */
      pt->v[0] = tetraVertices_ToRecv_inIdx[4*i];
      pt->v[1] = tetraVertices_ToRecv_inIdx[4*i+1];
      pt->v[2] = tetraVertices_ToRecv_inIdx[4*i+2];
      pt->v[3] = tetraVertices_ToRecv_inIdx[4*i+3];
      pt->ref  = tetraRef_ToRecv[i];
      pt->tag |= MG_OVERLAP;
    }

    if ( parmesh->info.imprim > PMMG_VERB_VERSION )
      fprintf(stdout, "        OVERLAP - part %d sends %d pts and %d tetra to part %d\n",
                      color_in,npTot_in2out,ntTot_in2out,color_out);

    /* Reset memory*/
    memset(n_ToSend,0x00,6*sizeof(int));
    PMMG_DEL_MEM(parmesh,n_ToRecv,int,"n_ToRecv");

    memset(pointCoordInterior_ToSend,0x00,3*npInterior_in2out*sizeof(double));
    PMMG_DEL_MEM(parmesh,pointCoordInterior_ToRecv,double,"pointCoordInterior_ToRecv");

    memset(pointCoordPBDY_ToSend,0x00,3*npPBDY_in2out*sizeof(double));
    PMMG_DEL_MEM(parmesh,pointCoordPBDY_ToRecv,double,"pointCoordPBDY_ToRecv");

    memset(pointTagInterior_ToSend,0x00,npInterior_in2out*sizeof(uint16_t));
    PMMG_DEL_MEM(parmesh,pointTagInterior_ToRecv,uint16_t,"pointTagInterior_ToRecv");

    memset(pointTagPBDY_ToSend,0x00,npPBDY_in2out*sizeof(uint16_t));
    PMMG_DEL_MEM(parmesh,pointTagPBDY_ToRecv,uint16_t,"pointTagPBDY_ToRecv");

    memset(pointRefInterior_ToSend,0x00,npInterior_in2out*sizeof(int));
    PMMG_DEL_MEM(parmesh,pointRefInterior_ToRecv,int,"pointRefInterior_ToRecv");

    memset(pointRefPBDY_ToSend,0x00,npPBDY_in2out*sizeof(int));
    PMMG_DEL_MEM(parmesh,pointRefPBDY_ToRecv,int,"pointRefPBDY_ToRecv");

    memset(pointIdxInterface_ToSend,0x00,nitem_ext*sizeof(int));
    PMMG_DEL_MEM(parmesh,pointIdxInterface_ToRecv,int,"pointIdxInterface_ToRecv");

    memset(pointIdxPBDY_ToSend,0x00,npPBDY_in2out*sizeof(int));
    PMMG_DEL_MEM(parmesh,pointIdxPBDY_ToRecv,int,"pointIdxPBDY_ToRecv");

    memset(tetraVertices_ToSend,0x00,4*ntTot_in2out*sizeof(int));
    PMMG_DEL_MEM(parmesh,tetraVertices_ToRecv_inIdx, int,"tetraVertices_ToRecv_inIdx");
    PMMG_DEL_MEM(parmesh,tetraVertices_ToRecv_outIdx,int,"tetraVertices_ToRecv_outIdx");

    memset(tetraVerticesSeen_ToSend,0x00,4*ntTot_in2out*sizeof(int));
    PMMG_DEL_MEM(parmesh,tetraVerticesSeen_ToRecv,int,"tetraVerticesSeen_ToRecv");

    memset(tetraRef_ToSend,0x00,ntTot_in2out*sizeof(int));
    PMMG_DEL_MEM(parmesh,tetraRef_ToRecv,int,"tetraRef_ToRecv");

    memset(lsInterior_ToSend,0x00,npInterior_in2out*sizeof(double));
    PMMG_DEL_MEM(parmesh,lsInterior_ToRecv,double,"lsInterior_ToRecv");

    memset(lsPBDY_ToSend,0x00,npPBDY_in2out*sizeof(double));
    PMMG_DEL_MEM(parmesh,lsPBDY_ToRecv,double,"lsPBDY_ToRecv");

    PMMG_DEL_MEM(parmesh,dataPBDY_ToSend,int,"dataPBDY_ToSend");
    PMMG_DEL_MEM(parmesh,dataPBDY_ToRecv,int,"dataPBDY_ToRecv");
  }

  /* Deallocate memory*/
  PMMG_DEL_MEM(parmesh,pointCoordInterior_ToSend,double,"pointCoordInterior_ToSend");
  PMMG_DEL_MEM(parmesh,pointCoordPBDY_ToSend,double,"pointCoordPBDY_ToSend");
  PMMG_DEL_MEM(parmesh,pointTagInterior_ToSend,uint16_t,"pointTagInterior_ToSend");
  PMMG_DEL_MEM(parmesh,pointTagPBDY_ToSend,uint16_t,"pointTagPBDY_ToSend");
  PMMG_DEL_MEM(parmesh,pointRefInterior_ToSend,int,"pointRefInterior_ToSend");
  PMMG_DEL_MEM(parmesh,pointRefPBDY_ToSend,int,"pointRefPBDY_ToSend");
  PMMG_DEL_MEM(parmesh,pointIdxInterface_ToSend,int,"pointIdxInterface_ToSend");
  PMMG_DEL_MEM(parmesh,pointIdxPBDY_ToSend,int,"pointIdxPBDY_ToSend");
  PMMG_DEL_MEM(parmesh,tetraVertices_ToSend,int,"tetraVertices_ToSend");
  PMMG_DEL_MEM(parmesh,tetraVerticesSeen_ToSend,int,"tetraVerticesSeen_ToSend");
  PMMG_DEL_MEM(parmesh,tetraRef_ToSend,int,"tetraRef_ToSend");
  PMMG_DEL_MEM(parmesh,lsInterior_ToSend,double,"lsInterior_ToSend");
  PMMG_DEL_MEM(parmesh,lsPBDY_ToSend,double,"lsPBDY_ToSend");

  PMMG_DEL_MEM(parmesh,dataPBDY_AlreadyAdded,int,"dataPBDY_AlreadyAdded");
  PMMG_DEL_MEM(parmesh,int_comm->intvalues,int,"intvalues");

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout, "        OVERLAP - part %d has %d pts and %d tetras after overlap creation\n",
                      color_in,mesh->np,mesh->ne);

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param comm MPI communicator for ParMmg
 *
 * \return 1 if success, 0 if fail.
 *
 * Delete the overlap points and tetras present in the mesh
 *
 */
int PMMG_delete_overlap(PMMG_pParMesh parmesh, MPI_Comm comm) {

  /* Local variables */
  MMG5_pMesh  mesh;
  MMG5_pTetra pt;
  MMG5_pPoint ppt;

  int i;

  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
      fprintf(stdout,"\n      ## Delete Overlap.\n");
  }

  /* Delete overlap works only when there is one group */
  /* Ensure only one group on each proc */
  assert (  (parmesh->ngrp == 1 || parmesh->ngrp == 0)
            && "Overlap not implemented for more than 1 group per rank");

  /* Global initialization */
  mesh =  parmesh->listgrp[0].mesh;

  /* Step 1 - Delete tetras with tag MG_OVERLAP */
  for (i=mesh->ne; i > 0; i--) {
    pt = &mesh->tetra[i];
    if ( !(pt->tag & MG_OVERLAP) ) continue;
    if ( !MMG3D_delElt(mesh,i) )   return 0;
  }

  /* Step 2 - Delete points with tag MG_OVERLAP */
  for (i=mesh->np; i > 0; i--) {
    ppt = &mesh->point[i];
    if ( !(ppt->tag & MG_OVERLAP) ) continue;
    MMG3D_delPt(mesh,i);
  }

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout, "        OVERLAP - part %d has %d pts and %d tetras after overlap deletion\n",
                      parmesh->myrank,mesh->np,mesh->ne);

  return 1;
}

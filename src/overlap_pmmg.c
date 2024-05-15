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
int PMMG_create_overlap(PMMG_pParMesh parmesh) {

  PMMG_pInt_comm int_node_comm;
  PMMG_pExt_comm ext_node_comm;
  PMMG_pGrp      grp,grp_init,grp_overlap;
  MMG5_pPoint    p0;
  MMG5_pMesh     mesh,meshOld;
  MMG5_pTetra    pt;
  MPI_Status     status;

  MMG5_int ip,ngrp,ne_init,n2inc_max,f2ifc_max,base_init;
  MMG5_int nt_overlap, np_overlap, s_init;

  int i, k;

  int inode;
  int i_commn;

  int nitem_int_node;
  int nitem_ext_node;
  int next_node_comm;

  int color_in_node,color_out_node;

  int idx_node_ext,idx_node_int,idx_node_mesh;

  double *rtosend,*rtorecv,*doublevalues;
  int    *itosend,*itorecv,*intvalues;

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
      fprintf(stdout,"\n      ## TODO:: PMMG_create_overlap.\n");

  /* For now - creationof overlap works only on packed mesh ,
      i.e. on one group */
  /* Ensure only one group on each proc */
  assert(parmesh->ngrp == 1);

  /* Initialization */
  grp_init  = &parmesh->listgrp[0];
  mesh = parmesh->listgrp[0].mesh;
  int_node_comm = parmesh->int_node_comm;

  doublevalues = int_node_comm->doublevalues;
  intvalues    = int_node_comm->intvalues;

  next_node_comm = parmesh->next_node_comm;  // Number of communicator for nodes

  /* STEP 1 - Identify nodes and tetra to communicate and MPI_Sendrecv them */

  /* Loop over the number of node communicator */
  for (i_commn=0; i_commn<next_node_comm; i_commn++) {

    /* Get external edge communicator information */
    ext_node_comm  = &parmesh->ext_node_comm[i_commn]; // External node communicator
    color_in_node  = ext_node_comm->color_in;          // Color of the hosting proc - this proc
    color_out_node = ext_node_comm->color_out;         // Color of the remote  proc - the proc to exchange with
    nitem_ext_node = ext_node_comm->nitem;             // Number of nodes in common between these 2 procs


    /* STEP 1.1 - First loop to estimate the number of nodes and tetra to communicate
                  for the CALLOC of itosend, itorecv, rtosend and rtorecv */
    nt_overlap = np_overlap = 0;
    s_init     = mesh->point[1].s;

    /* Loop over the nodes in the external edge communicator */
    for (inode=0; inode < nitem_ext_node; inode++) {

      /* Get the indices of the nodes in internal communicators */
      idx_node_ext  = ext_node_comm->int_comm_index[inode];
      idx_node_int  = grp->node2int_node_comm_index2[idx_node_ext];
      idx_node_mesh = grp->node2int_node_comm_index1[idx_node_int];

      /* Add the tag MG_OVERLAP to these nodes*/
      p0 = &mesh->point[idx_node_mesh];
      p0->flag=1;

    }

    /* Loop over number of tetra and assign MG_OVERLAP */
    for (k=1; k<=mesh->ne; k++) {
      /* Get the tetra k */
      pt  = &mesh->tetra[k];

      /* Loop over vertex of tetra and assign MG_OVERLAP to tetras */
      for (i=0; i<4; i++) {
        ip = pt->v[i];
        p0 = &mesh->point[ip];
        if ( p0->flag ) {
          pt->tag |= MG_OVERLAP;
          nt_overlap++;
          break;
        }
      }

      /* If tetra is now MG_OVERLAP, then assign MG_OVERLAP to nodes */
      if (pt->tag & MG_OVERLAP) {
        for (i=0; i<4; i++) {
          ip = pt->v[i];
          p0 = &mesh->point[ip];
          if ( !(p0->flag) && (p0->s != color_out_node+1) ) {
            p0->tag |= MG_OVERLAP;
            p0->s = color_out_node+1;
            np_overlap++;
          }
        }
      }
    }

    fprintf(stdout, "OVERLAP Proc=%d->Proc=%d :: npoint=%d, ntetra=%d \n",parmesh->myrank,color_out_node,np_overlap,nt_overlap);

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

    /* STEP 1.2 - Allocate itosend, itorecv, rtosend and rtorecv using CALLOC */

    /* STEP 1.3 - Second loop to fill node and tetras data to be communicated */

    for (k=1; k<=mesh->ne; k++) {

      /* Remove tetra tag - not now in the second loop */
      pt->tag =~ MG_OVERLAP;
    }

  /* STEP 3 -  */


  }


  // PMMG_CALLOC(parmesh,doublevalues,3*nitem_ext_node,double,"Alloc doublevalues",return 0);


    // PMMG_REALLOC(parmesh,ext_node_comm->rtosend,nitem_ext_node,3*np_overlap,double,"Realloc rtosend",return 0);
    // PMMG_REALLOC(parmesh,ext_node_comm->rtorecv,nitem_ext_node,3*np_overlap,double,"Realloc rtorecv",return 0);

    // PMMG_CALLOC(parmesh,ext_node_comm->rtosend,np_overlap,double,"Alloc rtosend",return 0);
    // PMMG_CALLOC(parmesh,ext_node_comm->rtorecv,np_overlap,double,"Alloc rtorecv",return 0);

    // itosend = ext_node_comm->itosend;
    // itorecv = ext_node_comm->itorecv;
    // rtosend = ext_node_comm->rtosend;
    // rtorecv = ext_node_comm->rtorecv;

    // MPI_CHECK(
    //   MPI_Sendrecv(itosend,nitem_ext_node,MPI_INT,color_out_node,MPI_OVERLAP_TAG,
    //                itorecv,nitem_ext_node,MPI_INT,color_out_node,MPI_OVERLAP_TAG,
    //                parmesh->info.read_comm,&status),return 0 );

  /* Create listgrp[1] */
  // PMMG_pGrp grpOld;
  // MMG5_pMesh meshOld;
  // int grpIdOld;
  // idx_t ngrp;

  // /* We are splitting group 0 */
  // grpIdOld = 0;

  // assert(parmesh->ngrp == 1);
  // grpOld = &parmesh->listgrp[grpIdOld];
  // meshOld = parmesh->listgrp[grpIdOld].mesh;

  // /* How many groups to split into */
  // ngrp = 2;

  // /* part array contains the groupID for each tetra */
  // PMMG_CALLOC(parmesh,part,meshOld->ne,idx_t,"metis buffer ", return 0);
  // for (k=0;k<meshOld->ne,k++) {
  //   part[k]=1;
  // }



  /* Loop over the number of node communicator */
  // for (i_commn=0; i_commn<next_node_comm; i_commn++) {

  //     /* Get external edge communicator information */
  //     ext_node_comm  = &parmesh->ext_node_comm[i_commn]; // External node communicator
  //     color_in_node  = ext_node_comm->color_in;          // Color of the hosting proc - this proc
  //     color_out_node = ext_node_comm->color_out;         // Color of the remote  proc - the proc to exchange with
  //     nitem_ext_node = ext_node_comm->nitem;             // Number of nodes in common between these 2 procs

  // }


  return 1;
}
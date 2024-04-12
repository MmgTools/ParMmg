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

  PMMG_pExt_comm ext_node_comm;
  PMMG_pGrp      grp;
  MMG5_pPoint    p0;
  MMG5_pMesh     mesh;
  MMG5_pTetra    pt;

  int i, k;

  int inode;
  int i_commn;

  int nitem_int_node;
  int nitem_ext_node;
  int next_node_comm;

  int color_in_node,color_out_node;

  int idx_node_ext,idx_node_int,idx_node_mesh;

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
      fprintf(stdout,"\n      ## TODO:: PMMG_create_overlap.\n");

  /* For now - creationof overlap works only on packed mesh ,
      i.e. on one group */
  /* Ensure only one group on each proc */
  assert(parmesh->ngrp == 1);

  /* Initialization */
  grp  = &parmesh->listgrp[0];
  mesh = parmesh->listgrp[0].mesh;

  /* Loop over the number of node communicator */
  for (i_commn=0; i_commn<next_node_comm; i_commn++) {

      /* Get external edge communicator information */
      ext_node_comm  = &parmesh->ext_node_comm[i_commn]; // External node communicator
      color_in_node  = ext_node_comm->color_in;          // Color of the hosting proc - this proc
      color_out_node = ext_node_comm->color_out;         // Color of the remote  proc - the proc to exchange with
      nitem_ext_node = ext_node_comm->nitem;             // Number of nodes in common between these 2 procs

    /* Loop over the nodes in the external edge communicator */
    for (inode=0; inode < nitem_ext_node; inode++) {

      /* Get the indices of the nodes in internal communicators */
      idx_node_ext  = ext_node_comm->int_comm_index[inode];
      idx_node_int  = grp->node2int_node_comm_index2[idx_node_ext];
      idx_node_mesh = grp->node2int_node_comm_index1[idx_node_int];

      /* Add the tag MG_OVERLAP to these nodes*/
      p0 = &mesh->point[idx_node_mesh];
      p0->tag |= MG_OVERLAP;

    }
  }

  /* Loop over number of tetra and assign MG_OVERLAP */
  for (k=1; k<=mesh->ne; k++) {

    /* Get the tetra k */
    pt  = &mesh->tetra[k];

    /* Loop over vertex of tetra and assign MG_OVERLAP */
    for (i=0; i<=3; i++) {
      if ( pt->v[i] & MG_OVERLAP) {
        pt->tag |= MG_OVERLAP;
        break;
      }
    }

    /* If tetra is now MG_OVERLAP, then assign all the nodes to MG_OVERLAP */
    if (pt->tag & MG_OVERLAP) {
      for (i=0; i<=3; i++) {
        pt->v[i] |= MG_OVERLAP;
      }
    }
 }

  return 1;
}
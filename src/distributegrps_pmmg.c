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
 * \file distributegrps_pmmg.c
 * \brief Group distribution on the processors
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \author Luca Cirrotola (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * \remark to not deaclock when 1 process fail, we continue the process even
 * when allocations fail. To avoid memory error in this case, we check each time
 * if the arrays are allocated. Try to remove this by dedicating one thread to
 * the error management.
 *
 */
#include <mpi.h>
#include "parmmg.h"
#include "mpitypes_pmmg.h"
#include "metis_pmmg.h"
#include "mpipack_pmmg.h"
#include "mpiunpack_pmmg.h"

/**
 * \param group pointer toward group to assign into another group structure
 *
 * \return The group \a out filled by the group \a group fields.
 *
 * Assign the group \a group into another group structure.
 *
 */
static inline
PMMG_Grp PMMG_assign_grp(PMMG_pGrp group) {
  PMMG_Grp out;

  out.mesh = group->mesh;
  out.field = group->field;
  out.met  = group->met;
  out.disp = group->disp;
  out.ls    = group->ls;

  out.nitem_int_node_comm = group->nitem_int_node_comm;
  out.nitem_int_edge_comm = group->nitem_int_edge_comm;
  out.nitem_int_face_comm = group->nitem_int_face_comm;

  out.node2int_node_comm_index1 = group->node2int_node_comm_index1;
  out.node2int_node_comm_index2 = group->node2int_node_comm_index2;
  out.edge2int_edge_comm_index1 = group->edge2int_edge_comm_index1;
  out.edge2int_edge_comm_index2 = group->edge2int_edge_comm_index2;
  out.face2int_face_comm_index1 = group->face2int_face_comm_index1;
  out.face2int_face_comm_index2 = group->face2int_face_comm_index2;

  out.flag = group->flag;

  return out;
}

/**
 * \param group pointer toward group to move into another group structure
 *
 * \return The group \a out filled by the group \a group fields.
 *
 * Move the group \a group into another group structure.
 *
 */
static inline
PMMG_Grp PMMG_move_grp(PMMG_pGrp group) {
  PMMG_Grp out;

  out = PMMG_assign_grp( group );

  group->mesh  = NULL;
  group->field = NULL;
  group->met   = NULL;
  group->disp  = NULL;
  group->ls    = NULL;

  group->nitem_int_node_comm = 0;
  group->nitem_int_edge_comm = 0;
  group->nitem_int_face_comm = 0;

  group->node2int_node_comm_index1 = NULL;
  group->node2int_node_comm_index2 = NULL;
  group->edge2int_edge_comm_index1 = NULL;
  group->edge2int_edge_comm_index2 = NULL;
  group->face2int_face_comm_index1 = NULL;
  group->face2int_face_comm_index2 = NULL;

  return out;
}

/**
 * \param parmesh pointer toward the parmesh structure.  \param grpI pointer
 * toward the group in which we want to merge.  \param grpJ pointer toward the
 * group that we want to merge.  \param grps pointer toward the new list of
 * groups \param first_idx index of the first group that has not been treated in
 * \a listgrp
 *
 * \return 0 if fail, 1 otherwise
 *
 * Update the node communicators when merging the group \a grpJ into the group
 * \a grpI.
 *
 */
static inline
int PMMG_mergeGrpJinI_nodeCommunicators( PMMG_pParMesh parmesh,PMMG_pGrp grpI,
                                         PMMG_pGrp grpJ,PMMG_pGrp grps,
                                         int first_idx ) {
  PMMG_pExt_comm ext_node_comm;
  int            nitem_int_node_commI,nitem_int_node_commJ,*intvalues;
  int           *node2int_node_commI_index1;
  int           *node2int_node_commI_index2,*node2int_node_commJ_index2;
  int            idx,k,i,ip;
  int            new_nitem;

  intvalues = parmesh->int_node_comm->intvalues;

  /** Step 1: Travel through the external communicators and mark their points as
   * used in the intvalues array (set a poisitive value even if the point is
   * shared between the merged groups (negative)) */
  for ( k=0; k<parmesh->next_node_comm; ++k ) {

    /* currently working external communicator */
    ext_node_comm = &parmesh->ext_node_comm[k];
    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx = ext_node_comm->int_comm_index[i];
      assert( (0<=idx ) && (idx<parmesh->int_node_comm->nitem) &&
              "check ext_comm index" );
      intvalues[idx] = abs(intvalues[idx]);
    }
  }

  /** Step 2: Process all the other groups to detec if a point marked as
   * interface between J in I is at the interface of another group and must be
   * keeped in the internal communicator */

  /* Groups stored in the grps array */
  for ( k=0; k<parmesh->nprocs; ++k ) {
    if ( !grps[k].node2int_node_comm_index2 ) continue;

    if ( &grps[k] == grpI )  continue;

    for ( i=0; i<grps[k].nitem_int_node_comm; ++i ) {
      idx = grps[k].node2int_node_comm_index2[i];
      intvalues[idx] = abs(intvalues[idx]);
    }
  }

  /* Groups remaining in the listgroup array */
  for ( k=first_idx; k<parmesh->ngrp; ++k ) {
    /* Empty group */
    if ( !parmesh->listgrp[k].node2int_node_comm_index2 ) continue;

    if ( &parmesh->listgrp[k] == grpJ ) continue;

    for ( i=0; i<parmesh->listgrp[k].nitem_int_node_comm; ++i ) {
      idx = parmesh->listgrp[k].node2int_node_comm_index2[i];
      intvalues[idx] = abs(intvalues[idx]);
    }
  }

  /** Step 3: Process the interface points of the group \a grpI and delete the
   * points that will be removed from the internal communicator (points shared
   * with the group \a grpJ and that doesn't belong to the external
   * communicators) */
  nitem_int_node_commJ        = grpJ->nitem_int_node_comm;
  nitem_int_node_commI        = grpI->nitem_int_node_comm;

  /* Worst case allocation */
  PMMG_REALLOC(parmesh, grpI->node2int_node_comm_index1,
               nitem_int_node_commI+nitem_int_node_commJ,
               nitem_int_node_commI, int,"node2int_node_comm_index1",
               return 0);
  PMMG_REALLOC(parmesh, grpI->node2int_node_comm_index2,
               nitem_int_node_commI+nitem_int_node_commJ,
               nitem_int_node_commI, int,"node2int_node_comm_index2",
               return 0);
  node2int_node_commI_index1  = grpI->node2int_node_comm_index1;
  node2int_node_commI_index2  = grpI->node2int_node_comm_index2;

  new_nitem = 0;
  for ( k=0; k<nitem_int_node_commI; ++k ) {
    idx = node2int_node_commI_index2[k];
    if ( intvalues[idx]<0 ) continue;

    /* Mark the position of intvalue as already treated */
    intvalues[idx] *=-1;

    node2int_node_commI_index1[new_nitem] = node2int_node_commI_index1[k];
    node2int_node_commI_index2[new_nitem] = node2int_node_commI_index2[k];
    ++new_nitem;
  }

  /** Step 4: Process the interface points of the group \a grpJ and add the
   * points used by the internal communicator in the list of interface points of
   * the group \a grpI */
  node2int_node_commJ_index2  = grpJ->node2int_node_comm_index2;

  for ( k=0; k<nitem_int_node_commJ; ++k ) {
    idx = node2int_node_commJ_index2[k];
    if ( intvalues[idx]<0 ) continue;

    /* Mark the position of intvalue as already treated */
    intvalues[idx] *=-1;

    /* recover the point index in grpI from the point in grpJ */
    ip = grpJ->node2int_node_comm_index1[k];
    assert ( MG_VOK(&grpJ->mesh->point[ip]) );

    node2int_node_commI_index1[new_nitem] = grpJ->mesh->point[ip].tmp;
    node2int_node_commI_index2[new_nitem] = node2int_node_commJ_index2[k];
    ++new_nitem;
  }
  PMMG_REALLOC(parmesh,grpI->node2int_node_comm_index1,new_nitem,
               nitem_int_node_commI+nitem_int_node_commJ,
               int,"node2int_node_comm_index1",return 0);
  PMMG_REALLOC(parmesh, grpI->node2int_node_comm_index2,new_nitem,
               nitem_int_node_commI+nitem_int_node_commJ,
               int,"node2int_node_comm_index2",return 0);
  grpI->nitem_int_node_comm = new_nitem;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param grpI pointer toward the group in which we want to merge.
 * \param grpJ pointer toward the group that we want to merge.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Update the face communicators when merging the group \a grpJ into the group
 * \a grpI.
 *
 */
static inline
int PMMG_mergeGrpJinI_faceCommunicators( PMMG_pParMesh parmesh,PMMG_pGrp grpI,
                                         PMMG_pGrp grpJ ) {
  int            nitem_int_face_commI,nitem_int_face_commJ,*intvalues;
  int           *face2int_face_commI_index1;
  int           *face2int_face_commI_index2,*face2int_face_commJ_index2;
  int            idx,k,iel;
  int            new_nitem;

  intvalues = parmesh->int_face_comm->intvalues;

  /** Step 1: Process the interface points of the group \a grpI and delete the
   * points that will be removed from the internal communicator (points shared
   * with the group \a grpJ and that doesn't belong to the external
   * communicators) */
  nitem_int_face_commJ        = grpJ->nitem_int_face_comm;
  nitem_int_face_commI        = grpI->nitem_int_face_comm;

  /* Worst case allocation */
  PMMG_REALLOC(parmesh, grpI->face2int_face_comm_index1,
               nitem_int_face_commI+nitem_int_face_commJ,
               nitem_int_face_commI, int,"face2int_face_comm_index1",
               return 0);
  PMMG_REALLOC(parmesh, grpI->face2int_face_comm_index2,
               nitem_int_face_commI+nitem_int_face_commJ,
               nitem_int_face_commI, int,"face2int_face_comm_index2",
               return 0);

  face2int_face_commI_index1  = grpI->face2int_face_comm_index1;
  face2int_face_commI_index2  = grpI->face2int_face_comm_index2;
  new_nitem = 0;
  for ( k=0; k<nitem_int_face_commI; ++k ) {
    idx = face2int_face_commI_index2[k];
    if ( intvalues[idx]<0 ) continue;

    face2int_face_commI_index1[new_nitem] = face2int_face_commI_index1[k];
    face2int_face_commI_index2[new_nitem] = face2int_face_commI_index2[k];
    ++new_nitem;
  }

  /** Step 2: Process the interface points of the group \a grpJ and add the
   * points used by the internal communicator in the list of interface points of
   * the group \a grpI */
  face2int_face_commJ_index2  = grpJ->face2int_face_comm_index2;

  for ( k=0; k<nitem_int_face_commJ; ++k ) {
    idx = face2int_face_commJ_index2[k];
    iel = intvalues[idx];

    assert ( iel );
    if ( iel<0 ) continue;

    assert ( MG_EOK(&grpI->mesh->tetra[iel/12]) );
    face2int_face_commI_index1[new_nitem] = iel;
    face2int_face_commI_index2[new_nitem] = face2int_face_commJ_index2[k];
    ++new_nitem;
  }
  PMMG_REALLOC(parmesh,grpI->face2int_face_comm_index1,new_nitem,
               nitem_int_face_commI+nitem_int_face_commJ,
               int,"face2int_face_comm_index1",return 0);
  PMMG_REALLOC(parmesh,grpI->face2int_face_comm_index2,new_nitem,
               nitem_int_face_commI+nitem_int_face_commJ,
               int,"face2int_face_comm_index2",return 0);
  grpI->nitem_int_face_comm = new_nitem;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param grpI pointer toward the group in which we want to merge
 * \param grpJ pointer toward the group that we want to merge with group \a grpI
 *
 * \return 0 if fail, 1 otherwise
 *
 * Update the communicators when merging the group \a grpJ into the group \a
 * grpI.
 *
 */
static inline
int PMMG_mergeGrpJinI_communicators(PMMG_pParMesh parmesh,PMMG_pGrp grpI,
                                    PMMG_pGrp grpJ, PMMG_pGrp grps,
                                    int first_idx) {

  if ( !PMMG_mergeGrpJinI_faceCommunicators(parmesh,grpI,grpJ) ) return 0;

  if ( !PMMG_mergeGrpJinI_nodeCommunicators(parmesh,grpI,grpJ,grps,first_idx) )
    return 0;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param grpI pointer toward the group in which we want to merge.
 * \param grpJ pointer toward the group that we want to merge.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Merge the interface points of the group \a grpJ (listed inside the
 * internal communicator) into the \a grpI group.
 *
 */
static inline
int PMMG_mergeGrpJinI_interfacePoints( PMMG_pParMesh parmesh,PMMG_pGrp grpI,
                                       PMMG_pGrp grpJ ) {
  MMG5_pMesh     meshI,meshJ;
  int            *intvalues;
  int            poi_id_int,poi_id_glo,k;

  meshI     = grpI->mesh;

  meshJ     = grpJ->mesh;

  intvalues = parmesh->int_node_comm->intvalues;

  /** Use the tmp field of points in meshes to remember the id in the merged mesh
   * of points that have already been added to the merged mesh or 0 if they
   * haven't been merged yet */
  for ( k = 1; k <= meshI->np; ++k )
    meshI->point[k].tmp = 0;

  for ( k = 1; k <= meshJ->np; ++k )
    meshJ->point[k].tmp = 0;

  /** Step 1: store the indices of the interface entities of meshI into the
   * internal communicators */
  for ( k = 0; k < grpI->nitem_int_node_comm; ++k ) {
    poi_id_int = grpI->node2int_node_comm_index1[k];
    poi_id_glo = grpI->node2int_node_comm_index2[k];
    assert(   ( 0 <= poi_id_glo )
           && ( poi_id_glo < parmesh->int_node_comm->nitem )
           && "check intvalues indices" );
    intvalues[ poi_id_glo ] = poi_id_int;
  }

  /** Step 2: add the interface points of the group grpJ into the group grpI */
  if ( !PMMG_mergeGrpJinI_interfacePoints_addGrpJ(parmesh,grpI,grpJ ) )
    return 0;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param grpI pointer toward the group in which we want to merge
 * \param grpJ pointer toward the group that we want to merge with group \a grpI
 *
 * \return 0 if fail, 1 otherwise
 *
 * Merge \a grpJ into grpI.
 *
 */
static inline
int PMMG_merge_grpJinI(PMMG_pParMesh parmesh,PMMG_pGrp grpI, PMMG_pGrp grpJ) {

  PMMG_pInt_comm int_node_comm,int_face_comm;
  int            *face2int_face_comm_index1,*face2int_face_comm_index2;
  int            k,iel;

  /** Use the internal communicators to mark the interface entities indices: a
   * null entitie means that the point belongs to another group or to an
   * external communicator and must be keeped, a negative value means that the
   * entitie is shared between twe groups and must be removed. */
  int_node_comm = parmesh->int_node_comm;
  int_face_comm = parmesh->int_face_comm;
  memset( int_node_comm->intvalues,0,int_node_comm->nitem*sizeof(int) );
  memset( int_face_comm->intvalues,0,int_face_comm->nitem*sizeof(int) );

  /** Step 1: Merge interface points of grpJ into grpI */
  if ( !PMMG_mergeGrpJinI_interfacePoints(parmesh,grpI,grpJ) ) return 0;

  /** Step 2: Merge internal points of grpJ into grpI */
  if ( !PMMG_mergeGrpJinI_internalPoints(grpI,grpJ) ) return 0;

  /** Step 3: Merge interfaces tetras of the grpJ into the grpI mesh */
  /* 1) Store the indices of the interface faces of grpI into the internal
   * face communicator */
  face2int_face_comm_index1 = grpI->face2int_face_comm_index1;
  face2int_face_comm_index2 = grpI->face2int_face_comm_index2;

  for ( k=0; k<grpI->nitem_int_face_comm; ++k ) {
    iel = face2int_face_comm_index1[k];
    assert(   ( 0 <= face2int_face_comm_index2[k] )
              && ( face2int_face_comm_index2[k] < parmesh->int_face_comm->nitem )
              && "check intvalues indices" );
    parmesh->int_face_comm->intvalues[face2int_face_comm_index2[k]] = iel;
  }

  /* 2) Add the interfaces tetra of the grpJ mesh to the meshI mesh */
  if ( !PMMG_mergeGrpJinI_interfaceTetra(parmesh,grpI,grpJ) ) return 0;

  /** Step 4: Merge internal tetras of the grpJ mesh mesh into the meshI mesh */
  if ( !PMMG_mergeGrpJinI_internalTetra(grpI,grpJ) ) return 0;

  grpI->mesh->npi = grpI->mesh->np;
  grpI->mesh->nei = grpI->mesh->ne;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param grps pointer toward the list of groups to pack.
 * \param ngrp pointer toward the number of groups (to update)
 * \param part pointer toward the metis array containing the partitions to
 * update and pack.
 * \param partSize initial size of the part array.
 *
 * \return 0 if fail, 1 if success.
 *
 * Pack the group array \a grps filled such as grps[k] must be send to the proc
 * k and update the parmetis partition array.
 *
 */
static inline
int PMMG_pack_grpsAndPart( PMMG_pParMesh parmesh,PMMG_pGrp *grps,
                            int *ngrp,idx_t **part,int partSiz ) {
  int k,nbl;

  nbl   = 0;
  for ( k=0; k<*ngrp; ++k ) {
    if ( !(*grps)[k].mesh ) continue;

    if ( k!=nbl )
      (*grps)[nbl] = PMMG_assign_grp( &(*grps)[k] );

    (*part)[nbl] = k;

    ++nbl;
  }
  /* Here we should never fail (nbl <= ngrp) so we can ignore the fact that if
   * the part array realloc fail, we have different sizez for the grps and the
   * part array and the grps size is not stored. */
  PMMG_REALLOC( parmesh,*grps,nbl,*ngrp,PMMG_Grp,"Groups to send", return 0 );
  PMMG_REALLOC( parmesh,*part,nbl,partSiz,idx_t,"parmetis partition",return 0 );

  *ngrp = nbl;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param grps pointer toward the list of groups to pack.
 *
 * \return 1 if success
 *
 * Pack the group array \a grps. A group is unused if its mesh is NULL.
 *
 */
static inline
int PMMG_pack_grps( PMMG_pParMesh parmesh,PMMG_pGrp *grps ) {
  int k,nbl;

  nbl   = 0;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    if ( !(*grps)[k].mesh ) continue;

    if ( k!=nbl )
      (*grps)[nbl] = PMMG_move_grp( &(*grps)[k] );

    ++nbl;
  }
  PMMG_REALLOC( parmesh,*grps,nbl,parmesh->ngrp,PMMG_Grp,"grps array", return 0 );

  parmesh->ngrp = nbl;

  return 1;
}

/**
 * \param parmesh pointer toward the mesh structure.
 *
 * \return 0 if fail, 1 if success.
 *
 * Remove the unused positions in the internal node communicator and update
 * the other communicators (external communicator, group communicators).
 *
 */
int PMMG_pack_nodeCommunicators(PMMG_pParMesh parmesh) {
  PMMG_pGrp      grp;
  PMMG_pInt_comm int_node_comm;
  PMMG_pExt_comm ext_node_comm;
  int            *intvalues;
  int            *node2int_node_comm_index2;
  int            k,nitem_int,nitem_ext,idx_int,idx_ext,i;

  intvalues = parmesh->int_node_comm->intvalues;

  /** Step 1: initialization of the communicator */
  int_node_comm = parmesh->int_node_comm;
  for ( k=0; k<int_node_comm->nitem; ++k )
    intvalues[k] = PMMG_UNSET;

  /** Step 2: Process the interface points of each group and store in intvalue
   * their new position in the internal communicator */
  nitem_int = 0;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp = &parmesh->listgrp[k];

    for ( i=0; i<grp->nitem_int_node_comm; ++i )
      if ( intvalues[grp->node2int_node_comm_index2[i]]<0 )
        intvalues[grp->node2int_node_comm_index2[i]] = nitem_int++;
  }

  /** Step 3: Update the group communicators arrays with the new points
   * position */
  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp                       = &parmesh->listgrp[k];
    node2int_node_comm_index2 = grp->node2int_node_comm_index2;

    for ( i=0; i<grp->nitem_int_node_comm; ++i ) {
      idx_int = intvalues[node2int_node_comm_index2[i]];
      assert ( idx_int>=0 && "Deleted position in internal communicator" );

      node2int_node_comm_index2[i] = idx_int;
    }
  }

  /** Step 4: Update the external communicators arrays with the new points
   * position */
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];

    nitem_ext = 0;
    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx_ext = ext_node_comm->int_comm_index[i];
      if ( idx_ext<0 ) continue;

      idx_int = intvalues[idx_ext];
      assert ( idx_int>=0 && "Deleted position in internal communicator" );

      ext_node_comm->int_comm_index[nitem_ext++] = idx_int;
    }

    if ( !PMMG_resize_extComm(parmesh,ext_node_comm,nitem_ext,&ext_node_comm->nitem) )
      return 0;
  }

  /** Step 5: Remove the empty external communicators */
  i = 0;
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];

    if ( !ext_node_comm->nitem ) continue;

    if ( i!=k ) {
      parmesh->ext_node_comm[i].nitem          = ext_node_comm->nitem;
      parmesh->ext_node_comm[i].color_in       = ext_node_comm->color_in;
      parmesh->ext_node_comm[i].color_out      = ext_node_comm->color_out;
      parmesh->ext_node_comm[i].int_comm_index = ext_node_comm->int_comm_index;
    }
    ++i;
  }
  if ( !PMMG_resize_extCommArray(parmesh,&parmesh->ext_node_comm,i,
                                 &parmesh->next_node_comm) )
    return 0;

  /** Step 6: unallocate intvalues array and set the nitem field of the internal
   * communicator to the suitable value */
  PMMG_DEL_MEM( parmesh,int_node_comm->intvalues,int,"node communicator");
  int_node_comm->nitem = nitem_int;

  return 1;
}

/**
 * \param parmesh pointer toward the mesh structure.
 *
 * \return 0 if fail, 1 if success.
 *
 * Remove the unused positions in the internal face communicator and update
 * the other communicators (external communicator, group communicators).
 *
 */
int PMMG_pack_faceCommunicators(PMMG_pParMesh parmesh) {
  PMMG_pGrp      grp;
  PMMG_pInt_comm int_face_comm;
  PMMG_pExt_comm ext_face_comm;
  int            *intvalues;
  int            *face2int_face_comm_index2;
  int            k,nitem_int,nitem_ext,idx_int,idx_ext,i;

  /** Step 1: initialization of the communicator */
  int_face_comm = parmesh->int_face_comm;

  if ( !(int_face_comm->intvalues) )
    PMMG_MALLOC(parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,"intvalues",
                 assert(0); return 0;);
  intvalues = int_face_comm->intvalues;


  for ( k=0; k<int_face_comm->nitem; ++k )
    intvalues[k] = PMMG_UNSET;

  /** Step 2: Process the interface points of each group and store in intvalue
   * their new position in the internal communicator */
  nitem_int = 0;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp = &parmesh->listgrp[k];

    for ( i=0; i<grp->nitem_int_face_comm; ++i )
      if ( intvalues[grp->face2int_face_comm_index2[i]]<0 )
        intvalues[grp->face2int_face_comm_index2[i]] = nitem_int++;
  }

  /** Step 3: Update the group communicators arrays with the new points
   * position */
  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp                       = &parmesh->listgrp[k];
    face2int_face_comm_index2 = grp->face2int_face_comm_index2;

    for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
      idx_int = intvalues[face2int_face_comm_index2[i]];
      assert ( idx_int>=0 && "Deleted position in internal communicator" );

      face2int_face_comm_index2[i] = idx_int;
    }
  }

  /** Step 4: Update the external communicators arrays with the new points
   * position and delete the unused position in ext_face_comm */
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];

    nitem_ext = 0;
    for ( i=0; i<ext_face_comm->nitem; ++i ) {

      idx_ext = ext_face_comm->int_comm_index[i];
      if ( idx_ext < 0 ) continue;

      idx_int = intvalues[idx_ext];
      assert ( idx_int>=0 && "Deleted position in internal communicator" );

      ext_face_comm->int_comm_index[nitem_ext++] = idx_int;
    }
    if ( !PMMG_resize_extComm(parmesh,ext_face_comm,nitem_ext,&ext_face_comm->nitem) )
      return 0;
  }

  /** Step 5: Remove the empty external communicators */
  i = 0;
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];

    if ( !ext_face_comm->nitem ) {
      if ( ext_face_comm->nitem_to_share ) {
        if ( ext_face_comm->itosend )
          PMMG_DEL_MEM( parmesh,ext_face_comm->itosend,int,"itosend");
        if ( ext_face_comm->itorecv )
          PMMG_DEL_MEM( parmesh,ext_face_comm->itorecv,int,"itorecv");
        if ( ext_face_comm->rtosend )
          PMMG_DEL_MEM( parmesh,ext_face_comm->rtosend,double,"rtosend");
        if ( ext_face_comm->rtorecv )
          PMMG_DEL_MEM( parmesh,ext_face_comm->rtorecv,double,"rtorecv");
        ext_face_comm->nitem_to_share = 0;
      }
      continue;
    }

    if ( i!=k ) {
      parmesh->ext_face_comm[i].nitem          = ext_face_comm->nitem;
      parmesh->ext_face_comm[i].nitem_to_share = ext_face_comm->nitem_to_share;
      parmesh->ext_face_comm[i].color_out      = ext_face_comm->color_out;
      parmesh->ext_face_comm[i].int_comm_index = ext_face_comm->int_comm_index;
      parmesh->ext_face_comm[i].itosend        = ext_face_comm->itosend;
      parmesh->ext_face_comm[i].itorecv        = ext_face_comm->itorecv;
      parmesh->ext_face_comm[i].rtosend        = ext_face_comm->rtosend;
      parmesh->ext_face_comm[i].rtorecv        = ext_face_comm->rtorecv;

      /* Avoid double free of ito*  and rto* arrays */
      ext_face_comm->nitem_to_share = 0;
    }
    ++i;
  }

  PMMG_REALLOC(parmesh,parmesh->ext_face_comm,i,parmesh->next_face_comm,
               PMMG_Ext_comm,"ext_face_comm",return 0);
  parmesh->next_face_comm = i;

  /** Step 6: unallocate intvalues array and set the nitem field of the internal
   * communicator to the suitable value */
  PMMG_DEL_MEM( parmesh,parmesh->int_face_comm->intvalues,int,"face communicator");
  int_face_comm->nitem = nitem_int;

  return 1;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param grps poiner toward list of groups partially filled
 * \param part pointer toward the metis partition array
 * \param k index of the first unmerged group (thos one that has failed).
 *
 * \return 0 if fail, the new number of groups if success.
 *
 * Try to recover a valid list of group and partition array after the failure of
 * the merge of one group into another one.
 *
 */
static inline
int PMMG_mergeGrps2Send_errorHandler( PMMG_pParMesh parmesh,PMMG_pGrp *grps,
                                      idx_t **part, int k ) {
  PMMG_pGrp listgrp;
  idx_t     *tmpPart;
  int       ngrp,nprocs,ngrpOld,j;

  nprocs  = parmesh->nprocs;
  ngrp    = parmesh->ngrp;
  listgrp = parmesh->listgrp;

  /** Error handler: Try to recover all the existing groups to continue the
   * process. */
  PMMG_RECALLOC( parmesh,*grps,(ngrp-k)+nprocs,nprocs,PMMG_Grp,"Groups to send",
                 return 0;);
  ngrp    = (ngrp-k)+nprocs;
  PMMG_CALLOC( parmesh,tmpPart,ngrp,idx_t,"New parmetis partition",return 0);

  /* At this step, the groups of listgroups that have been treated, are
   * moved/merged in grps[part[k]]: update the partition array */
  for ( j=0; j<parmesh->ngrp; ++j ) {
    if ( listgrp[j].mesh ) continue;
    assert( (*grps)[(*part)[j]].mesh );
    tmpPart[(*part)[j]] = (*part)[j];
  }
  /* Concatenation of the remaining groups of listgrp into grps */
  j = 0;
  for ( ; k<parmesh->ngrp; ++k ) {
    if ( !listgrp[k].mesh )
      continue;

    /* Find the first available location in grps array */
    while ( (*grps)[j].mesh ) ++j;
    assert ( j<ngrp );

    tmpPart[j]= (*part)[k];
    (*grps)[j++] = PMMG_move_grp(&parmesh->listgrp[k]);
  }
  /* If j<nprocs, it is possible to still have holes in the group array: remove
   * it */
  k = j;
  for ( ; j<nprocs; ++j ) {
    if ( !(*grps)[j].mesh ) continue;
    if ( k != j ) (*grps)[k] = PMMG_assign_grp(&(*grps)[j]);
    tmpPart[k] = j;
    ++k;
  }
  ngrpOld = ngrp;
  ngrp    = k;
  PMMG_REALLOC( parmesh,*grps,ngrp,ngrpOld,PMMG_Grp,"Groups to send",
                return 0;);
  PMMG_REALLOC( parmesh,tmpPart,ngrp,ngrpOld,idx_t,"parmetis partition",
                return 0;);
  PMMG_DEL_MEM( parmesh,*part,idx_t,"parmetis partition" );
  (*part) = tmpPart;

  return ngrp;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param part pointer toward the metis array containing the partitions.
 *
 * \return 0 if fail but the process can continue, -1 if the process cannot
 * continue, 1 if success.
 *
 * Merge all the groups that must be send to a given proc into 1 group.
 * Update the \a part array to match with the new groups.
 *
 * \Warning To optimize...
 *
 */
static inline
int PMMG_merge_grps2send(PMMG_pParMesh parmesh,idx_t **part) {
  PMMG_pGrp     grps,listgrp,grpI,grpJ;
  PMMG_Int_comm *int_node_comm,*int_face_comm;
  MMG5_pMesh    meshI,meshJ;
  int           nprocs,ngrp,k,j,ier;

  nprocs = parmesh->nprocs;
  ngrp   = parmesh->ngrp;

  /** Step 1: New groups allocation and initialization: move the groups to have
   * a group that will be send to proc k stored in grps[k]. Free the adja
   * array. */

  /* If this step is too expensive in memory, we can count the number of procs
     with which we will communicate and fill directly the pack array. */
  PMMG_CALLOC( parmesh,grps,nprocs,PMMG_Grp,"Groups to send",return 0 );

  j = 0;
  for ( k=0; k<ngrp; ++k ) {
    /* Free the adja array */
    meshI = parmesh->listgrp[k].mesh;
    if ( meshI->adja ) {
      PMMG_DEL_MEM(meshI, meshI->adja, int, "adjacency table" );
    }

    /* Group initialization */
    if ( grps[(*part)[k]].mesh ) continue;

    grps[(*part)[k]] = PMMG_move_grp(&parmesh->listgrp[k]);
    if ( ++j==nprocs ) break;
  }

  /** Step 2: Merge the groups that will be send to the same proc into the
   * same group */
  listgrp = parmesh->listgrp;

  // If it is too expensive in memory, this array can be allocated/unallocated
  // inside the next loop.
  int_node_comm = parmesh->int_node_comm;
  PMMG_MALLOC(parmesh,int_node_comm->intvalues,int_node_comm->nitem,int,
              "node communicator",goto strong_fail1);

  int_face_comm = parmesh->int_face_comm;
  PMMG_MALLOC(parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,
              "face communicator",goto strong_fail2);

  for ( k=0; k<ngrp; ++k ) {

    /* Merge the group that have not been moved */
    if ( !listgrp[k].mesh ) continue;

    grpI = &grps[(*part)[k]];
    grpJ = &listgrp[k];

    /* Free the adja arrays */
    meshI = grpI->mesh;
    meshJ = grpJ->mesh;
    if ( meshI->adja ) {
      PMMG_DEL_MEM(meshI, meshI->adja,int, "adjacency table" );
    }
    if ( meshJ->adja ) {
      PMMG_DEL_MEM(meshJ, meshJ->adja,int, "adjacency table" );
    }

    if ( !PMMG_merge_grpJinI(parmesh,grpI,grpJ) ) goto low_fail;

    if ( !PMMG_mergeGrpJinI_communicators(parmesh,grpI,grpJ,grps,k) ) goto low_fail;

    /* Delete the useless group to gain memory space */
    PMMG_grp_free(parmesh,&listgrp[k]);
  }

  ngrp = nprocs;
  if ( !PMMG_pack_grpsAndPart( parmesh,&grps,&ngrp,part,parmesh->ngrp ) ) ier = -1;
  else ier  = 1;
  goto end;

strong_fail2:
  PMMG_DEL_MEM( parmesh,int_node_comm->intvalues,int,"node communicator");
strong_fail1:
  PMMG_DEL_MEM( parmesh,grps,PMMG_Grp,"Groups to send");
  return -1;

low_fail:
  ngrp = PMMG_mergeGrps2Send_errorHandler(parmesh,&grps,part,k);
  if ( !ngrp ) ier=-1;
  else ier = 0;

end:

  /** Step 3: Update the parmesh */
  if( parmesh->ngrp ) {
    PMMG_DEL_MEM(parmesh,parmesh->listgrp,PMMG_Grp,"listgrp");
    parmesh->listgrp = grps;
    parmesh->ngrp    = ngrp;
  }

  /* Pack the communicators */
  if ( !PMMG_pack_nodeCommunicators(parmesh) ) ier = -1;

  if ( !PMMG_pack_faceCommunicators(parmesh) ) ier = -1;

  assert ( PMMG_check_extFaceComm(parmesh,parmesh->comm) );
  assert ( PMMG_check_extNodeComm(parmesh,parmesh->comm) );

  /* Update tag on points, tetra */
  if ( !PMMG_updateTag(parmesh) ) return -1;

  if ( !PMMG_updateMeshSize(parmesh, 1) ) {
    fprintf(stderr,"\n  ## Error: %s: Unable to update the memory repartition"
            " between meshes and communicators.\n",__func__);
    return -1;
  }

  return ier;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param recv index of the proc that receive the groups
 * \param interaction_map map of interactions with the other processors
 * \param intcomm_flag flag of the boundary faces of the sended group
 * \param nitem_intcomm_flag size of the incomm_flag array
 * \param recv_ext_idx buffer to receive data
 * \param nitem_recv_ext_idx size of recv_ext_idx buffer
 * \param ext_recv_comm external communicator \a myrank - \a recv
 * \param irequest mpi request of the send of the integer buffer
 * \param drequest mpi request of the send of the double buffer
 * \param trequest array of mpi requests of the send of the external comm
 *
 * \return 0 if fail, 1 if we success
 *
 * Transfer and update the data that are modified due to the transfer of the
 * groups from the local processor (\myrank) toward the proc \a recv.
 *
 */
static inline
int PMMG_transfer_grps_fromMetoJ(PMMG_pParMesh parmesh,const int recv,
                                 int *interaction_map,int **intcomm_flag,
                                 int *nitem_intcomm_flag,int **recv_ext_idx,
                                 int *nitem_recv_ext_idx,
                                 PMMG_pExt_comm ext_recv_comm,char **grps2send,
                                 int *pack_size,MPI_Request *irequest,
                                 MPI_Request *drequest,MPI_Request **trequest ) {

  PMMG_pGrp      grp;
  PMMG_pInt_comm int_comm;
  PMMG_pExt_comm ext_face_comm;
  MPI_Status     status;
  int            offset,nitem_recv_intcomm;
  int            k,i,count,ier,ier0,old_nitem,idx;
  int            *send2recv_int_comm,old_offset,nitem,nextcomm;
  int            nitem_ext_recv_comm;
  char           *ptr;

  const int      myrank      = parmesh->myrank;
  const int      nprocs      = parmesh->nprocs;
  const int      ngrp        = parmesh->ngrp;
  const MPI_Comm comm        = parmesh->comm;

  ier = 1;

  int_comm = parmesh->int_face_comm;

  send2recv_int_comm = NULL;

  /** Step 1: fill the send2recv_int_comm array that contains the positions (in
   * the internal communicator of the proc \a recv) of the faces of the groups
   * that are transfered */

  /* Flag the faces of the internal communicator. A face is :
   *  - 1 if it is a // boundary of a group that doesn't move;
   *  - 2 if it is at the interface between 2 groups that doesn't move;
   *  - 3 if it is a // boundary of a group that is transfered;
   *  - 4 if it is at the interface between a group that move and another that doesn't;
   *  - 6 if it is at the interface of 2 groups that are transfered.
   */
  *nitem_intcomm_flag = int_comm->nitem;
  PMMG_CALLOC(parmesh,*intcomm_flag,*nitem_intcomm_flag,int,"intcomm_flag",ier=0);
  count = 0;

  if ( *intcomm_flag ) {
    for ( k=0; k<ngrp; ++k ) {
      grp = &parmesh->listgrp[k];

      if ( recv != grp->flag ) {
        for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
          idx = grp->face2int_face_comm_index2[i];
          (*intcomm_flag)[idx] += 1;
        }
      }
      else {
        ++count;
        for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
          idx = grp->face2int_face_comm_index2[i];
          (*intcomm_flag)[idx] += 3;
        }
      }
    }
  }

  /* Store the size of the internal communicator of \a recv and the index in
   * this communicator of the faces of the external communicator myrank-recv */
  nitem = ext_recv_comm ? ext_recv_comm->nitem : 0;

  *nitem_recv_ext_idx = nitem + 1;
  PMMG_MALLOC(parmesh,*recv_ext_idx,*nitem_recv_ext_idx,int,"recv_ext_idx",ier=0);
  /* MPI_Recv will fail if recv_ext_idx is not allocated */
  MPI_CHECK( MPI_Recv(*recv_ext_idx,*nitem_recv_ext_idx,MPI_INT,recv,
                      MPI_TRANSFER_GRP_TAG+1,comm,&status), ier=0 );

  /* Fill the send2recv_int_comm array */
  offset = 0;
  nitem_recv_intcomm = (*recv_ext_idx)[offset++];

  PMMG_MALLOC(parmesh,send2recv_int_comm,int_comm->nitem,int,"send2recv_int_comm",
              ier = MG_MIN (ier,0) );

  if ( send2recv_int_comm ) {
    for ( k=0; k<int_comm->nitem; ++k ) {
      send2recv_int_comm[k] = PMMG_UNSET;
    }

    /* Process the myrank-recv external communicator and store the position of
     * this faces in the internal communicator of \a recv */
    for ( i=0; i<nitem; ++i ) {
      idx = ext_recv_comm->int_comm_index[i];
      send2recv_int_comm[idx] = (*recv_ext_idx)[offset++];
    }

    /* Process the groups and store the position (in the internal comm of \a recv)
     * of the new faces of \a recv */
    for ( k=0; k<ngrp; ++k ) {
      grp = &parmesh->listgrp[k];
      if ( recv != grp->flag ) continue;

      for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
        idx = grp->face2int_face_comm_index2[i];
        if ( send2recv_int_comm[idx] >=0 ) continue;

        send2recv_int_comm[idx] = nitem_recv_intcomm++;
      }
    }
  }

  /** Step 2: Add the new faces of the external communicator myrank-recv to the
   * external communicator myrank-recv and to the recv_ext_idx array (due to the
   * transfer of the groups from myrank toward recv) */

  /* Count the number of faces to add to the external communicator (faces flagged 4 ) */
  idx = 0;
  if ( *intcomm_flag ) {
    for ( k=0; k<int_comm->nitem; ++k ) {
      if ( 4 == (*intcomm_flag)[k] ) ++idx;
    }
  }

  /* Store in recv_ext_idx the number of groups to send, the new size of the
   * internal communicator of \a recv and the number of new faces in the
   * external communicator myrank-recv */
  PMMG_REALLOC ( parmesh,*recv_ext_idx,idx+3,*nitem_recv_ext_idx,int,"recv_ext_idx",
                 ier = 0);
  *nitem_recv_ext_idx = idx + 3;

  offset = 0;
  (*recv_ext_idx)[offset++] = count;
  (*recv_ext_idx)[offset++] = nitem_recv_intcomm;
  (*recv_ext_idx)[offset++] = idx;

  if ( idx ) {
    old_nitem = nitem;

    ier0 = 1;
    if ( !ext_recv_comm ) {
      /* Enlarge the list of external communicators */
      if ( !PMMG_resize_extCommArray ( parmesh,&parmesh->ext_face_comm,
                                       parmesh->next_face_comm+1,
                                       &parmesh->next_face_comm ) ) {
        ier0 = 0;
        ier = MG_MIN(ier,ier0);
      }

      ext_recv_comm = &parmesh->ext_face_comm[parmesh->next_face_comm-1];
      ext_recv_comm->color_in  = myrank;
      ext_recv_comm->color_out = recv;
    }

    if ( ier0 ) {
      if ( !PMMG_resize_extComm ( parmesh,ext_recv_comm,old_nitem+idx,
                                  &ext_recv_comm->nitem ) ) {
        ier0 = 0;
        ier = MG_MIN(ier,ier0);
      }

      for ( k=0; k<int_comm->nitem; ++k ) {
        if ( 4 != (*intcomm_flag)[k] ) continue;

        /* Add the new faces to the external comm myrank-recv */
        if ( ier0 )
          ext_recv_comm->int_comm_index[old_nitem++] = k;

        /* Add the new faces to the external comm to send to recv */
        (*recv_ext_idx)[offset++] = send2recv_int_comm ? send2recv_int_comm[k] : 0;
      }
    }
    else {
      /* Unable to allocate ext_recv_comm */
      for ( k=0; k<int_comm->nitem; ++k ) {
        if ( 4 != (*intcomm_flag)[k] ) continue;
        offset++;
      }
    }
  }
  assert ( offset==*nitem_recv_ext_idx );

  /** Step 3: Delete the old faces of the external communicator myrank-recv (due
   * to the transfer of the groups from myrank toward recv) from the external
   * communicator myrank-recv and add at the end of the recv_ext_idx array the
   * mark of the faces of the external communicator (0 if a face is
   * keeped, 1 otherwise). */
  old_offset = *nitem_recv_ext_idx;

  nitem_ext_recv_comm = ext_recv_comm ? ext_recv_comm->nitem : 0;
  PMMG_REALLOC ( parmesh,*recv_ext_idx,*nitem_recv_ext_idx+nitem_ext_recv_comm+1,
                 *nitem_recv_ext_idx,int,"recv_ext_idx",ier=0);
  *nitem_recv_ext_idx += nitem_ext_recv_comm + 1;

  /* Add the mark of the faces of the external comm (1 if deleted, 0 if keeped)
   * to recv_ext_idx and count the number of faces to remove from the external
   * communicator (faces flagged 3 ) */
  count = 0;

  if ( *intcomm_flag && *recv_ext_idx ) {
    for ( k=0; k<nitem_ext_recv_comm; ++k ) {
      idx = ext_recv_comm->int_comm_index[k];

      if ( 3 == (*intcomm_flag)[idx] ) {
        (*recv_ext_idx)[offset++] = 1;
        ++count;
      }
      else {
        (*recv_ext_idx)[offset++] = 0;
      }
    }
  }

  if ( *recv_ext_idx )
    (*recv_ext_idx)[offset++] = count;

  /* Remove the faces from the external communicator */
  nitem = 0;
  if ( *recv_ext_idx ) {
    for ( k=0; k<nitem_ext_recv_comm; ++k ) {
      idx = ext_recv_comm->int_comm_index[k];

      if ( (*recv_ext_idx)[old_offset+k] ) {
        /* this face must be removed */
        continue;
      }
      else {
        /* Pack the external communicator */
        ext_recv_comm->int_comm_index[nitem++] = idx;
      }
    }
    if ( nitem_ext_recv_comm ) {
      if ( !PMMG_resize_extComm ( parmesh,ext_recv_comm,nitem,
                                  &ext_recv_comm->nitem ) ) {
        ier = MG_MIN ( ier, 0 );
      }
    }
  }

  /** Step 4: When transferring the groups from myrank toward recv, the groups
   *   faces of the external communicators myrank-procB (procB != recv) becames
   *   faces of the external communicators (recv-procB):
   *
   * - remove the groups faces from the ext comm procB-myrank on procB and add
   *   it in procB-recv (send their positions in the internal comm of procB);
   *
   * - add the groups faces to the ext comm recv-procB on recv: append to the
   *   recv_ext_idx array the number of external comm to update, for each
   *   external communicator, the number of faces to add and the list of faces.
   */
  ier0 = 1;
  PMMG_MALLOC ( parmesh,*trequest,nprocs,MPI_Request,"request_tab",
                ier0 = 0; ier = MG_MIN(ier,ier0); );
  if ( ier0 ) {
    for ( k=0; k<nprocs; ++k ) {
      (*trequest)[k] = MPI_REQUEST_NULL;
    }
  }

  nextcomm = 0;
  count    = 0;

  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    if ( ext_face_comm->color_out == recv || !ext_face_comm->nitem ) continue;

    ++nextcomm;

    ext_face_comm->nitem_to_share = ext_face_comm->nitem;

    PMMG_MALLOC ( parmesh,ext_face_comm->itosend,ext_face_comm->nitem_to_share,
                  int,"itosend",ier = MG_MIN(ier,0); );

    if ( *intcomm_flag && ext_face_comm->itosend ) {
      for ( i=0; i<ext_face_comm->nitem; ++i ) {
        idx = ext_face_comm->int_comm_index[i];

        if ( 3 != (*intcomm_flag)[idx] ) {
          ext_face_comm->itosend[i] = 0;
        }
        else {
          /* We want to remove this face */
          ++count;
          ext_face_comm->itosend[i] = PMMG_UNSET;
        }
      }
    }

    if ( interaction_map[ext_face_comm->color_out] ) {
      MPI_CHECK ( MPI_Isend(ext_face_comm->itosend,ext_face_comm->nitem_to_share,
                            MPI_INT,ext_face_comm->color_out,MPI_TRANSFER_GRP_TAG+30,
                            comm,&((*trequest)[ext_face_comm->color_out])), ier=0 );
    }
  }

  PMMG_REALLOC ( parmesh,*recv_ext_idx,*nitem_recv_ext_idx+1+2*nextcomm+count,
                 *nitem_recv_ext_idx,int,"recv_ext_idx",ier=0);
  *nitem_recv_ext_idx += 1 + 2*nextcomm + count;

  if ( *recv_ext_idx )
    (*recv_ext_idx)[offset++] = nextcomm;

  old_offset = offset;
  /* Store in [old_offset:old_offset+2*nextcomm] the number of faces to send for
   * each communicator and the color_out of the external comm */
  offset += 2*nextcomm;

  nextcomm = 0;
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    if ( ext_face_comm->color_out == recv || !ext_face_comm->nitem ) continue;

    count = 0;
    nitem = 0;
    if ( *intcomm_flag && *recv_ext_idx ) {
      for ( i=0; i<ext_face_comm->nitem; ++i ) {
        idx = ext_face_comm->int_comm_index[i];

        if ( 3 == (*intcomm_flag)[idx] ) {
          (*recv_ext_idx)[offset++] = send2recv_int_comm ? send2recv_int_comm[idx] : 0;
          ++count;
        }
        else {
          /* We want to add this face to the remote communicator and to remove it
           * from the local one */
          ext_face_comm->int_comm_index[nitem++] = idx;
        }
      }
    }
    if ( !PMMG_resize_extComm(parmesh,ext_face_comm,nitem,&ext_face_comm->nitem) )
      ier = 0;

    if ( *recv_ext_idx ) {
      (*recv_ext_idx)[old_offset + nextcomm++] = count;
      (*recv_ext_idx)[old_offset + nextcomm++] = ext_face_comm->color_out;
    }
  }

  /** Step 5: send the buffer to the proc recv */
  *irequest = MPI_REQUEST_NULL;
  assert ( *nitem_recv_ext_idx == offset );
  MPI_CHECK ( MPI_Isend(*recv_ext_idx,*nitem_recv_ext_idx,MPI_INT,recv,
                        MPI_TRANSFER_GRP_TAG+3, comm,irequest), ier = 0 );

  /** Step 6: send and receive the groups */
  *pack_size = 0;
  for ( k=0; k<ngrp; ++k ) {
    grp = &parmesh->listgrp[k];
    if ( grp->flag != recv ) continue;

    /* Update the face2int_index2 array */
    for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
      idx = grp->face2int_face_comm_index2[i];
      grp->face2int_face_comm_index2[i] = send2recv_int_comm ? send2recv_int_comm[idx] : 0;
    }

    *pack_size += PMMG_mpisizeof_grp(grp);
  }

  /* Pack the groups */
  PMMG_MALLOC ( parmesh,*grps2send,*pack_size,char,"grps2send",
                ier = MG_MIN(ier,0) );

  ptr = *grps2send;
  for ( k=0; k<ngrp; ++k ) {
    grp = &parmesh->listgrp[k];

    if ( grp->flag != recv ) continue;
    PMMG_mpipack_grp(grp,&ptr);
  }

  /* Send its */
  *drequest = MPI_REQUEST_NULL;
  MPI_CHECK ( MPI_Isend ( *grps2send,*pack_size,MPI_CHAR,recv,MPI_SENDGRP_TAG,
                           comm,drequest), ier = 0 );

  /** Free the memory */
  /* Group deletion */
  for ( k=0; k<parmesh->ngrp; ++k ) {
    if ( parmesh->listgrp[k].flag == recv ) {
      PMMG_grp_free ( parmesh,&parmesh->listgrp[k] );
    }
  }

  PMMG_DEL_MEM ( parmesh,send2recv_int_comm,int,"send2recv_int_comm" );

  return ier;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param sndr index of the proc that send the groups
 * \param interaction_map map of interactions with the other processors
 * \param intcomm_flag position in the internal comm of the faces of the
 * external comm sndr-recv
 * \param nitem_intcomm_flag size of the incomm_flag array
 * \param recv_ext_idx buffer to receive data
 * \param nitem_recv_ext_idx size of recv_ext_idx buffer
 * \param ext_send_comm external communicator \a myrank - \a sndr
 *
 * \return 0 if fail, 1 if we success
 *
 * Transfer and update the data that are modified due to the transfer of the
 * groups from the processor \a sndr toward the local proc (\a myrank).
 *
 */
static inline
int PMMG_transfer_grps_fromItoMe(PMMG_pParMesh parmesh,const int sndr,
                                 int *interaction_map,int **intcomm_flag,
                                 int *nitem_intcomm_flag,int **recv_ext_idx,
                                 int *nitem_recv_ext_idx,
                                 PMMG_pExt_comm ext_send_comm,
                                 MPI_Request *irequest) {

  PMMG_pExt_comm ext_face_comm;
  MPI_Status     status;
  int            pack_size;
  int            k,ier,ier0,recv_int_nitem,offset,old_nitem;
  int            *send2recv_int_comm,nitem,nextcomm;
  int            old_offset,grpscount,idx,color_out,n,err;
  char           *buffer,*ptr;

  const int      myrank      = parmesh->myrank;
  const int      ngrp        = parmesh->ngrp;
  const MPI_Comm comm        = parmesh->comm;

#ifndef NDEBUG
  MMG5_pMesh     mesh;
#endif

  ier = 1;

  send2recv_int_comm = NULL;

  /** Step 1: build the send2recv_int_comm array that gives the new
      face2int_face_comm_index pointer for the received group */
  *nitem_intcomm_flag = ext_send_comm ? ext_send_comm->nitem + 1 : 1;
  PMMG_MALLOC ( parmesh,(*intcomm_flag),*nitem_intcomm_flag,int,"intcomm_flag",
                *nitem_intcomm_flag = 1;ier = 0 );

  if ( *intcomm_flag ) {
    (*intcomm_flag)[0] = parmesh->int_face_comm->nitem;

    if ( ext_send_comm ) {
      for ( k=0; k<ext_send_comm->nitem; ++k ) {
        (*intcomm_flag)[k+1] = ext_send_comm->int_comm_index[k];
      }
    }
  }

  /* Send the size of the internal communicator and the position of the item of
   * the external communicator in the internal communicator */
  *irequest = MPI_REQUEST_NULL;
  MPI_CHECK ( MPI_Isend ( *intcomm_flag,*nitem_intcomm_flag,MPI_INT,sndr,
                          MPI_TRANSFER_GRP_TAG+1,comm,irequest ), ier = 0 );

  /* Receive the external communicator */
  MPI_CHECK ( MPI_Probe(sndr,MPI_TRANSFER_GRP_TAG+3,comm,&status), ier = 0 );
  MPI_CHECK ( MPI_Get_count(&status,MPI_INT,nitem_recv_ext_idx), ier = 0 );

  PMMG_MALLOC ( parmesh,*recv_ext_idx,*nitem_recv_ext_idx,int,"recv_ext_idx",
                ier = 0 );

  MPI_CHECK ( MPI_Recv(*recv_ext_idx,*nitem_recv_ext_idx,MPI_INT,sndr,MPI_TRANSFER_GRP_TAG+3,comm,&status),
              ier = 0 );
  offset = 0;

  /* Get the number of groups to receive */
  grpscount = (*recv_ext_idx)[offset++];

  /* Get the size of the new internal communicator */
  recv_int_nitem = (*recv_ext_idx)[offset++];

  PMMG_MALLOC ( parmesh,parmesh->int_face_comm->intvalues,recv_int_nitem,
                int,"int_face_comm",ier = MG_MIN (ier,0) );
  parmesh->int_face_comm->nitem = recv_int_nitem;

  /** Step 2: Add the new faces of the external communicator myrank-sndr to the
   * external communicator myrank-sndr (due to the transfer of the groups from
   * sndr toward myrank) (intcomm_flag = 4 ) */
  n = (*recv_ext_idx)[offset++];

  if ( n ) {

    ier0 = 1;
    if ( !ext_send_comm ) {
      /* Enlarge external communicator */
      if ( !PMMG_resize_extCommArray ( parmesh,&parmesh->ext_face_comm,
                                       parmesh->next_face_comm+1,
                                       &parmesh->next_face_comm ) ) {
        ier0 = 0;
        ier = MG_MIN(ier,ier0);
      }
      else {
        ext_send_comm = &parmesh->ext_face_comm[parmesh->next_face_comm-1];
        ext_send_comm->color_in  = myrank;
        ext_send_comm->color_out = sndr;
      }
    }

    if ( ier0 ) {
      old_nitem = ext_send_comm->nitem;

      if ( !PMMG_resize_extComm(parmesh,ext_send_comm,n+ext_send_comm->nitem,
                                &ext_send_comm->nitem) ) {
        ier = MG_MIN(ier,ier0);
        /* Fail to alloc the external communicator sndr - myrank */
        for ( k=0; k<n; ++k ) offset++;
      }
      else {
        for ( k=0; k<n; ++k ) {
          ext_send_comm->int_comm_index[old_nitem+k] = (*recv_ext_idx)[offset++];
        }
      }
    }
    else {
      /* Fail to alloc the external communicator sndr - myrank */
      for ( k=0; k<n; ++k ) offset++;
    }
  }

  /** Step 3: Delete the old faces of the external communicator myrank-sndr (due
   * to the transfer of the groups from sndr toward myrank) from the external
   * communicator myrank-sndr (intcomm_flag = 3) */

  /* New number of faces */
  n = 0;
  if ( ext_send_comm ) {
    for ( k=0; k<ext_send_comm->nitem; ++k ) {
      idx = ext_send_comm->int_comm_index[k];

      if ( (*recv_ext_idx)[offset++] ) {
        /* Face to remove from the external communicator */
        continue;
      }
      else {
        ext_send_comm->int_comm_index[n++] = idx;
      }
    }

    assert ( ext_send_comm->nitem-n==(*recv_ext_idx)[offset] );

    if ( !PMMG_resize_extComm(parmesh,ext_send_comm,n,&ext_send_comm->nitem) ) {
      ier  = MG_MIN ( ier, 0 );
    }
  }
  else
    assert ( !(*recv_ext_idx)[offset] );

  offset++;

  /** Step 4: Update the faces of the external comms between sndr and another
   * proc than myrank (flag 3) */
  /* number of external communicators */
  nextcomm   = (*recv_ext_idx)[offset++];

  /* position of the number of faces to update in each ext_comm */
  old_offset = offset;
  /* position of the list of faces to update */
  offset    += 2*nextcomm;

  for ( n=0; n<nextcomm; ++n ) {

    nitem     = (*recv_ext_idx)[old_offset++];
    color_out = (*recv_ext_idx)[old_offset++];

    for ( k=0; k<parmesh->next_face_comm; ++k ) {
      ext_face_comm = &parmesh->ext_face_comm[k];
      if ( ext_face_comm->color_out == color_out ) break;
    }

    ier0 = 1;
    if ( k == parmesh->next_face_comm ) {
      if ( !PMMG_resize_extCommArray ( parmesh,&parmesh->ext_face_comm,
                                       parmesh->next_face_comm+1,
                                       &parmesh->next_face_comm ) ) {
        ier0 = 0;
        ier = MG_MIN(ier,0);
      }

      if ( ier0 ) {
        ext_face_comm = &parmesh->ext_face_comm[parmesh->next_face_comm-1];
        ext_face_comm->color_in  = myrank;
        ext_face_comm->color_out = color_out;
      }
    }

    if ( ier0 ) {
      old_nitem = ext_face_comm->nitem;

      if ( !PMMG_resize_extComm(parmesh,ext_face_comm,old_nitem+nitem,&ext_face_comm->nitem) ) {
        ier0 = 0;
        ier  = MG_MIN(ier,ier0);
      }
    }

    if ( ier0 ) {
      for ( idx=0; idx<nitem; ++idx ) {
        ext_face_comm->int_comm_index[old_nitem + idx] = (*recv_ext_idx)[offset++];
      }
    }
  }


  /** Step 5: Receive the new groups */
  MPI_CHECK ( MPI_Probe(sndr,MPI_SENDGRP_TAG,comm,&status), ier = 0 );
  MPI_CHECK ( MPI_Get_count(&status,MPI_CHAR,&pack_size), ier = 0 );

  PMMG_MALLOC ( parmesh,buffer,pack_size,char,"buffer", ier = 0 );

  MPI_CHECK ( MPI_Recv(buffer,pack_size,MPI_CHAR,sndr,MPI_SENDGRP_TAG,comm,&status),
              ier = 0 );

  ier0 = 1;
  if( ngrp ) {
    PMMG_RECALLOC ( parmesh,parmesh->listgrp,ngrp+grpscount,ngrp,PMMG_Grp,"listgrp",
                    ier0 = 0;ier = 0 );
  }
  else {
    if ( !parmesh->listgrp )
      PMMG_CALLOC ( parmesh,parmesh->listgrp,grpscount,PMMG_Grp,"listgrp",
                    ier0 = 0;ier = 0 );
    else
      PMMG_RECALLOC( parmesh,parmesh->listgrp,grpscount,1,PMMG_Grp,"listgrp",
                     ier0 = 0;ier = 0 );
  }

  if ( ier0 )
    parmesh->ngrp += grpscount;

    /* Count the available memory and update it while unpacking the received groups */
#ifndef NDEBUG
  for ( k=0; k<ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    assert ( mesh->npmax == mesh->np );
    assert ( mesh->xpmax == mesh->xp );
    assert ( mesh->nemax == mesh->ne );
    assert ( mesh->xtmax == mesh->xt );
  }
#endif

  ptr = buffer;
  if ( ier0 ) {
    for ( k=0; k<grpscount; ++k ) {
      err = PMMG_mpiunpack_grp(parmesh,parmesh->listgrp,ngrp+k,&ptr);
      ier = MG_MIN(ier,err);
      parmesh->listgrp[ngrp+k].flag = PMMG_UNSET;
    }

  }

  PMMG_DEL_MEM ( parmesh,buffer,char,"buffer" );
  return ier;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param sndr index of the proc that send the groups
 * \param recv index of the proc that receive the groups
 * \param interaction_map map of interactions with the other processors
 * \param called_from_distrib_mesh 1 if called for initial mesh distrib.
 * In this case do not print warnings about empty procs.
 *
 * \return 0 if fail, 1 if we success
 *
 * Transfer and update the data that are modified due to the transfer of the
 * groups from the proc \a sndr toward the proc \a recv.
 *
 */
static inline
int PMMG_transfer_grps_fromItoJ(PMMG_pParMesh parmesh,const int sndr,
                                const int recv,int *interaction_map,
                                int called_from_distrib_mesh) {

  PMMG_pExt_comm ext_face_comm,ext_send_comm,ext_recv_comm;
  MPI_Status     status;
  MPI_Request    irequest,drequest;
  MPI_Request    *trequest;
  int            k,count,ier,ier0,*recv_ext_idx,old_nitem,idx,pack_size,err;
  int            *intcomm_flag,nitem_intcomm_flag,nitem_recv_ext_idx;
  char           *grps2send;
  static int8_t  pmmgWarn = 0;

  const int      myrank      = parmesh->myrank;
  const int      nprocs      = parmesh->nprocs;

  assert ( sndr != recv );

  /* I am the receiver and I don't have any interactions with the sender */
  if ( ( recv == myrank ) && !interaction_map[sndr] ) return 1;

  /* I am the sender and I don't have any interactions with the receiver */
  if ( ( sndr == myrank ) && !interaction_map[sndr] ) return 1;

  /* I am neither the sender nor the receiver and I don't have any interactions
   * with the sender */
  if ( ( sndr != myrank ) && ( recv != myrank )  && !interaction_map[sndr] ) return 1;

  ier = 1;

  /** Step 1: find the myrank-sndr and myrank-recv external communicators */
  ext_send_comm = NULL;
  ext_recv_comm = NULL;

  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];

    if ( ext_face_comm->color_out == sndr )
      ext_send_comm = ext_face_comm;
    else if ( ext_face_comm->color_out == recv )
      ext_recv_comm = ext_face_comm;
  }


  /** Step 2: transfer and update the data */
  nitem_intcomm_flag = 0;
  nitem_recv_ext_idx = 0;
  intcomm_flag = NULL;
  recv_ext_idx = NULL;
  trequest     = NULL;

  if ( myrank == sndr ) {
    /* j = recv */
    ier = PMMG_transfer_grps_fromMetoJ(parmesh,recv,interaction_map,
                                       &intcomm_flag,&nitem_intcomm_flag,
                                       &recv_ext_idx,&nitem_recv_ext_idx,
                                       ext_recv_comm,&grps2send,&pack_size,
                                       &irequest,&drequest,&trequest);
  }
  else if ( myrank == recv ) {
    /* i = sndr */
    ier = PMMG_transfer_grps_fromItoMe(parmesh,sndr,interaction_map,
                                       &intcomm_flag,&nitem_intcomm_flag,
                                       &recv_ext_idx,&nitem_recv_ext_idx,
                                       ext_send_comm,&irequest);
  }
  else {
    /* Transfer the faces of external communicators between the sender and a
     * processor other than the receiver (intcomm_flag = 3): When proc sndr
     * send groups to proc recv, the groups faces must be moved from the
     * external communicator myrank-sndr toward the external communicator
     * myrank-recv. */
    nitem_recv_ext_idx = 0;
    if ( ext_send_comm && ext_send_comm->nitem > 0 ) {
      nitem_recv_ext_idx = ext_send_comm->nitem;
      PMMG_MALLOC(parmesh,recv_ext_idx,nitem_recv_ext_idx,int,"recv_ext_idx",
                  ier=0);
      MPI_CHECK( MPI_Recv( recv_ext_idx,nitem_recv_ext_idx,MPI_INT,sndr,
                           MPI_TRANSFER_GRP_TAG+30,parmesh->comm,&status ),
                 ier = MG_MIN(ier,0) );

      count = 0;
      for ( k=0; k<ext_send_comm->nitem; ++k ) {
        if ( PMMG_UNSET == recv_ext_idx[k] ) {
          ++count;
        }
      }

      ier0 = 1;
      if ( !ext_recv_comm ) {
        /* Enlarge the list of external communicators */
        if ( !PMMG_resize_extCommArray ( parmesh,&parmesh->ext_face_comm,
                                         parmesh->next_face_comm+1,
                                         &parmesh->next_face_comm ) ) {
          if ( !pmmgWarn ) {
            pmmgWarn = 1;
            fprintf(stderr,"  ## Error: %s: unable to resize the array of external"
                    " communicators\n",__func__);
          }
          ier0 = 0;
          ier  = MG_MIN(ier,ier0);
        }

        if ( ier0 ) {
          ext_recv_comm = &parmesh->ext_face_comm[parmesh->next_face_comm-1];
          ext_recv_comm->color_in  = myrank;
          ext_recv_comm->color_out = recv;

          /* The pointer toward ext_send_comm may be corrupted */
          for ( k=0; k<parmesh->next_face_comm; ++k ) {
            ext_face_comm = &parmesh->ext_face_comm[k];

            if ( ext_face_comm->color_out == sndr )
              ext_send_comm = ext_face_comm;
          }
        }
      }

      if ( ier0 ) {
        old_nitem = ext_recv_comm->nitem;
        if ( !PMMG_resize_extComm(parmesh,ext_recv_comm,old_nitem+count,
                                  &ext_recv_comm->nitem) ) {
          if ( !pmmgWarn ) {
            pmmgWarn = 1;
            fprintf(stderr,"  ## Error: %s: unable to resize at least one external"
                    " communicator\n",__func__);
          }
          ier0 = 0;
          ier  = MG_MIN(ier,ier0);
        }

        count = 0;
        if ( ier0 ) {
          for ( k=0; k<ext_send_comm->nitem; ++k ) {
            idx = ext_send_comm->int_comm_index[k];
            if ( PMMG_UNSET == recv_ext_idx[k] ) {
              /* Add the face to the external communicator myrank-recv */
              ext_recv_comm->int_comm_index[old_nitem++] = idx;
            }
            else {
              /* Remove the face from the external communicator myrank-sndr */
              ext_send_comm->int_comm_index[count++] = idx;
            }
          }
        }
        if ( !PMMG_resize_extComm(parmesh,ext_send_comm,count,
                                  &ext_send_comm->nitem) ) {
          ier0 = 0;
          ier  = MG_MIN(ier,ier0);
        }
      }
      PMMG_DEL_MEM(parmesh,recv_ext_idx,int,"recv_ext_idx");
    }
  }


  /** Step 3: Wait for the end of the MPI communications */
  if ( myrank == sndr) {
    MPI_CHECK( MPI_Waitall(nprocs,trequest,MPI_STATUSES_IGNORE), return 0 );
    PMMG_DEL_MEM ( parmesh, trequest,MPI_Request,"request_tab" );

    MPI_CHECK( MPI_Wait(&irequest,&status), return 0 );
    MPI_CHECK( MPI_Wait(&drequest,&status), return 0 );

    /* Free the memory */
    PMMG_DEL_MEM ( parmesh,grps2send,char,"grps2send" );
  }
  else if ( myrank == recv ) {
    MPI_CHECK( MPI_Wait(&irequest,&status), return 0 );
  }

  /** Step 4: pack the groups and face communicators */
  err = PMMG_pack_grps( parmesh,&parmesh->listgrp);
  ier = MG_MIN ( ier, err );
  if ( (!parmesh->ngrp) && parmesh->ddebug && (!called_from_distrib_mesh) ) {
    fprintf(stderr,"  ## Warning: %s: rank %d: processor without any groups.\n",
            __func__,myrank);
  }

  if ( ier ) {
    /* Pack the face communicators */
    err = PMMG_pack_faceCommunicators(parmesh);
    ier = MG_MIN ( ier, err );
  }

  PMMG_DEL_MEM ( parmesh,intcomm_flag,int,"intcomm_flag" );
  PMMG_DEL_MEM ( parmesh,recv_ext_idx,int,"recv_ext_idx" );

  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    if ( ext_face_comm->itosend ) {
      PMMG_DEL_MEM ( parmesh,ext_face_comm->itosend,int,"itosend" );
    }
  }

  return ier;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param part pointer toward the metis array containing the partitions.
 * \param called_from_distrib_mesh 1 if function is called from the
 * distributedmesh one. In this case we don't want to build the interaction map
 * (no interactions).
 *
 * \return 0 if fail but we can try to save a mesh, -1 if we fail and are unable
 * to save the mesh, 1 if we success
 *
 * Send the suitable groups to other procs and recieve their groups.
 * Deallocate the \a part array.
 *
 */
int PMMG_transfer_all_grps(PMMG_pParMesh parmesh,idx_t *part,int called_from_distrib_mesh) {
  PMMG_pGrp      grp;
  MPI_Comm       comm;
  int            myrank,nprocs;
  int            ninteractions,*interaction_map,*interactions;
  int            *send_grps,*recv_grps,*nfaces2send,*nfaces2recv;
  int            *next_comm2send,*ext_comms_next_idx,*nitems2send;
  int            *extComm_next_idx,*items_next_idx,*recv_array;
  int            *extComm_grpFaces2extComm,*extComm_grpFaces2face2int;
  int            max_ngrp;
  int            ier,k,i,j,err;

  myrank    = parmesh->myrank;
  nprocs    = parmesh->nprocs;
  comm      = parmesh->comm;
  max_ngrp  = 0;
  ier = 1;

  send_grps                 = NULL;
  recv_grps                 = NULL;
  nfaces2send               = NULL;
  nfaces2recv               = NULL;
  next_comm2send            = NULL;
  ext_comms_next_idx        = NULL;
  nitems2send               = NULL;
  items_next_idx            = NULL;
  extComm_next_idx          = NULL;
  extComm_grpFaces2extComm  = NULL;
  extComm_grpFaces2face2int = NULL;
  recv_array                = NULL;
  interaction_map           = NULL;
  interactions              = NULL;

  /** Step 1: Merge all the groups that must be sent to a given proc into 1
   * group */
  ier = PMMG_merge_grps2send(parmesh,&part);

  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp = &parmesh->listgrp[k];
    err = PMMG_fitMeshSize( parmesh,grp );

    ier = MG_MIN( err, ier );
  }

  /* Store the group destination in the group flag */
  for ( k=0; k<parmesh->ngrp; ++k ) {
    parmesh->listgrp[k].flag = part[k];
  }
  PMMG_DEL_MEM(parmesh,part,idx_t,"parmetis partition");

  MPI_Allreduce( MPI_IN_PLACE, &ier, 1, MPI_INT, MPI_MIN, comm);

  if ( ier < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to compute the new group"
            " partition.\n",__func__);
    goto end;
  }
  else if ( !ier ) {
    fprintf(stderr,"\n  ## Warning: %s: unable to compute the new group"
            " partition. Try to send it nevertheless.\n",__func__);
  }

  /** Step 2: Nodal communicators deletion (unallocation to remove, at least for
   * the external communicators) */
  PMMG_node_comm_free(parmesh);

  /** Step 3: Compute the map of interactions, if fail, consider that each proc
   * interacts with all the other procs (worst case)
   */
  if ( (!called_from_distrib_mesh) && parmesh->nprocs > 2 ) {
    ninteractions = PMMG_interactionMap(parmesh,&interactions,&interaction_map);
  }
  else {
    ninteractions = -1;
  }

  if ( ninteractions <= 0 ) {
    PMMG_DEL_MEM ( parmesh, interactions, int, "interactions");
    if ( !interaction_map ) {
      PMMG_CALLOC ( parmesh,interaction_map,nprocs,int,"interaction_map" ,ier=0 );
    }
    for ( k=0; k<nprocs; ++k ) {
      interaction_map[k] = 1;
    }
  }

  /** Step 4: proc k send its data (group and/or communicators), proc j receive
   * data */
  ier = 1;

  if ( interactions ) {
    for ( k=0; k<ninteractions; ++k ) {
      i = interactions[2*k];
      j = interactions[2*k+1];
      if ( i==j ) {
        continue;
      }
      err =  PMMG_transfer_grps_fromItoJ(parmesh,i,j,interaction_map,called_from_distrib_mesh);
      ier = MG_MIN ( ier,err );
    }
  }
  else {
    for ( k=0; k<nprocs; ++k ) {
      for ( j=0; j<nprocs; ++j ) {
        if ( j==k ) {
          continue;
        }
        err =  PMMG_transfer_grps_fromItoJ(parmesh,k,j,interaction_map,called_from_distrib_mesh);
        ier = MG_MIN ( ier,err );
      }
    }
  }
  MPI_Allreduce( MPI_IN_PLACE, &ier, 1, MPI_INT, MPI_MIN, comm);

  if ( ier <= 0 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to transfer the groups toward"
            " their new processor.\n",__func__);
    ier = -1;
    goto end;
  }

  /** Step 6: Node communicators reconstruction from the face ones */
  if ( !PMMG_build_nodeCommFromFaces(parmesh,parmesh->comm) ) {
    fprintf(stderr,"\n  ## Unable to build the new node communicators from"
            " the face ones.\n");
    ier = -1;
    goto end;
  }

#ifndef NDEBUG
  for ( k=0; k<parmesh->ngrp; ++k ) {
    for ( j=0; j<parmesh->listgrp[k].nitem_int_face_comm; ++j )
      assert( parmesh->listgrp[k].face2int_face_comm_index2[j]>=0 &&
              "face2int_face_comm not filled.");
  }
#endif

  /** Success */
  if ( !err )  ier = 1;

end:
  if ( interaction_map )
    PMMG_DEL_MEM(parmesh,interaction_map,int,"interaction_map");
  if ( interactions )
    PMMG_DEL_MEM(parmesh,interactions,int,"interactions");

  if ( send_grps )
    PMMG_DEL_MEM(parmesh,send_grps,int,"send_grps");
  if ( recv_grps )
    PMMG_DEL_MEM(parmesh,recv_grps,int,"recv_grps");
  if ( parmesh->int_face_comm->intvalues ) {
    PMMG_DEL_MEM(parmesh,parmesh->int_face_comm->intvalues,
                 int,"intvalues");
  }
  if ( recv_array )
    PMMG_DEL_MEM(parmesh,recv_array,int,"recv_array");

  /* Arrays allocated in fill_extFaceComm2send */
  if ( nfaces2send )
    PMMG_DEL_MEM(parmesh,nfaces2send,int,"nfaces2send");
  if ( nfaces2recv )
    PMMG_DEL_MEM(parmesh,nfaces2recv,int,"nfaces2recv");

  if ( extComm_next_idx )
    PMMG_DEL_MEM(parmesh,extComm_next_idx,int,
                 "positions of the external comms in pos_grp_faces* arrays");
  if ( extComm_grpFaces2extComm )
    PMMG_DEL_MEM(parmesh,extComm_grpFaces2extComm,
                 int,"extComm_grpFaces2extComm");
  if ( extComm_grpFaces2face2int )
    PMMG_DEL_MEM(parmesh,extComm_grpFaces2face2int,
                 int,"extComm_grpFaces2face2int");

  if ( next_comm2send )
    PMMG_DEL_MEM(parmesh,next_comm2send,int,"next_comm2send");
  if ( nitems2send )
    PMMG_DEL_MEM(parmesh,nitems2send,int,"nitems2send");
  if ( ext_comms_next_idx )
    PMMG_DEL_MEM(parmesh,ext_comms_next_idx,int,"ext_comms_next_idx");
  if ( items_next_idx )
    PMMG_DEL_MEM(parmesh,items_next_idx,int,"items_next_idx");

  return ier;
}


/**
 * \param parmesh pointer toward the mesh structure.
 * \param repartitionning_mode strategy to use for repartitionning (we want to balance
 * the graph if the function is called at preprocessing stage due to inputs from
 * a different number of parititions than the number used during the run (before
 * libparmmg1 call) or at the end of the iterations to return a balanced mesh
 * but we want to perform interfaces migration during internal iterations )
 *
 * \return -1 if we fail and can not save the meshes, 0 if we fail but can save
 * the meshes, 1 otherwise
 *
 * Compute the group distribution that load balance the groups over the
 * processors and send and recieve the groups from the other processors.
 *
 */
int PMMG_distribute_grps( PMMG_pParMesh parmesh, int repartitioning_mode ) {
  idx_t *part;
  int   ngrp,ier;

  MPI_Allreduce( &parmesh->ngrp, &ngrp, 1, MPI_INT, MPI_MIN, parmesh->comm);

  /** Get the new partition of groups (1 group = 1 metis node) */
  part = NULL;
  PMMG_CALLOC(parmesh,part,parmesh->ngrp,idx_t,"allocate parmetis buffer",
              return 0);

  if( repartitioning_mode == PMMG_REDISTRIBUTION_ifc_displacement ) {

    ier = PMMG_part_getProcs( parmesh, part );

  } else {

    switch ( parmesh->info.loadbalancing_mode ) {

#ifdef USE_PARMETIS
    case PMMG_LOADBALANCING_parmetis:
      ier = PMMG_part_parmeshGrps2parmetis(parmesh,part,parmesh->nprocs);
      break;
#endif

    case PMMG_LOADBALANCING_metis:
    default:

      ier = PMMG_part_parmeshGrps2metis(parmesh,part,parmesh->nprocs);
      break;
    }
  }

  if ( !ier )
    fprintf(stderr,"\n  ## Unable to compute the new group partition.\n");

  MPI_Allreduce( MPI_IN_PLACE, &ier, 1, MPI_INT, MPI_MIN, parmesh->comm);

  if ( !ier ) {
    PMMG_DEL_MEM(parmesh,part,idx_t,"deallocate parmetis partition");
    return ier;
  }

  /** Send the suitable groups to other procs and recieve their groups */
  ier = PMMG_transfer_all_grps(parmesh,part,0);
  if ( ier <= 0 )
    fprintf(stderr,"\n  ## Unable to communicate groups through processors.\n");

  /** Check grps contiguity */
  ier = PMMG_checkAndReset_grps_contiguity( parmesh );

  return ier;
}

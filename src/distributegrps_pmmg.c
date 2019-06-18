/**
 * \file distributegrps_pmmg.c
 * \brief Group distribution on the processors
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \version 5
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
  out.sol  = group->sol;
  out.met  = group->met;
  out.disp = group->disp;

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

  group->mesh = NULL;
  group->sol  = NULL;
  group->met  = NULL;
  group->disp = NULL;

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
  size_t        memAv,oldMemMax;
  int           nprocs,ngrp,k,j,ier;

  nprocs = parmesh->nprocs;
  ngrp   = parmesh->ngrp;

  /** Step 0: Set mesh sizes to their minimal size and count the available
   * memory */
  parmesh->memMax = parmesh->memCur;
  assert ( parmesh->memGloMax > parmesh->memMax );
  memAv = parmesh->memGloMax - parmesh->memMax;
  for ( k=0; k<ngrp; ++k ) {
    parmesh->listgrp[k].mesh->memMax = parmesh->listgrp[k].mesh->memCur;
    assert ( memAv > parmesh->listgrp[k].mesh->memMax );
    memAv -= parmesh->listgrp[k].mesh->memMax;
  }

  /** Step 1: New groups allocation and initialization: move the groups to have
   * a group that will be send to proc k stored in grps[k]. Free the adja
   * array. */

  /* If this step is too expensive in memory, we can count the number of procs
     with which we will communicate and fill directly the pack array. */
  oldMemMax        = parmesh->memMax;
  parmesh->memMax += memAv;
  PMMG_CALLOC( parmesh,grps,nprocs,PMMG_Grp,"Groups to send",return 0 );

  j = 0;
  for ( k=0; k<ngrp; ++k ) {
    /* Free the adja array */
    meshI = parmesh->listgrp[k].mesh;
    if ( meshI->adja )
      PMMG_DEL_MEM(meshI, meshI->adja, int, "adjacency table" );

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

    if ( !listgrp[k].mesh ) continue;

    grpI = &grps[(*part)[k]];
    grpJ = &listgrp[k];

    /* Free the adja arrays */
    meshI = grpI->mesh;
    meshJ = grpJ->mesh;
    if ( meshI->adja )
      PMMG_DEL_MEM(meshI, meshI->adja,int, "adjacency table" );
    if ( meshJ->adja )
      PMMG_DEL_MEM(meshJ, meshJ->adja,int, "adjacency table" );

    /* Update the available memory and give it to the mesh */
    assert ( meshI->memMax >= meshI->memCur );
    assert ( meshJ->memMax >= meshJ->memCur );
    memAv += ( meshI->memMax-meshI->memCur );
    memAv += ( meshJ->memMax-meshJ->memCur );
    meshI->memMax = meshI->memCur;
    meshJ->memMax = meshJ->memCur;

    PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,meshI,memAv,oldMemMax);

    if ( !PMMG_merge_grpJinI(parmesh,grpI,grpJ) ) goto low_fail;

    /* Update the communicators: WARNING TO IMPROVE: INEFFICIENT */
    PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,meshI,memAv,oldMemMax);

    if ( !PMMG_mergeGrpJinI_communicators(parmesh,grpI,grpJ,grps,k) ) goto low_fail;

    /* Delete the useless group to gain memory space */
    memAv += meshJ->memCur;
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
  PMMG_DEL_MEM(parmesh,parmesh->listgrp,PMMG_Grp,"listgrp");
  parmesh->listgrp = grps;
  parmesh->ngrp    = ngrp;

  /* Pack the communicators */
  if ( !PMMG_pack_nodeCommunicators(parmesh) ) ier = -1;

  if ( !PMMG_pack_faceCommunicators(parmesh) ) ier = -1;

  assert ( PMMG_check_extFaceComm(parmesh) );
  assert ( PMMG_check_extNodeComm(parmesh) );

  /* Update tag on points, tetra */
  if ( !PMMG_updateTag(parmesh) ) return -1;

  if ( !PMMG_parmesh_updateMemMax(parmesh, 5, 1) ) {
    fprintf(stderr,"\n  ## Error: %s: Unable to update the memory repartition"
            " between meshes and communicators.\n",__func__);
    return -1;
  }

  return ier;
}

/**
 * \param grp pointer toward a PMMG_Grp structure.
 * \param buffer pointer toward the buffer in which we pack the group
 * \param idx index from which we start to write in \a buffer
 *
 * \warning the mesh prisms are not treated.
 *
 * \warning IMPROVEMENT: do not use an mpi type, just pack the
 * point+xpoint+tetra+xtetra integer and char fields into one unique array and
 * the point+xpoint double fields into another one.
 *
 * Pack a group into a double buffer to allow mpi communication.
 *
 */
static
int PMMG_mpisizeof_grp ( PMMG_pGrp grp ) {
  const MMG5_pMesh mesh = grp->mesh;
  const MMG5_pSol  met  = grp->met;

  int k,i,idx;

  idx = 0;
  if ( !grp->mesh ) {
    /* unused group */
    return sizeof(int);
  }

  /* used group */
  idx += sizeof(int);

  /** Mesh size */
  idx += sizeof(int); // mesh->np
  idx += sizeof(int); // mesh->xp
  idx += sizeof(int); // mesh->ne
  idx += sizeof(int); // mesh->xt

  /** Mesh names */
  idx += sizeof(int); // meshin
  idx += sizeof(int); // metin
  idx += sizeof(int); // meshout
  idx += sizeof(int); // metout

  idx += (strlen(mesh->namein) + 1) * sizeof(char);
  idx += (strlen(met->namein) + 1) * sizeof(char);
  idx += (strlen(mesh->nameout) + 1) * sizeof(char);
  idx += (strlen(met->nameout) + 1) * sizeof(char);

  /** Mesh infos: warning, some "useless" info are not sended */
  idx += sizeof(double); // mesh->info.dhd
  idx += sizeof(double); // mesh->info.hmin
  idx += sizeof(double); // mesh->info.hmax
  idx += sizeof(double); // mesh->info.hsiz
  idx += sizeof(double); // mesh->info.hgrad
  //idx += sizeof(double); // mesh->info.hgradreq
  idx += sizeof(double); // mesh->info.hausd

  idx += sizeof(double); // mesh->info.delta
  idx += sizeof(double); // mesh->info.min[0]
  idx += sizeof(double); // mesh->info.min[1]
  idx += sizeof(double); // mesh->info.min[2]

  idx += sizeof(double); // mesh->info.ls

  idx += sizeof(int); // npar
  idx += sizeof(int); // openbdy
  idx += sizeof(int); // renum
  idx += sizeof(int); // PROctree
  idx += sizeof(int); // nmat

  idx += sizeof(char); // nreg
  idx += sizeof(char); // imprim
  idx += sizeof(char); // ddebug
  idx += sizeof(char); // iso
  idx += sizeof(char); // lag
  idx += sizeof(char); // parTyp
  idx += sizeof(char); // optim
  idx += sizeof(char); // optimLES
  idx += sizeof(char); // noinsert
  idx += sizeof(char); // noswap
  idx += sizeof(char); // nomove
  idx += sizeof(char); // nosurf
  idx += sizeof(char); // inputMet

  /* affectation of old refs in ls-mode */
  if ( mesh->info.nmat ) {
    assert( mesh->info.mat );
    for ( k=0; k<mesh->info.nmat; ++k ) {
      idx += sizeof(char); // mat->dospl
      idx += sizeof(int); //  mat->ref
      idx += sizeof(int); //  mat->rin
      idx += sizeof(int); //  mat->rex
    }
  }

  /* local parameters */
  if ( mesh->info.npar ) {
    assert( mesh->info.par );
    for ( k=0; k<mesh->info.npar; ++k ) {
      idx += sizeof(double); // par->hmin
      idx += sizeof(double); // par->hmax
      idx += sizeof(double); // par->hausd
      idx += sizeof(int); //  par->ref
      idx += sizeof(char); // par->elt
    }
  }

  /** Metric info and sizes */
  idx += sizeof(int); // (met->m ? 1 : 0 );
  idx += sizeof(int); // met->size;
  idx += sizeof(int); // met->type;

  /** Communicator sizes */
  idx += sizeof(int); // grp->nitem_int_node_comm;
  idx += sizeof(int); // grp->nitem_int_face_comm;

  /** Pack mesh points */
  for ( k=1; k<=mesh->np; ++k ) {
    /* Coordinates */
    idx += sizeof(double); // mesh->point[k].c[0];
    idx += sizeof(double); // mesh->point[k].c[1];
    idx += sizeof(double); // mesh->point[k].c[2];
    /* Tangent */
    idx += sizeof(double); // mesh->point[k].n[0];
    idx += sizeof(double); // mesh->point[k].n[1];
    idx += sizeof(double); // mesh->point[k].n[2];
    /* Pointer toward the boundary entity */
    idx += sizeof(int); // mesh->point[k].xp;
    /* Ref */
    idx += sizeof(int); // mesh->point[k].ref;
    /* Tag */
    idx += sizeof(int16_t); // mesh->point[k].tag;
  }

  /** Pack mesh boundary points */
  for ( k=1; k<=mesh->xp; ++k ) {
    /* First normal */
    idx += sizeof(double); // mesh->xpoint[k].n1[0];
    idx += sizeof(double); // mesh->xpoint[k].n1[1];
    idx += sizeof(double); // mesh->xpoint[k].n1[2];
    /* Second normal */
    idx += sizeof(double); // mesh->xpoint[k].n2[0];
    idx += sizeof(double); // mesh->xpoint[k].n2[1];
    idx += sizeof(double); // mesh->xpoint[k].n2[2];
  }

  /** Pack mesh elements */
  for ( k=1; k<=mesh->ne; ++k ) {
    /* Tetra vertices */
    idx += sizeof(int); // mesh->tetra[k].v[0];
    idx += sizeof(int); // mesh->tetra[k].v[1];
    idx += sizeof(int); // mesh->tetra[k].v[2];
    idx += sizeof(int); // mesh->tetra[k].v[3];
    /* Pointer toward the boundary entity */
    idx += sizeof(int); // mesh->tetra[k].xt;
    /* Ref */
    idx += sizeof(int); // mesh->tetra[k].ref;
    /* Mark */
    idx += sizeof(int); // mesh->point[k].mark;
    /* Tag */
    idx += sizeof(int16_t); // mesh->tetra[k].tag;
    /* Quality */
    idx += sizeof(double); // mesh->tetra[k].qual;
  }

  /** Pack mesh boundary tetra */
  for ( k=1; k<=mesh->xt; ++k ) {
    /* Faces references  */
    idx += sizeof(int); // mesh->xtetra[k].ref[0];
    idx += sizeof(int); // mesh->xtetra[k].ref[1];
    idx += sizeof(int); // mesh->xtetra[k].ref[2];
    idx += sizeof(int); // mesh->xtetra[k].ref[3];
    /* Faces tags */
    idx += sizeof(int16_t); // mesh->xtetra[k].ftag[0];
    idx += sizeof(int16_t); // mesh->xtetra[k].ftag[1];
    idx += sizeof(int16_t); // mesh->xtetra[k].ftag[2];
    idx += sizeof(int16_t); // mesh->xtetra[k].ftag[3];
    /* Edges references  */
    idx += sizeof(int); // mesh->xtetra[k].edg[0];
    idx += sizeof(int); // mesh->xtetra[k].edg[1];
    idx += sizeof(int); // mesh->xtetra[k].edg[2];
    idx += sizeof(int); // mesh->xtetra[k].edg[3];
    idx += sizeof(int); // mesh->xtetra[k].edg[4];
    idx += sizeof(int); // mesh->xtetra[k].edg[5];
    /* Edges tags */
    idx += sizeof(int16_t); // mesh->xtetra[k].tag[0];
    idx += sizeof(int16_t); // mesh->xtetra[k].tag[1];
    idx += sizeof(int16_t); // mesh->xtetra[k].tag[2];
    idx += sizeof(int16_t); // mesh->xtetra[k].tag[3];
    idx += sizeof(int16_t); // mesh->xtetra[k].tag[4];
    idx += sizeof(int16_t); // mesh->xtetra[k].tag[5];
  }

  /** Pack metric */
  if ( met->m ) {
    for ( k=1; k<=met->np; ++k ) {
      for ( i=0; i<met->size; ++i ) {
        idx += sizeof(double); // met->m[met->size*k + i];
      }
    }
  }

  /** Pack communicators */
  /* Node communicator */
  for ( k=0; k<grp->nitem_int_node_comm; ++k ) {
    idx += sizeof(int); // grp->node2int_node_comm_index1[k];
  }
  for ( k=0; k<grp->nitem_int_node_comm; ++k ) {
    idx += sizeof(int); // grp->node2int_node_comm_index2[k];
  }
  /* Face communicator */
  for ( k=0; k<grp->nitem_int_face_comm; ++k ) {
    idx += sizeof(int); // grp->face2int_face_comm_index1[k];
  }
  for ( k=0; k<grp->nitem_int_face_comm; ++k ) {
    idx += sizeof(int); // grp->face2int_face_comm_index2[k];
  }

  return idx;
}

/**
 * \param grp pointer toward a PMMG_Grp structure.
 * \param buffer pointer toward the buffer in which we pack the group
 *
 * \return 1 if success, 0 if fail
 *
 * \warning the mesh prisms are not treated.
 *
 * \remark To store the correct value of a variable of type "x" into buffer, we
 * need to cast the void pointer toward the buffer into a pointer toward a
 * buffer of type "x". Then we can store the variable value by dereferencing the
 * adress of the buffer.
 *
 * Pack a group into a buffer (to allow mpi communication) and shift the buffer
 * pointer at the end of the written area.
 *
 */
static
int PMMG_mpipack_grp ( PMMG_pGrp grp,char **buffer ) {
  const MMG5_pMesh mesh = grp->mesh;
  const MMG5_pSol  met  = grp->met;

  int   k,i,ier,meshin_s,meshout_s,metin_s,metout_s;
  char  *tmp;

  ier = 1;

  tmp = *buffer;

  if ( !grp->mesh ) {
    /* unused group */
    *( (int *) tmp ) = 0; tmp += sizeof(int);

    return ier;
  }

  /* used group */
  *( (int *) tmp) = 1; tmp += sizeof(int);

  /** Mesh size */
  *( (int *) tmp) = mesh->np; tmp += sizeof(int);
  *( (int *) tmp) = mesh->xp; tmp += sizeof(int);
  *( (int *) tmp) = mesh->ne; tmp += sizeof(int);
  *( (int *) tmp) = mesh->xt; tmp += sizeof(int);

 /** Mesh names */
  meshin_s  = (strlen(mesh->namein) + 1);
  metin_s   = (strlen(met->namein) + 1);
  meshout_s = (strlen(mesh->nameout) + 1);
  metout_s  = (strlen(met->nameout) + 1);

  if ( meshin_s > 255 || metin_s > 255 || meshout_s > 255 || metout_s > 255) {
    printf("  ## Error: input filenames too long.");
    return 0;
  }

  *( (int *) tmp) = meshin_s; tmp += sizeof(int);
  *( (int *) tmp) = metin_s; tmp += sizeof(int);
  *( (int *) tmp) = meshout_s; tmp += sizeof(int);
  *( (int *) tmp) = metout_s; tmp += sizeof(int);

  strncpy( ( (char *) tmp), mesh->namein,meshin_s);
  tmp += meshin_s * sizeof(char);
  strncpy( ( (char *) tmp), met->namein,metin_s);
  tmp += metin_s * sizeof(char);
  strncpy( ( (char *) tmp), mesh->nameout,meshout_s);
  tmp += meshout_s * sizeof(char);
  strncpy( ( (char *) tmp), met->nameout, metout_s);
  tmp += metout_s * sizeof(char);

  /** Mesh infos */
  *( (double *) tmp) = mesh->info.dhd;      tmp += sizeof(double);
  *( (double *) tmp) = mesh->info.hmin;     tmp += sizeof(double);
  *( (double *) tmp) = mesh->info.hmax;     tmp += sizeof(double);
  *( (double *) tmp) = mesh->info.hsiz;     tmp += sizeof(double);
  *( (double *) tmp) = mesh->info.hgrad;    tmp += sizeof(double);
  //*( (double *) tmp) = mesh->info.hgradreq; tmp += sizeof(double);
  *( (double *) tmp) = mesh->info.hausd;    tmp += sizeof(double);

  *( (double *) tmp) = mesh->info.delta;    tmp += sizeof(double);
  *( (double *) tmp) = mesh->info.min[0];   tmp += sizeof(double);
  *( (double *) tmp) = mesh->info.min[1];   tmp += sizeof(double);
  *( (double *) tmp) = mesh->info.min[2];   tmp += sizeof(double);

  *( (double *) tmp) = mesh->info.ls;       tmp += sizeof(double);

  *( (int *) tmp) = mesh->info.npar;      tmp += sizeof(int);
  *( (int *) tmp) = mesh->info.opnbdy;    tmp += sizeof(int);
  *( (int *) tmp) = mesh->info.renum;     tmp += sizeof(int);
  *( (int *) tmp) = mesh->info.PROctree;  tmp += sizeof(int);
  *( (int *) tmp) = mesh->info.nmat;      tmp += sizeof(int);

  *( (char *) tmp) = mesh->info.nreg;     tmp += sizeof(char);
  *( (char *) tmp) = mesh->info.imprim;   tmp += sizeof(char);
  *( (char *) tmp) = mesh->info.ddebug;   tmp += sizeof(char);
  *( (char *) tmp) = mesh->info.iso;      tmp += sizeof(char);
  *( (char *) tmp) = mesh->info.lag;      tmp += sizeof(char);
  *( (char *) tmp) = mesh->info.parTyp;   tmp += sizeof(char);
  *( (char *) tmp) = mesh->info.optim;    tmp += sizeof(char);
  *( (char *) tmp) = mesh->info.optimLES; tmp += sizeof(char);
  *( (char *) tmp) = mesh->info.noinsert; tmp += sizeof(char);
  *( (char *) tmp) = mesh->info.noswap;   tmp += sizeof(char);
  *( (char *) tmp) = mesh->info.nomove;   tmp += sizeof(char);
  *( (char *) tmp) = mesh->info.nosurf;   tmp += sizeof(char);
  *( (char *) tmp) = mesh->info.inputMet; tmp += sizeof(char);

  /* affectation of old refs in ls-mode */
  if ( mesh->info.nmat ) {
    assert( mesh->info.mat );
    for ( k=0; k<mesh->info.nmat; ++k ) {
      *( (char *) tmp) = mesh->info.mat[k].dospl; tmp += sizeof(char);
      *( (int *) tmp)  = mesh->info.mat[k].ref; tmp += sizeof(int);
      *( (int *) tmp)  = mesh->info.mat[k].rin; tmp += sizeof(int);
      *( (int *) tmp)  = mesh->info.mat[k].rex; tmp += sizeof(int);
    }
  }

  /* local parameters */
  if ( mesh->info.npar ) {
    assert( mesh->info.par );
    for ( k=0; k<mesh->info.npar; ++k ) {
      *( (double *) tmp) =  mesh->info.par[k].hmin; tmp += sizeof(double);
      *( (double *) tmp) =  mesh->info.par[k].hmax; tmp += sizeof(double);
      *( (double *) tmp) =  mesh->info.par[k].hausd; tmp += sizeof(double);
      *( (int *)    tmp) =  mesh->info.par[k].ref; tmp += sizeof(int);
      *( (char *)   tmp) =  mesh->info.par[k].elt; tmp += sizeof(char);
    }
  }

  /** Metric info and sizes */
  assert ( mesh->npmax == met->npmax );
  assert ( mesh->np    == met->np    );

  *( (int *) tmp) = (met->m ? 1 : 0 );  tmp += sizeof(int);
  *( (int *) tmp) = met->size;          tmp += sizeof(int);
  *( (int *) tmp) = met->type;          tmp += sizeof(int);

  /** Communicator sizes */
  *( (int *) tmp) = grp->nitem_int_node_comm; tmp += sizeof(int);
  *( (int *) tmp) = grp->nitem_int_face_comm; tmp += sizeof(int);

  /** Pack mesh points */
  for ( k=1; k<=mesh->np; ++k ) {
    /* Coordinates */
    *( (double *) tmp) = mesh->point[k].c[0]; tmp += sizeof(double);
    *( (double *) tmp) = mesh->point[k].c[1]; tmp += sizeof(double);
    *( (double *) tmp) = mesh->point[k].c[2]; tmp += sizeof(double);
    /* Tangent */
    *( (double *) tmp) = mesh->point[k].n[0]; tmp += sizeof(double);
    *( (double *) tmp) = mesh->point[k].n[1]; tmp += sizeof(double);
    *( (double *) tmp) = mesh->point[k].n[2]; tmp += sizeof(double);
    /* Pointer toward the boundary entity */
    *( (int *) tmp) = mesh->point[k].xp; tmp += sizeof(int);
    /* Ref */
    *( (int *) tmp) = mesh->point[k].ref; tmp += sizeof(int);
    /* Tag */
    *( (int16_t *) tmp) = mesh->point[k].tag; tmp += sizeof(int16_t);
  }

  /** Pack mesh boundary points */
  for ( k=1; k<=mesh->xp; ++k ) {
    /* First normal */
    *( (double *) tmp) = mesh->xpoint[k].n1[0]; tmp += sizeof(double);
    *( (double *) tmp) = mesh->xpoint[k].n1[1]; tmp += sizeof(double);
    *( (double *) tmp) = mesh->xpoint[k].n1[2]; tmp += sizeof(double);
    /* Second normal */
    *( (double *) tmp) = mesh->xpoint[k].n2[0]; tmp += sizeof(double);
    *( (double *) tmp) = mesh->xpoint[k].n2[1]; tmp += sizeof(double);
    *( (double *) tmp) = mesh->xpoint[k].n2[2]; tmp += sizeof(double);
  }

  /** Pack mesh elements */
  for ( k=1; k<=mesh->ne; ++k ) {
    /* Tetra vertices */
    *( (int *) tmp) = mesh->tetra[k].v[0]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->tetra[k].v[1]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->tetra[k].v[2]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->tetra[k].v[3]; tmp += sizeof(int);
    /* Pointer toward the boundary entity */
    *( (int *) tmp) = mesh->tetra[k].xt; tmp += sizeof(int);
    /* Ref */
    *( (int *) tmp) = mesh->tetra[k].ref; tmp += sizeof(int);
    /* Mark */
    *( (int *) tmp) = mesh->tetra[k].mark; tmp += sizeof(int);
    /* Tag */
    *( (int16_t *) tmp) = mesh->tetra[k].tag; tmp += sizeof(int16_t);
    /* Quality */
    *( (double *) tmp) = mesh->tetra[k].qual; tmp += sizeof(double);
  }

  /** Pack mesh boundary tetra */
  for ( k=1; k<=mesh->xt; ++k ) {
    /* Faces references  */
    *( (int *) tmp) = mesh->xtetra[k].ref[0]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->xtetra[k].ref[1]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->xtetra[k].ref[2]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->xtetra[k].ref[3]; tmp += sizeof(int);
    /* Faces tags */
    *( (int16_t *) tmp) = mesh->xtetra[k].ftag[0]; tmp += sizeof(int16_t);
    *( (int16_t *) tmp) = mesh->xtetra[k].ftag[1]; tmp += sizeof(int16_t);
    *( (int16_t *) tmp) = mesh->xtetra[k].ftag[2]; tmp += sizeof(int16_t);
    *( (int16_t *) tmp) = mesh->xtetra[k].ftag[3]; tmp += sizeof(int16_t);
    /* Edges references  */
    *( (int *) tmp) = mesh->xtetra[k].edg[0]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->xtetra[k].edg[1]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->xtetra[k].edg[2]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->xtetra[k].edg[3]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->xtetra[k].edg[4]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->xtetra[k].edg[5]; tmp += sizeof(int);
    /* Edges tags */
    *( (int16_t *) tmp) = mesh->xtetra[k].tag[0]; tmp += sizeof(int16_t);
    *( (int16_t *) tmp) = mesh->xtetra[k].tag[1]; tmp += sizeof(int16_t);
    *( (int16_t *) tmp) = mesh->xtetra[k].tag[2]; tmp += sizeof(int16_t);
    *( (int16_t *) tmp) = mesh->xtetra[k].tag[3]; tmp += sizeof(int16_t);
    *( (int16_t *) tmp) = mesh->xtetra[k].tag[4]; tmp += sizeof(int16_t);
    *( (int16_t *) tmp) = mesh->xtetra[k].tag[5]; tmp += sizeof(int16_t);
  }

  /** Pack metric */
  if ( met->m ) {
    for ( k=1; k<=met->np; ++k ) {
      for ( i=0; i<met->size; ++i ) {
        *( (double *) tmp) = met->m[met->size*k + i]; tmp += sizeof(double);
      }
    }
  }

  /** Pack communicators */
  /* Node communicator */
  for ( k=0; k<grp->nitem_int_node_comm; ++k ) {
    *( (int *) tmp) = grp->node2int_node_comm_index1[k]; tmp += sizeof(int);
  }
  for ( k=0; k<grp->nitem_int_node_comm; ++k ) {
    *( (int *) tmp) = grp->node2int_node_comm_index2[k]; tmp += sizeof(int);
  }
  /* Face communicator */
  for ( k=0; k<grp->nitem_int_face_comm; ++k ) {
    *( (int *) tmp) = grp->face2int_face_comm_index1[k]; tmp += sizeof(int);
  }
  for ( k=0; k<grp->nitem_int_face_comm; ++k ) {
    *( (int *) tmp) = grp->face2int_face_comm_index2[k]; tmp += sizeof(int);
  }

  *buffer = tmp;

  return ier;
}


/**
 * \param parmesh pointer toward a parmesh structure.
 * \param grp pointer toward a PMMG_Grp structure.
 * \param buffer pointer toward the buffer in which we unpack the group
 * \param memAv pointer toward the available memory whose value is updated.
 *
 * \return 0 if fail, 1 otherwise
 *
 * \warning the mesh prisms are not treated.
 *
 * Upack a group from a double buffer and shift the pointer toward the buffer to
 * point to the next group stored in the buffer.
 *
 * \remark  To get the correct value of a variable
 * of type "x", we need to cast the void pointer toward the buffer into a
 * pointer toward a buffer of type "x". Then we can get the variable value by
 * dereferencing the adress of the buffer.
 *
 */
static
int PMMG_mpiunpack_grp ( PMMG_pParMesh parmesh,PMMG_pGrp grp,char **buffer,
                         size_t *memAv) {
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  double     ddummy;
  int        k,i,ier,ier_grp,ier_mesh,ier_sol,ier_comm,np,xp,ne,xt;
  int        size,ismet,used,idummy,meshin_s,metin_s,meshout_s,metout_s;
  int16_t    i16dummy;
  char       cdummy,chaine[256];

  ier = 1;

  used = *( (int *) *buffer ); *buffer += sizeof(int);
  if ( !used ) {
    /* unused group */
    return ier;
  }

  /* used group */
  /** Mesh size */
  ier_grp = MMG3D_Init_mesh(MMG5_ARG_start,
                            MMG5_ARG_ppMesh,&(grp->mesh),
                            MMG5_ARG_ppMet ,&(grp->met),
                            MMG5_ARG_end);

  mesh = grp->mesh;
  met  = grp->met;

  /** Get the mesh maximal authorized memory */
  np = *( (int *) *buffer); *buffer += sizeof(int);
  xp = *( (int *) *buffer); *buffer += sizeof(int);
  ne = *( (int *) *buffer); *buffer += sizeof(int);
  xt = *( (int *) *buffer); *buffer += sizeof(int);

  if ( ier_grp ) {
    /* Give all the available memory to the mesh */
    mesh->memMax = *memAv;

    /** Set the mesh size */
    ier_mesh = PMMG_grpSplit_setMeshSize( mesh,np,ne,0,xp,xt );
  }
  else ier = ier_mesh = 0;

  if ( ier_mesh ) {

    /* File names */
    meshin_s = *( (int *) *buffer); *buffer += sizeof(int);
    metin_s = *( (int *) *buffer); *buffer += sizeof(int);
    meshout_s = *( (int *) *buffer); *buffer += sizeof(int);
    metout_s = *( (int *) *buffer); *buffer += sizeof(int);

    if ( meshin_s > 255 || metin_s > 255 || meshout_s > 255 || metout_s > 255) {
      printf("  ## Error: input filenames too long.");
      return 0;
    }
    else {
      strncpy(chaine,( (char *) *buffer),meshin_s); *buffer += meshin_s * sizeof(char);
      ier = MMG5_Set_inputMeshName( mesh, chaine );
      strncpy(chaine,( (char *) *buffer),metin_s); *buffer += metin_s * sizeof(char);
      if ( !MMG5_Set_inputSolName( mesh, met, chaine ) ) { ier = 0; }
      strncpy(chaine,( (char *) *buffer),meshout_s); *buffer += meshout_s * sizeof(char);
      if ( !MMG5_Set_outputMeshName( mesh, chaine ) ) { ier = 0; }
      strncpy(chaine,( (char *) *buffer),metout_s); *buffer += metout_s * sizeof(char);
      if ( !MMG5_Set_outputSolName( mesh,met, chaine ) ) { ier = 0; }
    }

    /** Mesh infos */
    mesh->info.dhd       = *( (double *) *buffer); *buffer += sizeof(double);
    mesh->info.hmin      = *( (double *) *buffer); *buffer += sizeof(double);
    mesh->info.hmax      = *( (double *) *buffer); *buffer += sizeof(double);
    mesh->info.hsiz      = *( (double *) *buffer); *buffer += sizeof(double);
    mesh->info.hgrad     = *( (double *) *buffer); *buffer += sizeof(double);
    //mesh->info.hgradreq  = *( (double *) *buffer); *buffer += sizeof(double);
    mesh->info.hausd     = *( (double *) *buffer); *buffer += sizeof(double);

    mesh->info.delta  = *( (double *) *buffer); *buffer += sizeof(double);
    mesh->info.min[0] = *( (double *) *buffer); *buffer += sizeof(double);
    mesh->info.min[1] = *( (double *) *buffer); *buffer += sizeof(double);
    mesh->info.min[2] = *( (double *) *buffer); *buffer += sizeof(double);

    mesh->info.ls     = *( (double *) *buffer); *buffer += sizeof(double);

    mesh->info.npar      = *( (int *) *buffer); *buffer += sizeof(int);
    mesh->info.opnbdy    = *( (int *) *buffer); *buffer += sizeof(int);
    mesh->info.renum     = *( (int *) *buffer); *buffer += sizeof(int);
    mesh->info.PROctree  = *( (int *) *buffer); *buffer += sizeof(int);
    mesh->info.nmat      = *( (int *) *buffer); *buffer += sizeof(int);

    mesh->info.nreg      = *( (char *) *buffer); *buffer += sizeof(char);
    mesh->info.imprim    = *( (char *) *buffer); *buffer += sizeof(char);
    mesh->info.ddebug    = *( (char *) *buffer); *buffer += sizeof(char);
    mesh->info.iso       = *( (char *) *buffer); *buffer += sizeof(char);
    mesh->info.lag       = *( (char *) *buffer); *buffer += sizeof(char);
    mesh->info.parTyp    = *( (char *) *buffer); *buffer += sizeof(char);
    mesh->info.optim     = *( (char *) *buffer); *buffer += sizeof(char);
    mesh->info.optimLES  = *( (char *) *buffer); *buffer += sizeof(char);
    mesh->info.noinsert  = *( (char *) *buffer); *buffer += sizeof(char);
    mesh->info.noswap    = *( (char *) *buffer); *buffer += sizeof(char);
    mesh->info.nomove    = *( (char *) *buffer); *buffer += sizeof(char);
    mesh->info.nosurf    = *( (char *) *buffer); *buffer += sizeof(char);
    mesh->info.inputMet  = *( (char *) *buffer); *buffer += sizeof(char);

    /* affectation of old refs in ls-mode */
    if ( mesh->info.nmat ) {

      MMG5_SAFE_CALLOC(mesh->info.mat,mesh->info.nmat,MMG5_Mat, ier = 0);

      if ( ier ) {
        for ( k=0; k<mesh->info.nmat; ++k ) {
          mesh->info.mat[k].dospl = *( (char *) *buffer); *buffer += sizeof(char);
          mesh->info.mat[k].ref   = *( (int *) *buffer); *buffer += sizeof(int);
          mesh->info.mat[k].rin   = *( (int *) *buffer); *buffer += sizeof(int);
          mesh->info.mat[k].rex   = *( (int *) *buffer); *buffer += sizeof(int);
        }
      }
    }

    /* local parameters */
    if ( mesh->info.npar ) {

      MMG5_SAFE_CALLOC(mesh->info.par,mesh->info.npar,MMG5_Par, ier = 0);

      if ( ier ) {
        for ( k=0; k<mesh->info.npar; ++k ) {
          mesh->info.par[k].hmin = *( (double *) *buffer); *buffer += sizeof(double);
          mesh->info.par[k].hmax = *( (double *) *buffer); *buffer += sizeof(double);
          mesh->info.par[k].hausd = *( (double *) *buffer); *buffer += sizeof(double);
          mesh->info.par[k].ref = *( (int *) *buffer); *buffer += sizeof(int);
          mesh->info.par[k].elt = *( (char *) *buffer); *buffer += sizeof(char);
        }
      }
    }
  }
  else {
    ier = 0;

    /** Mesh infos */
    ddummy = *( (double *) *buffer); *buffer += sizeof(double);
    ddummy = *( (double *) *buffer); *buffer += sizeof(double);
    ddummy = *( (double *) *buffer); *buffer += sizeof(double);
    ddummy = *( (double *) *buffer); *buffer += sizeof(double);
    ddummy = *( (double *) *buffer); *buffer += sizeof(double);
    ddummy = *( (double *) *buffer); *buffer += sizeof(double);
    ddummy = *( (double *) *buffer); *buffer += sizeof(double);

    ddummy = *( (double *) *buffer); *buffer += sizeof(double);
    ddummy = *( (double *) *buffer); *buffer += sizeof(double);
    ddummy = *( (double *) *buffer); *buffer += sizeof(double);
    ddummy = *( (double *) *buffer); *buffer += sizeof(double);

    ddummy = *( (double *) *buffer); *buffer += sizeof(double);

    idummy = *( (int *) *buffer); *buffer += sizeof(int);
    idummy = *( (int *) *buffer); *buffer += sizeof(int);
    idummy = *( (int *) *buffer); *buffer += sizeof(int);
    idummy = *( (int *) *buffer); *buffer += sizeof(int);
    idummy = *( (int *) *buffer); *buffer += sizeof(int);

    cdummy = *( (char *) *buffer); *buffer += sizeof(char);
    cdummy = *( (char *) *buffer); *buffer += sizeof(char);
    cdummy = *( (char *) *buffer); *buffer += sizeof(char);
    cdummy = *( (char *) *buffer); *buffer += sizeof(char);
    cdummy = *( (char *) *buffer); *buffer += sizeof(char);
    cdummy = *( (char *) *buffer); *buffer += sizeof(char);
    cdummy = *( (char *) *buffer); *buffer += sizeof(char);
    cdummy = *( (char *) *buffer); *buffer += sizeof(char);
    cdummy = *( (char *) *buffer); *buffer += sizeof(char);
    cdummy = *( (char *) *buffer); *buffer += sizeof(char);
    cdummy = *( (char *) *buffer); *buffer += sizeof(char);
    cdummy = *( (char *) *buffer); *buffer += sizeof(char);
    cdummy = *( (char *) *buffer); *buffer += sizeof(char);

    if ( mesh->info.nmat ) {
      for ( k=0; k<mesh->info.nmat; ++k ) {
        cdummy = *( (char *) *buffer); *buffer += sizeof(char);
        idummy = *( (int *) *buffer); *buffer += sizeof(int);
        idummy = *( (int *) *buffer); *buffer += sizeof(int);
        idummy = *( (int *) *buffer); *buffer += sizeof(int);
      }
    }

    /* local parameters */
    if ( mesh->info.npar ) {
      for ( k=0; k<mesh->info.npar; ++k ) {
        ddummy = *( (double *) *buffer); *buffer += sizeof(double);
        ddummy = *( (double *) *buffer); *buffer += sizeof(double);
        ddummy = *( (double *) *buffer); *buffer += sizeof(double);
        idummy = *( (int *) *buffer);    *buffer += sizeof(int);
        cdummy = *( (char *) *buffer);   *buffer += sizeof(char);
      }
    }
  }

  /** Metric info and sizes */
  ismet     = *( (int *) *buffer); *buffer += sizeof(int);
  size      = *( (int *) *buffer); *buffer += sizeof(int);

  np = ismet ? mesh->np : 0;

  if ( ier_grp ) {
    met->size = size;
    /** Metric type */
    met->type = *( (int *) *buffer); *buffer += sizeof(int);

    /** Set the metric size */
    ier_sol = MMG3D_Set_solSize(mesh,met,MMG5_Vertex,np,met->type);

    /* Use exactly the amount of needed memory for this mesh and metric */
    mesh->memMax = mesh->memCur;

    /* Update the available memory count */
    *memAv -= mesh->memMax;

    ier = MG_MIN ( ier, ier_sol );
  }
  else {
    ier = ier_sol = 0;
    /** Metric type */
    idummy = *( (int *) *buffer); *buffer += sizeof(int);
  }

  /** Communicator sizes */
  grp->nitem_int_node_comm = *( (int *) *buffer); *buffer += sizeof(int);
  grp->nitem_int_face_comm = *( (int *) *buffer); *buffer += sizeof(int);

  if ( ier_mesh ) {
    /** Get mesh points */
    for ( k=1; k<=mesh->np; ++k ) {
      /* Coordinates */
      mesh->point[k].c[0] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->point[k].c[1] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->point[k].c[2] = *( (double *) *buffer); *buffer += sizeof(double);
      /* Tangent */
      mesh->point[k].n[0] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->point[k].n[1] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->point[k].n[2] = *( (double *) *buffer); *buffer += sizeof(double);
      /* Pointer toward the boundary entity */
      mesh->point[k].xp = *( (int *) *buffer); *buffer += sizeof(int);
      /* Ref */
      mesh->point[k].ref = *( (int *) *buffer); *buffer += sizeof(int);
      /* Tag */
      mesh->point[k].tag = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
    }

    /** Pack mesh boundary points */
    for ( k=1; k<=mesh->xp; ++k ) {
      /* First normal */
      mesh->xpoint[k].n1[0] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->xpoint[k].n1[1] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->xpoint[k].n1[2] = *( (double *) *buffer); *buffer += sizeof(double);
      /* Second normal */
      mesh->xpoint[k].n2[0] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->xpoint[k].n2[1] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->xpoint[k].n2[2] = *( (double *) *buffer); *buffer += sizeof(double);
    }

    /** Pack mesh elements */
    for ( k=1; k<=mesh->ne; ++k ) {
      /* Tetra vertices */
      mesh->tetra[k].v[0] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->tetra[k].v[1] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->tetra[k].v[2] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->tetra[k].v[3] = *( (int *) *buffer); *buffer += sizeof(int);
      /* Pointer toward the boundary entity */
      mesh->tetra[k].xt = *( (int *) *buffer); *buffer += sizeof(int);
      /* Ref */
      mesh->tetra[k].ref = *( (int *) *buffer); *buffer += sizeof(int);
      /* Mark */
      mesh->tetra[k].mark = *( (int *) *buffer); *buffer += sizeof(int);
      /* Tag */
      mesh->tetra[k].tag = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      /* Quality */
      mesh->tetra[k].qual = *( (double *) *buffer); *buffer += sizeof(double);
    }

    /** Pack mesh boundary tetra */
    for ( k=1; k<=mesh->xt; ++k ) {
      /* Faces references  */
      mesh->xtetra[k].ref[0] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->xtetra[k].ref[1] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->xtetra[k].ref[2] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->xtetra[k].ref[3] = *( (int *) *buffer); *buffer += sizeof(int);
      /* Faces tags */
      mesh->xtetra[k].ftag[0] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      mesh->xtetra[k].ftag[1] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      mesh->xtetra[k].ftag[2] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      mesh->xtetra[k].ftag[3] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      /* Edges references  */
      mesh->xtetra[k].edg[0] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->xtetra[k].edg[1] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->xtetra[k].edg[2] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->xtetra[k].edg[3] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->xtetra[k].edg[4] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->xtetra[k].edg[5] = *( (int *) *buffer); *buffer += sizeof(int);
      /* Edges tags */
      mesh->xtetra[k].tag[0] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      mesh->xtetra[k].tag[1] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      mesh->xtetra[k].tag[2] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      mesh->xtetra[k].tag[3] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      mesh->xtetra[k].tag[4] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      mesh->xtetra[k].tag[5] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
    }
  }
  else {
    /* The mesh can't be allocated */
   /** Get mesh points */
    for ( k=1; k<=mesh->np; ++k ) {
      /* Coordinates */
      ddummy = *( (double *) *buffer); *buffer += sizeof(double);
      ddummy = *( (double *) *buffer); *buffer += sizeof(double);
      ddummy = *( (double *) *buffer); *buffer += sizeof(double);
      /* Tangent */
      ddummy = *( (double *) *buffer); *buffer += sizeof(double);
      ddummy = *( (double *) *buffer); *buffer += sizeof(double);
      ddummy = *( (double *) *buffer); *buffer += sizeof(double);
      /* Pointer toward the boundary entity */
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      /* Ref */
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      /* Tag */
      i16dummy = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
    }

    /** Pack mesh boundary points */
    for ( k=1; k<=mesh->xp; ++k ) {
      /* First normal */
      ddummy = *( (double *) *buffer); *buffer += sizeof(double);
      ddummy = *( (double *) *buffer); *buffer += sizeof(double);
      ddummy = *( (double *) *buffer); *buffer += sizeof(double);
      /* Second normal */
      ddummy = *( (double *) *buffer); *buffer += sizeof(double);
      ddummy = *( (double *) *buffer); *buffer += sizeof(double);
      ddummy = *( (double *) *buffer); *buffer += sizeof(double);
    }

    /** Pack mesh elements */
    for ( k=1; k<=mesh->ne; ++k ) {
      /* Tetra vertices */
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      /* Pointer toward the boundary entity */
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      /* Ref */
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      /* Tag */
      i16dummy = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
    }

    /** Pack mesh boundary tetra */
    for ( k=1; k<=mesh->xt; ++k ) {
      /* Faces references  */
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      /* Faces tags */
      i16dummy = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      i16dummy = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      i16dummy = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      i16dummy = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      /* Edges references  */
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      idummy = *( (int *) *buffer); *buffer += sizeof(int);
      /* Edges tags */
      i16dummy = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      i16dummy = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      i16dummy = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      i16dummy = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      i16dummy = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      i16dummy = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
    }
  }

  /** Pack metric */
  if ( ier_sol ) {
    for ( k=1; k<=mesh->np; ++k ) {
      for ( i=0; i<size; ++i ) {
        met->m[size*k + i] = *( (double *) *buffer);
        *buffer += sizeof(double);
      }
    }
  }
  else {
    /* The solution array can't be allocated */
    for ( k=1; k<=mesh->np; ++k ) {
      for ( i=0; i<size; ++i ) {
        ddummy = *( (double *) *buffer);
        *buffer += sizeof(double);
      }
    }
  }

  /** Pack communicators */
  ier_comm = 1;

  /* Give all the available mem to the communicators */
  parmesh->memMax = *memAv;

  /* Node communicator */
  PMMG_MALLOC(parmesh,grp->node2int_node_comm_index1,grp->nitem_int_node_comm,
              int,"node2int_node_comm_index1",ier = ier_comm = 0);
  PMMG_MALLOC(parmesh,grp->node2int_node_comm_index2,grp->nitem_int_node_comm,
              int,"node2int_node_comm_index1",ier = ier_comm = 0);

  if ( ier_comm ) {
    for ( k=0; k<grp->nitem_int_node_comm; ++k ) {
      grp->node2int_node_comm_index1[k] = *( (int *) *buffer);
      *buffer += sizeof(int);
    }
    for ( k=0; k<grp->nitem_int_node_comm; ++k ) {
      grp->node2int_node_comm_index2[k] = *( (int *) *buffer);
      *buffer += sizeof(int);
    }
  }
  else {
    /* The communicators arrays can't be allocated */
    for ( k=0; k<2*grp->nitem_int_node_comm; ++k ) {
      idummy = *( (int *) *buffer);
      *buffer += sizeof(int);
    }
  }

  /* Face communicator */
  PMMG_MALLOC(parmesh,grp->face2int_face_comm_index1,grp->nitem_int_face_comm,
              int,"face2int_face_comm_index1",ier = ier_comm = 0);
  PMMG_MALLOC(parmesh,grp->face2int_face_comm_index2,grp->nitem_int_face_comm,
              int,"face2int_face_comm_index1",ier = ier_comm = 0);

  /* Use the minimal memory needed */
  parmesh->memMax = parmesh->memCur;
  *memAv         -= parmesh->memMax;

  if ( ier_comm ) {
    for ( k=0; k<grp->nitem_int_face_comm; ++k ) {
      grp->face2int_face_comm_index1[k] = *( (int *) *buffer);
      *buffer += sizeof(int);
    }
    for ( k=0; k<grp->nitem_int_face_comm; ++k ) {
      grp->face2int_face_comm_index2[k] = *( (int *) *buffer);
      *buffer += sizeof(int);
    }
  }
  else {
    /* The communicators arrays can't be allocated */
    for ( k=0; k<2*grp->nitem_int_face_comm; ++k ) {
      idummy = *( (int *) *buffer);
      *buffer += sizeof(int);
    }
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
  if ( ier0 )
    for ( k=0; k<nprocs; ++k ) {
      (*trequest)[k] = MPI_REQUEST_NULL;
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
  size_t         available;
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
  PMMG_RECALLOC ( parmesh,parmesh->listgrp,ngrp+grpscount,ngrp,PMMG_Grp,"listgrp",
                  ier0 = 0;ier = 0 );

  if ( ier0 )
    parmesh->ngrp += grpscount;

    /* Count the available memory and update it while unpacking the received groups */
#ifndef NDEBUG
  for ( k=0; k<ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    assert ( mesh->memCur == mesh->memMax );
    assert ( mesh->npmax == mesh->np );
    assert ( mesh->xpmax == mesh->xp );
    assert ( mesh->nemax == mesh->ne );
    assert ( mesh->xtmax == mesh->xt );
  }
#endif

  parmesh->memMax = parmesh->memCur;
  available = parmesh->memGloMax - parmesh->memMax;
  for ( k=0; k<ngrp; ++k ) {
    available -= parmesh->listgrp[k].mesh->memMax;
  }
  assert ( available >= 0 );

  ptr = buffer;
  if ( ier0 ) {
    for ( k=0; k<grpscount; ++k ) {
      err = PMMG_mpiunpack_grp(parmesh,&parmesh->listgrp[ngrp+k],&ptr,&available);
      ier = MG_MIN(ier,err);
      parmesh->listgrp[ngrp+k].flag = PMMG_UNSET;
    }

    /* Give the remaining memory to the parmesh */
    parmesh->memMax += available;
  }

  PMMG_DEL_MEM ( parmesh,buffer,char,"buffer" );
  return ier;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param sndr index of the proc that send the groups
 * \param recv index of the proc that receive the groups
 * \param interaction_map map of interactions with the other processors
 *
 * \return 0 if fail, 1 if we success
 *
 * Transfer and update the data that are modified due to the transfer of the
 * groups from the proc \a sndr toward the proc \a recv.
 *
 */
static inline
int PMMG_transfer_grps_fromItoJ(PMMG_pParMesh parmesh,const int sndr,
                                const int recv,int *interaction_map) {

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

  /** Step 3: pack the groups and face communicators */
  err = PMMG_pack_grps( parmesh,&parmesh->listgrp);
  ier = MG_MIN ( ier, err );
  if ( (!parmesh->ngrp) && parmesh->ddebug ) {
    fprintf(stderr,"  ## Warning: %s: rank %d: processor without any groups.\n",
            __func__,myrank);
  }

  if ( ier ) {
    /* Pack the face communicators */
    err = PMMG_pack_faceCommunicators(parmesh);
    ier = MG_MIN ( ier, err );
  }

  /** Step 4: Wait for the end of the MPI communications */
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
 *
 * \return 0 if fail but we can try to save a mesh, -1 if we fail and are unable
 * to save the mesh, 1 if we success
 *
 * Send the suitable groups to other procs and recieve their groups.
 * Deallocate the \a part array.
 *
 */
static inline
int PMMG_transfer_all_grps(PMMG_pParMesh parmesh,idx_t *part) {
  MPI_Comm       comm;
  int            myrank,nprocs;
  int            *interaction_map;
  int            *send_grps,*recv_grps,*nfaces2send,*nfaces2recv;
  int            *next_comm2send,*ext_comms_next_idx,*nitems2send;
  int            *extComm_next_idx,*items_next_idx,*recv_array;
  int            *extComm_grpFaces2extComm,*extComm_grpFaces2face2int;
  int            max_ngrp;
  int            ier,ier_glob,k,j,err;

  myrank    = parmesh->myrank;
  nprocs    = parmesh->nprocs;
  comm      = parmesh->comm;
  max_ngrp  = 0;

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

  /** Step 1: Merge all the groups that must be sended to a given proc into 1
   * group */
  ier = PMMG_merge_grps2send(parmesh,&part);

  /* Give all the available memory to the parmesh */
  assert ( parmesh->memGloMax > parmesh->memCur );

  parmesh->memMax = parmesh->memGloMax;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    err = PMMG_parmesh_fitMesh( parmesh,parmesh->listgrp[k].mesh,
                                parmesh->listgrp[k].met );
    ier = MG_MIN( err, ier );

    if ( parmesh->memMax < parmesh->listgrp[k].mesh->memMax ) {
      fprintf(stderr,"\n  ## Error: %s: not enough memory\n",__func__);
      ier = MG_MIN(ier,-1);
    }
    parmesh->memMax -= parmesh->listgrp[k].mesh->memMax;
  }

  /* Store the group destination in the group flag */
  for ( k=0; k<parmesh->ngrp; ++k ) {
    parmesh->listgrp[k].flag = part[k];
  }
  PMMG_DEL_MEM(parmesh,part,idx_t,"parmetis partition");

  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, comm);

  if ( ier_glob < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to compute the new group"
            " partition.\n",__func__);
    ier = -1;
    goto end;
  }
  else if ( !ier_glob ) {
    fprintf(stderr,"\n  ## Warning: %s: unable to compute the new group"
            " partition. Try to send it nevertheless.\n",__func__);
  }

  /** Step 2: Nodal communicators deletion (unallocation to remove, at least for
   * the external communicators) */
  PMMG_node_comm_free(parmesh);

  /** Step 3:
   *
   * Compute the map of interactions: each proc interacts with all the other
   * procs (worst case)
   *
   */
#warning add interaction map
  ier = 0;
  PMMG_CALLOC(parmesh,interaction_map,nprocs,int,"interaction_map",goto end);
  for ( k=0; k<nprocs; ++k ) {
    interaction_map[k] = 1;
  }

  /** Step 4: proc k send its data (group and/or communicators), proc j receive
   * data */
  ier = 1;
  for ( k=0; k<nprocs; ++k ) {
    for ( j=0; j<nprocs; ++j ) {
      if ( j==k ) continue;
      err =  PMMG_transfer_grps_fromItoJ(parmesh,k,j,interaction_map);
      ier = MG_MIN ( ier,err );
    }
  }
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, comm);

  if ( ier_glob <= 0 ) {
    fprintf(stderr,"\n  ## Error: %s: unable to transfer the groups toward"
            " their new processor.\n",__func__);
    ier = -1;
    goto end;
  }

  /** Step 6: Node communicators reconstruction from the face ones */
  if ( !PMMG_build_nodeCommFromFaces(parmesh) ) {
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
 * \param moveIfcs 0 for metis, 1 for advancing-front method.
 *
 * \return -1 if we fail and can not save the meshes, 0 if we fail but can save
 * the meshes, 1 otherwise
 *
 * Compute the group distribution that load balance the groups over the
 * processors and send and recieve the groups from the other processors.
 *
 */
int PMMG_distribute_grps( PMMG_pParMesh parmesh,int moveIfcs ) {
  idx_t *part;
  int   ngrp,ier,ier_glob;

  MPI_Allreduce( &parmesh->ngrp, &ngrp, 1, MPI_INT, MPI_MIN, parmesh->comm);

  if ( !ngrp ) {
    fprintf(stderr,"Error:%s:%d: Empty partition. Not yet implemented\n",
            __func__,__LINE__);
    return 0;
  }

  /** Get the new partition of groups (1 group = 1 metis node) */
  part = NULL;
  PMMG_CALLOC(parmesh,part,parmesh->ngrp,idx_t,"allocate parmetis buffer",
              return 0);

  if( moveIfcs ) {

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

  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);

  if ( !ier_glob ) {
    PMMG_DEL_MEM(parmesh,part,idx_t,"deallocate parmetis partition");
    return ier_glob;
  }

  /** Send the suitable groups to other procs and recieve their groups */
  ier = PMMG_transfer_all_grps(parmesh,part);
  if ( ier <= 0 )
    fprintf(stderr,"\n  ## Unable to communicate groups through processors.\n");

  /** Check grps contiguity */
  ier = PMMG_checkAndReset_grps_contiguity( parmesh );

  return ier;
}

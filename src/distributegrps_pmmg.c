/**
 * \file distributegrps_pmmg.c
 * \brief Group distribution on the processors
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \author Nikos Pattakos (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
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
 * \param parmesh pointer toward the parmesh structure.
 * \param grpI pointer toward the group in which we want to merge.
 * \param grpJ pointer toward the group that we want to merge.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Update the node communicators when merging the group \a grpJ into the group
 * \a grpI.
 *
 */
static inline
int PMMG_mergeGrpJinI_nodeCommunicators( PMMG_pParMesh parmesh,PMMG_pGrp grpI,
                                         PMMG_pGrp grpJ ) {
  PMMG_pext_comm ext_node_comm;
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

  /** Step 2: Process the interface points of the group \a grpI and delete the
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

  /** Step 3: Process the interface points of the group \a grpJ and add the
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

    assert ( MG_EOK(&grpI->mesh->tetra[iel/4]) );
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
                                    PMMG_pGrp grpJ) {

  if ( !PMMG_mergeGrpJinI_nodeCommunicators(parmesh,grpI,grpJ) ) return 0;

  if ( !PMMG_mergeGrpJinI_faceCommunicators(parmesh,grpI,grpJ) ) return 0;

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
  MMG5_pSol      metI,metJ;
  int            *intvalues;
  int            poi_id_int,poi_id_glo,k;

  meshI     = grpI->mesh;
  metI      = grpI->met;

  meshJ     = grpJ->mesh;
  metJ      = grpJ->met;

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

  PMMG_pint_comm int_node_comm,int_face_comm;
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
  if ( !PMMG_mergeGrpJinI_internalPoints(parmesh,grpI,grpJ) ) return 0;

  /** Step 3: Merge interfaces tetras of the grpJ into the grpI mesh */
  /* 1) Store the indices of the interface faces of grpI into the internal
   * face communicator */
  face2int_face_comm_index1 = grpI->face2int_face_comm_index1;
  face2int_face_comm_index2 = grpI->face2int_face_comm_index2;

  for ( k=0; k<grpI->nitem_int_face_comm; ++k ) {
    iel = face2int_face_comm_index1[k];
    assert(   ( 0 <= face2int_face_comm_index2[k] )
              && ( face2int_face_comm_index2[k] < parmesh->int_node_comm->nitem )
              && "check intvalues indices" );
    parmesh->int_face_comm->intvalues[face2int_face_comm_index2[k]] = iel;
  }

  /* 2) Add the interfaces tetra of the grpJ mesh to the meshI mesh */
  if ( !PMMG_mergeGrpJinI_interfaceTetra(parmesh,grpI,grpJ) ) return 0;

  /** Step 4: Merge internal tetras of the grpJ mesh mesh into the meshI mesh */
  if ( !PMMG_mergeGrpJinI_internalTetra(parmesh,grpI,grpJ) ) return 0;

  /** Step 5: Update the communicators */
  if ( !PMMG_mergeGrpJinI_communicators(parmesh,grpI,grpJ) ) return 0;

  return 1;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param grps pointer toward the list of groups to pack.
 * \param ngrp pointer toward the number of groups (to update)
 * \param part pointer toward the metis array containing the partitions to
 * update and pack.
 *
 * \return 0 if fail, 1 if success.
 *
 * Pack the group array \a grps filled such as grps[k] must be send to the proc
 * k and update the parmetis partition array.
 *
 */
static inline
int PMMG_pack_grpsAndPart( PMMG_pParMesh parmesh,PMMG_pGrp *grps,
                            int *ngrp,idx_t **part ) {
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
  PMMG_REALLOC( parmesh,*part,nbl,*ngrp,idx_t,"parmetis partition",return 0 );

  *ngrp = nbl;

  return 1;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param grps pointer toward the list of groups to pack.
 * \param ngrp pointer toward the number of groups (to update)
 *
 * \return the new number of groups
 *
 * Pack the group array \a grps. A group is unused if its mesh is NULL.
 *
 */
static inline
int PMMG_pack_grps_norealloc( PMMG_pParMesh parmesh,PMMG_pGrp *grps,
                              int ngrp ) {
  int k,nbl;

  nbl   = 0;
  for ( k=0; k<ngrp; ++k ) {
    if ( !(*grps)[k].mesh ) continue;

    if ( k!=nbl )
      (*grps)[nbl] = PMMG_move_grp( &(*grps)[k] );

    ++nbl;
  }

  return nbl;
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
static inline
int PMMG_pack_nodeCommunicators(PMMG_pParMesh parmesh) {
  PMMG_pGrp      grp;
  PMMG_pint_comm int_node_comm;
  PMMG_pext_comm ext_node_comm;
  int            *intvalues;
  int            *node2int_node_comm_index2;
  int            k,nitem,idx,i;

  intvalues = parmesh->int_node_comm->intvalues;

  /** Step 1: initialization of the communicator */
  int_node_comm = parmesh->int_node_comm;
  for ( k=0; k<int_node_comm->nitem; ++k )
    intvalues[k] = PMMG_UNSET;

  /** Step 2: Process the interface points of each group and store in intvalue
   * their new position in the internal communicator */
  nitem = 0;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp = &parmesh->listgrp[k];

    for ( i=0; i<grp->nitem_int_node_comm; ++i )
      if ( intvalues[grp->node2int_node_comm_index2[i]]<0 )
        intvalues[grp->node2int_node_comm_index2[i]] = nitem++;
  }

  /** Step 3: Update the group communicators arrays with the new points
   * position */
  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp                       = &parmesh->listgrp[k];
    node2int_node_comm_index2 = grp->node2int_node_comm_index2;

    for ( i=0; i<grp->nitem_int_node_comm; ++i ) {
      idx = intvalues[node2int_node_comm_index2[i]];
      assert ( idx>=0 && "Deleted position in internal communicator" );

      node2int_node_comm_index2[i] = idx;
    }
  }

  /** Step 4: Update the external communicators arrays with the new points
   * position */
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];

    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx = intvalues[ext_node_comm->int_comm_index[i]];
      assert ( idx>=0 && "Deleted position in internal communicator" );

      ext_node_comm->int_comm_index[i] = idx;
    }
  }

  /** Step 5: unallocate intvalues array and set the nitem field of the internal
   * communicator to the suitable value */
  PMMG_DEL_MEM( parmesh,int_node_comm->intvalues,int_node_comm->nitem,int,
                "node communicator");
  int_node_comm->nitem = nitem;

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
static inline
int PMMG_pack_faceCommunicators(PMMG_pParMesh parmesh) {
  PMMG_pGrp      grp;
  PMMG_pint_comm int_face_comm;
  PMMG_pext_comm ext_face_comm;
  int            *intvalues;
  int            *face2int_face_comm_index2;
  int            k,nitem,idx,i;

  intvalues = parmesh->int_face_comm->intvalues;

  /** Step 1: initialization of the communicator */
  int_face_comm = parmesh->int_face_comm;
  for ( k=0; k<int_face_comm->nitem; ++k )
    intvalues[k] = PMMG_UNSET;

  /** Step 2: Process the interface points of each group and store in intvalue
   * their new position in the internal communicator */
  nitem = 0;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp = &parmesh->listgrp[k];

    for ( i=0; i<grp->nitem_int_face_comm; ++i )
      if ( intvalues[grp->face2int_face_comm_index2[i]]<0 )
        intvalues[grp->face2int_face_comm_index2[i]] = nitem++;
  }

  /** Step 3: Update the group communicators arrays with the new points
   * position */
  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp                       = &parmesh->listgrp[k];
    face2int_face_comm_index2 = grp->face2int_face_comm_index2;

    for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
      idx = intvalues[face2int_face_comm_index2[i]];
      assert ( idx>=0 && "Deleted position in internal communicator" );

      face2int_face_comm_index2[i] = idx;
    }
  }

  /** Step 4: Update the external communicators arrays with the new points
   * position */
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];

    for ( i=0; i<ext_face_comm->nitem; ++i ) {
      idx = intvalues[ext_face_comm->int_comm_index[i]];
      assert ( idx>=0 && "Deleted position in internal communicator" );

      ext_face_comm->int_comm_index[i] = idx;
    }
  }

  /** Step 5: unallocate intvalues array and set the nitem field of the internal
   * communicator to the suitable value */
  PMMG_DEL_MEM( parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,
                "face communicator");
  int_face_comm->nitem = nitem;

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
  ngrpOld  = ngrp;
  ngrp    = k;
  PMMG_REALLOC( parmesh,*grps,ngrp,ngrpOld,PMMG_Grp,"Groups to send",
                return 0;);
  PMMG_REALLOC( parmesh,tmpPart,ngrp,ngrpOld,idx_t,"parmetis partition",
                return 0;);
  PMMG_DEL_MEM( parmesh,*part,nprocs,idx_t,"parmetis partition" );
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
  PMMG_int_comm *int_node_comm,*int_face_comm;
  MMG5_pMesh    meshI,meshJ;
  long long     memAv;
  int           nprocs,ngrp,k,j,ier;

  nprocs = parmesh->nprocs;
  ngrp   = parmesh->ngrp;

  /** Step 0: Set mesh sizes to their minimal size and count the available
   * memory */
  memAv = parmesh->memGloMax - parmesh->memMax;
  for ( k=0; k<ngrp; ++k ) {
    parmesh->listgrp[k].mesh->memMax = parmesh->listgrp[k].mesh->memCur;
    memAv -= parmesh->listgrp[k].mesh->memMax;
  }
  assert ( memAv >=0 );

  /** Step 1: New groups allocation and initialization: move the groups to have
   * a group that will be send to proc k stored in grps[k]. Free the adja
   * array. */

  // If this step is too expensive in memory, we can count the number of procs
  // with which we will communicate and fill directly the pack array.
  PMMG_CALLOC( parmesh,grps,nprocs,PMMG_Grp,"Groups to send",return 0 );

  j = 0;
  for ( k=0; k<ngrp; ++k ) {
    /* Free the adja array */
    meshI = parmesh->listgrp[k].mesh;
    if ( meshI->adja )
      PMMG_DEL_MEM(meshI, meshI->adja, 4*meshI->nemax+5, int, "adjacency table" );

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
      PMMG_DEL_MEM(meshI, meshI->adja, 4*meshI->nemax+5, int, "adjacency table" );
    if ( meshJ->adja )
      PMMG_DEL_MEM(meshJ, meshJ->adja, 4*meshJ->nemax+5, int, "adjacency table" );

    meshI->memMax = memAv;
    if ( !PMMG_merge_grpJinI(parmesh,grpI,grpJ) ) goto low_fail;

    memAv -= meshI->memCur;
    meshI->memMax = meshI->memCur;

    /* Delete the useless group to gain memory space */
    memAv += meshJ->memCur;
    PMMG_grp_free(parmesh,&listgrp[k]);
  }
  ngrp = nprocs;

  if ( !PMMG_pack_grpsAndPart( parmesh,&grps,&ngrp,part ) ) ier = -1;
  else ier  = 1;
  goto end;

strong_fail2:
  PMMG_DEL_MEM( parmesh,int_node_comm->intvalues,int_node_comm->nitem,int,
                "node communicator");
strong_fail1:
  PMMG_DEL_MEM( parmesh,grps,nprocs,PMMG_Grp,"Groups to send");
  return -1;

low_fail:
  ngrp = PMMG_mergeGrps2Send_errorHandler(parmesh,&grps,part,k);
  if ( !ngrp ) ier=-1;
  else ier = 0;

end:

  /** Step 3: Update the parmesh */
  PMMG_DEL_MEM(parmesh,parmesh->listgrp,parmesh->ngrp,PMMG_Grp,"listgrp");
  parmesh->listgrp = grps;
  parmesh->ngrp    = ngrp;

  /* Pack the communicators */
  if ( !PMMG_pack_nodeCommunicators(parmesh) ) ier = -1;

  if ( !PMMG_pack_faceCommunicators(parmesh) ) ier = -1;

  if ( !PMMG_parmesh_updateMemMax(parmesh, 105, 1) ) {
    fprintf(stderr,"\n  ## Error: %s: Unable to update the memory repartition"
            " between meshes and communicators.\n",__func__);
    return -1;
  }

  return ier;
}

/**
 * \param mesh pointer toward a mesh structure.
 * \param met pointer toward a metric structure.
 * \param dest rank of destination.
 * \param tag message tag
 * \param comm communicator
 *
 * \return 0 if fail, 1 if success.
 *
 * Send mesh and metric to proc \a dest.
 *
 * \remark the prisms are not treated.
 *
 */
int PMMG_send_mesh( MMG5_pMesh mesh,MMG5_pSol met,int dest,int tag,MPI_Comm comm)
{
  MPI_Datatype   mpi_point,mpi_xpoint,mpi_tetra,mpi_xtetra;
  int            ier,ismet;

  /* Fit the mesh and metric at their minimal sizes */
  PMMG_REALLOC(mesh,mesh->point,mesh->np+1,mesh->npmax+1,MMG5_Point,
               "points",return 0);
  mesh->npmax = mesh->np;
  PMMG_REALLOC(mesh,mesh->xpoint,mesh->xp+1,mesh->xpmax+1,MMG5_xPoint,
               "xpoints",return 0);
  mesh->xpmax = mesh->xp;
  PMMG_REALLOC(mesh,mesh->tetra,mesh->ne+1,mesh->nemax+1,MMG5_Tetra,
               "tetra",return 0);
  mesh->nemax = mesh->ne;
  PMMG_REALLOC(mesh,mesh->xtetra,mesh->xt+1,mesh->xpmax+1,MMG5_xTetra,
               "xtetra",return 0);
  mesh->xtmax = mesh->xt;

  assert(met->npmax==mesh->npmax);
  assert(met->np   == met->np);
  PMMG_REALLOC(mesh,met->m,met->size*met->np+1,met->size*met->npmax+1,double,
               "metric array",return 0);
  met->npmax = met->np;

  /* Send the mesh size */
  mesh->memMax = mesh->memCur;
  MPI_CHECK( MPI_Send( &mesh->memMax, 1, MPI_LONG_LONG, dest, tag, comm ),
             return 0 );
  MPI_CHECK( MPI_Send( &mesh->np, 1, MPI_INT, dest, tag, comm ), return 0 );
  MPI_CHECK( MPI_Send( &mesh->ne, 1, MPI_INT, dest, tag, comm ), return 0 );
  MPI_CHECK( MPI_Send( &mesh->xp, 1, MPI_INT, dest, tag, comm ), return 0 );
  MPI_CHECK( MPI_Send( &mesh->xt, 1, MPI_INT, dest, tag, comm ), return 0 );

  ier = 0;
  /* Create mesh type */
  if ( !PMMG_create_MPI_Point( mesh->point,  &mpi_point ) ) goto end_5;
  if ( !PMMG_create_MPI_xPoint(mesh->xpoint, &mpi_xpoint) ) goto end_4;
  if ( !PMMG_create_MPI_Tetra( mesh->tetra,  &mpi_tetra ) ) goto end_3;
  if ( !PMMG_create_MPI_xTetra(mesh->xtetra, &mpi_xtetra) ) goto end_2;

  /* Send mesh */
  MPI_CHECK( MPI_Send( &mesh->point[1] ,mesh->np,mpi_point, dest,tag,comm),
             goto end_1;);
  MPI_CHECK( MPI_Send( &mesh->xpoint[1],mesh->xp,mpi_xpoint,dest,tag,comm),
             goto end_1;);
  MPI_CHECK( MPI_Send( &mesh->tetra[1] ,mesh->ne,mpi_tetra, dest,tag,comm),
             goto end_1;);
  MPI_CHECK( MPI_Send( &mesh->xtetra[1],mesh->xt,mpi_xtetra,dest,tag,comm),
             goto end_1;);

  /* Met size */
  assert ( met->np == mesh->np );
  MPI_CHECK( MPI_Send( &met->size, 1, MPI_INT, dest, tag, comm ),goto end_1 );
  MPI_CHECK( MPI_Send( &met->type, 1, MPI_INT, dest, tag, comm ),goto end_1 );

  /* Send metric */
  ismet = met->m ? 1:0;
  MPI_CHECK( MPI_Send( &ismet, 1, MPI_INT, dest, tag, comm ),
             goto end_1 );

  if ( ismet )
    MPI_CHECK( MPI_Send( met->m, met->size*met->np, MPI_DOUBLE, dest, tag, comm ),
               goto end_1 );

  ier = 1;

  /* Free useless mesh */
  /* 1: mesh */
  PMMG_DEL_MEM(mesh,mesh->point, mesh->npmax+1,sizeof(MMG5_Point),"point");
  PMMG_DEL_MEM(mesh,mesh->tetra, mesh->nemax+1,sizeof(MMG5_Tetra),"tetra");
  PMMG_DEL_MEM(mesh,mesh->xpoint,mesh->xpmax+1,sizeof(MMG5_xPoint),"xpoint");
  PMMG_DEL_MEM(mesh,mesh->xtetra,mesh->xtmax+1,sizeof(MMG5_xTetra),"xtetra");
  if ( mesh->adja ) {
    PMMG_DEL_MEM( mesh,mesh->adja,4*mesh->nemax+5,int,"adjacency table");
  }
  if ( ismet )
    PMMG_DEL_MEM(mesh,met->m,(met->npmax+1)*met->size,double,"metric");

  mesh->np = mesh->npmax = 0;
  mesh->xp = mesh->xpmax = 0;
  mesh->ne = mesh->nemax = 0;
  mesh->xt = mesh->xtmax = 0;

end_1:
  MPI_Type_free( &mpi_xtetra );
end_2:
  MPI_Type_free( &mpi_tetra );
end_3:
  MPI_Type_free( &mpi_xpoint );
end_4:
  MPI_Type_free( &mpi_point );
end_5:
  return ier;
}

/**
 * \param mesh pointer toward a mesh structure.
 * \param source rank of source.
 * \param tag message tag
 * \param comm communicator
 *
 * \return 0 if fail, 1 if success.
 *
 * Receive mesh and metric from proc \a source.
 *
 * \remark the prism communication is not yet implemented
 *
 */
int PMMG_recv_mesh( MMG5_pMesh mesh,MMG5_pSol met,int source,int tag,MPI_Comm comm)
{
  MPI_Datatype   mpi_point,mpi_xpoint,mpi_tetra,mpi_xtetra;
  MPI_Status     *status;
  int            ier,ismet,np,ne,xp,xt;

  /* Get the mesh size */
  MPI_CHECK( MPI_Recv( &mesh->memMax,1,MPI_LONG_LONG,source,tag,comm,
                       status),return 0 );
  mesh->memCur = 0;

  /* Mesh size */
  MPI_CHECK( MPI_Recv( &np,1,MPI_INT,source,tag,comm,status ), return 0 );
  MPI_CHECK( MPI_Recv( &ne,1,MPI_INT,source,tag,comm,status ), return 0 );
  MPI_CHECK( MPI_Recv( &xp,1,MPI_INT,source,tag,comm,status ), return 0 );
  MPI_CHECK( MPI_Recv( &xt,1,MPI_INT,source,tag,comm,status ), return 0 );

  if ( !PMMG_grpSplit_setMeshSize( mesh, np, ne, 0, xp, xt ) ) return 0;

  /* Create mesh type */
  if ( !PMMG_create_MPI_Point( mesh->point,  &mpi_point ) ) goto end_5;
  if ( !PMMG_create_MPI_xPoint(mesh->xpoint, &mpi_xpoint) ) goto end_4;
  if ( !PMMG_create_MPI_Tetra( mesh->tetra,  &mpi_tetra ) ) goto end_3;
  if ( !PMMG_create_MPI_xTetra(mesh->xtetra, &mpi_xtetra) ) goto end_2;

  /* Recv mesh */
  MPI_CHECK( MPI_Recv(&mesh->point[1],mesh->np,mpi_point,source,tag,comm,
                      status),goto end_1;);
  MPI_CHECK( MPI_Recv(&mesh->xpoint[1],mesh->xp,mpi_xpoint,source,tag,comm,
                      status),goto end_1;);
  MPI_CHECK( MPI_Recv( &mesh->tetra[1] ,mesh->ne,mpi_tetra, source,tag,comm,
                       status),goto end_1;);
  MPI_CHECK( MPI_Recv( &mesh->xtetra[1],mesh->xt,mpi_xtetra,source,tag,comm,
                       status),goto end_1;);

  /* Met size */
  MPI_CHECK( MPI_Recv( &met->size,1,MPI_INT,source,tag,comm,
                       status ),goto end_1 );
  MPI_CHECK( MPI_Recv( &met->type, 1, MPI_INT, source, tag, comm,
                       status ),goto end_1 );
  if ( !MMG3D_Set_solSize(mesh,met,MMG5_Vertex,np,met->type) ) goto end_1;

  /* Recv metric */
  MPI_CHECK( MPI_Recv(&ismet,1,MPI_INT,source,tag,comm,status ),
             goto end_1 );

  if ( ismet )
    MPI_CHECK( MPI_Recv( met->m,met->size*met->np,MPI_DOUBLE,source,tag,comm,
                         status ),goto end_1 );

  ier = 1;

end_1:
  MPI_Type_free( &mpi_xtetra );
end_2:
  MPI_Type_free( &mpi_tetra );
end_3:
  MPI_Type_free( &mpi_xpoint );
end_4:
  MPI_Type_free( &mpi_point );
end_5:
  return ier;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param part pointer toward the metis array containing the partitions.
 * \param recv_comm number of comms to wait from another proc.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Process the parmesh groups and send and receive the mesh and metrics over
 * suitable processor.
 *
 */
static inline
int PMMG_mpiexchange_mesh(PMMG_pParMesh parmesh,idx_t *part,int *recv_comm) {
  MMG5_pMesh     mesh;
  MMG5_pSol      met;
  int            tag,ier,nprocs,ngrp,new_ngrp,k,i;

  ier    = 0;
  nprocs = parmesh->nprocs;
  ngrp   = parmesh->ngrp;

  /** Process the groups and send/receive the mesh and metrics that need
   * to be exchanged */
  for ( i=0; i<ngrp; ++i ) {
    if ( part[i] == parmesh->myrank ) continue;

    mesh = parmesh->listgrp[i].mesh;
    met  = parmesh->listgrp[i].met;

    tag = part[i]*nprocs + parmesh->myrank;
    if ( !PMMG_send_mesh(mesh,met,part[i],tag,parmesh->comm) ) goto end;
  }

  /** Pack the used groups on the mesh */
  new_ngrp = PMMG_pack_grps_norealloc( parmesh,&parmesh->listgrp,parmesh->ngrp);

  /** Count the final number of groups and reallocate the listgrp array */
  k = 0;
  for ( i=0; i<nprocs; ++i )
    k += recv_comm[i];

  PMMG_RECALLOC(parmesh,parmesh->listgrp,k+new_ngrp,parmesh->ngrp,PMMG_Grp,
                "listgrps",goto end);
  parmesh->ngrp = k+new_ngrp;

  for ( k=new_ngrp; k<parmesh->ngrp; ++k ) {
    MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&parmesh->listgrp[k].mesh,
                    MMG5_ARG_ppMet,&parmesh->listgrp[k].met,MMG5_ARG_end);
  }

  /** Receive the groups */
  for ( i=0; i<nprocs; ++i ) {
    /* If !recv_comm[i], we have nothing to wait from this proc */
    for ( k=0; k<recv_comm[i]; ++k ) {

      mesh = parmesh->listgrp[new_ngrp].mesh;
      met  = parmesh->listgrp[new_ngrp].met;

      tag = parmesh->myrank*nprocs + i;
      if ( !PMMG_recv_mesh(mesh,met,i,tag,parmesh->comm) ) goto end;
    }
  }

  /* Success */
  ier = 1;

end:
  return ier;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param idgrp index of the group that we process.
 * \param dest rank of destination.
 * \param tag message tag.
 * \param comm communicator.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Send the wanted part of the node communicator to the suitable processor.
 *
 */
static inline
int PMMG_send_nodeCommunicators(PMMG_pParMesh parmesh,int idgrp,int dest,int tag,
                                MPI_Comm comm) {
  PMMG_pint_comm int_comm;
  PMMG_pext_comm ext_comm;
  PMMG_ext_comm  ext_comm_2send;
  int            *intvalues,nitem;
  int            *node2int_node_comm_idx1,*node2int_node_comm_idx2;
  int            nitem_int_comm;
  int            next_comm_2send,idx,max_idx;
  int            ier,k,i,new_next_comm,old_next_comm,old_nitem;

return 1;
/*   ier    = 0; */

/*   /\** Process the list of external node communicators to create a new external */
/*    * communicator: if the node belong to this group (position in intvalue is */
/*    * filled), stores the position in the external communicator in a pack */
/*    * array and the intvalues value in another array. Both arrays will be */
/*    * sended by the external communicator. *\/ */
/*   int_comm  = parmesh->int_node_comm; */
/*   nitem     = int_comm->nitem; */
/*   intvalues = int_comm->intvalues; */

/*   /\** Step 1: intvalues array initialisation with the position of the interface */
/*    * points of the group in a packed intvalue array *\/ */
/*   for ( k=0; k<nitem; ++k ) */
/*     intvalues[k] = PMMG_UNSET; */

/*   node2int_node_comm_idx1 = parmesh->listgrp[idgrp].node2int_node_comm_index1; */
/*   node2int_node_comm_idx2 = parmesh->listgrp[idgrp].node2int_node_comm_index2; */
/*   nitem_int_comm          = parmesh->listgrp[idgrp].nitem_int_node_comm; */

/*   /\** Step 2: send the list of interface nodes of this group *\/ */
/*   MPI_CHECK( MPI_Send(&nitem_int_comm,1,MPI_INT,dest,tag, */
/*                       parmesh->comm),goto end ); */

/*   MPI_CHECK( MPI_Send(node2int_node_comm_idx1,nitem_int_comm,MPI_INT,dest,tag, */
/*                       parmesh->comm),goto end ); */

/*   for ( k=0; k<nitem_int_comm; ++k ) */
/*     intvalues[node2int_node_comm_idx2[k]] = k; */

/*   /\** Step 3: Process all the groups and find the nodes at the interface between */
/*    * the group that will move and other groups. Add this nodes to the external */
/*    * communicator (myrank-dest) *\/ */
/*   idx = 0; */
/*   new_next_comm = parmesh->next_node_comm; */
/*   old_next_comm = parmesh->next_node_comm; */
/*   for ( i=0; i<parmesh->ngrp; ++i ) { */
/*     if ( i==idgrp ) continue; */

/*     node2int_node_comm_idx1 = parmesh->listgrp[i].node2int_node_comm_index1; */
/*     node2int_node_comm_idx2 = parmesh->listgrp[i].node2int_node_comm_index2; */
/*     nitem_int_comm          = parmesh->listgrp[i].nitem_int_node_comm; */

/*     /\* Mark the nodes that already belong to the external comm to avoid */
/*      * doublons *\/ */
/*     for ( k=0; k<old_next_node_comm; ++k ) { */
/*       ext_comm = &parmesh->ext_node_comm[k]; */
/*       if ( ext_comm->color_out == dest ) { */
/*         for ( j=0; j<ext_comm->nitem; ++j ) */
/*           intvalues[ext_comm->int_comm_index[j]] *=-1; */
/*         break; */
/*       } */
/*     } */

/*     for ( j=0; j<nitem_int_comm; ++j ) */
/*       if ( intvalues[node2int_node_comm_idx2[j]] >= 0 ) ++idx; */

/*     /\* No common interfaces between the two groups: unmark the nodes and pass to */
/*      * the next group *\/ */
/*     if ( !idx ) { */
/*       for ( j=0; j<ext_comm->nitem; ++j ) */
/*         intvalues[ext_comm->int_comm_index[j]] *=-1; */
/*       continue; */
/*     } */

/*     /\* Common interfaces: creation or completion of the external communicator *\/ */
/*     if ( k==parmesh->next_node_comm ) { */
/*       /\* Creation of a new external communicator *\/ */
/*       ++new_next_comm; */
/*       if ( new_next_comm >= parmesh->next_node_comm ) { */
/*         PMMG_REALLOC(parmesh,parmesh->ext_node_comm, */
/*                      (1+PMMG_GAP)*parmesh->next_node_comm,parmesh->next_node_comm, */
/*                      PMMG_ext_comm,"list of external communicators",goto end); */
/*         parmesh->next_node_comm = (1+PMMG_GAP)*parmesh->next_node_comm; */
/*       } */
/*       ext_comm = parmesh->ext_node_comm[new_next_comm-1]; */
/*       ext_comm->color_in  = parmesh->myrank; */
/*       ext_comm->color_out = dest; */
/*     } */

/*     PMMG_REALLOC(parmesh,ext_comm->int_comm_index,ext_comm->nitem+idx, */
/*                  ext_comm->nitem,int,"int_comm_index",goto end); */
/*     old_nitem = ext_comm->nitem; */
/*     ext_comm->nitem += idx; */

/*     idx = 0; */
/*     for ( j=0; j<nitem_int_comm; ++j ) { */
/*       if ( intvalues[node2int_node_comm_idx2[j]] < 0 ) continue; */
/*       ext_comm->int_comm_index[old_nitem+idx] = node2int_node_comm_idx2[j]; */
/*     } */

/*     if ( k<parmesh->next_node_comm ) { */
/*       /\* Unmark the nodes that were in the external comm before *\/ */
/*       for ( j=0; j<old_nitem; ++j ) */
/*           intvalues[ext_comm->int_comm_index[j]] *=-1; */
/*     } */
/*   } */
/*   PMMG_REALLOC(parmesh,parmesh->ext_node_comm,new_next_comm, */
/*                parmesh->next_node_comm,PMMG_ext_comm, */
/*                "list of external communicators",goto end); */
/*   parmesh->next_node_comm = new_next_comm; */

/*   /\** Step 4: Count the number of external communicators to send and the max */
/*    * number of items in this communicators to allocate a new communicator */
/*    * sufficiently large to store the wanted data *\/ */
/*   for ( k=0; k<parmesh->next_node_comm; ++k ) { */
/*     ext_comm = &parmesh->ext_node_comm[k]; */

/*     /\* Count the number of items to send for this external comm *\/ */
/*     idx = 0; */
/*     for ( i=0; i<ext_comm->nitem; ++i ) ++idx; */

/*     if ( !idx )  continue; */

/*     max_idx = MG_MAX(idx,max_idx); */
/*     ++next_comm_2send; */
/*   } */

/*   PMMG_MALLOC(parmesh,ext_comm_2send.int_comm_index,max_idx,int, */
/*               "int_comm_index to send",goto end); */
/*   PMMG_MALLOC(parmesh,ext_comm_2send.int_comm_index,max_idx,int,"i2send", */
/*               goto end); */

/*   MPI_CHECK( MPI_Send( &next_comm_2send, 1,MPI_INT,dest,tag,parmesh->comm ), */
/*              goto end ); */

/*   /\** Step 5: Fill the new external communicator with the data that we want to */
/*    * communicate: the list of the position of the nodes of this group (in the */
/*    * external communicator) in the packed internal communicator as well as their */
/*    * positions positions int the external communicator *\/ */
/*   next_comm_2send = 0; */
/*   for ( k=0; k<parmesh->next_node_comm; ++k ) { */
/*     ext_comm = &parmesh->ext_node_comm[k]; */

/*     /\* Count the number of items to send for this external comm *\/ */
/*     idx = 0; */
/*     for ( i=0; i<ext_comm->nitem; ++i ) ++idx; */

/*     if ( !idx ) continue; */

/*     /\* Fill the ext_comm_2send array *\/ */
/*     ext_comm_2send.nitem     = idx; */
/*     ext_comm_2send.color_in  = ext_comm->color_in; */
/*     ext_comm_2send.color_out = ext_comm->color_out; */

/*     idx = 0; */
/*     for ( i=0; i<ext_comm->nitem; ++i ) { */
/*       if ( PMMG_UNSET != intvalues[ext_comm->int_comm_index[i]] ) { */
/*         ext_comm_2send.int_comm_index[idx] = i; */
/*         ext_comm_2send.itosend[idx] = */
/*           intvalues[ext_comm->int_comm_index[i]]; */
/*         ++idx; */
/*       } */
/*     } */
/* #warning add the creation of an MPI_type to send all this data efficiently */
/*     MPI_CHECK( MPI_Send(&ext_comm_2send.nitem,1,MPI_INT,dest,tag,parmesh->comm), */
/*                goto end ); */
/*     MPI_CHECK( MPI_Send(&ext_comm_2send.color_out,1,MPI_INT,dest,tag,parmesh->comm), */
/*                goto end ); */
/*     MPI_CHECK( MPI_Send(ext_comm_2send.int_comm_index,ext_comm_2send.nitem, */
/*                         MPI_INT,dest,tag,parmesh->comm ),goto end ); */
/*     MPI_CHECK( MPI_Send(ext_comm_2send.itosend,ext_comm_2send.nitem,MPI_INT, */
/*                         dest,tag,parmesh->comm ),goto end ); */
/*   } */
/*   PMMG_DEL_MEM(parmesh,ext_comm_2send.int_comm_index,max_idx,int, */
/*                "int_comm_index to send"); */
/*   PMMG_DEL_MEM(parmesh,ext_comm_2send.int_comm_index,max_idx,int,"i2send"); */

/*   /\* Success *\/ */
/*   ier = 1; */

/* end: */
/*   if ( ext_comm_2send.int_comm_index ) */
/*     PMMG_DEL_MEM(parmesh,ext_comm_2send.int_comm_index, */
/*                  ext_comm_2send.nitem,int,"int_comm_index to send"); */

/*   if ( ext_comm_2send.itosend ) */
/*     PMMG_DEL_MEM(parmesh,ext_comm_2send.itosend, */
/*                  ext_comm_2send.nitem,int,"i2send array to send"); */

/*   return ier; */
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param source rank of source.
 * \param tag message tag
 * \param comm communicator
 *
 * \return 0 if fail, 1 otherwise
 *
 * Send the wanted part of the node communicator to the suitable processor.
 *
 */
static inline
int PMMG_recv_nodeCommunicators(PMMG_pParMesh parmesh,int source,int tag,
                                MPI_Comm comm) {
  PMMG_pint_comm int_comm;
  PMMG_pext_comm ext_comm;
  PMMG_ext_comm  ext_comm_2send;
  int            *intvalues,nitem;
  int            *node2int_node_comm_idx1,*node2int_node_comm_idx2;
  int            nitem_int_comm;
  int            next_comm_2recv,idx,max_idx;
  int            ier,k,i,new_next_comm,old_next_comm;

  return 1;

/*   ier    = 0; */

/*   /\** Step 1: recv the list of interface nodes of this group *\/ */
/*   MPI_CHECK( MPI_Recv(&nitem_int_comm,1,MPI_INT,source,tag, */
/*                       parmesh->comm,MPI_STATUS_IGNORE),goto end ); */
/*   MPI_CHECK( MPI_Recv(node2int_node_comm_idx1,nitem_int_comm,MPI_INT,source,tag, */
/*                       parmesh->comm,MPI_STATUS_IGNORE),goto end ); */

/*   MPI_CHECK( MPI_Recv(&next_comm_2recv,1,MPI_INT,source,tag,parmesh->comm, */
/*                       MPI_STATUS_IGNORE), goto end ); */

/*   new_next_comm = parmesh->next_node_comm; */
/*   old_next_comm = parmesh->next_node_comm; */

/*   for ( i=0; i<next_comm_2recv; ++i ) { */
/*     /\** Step 2: Receive the external communicator *\/ */
/*     MPI_CHECK( MPI_Recv(&ext_comm_2send.nitem,1,MPI_INT,source,tag,parmesh->comm, */
/*                         MPI_STATUS_IGNORE),goto end ); */
/*     MPI_CHECK( MPI_Recv(&ext_comm_2send.color_out,1,MPI_INT,source,tag, */
/*                         parmesh->comm,MPI_STATUS_IGNORE),goto end ); */
/*     MPI_CHECK( MPI_Recv(ext_comm_2send.int_comm_index,ext_comm_2send.nitem, */
/*                         MPI_INT,source,tag,parmesh->comm, */
/*                         MPI_STATUS_IGNORE),goto end ); */
/*     MPI_CHECK( MPI_Recv(ext_comm_2send.itosend,ext_comm_2send.nitem, */
/*                         MPI_INT,source,tag,parmesh->comm, */
/*                         MPI_STATUS_IGNORE),goto end ); */

/*     /\** Step 2: Merge the external communicators *\/ */
/*     if ( ext_comm_2send->color_out==parmesh->myrank ) { */
/*       /\* External communicator for points between me and the proc that send me */
/*        * the group: remove the points that are no longer in the external */
/*        * communicator *\/ */
/*       for ( k=0; k<old_next_node_comm; ++k ) { */
/*         ext_comm = parmesh->ext_node_comm[k]; */
/*         if ( ext_comm->color_out == source ) break; */
/*       } */
/*       if ( k==old_next_node_comm ) { */
/*         printf("  ## Error: %s: unable to find the external communicator %d-%d.\n", */
/*           __func__,parmesh->myrank,source); */
/*         goto end; */
/*       } */

/*       for ( k=0; k<ext_comm_2send.nitem; ++k ) { */
/*         idx = ext_comm_2send.i2send[k]; */
        
/*       } */

/*     } */
/*     else { */
/*       /\* External communicator for points between me and an arbitrary proc *\/ */


/*     } */


/*     /\** Step 2: Search if this communicator already exists in the list of */
/*      * external communicators or if we need to add it *\/ */
/*     for ( k=0; k<old_next_node_comm; ++k ) { */
/*       ext_comm = parmesh->ext_node_comm[k]; */
/*       if ( ext_comm->color_out == source ) break; */
/*     } */

/*     if ( k==old_next_node_comm ) { */
/*       ++new_next_comm; */
/*       if ( new_next_comm >= parmesh->next_node_comm ) { */
/*         PMMG_REALLOC(parmesh,parmesh->ext_node_comm, */
/*                      (1+PMMG_GAP)*parmesh->next_node_comm,parmesh->next_node_comm, */
/*                      PMMG_ext_comm,"list of external communicators",goto end); */
/*         parmesh->next_node_comm = (1+PMMG_GAP)*parmesh->next_node_comm; */
/*       } */
/*       ext_comm = parmesh->ext_node_comm[new_next_comm-1]; */
/*     } */

/*     /\** Step 3: *\/ */

/*   } */


/*   /\** Process the list of external node communicators to create a new external */
/*    * communicator: if the node belong to this group (position in intvalue is */
/*    * filled), stores the position in the external communicator in a pack */
/*    * array and the intvalues value in another array. Both arrays will be */
/*    * sended by the external communicator. *\/ */
/*   int_comm  = parmesh->int_node_comm; */
/*   nitem     = int_comm->nitem; */
/*   intvalues = int_comm->intvalues; */

/*   /\** Step 1: intvalues array initialisation with the index of the interface */
/*    * points of the group *\/ */
/*   for ( k=0; k<nitem; ++k ) */
/*     intvalues[k] = PMMG_UNSET; */

/*   node2int_node_comm_idx1 = parmesh->listgrp[i].node2int_node_comm_index1; */
/*   node2int_node_comm_idx2 = parmesh->listgrp[i].node2int_node_comm_index2; */
/*   nitem_int_comm          = parmesh->listgrp[i].nitem_int_node_comm; */

/*   for ( k=0; k<nitem_int_comm; ++k ) */
/*     intvalues[node2int_node_comm_idx2[k]] = node2int_node_comm_idx1[k]; */

/*   /\** Step 2: Count the number of external communicators to send and the max */
/*    * number of items in this communicators to allocate a new communicator */
/*    * sufficiently large to store the wanted data *\/ */
/*   max_idx = 0; */
/*   next_comm_2send = 0; */
/*   for ( k=0; k<parmesh->next_node_comm; ++k ) { */
/*     ext_comm = &parmesh->ext_node_comm[k]; */

/*     /\* Count the number of items to send for this external comm *\/ */
/*     idx = 0; */
/*     for ( i=0; i<ext_comm->nitem; ++i ) ++idx; */

/*     if ( !idx )  continue; */

/*     max_idx = MG_MAX(idx,max_idx); */
/*     ++next_comm_2send; */
/*   } */

/*   PMMG_MALLOC(parmesh,ext_comm_2send.int_comm_index,max_idx,int, */
/*               "int_comm_index to send",goto end); */
/*   PMMG_MALLOC(parmesh,ext_comm_2send.int_comm_index,max_idx,int,"i2send", */
/*               goto end); */


/*   /\** Step 3: Fill the new external communicator with the data that we want to */
/*    * communicate: the list of nodes of this group that are in the external */
/*    * communicator and their positions int this communicator *\/ */
/*   next_comm_2send = 0; */
/*   for ( k=0; k<parmesh->next_node_comm; ++k ) { */
/*     ext_comm = &parmesh->ext_node_comm[k]; */

/*     /\* Count the number of items to send for this external comm *\/ */
/*     idx = 0; */
/*     for ( i=0; i<ext_comm->nitem; ++i ) ++idx; */

/*     if ( !idx ) continue; */

/*     /\* Fill the ext_comm_2send array *\/ */
/*     ext_comm_2send.nitem     = idx; */
/*     ext_comm_2send.color_in  = ext_comm->color_in; */
/*     ext_comm_2send.color_out = ext_comm->color_out; */

/*     idx = 0; */
/*     for ( i=0; i<ext_comm->nitem; ++i ) { */
/*       if ( PMMG_UNSET != intvalues[ext_comm->int_comm_index[i]] ) { */
/*         ext_comm_2send.int_comm_index[idx] = i; */
/*         ext_comm_2send.itosend[idx] = */
/*           intvalues[ext_comm->int_comm_index[i]]; */
/*         ++idx; */
/*       } */
/*     } */
/*     MPI_CHECK( MPI_Recv(ext_comm_2send.int_comm_index,idx,MPI_INT, */
/*                         part[i],tag,parmesh->comm ),goto end ); */
/*     MPI_CHECK( MPI_Recv(ext_comm_2send.itosend,idx,MPI_INT, */
/*                         part[i],tag,parmesh->comm ),goto end ); */
/*   } */
/*   PMMG_DEL_MEM(parmesh,ext_comm_2send.int_comm_index,max_idx,int, */
/*                "int_comm_index to send"); */
/*   PMMG_DEL_MEM(parmesh,ext_comm_2send.int_comm_index,max_idx,int,"i2send"); */

/*   /\* Success *\/ */
/*   ier = 1; */

/* end: */
/*   if ( ext_comm_2send.int_comm_index ) */
/*     PMMG_DEL_MEM(parmesh,ext_comm_2send.int_comm_index, */
/*                  ext_comm_2send.nitem,int,"int_comm_index to send"); */

/*   if ( ext_comm_2send.itosend ) */
/*     PMMG_DEL_MEM(parmesh,ext_comm_2send.itosend, */
/*                  ext_comm_2send.nitem,int,"i2send array to send"); */

/*   return ier; */
}


/**
 * \param parmesh pointer toward the mesh structure.
 * \param part pointer toward the metis array containing the partitions.
 * \param recv_comm number of comms to wait from another proc.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Process the parmesh groups and send and receive the wanted part of the node
 * communicators over the suitable processor.
 *
 */
static inline
int PMMG_mpiexchange_nodeCommunicators(PMMG_pParMesh parmesh,idx_t *part,
                                       int* recv_comm) {
  PMMG_pint_comm int_comm;
  int            nitem,*intvalues;
  int            tag,ier,nprocs,ngrp,k,i;

  ier    = 0;
  nprocs = parmesh->nprocs;
  ngrp   = parmesh->ngrp;

  /** Step 1: Process the groups and send the node communicators that
   * need to be sended */
  int_comm = parmesh->int_node_comm;
  nitem    = int_comm->nitem;

  PMMG_MALLOC(parmesh,int_comm->intvalues,nitem,int,"intvalues",goto end);

  for ( i=0; i<ngrp; ++i ) {
    if ( part[i] == parmesh->myrank ) continue;

    /* Send the wanted data of the external node communicator of this group */
    tag = part[i]*nprocs + parmesh->myrank;
    if ( !PMMG_send_nodeCommunicators(parmesh,i,part[i],tag,parmesh->comm) )
      goto end;
  }
  PMMG_DEL_MEM(parmesh,intvalues,nitem,int,"intvalues");

  /** Step 2: Process the groups and received the node communicators that
   * need to be received */
  for ( i=0; i<nprocs; ++i ) {
    /* If !recv_comm[i], we have nothing to wait from this proc */
    for ( k=0; k<recv_comm[i]; ++k ) {

      tag = parmesh->myrank*nprocs + i;
      if ( !PMMG_recv_nodeCommunicators(parmesh,i,tag,parmesh->comm) )
        goto end;
    }
  }

  /* Success */
  ier = 1;

end:
  if ( intvalues )
    PMMG_DEL_MEM(parmesh,intvalues,nitem,int,"intvalues");

  return ier;

}


/**
 * \param parmesh pointer toward the mesh structure.
 * \param part pointer toward the metis array containing the partitions.
 * \param recv_comm number of comms to wait from another proc.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Process the parmesh groups and send and receive the wanted part of the face
 * communicators over the suitable processor.
 *
 */
static inline
int PMMG_mpiexchange_faceCommunicators(PMMG_pParMesh parmesh,idx_t *part,
                                       int* recv_comm) {

  return 1;
}


/**
 * \param parmesh pointer toward the mesh structure.
 * \param part pointer toward the metis array containing the partitions.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Send the suitable group to other procs and recieve their groups.
 * Deallocate the \a part array.
 *
 */
static inline
int PMMG_mpiexchange_grps(PMMG_pParMesh parmesh,idx_t *part) {
  int *send_comm,*recv_comm,nprocs;
  int ier,ier_glob,k;

  /** Merge all the groups that must be send to a given proc into 1 group */
  ier = PMMG_merge_grps2send(parmesh,&part);
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);

  if ( ier_glob < 0 ) {
    fprintf(stderr,"\n  ## Unable to compute the new group partition.\n");
    ier = 0;
    goto end;
  }
  else if ( !ier_glob ) {
    fprintf(stderr,"\n  ## Unable to compute the new group partition."
            " Try to send it nevertheless\n");
  }

  /** Send each group to the suitable processor */
  /* Send a logical to say the number of comms that we will receive from this
   * proc */
  nprocs = parmesh->nprocs;
  PMMG_CALLOC(parmesh,send_comm,nprocs+1,int,"sended comm with other procs",
              goto end);
  PMMG_CALLOC(parmesh,recv_comm,nprocs+1,int,"received comm with other procs",
              goto end);

  for ( k=0; k<parmesh->ngrp; ++k )
    if ( part[k] != parmesh->myrank ) ++send_comm[part[k]];

  MPI_CHECK( MPI_Alltoall(send_comm,1,MPI_INT,recv_comm,1,MPI_INT,parmesh->comm),
             goto end);

  PMMG_DEL_MEM(parmesh,send_comm,nprocs+1,int,"send_comm");

  /* Exchange the meshes and metrics */
  if ( !PMMG_mpiexchange_mesh(parmesh,part,recv_comm) ) {
    fprintf(stderr,"\n  ## Unable to move meshes and metrics over"
            " the suitable proc.\n");
    ier = 0;
    goto end;
  }

  /* Exchange the needed part of the node communicators */
  /* if ( !PMMG_mpiexchange_nodeCommunicators(parmesh,part,recv_comm) ) { */
  /*   fprintf(stderr,"\n  ## Unable to move node communicators over" */
  /*           " the suitable proc.\n"); */
  /*   ier = 0; */
  /*   goto end; */
  /* } */

  /* Exchange the needed part of the face communicators */
  /* if ( !PMMG_mpiexchange_faceCommunicators(parmesh,part,recv_comm) ) { */
  /*   fprintf(stderr,"\n  ## Unable to move face communicators over" */
  /*           " the suitable proc.\n"); */
  /*   ier = 0; */
  /*   goto end; */
  /* } */

  PMMG_DEL_MEM(parmesh,part,parmesh->ngrp,idx_t,"deallocate parmetis partition");

  /** Success */
  ier = 1;

end:
  if ( part )
    PMMG_DEL_MEM(parmesh,part,parmesh->ngrp,idx_t,"deallocate parmetis partition");
  if ( send_comm )
    PMMG_DEL_MEM(parmesh,send_comm,nprocs+1,int,"send_comm");
  if ( recv_comm )
    PMMG_DEL_MEM(parmesh,recv_comm,nprocs+1,int,"recv_comm");

  return ier;
}

/**
 * \param parmesh pointer toward the mesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Compute the group distribution that load balance the groups over the
 * processors and send and recieve the groups from the other processors.
 *
 */
int PMMG_distribute_grps( PMMG_pParMesh parmesh ) {
  idx_t *part;

  /** Get the new partition of groups (1 group = 1 metis node) */
  PMMG_CALLOC(parmesh,part,parmesh->ngrp,idx_t,"allocate parmetis buffer",
              return 0);

  if ( !PMMG_part_parmeshGrps2parmetis(parmesh,part,parmesh->nprocs) ) {
    fprintf(stderr,"\n  ## Unable to compute the new group partition.\n");
    PMMG_DEL_MEM(parmesh,part,parmesh->ngrp,idx_t,"deallocate parmetis partition");
    return 0;
  }

  /** Send the suitable groups to other procs and recieve their groups */
  if ( !PMMG_mpiexchange_grps(parmesh,part) ) {
    fprintf(stderr,"\n  ## Unable to communicate groups through processors.\n");
    return 0;
  }

  return 1;
}

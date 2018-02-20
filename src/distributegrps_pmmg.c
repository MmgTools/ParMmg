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
    PMMG_REALLOC(parmesh,ext_node_comm->int_comm_index,nitem_ext,
                 ext_node_comm->nitem,int,"int_comm_index",return 0);
    ext_node_comm->nitem = nitem_ext;
  }

  /** Step 5: unallocate intvalues array and set the nitem field of the internal
   * communicator to the suitable value */
  PMMG_DEL_MEM( parmesh,int_node_comm->intvalues,int_node_comm->nitem,int,
                "node communicator");
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
static inline
int PMMG_pack_faceCommunicators(PMMG_pParMesh parmesh) {
  PMMG_pGrp      grp;
  PMMG_pint_comm int_face_comm;
  PMMG_pext_comm ext_face_comm;
  int            *intvalues;
  int            *face2int_face_comm_index2;
  int            k,nitem_int,nitem_ext,idx_int,idx_ext,i;

  intvalues = parmesh->int_face_comm->intvalues;

  /** Step 1: initialization of the communicator */
  int_face_comm = parmesh->int_face_comm;
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
    PMMG_REALLOC(parmesh,ext_face_comm->int_comm_index,nitem_ext,
                 ext_face_comm->nitem,int,"int_comm_index",return 0);
    ext_face_comm->nitem = nitem_ext;
  }

  /** Step 5: unallocate intvalues array and set the nitem field of the internal
   * communicator to the suitable value */
  PMMG_DEL_MEM( parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,
                "face communicator");
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

    /* Update the communicators: WARNING TO IMPROVE: INEFFICIENT */
    if ( !PMMG_mergeGrpJinI_communicators(parmesh,grpI,grpJ,grps,k) ) goto low_fail;

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
 * \param parmesh pointer toward a parmesh structure
 * \param grp_id id of the current group
 * \param source id or the proc from which we receive
 * \param tag communication tag
 * \param recv_meshsize array to store the mesh size data
 * \param req mpi request array
 *
 * \return 1 if success, 0 if fail.
 *
 * Receive the needed info to allocate a mesh and a metric.
 *
 */
int PMMG_mpirecv_meshSize( PMMG_pParMesh parmesh,int grp_id,int source,int tag,
                           int *recv_meshsize,MPI_Request *req )
{
  MMG5_pMesh     mesh;

  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&parmesh->listgrp[grp_id].mesh,
                  MMG5_ARG_ppMet ,&parmesh->listgrp[grp_id].met,
                  MMG5_ARG_end);

  mesh = parmesh->listgrp[grp_id].mesh;

  /* Get the mesh size */
  MPI_CHECK( MPI_Irecv( &mesh->memMax,1,MPI_LONG_LONG,source,tag,parmesh->comm,
                        &req[0] ),return 0 );

  /* Mesh and met size:
   * meshsize[]=[mesh->np,mesh->ne,mesh->xp,mesh->xt,
   *             met->size,met->type,ismet] */
  MPI_CHECK( MPI_Irecv( recv_meshsize,7,MPI_INT,source,tag,parmesh->comm,&req[1] ),
             return 0 );

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param grp_id id of the current group
 * \param dest id of the proc to which we send data
 * \param send_meshsize array to store the mesh size data
 * \param req mpi request array
 *
 * \return 1 if success, 0 if fail.
 *
 * Fit the mesh and metric to its minimal size and send the needed info to
 * allocate its.
 *
 */
int PMMG_mpisend_meshSize( PMMG_pParMesh parmesh,int grp_id,int dest,int tag,
                           int *send_meshsize,MPI_Request *req )
{
  MMG5_pMesh     mesh;
  MMG5_pSol      met;
  MPI_Comm       comm;

  mesh = parmesh->listgrp[grp_id].mesh;
  met  = parmesh->listgrp[grp_id].met;
  comm = parmesh->comm;

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

  /* Send the mesh/met size:
   * meshsize[]=[mesh->np,mesh->ne,mesh->xp,mesh->xt,
   *             met->size,met->type,ismet] */
  mesh->memMax = mesh->memCur;
  MPI_CHECK( MPI_Isend( &mesh->memMax,1,MPI_LONG_LONG,dest,tag,comm,&req[0] ),
             return 0 );

  send_meshsize[0] = mesh->np;
  send_meshsize[1] = mesh->ne;
  send_meshsize[2] = mesh->xp;
  send_meshsize[3] = mesh->xt;
  send_meshsize[4] = met->size;
  send_meshsize[5] = met->type;
  send_meshsize[6] = met->m ? 1:0;

  MPI_CHECK( MPI_Isend( send_meshsize,7,MPI_INT,dest,tag,comm,&req[1] ),
             return 0 );

  return 1;
}


/**
 * \param parmesh pointer toward a parmesh structure
 * \param first_grp_id id of the first available group
 * \param part pointer toward the metis array containing the partitions
 * \param send_grps array of the groups that we send to each proc
 * \param recv_grps array of the groups waited from each proc
 *
 * \return the number of groups for which we have successfully send and receive
 * the array sizes
 *
 * Resizes the meshes and send and receive the needed sizes to communicate
 * groups.
 *
 */
/* int PMMG_mpiexchange_grpSize( PMMG_pParMesh parmesh,int first_grp_id, */
/*                               int *part,int *send_grps,int *recv_grps) */
/* { */
/*   int            ier,nprocs,k,i,pos,grp_idx,nrecv_success,nsend_success; */
/*   int            recv_meshsize[7],send_meshsize[7],nsend,nrecv,tag; */
/*   MPI_Request    *send_req ,*recv_req; */
/*   MPI_Status     *send_stat,*recv_stat; */
/*   MPI_Comm       comm; */

/*   ier = 0; */

/*   nprocs = parmesh->nprocs; */
/*   comm   = parmesh->comm; */

/*   /\** Step 1: Mesh initialization and reception of all the size (mesh+metric) *\/ */
/*   nrecv = 0; */
/*   recv_req  = NULL; */
/*   recv_stat = NULL; */

/*   nsend = 0; */
/*   send_req  = NULL; */
/*   send_stat = NULL; */

/*   for ( i=0; i<nprocs; ++i ) nrecv += recv_grps[i]; */

/*   PMMG_MALLOC(parmesh,recv_req, 2*nrecv,MPI_Request,"recv_req" ,goto end); */
/*   PMMG_MALLOC(parmesh,recv_stat,2*nrecv,MPI_Status ,"recv_stat",goto end); */

/*   nrecv_success = 0; */
/*   for ( i=0; i<nprocs; ++i ) { */
/*     nrecv = recv_grps[i]; */

/*     for ( k=0; k<recv_grps[i]; ++k ) { */
/*       grp_idx = first_grp_id+nrecv_success; */
/*       pos     = 2*nrecv_success; */

/*       tag = MPI_DISTRIBUTEGRPS_MESHSIZE_TAG; */
/*       if ( !PMMG_mpirecv_meshSize( parmesh,grp_idx,i,tag, */
/*                                    recv_meshsize,&recv_req[pos]) ) goto end; */

/*     } */
/*     ++nrecv_success; */
/*   } */

/*   for ( k=0; k<parmesh->ngrp; ++k ) */
/*     if ( part[k] != parmesh->myrank ) ++nsend; */

/*   PMMG_MALLOC(parmesh,send_req ,2*nsend,MPI_Request,"send_req" ,goto end); */
/*   PMMG_MALLOC(parmesh,send_stat,2*nsend,MPI_Status ,"send_stat",goto end); */


/*   for ( k=0; k<parmesh->ngrp; ++k ) { */
/*     if ( part[k] == parmesh->myrank ) continue; */

/*     tag = MPI_DISTRIBUTEGRPS_MESHSIZE_TAG; */
/*     pos = 2*nsend_success; */

/*     if ( !PMMG_mpisend_meshSize( parmesh,k,part[k],tag,send_meshsize, */
/*                                  &send_req[pos]) )  goto end; */

/*     ++nsend_success; */
/*   } */

/*   MPI_Waitall(2*nrecv,recv_req,recv_stat); */
/*   MPI_Waitall(2*nsend,send_req,send_stat); */

/*   /\* Success *\/ */
/*   ier = 1; */

/*   end: */
/*   if ( recv_req ) */
/*     PMMG_DEL_MEM(parmesh,recv_req, 2*nrecv,MPI_Request,"recv_req"); */
/*   if ( recv_stat ) */
/*     PMMG_DEL_MEM(parmesh,recv_stat,2*nrecv,MPI_Status, "recv_stat"); */
/*   if ( send_req ) */
/*     PMMG_DEL_MEM(parmesh,send_req, 2*nsend,MPI_Request,"send_req"); */
/*   if ( send_stat ) */
/*     PMMG_DEL_MEM(parmesh,send_stat,2*nsend,MPI_Status, "send_stat"); */

/*   return ier; */
/* } */

/**
 * \param parmesh pointer toward a parmesh structure.
 * \param idgrp index of the group to treat.
 * \param dest rank of destination.
 * \param tag message tag
 * \param comm communicator
 *
 * \return 0 if fail, 1 if success.
 *
 * Send mesh and metric of the group \a grp to proc \a dest.
 *
 * \remark the prisms are not treated.
 *
 */
int PMMG_send_mesh( PMMG_pParMesh parmesh,int idgrp,int dest,int tag,MPI_Comm comm)
{
  MMG5_pMesh     mesh;
  MMG5_pSol      met;
  MPI_Datatype   mpi_point,mpi_xpoint,mpi_tetra,mpi_xtetra;
  int            ier,ismet,k;
  MPI_Request    req[12];
  MPI_Status     status;

  mesh = parmesh->listgrp[idgrp].mesh;
  met  = parmesh->listgrp[idgrp].met;

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
  MPI_CHECK( MPI_Isend( &mesh->memMax,1,MPI_LONG_LONG,dest,tag,comm,&req[0] ),
             return 0 );
  MPI_CHECK( MPI_Isend( &mesh->np,1,MPI_INT,dest,tag,comm,&req[1] ),return 0 );
  MPI_CHECK( MPI_Isend( &mesh->ne,1,MPI_INT,dest,tag,comm,&req[2] ),return 0 );
  MPI_CHECK( MPI_Isend( &mesh->xp,1,MPI_INT,dest,tag,comm,&req[3] ),return 0 );
  MPI_CHECK( MPI_Isend( &mesh->xt,1,MPI_INT,dest,tag,comm,&req[4] ),return 0 );
  MPI_Wait(&req[0], &status);
  MPI_Wait(&req[1], &status);
  MPI_Wait(&req[2], &status);
  MPI_Wait(&req[3], &status);
  MPI_Wait(&req[4], &status);

  ier = 0;
  /* Create mesh type */
  if ( !PMMG_create_MPI_Point( mesh->point,  &mpi_point ) ) goto end_5;
  if ( !PMMG_create_MPI_xPoint(mesh->xpoint, &mpi_xpoint) ) goto end_4;
  if ( !PMMG_create_MPI_Tetra( mesh->tetra,  &mpi_tetra ) ) goto end_3;
  if ( !PMMG_create_MPI_xTetra(mesh->xtetra, &mpi_xtetra) ) goto end_2;

  /* Send mesh */
  MPI_CHECK( MPI_Isend( &mesh->point[1] ,mesh->np,mpi_point, dest,tag,comm,
                        &req[5]),goto end_1;);
  MPI_Wait(&req[5], &status);

  MPI_CHECK( MPI_Isend( &mesh->xpoint[1],mesh->xp,mpi_xpoint,dest,tag,comm,
                        &req[6]),goto end_1;);
  MPI_Wait(&req[6], &status);

  MPI_CHECK( MPI_Isend( &mesh->tetra[1] ,mesh->ne,mpi_tetra, dest,tag,comm,
                        &req[7]),goto end_1;);
  MPI_Wait(&req[7], &status);

  MPI_CHECK( MPI_Isend( &mesh->xtetra[1],mesh->xt,mpi_xtetra,dest,tag,comm,
                        &req[8]),goto end_1;);
  MPI_Wait(&req[8], &status);


  /* Met size */
  assert ( met->np == mesh->np );
  MPI_CHECK( MPI_Isend( &met->size, 1, MPI_INT, dest, tag, comm, &req[9] ),
             goto end_1 );
  MPI_CHECK( MPI_Isend( &met->type, 1, MPI_INT, dest, tag, comm, &req[10]),
             goto end_1 );

  /* Send metric */
  ismet = met->m ? 1:0;
  MPI_CHECK( MPI_Isend( &ismet, 1, MPI_INT, dest, tag, comm, &req[11] ),
             goto end_1 );

  if ( ismet )
    MPI_CHECK( MPI_Isend( met->m,met->size*met->np,MPI_DOUBLE,dest,tag,comm,
                          &req[12] ),goto end_1 );

  for ( k=0; k<13; ++k )
    MPI_Wait(&req[k], &status);

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
 * \param parmesh pointer toward a parmesh structure
 * \param first_grp_id id of the first available group
 * \param source rank of source.
 * \param tag message tag
 * \param comm communicator
 * \param recv_grps number of groups coming from the proc source.
 *
 * \return the number of groups successfully received (so recv_grps if success)
 *
 * Receive mesh and metric from proc \a source.
 *
 * \remark the prism communication is not yet implemented
 *
 */
int PMMG_recv_mesh( PMMG_pParMesh parmesh,int first_grp_id,int source,int tag,
                    MPI_Comm comm,int recv_grps)
{
  MMG5_pMesh     mesh;
  MMG5_pSol      met;
  MPI_Datatype   mpi_point,mpi_xpoint,mpi_tetra,mpi_xtetra;
  int            k,l,ismet,np,ne,xp,xt;
  MPI_Request    req[12];
  MPI_Status     status;

  mpi_xtetra = NULL;
  mpi_tetra  = NULL;
  mpi_xpoint = NULL;
  mpi_point  = NULL;

  for ( k=0; k<recv_grps; ++k ) {

    MMG3D_Init_mesh(MMG5_ARG_start,
                    MMG5_ARG_ppMesh,&parmesh->listgrp[first_grp_id+k].mesh,
                    MMG5_ARG_ppMet,&parmesh->listgrp[first_grp_id+k].met,
                    MMG5_ARG_end);

    mesh = parmesh->listgrp[first_grp_id+k].mesh;
    met  = parmesh->listgrp[first_grp_id+k].met;

    /* Get the mesh size */
    MPI_CHECK( MPI_Irecv( &mesh->memMax,1,MPI_LONG_LONG,source,tag,comm,
                          &req[0]),return k );
    mesh->memCur = 0;

    /* Mesh size */
    MPI_CHECK( MPI_Irecv( &np,1,MPI_INT,source,tag,comm,&req[1] ), return k );
    MPI_CHECK( MPI_Irecv( &ne,1,MPI_INT,source,tag,comm,&req[2] ), return k );
    MPI_CHECK( MPI_Irecv( &xp,1,MPI_INT,source,tag,comm,&req[3] ), return k );
    MPI_CHECK( MPI_Irecv( &xt,1,MPI_INT,source,tag,comm,&req[4] ), return k );

    if ( !PMMG_grpSplit_setMeshSize( mesh, np, ne, 0, xp, xt ) ) return k;

    /* Create mesh type */
    // To test : creation of the mpi types outside this loop
    if ( !PMMG_create_MPI_Point( mesh->point,  &mpi_point ) ) goto end_5;
    if ( !PMMG_create_MPI_xPoint(mesh->xpoint, &mpi_xpoint) ) goto end_4;
    if ( !PMMG_create_MPI_Tetra( mesh->tetra,  &mpi_tetra ) ) goto end_3;
    if ( !PMMG_create_MPI_xTetra(mesh->xtetra, &mpi_xtetra) ) goto end_2;

    /* Recv mesh */
    MPI_CHECK( MPI_Irecv(&mesh->point[1],mesh->np,mpi_point,source,tag,comm,
                 &req[5]),goto end_1;);
    MPI_CHECK( MPI_Irecv(&mesh->xpoint[1],mesh->xp,mpi_xpoint,source,tag,comm,
                         &req[6]),goto end_1;);
    MPI_CHECK( MPI_Irecv( &mesh->tetra[1] ,mesh->ne,mpi_tetra, source,tag,comm,
                          &req[7]),goto end_1;);
    MPI_CHECK( MPI_Irecv( &mesh->xtetra[1],mesh->xt,mpi_xtetra,source,tag,comm,
                          &req[8]),goto end_1;);

    /* Met size */
    MPI_CHECK( MPI_Irecv( &met->size,1,MPI_INT,source,tag,comm,
                          &req[9] ),goto end_1 );
    MPI_CHECK( MPI_Irecv( &met->type, 1, MPI_INT, source, tag, comm,
                          &req[10] ),goto end_1 );
    if ( !MMG3D_Set_solSize(mesh,met,MMG5_Vertex,np,met->type) ) goto end_1;

    /* Recv metric */
    MPI_CHECK( MPI_Irecv(&ismet,1,MPI_INT,source,tag,comm,&req[11] ),
               goto end_1 );

    if ( ismet )
      MPI_CHECK( MPI_Irecv( met->m,met->size*met->np,MPI_DOUBLE,source,tag,comm,
                            &req[12] ),goto end_1 );

    for ( l=0; l<13; ++l )
      MPI_Wait(&req[l], &status);

    MPI_Type_free( &mpi_xtetra );
    MPI_Type_free( &mpi_tetra );
    MPI_Type_free( &mpi_xpoint );
    MPI_Type_free( &mpi_point );
  }

end_1:
  if ( mpi_xtetra )
    MPI_Type_free( &mpi_xtetra );
end_2:
  if ( mpi_tetra )
    MPI_Type_free( &mpi_tetra );
end_3:
  if ( mpi_xpoint )
    MPI_Type_free( &mpi_xpoint );
end_4:
  if ( mpi_point )
    MPI_Type_free( &mpi_point );
end_5:

  return k;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param idgrp index of the group that we process.
 * \param dest rank of destination.
 * \param send_ext_comm \a send_ext_comm[k] is incremented if we must send an
 * external face communicator to the proc k.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Count the number of external face communicators that this will send to each
 * other processor and increment the \a send_ext_comm array.
 *
 */
static inline
int PMMG_count_extFaceComm2send(PMMG_pParMesh parmesh,int idgrp,int dest,
                                int *send_ext_comm) {
  PMMG_pint_comm int_comm;
  PMMG_pext_comm ext_comm;
  int            *intvalues,*face2int_face_comm_idx2,nitem,nitem_int_comm;
  int            k,j;

  int_comm  = parmesh->int_face_comm;
  nitem     = int_comm->nitem;
  intvalues = int_comm->intvalues;

  /** Step 1: intvalues array initialisation by 1 for the faces of the current
   * group (and PMMG_UNSET for the other faces) */
  for ( k=0; k<nitem; ++k ) intvalues[k] = PMMG_UNSET;

  face2int_face_comm_idx2 = parmesh->listgrp[idgrp].face2int_face_comm_index2;
  nitem_int_comm          = parmesh->listgrp[idgrp].nitem_int_face_comm;

  for ( k=0; k<nitem_int_comm; ++k )
    intvalues[face2int_face_comm_idx2[k]] = 1;

  /** Step 3: Process the external communicators and increment \a
   * send_ext_comm[color_out] if we must send an external face communicator to
   * the color_out processor.
   */
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_comm = &parmesh->ext_face_comm[k];

    for ( j=0; j<ext_comm->nitem; ++j ) {
      if ( intvalues[ext_comm->int_comm_index[j]]>0 ) {
        ++send_ext_comm[ext_comm->color_out];
        break;
      }
    }
  }
  return 1;
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
 * Send the wanted part of the face communicator to the suitable processor.
 *
 */
static inline
int PMMG_send_faceCommunicators(PMMG_pParMesh parmesh,int idgrp,int dest,int tag,
                                MPI_Comm comm) {
  PMMG_pint_comm int_comm;
  PMMG_pext_comm ext_comm;
  int            *pos_grp_faces_in_ext_comm,*pos_grp_faces_in_face2int;
  int            *ext_comm_first_idx,*intvalues,intval,nitem;
  int            *face2int_face_comm_idx1,*face2int_face_comm_idx2;
  int            myrank,color_out,nitem_int_comm,nfaces_in_grp;
  int            next_comm,nitem_tot,idx,idx_glob;
  int            idx_next,id_myrank_dest,nitem_ext_comm_myrank_dest;
  int            ier,k,i,j;
  MPI_Request    req[7];
  MPI_Status     status;

  ier    = 0;

  ext_comm_first_idx        = NULL;
  pos_grp_faces_in_ext_comm = NULL;
  pos_grp_faces_in_face2int = NULL;

  /** Process the list of external face communicators to create a new external
   * communicator: if the face belong to this group (position in intvalue is
   * filled), stores the position in the external communicator in a pack
   * array and the intvalues value in another array. Both arrays will be
   * sended by the external communicator. */
  int_comm  = parmesh->int_face_comm;
  nitem     = int_comm->nitem;
  intvalues = int_comm->intvalues;

  /** Step 1: intvalues array initialisation with the face index (ensure that
   * intvalues is strictly greater than 0 for a face of this group and strictly
   * lower for the other faces */
  for ( k=0; k<nitem; ++k )
    intvalues[k] = PMMG_UNSET;

  face2int_face_comm_idx1 = parmesh->listgrp[idgrp].face2int_face_comm_index1;
  face2int_face_comm_idx2 = parmesh->listgrp[idgrp].face2int_face_comm_index2;
  nitem_int_comm          = parmesh->listgrp[idgrp].nitem_int_face_comm;

  /** Step 2: send the list of interface faces of this group */
  MPI_CHECK( MPI_Isend(&nitem_int_comm,1,MPI_INT,dest,tag,
                       parmesh->comm,&req[0]),goto end );

  MPI_CHECK( MPI_Isend(face2int_face_comm_idx1,nitem_int_comm,MPI_INT,dest,tag,
                       parmesh->comm,&req[1]),goto end );

  for ( k=0; k<nitem_int_comm; ++k )
    intvalues[face2int_face_comm_idx2[k]] = face2int_face_comm_idx1[k];

  /** Step 3: Process the external communicators and:
   *
   * - fill \a next_comm by the number of external communicators that need to be
       communicated to the proc dest (because they contains at least 1 face of
       the group to move);
   *
   * - fill \a nitem_tot by the sum of the number of faces of the group in each
       external communicator;
   *
   * - mark the faces that are in the external communicators to be able to
       detect the list of the faces at the interface between two groups of the
       local proc.
   */
  id_myrank_dest             = PMMG_UNSET;
  nitem_ext_comm_myrank_dest = 0;
  next_comm                  = 0;
  nitem_tot                  = 0;
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_comm = &parmesh->ext_face_comm[k];

    nfaces_in_grp = 0;
    if ( ext_comm->color_out == dest ) {
      id_myrank_dest = k;
      nitem_ext_comm_myrank_dest = ext_comm->nitem;
    }

    for ( j=0; j<ext_comm->nitem; ++j ) {
      if ( intvalues[ext_comm->int_comm_index[j]] > 0 ) {
        intvalues[ext_comm->int_comm_index[j]] *= -1;
        ++nfaces_in_grp;
      }
    }
    if ( nfaces_in_grp ) {
      ++next_comm;
      nitem_tot += nfaces_in_grp;
    }
  }

  /** Step 4: Process the face2int_face_comm_index2 array and count the number
   * of faces at the interface of 2 local groups */
  nfaces_in_grp = 0;
  for ( k=0; k<nitem_int_comm; ++k ) {
    assert(intvalues[face2int_face_comm_idx2[k]]);
    if ( intvalues[face2int_face_comm_idx2[k]] > 0 ) {
      ++nfaces_in_grp;
    }
  }

  if ( nfaces_in_grp ) {
    ++next_comm; nitem_tot += nfaces_in_grp;
  }

  /** Step 5:
   *
   * - process the internal faces of the group (if needed) and complete the
   *   external communicator lists by a new communicator containing the list of
   *   the index in \a face2int_face_comm_index1 of the faces that must be added
   *   to the external communicator myrank-dest (and dest-myrank) after the
   *   group displacement;
   *
   * - reset the mark in the intvalues array;
   *
   * - process the external communicators and fill \a ext_comm_first_idx by the
   *   position of the next external communicator in the \a pos_grp_faces_*
   *   arrays, \a pos_grp_faces_in_ext_comm by the position of the group faces
   *   in the external communicator and the destination of the group and \a
   *   pos_grp_faces_in_face2int by the position of the group faces in the
   *   \a face2int_face_comm arrays and the \a color_out field of the matching
   *   external communicator.
   *
   */
  PMMG_MALLOC(parmesh,ext_comm_first_idx,next_comm,int,
              "positions of the external comms in pos_grp_faces arrays",
              goto end);
  PMMG_MALLOC(parmesh,pos_grp_faces_in_ext_comm,nitem_tot+next_comm,int,
              "pos of group faces in the external communicator",
              goto end);
  PMMG_MALLOC(parmesh,pos_grp_faces_in_face2int,nitem_tot+next_comm,int,
              "pos of group faces in the face2int_face_comm arrays",
              goto end);

  idx      = 0;
  idx_glob = 0;
  myrank   = parmesh->myrank;

  if ( nfaces_in_grp ) {
    /* Process the internal faces first: store the position of the faces in the
    external comm and in the face2intface_comm_index array */
    for ( k=0; k<nitem_int_comm; ++k ) {
      assert(intvalues[face2int_face_comm_idx2[k]]);
      if ( intvalues[face2int_face_comm_idx2[k]] > 0 ) {
        pos_grp_faces_in_face2int[idx_glob] = k;

        /* Remark: filling \a pos_grp_faces_in_ext_comm is not mandatory but
         * it allows to have the same size than the array
         * \a pos_grp_faces_in_face2int. */
        pos_grp_faces_in_ext_comm[idx_glob] = nitem_ext_comm_myrank_dest+idx_glob;
        ++idx_glob;
      }
    }
    pos_grp_faces_in_face2int[idx_glob] = myrank;
    pos_grp_faces_in_ext_comm[idx_glob] = myrank;
    ext_comm_first_idx[idx++] = ++idx_glob;
  }

  /* Process the external faces: fill intvalues by the position of the face in
   * the face2int_face_comm arrays, then fill \a pos_grp_faces_in_face2int by
   * this position */
  for ( k=0; k<nitem_int_comm; ++k )
    intvalues[face2int_face_comm_idx2[k]] = k;

  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_comm = &parmesh->ext_face_comm[k];

    nfaces_in_grp = 0;
    for ( j=0; j<ext_comm->nitem; ++j ) {
      intval = intvalues[ext_comm->int_comm_index[j]];
      if ( intval >= 0 ) {
        pos_grp_faces_in_face2int[idx_glob  ]  = intval;
        pos_grp_faces_in_ext_comm[idx_glob++]  = j;
        ++nfaces_in_grp;
      }
    }
    if ( nfaces_in_grp ) {
      pos_grp_faces_in_face2int[idx_glob] = ext_comm->color_out;
      pos_grp_faces_in_ext_comm[idx_glob] = dest;
      ext_comm_first_idx[idx++] = ++idx_glob;
    }
  }
  assert ( idx_glob == nitem_tot+next_comm );
  assert ( idx == next_comm+1 );

  /** Step 6: Send the data to the other processors */
  MPI_CHECK( MPI_Isend(&next_comm,1,MPI_INT,dest,tag,parmesh->comm,&req[2]),
             goto end );
  MPI_CHECK( MPI_Isend(ext_comm_first_idx,next_comm,MPI_INT,dest,tag,
                       parmesh->comm,&req[3]),goto end );
  MPI_CHECK( MPI_Isend(&nitem_tot,1,MPI_INT,dest,tag,parmesh->comm,&req[4]),
             goto end );
  MPI_CHECK( MPI_Isend(pos_grp_faces_in_face2int,nitem_tot+next_comm,MPI_INT,
                       dest,tag,parmesh->comm,&req[5]),goto end );

  /* If we have old internal faces that became external, they are stored in
   * the first position of the array and they have a
   * color_out==myrank. In this case, no need to send it */
  idx_next  = ext_comm_first_idx[0];
  color_out = pos_grp_faces_in_ext_comm[idx_next]-1;
  if ( color_out==myrank ) {
    k=1;
    idx = idx_next;
  }
  else
    k = idx = 0;

  for ( ; k<next_comm; ++k ) {
    idx_next  = ext_comm_first_idx[k];
    nitem     = idx_next-idx;
    color_out = pos_grp_faces_in_face2int[idx_next]-1;

    MPI_CHECK( MPI_Isend(&nitem,1,MPI_INT,color_out,tag,parmesh->comm,&req[6]),
               goto end );
    MPI_CHECK( MPI_Isend(&pos_grp_faces_in_ext_comm[idx],nitem,MPI_INT,color_out,
                        tag,parmesh->comm,&req[7]),goto end );
    idx = idx_next;
  }

  /** Step 7: Update the local external communicators:
   *
   * - Remove the external group faces from the suitable external comms
   *
   * - Add the internal group faces to the suitable external communicator
   */
  idx = 0;
  for ( k=0; k<next_comm; ++k ) {
    idx_next  = ext_comm_first_idx[k];
    nitem     = idx_next-1-idx;
    color_out = pos_grp_faces_in_ext_comm[idx_next]-1;

    if ( color_out == myrank ) {
      /** Add the new external faces */
      if ( id_myrank_dest==PMMG_UNSET ) {
        /* Creation of a new external communicator */
        PMMG_REALLOC(parmesh,parmesh->ext_face_comm,
                     parmesh->next_face_comm+1,parmesh->next_face_comm,
                     PMMG_ext_comm,"list of external communicators",goto end);
        parmesh->next_face_comm = parmesh->next_face_comm+1;

        ext_comm = &parmesh->ext_face_comm[parmesh->next_face_comm];
        ext_comm->color_in  = myrank;
        ext_comm->color_out = dest;
        ext_comm->nitem     = 0;
      }
      else {
        /* Completion of the existing external communicator */
        ext_comm = &parmesh->ext_face_comm[id_myrank_dest];
      }

      /* New interface faces adding */
      PMMG_REALLOC(parmesh,ext_comm->int_comm_index,ext_comm->nitem+nitem,
                   ext_comm->nitem,int,"int_comm_index",goto end);
      ext_comm->nitem += nitem;

      for ( j=0; j<nitem-1; ++j ) {
        i = face2int_face_comm_idx2[pos_grp_faces_in_ext_comm[j+idx]];
        ext_comm->int_comm_index[ext_comm->nitem+j] = i;
        assert ( pos_grp_faces_in_ext_comm[j]==ext_comm->nitem+j );
      }
    }
    else {
      /** Remove the old external faces of the group */
      for ( i=0; i<parmesh->next_face_comm; ++i ) {
        ext_comm = &parmesh->ext_face_comm[i];
        if ( ext_comm->color_out == color_out ) break;
      }
      for ( j=idx; j<idx_next-1; ++j ) {
        ext_comm->int_comm_index[pos_grp_faces_in_ext_comm[j]] = PMMG_UNSET;
      }
    }
    idx = idx_next;
  }

  for ( k=0; k<8; ++k )
    MPI_Wait(&req[k], &status);

  /* Success */
  ier = 1;

end:
  if ( ext_comm_first_idx )
    PMMG_DEL_MEM(parmesh,ext_comm_first_idx,next_comm,int,
                 "positions of the external comms in pos_grp_faces* arrays");

  if ( pos_grp_faces_in_ext_comm )
    PMMG_DEL_MEM(parmesh,pos_grp_faces_in_ext_comm,nitem_tot+next_comm,int,
                 "pos_grp_faces_in_ext_comm");

  if ( pos_grp_faces_in_face2int )
    PMMG_DEL_MEM(parmesh,pos_grp_faces_in_face2int,nitem_tot+next_comm,int,
                 "pos_grp_faces_in_face2int");

  return ier;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param idgrp index of the group that we process.
 * \param source rank of source.
 * \param tag message tag
 * \param comm communicator
 * \param recv_grps number of groups coming from the proc source.
 * \param recv_ext_comm number of external face communicators coming from the
 * proc source
 *
 * \return 0 if fail, 1 otherwise
 *
 * Receive the wanted part of the face communicator from the suitable processor.
 *
 */
static inline
int PMMG_recv_faceCommunicators(PMMG_pParMesh parmesh,int idgrp,int source,
                                int tag,MPI_Comm comm,int recv_grps,
                                int recv_ext_comm)
{
  PMMG_pGrp      grp;
  PMMG_pint_comm int_comm;
  PMMG_pext_comm ext_comm,ext_comm2;
  int            *face2int_face_comm_idx1,*face2int_face_comm_idx2;
  int            *ext_comm_first_idx;
  int            *pos_grp_faces_in_face2int,*pos_grp_faces_in_ext_comm;
  int            nitem_int_comm,next_comm,nitem_tot;
  int            next_comm_prev,next_comm_max,sum_nitem_prev,sum_nitem_max;
  int            nitem,nitem_chk,nitem_max,nitem_prev;
  int            color_out,id_in_out,id_in_source;
  int            ier,k,j,l,idx,idx_next,idx_in_face2int,pos_in_int_comm;
  MPI_Request    req[9];
  MPI_Status     status;


  ier    = 0;

  /** Step 1: Receive and update the communicators of the groups that we receive
   * from the proc source */
  int_comm  = parmesh->int_face_comm;

  id_in_source  = PMMG_UNSET;
  next_comm_max = 0;
  sum_nitem_max = 0;

  ext_comm_first_idx        = NULL;
  pos_grp_faces_in_face2int = NULL;
  pos_grp_faces_in_ext_comm = NULL;

  for ( k=0; k<recv_grps; ++k ) {
    grp = &parmesh->listgrp[idgrp+k];

    /** Step 1-a: receive the list of the interface faces of this group */
    MPI_CHECK( MPI_Irecv(&nitem_int_comm,1,MPI_INT,source,tag,parmesh->comm,
                         &req[0]),goto end );

    assert ( !grp->face2int_face_comm_index1 );
    assert ( !grp->face2int_face_comm_index2 );
    PMMG_MALLOC(parmesh,grp->face2int_face_comm_index1,nitem_int_comm,int,
                "alloc face2int_face_comm_index1 ",goto end);
    PMMG_MALLOC(parmesh,grp->face2int_face_comm_index2,nitem_int_comm,int,
                "alloc face2int_face_comm_index2 ",goto end);

    face2int_face_comm_idx1 = parmesh->listgrp[idgrp].face2int_face_comm_index1;
    face2int_face_comm_idx2 = parmesh->listgrp[idgrp].face2int_face_comm_index2;

    MPI_CHECK( MPI_Irecv(face2int_face_comm_idx1,nitem_int_comm,MPI_INT,
                         source,tag,parmesh->comm,&req[1]),goto end );

#ifndef NDEBUG
    for ( j=0; j<nitem_int_comm; ++j )
      face2int_face_comm_idx2[j] = PMMG_UNSET;
#endif

    /** Step 1-b: receive the list of the external communicators of this group */
    MPI_CHECK( MPI_Irecv(&next_comm,1,MPI_INT,source,tag,parmesh->comm,&req[2]),
               goto end );

    next_comm_prev = next_comm_max;
    next_comm_max  = MG_MAX(next_comm,next_comm_max);
    PMMG_REALLOC(parmesh,ext_comm_first_idx,next_comm_max,next_comm_prev,int,
                "positions of the external comms in grp_ext_comm",
                 next_comm_max=next_comm_prev;goto end);

    MPI_CHECK( MPI_Irecv(ext_comm_first_idx,next_comm,MPI_INT,source,
                        tag,parmesh->comm,&req[3]),goto end );

    MPI_CHECK( MPI_Irecv(&nitem_tot,1,MPI_INT,source,tag,parmesh->comm,&req[4]),
               goto end );

    sum_nitem_prev = sum_nitem_max;
    sum_nitem_max  = MG_MAX(next_comm+nitem_tot,sum_nitem_max);
    PMMG_REALLOC(parmesh,pos_grp_faces_in_face2int,sum_nitem_max,sum_nitem_prev,
                int,"position of the group faces in the face2int* arrays",
                sum_nitem_max=sum_nitem_prev;goto end);

    MPI_CHECK( MPI_Irecv(pos_grp_faces_in_face2int,nitem_tot+next_comm,MPI_INT,
                         source,tag,parmesh->comm,&req[5]),goto end );

    /** Step 1-c: update the list of the external communicators of this group */
    idx       = 0;
    nitem_max = 0;
    for ( j=0; j<next_comm; ++j ) {
      idx_next  = ext_comm_first_idx[j];
      nitem     = idx_next-idx;
      color_out = pos_grp_faces_in_face2int[idx_next-1];

      /* Search the index of the source(=color_in)-color_out communicator */
      for ( id_in_out=0;id_in_out<parmesh->next_face_comm;++id_in_out ) {
        if ( parmesh->ext_face_comm[id_in_out].color_out == color_out ) break;
      }

      if ( color_out == parmesh->myrank ) {
        /* Delete the old external faces between source and myrank (both groups
         * are now on the myrank processor) */

        /* Receive the position of the faces in the external communicator */
        assert ( recv_ext_comm && "Unexpected external communicator.");

        MPI_CHECK( MPI_Irecv(&nitem_chk,1,MPI_INT,color_out,tag,parmesh->comm,
                            &req[6]), goto end );
        assert(nitem_chk==nitem && "Mismatch sizes");

        nitem_prev = nitem_max;
        nitem_max  = MG_MAX(nitem,nitem_max);
        PMMG_REALLOC(parmesh,pos_grp_faces_in_ext_comm,nitem_max,nitem_prev,int,
                     "position of group faces in the ext comm",
                     nitem_max=nitem_prev;goto end);

        MPI_CHECK( MPI_Irecv(&pos_grp_faces_in_ext_comm,nitem,MPI_INT,
                            source,tag,parmesh->comm,&req[7]),goto end );
        assert ( pos_grp_faces_in_ext_comm[nitem-1]==color_out );

        /* External comm in which we will work */
        assert ( id_in_out<parmesh->next_face_comm &&
                 "Received (k-myrank) interfaces from k but (myrank-k)"
                 " interfaces not found on myrank"  );

        ext_comm = &parmesh->ext_face_comm[id_in_out];

        /* Get the position of the face in the internal comm and update \a
         * face2int_face_comm_index2. Delete the face from the external comm. */
        for ( l=0; l<nitem-1; ++l ) {
          idx_in_face2int = pos_grp_faces_in_face2int[l+idx];
          pos_in_int_comm = ext_comm->int_comm_index[pos_grp_faces_in_ext_comm[l]];
          face2int_face_comm_idx2[idx_in_face2int] = pos_in_int_comm;
          ext_comm->int_comm_index[pos_grp_faces_in_ext_comm[l]] = PMMG_UNSET;
        }
      }
      else {
        if ( color_out==source ) id_in_source = id_in_out;

        /* Add the new external faces between source and myrank (the group that
         * was on source is now on myrank) */
        if ( id_in_out==parmesh->next_face_comm ) {
          /* Creation of a new external communicator */
          PMMG_REALLOC(parmesh,parmesh->ext_face_comm,
                       parmesh->next_face_comm+1,parmesh->next_face_comm,
                       PMMG_ext_comm,"list of external communicators",goto end);
          parmesh->next_face_comm = parmesh->next_face_comm+1;

          ext_comm = &parmesh->ext_face_comm[parmesh->next_face_comm];
          ext_comm->color_in  = parmesh->myrank;
          ext_comm->color_out = k;
          ext_comm->nitem     = 0;
        }
        else {
          /* Completion of the existing external communicator */
          ext_comm = &parmesh->ext_face_comm[id_in_out];
        }

        /* New interface faces adding */
        assert(nitem && "nothing to add... weird");
        PMMG_REALLOC(parmesh,ext_comm->int_comm_index,ext_comm->nitem+nitem-1,
                     ext_comm->nitem,int,"int_comm_index",goto end);
        ext_comm->nitem += nitem-1;

        /* Assign a position in the internal comm to the face */
        for ( l=0; l<nitem-1; ++l ) {
          idx_in_face2int = pos_grp_faces_in_face2int[l+idx];
          face2int_face_comm_idx2[idx_in_face2int]    = int_comm->nitem;
          ext_comm->int_comm_index[ext_comm->nitem+l] = int_comm->nitem++;
        }
      }
      idx = idx_next;
    }

#ifndef NDEBUG
    for ( j=0; j<nitem_int_comm; ++j )
      assert( face2int_face_comm_idx2[j]>=0 && "face2int_face_comm not filled.");
#endif

    for ( j=0; j<8; ++j )
      MPI_Wait(&req[j], &status);
  }
  PMMG_DEL_MEM(parmesh,ext_comm_first_idx,next_comm_max,int,
               "positions of the external comms in pos_grp_faces* arrays");

  PMMG_DEL_MEM(parmesh,pos_grp_faces_in_face2int,sum_nitem_max,int,
               "pos_grp_faces_in_face2int");


  /** Step 2: receive the external communicators that we must update because the
   * distant group has moved from proc source to another proc (not myrank). */
  if ( id_in_source==PMMG_UNSET )
    for ( id_in_source=0; id_in_source<parmesh->next_face_comm; ++id_in_source )
      if ( parmesh->ext_face_comm[k].color_out == source ) break;

  assert ( id_in_source!=PMMG_UNSET && "Communicator to update not found" );

  for ( k=0; k<recv_ext_comm; ++k ) {

    /** Step 2-a: receive the external communicator to update */
    MPI_CHECK( MPI_Irecv(&nitem,1,MPI_INT,source,tag,parmesh->comm,&req[8]),
               goto end );
    assert(nitem && "Empty external comm");

    nitem_max = MG_MAX(nitem,nitem_max);
    PMMG_REALLOC(parmesh,pos_grp_faces_in_ext_comm,nitem_max,nitem,int,
                 "position of group faces in the ext comm",goto end);

    MPI_CHECK( MPI_Irecv(&pos_grp_faces_in_ext_comm,nitem,MPI_INT,source,tag,
                        parmesh->comm,&req[9]),goto end );

    /** Step 2-b: update the matching external communicator: delete the faces
     * from the myrank-source communicator and add it in the myrank-color_out
     * one */
    color_out = pos_grp_faces_in_ext_comm[nitem-1];

    /* Find the communicator myrank-color_out */
    for ( id_in_out=0; id_in_out<parmesh->next_face_comm; ++id_in_out )
      if ( parmesh->ext_face_comm[k].color_out == color_out ) break;

    /* Add the new external faces between source and myrank (the group that
     * was on source is now on myrank) */
    if ( id_in_out==parmesh->next_face_comm ) {
      /* Creation of a new external communicator */
      PMMG_REALLOC(parmesh,parmesh->ext_face_comm,
                   parmesh->next_face_comm+1,parmesh->next_face_comm,
                   PMMG_ext_comm,"list of external communicators",goto end);
      parmesh->next_face_comm = parmesh->next_face_comm+1;

      ext_comm = &parmesh->ext_face_comm[parmesh->next_face_comm];
      ext_comm->color_in  = parmesh->myrank;
      ext_comm->color_out = color_out;
      ext_comm->nitem     = 0;
    }
    else {
      /* Completion of the existing external communicator */
      ext_comm = &parmesh->ext_face_comm[id_in_out];
    }

    ext_comm2  = &parmesh->ext_face_comm[id_in_source];
    PMMG_REALLOC(parmesh,ext_comm->int_comm_index,ext_comm->nitem+nitem-1,
                 ext_comm->nitem,int,"int_comm_index",goto end);
    ext_comm->nitem += nitem-1;

    /* Move the faces from ext_comm2 toward ext_comm */
    for ( l=0; l<nitem-1; ++l ) {
      idx = pos_grp_faces_in_ext_comm[l];
      ext_comm->int_comm_index[ext_comm->nitem+l] = ext_comm2->int_comm_index[idx];
      ext_comm2->int_comm_index[idx] = PMMG_UNSET;
    }

    for ( l=8; l<10; ++l )
      MPI_Wait(&req[l], &status);
  }

  /* Success */
  ier = 1;

end:
  if ( ext_comm_first_idx )
    PMMG_DEL_MEM(parmesh,ext_comm_first_idx,next_comm_max,int,
                 "positions of the external comms in pos_grp_faces* arrays");

  if ( pos_grp_faces_in_ext_comm )
    PMMG_DEL_MEM(parmesh,pos_grp_faces_in_ext_comm,nitem_max,int,
                 "pos_grp_faces_in_ext_comm");

  if ( pos_grp_faces_in_face2int )
    PMMG_DEL_MEM(parmesh,pos_grp_faces_in_face2int,sum_nitem_max,int,
                 "pos_grp_faces_in_face2int");
  return ier;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param part pointer toward the metis array containing the partitions.
 *
 * \return 0 if fail but we can try to save a mesh, -1 if we fail and are unable
 * to save the mesh, 1 if we success
 *
 * Send the suitable group to other procs and recieve their groups.
 * Deallocate the \a part array.
 *
 */
/* static inline */
/* int PMMG_mpiexchange_grps(PMMG_pParMesh parmesh,idx_t *part) { */
/*   PMMG_pint_comm int_comm; */
/*   int            nitem; */
/*   int            *send_grps,*recv_grps,*send_ext_comm,*recv_ext_comm,nprocs; */
/*   int            nrecv,ngrps_added,old_ngrp,new_ngrp; */
/*   int            ier,ier_glob,k,tag; */

/*   /\** Merge all the groups that must be send to a given proc into 1 group *\/ */
/*   ier = PMMG_merge_grps2send(parmesh,&part); */
/*   MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm); */

/*   if ( ier_glob < 0 ) { */
/*     fprintf(stderr,"\n  ## Unable to compute the new group partition.\n"); */
/*     ier = -1; */
/*     goto end; */
/*   } */
/*   else if ( !ier_glob ) { */
/*     fprintf(stderr,"\n  ## Unable to compute the new group partition." */
/*             " Try to send it nevertheless\n"); */
/*   } */

/*   /\** Step 1: */
/*    * */
/*    * - Fill the array of the number of groups that we will send to the other */
/*    * processors (\a send_grps) and receive the array of the number of group that */
/*    * we will receive from the other procs (\a recv_grps); */
/*    * */
/*    * - Fill the array of the number of external communicators that we will send */
/*    * to the other processors (\a send_ext_comm) and receive the array of the */
/*    * number of external communicators that we will receive from each other procs */
/*    * (\a recv_ext_comm). */
/*    * */
/*    *\/ */
/*   ier = 0; */

/*   nprocs = parmesh->nprocs; */
/*   PMMG_CALLOC(parmesh,send_grps,nprocs,int,"sended groups to other procs", */
/*               goto end); */
/*   PMMG_CALLOC(parmesh,recv_grps,nprocs,int,"received groups from other procs", */
/*               goto end); */
/*   PMMG_CALLOC(parmesh,send_ext_comm,nprocs*nprocs,int, */
/*               "sended external comm to other procs",goto end); */
/*   PMMG_CALLOC(parmesh,recv_ext_comm,nprocs*nprocs,int, */
/*               "received external comm from other procs",goto end); */

/*   int_comm = parmesh->int_face_comm; */
/*   nitem    = int_comm->nitem; */

/*   PMMG_MALLOC(parmesh,int_comm->intvalues,nitem,int,"intvalues",goto end); */

/*   for ( k=0; k<parmesh->ngrp; ++k ) { */
/*     if ( part[k] == parmesh->myrank ) continue; */
/*     ++send_grps[part[k]]; */

/*     if ( !PMMG_count_extFaceComm2send(parmesh,k,part[k],send_ext_comm) ) { */
/*       fprintf(stderr,"\n  ## Unable to compute the map of processors to which" */
/*               " the group k nees to send data."); */
/*       goto end; */
/*     } */
/*   } */

/*   MPI_CHECK( MPI_Alltoall(send_grps,1,MPI_INT,recv_grps,1,MPI_INT,parmesh->comm), */
/*              goto end); */
/*   MPI_CHECK( MPI_Alltoall(send_ext_comm,1,MPI_INT,recv_ext_comm,1,MPI_INT, */
/*                           parmesh->comm),goto end); */


/*   /\** Count the final number of groups and reallocate the listgrp array *\/ */
/*   ngrps_added = 0; */
/*   for ( k=0; k<nprocs; ++k ) */
/*     ngrps_added += recv_grps[k]; */

/*   PMMG_RECALLOC(parmesh,parmesh->listgrp,ngrps_added+parmesh->ngrp,parmesh->ngrp, */
/*                 PMMG_Grp,"listgrps",goto end); */
/*   old_ngrp      = parmesh->ngrp; */
/*   parmesh->ngrp = ngrps_added+parmesh->ngrp+ngrps_added; */

/*   /\** Step 4: Receive the groups *\/ */
/*   for ( k=0; k<nprocs; ++k ) { */
/*       tag   = parmesh->myrank*nprocs + k; */
/*       nrecv = recv_grps[k]; */

/*       if ( nrecv != PMMG_recv_mesh(parmesh,old_ngrp,k,tag,parmesh->comm,nrecv) ) */
/*         goto end; */

/* #warning if we fail here, we would like to update the communicators of the moved meshes */
/*       /\* if ( !PMMG_recv_faceCommunicators(parmesh,new_ngrp,k,tag,parmesh->comm, *\/ */
/*       /\*                                   nrecv,recv_ext_comm[k]) ) *\/ */
/*       /\*   goto end; *\/ */
/*   } */

/*   /\** Step 2: Send each group to the suitable processor *\/ */
/*   for ( k=0; k<parmesh->ngrp; ++k ) { */
/*     if ( part[k] == parmesh->myrank ) continue; */

/*     tag = part[k]*nprocs + parmesh->myrank; */

/*     if ( !PMMG_send_mesh(parmesh,k,part[k],tag,parmesh->comm) ) */
/*       goto end; */

/*     /\* if ( !PMMG_send_faceCommunicators(parmesh,k,part[k],tag,parmesh->comm) ) *\/ */
/*     /\*   goto end; *\/ */
/*   } */
/*   PMMG_DEL_MEM(parmesh,part,parmesh->ngrp,idx_t,"deallocate parmetis partition"); */

/*   /\** Step 3: Pack the used groups on the mesh *\/ */
/*   new_ngrp = PMMG_pack_grps_norealloc( parmesh,&parmesh->listgrp,parmesh->ngrp ); */

/* /\*   /\\** Count the final number of groups and reallocate the listgrp array *\\/ *\/ */
/* /\*   ngrps_added = 0; *\/ */
/* /\*   for ( k=0; k<nprocs; ++k ) *\/ */
/* /\*     ngrps_added += recv_grps[k]; *\/ */

/* /\*   PMMG_RECALLOC(parmesh,parmesh->listgrp,ngrps_added+new_ngrp,parmesh->ngrp, *\/ */
/* /\*                 PMMG_Grp,"listgrps",goto end); *\/ */
/* /\*   parmesh->ngrp = ngrps_added+new_ngrp; *\/ */

/* /\*   /\\** Step 4: Receive the groups *\\/ *\/ */
/* /\*   for ( k=0; k<nprocs; ++k ) { *\/ */
/* /\*       tag   = parmesh->myrank*nprocs + k; *\/ */
/* /\*       nrecv = recv_grps[k]; *\/ */

/* /\*       if ( nrecv != PMMG_recv_mesh(parmesh,new_ngrp,k,tag,parmesh->comm,nrecv) ) *\/ */
/* /\*         goto end; *\/ */

/* /\* #warning if we fail here, we would like to update the communicators of the moved meshes *\/ */
/* /\*       /\\* if ( !PMMG_recv_faceCommunicators(parmesh,new_ngrp,k,tag,parmesh->comm, *\\/ *\/ */
/* /\*       /\\*                                   nrecv,recv_ext_comm[k]) ) *\\/ *\/ */
/* /\*       /\\*   goto end; *\\/ *\/ */
/* /\*   } *\/ */

/*   /\* Pack the face communicators *\/ */
/*   if ( !PMMG_pack_faceCommunicators(parmesh) ) goto end; */

/*   /\* Node communicators reconstruction from the face ones *\/ */
/*   if ( !PMMG_build_nodeCommFromFaces(parmesh) ) { */
/*     fprintf(stderr,"\n  ## Unable to build the new node communicators from" */
/*             " the face ones.\n"); */
/*     ier = -1; */
/*     goto end; */
/*   } */

/*   /\* Pack the communicators *\/ */
/*   // if ( !PMMG_pack_nodeCommunicators(parmesh) ) goto end; */

/*   if ( !PMMG_parmesh_updateMemMax(parmesh, 105, 1) ) { */
/*     fprintf(stderr,"\n  ## Error: %s: Unable to update the memory repartition" */
/*             " between meshes and communicators.\n",__func__); */
/*     goto end; */
/*   } */

/*   /\** Success *\/ */
/*   ier = 1; */
/*   puts("NO DEADLOCK!!!!"); */

/* end: */
/*   if ( part ) */
/*     PMMG_DEL_MEM(parmesh,part,parmesh->ngrp,idx_t,"deallocate parmetis partition"); */
/*   if ( send_grps ) */
/*     PMMG_DEL_MEM(parmesh,send_grps,nprocs,int,"send_grps"); */
/*   if ( recv_grps ) */
/*     PMMG_DEL_MEM(parmesh,recv_grps,nprocs,int,"recv_grps"); */
/*   if ( send_ext_comm ) */
/*     PMMG_DEL_MEM(parmesh,send_ext_comm,nprocs*nprocs,int,"send_ext_comm"); */
/*   if ( recv_ext_comm ) */
/*     PMMG_DEL_MEM(parmesh,recv_ext_comm,nprocs*nprocs,int,"recv_ext_comm"); */
/*   if ( parmesh->int_face_comm->intvalues ) */
/*     PMMG_DEL_MEM(parmesh,parmesh->int_face_comm->intvalues,nitem,int,"intvalues"); */

/*   return ier; */
/* } */

/**
 * \param parmesh pointer toward the mesh structure.
 *
 * \return -1 if we fail and can not save the meshes, 0 if we fail but can save
 * the meshes, 1 otherwise
 *
 * Compute the group distribution that load balance the groups over the
 * processors and send and recieve the groups from the other processors.
 *
 */
int PMMG_distribute_grps( PMMG_pParMesh parmesh ) {
  idx_t *part;
  int    ier,ier_glob;

  /** Get the new partition of groups (1 group = 1 metis node) */
  part = NULL;
  PMMG_CALLOC(parmesh,part,parmesh->ngrp,idx_t,"allocate parmetis buffer",
              return 0);

  ier = PMMG_part_parmeshGrps2parmetis(parmesh,part,parmesh->nprocs);
  if ( !ier )
    fprintf(stderr,"\n  ## Unable to compute the new group partition.\n");

  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);

  if ( !ier_glob ) {
    PMMG_DEL_MEM(parmesh,part,parmesh->ngrp,idx_t,"deallocate parmetis partition");
    return ier_glob;
  }

  /** Send the suitable groups to other procs and recieve their groups */
  ier = 1;//PMMG_mpiexchange_grps(parmesh,part);
  if ( ier <= 0 )
    fprintf(stderr,"\n  ## Unable to communicate groups through processors.\n");

  return ier;
}

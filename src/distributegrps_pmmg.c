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
 * \param parmesh pointer toward the mesh structure.
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
    PMMG_REALLOC(parmesh,ext_node_comm->int_comm_index,nitem_ext,
                 ext_node_comm->nitem,int,"int_comm_index",return 0);
    ext_node_comm->nitem = nitem_ext;
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
  PMMG_REALLOC(parmesh,parmesh->ext_node_comm,i,
               parmesh->next_node_comm,PMMG_Ext_comm,"ext_node_comm",return 0);
  parmesh->next_node_comm = i;

  /** Step 6: unallocate intvalues array and set the nitem field of the internal
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
int PMMG_pack_faceCommunicators(PMMG_pParMesh parmesh) {
  PMMG_pGrp      grp;
  PMMG_pInt_comm int_face_comm;
  PMMG_pExt_comm ext_face_comm;
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

  /** Step 5: Remove the empty external communicators */
  i = 0;
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];

    if ( !ext_face_comm->nitem ) continue;

    if ( i!=k ) {
      parmesh->ext_face_comm[i].nitem          = ext_face_comm->nitem;
      parmesh->ext_face_comm[i].color_out      = ext_face_comm->color_out;
      parmesh->ext_face_comm[i].int_comm_index = ext_face_comm->int_comm_index;
    }
    ++i;
  }
  PMMG_REALLOC(parmesh,parmesh->ext_face_comm,i,
               parmesh->next_face_comm,PMMG_Ext_comm,"ext_face_comm",return 0);
  parmesh->next_face_comm = i;

  /** Step 6: unallocate intvalues array and set the nitem field of the internal
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
  ngrpOld = ngrp;
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
  PMMG_Int_comm *int_node_comm,*int_face_comm;
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
  if ( !PMMG_pack_grpsAndPart( parmesh,&grps,&ngrp,part,parmesh->ngrp ) ) ier = -1;
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

  assert ( PMMG_check_extFaceComm(parmesh) );
  assert ( PMMG_check_extNodeComm(parmesh) );

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
 * \param tag communication tag (updated by the number of comms
 * \param size array to store the mesh size data
 * \param ismet 1 if we must received a metric field
 *
 * \return 1 if success, 0 if fail.
 *
 * Receive the needed info to allocate a mesh and a metric and perform the
 * wanted allocations.
 *
 */
static
int PMMG_mpirecv_meshSize( PMMG_pParMesh parmesh,int grp_id,int source,int *tag,
                           int size[7],int *ismet )
{
  MMG5_pMesh     mesh;
  MMG5_pSol      met;
  MPI_Status     status;

  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&(parmesh->listgrp[grp_id].mesh),
                  MMG5_ARG_ppMet ,&(parmesh->listgrp[grp_id].met),
                  MMG5_ARG_end);

  mesh = parmesh->listgrp[grp_id].mesh;

  /** Get the mesh maximal authorized memory */
  MPI_CHECK( MPI_Recv( &mesh->memMax,1,MPI_LONG_LONG,source,*tag,parmesh->comm,
               &status),return 0 );
  ++(*tag);

  /** Get the mesh and met size:
   * size[]=[mesh->np,mesh->ne,mesh->xp,mesh->xt,met->size,met->type,ismet] */
  MPI_CHECK( MPI_Recv( size,7,MPI_INT,source,*tag,parmesh->comm,&status),
             return 0 );
  ++(*tag);

  /** Mesh and metric allocations */
  mesh->memCur = 0;
  if ( !PMMG_grpSplit_setMeshSize( mesh,size[0],size[1],0,size[2],size[3] ) )
    return 0;

  met       = parmesh->listgrp[grp_id].met;
  met->size = size[4];
  met->type = size[5];
  *ismet    = size[6];

  if ( !MMG3D_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,met->type) ) return 0;

  /** Mark the mesh as not scaled */
  mesh->info.delta = 1.;

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param grp_id id of the current group
 * \param dest id of the proc to which we send data
 * \param tag tag of the first communication (updated by the number of comms)
 * \param size array to store the mesh size data
 *
 * \return 1 if success, 0 if fail.
 *
 * Fit the mesh and metric to its minimal size and send the needed info to
 * allocate its.
 *
 */
static
int PMMG_mpisend_meshSize( PMMG_pParMesh parmesh,int grp_id,int dest,int *tag,
                           int size[7] )
{
  MMG5_pMesh     mesh;
  MMG5_pSol      met;
  MPI_Comm       comm;

  mesh = parmesh->listgrp[grp_id].mesh;
  met  = parmesh->listgrp[grp_id].met;
  comm = parmesh->comm;

  /** Fit the mesh and metric at their minimal sizes */
  PMMG_REALLOC(mesh,mesh->point,mesh->np+1,mesh->npmax+1,MMG5_Point,
               "points",return 0);
  mesh->npmax = mesh->np;
  PMMG_REALLOC(mesh,mesh->xpoint,mesh->xp+1,mesh->xpmax+1,MMG5_xPoint,
               "xpoints",return 0);
  mesh->xpmax = mesh->xp;
  PMMG_REALLOC(mesh,mesh->tetra,mesh->ne+1,mesh->nemax+1,MMG5_Tetra,
               "tetra",return 0);
  mesh->nemax = mesh->ne;
  PMMG_REALLOC(mesh,mesh->xtetra,mesh->xt+1,mesh->xtmax+1,MMG5_xTetra,
               "xtetra",return 0);
  mesh->xtmax = mesh->xt;

  assert(met->npmax==mesh->npmax);
  assert(met->np   == met->np);
  PMMG_REALLOC(mesh,met->m,met->size*met->np+1,met->size*met->npmax+1,double,
               "metric array",return 0);
  met->npmax = met->np;

  /** Send the mesh/met size:
   * meshsize[]=[mesh->np,mesh->ne,mesh->xp,mesh->xt,
   *             met->size,met->type,ismet] */
  mesh->memMax = mesh->memCur;
  MPI_CHECK( MPI_Send( &mesh->memMax,1,MPI_LONG_LONG,dest,*tag,comm ),
             return 0 );
  ++(*tag);

  size[0] = mesh->np;
  size[1] = mesh->ne;
  size[2] = mesh->xp;
  size[3] = mesh->xt;
  size[4] = met->size;
  size[5] = met->type;
  size[6] = met->m ? 1:0;

  MPI_CHECK( MPI_Send( size,7,MPI_INT,dest,*tag,comm ), return 0 );
  ++(*tag);

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure.
 * \param dest rank of destination.
 * \param nsend number of mesh to send from \a source to \a dest
 * \param part pointer toward the metis array containing the partitions.
 *
 * \param next_comm2send number of external communicator to send to each proc
 *
 * \param ext_comms_next_idx[k] position of the external communicators to send
 * to the proc \a k+1 in the \a extComm_next_idx array.
 *
 * \param nitems2send number of items to send to each proc
 *
 * \param items_next_idx positions of the items to send to each proc in the
 * \a extComm_grpFaces* arrays
 *
 * \param extComm_next_idx \a extComm_next_idx[k] position of the external
 * communicator number \a k+1 in the part of the \a extComm_grpFaces* arrays
 * that goes to a given proc (so for each proc, the communicator number 0 is at
 * position 0)
 *
 * \param extComm_grpFaces2face2int for each group, position of the group faces
 * that belongs to the external communicators in the face2int* arrays
 *
 * \param extComm_grpFaces2extComm for each group, position of the group faces
 * that belongs to the external communicators in the external communicators.
 *
 * \return -1 if fail, the number of groups successfully sended otherwise
 *
 * Send mesh, metric and internal face communicator of the groups that must go
 * to proc \a dest. Last, send the list of external communicators that must go
 * to \a dest
 *
 * \warning the mesh prisms are not treated.
 *
 * \warning IMPROVEMENT: do not use an mpi type, just pack the
 * point+xpoint+tetra+xtetra integer and char fields into one unique array and
 * the point+xpoint double fields into another one.
 *
 * \warning IMPROVEMENT: is it faster to communicate the xpoint+xtetra or the
 * loose it and build it again with a light analysis step?
 *
 */
static
int PMMG_send_grp( PMMG_pParMesh parmesh,int dest,
                   int nsend,idx_t *part,int *next_comm2send,
                   int *ext_comms_next_idx,int *nitems2send,int *items_next_idx,
                   int *extComm_next_idx,int *extComm_grpFaces2face2int,
                   int *extComm_grpFaces2extComm)
{
  PMMG_pExt_comm ext_comm,ext_comm_dest;
  MMG5_pMesh     mesh;
  MMG5_pSol      met;
  MPI_Datatype   mpi_point,mpi_xpoint,mpi_tetra,mpi_xtetra;
  MPI_Comm       comm;
  int            *face2int_face_comm_idx1,*face2int_face_comm_idx2,nitem;
  int            size[7],send_count,send_count2,tag,myrank,color_out;
  int            grp_id_prev,grp_id,k,i,j,idx,idx_next,ext_comm_idx,item_idx;
  int            ier,ismet,old_nitem;

  if ( !nsend ) return 0;

  ier = 1;
  send_count = 0;

  comm = parmesh->comm;

  mpi_point  = NULL;
  mpi_xpoint = NULL;
  mpi_tetra  = NULL;
  mpi_xtetra = NULL;

  tag = MPI_SENDGRP_TAG;

  grp_id = 0;
  while ( send_count<nsend ) {

    /* Id of the group to send */
    while ( part[grp_id] != dest ) {
      ++grp_id;
    }

    mesh = parmesh->listgrp[grp_id].mesh;
    met  = parmesh->listgrp[grp_id].met;

    /** Fit the mesh and metric at their minimal sizes and send the mesh sizes */
    if ( !PMMG_mpisend_meshSize( parmesh,grp_id,dest,&tag,size ) ) ier=0;

    /** Send mesh */
    /* Create mesh type */
    if ( !PMMG_create_MPI_Point( mesh->point,  &mpi_point ) ) ier=0;
    if ( !PMMG_create_MPI_xPoint(mesh->xpoint, &mpi_xpoint) ) ier=0;
    if ( !PMMG_create_MPI_Tetra( mesh->tetra,  &mpi_tetra ) ) ier=0;
    if ( !PMMG_create_MPI_xTetra(mesh->xtetra, &mpi_xtetra) ) ier=0;

    MPI_CHECK( MPI_Send( &mesh->point[1] ,mesh->np,mpi_point,dest,tag,comm),
               ier=0 );

    ++tag;
    MPI_CHECK( MPI_Send( &mesh->xpoint[1],mesh->xp,mpi_xpoint,dest,tag,comm),
               ier=0 );

    ++tag;
    MPI_CHECK( MPI_Send( &mesh->tetra[1] ,mesh->ne,mpi_tetra,dest,tag,comm),
               ier=0 );

    ++tag;
    MPI_CHECK( MPI_Send( &mesh->xtetra[1],mesh->xt,mpi_xtetra,dest,tag,comm),
               ier=0 );
    ++tag;

    ismet = met->m ? 1 : 0;
    if ( ismet ) {
      MPI_CHECK( MPI_Send( met->m,met->size*met->np,MPI_DOUBLE,dest,tag,comm),
                 ier=0 );
      ++tag;
    }

    PMMG_Free_MPI_meshDatatype( &mpi_point,&mpi_xpoint,&mpi_tetra,&mpi_xtetra );

    /** Free useless mesh */
    /* 1: mesh */
    PMMG_DEL_MEM(mesh,mesh->point, mesh->npmax+1,MMG5_Point,"point");
    PMMG_DEL_MEM(mesh,mesh->tetra, mesh->nemax+1,MMG5_Tetra,"tetra");
    PMMG_DEL_MEM(mesh,mesh->xpoint,mesh->xpmax+1,MMG5_xPoint,"xpoint");
    PMMG_DEL_MEM(mesh,mesh->xtetra,mesh->xtmax+1,MMG5_xTetra,"xtetra");
    if ( mesh->adja ) {
      PMMG_DEL_MEM( mesh,mesh->adja,4*mesh->nemax+5,int,"adjacency table");
    }
    if ( ismet )
      PMMG_DEL_MEM(mesh,met->m,(met->npmax+1)*met->size,double,"metric");

    mesh->np = mesh->npmax = 0;
    mesh->xp = mesh->xpmax = 0;
    mesh->ne = mesh->nemax = 0;
    mesh->xt = mesh->xtmax = 0;

    MMG3D_Free_all( MMG5_ARG_start,
                    MMG5_ARG_ppMesh, &parmesh->listgrp[grp_id].mesh,
                    MMG5_ARG_ppMet, &parmesh->listgrp[grp_id].met,
                    MMG5_ARG_end );

    /** Send the internal face communicator of the group */
    face2int_face_comm_idx1 = parmesh->listgrp[grp_id].face2int_face_comm_index1;
    nitem                   = parmesh->listgrp[grp_id].nitem_int_face_comm;
    MPI_CHECK( MPI_Send(face2int_face_comm_idx1,nitem,MPI_INT,dest,tag,
                        parmesh->comm),ier = 0 );
    ++tag;

    /* If we fail to send this group, try to keep things clean to be able to
     * save something */
    if ( !ier ) break;

    ++send_count;
  }

  /** Send the list of external commnicators that must go on dest because the
   * groups are going on dest and update the local external communicators:
   *
   * - Remove the external group faces from the suitable external comms
   *
   * - Add the internal group faces to the suitable external communicator
   */
  ier = 1;
  send_count2 = 0;

  /* Position of the external communicators in the global arrays */
  if ( !dest ) {
    ext_comm_idx = 0;
    item_idx     = 0;
  }
  else {
    ext_comm_idx = ext_comms_next_idx[dest-1];
    item_idx     =     items_next_idx[dest-1];
  }

  MPI_CHECK( MPI_Send(&extComm_next_idx[ext_comm_idx],next_comm2send[dest],
                      MPI_INT,dest,tag,comm),ier = 0 );
  ++tag;

  MPI_CHECK( MPI_Send(&extComm_grpFaces2face2int[item_idx],nitems2send[dest],
                      MPI_INT,dest,tag,comm),ier = 0 );
  ++tag;

  /* Update of the communicator */
  /* Search the external comm myrank-dest */
  ext_comm_dest = NULL;
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    if ( parmesh->ext_face_comm[k].color_out == dest ) {
      ext_comm_dest = &parmesh->ext_face_comm[k];
      break;
    }
  }

  /* Process the external communicators */
  myrank   = parmesh->myrank;
  idx_next = item_idx;
  grp_id   = PMMG_UNSET;

  for ( k=0; k<next_comm2send[dest]; ++k ) {
    idx         = idx_next;
    idx_next    = item_idx + extComm_next_idx[ext_comm_idx+k];
    nitem       = idx_next-idx-2;

    grp_id_prev = grp_id;
    grp_id      = extComm_grpFaces2extComm[idx+nitem];
    /* Color_out of the external communicator to be able to detect the faces
     * to remove (myrank-dest) */
    color_out   = extComm_grpFaces2face2int[idx+nitem+1];

    if ( grp_id != grp_id_prev ) {
      /* Avoid the update of the communicator for a group that has not been
       * sended */
      if ( send_count2==send_count ) break;
    }

    if ( color_out == myrank ) {

      /* Add the new external faces */
      if ( !ext_comm_dest ) {
        /* Creation of a new external communicator */
        PMMG_REALLOC(parmesh,parmesh->ext_face_comm,
                     parmesh->next_face_comm+1,parmesh->next_face_comm,
                     PMMG_Ext_comm,"list of external communicators",ier=-1);
        parmesh->next_face_comm = parmesh->next_face_comm+1;

        ext_comm = &parmesh->ext_face_comm[parmesh->next_face_comm-1];
        ext_comm->color_in       = myrank;
        ext_comm->color_out      = dest;
        ext_comm->nitem          = 0;
        ext_comm->int_comm_index = NULL;
        ext_comm_dest            = ext_comm;
      }
      else ext_comm = ext_comm_dest;

      /* New interface faces adding */
      face2int_face_comm_idx2 = parmesh->listgrp[grp_id].face2int_face_comm_index2;

      PMMG_REALLOC(parmesh,ext_comm->int_comm_index,ext_comm->nitem+nitem,
                   ext_comm->nitem,int,"int_comm_index",ier=-1);
      old_nitem = ext_comm->nitem;
      ext_comm->nitem += nitem;

      for ( j=0; j<nitem; ++j ) {
        i = face2int_face_comm_idx2[extComm_grpFaces2face2int[j+idx-item_idx]];
        ext_comm->int_comm_index[old_nitem+j] = i;
      }
    }
    else {
      if ( color_out == dest ) {
        /* We need to send the extComm_grpFaces2extComm array to allow dest to
         * remove the old external faces from its communicator */
        MPI_CHECK( MPI_Send(&extComm_grpFaces2extComm[idx],nitem+2,
                            MPI_INT,dest,tag,comm),ier = 0 );
        ++tag;
      }

      /* Remove the old external faces of the group */
      for ( i=0; i<parmesh->next_face_comm; ++i ) {
        ext_comm = &parmesh->ext_face_comm[i];
        if ( ext_comm->color_out == color_out ) break;
      }
      for ( j=idx; j<idx_next-2; ++j ) {
        ext_comm->int_comm_index[extComm_grpFaces2extComm[j]] = PMMG_UNSET;
      }
    }

    if ( grp_id != grp_id_prev ) ++send_count2;

    if ( ier<=0 ) goto end;
  }

  /* Success */
  ier = 1;

end:
  PMMG_Free_MPI_meshDatatype( &mpi_point,&mpi_xpoint,&mpi_tetra,&mpi_xtetra );

  if ( ier<0 || send_count2!=send_count )
    return -1;
  else
    return send_count;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param first_grp_id id of the first available group (updated)
 * \param source rank of source.
 * \param nrecv number of groups to receive from the proc source.
 *
 * \return -1 if fail, the number of groups successfully received (so nrecv if
 * success)
 *
 * Receive mesh, metric and internal face communicator from proc \a source.
 * Last, recv the list of external communicators that comes from \a dest.
 *
 * \remark the prism communication is not yet implemented
 *
 * \warning IMPROVEMENT: try to move the creation of mpi_types outside the loop
 *
 */
int PMMG_recv_grp( PMMG_pParMesh parmesh,int *first_grp_id,int source,int nrecv)
{
  PMMG_pExt_comm ext_comm_col_out_ptr,ext_comm;
  PMMG_pInt_comm int_comm;
  MMG5_pMesh     mesh;
  MMG5_pSol      met;
  MPI_Comm       comm;
  MPI_Datatype   mpi_point,mpi_xpoint,mpi_tetra,mpi_xtetra;
  int            myrank;
  int            *face2int_face_comm_idx1,*face2int_face_comm_idx2,nitem;
  int            *extComm_next_idx_recv,*extComm_grpFaces2face2int_recv;
  int            *extComm_grpFaces2extComm_recv;
  int            size[7],tag,next_comm2recv,nitems2recv,nitem_recv,color_out;
  int            grp_id,ier,ismet,k,j,irecv,cur_grp_id,idx,idx_next,grp_idx;
  int            idx_in_face2int,pos_in_int_comm,old_nitem;
  MPI_Status     status;

  if ( !nrecv ) return 0;

  ier = -1;

  comm   = parmesh->comm;
  myrank = parmesh->myrank;

  mpi_xtetra = NULL;
  mpi_tetra  = NULL;
  mpi_xpoint = NULL;
  mpi_point  = NULL;

  extComm_next_idx_recv          = NULL;
  extComm_grpFaces2face2int_recv = NULL;

  tag = MPI_SENDGRP_TAG;

  for ( irecv=0; irecv<nrecv; ++irecv ) {

    grp_id = (*first_grp_id)+irecv;

    /** Get the mesh size */
    if ( !PMMG_mpirecv_meshSize( parmesh,grp_id,source,&tag,size,&ismet ) )
      goto end;

    mesh = parmesh->listgrp[grp_id].mesh;
    met  = parmesh->listgrp[grp_id].met;

    /* Create mesh type */
    if ( !PMMG_create_MPI_Point( mesh->point,  &mpi_point ) ) goto end;
    if ( !PMMG_create_MPI_xPoint(mesh->xpoint, &mpi_xpoint) ) goto end;
    if ( !PMMG_create_MPI_Tetra( mesh->tetra,  &mpi_tetra ) ) goto end;
    if ( !PMMG_create_MPI_xTetra(mesh->xtetra, &mpi_xtetra) ) goto end;

    /** Recv mesh */
    MPI_CHECK( MPI_Recv(&mesh->point[1],mesh->np,mpi_point,source,tag,comm,
                        &status),goto end);
    ++tag;
    MPI_CHECK( MPI_Recv(&mesh->xpoint[1],mesh->xp,mpi_xpoint,source,tag,comm,
                        &status),goto end);
    ++tag;
    MPI_CHECK( MPI_Recv( &mesh->tetra[1] ,mesh->ne,mpi_tetra, source,tag,comm,
                         &status),goto end);
    ++tag;
    MPI_CHECK( MPI_Recv( &mesh->xtetra[1],mesh->xt,mpi_xtetra,source,tag,comm,
                         &status),goto end);
    ++tag;

    /** Recv metric */
    if ( ismet ) {
      MPI_CHECK( MPI_Recv( met->m,met->size*met->np,MPI_DOUBLE,source,tag,comm,
                   &status),goto end );
      ++tag;
    }

    PMMG_Free_MPI_meshDatatype( &mpi_point,&mpi_xpoint,&mpi_tetra,&mpi_xtetra );

    /** Recv the internal face communicator of the group */
    MPI_Probe(source,tag,comm,&status);
    MPI_Get_count(&status,MPI_INT,&nitem);

    PMMG_MALLOC(parmesh,parmesh->listgrp[grp_id].face2int_face_comm_index1,nitem,
                int,"alloc face2int_face_comm_index1 ",goto end);
    PMMG_MALLOC(parmesh,parmesh->listgrp[grp_id].face2int_face_comm_index2,nitem,
                int,"alloc face2int_face_comm_index2 ",goto end);
    parmesh->listgrp[grp_id].nitem_int_face_comm = nitem;

    face2int_face_comm_idx1 = parmesh->listgrp[grp_id].face2int_face_comm_index1;

    MPI_CHECK( MPI_Recv(face2int_face_comm_idx1,nitem,MPI_INT,source,tag,comm,
                        &status),goto end );

#ifndef NDEBUG
    face2int_face_comm_idx2 = parmesh->listgrp[grp_id].face2int_face_comm_index2;
    for ( j=0; j<parmesh->listgrp[grp_id].nitem_int_face_comm; ++j )
      face2int_face_comm_idx2[j] = PMMG_UNSET;
#endif

    ++tag;
  }

  /** Recv the list of external commnicators that comes from source because
   * source send us groups */
  MPI_Probe(source,tag,comm,&status);
  MPI_Get_count(&status,MPI_INT,&next_comm2recv);

  PMMG_MALLOC(parmesh,extComm_next_idx_recv,next_comm2recv,int,
              "extComm_next_idx",goto end);

  MPI_CHECK( MPI_Recv(extComm_next_idx_recv,next_comm2recv,
                      MPI_INT,source,tag,comm,&status),ier = 0 );
  ++tag;

  MPI_Probe(source,tag,parmesh->comm,&status);
  MPI_Get_count(&status,MPI_INT,&nitems2recv);

  PMMG_MALLOC(parmesh,extComm_grpFaces2face2int_recv,nitems2recv,int,
              "extComm_grpFaces2face2int",goto end);

  MPI_CHECK( MPI_Recv(extComm_grpFaces2face2int_recv,nitems2recv,
                      MPI_INT,source,tag,comm,&status),ier = 0 );
  ++tag;

  grp_id      = (*first_grp_id);
  cur_grp_id  = extComm_grpFaces2face2int_recv[extComm_next_idx_recv[0]-2];

  int_comm  = parmesh->int_face_comm;
  idx       = 0;
  for ( j=0; j<next_comm2recv; ++j ) {
    idx_next  = extComm_next_idx_recv[j];
    nitem     = idx_next-idx;
    grp_idx   = extComm_grpFaces2face2int_recv[idx_next-2];
    color_out = extComm_grpFaces2face2int_recv[idx_next-1];

    /* Group that we want to update */
    if ( grp_idx != cur_grp_id ) ++grp_id;
    cur_grp_id = grp_idx;

    face2int_face_comm_idx2 = parmesh->listgrp[grp_id].face2int_face_comm_index2;
    face2int_face_comm_idx1 = parmesh->listgrp[grp_id].face2int_face_comm_index1;

    if ( color_out == myrank ) {
      /* Delete the old external faces between source and myrank (both groups
       * are now on the myrank processor) */

      /* Search the myrank-source communicator */
      ext_comm = NULL;
      for ( k=0; k<parmesh->next_face_comm; ++k ) {
        if ( parmesh->ext_face_comm[k].color_out == source ) {
          ext_comm = &parmesh->ext_face_comm[k];
          break;
        }
      }

      assert ( ext_comm &&
               "Faces to delete in a not found communicator" );

      /* Receive the additional data nedded to treat this case : the position of
       * the old faces in the external communicator */
      MPI_Probe(source,tag,parmesh->comm,&status);
      MPI_Get_count(&status,MPI_INT,&nitem_recv);

      PMMG_MALLOC(parmesh,extComm_grpFaces2extComm_recv,nitem_recv,int,
                  "extComm_grpFaces2extComm",goto end);

      MPI_CHECK( MPI_Recv(extComm_grpFaces2extComm_recv,nitem_recv,
                          MPI_INT,source,tag,comm,&status),ier = 0 );
      ++tag;

      assert ( extComm_grpFaces2extComm_recv[nitem_recv-1]==color_out );

      /* Get the position of the face in the internal comm and update \a
       * face2int_face_comm_index2. Delete the face from the external comm. */
      for ( k=0; k<nitem_recv-2; ++k ) {
        idx_in_face2int = extComm_grpFaces2face2int_recv[k+idx];
        pos_in_int_comm = ext_comm->int_comm_index[extComm_grpFaces2extComm_recv[k]];

        if ( pos_in_int_comm < 0 ) {
          /* This position has been deleted by the migration of another group */
          pos_in_int_comm = int_comm->nitem++;
        }

        face2int_face_comm_idx2[idx_in_face2int] = pos_in_int_comm;
        ext_comm->int_comm_index[extComm_grpFaces2extComm_recv[k]] = PMMG_UNSET;
      }
      PMMG_DEL_MEM(parmesh,extComm_grpFaces2extComm_recv,nitem_recv,int,
                   "extComm_grpFaces2extComm");

    }
    else {
      /* Search the index of the myrank-color_out communicator */
      for ( k=0; k<parmesh->next_face_comm; ++k ) {
        if ( parmesh->ext_face_comm[k].color_out == color_out ) {
          ext_comm_col_out_ptr = &parmesh->ext_face_comm[k];
          break;
        }
      }

      /* Add the new external faces in the external communicator
       * myrank-color_out */
      if ( !ext_comm_col_out_ptr ) {
        /* Creation of a new external communicator */
        PMMG_REALLOC(parmesh,parmesh->ext_face_comm,
                     parmesh->next_face_comm+1,parmesh->next_face_comm,
                     PMMG_Ext_comm,"list of external communicators",goto end);
        parmesh->next_face_comm = parmesh->next_face_comm+1;
        ext_comm = &parmesh->ext_face_comm[parmesh->next_face_comm];
        ext_comm->color_in  = parmesh->myrank;
        ext_comm->color_out = k;
        ext_comm->nitem     = 0;
        ext_comm_col_out_ptr = ext_comm;
      }
      else {
        /* Completion of the existing external communicator */
        ext_comm = ext_comm_col_out_ptr;
      }

      /* New interface faces adding */
      assert(nitem && "nothing to add... weird");
      PMMG_REALLOC(parmesh,ext_comm->int_comm_index,ext_comm->nitem+nitem-2,
                   ext_comm->nitem,int,"int_comm_index",goto end);
      old_nitem = ext_comm->nitem;
      ext_comm->nitem += nitem-2;

      /* Assign a position in the internal comm to the face */
      for ( k=0; k<nitem-2; ++k ) {
        idx_in_face2int = extComm_grpFaces2face2int_recv[k+idx];
        face2int_face_comm_idx2[idx_in_face2int]    = int_comm->nitem;
        ext_comm->int_comm_index[old_nitem+k]       = int_comm->nitem++;
      }
    }
    idx = idx_next;
  }
  assert ( grp_id==(*first_grp_id)+nrecv-1 );

end:
  if ( extComm_next_idx_recv )
    PMMG_DEL_MEM(parmesh,extComm_next_idx_recv,next_comm2recv,int,
                 "extComm_next_idx");
  if ( extComm_grpFaces2face2int_recv )
    PMMG_DEL_MEM(parmesh,extComm_grpFaces2face2int_recv,nitems2recv,int,
                 "extComm_grpFaces2face2int");

  PMMG_Free_MPI_meshDatatype( &mpi_point,&mpi_xpoint,&mpi_tetra,&mpi_xtetra );

  (*first_grp_id) += irecv;

  return irecv;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param part pointer toward the metis partitions.
 * \param max_ngrp maximum number of groups per proc
 *
 * \param nfaces2send \a nfaces2send[k*max_ngrp+id_grp] contains the number of
 * faces of the group \a group_id that belongs to the external communicator
 * myrank-k (so the number of faces of \a grp_id that must be sended to \a k)
 *
 * \param nfaces2recv \a nfaces2recv[k*max_ngrp+id_grp] contains the number of
 * faces that we must receive from the proc \a and for a given group
 *
 * \param next_comm2send number of external communicator to send to each proc
 *
 * \param ext_comms_next_idx[k] position of the external communicators to send
 * to the proc \a k+1 in the \a extComm_next_idx array.
 *
 * \param nitems2send number of items to send to each proc
 *
 * \param items_next_idx positions of the items to send to each proc in the
 * \a extComm_grpFaces* arrays
 *
 * \param extComm_next_idx \a extComm_next_idx[k] position of the external
 * communicator number \a k+1 in the part of the \a extComm_grpFaces* arrays
 * that goes to a given proc (so for each proc, the communicator number 0 is at
 * position 0)
 *
 * \param extComm_grpFaces2face2int for each group, position of the group faces
 * that belongs to the external communicators in the face2int* arrays
 *
 * \param extComm_grpFaces2extComm for each group, position of the group faces
 * that belongs to the external communicators in the external communicators.
 *
 * \return 0 if fail, 1 otherwise
 *
 * For each group, count the number of external faces that the group will send
 * to each processor and fill the \a nfaces2send array.
 * Then fill the array of the list of the external communicators to send to
 * each proc.
 *
 */
static inline
int PMMG_fill_extFaceCommData(PMMG_pParMesh parmesh,idx_t *part,int *max_ngrp,
                              int **nfaces2send,int **nfaces2recv,
                              int **next_comm2send,int **ext_comms_next_idx,
                              int **nitems2send,int **items_next_idx,
                              int **extComm_next_idx,
                              int **extComm_grpFaces2face2int,
                              int **extComm_grpFaces2extComm) {
  PMMG_pInt_comm int_comm;
  PMMG_pExt_comm ext_comm;
  MPI_Comm       comm;
  int            *intvalues,intval,nitem,nprocs;
  int            *face2int_face_comm_idx1,*face2int_face_comm_idx2;
  int            nitem_int_comm,nitem_myrank_dest;
  int            myrank,*nfaces_in_grp;
  int            next_comm,idx,idx_glob,new,idx_per_dest_proc;
  int            ier,proc,k,j,grp_id,dest,next_comm_tot,nitems_tot;

  ier = 0;

  comm   = parmesh->comm;
  nprocs = parmesh->nprocs;
  myrank = parmesh->myrank;

  /** Step 1: Process the external communicators and:
   *
   * - fill \a next_comm2send[k] by the number of external communicators that
   need to be communicated to the processor k (because they contains at least 1
   face of the group to move) and \a ext_comms_next_idx by the last index of
   each list of communicators in the array that list all the communicators;
   *
   * - fill \a nitems2send[k] by the number of integers to send to the proc k
   to communicate the needed data of each external communicators and the
   group ids. Fill \a items_next_idx by the last index of
   each list of items in the array that list all the items;
   *
   * - fill \a nfaces2send with the number of external communicator faces that
   each group send to each proc
   *
   * - fill \a nfaces2recv with the number of external communicator faces that
   each proc receive from each proc
   *
   * - mark the faces that are in the external communicators to be able to
   detect the list of the faces at the interface between two groups of the
   local proc.
  */

  MPI_CHECK( MPI_Allreduce(&parmesh->ngrp,max_ngrp,1,MPI_INT,MPI_MAX,comm),
             goto end);

  PMMG_CALLOC(parmesh,*nfaces2send,(*max_ngrp)*nprocs,int,
              "sended external comm to other procs",goto end);
  PMMG_CALLOC(parmesh,*nfaces2recv,(*max_ngrp)*nprocs,int,
              "received external comm from other procs",goto end);
  PMMG_CALLOC(parmesh,*next_comm2send,nprocs,int,"next_comm2send",goto end);
  PMMG_CALLOC(parmesh,*nitems2send,nprocs,int,"nitems2send",goto end);
  PMMG_CALLOC(parmesh,*ext_comms_next_idx,nprocs,int,"ext_comms_next_idx",goto end);
  PMMG_CALLOC(parmesh,*items_next_idx,nprocs,int,"items_next_idx",goto end);

  PMMG_CALLOC(parmesh,nfaces_in_grp,parmesh->ngrp,int,"nfaces_in_grp",goto end);

  int_comm  = parmesh->int_face_comm;
  nitem     = int_comm->nitem;
  intvalues = int_comm->intvalues;

  for ( grp_id=0; grp_id<parmesh->ngrp; ++grp_id ) {

    dest = part[grp_id];
    if ( dest==myrank ) continue;

    /* intvalues array initialisation with the face index (ensure that intvalues
     * is strictly greater than 0 for a face of this group and strictly lower
     * for the other faces */
    for ( k=0; k<nitem; ++k ) intvalues[k] = PMMG_UNSET;

    face2int_face_comm_idx1 = parmesh->listgrp[grp_id].face2int_face_comm_index1;
    face2int_face_comm_idx2 = parmesh->listgrp[grp_id].face2int_face_comm_index2;
    nitem_int_comm          = parmesh->listgrp[grp_id].nitem_int_face_comm;

    for ( k=0; k<nitem_int_comm; ++k )
      intvalues[face2int_face_comm_idx2[k]] = face2int_face_comm_idx1[k];

    next_comm                  = 0;
    nitems_tot                 = 0;
    for ( k=0; k<parmesh->next_face_comm; ++k ) {
      ext_comm = &parmesh->ext_face_comm[k];
      idx      = (*max_ngrp) * ext_comm->color_out + grp_id;

      for ( j=0; j<ext_comm->nitem; ++j ) {
        if ( intvalues[ext_comm->int_comm_index[j]] > 0 ) {
          intvalues[ext_comm->int_comm_index[j]] *= -1;
          ++(*nfaces2send)[idx];
        }
      }
      if ( (*nfaces2send)[idx] ) {
        ++next_comm;
        nitems_tot += (*nfaces2send)[idx];
      }
    }

    /** Step 3: Process the face2int_face_comm_index2 array and count the number
     * of faces at the interface of 2 local groups */
    for ( k=0; k<nitem_int_comm; ++k ) {
      assert(intvalues[face2int_face_comm_idx2[k]]);
      if ( intvalues[face2int_face_comm_idx2[k]] > 0 ) {
        ++nfaces_in_grp[grp_id];
      }
    }

    if ( nfaces_in_grp[grp_id] ) {
      ++next_comm;
      nitems_tot += nfaces_in_grp[grp_id];
    }

    /* Number of external comm to send to the proc dest */
    (*next_comm2send)[dest] += next_comm;
    /* Number of items (integers) to send to dest (list of external comms
     * (int_comm_index and colour_out + group id) */
    (*nitems2send)[dest]    += nitems_tot+2*next_comm;
  }

  /* Fill ext_comms_next_idx[k] with the position of the first external
   * communicator going to proc k+1 in extComms_next_idx
   *
   * Fill items_next_idx[k] with the position of the first item to communicate to
   * the proc k+1 (position in extComms_grpFaces arrayx) */
  (*ext_comms_next_idx)[0] = (*next_comm2send)[0];
  (*items_next_idx)[0]     =    (*nitems2send)[0];
  for ( k=1; k<nprocs; ++k ) {
    (*ext_comms_next_idx)[k] = (*ext_comms_next_idx)[k-1] +(*next_comm2send)[k];
    (*items_next_idx)[k]     =     (*items_next_idx)[k-1] +   (*nitems2send)[k];
  }

  MPI_CHECK( MPI_Alltoall((*nfaces2send),(*max_ngrp),MPI_INT,
                          (*nfaces2recv),(*max_ngrp),MPI_INT, comm),goto end);

  /** Step 2: Process the groups in growing destination order. For each group:
   *
   * - process the internal faces of the group (if needed) and complete the
   *   external communicator lists by a new communicator containing the list of
   *   the index in \a face2int_face_comm_index1 of the faces that must be added
   *   to the external communicator myrank-dest (and dest-myrank) after the
   *   group displacement;
   *
   * - reset the mark in the intvalues array;
   *
   * - process the external communicators and fill \a extComm_next_idx by the
   *   position of the last item of each external communicator in the \a
   *   extComm_grpFaces* arrays (starting back from the 0 position for each new
   *   proc), fill \a extComm_grpFaces2extComm by the position of the group
   *   faces in the external comm and the id and destination of the group. Last,
   *   fill extComm_grpFaces2face2int by the position of the group faces in the
   *   \a face2int_face_comm arrays, the id of the group and the \a color_out
   *   field of the matching external communicator.
   *
   */
  next_comm_tot = (*ext_comms_next_idx)[nprocs-1];
  nitems_tot    =     (*items_next_idx)[nprocs-1];

  PMMG_MALLOC(parmesh,(*extComm_next_idx),next_comm_tot,int,
              "extComm_next_idx",
              goto end);
  PMMG_MALLOC(parmesh,(*extComm_grpFaces2extComm),nitems_tot,int,
              "extComm_grpFaces2extComm",
              goto end);
  PMMG_MALLOC(parmesh,(*extComm_grpFaces2face2int),nitems_tot,int,
              "extComm_grpFaces2face2int",
              goto end);

  idx      = 0;
  idx_glob = 0;

  for ( proc=0; proc<parmesh->nprocs; ++proc ) {

    idx_per_dest_proc = 0;
    for ( grp_id=0; grp_id<parmesh->ngrp; ++grp_id ) {

      dest = part[grp_id];

      if ( dest != proc || dest==myrank ) continue;

      /* Search the external communicator myrank-dest */
      nitem_myrank_dest = 0;
      for ( k=0; k<parmesh->next_face_comm; ++k ) {
        if ( parmesh->ext_face_comm[k].color_out==dest ) {
          nitem_myrank_dest = parmesh->ext_face_comm[k].nitem;
        }
      }

      face2int_face_comm_idx1 = parmesh->listgrp[grp_id].face2int_face_comm_index1;
      face2int_face_comm_idx2 = parmesh->listgrp[grp_id].face2int_face_comm_index2;
      nitem_int_comm          = parmesh->listgrp[grp_id].nitem_int_face_comm;

      new = 0;
      if ( nfaces_in_grp[grp_id] ) {
        for ( k=0; k<nitem_int_comm; ++k )
          intvalues[face2int_face_comm_idx2[k]] = face2int_face_comm_idx1[k];

        for ( k=0; k<parmesh->next_face_comm; ++k ) {
          ext_comm = &parmesh->ext_face_comm[k];

          for ( j=0; j<ext_comm->nitem; ++j ) {
            if ( intvalues[ext_comm->int_comm_index[j]] > 0 ) {
              intvalues[ext_comm->int_comm_index[j]] *= -1;
            }
          }
        }

        /* Process the internal faces first: store the position of the group
           interfaces in the external comm and in the face2intface_comm_index
           array: It will allows the creation of the new external faces between
           myrank and the destination proc of the group */
        for ( k=0; k<nitem_int_comm; ++k ) {

          assert(intvalues[face2int_face_comm_idx2[k]]);

          if ( intvalues[face2int_face_comm_idx2[k]] > 0 ) {
            (*extComm_grpFaces2face2int)[idx_glob] = k;
            (*extComm_grpFaces2extComm)[idx_glob]  = nitem_myrank_dest+new;
            ++idx_glob;
            ++new;
            ++idx_per_dest_proc;
          }
        }
        (*extComm_grpFaces2face2int)[idx_glob] = grp_id;
        (*extComm_grpFaces2extComm )[idx_glob] = grp_id;
        ++idx_glob;
        ++idx_per_dest_proc;
        (*extComm_grpFaces2face2int)[idx_glob] = myrank;
        (*extComm_grpFaces2extComm )[idx_glob] = myrank;
        (*extComm_next_idx)[idx++] = ++idx_per_dest_proc ;
        ++idx_glob;
      }

      /* Process the external faces: fill intvalues by the position of the face
       * in the face2int_face_comm arrays, then fill \a
       * (*extComm_grpFaces2face2int) by this position */
      for ( k=0; k<nitem_int_comm; ++k )
        intvalues[face2int_face_comm_idx2[k]] = k;

      for ( k=0; k<parmesh->next_face_comm; ++k ) {
        ext_comm = &parmesh->ext_face_comm[k];

        nfaces_in_grp[grp_id] = 0;
        for ( j=0; j<ext_comm->nitem; ++j ) {
          intval = intvalues[ext_comm->int_comm_index[j]];
          if ( intval >= 0 ) {
            (*extComm_grpFaces2face2int)[idx_glob]  = intval;
            (*extComm_grpFaces2extComm )[idx_glob]  = j;
            ++idx_glob;
            ++idx_per_dest_proc;
            ++nfaces_in_grp[grp_id];
          }
        }
        if ( nfaces_in_grp[grp_id] ) {
          (*extComm_grpFaces2face2int)[idx_glob] = grp_id;
          (*extComm_grpFaces2extComm )[idx_glob] = grp_id;
          ++idx_glob;
          ++idx_per_dest_proc;
          /* Store the id of the proc with which we list the interfaces */
          (*extComm_grpFaces2face2int)[idx_glob] = ext_comm->color_out;
          /* Store the id of the proc on which the group move */
          (*extComm_grpFaces2extComm )[idx_glob] = dest;
          (*extComm_next_idx)[idx++] = ++idx_per_dest_proc;
          ++idx_glob;
        }
      }
    }
    assert ( idx_glob == (*items_next_idx)[proc] );
    assert ( idx == (*ext_comms_next_idx)[proc] );
  }

  /* Success */
  ier = 1;

end:
  PMMG_DEL_MEM(parmesh,nfaces_in_grp,parmesh->ngrp,int,"nfaces_in_grp");

  return ier;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param dest rank of destination.
 * \param max_ngrp maximum number of groups per proc
 *
 * \param nfaces2send \a nfaces2send[k*max_ngrp+id_grp] contains the number of
 * faces of the group \a group_id that belongs to the external communicator
 * myrank-k (so the number of faces of \a grp_id that must be sended to \a k)
 *
 * \param next_comm2send number of external communicator to send to each proc
 *
 * \param ext_comms_next_idx[k] position of the external communicators to send
 * to the proc \a k+1 in the \a extComm_next_idx array.
 *
 * \param extComm_next_idx \a extComm_next_idx[k] position of the external
 * communicator number \a k+1 in the part of the \a extComm_grpFaces* arrays
 * that goes to a given proc (so for each proc, the communicator number 0 is at
 * position 0)
 *
 * \param extComm_grpFaces2face2int for each group, position of the group faces
 * that belongs to the external communicators in the face2int* arrays
 *
 * \param extComm_grpFaces2extComm for each group, position of the group faces
 * that belongs to the external communicators in the external communicators.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Send the all the external face communicators that must be updated on the proc
 * \a dest due to the migrations of groups (for example, if the a group move
 * from proc i to proc j and the group has interface faces with the proc \a
 * dest, we must update the external face communicator on \a dest).
 *
 *
 * \warning IMPROVEMENT: change the extComm_grpFaces2extComm storage to have all
 * the external communicators going to dest contiguous so we can send it with
 * only 1 MPI_Send.
 *
 */
static inline
int PMMG_send_extFaceComm(PMMG_pParMesh parmesh,int dest,int max_ngrp,
                          int *nfaces2send,int *next_comm2send,
                          int *ext_comms_next_idx,
                          int *extComm_next_idx,int *extComm_grpFaces2face2int,
                          int *extComm_grpFaces2extComm) {

  PMMG_pExt_comm ext_comm_color_out;
  MPI_Comm       comm;
  int            myrank,tag,grp_id,nprocs,color_out;
  int            k,i,j,idx,idx_next,item_idx,nitem,ier,ext_comm_idx,proc;

  ier = 1;

  comm   = parmesh->comm;
  myrank = parmesh->myrank;
  nprocs = parmesh->nprocs;
  tag    = MPI_SENDEXTFACECOMM_TAG;

  if ( dest == myrank ) return 1;

  /** Process all the external communicators that are listed and send those ones
   * that must go on dest */
  idx_next     = 0;
  item_idx     = 0;
  ext_comm_idx = 0;
  proc         = 0;
  while ( proc<nprocs ) {

    if ( next_comm2send[proc] )  item_idx += idx_next;

    for ( k=0; k<next_comm2send[proc]; ++k ) {
      idx         = idx_next;
      idx_next    = item_idx + extComm_next_idx[ext_comm_idx+k];
      nitem       = idx_next-idx-2;

      grp_id      = extComm_grpFaces2extComm[idx_next-2];
      color_out   = extComm_grpFaces2face2int[idx_next-1];

      if ( color_out != dest ) continue;

      assert ( nitem == nfaces2send [ max_ngrp * dest + grp_id ] );

      MPI_CHECK( MPI_Send(&extComm_grpFaces2extComm[idx],nitem+2,MPI_INT,dest,
                          tag,comm),ier = 0 );
      ++tag;

      /** Remove the old interface faces from the external communicator */
      /* Find the communicator myrank-color_out */
      ext_comm_color_out = NULL;
      for ( i=0; i<parmesh->next_face_comm; ++i ) {
        ext_comm_color_out = &parmesh->ext_face_comm[i];
        if ( ext_comm_color_out->color_out==color_out ) break;
      }
      assert ( ext_comm_color_out );

      for ( i=0; i<nitem; ++i ) {
        j = extComm_grpFaces2extComm[idx+i];
        ext_comm_color_out->int_comm_index[j] = PMMG_UNSET;
      }
    }
    ext_comm_idx = ext_comms_next_idx[proc];
    ++proc;
  }

  return ier;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param source source rank.
 * \param nfaces2recv array that stores the number of communicators to receive
 * from each group and from each proc
 * \param recv_array buffer to receive the external communicators
 *
 * \return 0 if fail, 1 otherwise
 *
 * Send the all the external face communicators that must be updated on the proc
 * \a dest due to the migrations of groups (for example, if the a group move
 * from proc i to proc j and the group has interface faces with the proc \a
 * dest, we must update the external face communicator on \a dest).
 *
 */
static inline
int PMMG_recv_extFaceComm(PMMG_pParMesh parmesh,int source,int max_ngrp,
                          int *nfaces2recv,int *recv_array) {

  PMMG_pExt_comm ext_comm_color_out,ext_comm_source,ext_comm;
  MPI_Comm       comm;
  MPI_Status     status;
  int            myrank,tag,color_out,old_nitem;
  int            k,i,idx,nitem,ier;

  ier = 0;

  comm   = parmesh->comm;
  myrank = parmesh->myrank;
  tag    = MPI_SENDEXTFACECOMM_TAG;

  if ( source == myrank ) return 1;

  /** Process the nfaces2recv array and receive the external communicator of the
   * groups on proc source that send something to us */

  /* Find the communicator myrank-source */
  ext_comm_source = NULL;
  for ( i=0; i<parmesh->next_face_comm; ++i ) {
    ext_comm_source = &parmesh->ext_face_comm[i];
    if ( ext_comm_source->color_out==source ) break;
  }

  for ( k=0; k<max_ngrp; ++k ) {
    if ( !nfaces2recv[ source * max_ngrp + k ] ) continue;

    /** Step 1: reception */
    MPI_Probe(source,tag,parmesh->comm,&status);
    MPI_Get_count(&status,MPI_INT,&nitem);

    assert ( nitem == nfaces2recv[source*max_ngrp+k]+2 );

    MPI_CHECK( MPI_Recv(recv_array,nitem,MPI_INT,source,tag,comm,&status),goto end );
    ++tag;

    /** Step 2: update the matching external communicator: delete the faces from
     * the myrank-source communicator and add it in the myrank-color_out one (with
     * color_out the proc on which move the group that has send us the
     * communicator) */
    color_out = recv_array[nitem-1];

    /* Find the communicator myrank-color_out */
    ext_comm_color_out = NULL;
    for ( i=0; i<parmesh->next_face_comm; ++i ) {
      ext_comm_color_out = &parmesh->ext_face_comm[i];
      if ( ext_comm_color_out->color_out==color_out ) break;
    }

    /* Add the new external faces between source and myrank (the group that
     * was on source is now on myrank) */
    if ( !ext_comm_color_out ) {
      /* Creation of a new external communicator */
      PMMG_REALLOC(parmesh,parmesh->ext_face_comm,
                   parmesh->next_face_comm+1,parmesh->next_face_comm,
                   PMMG_Ext_comm,"list of external communicators",goto end);
      parmesh->next_face_comm = parmesh->next_face_comm+1;

      ext_comm = &parmesh->ext_face_comm[parmesh->next_face_comm];
      ext_comm->color_in  = parmesh->myrank;
      ext_comm->color_out = color_out;
      ext_comm->nitem     = 0;
    }
    else
      ext_comm = ext_comm_color_out;

    PMMG_REALLOC(parmesh,ext_comm->int_comm_index,ext_comm->nitem+nitem-2,
                 ext_comm->nitem,int,"int_comm_index",goto end);
    old_nitem = ext_comm->nitem;
    ext_comm->nitem += nitem-2;

    /* Move the faces from ext_comm_source toward ext_comm */
    assert ( ext_comm_source );
    for ( i=0; i<nitem-2; ++i ) {
      idx = recv_array[i];
      ext_comm->int_comm_index[old_nitem+i] =
        ext_comm_source->int_comm_index[idx];
      ext_comm_source->int_comm_index[idx] = PMMG_UNSET;
    }
  }

  /* Success */
  ier = 1;

end:
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
static inline
int PMMG_mpiexchange_grps(PMMG_pParMesh parmesh,idx_t *part) {
  PMMG_pInt_comm int_comm;
  MPI_Comm       comm;
  int            nitem,myrank,nprocs;
  int            *send_grps,*recv_grps,*nfaces2send,*nfaces2recv;
  int            *next_comm2send,*ext_comms_next_idx,*nitems2send;
  int            *extComm_next_idx,*items_next_idx,*recv_array,nfaces2recv_max;
  int            *extComm_grpFaces2extComm,*extComm_grpFaces2face2int;
  int            nrecv,nsend,max_ngrp,prev_ngrp,new_ngrp,pleft,pright,comm_order;
  int            ier,ier_glob,k,err,part_size;

  myrank    = parmesh->myrank;
  nprocs    = parmesh->nprocs;
  comm      = parmesh->comm;
  err       = 0;
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

  /** Step 1: Merge all the groups that must be send to a given proc into 1
   * group */
  ier = PMMG_merge_grps2send(parmesh,&part);
  part_size = parmesh->ngrp;
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
   * Fill the array of the number of groups that we will send to the other
   * processors (\a send_grps) and receive the array of the number of group that
   * we will receive from the other procs (\a recv_grps);
   *
   */
  ier = 0;

  PMMG_CALLOC(parmesh,send_grps,nprocs,int,"sended groups to other procs",
              goto end);
  PMMG_CALLOC(parmesh,recv_grps,nprocs,int,"received groups from other procs",
              goto end);

  int_comm = parmesh->int_face_comm;
  nitem    = int_comm->nitem;

  PMMG_MALLOC(parmesh,int_comm->intvalues,nitem,int,"intvalues",goto end);

  for ( k=0; k<parmesh->ngrp; ++k ) {
    if ( part[k] == parmesh->myrank ) continue;
    ++send_grps[part[k]];
  }
  MPI_CHECK( MPI_Alltoall(send_grps,1,MPI_INT,recv_grps,1,MPI_INT,comm),
             goto end);

  ier = PMMG_fill_extFaceCommData(parmesh,part,&max_ngrp,
                                  &nfaces2send,&nfaces2recv,
                                  &next_comm2send,&ext_comms_next_idx,
                                  &nitems2send,&items_next_idx,
                                  &extComm_next_idx,&extComm_grpFaces2face2int,
                                  &extComm_grpFaces2extComm);

  if ( !ier ) {
    fprintf(stderr,"\n  ## Error: %s: Unable to compute the sended data of"
            " each group to each proc.",__func__);
    goto end;
  }
  ier = 0;

  /** Step 4: Count the number of groups that we will received and reallocate
   * the listgrp array */
  prev_ngrp = parmesh->ngrp;

  nrecv = 0;
  for ( k=0; k<nprocs; ++k )
    nrecv += recv_grps[k];

  PMMG_RECALLOC(parmesh,parmesh->listgrp,nrecv+parmesh->ngrp,parmesh->ngrp,
                PMMG_Grp,"listgrps",goto end);
  parmesh->ngrp = nrecv+parmesh->ngrp;

  /* Allocation of the array that will receive the external communicators at the
   * largest needed size */
  nfaces2recv_max = 0;
  for ( k=0; k<max_ngrp*nprocs; ++k ) {
    nfaces2recv_max = MG_MAX(nfaces2recv_max,nfaces2recv[k]);
  }
  if ( nfaces2recv_max )
    PMMG_MALLOC(parmesh,recv_array,nfaces2recv_max+2,int,"recv_array",goto end);

  /** Step 5: exchange the groups */
  for ( k=1; k<nprocs; ++k ) {
    /* Index of the right proc to which we send data and of the left proc from
     * which we receive data */
    pleft  = ( myrank+nprocs-k )%nprocs;
    pright = ( myrank+k )%nprocs;

    /* Order of the communication: if 0, we send then receive, if 1, we receive
     * then send */
    comm_order = (myrank/k)%2;

    if ( !comm_order ) {
      nsend = PMMG_send_grp(parmesh,pright,send_grps[pright],part,
                            next_comm2send,ext_comms_next_idx,nitems2send,
                            items_next_idx,extComm_next_idx,
                            extComm_grpFaces2face2int,extComm_grpFaces2extComm);

      nrecv = PMMG_recv_grp(parmesh,&prev_ngrp,pleft,recv_grps[pleft]);

      if ( nsend != send_grps[pright] ) {
        fprintf(stderr,"\n  ## Error: %s: unable to send the %d th mesh from"
                " proc %d to proc %d.\n",__func__,nsend+1,myrank,pright);
        if ( nsend<0 ) {
          ier = -1;
          goto end;
        }
        err = 1;
      }
      if ( nrecv != recv_grps[pleft] ) {
        fprintf(stderr,"\n  ## Error: %s: unable to recv the %d th mesh from"
                " proc %d on proc %d.\n",__func__,nrecv+1,pleft,myrank);
        if ( nrecv<0 ) {
          ier = -1;
          goto end;
        }
        err = 1;
      }

      if ( !PMMG_send_extFaceComm(parmesh,pright,max_ngrp,nfaces2send,
                                  next_comm2send,ext_comms_next_idx,
                                  extComm_next_idx,
                                  extComm_grpFaces2face2int,
                                  extComm_grpFaces2extComm) )
        goto end;

      if ( !PMMG_recv_extFaceComm(parmesh,pleft,max_ngrp,nfaces2recv,
                                  recv_array) )
        goto end;

    }
    else {
      nrecv = PMMG_recv_grp(parmesh,&prev_ngrp,pleft,recv_grps[pleft]);

      nsend = PMMG_send_grp(parmesh,pright,send_grps[pright],part,
                            next_comm2send,ext_comms_next_idx,nitems2send,
                            items_next_idx,extComm_next_idx,
                            extComm_grpFaces2face2int,extComm_grpFaces2extComm);

      if ( nsend != send_grps[pright] ) {
        fprintf(stderr,"\n  ## Error: %s: unable to send the %d th mesh from"
                " proc %d to proc %d.\n",__func__,nsend+1,myrank,pright);
        if ( nsend<0 ) {
          ier = -1;
          goto end;
        }
        err = 1;
      }
      if ( nrecv != recv_grps[pleft] ) {
        fprintf(stderr,"\n  ## Error: %s: unable to recv the %d th mesh from"
                " proc %d on proc %d.\n",__func__,nrecv+1,pleft,myrank);
        if ( nrecv<0 ) {
          ier = -1;
          goto end;
        }
        err = 1;
      }

      if ( !PMMG_recv_extFaceComm(parmesh,pleft,max_ngrp,nfaces2recv,
                                  recv_array) )
        goto end;

      if ( !PMMG_send_extFaceComm(parmesh,pright,max_ngrp,nfaces2send,
                                  next_comm2send,ext_comms_next_idx,
                                  extComm_next_idx,
                                  extComm_grpFaces2face2int,
                                  extComm_grpFaces2extComm) )
        goto end;

    }
  }

  /** Step 6: pack the groups and face communicators */
  new_ngrp = PMMG_pack_grps_norealloc( parmesh,&parmesh->listgrp,parmesh->ngrp );
  PMMG_RECALLOC(parmesh,parmesh->listgrp,new_ngrp,parmesh->ngrp,PMMG_Grp,
                "listgrps",goto end);
  parmesh->ngrp = new_ngrp;

  if ( !parmesh->ngrp ) {
    fprintf(stderr,"  ## Warning: %s: rank %d: processor without any groups.\n",
            __func__,myrank);

  }

  /* Pack the face communicators */
  /* Reallocation because we have added faces in the communicator but not yet
   * deleted it */
  PMMG_REALLOC(parmesh,int_comm->intvalues,int_comm->nitem,nitem,int,
               "intvalues",goto end);

  if ( !PMMG_pack_faceCommunicators(parmesh) ) goto end;

  /** Step 7: Node communicators reconstruction from the face ones */
  if ( !PMMG_build_nodeCommFromFaces(parmesh) ) {
    fprintf(stderr,"\n  ## Unable to build the new node communicators from"
            " the face ones.\n");
    ier = -1;
    goto end;
  }

  /** Step 8: Update the memory repartition */
  if ( !PMMG_parmesh_updateMemMax(parmesh, 105, 1) ) {
    fprintf(stderr,"\n  ## Error: %s: Unable to update the memory repartition"
            " between meshes and communicators.\n",__func__);
    goto end;
  }

#ifndef NDEBUG
  for ( k=0; k<parmesh->ngrp; ++k ) {
    for ( int j=0; j<parmesh->listgrp[k].nitem_int_face_comm; ++j )
      assert( parmesh->listgrp[k].face2int_face_comm_index2[j]>=0 &&
              "face2int_face_comm not filled.");
  }
#endif

  /** Success */
  if ( !err )  ier = 1;

end:
  if ( part )
    PMMG_DEL_MEM(parmesh,part,part_size,idx_t,"deallocate parmetis partition");
  if ( send_grps )
    PMMG_DEL_MEM(parmesh,send_grps,nprocs,int,"send_grps");
  if ( recv_grps )
    PMMG_DEL_MEM(parmesh,recv_grps,nprocs,int,"recv_grps");
  if ( parmesh->int_face_comm->intvalues )
    PMMG_DEL_MEM(parmesh,parmesh->int_face_comm->intvalues,nitem,int,"intvalues");
  if ( recv_array )
    PMMG_DEL_MEM(parmesh,recv_array,nfaces2recv_max+2,int,"recv_array");

  /* Arrays allocated in fill_extFaceComm2send */
  if ( nfaces2send )
    PMMG_DEL_MEM(parmesh,nfaces2send,max_ngrp*nprocs,int,"nfaces2send");
  if ( nfaces2recv )
    PMMG_DEL_MEM(parmesh,nfaces2recv,max_ngrp*nprocs,int,"nfaces2recv");

  if ( extComm_next_idx )
    PMMG_DEL_MEM(parmesh,extComm_next_idx,ext_comms_next_idx[nprocs-1],int,
                 "positions of the external comms in pos_grp_faces* arrays");
  if ( extComm_grpFaces2extComm )
    PMMG_DEL_MEM(parmesh,extComm_grpFaces2extComm,items_next_idx[nprocs-1],
                 int,"extComm_grpFaces2extComm");
  if ( extComm_grpFaces2face2int )
    PMMG_DEL_MEM(parmesh,extComm_grpFaces2face2int,items_next_idx[nprocs-1],
                 int,"extComm_grpFaces2face2int");

  if ( next_comm2send )
    PMMG_DEL_MEM(parmesh,next_comm2send,nprocs,int,"next_comm2send");
  if ( nitems2send )
    PMMG_DEL_MEM(parmesh,nitems2send,nprocs,int,"nitems2send");
  if ( ext_comms_next_idx )
    PMMG_DEL_MEM(parmesh,ext_comms_next_idx,nprocs,int,"ext_comms_next_idx");
  if ( items_next_idx )
    PMMG_DEL_MEM(parmesh,items_next_idx,nprocs,int,"items_next_idx");

  return ier;
}

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

  ier = PMMG_part_parmeshGrps2parmetis(parmesh,part,parmesh->nprocs);
  if ( !ier )
    fprintf(stderr,"\n  ## Unable to compute the new group partition.\n");

  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);

  if ( !ier_glob ) {
    PMMG_DEL_MEM(parmesh,part,parmesh->ngrp,idx_t,"deallocate parmetis partition");
    return ier_glob;
  }

  /** Send the suitable groups to other procs and recieve their groups */
  ier = PMMG_mpiexchange_grps(parmesh,part);
  if ( ier <= 0 )
    fprintf(stderr,"\n  ## Unable to communicate groups through processors.\n");

  return ier;
}

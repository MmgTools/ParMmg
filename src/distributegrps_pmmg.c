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
 * \param grpI group in which we want to merge
 * \param grpJ group that we want to merge with group \a grpI
 *
 * \return 0 if fail, 1 otherwise
 *
 * Merge \a grpJ into grpI.
 *
 */
static inline
int PMMG_merge_grpJinI(PMMG_pGrp grpI, PMMG_pGrp grpJ) {

  return 0;
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
  PMMG_pGrp  grps,grp,listgrp,grpI,grpJ;
  MMG5_pMesh meshI,meshJ;
  long long  memAv;
  idx_t      *tmpPart;
  int        nprocs,ngrp,ngrpOld,ncomm,k,j,ier;

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

  /** Step 1: New groups allocation and initialization: store a group that will
   * be send to proc k in grps[k] */

  // If this step is too expensive in memory, we can count the number of procs
  // with which we will communicate and fill directly the pack array.
  PMMG_CALLOC( parmesh,grps,nprocs,PMMG_Grp,"Groups to send",return 0 );

  j = 0;
  for ( k=0; k<ngrp; ++k ) {
    if ( grps[(*part)[k]].mesh ) continue;

    grps[(*part)[k]] = PMMG_move_grp(&parmesh->listgrp[k]);
    if ( ++j==nprocs ) break;
  }

  /** Step 2: Merge the groups that will be send to the same proc into the
   * same group */
  listgrp = parmesh->listgrp;
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
    if ( !PMMG_merge_grpJinI(grpI,grpJ) ) goto fail;

    memAv -= meshI->memCur;
    meshI->memMax = meshI->memCur;

    /* Delete the useless group to gain memory space */
    memAv += meshJ->memCur;
    PMMG_grp_free(parmesh,&listgrp[k]);
  }
  ngrp = nprocs;

  if ( !PMMG_pack_grpsAndPart( parmesh,&grps,&ngrp,part ) ) {
    ier = -1;
    goto end;
  }

  ier  = 1;
  goto end;

fail:
  /** Error handler: Try to recover all the existing groups to continue the
   * process. */
  ier = 0;
  PMMG_RECALLOC( parmesh,grps,(ngrp-k)+nprocs,nprocs,PMMG_Grp,"Groups to send",
                 ier=-1;goto end);
  ngrp    = (ngrp-k)+nprocs;
  PMMG_CALLOC( parmesh,tmpPart,ngrp,idx_t,"New parmetis partition",
               ier=-1;goto end);

  /* At this step, the groups of listgroups that have been treated, are
   * moved/merged in grps[part[k]]: update the partition array */
  for ( j=0; j<parmesh->ngrp; ++j ) {
    if ( listgrp[j].mesh ) continue;
    assert( grps[(*part)[j]].mesh );
    tmpPart[(*part)[j]] = (*part)[j];
  }
  /* Concatenation of the remaining groups of listgrp into grps */
  j = 0;
  for ( ; k<parmesh->ngrp; ++k ) {
    if ( !listgrp[k].mesh )
      continue;

    /* Find the first available location in grps array */
    while ( grps[j].mesh ) ++j;
    assert ( j<ngrp );

    tmpPart[j]= (*part)[k];
    grps[j++] = PMMG_move_grp(&parmesh->listgrp[k]);
  }
  /* If j<nprocs, it is possible to still have holes in the group array: remove
   * it */
  k = j;
  for ( ; j<nprocs; ++j ) {
    if ( !grps[j].mesh ) continue;
    if ( k != j ) grps[k] = PMMG_assign_grp(&grps[j]);
    tmpPart[k] = j;
    ++k;
  }
  ngrpOld = ngrp;
  ngrp    = k;
  PMMG_REALLOC( parmesh,grps,ngrp,ngrpOld,PMMG_Grp,"Groups to send",
                ier=-1;goto end );
  PMMG_REALLOC( parmesh,tmpPart,ngrp,ngrpOld,idx_t,"parmetis partition",
                ier=-1;goto end );
  PMMG_DEL_MEM( parmesh,*part,nprocs,idx_t,"parmetis partition" );
  (*part) = tmpPart;

end:
  /** Step 3: Update listgrp */
  PMMG_DEL_MEM(parmesh,parmesh->listgrp,parmesh->ngrp,PMMG_Grp,"listgrp");
  parmesh->listgrp = grps;
  parmesh->ngrp    = ngrp;

  if ( !PMMG_parmesh_updateMemMax(parmesh, 105, 1) ) {
    fprintf(stderr,"\n  ## Error: %s: Unable to update the memory repartition"
            " between meshes and communicators.\n",__func__);
    return -1;
  }

  return ier;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param part pointer toward the metis array containing the partitions.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Send each group to the suitable processor.
 *
 */
static inline
int PMMG_send_grps(PMMG_pParMesh parmesh,idx_t *part) {

  return 1;
}

/**
 * \param parmesh pointer toward the mesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Recieve groups from other processors and update \a parmesh
 * (groups and communicators).
 *
 */
static inline
int PMMG_recv_grps(PMMG_pParMesh parmesh) {

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
  int ier,ier_glob;

  /** Merge all the groups that must be send to a given proc into 1 group */
  ier = PMMG_merge_grps2send(parmesh,&part);
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);

  if ( ier_glob < 0 ) {
    fprintf(stderr,"\n  ## Unable to compute the new group partition.\n");
    ier = 0;
    goto fail;
  }
  else if ( !ier_glob ) {
    fprintf(stderr,"\n  ## Unable to compute the new group partition."
            " Try to send it nevertheless\n");
  }

  /** Send each group to the suitable processor */
  if ( !PMMG_send_grps(parmesh,part) ) {
    fprintf(stderr,"\n  ## Unable to compute the new group partition.\n");
    ier = 0;
    goto fail;
  }
  PMMG_DEL_MEM(parmesh,part,parmesh->ngrp,idx_t,"deallocate parmetis partition");

  /** Recieve the groups */
  if ( !PMMG_recv_grps(parmesh) ) {
    fprintf(stderr,"\n  ## Unable to compute the new group partition.\n");
    ier = 0;
    goto fail;
  }

  /** Success */
  ier = 1;

fail:
  if ( part )
    PMMG_DEL_MEM(parmesh,part,parmesh->ngrp,idx_t,"deallocate parmetis partition");

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
int PMMG_distribute_grps( PMMG_pParMesh parmesh )
{
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

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
 * \param parmesh pointer toward the mesh structure.
 * \param part pointer toward the metis array containing the partitions.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Merge all the groups that must be send to a given proc into 1 group.
 * Update the \a part array to match with the new groups.
 *
 */
static inline
int PMMG_merge_grps2send(PMMG_pParMesh parmesh,idx_t *part) {

  return 1;
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
  int ier;

  /** Merge all the groups that must be send to a given proc into 1 group */
  if ( !PMMG_merge_grps2send(parmesh,part) ) {
    fprintf(stderr,"\n  ## Unable to compute the new group partition.\n");
    ier = 0;
    goto fail;
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

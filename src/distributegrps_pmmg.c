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
 * Send the suitable group to other procs and recieve their groups.
 *
 */
static inline
int PMMG_mpiexchange_groups(PMMG_pParMesh parmesh,idx_t *part) {

  return 1;
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
int PMMG_distribute_groups( PMMG_pParMesh parmesh )
{
  idx_t *part;
  int   ier;

  ier = 1;

  /** Get the new partition of groups (1 group = 1 metis node) */
  PMMG_CALLOC(parmesh,part,parmesh->ngrp,idx_t,"allocate parmetis buffer",
              return 0);

  if ( !PMMG_part_parmeshGrps2parmetis(parmesh,part,parmesh->nprocs) ) {
    fprintf(stderr,"\n  ## Unable to compute the new partition.\n"
            " Exit program.\n");
    ier = 0;
  }

  /** Send the suitable groups to other procs and recieve their groups */
  if ( !PMMG_mpiexchange_groups(parmesh,part) ) {
    fprintf(stderr,"\n  ## Unable to communicate groups through processors.\n"
            " Exit program.\n");
    ier = 0;
  }

  PMMG_DEL_MEM(parmesh,part,parmesh->ngrp,idx_t,"deallocate parmetis buffer");

  return ier;
}

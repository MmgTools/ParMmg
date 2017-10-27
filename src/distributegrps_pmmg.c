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
 * Delete parts of the mesh not on the processor.
 */
int PMMG_distribute_groups( PMMG_pParMesh parmesh )
{
  idx_t *part;

  /** Get the new partition of groups (1 group=1 metis node) */
  part = NULL;
  if ( !PMMG_part_parmeshGrps2parmetis(parmesh,part,parmesh->nprocs) ) return 0;

  /** Send the suitable group to other procs and recieve their groups */

  return 1;
}

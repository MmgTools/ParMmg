/**
 * \file loadbalancing_pmmg.c
 * \brief Load balancing after a remeshing step
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \author Nikos Pattakos (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 */
#include "parmmg.h"

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * \return 1 if success, 0 if fail.
 *
 * Count the number of parallel triangle in each tetra and store the info into
 * the mark field of the tetra.
 *
 */
static inline
int PMMG_count_parBdy(PMMG_pParMesh parmesh) {

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * \return 1 if success, 0 if fail.
 *
 * Load balancing of the mesh groups over the processors.
 *
 */
int PMMG_loadBalancing(PMMG_pParMesh parmesh) {

  /** Count the number of interface faces per tetra and store it in mark */
  if ( !PMMG_count_parBdy(parmesh) ) return 0;

  /** Split the ngrp groups of listgrp into a higher number of groups */
  if ( !PMMG_split_n2mGroups(parmesh,METIS_TARGET_MESH_SIZE) ) return 0;

  /** Distribute the groups over the processor to load balance the meshes */
  if ( !PMMG_distribute_groups(parmesh) ) return 0;

  /** Redistribute the ngrp groups of listgrp into a higher number of groups */
  if ( !PMMG_split_n2mGroups(parmesh,REMESHER_TARGET_MESH_SIZE) ) return 0;

  return 1;
}

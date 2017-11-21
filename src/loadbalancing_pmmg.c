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
  MMG5_pMesh   mesh;
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  int          k,i,j;

  for ( i=0; i<parmesh->ngrp; ++i ) {
    mesh = parmesh->listgrp[i].mesh;

    for ( k=1; k<=mesh->ne; ++k ) {
      pt       = &mesh->tetra[k];
      pt->mark = 1;

      if ( (!MG_EOK(pt)) || (!pt->xt) ) continue;
      pxt = &mesh->xtetra[pt->xt];

      for ( j=0; j<4; ++j ) {
        if ( pxt->ftag[j] & MG_PARBDY ) ++pt->mark;
      }
    }
  }

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
  int ier,ier_glob;

  ier = 1;

  /** Count the number of interface faces per tetra and store it in mark */
  if ( !PMMG_count_parBdy(parmesh) ) {
    fprintf(stderr,"\n  ## Problem when counting the number of interface faces."
            " Exit program.\n");
    ier=0;
  }

  /** Split the ngrp groups of listgrp into a higher number of groups */
  if ( !PMMG_split_n2mGroups(parmesh,METIS_TARGET_MESH_SIZE,1) ) {
    fprintf(stderr,"\n  ## Problem when splitting into a higher number of groups."
            " Exit program.\n");
    ier=0;
  }

  /** Distribute the groups over the processor to load balance the meshes */
  if ( !PMMG_distribute_groups(parmesh) ) {
    fprintf(stderr,"\n  ## Group distribution problem. Exit program.\n");
    ier = 0;
  }

  /** Redistribute the ngrp groups of listgrp into a higher number of groups */
  if ( !PMMG_split_n2mGroups(parmesh,REMESHER_TARGET_MESH_SIZE,0) ) {
    fprintf(stderr,"\n  ## Problem when splitting into a lower number of groups."
            " Exit program.\n");
    ier=0;
  }

  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);

  return ier_glob;
}

/* =============================================================================
**  This file is part of the parmmg software package for parallel tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux, 2017-
**
**  parmmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  parmmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with parmmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the parmmg distribution only if you accept them.
** =============================================================================
*/

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
 * \return 1 if success, 0 if fail but we can save the meshes, -1 if we cannot.
 *
 * Load balancing of the mesh groups over the processors.
 *
 */
int PMMG_loadBalancing(PMMG_pParMesh parmesh) {
  MMG5_pMesh mesh;
  int        ier,ier_glob,igrp,ne;
  mytime     ctim[5];
  int8_t     tim;
  char       stim[32];


  tminit(ctim,5);

  /** Count the number of interface faces per tetra and store it in mark,
   * retag old parallel faces*/
  if ( parmesh->info.imprim > PMMG_VERB_DETQUAL ) {
    tim = 0;
    chrono(ON,&(ctim[tim]));
  }

  ier = PMMG_resetOldTag(parmesh);
  if ( !ier ) {
    fprintf(stderr,"\n  ## Problem when counting the number of interface faces.\n");
  }
#ifndef NDEBUG
    /* In debug mode we have mpi comm in split_n2mGrps, thus, if 1 proc fails
     * and the other not we will deadlock */
    MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);
    if ( ier_glob <=0 ) return ier_glob;
#endif
  if ( parmesh->info.imprim > PMMG_VERB_DETQUAL ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"               count para. interfaces    %s\n",stim);
  }

  if ( parmesh->info.imprim > PMMG_VERB_DETQUAL ) {
    tim = 1;
    chrono(ON,&(ctim[tim]));
  }

  ne = 0;
  for ( igrp=0; igrp < parmesh->ngrp; igrp++ )
    ne += parmesh->listgrp[igrp].mesh->ne;

  if ( ier ) {
    /** Split the ngrp groups of listgrp into a higher number of groups */
    ier = PMMG_split_n2mGrps(parmesh,PMMG_GRPSPL_DISTR_TARGET,1);
  }

  /* There is mpi comms in distribute_grps thus we don't want that one proc
   * enters the function and not the other proc */
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);
  if ( parmesh->info.imprim > PMMG_VERB_DETQUAL ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"               group split for metis     %s\n",stim);
  }

  if ( ier_glob <= 0 && parmesh->myrank == parmesh->info.root ) {
    fprintf(stderr,"\n  ## Problem when splitting into a higher number of groups.\n");
    return ier_glob;
  }

  /** Distribute the groups over the processor to load balance the meshes */
  if ( parmesh->info.imprim > PMMG_VERB_DETQUAL ) {
    tim = 2;
    chrono(ON,&(ctim[tim]));
  }

  ier = PMMG_distribute_grps(parmesh);
  if ( ier <= 0 ) {
    fprintf(stderr,"\n  ## Group distribution problem.\n");
  }
#ifndef NDEBUG
  /* In debug mode we have mpi comm in split_n2mGrps, thus, if 1 proc fails
   * and the other not we will deadlock */
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);
  if ( ier_glob <=0 ) return ier_glob;
#endif
  if ( parmesh->info.imprim > PMMG_VERB_DETQUAL ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"               group distribution        %s\n",stim);
  }


  if ( parmesh->info.imprim > PMMG_VERB_DETQUAL ) {
    tim = 3;
    chrono(ON,&(ctim[tim]));
  }

  if ( ier ) {
    /** Redistribute the ngrp groups of listgrp into a higher number of groups */
    ier = PMMG_split_n2mGrps(parmesh,PMMG_GRPSPL_MMG_TARGET,0);
    if ( ier<=0 )
      fprintf(stderr,"\n  ## Problem when splitting into a lower number of groups.\n");
    }

  // Algiane: Optim: is this reduce needed?
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);

  /* Rebuild mesh adjacency for the next adaptation iteration */
  PMMG_TRANSFER_AVMEM_TO_MESHES(parmesh);
  for( igrp = 0; igrp < parmesh->ngrp; igrp++ ) {
    mesh = parmesh->listgrp[igrp].mesh;
    if ( !mesh->adja ) {
      if ( !MMG3D_hashTetra(mesh,0) ) {
        fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
        ier_glob = 0;
      }
    }
  }
  PMMG_TRANSFER_AVMEM_TO_PARMESH(parmesh);

  if ( parmesh->info.imprim > PMMG_VERB_DETQUAL ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"               group split for mmg       %s\n",stim);
  }

  return ier_glob;
}

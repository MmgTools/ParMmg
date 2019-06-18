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
 * Reset the OLDPARBDY tag after mesh adaptation in order to track previous
 * parallel faces during load balancing.
 *
 */
static inline
int PMMG_resetOldTag(PMMG_pParMesh parmesh) {
  MMG5_pMesh   mesh;
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  int          k,i,j;

  for ( i=0; i<parmesh->ngrp; ++i ) {
    mesh = parmesh->listgrp[i].mesh;

    if ( !mesh ) continue;

    for ( k=1; k<=mesh->ne; ++k ) {
      pt       = &mesh->tetra[k];
      pt->mark = 1;

      if ( (!MG_EOK(pt)) || (!pt->xt) ) continue;
      pxt = &mesh->xtetra[pt->xt];

      for ( j=0; j<4; ++j ) {
        if ( pxt->ftag[j] & MG_PARBDY ) {
          /* Increase counter */
          ++pt->mark;
          /* Mark face as a previously parallel one */
          pxt->ftag[j] |= MG_OLDPARBDY;
        } else if ( pxt->ftag[j] & MG_OLDPARBDY ) {
          /* Untag faces which are not parallel anymore */
          pxt->ftag[j] &= ~MG_OLDPARBDY;
        }
      }
    }
  }

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * \return 1 if success, 0 if fail but we can save the meshes, -1 if we cannot.
 *
 * Load balancing of the mesh groups over the processors.
 *
 */
int PMMG_loadBalancing(PMMG_pParMesh parmesh,int moveIfcs) {
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
    ier = PMMG_split_n2mGrps(parmesh,PMMG_GRPSPL_DISTR_TARGET,1,moveIfcs);
  }

  /* There is mpi comms in distribute_grps thus we don't want that one proc
   * enters the function and not the other proc */
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);
  if ( parmesh->info.imprim > PMMG_VERB_DETQUAL ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"               group split for metis     %s\n",stim);
  }

  if ( ier_glob <= 0) {
    fprintf(stderr,"\n  ## Problem when splitting into a higher number of groups.\n");
    return ier_glob;
  }

  /** Distribute the groups over the processor to load balance the meshes */
  if ( parmesh->info.imprim > PMMG_VERB_DETQUAL ) {
    tim = 2;
    chrono(ON,&(ctim[tim]));
  }

  ier = PMMG_distribute_grps(parmesh,moveIfcs);
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
    ier = PMMG_split_n2mGrps(parmesh,PMMG_GRPSPL_MMG_TARGET,0,moveIfcs);
    if ( ier<=0 )
      fprintf(stderr,"\n  ## Problem when splitting into a lower number of groups.\n");
  }

  // Algiane: Optim: is this reduce needed?
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);

  if ( parmesh->info.imprim > PMMG_VERB_DETQUAL ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"               group split for mmg       %s\n",stim);
  }

  return ier_glob;
}

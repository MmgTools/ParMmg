/**
 * \file libparmmg1.c
 * \brief Parallel remeshing library
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version
 * \copyright
 *
 */

#include "parmmg.h"

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 1 if success, 0 otherwise.
 *
 * Check the validity of the input mesh data (tetra orientation, solution
 * compatibility with respect to the provided mesh, Mmg options).
 *
 */
int _PMMG_check_inputData(PMMG_pParMesh parmesh) {
  MMG5_pMesh mesh;
  MMG5_pSol  sol;
  int        k;

  if ( !parmesh->myrank && parmesh->listgrp[0].mesh->info.imprim )
    fprintf(stdout,"\n  -- PMMG: CHECK INPUT DATA\n");

  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    sol  = parmesh->listgrp[k].sol;

    /* Check options */
    if ( mesh->info.lag > -1 ) {
      fprintf(stderr,"  ## Error: lagrangian mode unavailable (MMG3D_IPARAM_lag):\n");
      return 0;
    }
    else if ( mesh->info.iso ) {
      fprintf(stderr,"  ## Error: level-set discretisation unavailable"
              " (MMG3D_IPARAM_iso):\n");
      return 0;
    }
    else if ( mesh->info.optimLES && sol->size==6 ) {
      fprintf(stdout,"  ## Error: strong mesh optimization for LES methods"
              " unavailable (MMG3D_IPARAM_optimLES) with an anisotropic metric.\n");
      return 0;
    }

    /* load data */
    _MMG5_warnOrientation(mesh);

    if ( sol->np && (sol->np != mesh->np) ) {
      fprintf(stdout,"  ## WARNING: WRONG SOLUTION NUMBER. IGNORED\n");
      _MMG5_DEL_MEM(mesh,sol->m,(sol->size*(sol->npmax+1))*sizeof(double));
      sol->np = 0;
    }
    else if ( sol->size!=1 && sol->size!=6 ) {
      fprintf(stderr,"  ## ERROR: WRONG DATA TYPE.\n");
      return 0;
    }
  }
  if ( !parmesh->myrank && parmesh->listgrp[0].mesh->info.imprim )
    fprintf(stdout,"  -- CHECK INPUT DATA COMPLETED\n");

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure where the boundary entities
 * are stored into xtetra and xpoint strucutres
 *
 * Main program of the parallel remeshing library: split the meshes over each
 * proc into groups, then perform niter of sequential remeshing of each group
 * (with moving of the proc boundaries between two iterations) and last, merge
 * the groups over each proc.
 *
 * \return PMMG_STRONGFAILURE if fail and we can't save the mesh (non-conform),
 * PMMG_LOWFAILURE if fail but we can save the mesh, PMMG_SUCCESS if success.
 *
 */
int _PMMG_parmmglib1(PMMG_pParMesh parmesh) {
  MMG5_pMesh       mesh;
  MMG5_pSol        sol;
  int              it,i,niter;

#warning niter must be a param setted by the user
  niter = 1;

  /** Groups creation */
  // if ( !PMMG_splitGrps(parmesh) ) return PMMG_STRONGFAILURE;

  /** Mesh adaptation */
  for ( it=0; it<niter; ++it ) {
    for ( i=0; i<parmesh->ngrp; ++i ) {
      mesh = parmesh->listgrp[i].mesh;
      sol  = parmesh->listgrp[i].sol;

      if ( !MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_nosurf,1 ) )
        return PMMG_STRONGFAILURE;

      if ( !_MMG5_mmg3d1_delone(mesh,sol) ) {
        if ( (!mesh->adja) && !MMG3D_hashTetra(mesh,1) ) {
          fprintf(stderr,"  ## Hashing problem. Invalid mesh.\n");
          return PMMG_STRONGFAILURE;
        }

        if ( (!mesh->adja) && (!MMG3D_hashTetra(mesh,1)) ) {
          fprintf(stderr,"  ## Hashing problem. Invalid mesh.\n");
          return PMMG_STRONGFAILURE;
        }

        if ( !_MMG3D_packMesh(mesh,sol,NULL) )
          return PMMG_STRONGFAILURE;

        if ( !PMMG_mergeGrps(parmesh) ) return PMMG_STRONGFAILURE;

        return PMMG_LOWFAILURE;
      }

#warning Do we need to update the communicators? Does Mmg renum the boundary nodes with -nosurf option?
      /** load Balancing at group scale and communicators reconstruction */

    }
  }

  for ( i=0; i<parmesh->ngrp; ++i ) {
    mesh = parmesh->listgrp[i].mesh;
    sol  = parmesh->listgrp[i].sol;

    /** Merge the groups over each procs */
    if ( !_MMG3D_packMesh(mesh,sol,NULL) )
      return PMMG_STRONGFAILURE;

    _MMG5_DEL_MEM(mesh,mesh->adja,(4*mesh->nemax+5)*sizeof(int));
  }

  if ( !PMMG_mergeGrps(parmesh) ) return PMMG_STRONGFAILURE;

  return(PMMG_SUCCESS);
}

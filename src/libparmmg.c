/**
 * \file libparmmg.c
 * \brief Wrapper for the parallel remeshing library.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version
 * \copyright
 *
 */

#include "parmmg.h"

int PMMG_parmmglib(PMMG_pParMesh parmesh) {
  MMG5_pMesh       mesh;
  MMG5_pSol        sol;
  int              k,ier;

  if ( !parmesh->myrank && parmesh->listgrp[0].mesh->info.imprim ) {
    fprintf(stdout,"  -- PARMMG3d, Release %s (%s) \n",PMMG_VER,PMMG_REL);
    fprintf(stdout,"  -- MMG3d,    Release %s (%s) \n",MG_VER,MG_REL);
    fprintf(stdout,"     %s\n",PMMG_CPY);
    fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);
  }

  /** Check input data */
  if ( !_PMMG_check_inputData(parmesh) ) return PMMG_STRONGFAILURE;

  if ( !parmesh->myrank &&  parmesh->listgrp[0].mesh->info.imprim )
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");


  for ( k=0; k<parmesh->ngrp; ++k ) {

    mesh = parmesh->listgrp[k].mesh;
    sol  = parmesh->listgrp[k].sol;

    /** Function setters (must be assigned before qualisty computation) */
    _MMG3D_Set_commonFunc();
    MMG3D_setfunc(mesh,sol);

    /** Mesh scaling and qualisty histogram */
    if ( !_MMG5_scaleMesh(mesh,sol) ) return PMMG_STRONGFAILURE;

    if ( !_MMG3D_tetraQual(mesh,sol) ) return PMMG_STRONGFAILURE;

#warning todo: send the data to proc 0 and depending on the imprim value, print either the qualisty of each mesh or the reduced qualisty.
    if ( abs(mesh->info.imprim) > 0 ) {
      if ( !_MMG3D_inqua(mesh,sol) ) {
        if ( !_MMG5_unscaleMesh(mesh,sol) ) return PMMG_STRONGFAILURE;
        return PMMG_LOWFAILURE;
      }
    }

    /* specific meshing */
    if ( mesh->info.optim && !sol->np ) {
      if ( !MMG3D_doSol(mesh,sol) ) {
        if ( !_MMG5_unscaleMesh(mesh,sol) )  return PMMG_STRONGFAILURE;
        return PMMG_LOWFAILURE;
      }
      _MMG3D_scalarSolTruncature(mesh,sol);
    }

    /** Mesh analysis */
    if ( !MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_nosurf,1 ) )
      return PMMG_STRONGFAILURE;

    if ( !_MMG3D_analys(mesh) ) {
      if ( !_MMG5_unscaleMesh(mesh,sol) ) return PMMG_STRONGFAILURE;
      return PMMG_LOWFAILURE;
    }

#warning todo: send the data over the proc 0 and print the lengths of each mesh or the global lengths (reduced)
    if ( mesh->info.imprim > 1 && sol->m ) _MMG3D_prilen(mesh,sol,0);

  }

  if ( !parmesh->myrank &&  parmesh->listgrp[0].mesh->info.imprim )
    fprintf(stdout,"\n  -- PHASE 1 COMPLETED\n");

  /** Remeshing */
  if ( !parmesh->myrank &&  parmesh->listgrp[0].mesh->info.imprim )
    fprintf(stdout,"\n  -- PHASE 2 : %s MESHING\n",
            parmesh->listgrp[0].sol->size < 6 ? "ISOTROPIC" : "ANISOTROPIC");

  ier = _PMMG_parmmglib1(parmesh);
  fprintf(stdout,"  -- PHASE 2 COMPLETED.\n");

  if ( ier == MMG5_STRONGFAILURE ) return PMMG_STRONGFAILURE;

  /** Boundaries reconstruction */
  if ( !parmesh->myrank &&  parmesh->listgrp[0].mesh->info.imprim )
    fprintf(stdout,"\n   -- PHASE 3 : MESH PACKED UP\n");

  if ( !MMG3D_hashTetra(parmesh->listgrp[0].mesh,0) ) return PMMG_STRONGFAILURE;
  if ( _MMG3D_bdryBuild(parmesh->listgrp[0].mesh)<0 ) {
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      sol  = parmesh->listgrp[k].sol;
      if ( !_MMG5_unscaleMesh(mesh,sol) )  return PMMG_STRONGFAILURE;
    }
    return PMMG_LOWFAILURE;
  }

  if ( !parmesh->myrank &&  parmesh->listgrp[0].mesh->info.imprim )
    fprintf(stdout,"\n   -- PHASE 3 COMPLETED.\n");

  /** Mesh unscaling */
  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    sol  = parmesh->listgrp[k].sol;
    if ( !_MMG5_unscaleMesh(mesh,sol) )  return PMMG_STRONGFAILURE;
  }

  return(ier);
}

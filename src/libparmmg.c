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

/**
 * \param  parmesh pointer toward a parmesh structure
 *
 * \return PMMG_STRONFAILURE error and a conformant mesh can not be saved
 *         PMMG_LOWFAILURE   error but a conformant mesh can be saved
 *         PMMG_SUCCESS      on success
 *
 *  main parmmh library call
 */
int PMMG_parmmglib(PMMG_pParMesh parmesh) {
  MMG5_pMesh       mesh;
  MMG5_pSol        met;
  int              k,ier;
  mytime     ctim[TIMEMAX];
  char       stim[32];


  if ( parmesh->info.imprim ) {
    fprintf(stdout,"  -- PARMMG3d, Release %s (%s) \n",PMMG_VER,PMMG_REL);
    fprintf(stdout,"     %s\n",PMMG_CPY);
    fprintf(stdout,"     %s %s\n\n",__DATE__,__TIME__);

    fprintf(stdout,"  -- MMG3d,    Release %s (%s) \n",MG_VER,MG_REL);
    fprintf(stdout,"     %s\n",MG_CPY);
  }

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /** Check input data */
  if ( !PMMG_check_inputData(parmesh) )
    return PMMG_STRONGFAILURE;

  chrono(ON,&(ctim[1]));
  if ( parmesh->info.imprim ) {
    fprintf(stdout,"\n  %s\n   MODULE PARMMG3D: IMB-LJLL : %s (%s)\n  %s\n",
            PMMG_STR,PMMG_VER,PMMG_REL,PMMG_STR);
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  }


  for ( k=0; k<parmesh->ngrp; ++k ) {

    mesh = parmesh->listgrp[k].mesh;
    met  = parmesh->listgrp[k].met;

    /** Function setters (must be assigned before quality computation) */
    _MMG3D_Set_commonFunc();

    /** Mesh scaling and quality histogram */
    if ( !_MMG5_scaleMesh(mesh,met) )
      return PMMG_STRONGFAILURE;

    /* specific meshing */
    if ( mesh->info.optim && !met->np ) {
      if ( !MMG3D_doSol(mesh,met) ) {
        if ( !_MMG5_unscaleMesh(mesh,met) )
          return PMMG_STRONGFAILURE;
        return PMMG_LOWFAILURE;
      }
      _MMG3D_solTruncatureForOptim(mesh,met);
    }
    if ( mesh->info.hsiz > 0. ) {
      if ( !MMG3D_Set_constantSize(mesh,met) ) {
        if ( !_MMG5_unscaleMesh(mesh,met) )
          return PMMG_STRONGFAILURE;
        return PMMG_LOWFAILURE;
      }
    }

    MMG3D_setfunc(mesh,met);
    if ( !_MMG3D_tetraQual( mesh,met, 0 ) )
      return PMMG_STRONGFAILURE;

    if ( abs(mesh->info.imprim) > 0 ) {
      if ( !_MMG3D_inqua(mesh,met) ) {
        if ( !_MMG5_unscaleMesh(mesh,met) )
          return PMMG_STRONGFAILURE;
        return PMMG_LOWFAILURE;
      }
    }

    /** Mesh analysis */
    if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_nosurf,1 ) )
      return PMMG_STRONGFAILURE;

    if ( !_MMG3D_analys(mesh) ) {
      if ( !_MMG5_unscaleMesh(mesh,met) ) return PMMG_STRONGFAILURE;
      return PMMG_LOWFAILURE;
    }

#warning todo: send the data over the proc 0 and print the lengths of each mesh or the global lengths (reduced)
    if ( mesh->info.imprim > 1 && met->m )
      _MMG3D_prilen(mesh,met,0);
  }

  chrono(OFF,&(ctim[1]));
  if ( parmesh->info.imprim )
    fprintf(stdout,"\n  -- PHASE 1 COMPLETED.     %s\n",stim);

  /** Remeshing */
  chrono(ON,&(ctim[2]));
  if ( parmesh->info.imprim )
    fprintf(stdout,"\n  -- PHASE 2 : %s MESHING\n",
            parmesh->listgrp[0].met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC");

  ier = PMMG_parmmglib1(parmesh);

  chrono(OFF,&(ctim[2]));
  printim(ctim[2].gdif,stim);
  if ( parmesh->info.imprim ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE PARMMG3d: IMB-LJLL \n  %s\n",
            PMMG_STR,PMMG_STR);
  }


  if ( ier == MMG5_STRONGFAILURE )
    return PMMG_STRONGFAILURE;

  /** Boundaries reconstruction */
  chrono(ON,&(ctim[1]));
  if ( parmesh->info.imprim )
    fprintf(stdout,"\n   -- PHASE 3 : MESH PACKED UP\n");

  if ( MMG3D_bdryBuild(parmesh->listgrp[0].mesh)<0 ) {
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      met  = parmesh->listgrp[k].met;
      if ( !_MMG5_unscaleMesh(mesh,met) )
        return PMMG_STRONGFAILURE;
    }
    return PMMG_LOWFAILURE;
  }

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);

  if ( parmesh->info.imprim )
    fprintf(stdout,"\n   -- PHASE 3 COMPLETED.     %s\n",stim);

  /** Mesh unscaling */
  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    met  = parmesh->listgrp[k].met;
    if ( !_MMG5_unscaleMesh(mesh,met) ) {
      ier = MG_MAX(PMMG_STRONGFAILURE,ier);
    }
  }

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( parmesh->info.imprim )
    fprintf(stdout,"\n   PARMMG3DLIB: ELAPSED TIME  %s\n",stim);

  return(ier);
}

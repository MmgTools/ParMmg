/**
 * \file inout_pmmg.c
 * \brief io for the parmmg software
 * \author Algiane Froehly (InriaSoft)
 * \version 1
 * \date 07 2018
 * \copyright GNU Lesser General Public License.
 *
 * input/outputs for parmmg.
 *
 */

#include "parmmg.h"

int PMMG_loadMesh_centralized(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }
  mesh = parmesh->listgrp[0].mesh;

  return MMG3D_loadMesh(mesh,filename);
}

int PMMG_loadMet_centralized(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;
  MMG5_pSol  met;

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }
  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  return MMG3D_loadSol(mesh,met,filename);
}

int PMMG_loadAllSols_centralized(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;
  MMG5_pSol  sol;

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }
  mesh = parmesh->listgrp[0].mesh;
  sol  = parmesh->listgrp[0].sol;

  return MMG3D_loadAllSols(mesh,&sol,filename);
}

int PMMG_saveMesh_centralized(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }
  mesh = parmesh->listgrp[0].mesh;

  return MMG3D_saveMesh(mesh,filename);
}

int PMMG_saveMet_centralized(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;
  MMG5_pSol  met;

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }
  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  return MMG3D_saveSol(mesh,met,filename);
}

int PMMG_saveAllSols_centralized(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;
  MMG5_pSol  sol;

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }
  mesh = parmesh->listgrp[0].mesh;
  sol  = parmesh->listgrp[0].sol;

  return MMG3D_saveAllSols(mesh,&sol,filename);
}

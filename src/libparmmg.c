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
  int              ier;

  /** Mesh analysis: to implement */
  // if ( !PMMG_analys(parmesh) ) return PMMG_STRONGFAILURE;

  ier = _PMMG_parmmglib1(parmesh);
  if ( ier == MMG5_STRONGFAILURE ) return PMMG_STRONGFAILURE;

  /** Boundaries reconstruction */
  if ( !MMG3D_hashTetra(parmesh->listgrp[0].mesh,0) ) return PMMG_STRONGFAILURE;
  if ( _MMG3D_bdryBuild(parmesh->listgrp[0].mesh)<0 ) return PMMG_LOWFAILURE;

  return(ier);
}

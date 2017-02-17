/**
 * \file parmmg.h
 * \brief internal functions headers for parmmg
 * \author
 * \version
 * \date 11 2016
 * \copyright
 */

#ifndef _PARMMG_H
#define _PARMMG_H

#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>
#include <float.h>
#include <math.h>

#include "libparmmg.h"

#ifdef __cplusplus
extern "C" {
#endif

  /** Free allocated pointers of parmesh, call MPI_Finalize and return value val */
#define _PMMG_RETURN_AND_FREE(parmesh,val)do                  \
  {                                                           \
  int k;                                                      \
                                                              \
  for ( k=0; k<parmesh->ngrp; ++k ) {                         \
    grp  = &parmesh->listgrp[k];                              \
    mesh = grp->mesh;                                         \
    sol  = grp->sol;                                          \
                                                              \
    MMG3D_Free_all(MMG5_ARG_start,                            \
                   MMG5_ARG_ppMesh,&parmesh->listgrp[k].mesh,  \
                   MMG5_ARG_ppMet,&parmesh->listgrp[k].sol,    \
                   MMG5_ARG_end);                             \
  }                                                           \
                                                              \
  _MMG5_SAFE_FREE(parmesh);                                   \
                                                              \
  MPI_Finalize();                                             \
                                                              \
  return(val);                                                \
                                                              \
  }while(0)


  int _PMMG_Init_parMesh_var( va_list argptr );
  int _PMMG_check_inputData(PMMG_pParMesh parmesh);

  int _PMMG_parmmglib1(PMMG_pParMesh parmesh);

  int  PMMG_mergeGrps(PMMG_pParMesh parmesh);
  int  PMMG_bdryStore( MMG5_pMesh mesh);
  int  PMMG_bcastMesh(PMMG_pParMesh parmesh);

#ifdef __cplusplus
}
#endif

#endif

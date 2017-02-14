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

  int _PMMG_Init_parMesh_var( va_list argptr );
  int PMMG_mergeGrps(PMMG_pParMesh parmesh);
  int PMMG_bdryStore( MMG5_pMesh mesh);

  int _PMMG_parmmglib1(PMMG_pParMesh parmesh);

#ifdef __cplusplus
}
#endif

#endif

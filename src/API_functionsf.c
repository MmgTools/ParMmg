/**
 * \file API_functionsf.c
 * \brief Fortran API functions for PARMMG library.
 *
 * Define the Fortran API functions for PARMMG library: adds function
 * definitions with upcase, underscore and double underscore to match
 * any fortran compiler.
 *
 */
#include "libparmmg.h"
#include "mmgcommon.h"

FORTRAN_NAME(PMMG_SET_GRPSIZE, pmmg_set_grpsize,
    (PMMG_pGrp *grp, int* np, int* ne, int* nprism, int* nt,
     int* nquad, int* na, int* typEntity, int* typSol, int* typMet,
     int* retval),
    (grp, np, ne, nprism, nt, nquad, na, typEntity, typSol, typMet,
     retval)) {
  *retval = PMMG_Set_grpSize(*grp, *np, *ne, *nprism, *nt, *nquad,
                             *na, *typEntity, *typSol, *typMet);
  return;
}

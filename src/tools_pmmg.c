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
 * \file tools_pmmg.c
 * \brief Various tools for parMmg
 * \author Algiane Froehly (InriaSoft)
 * \version 1
 * \date 07 2018
 * \copyright GNU Lesser General Public License.
 *
 * various tools for parmmg.
 *
 */

#include "parmmg.h"

/**
 * \param typArg integer defined by the PMMG_ARG_* preprocessor commands
 *
 * \return the name of the \a PMMG_ARG_* variable
 *
 * Print the \a PMMG_ARG_* name associated to the \a typArg value.
 *
 *
 */
const char* PMMG_Get_pmmgArgName(int typArg)
{

  switch ( typArg )
  {
  case(PMMG_ARG_start):
    return "PMMG_ARG_start";
    break;

  case(PMMG_ARG_ppParMesh):
    return "PMMG_ARG_ppParMesh";
    break;

  case(PMMG_ARG_pMesh):
    return "PMMG_ARG_pMesh";
    break;

  case(PMMG_ARG_pMet):
    return "PMMG_ARG_pMet";
    break;

  case(PMMG_ARG_pSols):
    return "PMMG_ARG_pSols";
    break;

  case(PMMG_ARG_pDisp):
    return "PMMG_ARG_pDisp";
    break;

  case(PMMG_ARG_pLs):
    return "PMMG_ARG_pLs";
    break;

  case(PMMG_ARG_ngroups):
    return "PMMG_ARG_ngroups";
    break;

  case(PMMG_ARG_dim):
    return "PMMG_ARG_dim";
    break;

  case(PMMG_ARG_MPIComm):
    return "PMMG_ARG_MPIComm";
    break;

  case(PMMG_ARG_end):
    return "PMMG_ARG_end";
    break;

  default:
    return "PMMG_ARG_Unknown";
  }
}

/**
 * \param info pointer toward a MMG5 info structure that we want to copy
 * \param info_cpy pointer toward a MMG5 info structure into copy the data
 *
 * \return 1 if success, 0 if not
 *
 * Copy the info data into info_cpy, allocate the \a mat and \a par structures
 * if needed.
 *
 */
int PMMG_copy_mmgInfo ( MMG5_Info *info, MMG5_Info *info_cpy ) {
  MMG5_pMat mat_tmp;
  MMG5_pPar par_tmp;
  int       *lookup_tmp;
  int i;

  // assert to remove (we may authorize to have mat and par already allocated )
  assert ( (!info_cpy->mat) && (!info_cpy->par) );

  if ( info->nmat && (!info_cpy->mat) ) {
    MMG5_SAFE_CALLOC(mat_tmp,info->nmat,MMG5_Mat,return 0);
  }
  else {
    mat_tmp = info_cpy->mat;
  }
  if ( mat_tmp ) {
    *mat_tmp = *info->mat;
    for ( i=0; i<info->nmat; ++i ) {
      mat_tmp[i]=info->mat[i];
    }
  }

  if ( info->nmat && (!info_cpy->invmat.lookup) ) {
    MMG5_SAFE_CALLOC(lookup_tmp,info->invmat.size,int,return 0);
  }
  else {
    lookup_tmp = info_cpy->invmat.lookup;
  }
  if ( lookup_tmp ) {
    *lookup_tmp = *info->invmat.lookup;
    for ( i=0; i<info->invmat.size; ++i ) {
      lookup_tmp[i]=info->invmat.lookup[i];
    }
  }

  /* local parameters */
  if ( info->npar && !info_cpy->par ) {
    MMG5_SAFE_CALLOC(par_tmp,info->npar,MMG5_Par,return 0);
  }
  else {
    par_tmp = info_cpy->par;
  }
  if ( par_tmp ) {
    *par_tmp = *info->par;
    for ( i=0; i<info->npar; ++i ) {
      par_tmp[i]=info->par[i];
    }
  }

  *info_cpy = *info;

  info_cpy->mat = mat_tmp;
  info_cpy->par = par_tmp;
  info_cpy->invmat.lookup = lookup_tmp;

  return 1;
}

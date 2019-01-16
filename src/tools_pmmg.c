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
  }

  *info_cpy = *info;

  info_cpy->mat = mat_tmp;
  info_cpy->par = par_tmp;

  return 1;
}

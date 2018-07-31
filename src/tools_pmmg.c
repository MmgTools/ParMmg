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

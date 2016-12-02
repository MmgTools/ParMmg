/**
 * \file inout_3d.c
 * \brief Input / Output Functions.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */


#include "libparmmg.h"

/**
 * \param parmesh pointer toward the ParMesh structure.
 * \param filename name of file.
 * \return 0 if the file is not found, -1 if we detect mismatch parameters,
 * 1 otherwise.
 *
 * Read mesh data.
 *
 */
int PMMG_loadMesh(PMMG_pParMesh parmesh, const char *filename) {
  PMMG_pGrp  grp;

  /*listgrp allocation*/
  parmesh->ngrp = 1;
  _MMG5_SAFE_CALLOC(parmesh->listgrp,parmesh->ngrp,PMMG_Grp);

  grp = &parmesh->listgrp[0];

  grp->mesh = NULL;
  grp->sol  = NULL;

  MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&grp->mesh,
                  MMG5_ARG_ppMet,&grp->sol,
                  MMG5_ARG_end);
  if ( !MMG3D_Set_inputMeshName(grp->mesh,filename) ) return(0);

  if ( !parmesh->myrank )
    if ( !MMG3D_Set_iparameter(grp->mesh,grp->sol,MMG3D_IPARAM_verbose,0) )
      return(0);

  if ( !MMG3D_loadMesh(grp->mesh,filename) )  return(0);


  return(1);
}

/**
 * \param parmesh pointer toward the ParMesh structure.
 * \param filename pointer toward the name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Save mesh data.
 *
 *
 */
int PMMG_saveMesh(PMMG_pParMesh parmesh, const char *filename) {

  //_MMG3D_packMesh(parmesh->listgrp[0].mesh);
  if ( !MMG3D_saveMesh(parmesh->listgrp[0].mesh,filename) )  return(0);

  return(1);
}

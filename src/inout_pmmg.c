/**
 * \file inout_pmmg.c
 * \brief Input / Output Functions.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version
 * \copyright
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

  if ( parmesh->ngrp > 1 ) {
    printf("  ## Error: Can not call the PMMG_loadMesh function with more"
           " than 1 group per processor.\n Exit Programm.\n");
    return 0;
  }
  else if ( !parmesh->ngrp ) return 1;

  grp = &parmesh->listgrp[0];

  if ( !MMG3D_Set_inputMeshName(grp->mesh,filename) ) return(0);

  if ( !parmesh->myrank )
    if ( !MMG3D_Set_iparameter(grp->mesh,grp->sol,MMG3D_IPARAM_verbose,0) )
      return(0);

  if ( !MMG3D_loadMesh(grp->mesh,filename) )  return(0);


  return(1);
}

/**
 * \param parmesh pointer toward the ParMesh structure.
 * \param filename name of file.
 * \return 0 if the file is not found, -1 if we detect mismatch parameters,
 * 1 otherwise.
 *
 * Read mesh data.
 *
 */
int PMMG_loadSol(PMMG_pParMesh parmesh, const char *filename) {
  PMMG_pGrp  grp;

  if ( parmesh->ngrp > 1 ) {
    printf("  ## Error: Can not call the PMMG_loadSol function with more"
           " than 1 group per processor.\n Exit Programm.\n");
    return 0;
  }
  else if ( !parmesh->ngrp ) return 1;

  if ( !MMG3D_Set_inputSolName(grp->mesh,grp->sol,filename) ) return(0);

  if ( !MMG3D_loadSol(grp->mesh,grp->sol,filename) )  return(0);


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
  PMMG_pGrp  grp;
  MMG5_pMesh mesh;
  MMG5_pSol  sol;

  if ( parmesh->ngrp > 1 ) {
    printf("  ## Error: Groups must have been merged (PMMG_mergeGrps function)"
           " before calling the PMMG_saveMesh function.\n Exit Programm.\n");
    return 0;
  }
  else if ( !parmesh->ngrp ) return 1;

  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;
  sol  = grp->sol;


#warning I am not sur that we must hash the tetra and pack the mesh here...
  /** Tetra adjacency reconstruction */
  _MMG5_SAFE_FREE(mesh->adja);
  if ( !MMG3D_hashTetra(mesh,0) ) {
    fprintf(stderr,"  ## PMMG Hashing problem (1). Exit program.\n");
    return(0);
  }

  /* Pack the mesh */
  if ( !_MMG3D_packMesh(mesh,sol,0) ) {
    fprintf(stderr,"  ## PMMG Packing problem (1). Exit program.\n");
    return(0);
  }

  if ( !MMG3D_saveMesh(mesh,filename) )  return(0);

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
int PMMG_saveSol(PMMG_pParMesh parmesh, const char *filename) {

  if ( parmesh->ngrp > 1 ) {
    printf("  ## Error: Groups must have been merged (PMMG_mergeGrps function)"
           " before calling the PMMG_saveSol function.\n Exit Programm.\n");
    return 0;
  }
  else if ( !parmesh->ngrp ) return 1;

  return MMG3D_saveSol(parmesh->listgrp[0].mesh,
                       parmesh->listgrp[0].sol,
                       filename);
}

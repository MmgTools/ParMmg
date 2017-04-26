/**
 * \file inout_pmmg.c
 * \brief Input / Output Functions.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version
 * \copyright
 */


#include "parmmg.h"

/**
 * \param parmesh pointer to working ParMesh structure.
 * \return 0 file not found
 *        -1 parameters mismatch detected
 *         1 success
 *
 * Read mesh data from file in to the parmesh structure.
 *
 */
int PMMG_loadMesh( PMMG_pParMesh parmesh )
{
  assert ( parmesh && "please provide a non empty parmesh" );
  assert ( ( parmesh->listgrp ) && ( parmesh->ngrp != 0 ) && "please provide a none empty mesh" );
  PMMG_pGrp grp = &parmesh->listgrp[0];

  if ( parmesh->ngrp != 1 )
    printf("  ## WARNING: PMMG_loadMesh will only load the file specified in "
           "the FIRST mesh of the processor's group of meshes.\n");

  return ( MMG3D_loadMesh( grp->mesh, grp->mesh->namein ) );
}

/**
 * \param parmesh pointer to working ParMesh structure.
 * \param filename name of sol file to read from.
 * \return  0 file not found
 *         -1 parameters mismatch detected
 *          1 success
 *
 * Read Sol data from file in to the parmesh structure.
 *
 */
int PMMG_loadSol( PMMG_pParMesh parmesh )
{
  assert ( parmesh && "please provide a non empty parmesh" );
  assert ( ( parmesh->listgrp ) && ( parmesh->ngrp != 0 ) && "please provide a none empty mesh" );
  PMMG_pGrp grp = &parmesh->listgrp[0];

  if ( parmesh->ngrp != 1 )
    printf("  ## WARNING: PMMG_loadSol will only load the file specified in "
           "the FIRST mesh of the processor's group of meshes.\n");

  return ( MMG3D_loadSol( grp->mesh, grp->met, grp->met->namein ) );
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
//  MMG5_pSol  met;

  if ( parmesh->ngrp > 1 ) {
    printf("  ## Error: Groups must have been merged (PMMG_mergeGrps function)"
           " before calling the PMMG_saveMesh function.\n Exit Programm.\n");
    return 0;
  }
  else if ( !parmesh->ngrp ) return 1;

  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;
//  met  = grp->met;

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
                       parmesh->listgrp[0].met,
                       filename);
}

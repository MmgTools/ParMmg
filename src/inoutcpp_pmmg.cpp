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
 * \brief io functions that need cpp features.
 * \author Algiane Froehly (InriaSoft)
 * \version 1
 * \date 07 2018
 * \copyright GNU Lesser General Public License.
 *
 * input/outputs for parmmg.
 *
 */

#ifdef USE_VTK
#include "vtkparser.hpp"
#include <vtkMultiProcessController.h>
#include <vtkMPIController.h>
#endif

#include "parmmg.h"

int PMMG_loadVtuMesh_centralized(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;
  MMG5_pSol  met, sol;
  int        ier;

  if ( parmesh->myrank!=parmesh->info.root ) {
    return 1;
  }

#ifndef USE_VTK

  fprintf(stderr,"  ** VTK library not founded. Unavailable file format.\n");
  return -1;

#else

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }
  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  ier = MMG3D_loadVtuMesh(mesh,sol,met,filename);

  /* Restore the mmg verbosity to its initial value */
  mesh->info.imprim = parmesh->info.mmg_imprim;

  if ( 1 != ier ) return 0;

  return 1;
#endif
}


int PMMG_savePvtuMesh(PMMG_pParMesh parmesh, const char * filename) {

#ifndef USE_VTK
  if ( parmesh->myrank == parmesh->info.root ) {
    fprintf(stderr,"  ** VTK library not founded. Unavailable file format.\n");
  }
  return -1;

#else
  char* mdata=NULL; // master file name
  MMG5_SAFE_CALLOC(mdata,strlen(filename)+6,char,return 0);

  strcpy(mdata,filename);
  char *ptr = MMG5_Get_filenameExt(mdata);
  *ptr = '\0'; // get basename
  sprintf( mdata, "%s.pvtu",mdata);

  MMG5_pMesh mesh = parmesh->listgrp[0].mesh;
  MMG5_pSol  met  = parmesh->listgrp[0].met;

  // Warning : does it works with a communicator /= to MPI_COMM_WORLD?
  vtkMPIController *vtkController = vtkMPIController::New();
  vtkController->Initialize();
  vtkMultiProcessController::SetGlobalController(vtkController);

  return MMG5_saveVtkMesh_i<vtkUnstructuredGrid,vtkXMLUnstructuredGridWriter,vtkXMLPUnstructuredGridWriter>
    (mesh,&met,mdata,1,1,parmesh->nprocs,parmesh->myrank,parmesh->info.root);

  MMG5_SAFE_FREE(mdata);

#endif
  return 1;
}

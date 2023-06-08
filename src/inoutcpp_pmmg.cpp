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

  fprintf(stderr,"  ** VTK library not found. Unavailable file format.\n");
  return -1;

#else

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in your parmesh.",
            __func__);
    return 0;
  }
  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  ier = MMG3D_loadVtuMesh(mesh,met,sol,filename);

  /* Restore the mmg verbosity to its initial value */
  mesh->info.imprim = parmesh->info.mmg_imprim;

  if ( 1 != ier ) return 0;

  return 1;
#endif
}


int PMMG_savePvtuMesh(PMMG_pParMesh parmesh, const char * filename) {

#ifndef USE_VTK
  if ( parmesh->myrank == parmesh->info.root ) {
    fprintf(stderr,"  ** VTK library not found. Unavailable file format.\n");
  }
  return -1;

#else
  char* mdata=NULL; // master    file name
  char* sdata=NULL; // secondary file name
  int i;
  MMG5_SAFE_CALLOC(mdata,strlen(filename)+6,char,return 0);
  MMG5_SAFE_CALLOC(sdata,strlen(filename)+6,char,return 0);

  strcpy(mdata,filename);
  char *ptr = MMG5_Get_filenameExt(mdata);
  *ptr = '\0'; // get basename

  // If the output *.pvtu filename has dots "." in the basename,
  // replace them with dashes "-".
  // Why? In VTK function SetFileName(filename), the first dot "." in the
  // filename is interpreted as the extension start. So, whatever the
  // user specifies after the first dot "." will be ignored by VTK. To overcome
  // this, dots are replaced by dashes.
  for(i=0;mdata[i]!='\0';i++) {
    if(mdata[i]=='.') {
      mdata[i] = '-';
    }
  }

  snprintf( sdata,strlen(mdata)+6, "%s.pvtu",mdata);

  MMG5_pMesh mesh = parmesh->listgrp[0].mesh;
  MMG5_pSol  met  = parmesh->listgrp[0].met;

  // Warning : does it works with a communicator /= to MPI_COMM_WORLD?
  vtkSmartPointer<vtkMPIController> vtkController = vtkSmartPointer<vtkMPIController>::New();
  vtkController->Initialize();
  vtkMultiProcessController::SetGlobalController(vtkController);

  return MMG5_saveVtkMesh_i<vtkUnstructuredGrid,vtkXMLUnstructuredGridWriter,vtkXMLPUnstructuredGridWriter>
    (mesh,&met,sdata,1,0,parmesh->nprocs,parmesh->myrank,parmesh->info.root);

  MMG5_SAFE_FREE(mdata);
  MMG5_SAFE_FREE(sdata);

#endif
  return 1;
}

int PMMG_savePvtuMesh_and_allData(PMMG_pParMesh parmesh, const char * filename) {

  MMG5_pMesh mesh;
  MMG5_pSol  allSol[2];
  MMG5_pSol  met   = NULL;
  MMG5_pSol  field = NULL;
  int metricData = 0;

#ifndef USE_VTK
  if ( parmesh->myrank == parmesh->info.root ) {
    fprintf(stderr,"  ** VTK library not found. Unavailable file format.\n");
  }
  return -1;

#else
  char* mdata=NULL; // master    file name
  char* sdata=NULL; // secondary file name
  int i;
  MMG5_SAFE_CALLOC(mdata,strlen(filename)+6,char,return 0);
  MMG5_SAFE_CALLOC(sdata,strlen(filename)+6,char,return 0);

  strcpy(mdata,filename);
  char *ptr = MMG5_Get_filenameExt(mdata);
  *ptr = '\0'; // get basename

  // If the output *.pvtu filename has dots "." in the basename,
  // replace them with dashes "-".
  // Why? In VTK function SetFileName(filename), the first dot "." in the
  // filename is interpreted as the extension start. So, whatever the
  // user specifies after the first dot "." will be ignored by VTK. To overcome
  // this, dots are replaced by dashes.
  for(i=0;mdata[i]!='\0';i++) {
    if(mdata[i]=='.') {
      mdata[i] = '-';
    }
  }

  snprintf( sdata,strlen(mdata)+6, "%s.pvtu",mdata);

  mesh  = parmesh->listgrp[0].mesh;
  // Add met at the end of field to be able to save everything in the pvtu file
  if (parmesh->listgrp[0].met) {
    met = parmesh->listgrp[0].met;
    metricData = 1;
    allSol[0] = met;
  }
  else {
    allSol[0] = NULL;
  }
  if (parmesh->listgrp[0].field) {
    field = parmesh->listgrp[0].field;
    allSol[1] = field;
  }
  else {
    allSol[1] = NULL;
  }

  // Warning : does it works with a communicator /= to MPI_COMM_WORLD?
  vtkSmartPointer<vtkMPIController> vtkController = vtkSmartPointer<vtkMPIController>::New();
  vtkController->Initialize();
  vtkMultiProcessController::SetGlobalController(vtkController);

  return MMG5_saveVtkMesh_i<vtkUnstructuredGrid,vtkXMLUnstructuredGridWriter,vtkXMLPUnstructuredGridWriter>
    (mesh,allSol,sdata,metricData,0,parmesh->nprocs,parmesh->myrank,parmesh->info.root);

  MMG5_SAFE_FREE(mdata);
  MMG5_SAFE_FREE(sdata);

#endif
  return 1;
}

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
 * \file variadic_pmmg.c
 * \brief C variadic functions definitions for PMMG library.
 * \author Algiane Froehly (InriaSoft)
 * \version 1
 * \date 07 2018
 * \copyright GNU Lesser General Public License.
 *
 * \note This file contains some internal functions for the API, see
 * the \ref libparmmg.h header file for the documentation of all
 * the usefull user's API functions.
 *
 * variadic functions definitions for PMMG library.
 *
 */

#include "parmmg.h"

/**
 * \param argptr list of the type of structures that must be initialized inside
 * your parmesh and needed informations for ParMmg (mesh dimension and MPI
 * communicator).
 *
 * \param callFromC 1 if called from C API, 0 if called from the Fortran one.
 *
 * \a argptr contains at least a pointer toward a parmesh pointer preceeded by
 * the PMMG_ARG_pParMesh keyword and the mesh dimension preceeded by the
 * PMMG_ARG_dim keyword
 *
 * By default, ParMmg will initilized at least 1 mesh, 1 metric and the
 * MPI_COMM_WORLD_COMMUNICATOR inside your parmesh. Thus, the 2 following calls
 * are identicals:
 *
 * 1) MMG3D_Init_parmesh(PMMG_ARG_start,PMMG_ARG_ppParMesh,your_pParmesh_address,
 *    PMMG_ARG_pMesh,PMMG_ARG_pMet,
 *    PMMG_ARG_dim,mesh_dimension,PMMG_ARG_MPIComm,MPI_COMM_WORLD,PMMG_ARG_end)
 *
 * 2) MMG3D_Init_parmesh(PMMG_ARG_start,PMMG_ARG_ppParMesh,your_pParmesh_address,
 *    PMMG_ARG_dim,mesh_dimension,PMMG_ARG_end)
 *
 * \return 1 if success, 0 if fail
 *
 * Internal function for structure allocations (taking a va_list argument).
 *
 */

int PMMG_Init_parMesh_var_internal(va_list argptr, int callFromC ) {
  PMMG_pParMesh  *parmesh;
  MPI_Comm       comm;
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_pSol      met;
  size_t         memAv;
  int            typArg,dim,nsol,comm_f;
  int            parmeshCount,meshCount,metCount,dimCount,solCount,commCount;

  parmeshCount = 0;
  meshCount    = 0;
  metCount     = 0;
  solCount     = 0;
  dimCount     = 0;
  commCount    = 0;

  nsol = 0;
  dim  = 3;
  comm = MPI_COMM_WORLD;
  while ( (typArg = va_arg(argptr,int)) != PMMG_ARG_end )
  {
    switch ( typArg )
    {
    case(PMMG_ARG_ppParMesh):
      parmesh = va_arg(argptr,PMMG_pParMesh*);
      ++parmeshCount;
      break;
    case(PMMG_ARG_pMesh):
      ++meshCount;
      break;
    case(PMMG_ARG_pMet):
      ++metCount;
      break;
    case(PMMG_ARG_dim):
      ++dimCount;
      dim = va_arg(argptr,int);
      break;
    case(PMMG_ARG_MPIComm):
      ++commCount;
      if ( callFromC ) {
        comm = va_arg(argptr,MPI_Comm);
      }
      else {
        comm_f = va_arg(argptr,int);
        comm = MPI_Comm_f2c(comm_f);
      }
      break;
    default:
      fprintf(stderr,"\n  ## Error: PMMG_Init_parmesh:\n"
              " unexpected argument type: %s\n",PMMG_Get_pmmgArgName(typArg));
      return 0;
    }
  }

  if ( parmeshCount !=1 ) {
    fprintf(stderr,"\n  ## Error: PMMG_Init_parmesh:\n"
            " you need to initialize the parmesh structure that"
            " will contain your data (mesh, metric, communicator...\n");
    return 0;
  }

  if ( meshCount>1 ) {
    fprintf(stdout,"\n  ## Warning: PMMG_Init_parmesh:\n"
            " Only 1 mesh structure is allowed.\n");
  }
  if ( metCount>1 ) {
    fprintf(stdout,"\n  ## Warning: PMMG_Init_parmesh:\n"
            " Only 1 metric structure is allowed.\n");
  }
  if ( commCount>1 ) {
    fprintf(stdout,"\n  ## Warning: PMMG_Init_parmesh:\n"
            " More than 1 MPI communicator provided. Used the last one.\n");
  }
  if ( dimCount>1 ) {
    fprintf(stdout,"\n  ## Warning: PMMG_Init_parmesh:\n"
            " More than 1 dimension provided. Used the last one.\n");
  }
  else if ( !dimCount ) {
    fprintf(stderr,"\n  ## Error: PMMG_Init_parmesh:\n"
            " you need to provide the dimension of your mesh using the PMMG_dim"
            " keyword\n.");
    return 0;
  }

  if ( dim !=3 ) {
    fprintf(stderr,"\n  ## Error: PMMG_Init_parmesh:\n"
            " dimension other than 3D is not yet implemented.\n");
    return 0;
  }

  /* ParMesh allocation */
  assert ( (*parmesh == NULL) && "trying to initialize non empty parmesh" );
  MMG5_SAFE_CALLOC( *parmesh,1, PMMG_ParMesh,return 0 );

  if ( *parmesh == NULL ) {
    return 0;
  }

  /* Assign some values to memory related fields to begin working with */
  (*parmesh)->memGloMax = 4 * 1024L * 1024L;
  (*parmesh)->memMax = 4 * 1024L * 1024L;
  (*parmesh)->memCur = sizeof(PMMG_ParMesh);

  /** Init Group */
  grp = NULL;
  (*parmesh)->ngrp = 1;
  PMMG_CALLOC(*parmesh,(*parmesh)->listgrp,1,PMMG_Grp,
              "allocating groups container", goto fail_grplst );
  grp = &(*parmesh)->listgrp[0];
  grp->mesh  = NULL;
  grp->met   = NULL;
  grp->field = NULL;
  grp->disp  = NULL;
  grp->ls    = NULL;

  if ( 1 != MMG3D_Init_mesh( MMG5_ARG_start,
                             MMG5_ARG_ppMesh, &grp->mesh,
                             MMG5_ARG_ppMet, &grp->met,
                             MMG5_ARG_end ) )
    goto fail_mesh;

  PMMG_Init_parameters(*parmesh,comm);

  /* Check that Mmg does not take too much memory */
  mesh = (*parmesh)->listgrp[0].mesh;
  met  = (*parmesh)->listgrp[0].met;

  memAv = (*parmesh)->memGloMax - (*parmesh)->memMax;
  if ( (*parmesh)->listgrp[0].mesh->memMax > memAv ) {
    if ( (*parmesh)->listgrp[0].mesh->memCur > memAv ) {
      fprintf(stderr,"\n  ## Error: %s: %zu Mo of memory ",__func__,
              memAv/MMG5_MILLION);
      fprintf(stderr,"is not enough to initialize an empty mesh."
              " You need to ask %zu Mo minimum.\n",
              ((*parmesh)->listgrp[0].mesh->memCur + (*parmesh)->memMax) /MMG5_MILLION+1);
      fprintf(stderr,"\nTry to use the -m option to impose the maximal memory per process.\n");


      MMG3D_Free_all(MMG5_ARG_start,
                     MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,
                     MMG5_ARG_end);
      goto fail_mesh;
    }
    (*parmesh)->listgrp[0].mesh->memMax = memAv;
  }

  return 1;

fail_mesh:
    PMMG_DEL_MEM(*parmesh,(*parmesh)->listgrp,PMMG_Grp,
                 "deallocating groups container");

fail_grplst:
  (*parmesh)->ngrp = 0;
  (*parmesh)->memMax = 0;
  (*parmesh)->memCur = 0;
  MMG5_SAFE_FREE( *parmesh );

  return 0;
}

/**
 * \param argptr list of the parmmg structures that must be deallocated.
 *
 * \a argptr contains at least a pointer toward a \a PMMG_pParMesh structure
 * (that will contain the parmesh and identified by the PMMG_ARG_ppParMesh keyword)
 *
 * \return 0 if fail, 1 if success
 *
 * Deallocations of the parmmg structures before return
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 */
int PMMG_Free_all_var(va_list argptr)
{

  PMMG_pParMesh  *parmesh;
  int            typArg;
  int            parmeshCount;

  parmeshCount = 0;

  while ( (typArg = va_arg(argptr,int)) != PMMG_ARG_end )
  {
    switch ( typArg )
    {
    case(PMMG_ARG_ppParMesh):
      parmesh = va_arg(argptr,PMMG_pParMesh*);
      ++parmeshCount;
      break;
    default:
      fprintf(stdout,"\n  ## Warning: PMMG_Free_all:\n"
              " ignored argument: %s\n",PMMG_Get_pmmgArgName(typArg));
    }
  }

  PMMG_parmesh_Free_Comm( *parmesh );

  PMMG_parmesh_Free_Listgrp( *parmesh );

  (*parmesh)->memCur -= sizeof(PMMG_ParMesh);

  if ( (*parmesh)->info.imprim>5 || (*parmesh)->ddebug ) {
    printf("  MEMORY USED AT END (Bytes) %zu\n",(*parmesh)->memCur);
  }
  MMG5_SAFE_FREE(*parmesh);

  return 1;

}

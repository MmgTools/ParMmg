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
 * \file API_functions_pmmg.c
 * \brief C API functions definitions for PARMMG library.
 * \copyright GNU Lesser General Public License.
 *
 * C API for PARMMG library.
 *
 */
#include "parmmg.h"
#include "metis_pmmg.h"
#include "linkedlist_pmmg.h"

int PMMG_Init_parMesh(const int starter,...) {
  va_list argptr;
  int     ier;

  va_start(argptr, starter);

  ier = PMMG_Init_parMesh_var_internal(argptr,1);

  va_end(argptr);

  return ier;
}

/**
 * \param parmesh pointer toward a parmesh structure for memory count.
 * \param buffer pointer toward the buffer to store file name.
 * \param name string to store into \a buffer.
 * \param defname default string for \a buffer if \a name is empty.
 * \return 1 if success, 0 if fail
 *
 * Set the file name \a name into the \a buffer string. If \a name is empty, use
 * \a defname as default file name. If \a defname is empty too, the buffer
 * remains non allocated.
 *
 * \warning for internal use only.
 */
int PMMG_Set_name(PMMG_pParMesh parmesh,char **buffer,
                  const char* name, const char* defname) {

  if ( *buffer ) {
    PMMG_DEL_MEM(parmesh,*buffer,char,"buffer unalloc");
  }

  if ( name && strlen(name) ) {
    PMMG_MALLOC(parmesh,*buffer,strlen(name)+1,char,"name",return 0);
    strcpy(*buffer,name);
  }
  else if ( defname ) {
    PMMG_MALLOC(parmesh,*buffer,strlen(defname)+1,char,"defname",return 0);
    strcpy(*buffer,defname);
  }
  /* Remark: for solution fields: a non allocated buffer allows to detect that
   * user don't provide input fields */

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure for memory count.
 * \param buffer pointer toward the buffer to store file name.
 * \param name string to store into \a buffer.
 * \param defroot default root for the filename if unable to find it in
 * the mesh name.
 * \param in 1 if input file, 0 if output file.
 *
 * \return 1 if success, 0 if fail
 *
 * Set the solution name \a name into the \a buffer string. If \a name is empty
 * try to set a default name (depending on \a in):
 *  -# try get the mesh name (input one if in==1, output one otherwise)
 * to set a default filename.
 *  -# Second, if the mesh name is empty too, use \a defroot to set it.
 *
 * \remark If \a name is empty and if the input mesh is called <input>.mesh and
 * the output mesh <output>.mesh, we want:
 *  - an input metric, ls or displacement named <input>.sol;
 *  - an output metric named <ouput>.sol;
 *  - no ouput ls (no possibility to save it);
 *  - no ouput displacement (no possibility to save it);
 *  - no default input name for fields (it is mandatory to provide its name);
 *  - use the input field name <field>.sol to name the output field <field>.o.sol;
 * The current function allows to set default input name for met, ls, disp and
 * default ouput name for met.
 *
 * \warning for internal use only.
 */
static inline
int PMMG_Set_solName(PMMG_pParMesh parmesh,char **buffer,
                     const char* name,const char* defroot,const int8_t in) {
  int  ier,len;
  char *defname,*input;

  if ( *buffer ) {
    PMMG_DEL_MEM(parmesh,*buffer,char,"buffer unalloc");
  }

  /* Find default field name if solin isn't provided */
  defname = NULL;
  if ( in ) {
    input = parmesh->meshin;
  }
  else {
    input = parmesh->meshout;
  }

  if ( (!name) || (!*name) ) {
    /* Remove .mesh extension */
    defname = MMG5_Remove_ext ( input,".mesh" );

    if ( defname ) {
      /* Add .sol extension to the mesh name */
      /* Remark: multiple solution structures (metric, ls...) may have the same
       * default name but in practice it doesn't happens as it is mandatory that
       * at least a part of the filename must be provided (otherwise the
       * multiple input is not possible). */
      PMMG_REALLOC(parmesh,defname,strlen(defname)+5,strlen(defname)+1,char,"",return 0);
      strncat ( defname,".sol",5 );
    }
    else {
      /* Use the provided defroot name */
      assert ( defroot && *defroot );
      if ( in ) {
        len = strlen(defroot)+5; // defroot+".sol"
      } else {
        len = strlen(defroot)+7; // defroot+".o.sol"
      }
      PMMG_CALLOC(parmesh,defname,len,char,"",return 0);
      strncat ( defname, defroot,strlen(defroot)+1 );
      if ( in ) {
        strncat( defname,".sol",5);
      }
      else {
        strncat( defname,".o.sol",7);
      }
    }
  }

  ier = PMMG_Set_name(parmesh,buffer,name,defname);

  PMMG_DEL_MEM ( parmesh,defname,char,"" );

  return ier;
}

int PMMG_Set_inputMeshName(PMMG_pParMesh parmesh, const char* meshin) {
  MMG5_pMesh mesh;
  int        k,ier;

  ier = PMMG_Set_name(parmesh,&parmesh->meshin,meshin,"mesh.mesh");

  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    ier  = MG_MIN( ier, MMG3D_Set_inputMeshName(mesh,parmesh->meshin) );
  }
  return ier;
}

int PMMG_Set_inputSolsName(PMMG_pParMesh parmesh, const char* solin) {
  MMG5_pMesh mesh;
  MMG5_pSol  sol,psl;
  int        k,i,ier;

  /* No default name for solution fields as providing a file name is mandatory
   * to know that we have solutions. */
  ier = PMMG_Set_name ( parmesh,&parmesh->fieldin,solin,NULL);

  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    sol  = parmesh->listgrp[k].field;

    for ( i=0; i<mesh->nsols; ++i ) {
      psl = sol + i;
      ier  = MG_MIN( ier, MMG3D_Set_inputSolName(mesh,psl,parmesh->fieldin) );
    }
  }
  return ier;
}

int PMMG_Set_inputMetName(PMMG_pParMesh parmesh, const char* metin) {
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  int        k,ier;

  ier = PMMG_Set_solName ( parmesh,&parmesh->metin,metin,"met",1 );

  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    met  = parmesh->listgrp[k].met;

    /* Check that the metric structure is allocated */
    if ( !met ) {
      fprintf(stderr, "\n  ## Error: %s: metric structure must be initialized.\n" ,__func__);
      fprintf(stderr, "              Please add the PMMG_ARG_pMet argument to"
              " your call to the PMMG_Init_parMesh function.\n");
      return 0;
    }

    ier  = MG_MIN( ier, MMG3D_Set_inputSolName(mesh,met,parmesh->metin) );
  }
  return ier;
}

int PMMG_Set_inputLsName(PMMG_pParMesh parmesh, const char* lsin) {
  MMG5_pMesh mesh;
  MMG5_pSol  ls;
  int        k,ier;

  ier = PMMG_Set_solName ( parmesh,&parmesh->lsin,lsin,"ls",1 );

  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    ls   = parmesh->listgrp[k].ls;

    /* Check that the ls structure is allocated */
    if ( !ls ) {
      fprintf(stderr, "\n  ## Error: %s: level-set structure must be initialized.\n" ,__func__);
      fprintf(stderr, "              Please add the PMMG_ARG_pLs argument to"
              " your call to the PMMG_Init_parMesh function.\n");
      return 0;
    }

    ier  = MG_MIN( ier, MMG3D_Set_inputSolName(mesh,ls,parmesh->lsin) );
  }
  return ier;
}

int PMMG_Set_inputDispName(PMMG_pParMesh parmesh, const char* dispin) {
  MMG5_pMesh mesh;
  MMG5_pSol  disp;
  int        k,ier;

  ier = PMMG_Set_solName ( parmesh,&parmesh->dispin,dispin,"disp",1 );

  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    disp = parmesh->listgrp[k].disp;

    /* Check that the displacement structure is allocated */
    if ( !disp ) {
      fprintf(stderr, "\n  ## Error: %s: displacement structure must be initialized.\n" ,__func__);
      fprintf(stderr, "              Please add the PMMG_ARG_pDisp argument to"
              " your call to the PMMG_Init_parMesh function.\n");
      return 0;
    }

    ier  = MG_MIN( ier, MMG3D_Set_inputSolName(mesh,disp,parmesh->dispin) );
  }
  return ier;
}

int PMMG_Set_outputMeshName(PMMG_pParMesh parmesh, const char* meshout) {
  MMG5_pMesh mesh;
  int        k,ier;
  char       *defname;

  if ( parmesh->meshout ) {
    PMMG_DEL_MEM(parmesh,parmesh->meshout,char,"fieldout unalloc");
  }

  defname = NULL;
  /* Find default field name if meshout isn't provided */
  if ( (!meshout) || !(*meshout) ) {
    defname = MMG5_Remove_ext ( parmesh->meshin,".mesh" );
  }
  if ( defname && *defname ) {
    /* Add .o.mesh extension */
    PMMG_REALLOC(parmesh,defname,strlen(defname)+8,strlen(defname)+1,char,"",return 0);
    strncat ( defname,".o.mesh",8 );
  }
  else {
    PMMG_MALLOC(parmesh,defname,strlen("mesh.o.mesh")+1,char,"default mesh name",return 0);
    strncpy( defname,"mesh.o.mesh",strlen("mesh.o.mesh")+1);
  }

  ier = PMMG_Set_name(parmesh,&parmesh->meshout,meshout,defname);

  PMMG_DEL_MEM ( parmesh,defname,char,"" );

  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    ier  = MG_MIN( ier, MMG3D_Set_outputMeshName(mesh,parmesh->meshout) );
  }
  return ier;
}

int PMMG_Set_outputSolsName(PMMG_pParMesh parmesh, const char* solout) {
  MMG5_pMesh mesh;
  MMG5_pSol  sol,psl;
  int        k,i,ier,pathlen,baselen;
  char      *basename,*path,*nopath;

  /* If \a solout is not provided we want to use the basename of the input field
   * name and the path of the output mesh name */
  if ( parmesh->fieldout ) {
    PMMG_DEL_MEM(parmesh,parmesh->fieldout,char,"fieldout unalloc");
  }

  if ( (!solout) || (!*solout) ) {

    if ( (!parmesh->meshout) || (!*parmesh->meshout) ) {
      fprintf(stderr, "  ## Error: %s: please, provide an output mesh"
              " name before calling this function without string.\n",
              __func__);
      return 0;
    }

    /* Get input field base name and remove .mesh extension */
    if ( (!parmesh->fieldin) || (!*parmesh->fieldin) ) {
      fprintf(stderr, "  ## Error: %s: please, provide an input field"
              " name before calling this function without string.\n",
              __func__);
      return 0;
    }

    path   = MMG5_Get_path(parmesh->meshout);
    nopath = MMG5_Get_basename(parmesh->fieldin);
    basename = MMG5_Remove_ext ( nopath,".sol" );

    pathlen = baselen = 0;
    if ( path ) pathlen = strlen(path)+1;
    if ( basename ) baselen = strlen(basename);
    PMMG_MALLOC(parmesh,parmesh->fieldout,pathlen+baselen+1,char,"fieldout",return 0);
    if ( pathlen ) {
      strncpy(parmesh->fieldout,path,pathlen-1);
      parmesh->fieldout[pathlen-1] = MMG5_PATHSEP;
    }
    if ( baselen ) {
      strncpy(parmesh->fieldout+pathlen,basename,baselen);
      parmesh->fieldout[pathlen+baselen] = '\0';
    }

    if ( parmesh->fieldout ) {
      /* Add .o.sol extension */
      PMMG_REALLOC(parmesh,parmesh->fieldout,strlen(parmesh->fieldout)+7,
                   strlen(parmesh->fieldout)+1,char,"fieldout",return 0);
      strncat ( parmesh->fieldout,".o.sol",7 );
    }

    MMG5_SAFE_FREE ( path );
    free ( nopath ); nopath = NULL;
    MMG5_SAFE_FREE ( basename );
  }
  else {
    PMMG_MALLOC(parmesh,parmesh->fieldout,strlen(solout)+1,char,"fieldout",return 0);
    strcpy(parmesh->fieldout,solout);
  }

  ier = 1;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    sol  = parmesh->listgrp[k].field;
    for ( i=0; i<mesh->nsols; ++i ) {
      psl = sol + i;
      ier  = MG_MIN( ier, MMG3D_Set_outputSolName(mesh,psl,parmesh->fieldout) );
    }
  }
  return ier;
}

int PMMG_Set_outputMetName(PMMG_pParMesh parmesh, const char* metout) {
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  int        k,ier;

  ier = PMMG_Set_solName ( parmesh,&parmesh->metout,metout,"met",0 );

  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    met  = parmesh->listgrp[k].met;
    ier  = MG_MIN( ier, MMG3D_Set_outputSolName(mesh,met,parmesh->metout) );
  }
  return ier;
}

void PMMG_Init_parameters(PMMG_pParMesh parmesh,MPI_Comm comm) {
  MMG5_pMesh mesh;
  size_t     mem;
  int        k,flag;

  memset(&parmesh->info    ,0, sizeof(PMMG_Info));

  parmesh->info.mem                = PMMG_UNSET; /* [n/-1]   ,Set memory size to n Mbytes/keep the default value */
  parmesh->info.root               = PMMG_NUL;

  parmesh->ddebug                  = PMMG_NUL;
  parmesh->iter                    = PMMG_UNSET;
  parmesh->niter                   = PMMG_NITER;
  parmesh->info.fem                = MMG5_FEM;
  parmesh->info.repartitioning     = PMMG_REDISTRIBUTION_mode;
  parmesh->info.ifc_layers         = PMMG_MVIFCS_NLAYERS;
  parmesh->info.grps_ratio         = PMMG_GRPS_RATIO;
  parmesh->info.nobalancing        = MMG5_OFF;
  parmesh->info.loadbalancing_mode = PMMG_LOADBALANCING_metis;
  parmesh->info.contiguous_mode    = PMMG_CONTIG_DEF;
  parmesh->info.target_mesh_size   = PMMG_REMESHER_TARGET_MESH_SIZE;
  parmesh->info.metis_ratio        = PMMG_RATIO_MMG_METIS;
  parmesh->info.API_mode           = PMMG_APIDISTRIB_faces;
  parmesh->info.globalNum          = PMMG_NUL;
  parmesh->info.sethmin            = PMMG_NUL;
  parmesh->info.sethmax            = PMMG_NUL;
  parmesh->info.fmtout             = PMMG_FMT_Unknown;

  /* Init MPI data */
  parmesh->comm   = comm;

  MPI_Initialized(&flag);
  parmesh->size_shm = 1;
  if ( flag ) {
    MPI_Comm_size( parmesh->comm, &parmesh->nprocs );
    MPI_Comm_rank( parmesh->comm, &parmesh->myrank );
  }
  else {
    parmesh->nprocs = 1;
    parmesh->myrank = PMMG_NUL;
  }

  /* ParMmg verbosity */
  if ( parmesh->myrank==parmesh->info.root ) {
    parmesh->info.imprim = PMMG_IMPRIM;
  }
  else {
    parmesh->info.imprim = PMMG_UNSET;
  }
  parmesh->info.imprim0  = PMMG_IMPRIM;

  /* Set the mmg verbosity to -1 (no output) */
  parmesh->info.mmg_imprim = PMMG_MMG_IMPRIM;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    mesh->info.imprim = MG_MIN ( parmesh->info.imprim,PMMG_MMG_IMPRIM );
  }

  /* Default memory */
  PMMG_parmesh_SetMemGloMax( parmesh );
  PMMG_parmesh_SetMemMax( parmesh );

}

int PMMG_Set_meshSize(PMMG_pParMesh parmesh, int np, int ne, int nprism, int nt,
                      int nquad, int na){
  MMG5_pMesh mesh;
  int        ier;

  assert ( parmesh->ngrp == 1 );

  ier = 1;
  mesh = parmesh->listgrp[0].mesh;

  /* Check input data and set mesh->ne/na/np/nt to the suitable values */
  if ( !MMG3D_setMeshSize_initData(mesh,np,ne,nprism,nt,nquad,na) )
    return 0;

  /* Check the -m option */
  if( parmesh->info.mem > 0) {
    assert ( mesh->info.mem > 0 );
    if(mesh->info.mem < 39) {
      fprintf(stderr,"\n  ## Error: %s: not enough memory per mesh  %d\n",__func__,
              mesh->info.mem);
      return 0;
    }

    if((mesh->npmax < mesh->np || mesh->ntmax < mesh->nt || mesh->nemax < mesh->ne)) {
       if ( !MMG3D_memOption(mesh) )  return 0;
    }
  } else {
    if ( !MMG3D_memOption(mesh) )  return 0;
  }

  /* Mesh allocation and linkage */
  if ( !MMG3D_setMeshSize_alloc( mesh ) ) return 0;

  return ier;
}

int PMMG_Set_solsAtVerticesSize(PMMG_pParMesh parmesh, int nsols,int nentities,
                                int *typSol){
  MMG5_pMesh mesh;
  MMG5_pSol  *sol;
  int        ier;

  ier = 1;

  mesh = parmesh->listgrp[0].mesh;
  sol  = &parmesh->listgrp[0].field;

  ier  = MMG3D_Set_solsAtVerticesSize(mesh,sol,nsols,nentities,typSol);

  return ier;
}

int PMMG_Set_metSize(PMMG_pParMesh parmesh,int typEntity,int np,int typSol){
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  int        ier;

  ier = 1;

  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  ier  = MMG3D_Set_solSize(mesh,met,typEntity,np,typSol);

  return ier;
}

int PMMG_Set_iparameter(PMMG_pParMesh parmesh, int iparam,int val) {
  MMG5_pMesh  mesh;
  MMG5_pSol   met;
  size_t      mem;
  int         k;

  switch ( iparam ) {
  case PMMG_IPARAM_verbose :
    if ( parmesh->myrank==parmesh->info.root ) {
      parmesh->info.imprim = val;
    }
    parmesh->info.imprim0 = val;

    break;
  case PMMG_IPARAM_mmgVerbose :
    parmesh->info.mmg_imprim = val;

    /* Set the mmg verbosity */
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_verbose,val) ) return 0;
    }
    break;

  case PMMG_IPARAM_mem :
    if ( val <= 0 ) {
      fprintf( stdout,
        "  ## Warning: maximal memory per process must be strictly positive.\n");
      fprintf(stdout,"  Reset to default value.\n");
    } else {
      parmesh->info.mem = val;
    }
    PMMG_parmesh_SetMemGloMax(parmesh);
    parmesh->memMax = parmesh->memGloMax;
    mem = parmesh->memGloMax;

    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;

      /* Mesh reallocation if needed */
      if ( (mesh->memCur >> MMG5_BITWIZE_MB_TO_B) > mem ) {
        fprintf(stderr,"\n  ## Error: %s: Maximal memory must be setted "
                "before reading the mesh.\n",__func__);
        return 0;
      }

      if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_mem,(int)mem) ) return 0;
    }
    break;
  case PMMG_IPARAM_meshSize :
    parmesh->info.target_mesh_size = val;
    break;
  case PMMG_IPARAM_nobalancing :
    parmesh->info.nobalancing = val;
    break;
  case PMMG_IPARAM_metisRatio :
    parmesh->info.metis_ratio = val;
    break;
  case PMMG_IPARAM_ifcLayers :
    parmesh->info.ifc_layers = val;
    break;
  case PMMG_IPARAM_APImode :
    parmesh->info.API_mode = val;
    break;
  case PMMG_IPARAM_globalNum :
    parmesh->info.globalNum = val;
    break;
  case PMMG_IPARAM_niter :
    parmesh->niter = val;
    break;

#ifndef PATTERN
  case PMMG_IPARAM_octree :
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_octree,val) ) return 0;
    }
    break;
#endif
  case PMMG_IPARAM_debug :
    parmesh->ddebug = val;
    break;

  case PMMG_IPARAM_distributedOutput :

    if ( val == 1 ) {
      parmesh->info.fmtout = PMMG_FMT_Distributed;
    }
    else if ( val == 0 ) {
      parmesh->info.fmtout = PMMG_FMT_Centralized;
    }
    break;

  case PMMG_IPARAM_mmgDebug :

    /* Set the mmg debug mode */
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_debug,val) ) return 0;
    }
    break;

  case PMMG_IPARAM_angle :
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_angle,val) ) return 0;
    }
    break;
  case PMMG_IPARAM_iso :
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_iso,val) ) return 0;
    }
    break;
  case PMMG_IPARAM_lag :
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_lag,val) ) return 0;
    }
    break;

  case PMMG_IPARAM_nofem :
    parmesh->info.fem    = (val==1)? 0 : 1;
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      met  = parmesh->listgrp[k].met;
      if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_nofem,val) ) return 0;
    }
    break;
  case PMMG_IPARAM_opnbdy :
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_opnbdy,val) ) return 0;
    }
    if( val ) {
#warning opnbdy not supported with surface adaptation
      fprintf(stderr," ## Warning: Surface adaptation not supported with opnbdy."
          "\nSetting nosurf on.\n");
      for ( k=0; k<parmesh->ngrp; ++k ) {
        mesh = parmesh->listgrp[k].mesh;
        if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_nosurf,val) ) return 0;
      }
    }
    break;
  case PMMG_IPARAM_optim :
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_optim,val) ) return 0;
    }
    break;
  case PMMG_IPARAM_optimLES :
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_optimLES,val) ) return 0;
    }
    break;
  case PMMG_IPARAM_noinsert :
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_noinsert,val) ) return 0;
    }
    break;
  case PMMG_IPARAM_noswap :
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_noswap,val) ) return 0;
    }
    break;
  case PMMG_IPARAM_nomove :
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_nomove,val) ) return 0;
    }
    break;
  case PMMG_IPARAM_nosurf :
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
#warning opnbdy not supported with surface adaptation
      if( !val && mesh->info.opnbdy )
        fprintf(stderr," ## Warning: Surface adaptation not supported with opnbdy."
          "\nCannot set nosurf off.\n");
      else if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_nosurf,val) ) return 0;
    }
    break;
  case PMMG_IPARAM_numberOfLocalParam :
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_numberOfLocalParam,val) )
        return 0;
    }
    break;
  case PMMG_IPARAM_anisosize :
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      met  = parmesh->listgrp[k].met;
      if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_anisosize,val) ) return 0;
    }
    break;
  default :
    fprintf(stderr,"  ## Error: unknown type of parameter\n");
    return 0;
  }

  return 1;
}

int PMMG_Set_dparameter(PMMG_pParMesh parmesh, int dparam,double val){
  MMG5_pMesh  mesh;
  int         k;

  switch ( dparam ) {
  case PMMG_DPARAM_angleDetection :
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_dparameter(mesh,NULL,MMG3D_DPARAM_angleDetection,val) ) {
        return 0;
      }
    }
    break;
  case PMMG_DPARAM_hmin :
    parmesh->info.sethmin  = 1;
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_dparameter(mesh,NULL,MMG3D_DPARAM_hmin,val) ) {
        return 0;
      }
    }
    break;
  case PMMG_DPARAM_hmax :
    parmesh->info.sethmax  = 1;
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_dparameter(mesh,NULL,MMG3D_DPARAM_hmax,val) ) {
        return 0;
      }
    }
    break;
  case PMMG_DPARAM_hsiz :
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_dparameter(mesh,NULL,MMG3D_DPARAM_hsiz,val) ) {
        return 0;
      }
    }
    break;
  case PMMG_DPARAM_hgrad :
   for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_dparameter(mesh,NULL,MMG3D_DPARAM_hgrad,val) ) {
        return 0;
      }
    }
    break;
  case PMMG_DPARAM_hgradreq :
   for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_dparameter(mesh,NULL,MMG3D_DPARAM_hgradreq,val) ) {
        return 0;
      }
    }
    break;
  case PMMG_DPARAM_hausd :
    if ( val <=0 ) {
      fprintf(stderr,"\n  ## Error: %s: hausdorff number must be strictly"
              " positive.\n",__func__);
      return 0;
    }
    else {
      for ( k=0; k<parmesh->ngrp; ++k ) {
        mesh = parmesh->listgrp[k].mesh;
        if ( !MMG3D_Set_dparameter(mesh,NULL,MMG3D_DPARAM_hausd,val) ) {
          return 0;
        }
      }
    }
    break;
  case PMMG_DPARAM_ls :
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_dparameter(mesh,NULL,MMG3D_DPARAM_ls,val) ) {
        return 0;
      }
    }
    break;
  case PMMG_DPARAM_groupsRatio :
    parmesh->info.grps_ratio = val;
    break;
  default :
    fprintf(stderr,"  ## Error: unknown type of parameter\n");
    return 0;
  }

  return 1;
}

int PMMG_Set_vertex(PMMG_pParMesh parmesh, double c0, double c1, double c2,
                    int ref, int pos){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_vertex(parmesh->listgrp[0].mesh, c0, c1, c2, ref, pos));
}

int PMMG_Set_vertices(PMMG_pParMesh parmesh, double *vertices,int *refs){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_vertices(parmesh->listgrp[0].mesh, vertices, refs));
}

int PMMG_Set_tetrahedron(PMMG_pParMesh parmesh, int v0, int v1, int v2, int v3,
                         int ref, int pos){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_tetrahedron(parmesh->listgrp[0].mesh, v0, v1, v2, v3, ref, pos));
}

int PMMG_Set_tetrahedra(PMMG_pParMesh parmesh, int *tetra, int *refs){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_tetrahedra(parmesh->listgrp[0].mesh, tetra, refs));
}

int PMMG_Set_prism(PMMG_pParMesh parmesh, int v0, int v1, int v2,
                   int v3, int v4, int v5, int ref, int pos){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_prism(parmesh->listgrp[0].mesh, v0, v1, v2, v3, v4, v5, ref, pos));
}

int PMMG_Set_prisms(PMMG_pParMesh parmesh, int *prisms, int *refs){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_prisms(parmesh->listgrp[0].mesh, prisms, refs));
}

int PMMG_Set_triangle(PMMG_pParMesh parmesh, int v0, int v1, int v2,
                      int ref,int pos){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_triangle(parmesh->listgrp[0].mesh, v0, v1, v2, ref, pos));
}

int PMMG_Set_triangles(PMMG_pParMesh parmesh, int *tria, int *refs){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_triangles(parmesh->listgrp[0].mesh, tria, refs));
}

int PMMG_Set_quadrilateral(PMMG_pParMesh parmesh, int v0, int v1,
                           int v2, int v3, int ref,int pos){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_quadrilateral(parmesh->listgrp[0].mesh, v0, v1, v2, v3, ref, pos));
}

int PMMG_Set_quadrilaterals(PMMG_pParMesh parmesh, int *quads, int *refs){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_quadrilaterals(parmesh->listgrp[0].mesh, quads, refs));
}

int PMMG_Set_edge(PMMG_pParMesh parmesh, int v0, int v1, int ref, int pos){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_edge(parmesh->listgrp[0].mesh, v0, v1, ref, pos));
}
int PMMG_Set_edges(PMMG_pParMesh parmesh, int *edges, int *refs) {
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_edges(parmesh->listgrp[0].mesh, edges, refs));
}
int PMMG_Set_corner(PMMG_pParMesh parmesh, int k){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_corner(parmesh->listgrp[0].mesh, k));
}

int PMMG_Set_requiredVertex(PMMG_pParMesh parmesh, int k){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_requiredVertex(parmesh->listgrp[0].mesh, k));
}

int PMMG_Set_requiredTetrahedron(PMMG_pParMesh parmesh, int k){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_requiredTetrahedron(parmesh->listgrp[0].mesh, k));
}

int PMMG_Set_requiredTetrahedra(PMMG_pParMesh parmesh, int *reqIdx, int nreq){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_requiredTetrahedra(parmesh->listgrp[0].mesh, reqIdx, nreq));
}

int PMMG_Set_requiredTriangle(PMMG_pParMesh parmesh, int k){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_requiredTriangle(parmesh->listgrp[0].mesh, k));
}

int PMMG_Set_requiredTriangles(PMMG_pParMesh parmesh, int *reqIdx, int nreq){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_requiredTriangles(parmesh->listgrp[0].mesh, reqIdx, nreq));
}

int PMMG_Set_ridge(PMMG_pParMesh parmesh, int k){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_ridge(parmesh->listgrp[0].mesh, k));
}

int PMMG_Set_requiredEdge(PMMG_pParMesh parmesh, int k){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_requiredEdge(parmesh->listgrp[0].mesh, k));
}

int PMMG_Set_normalAtVertex(PMMG_pParMesh parmesh, int k, double n0, double n1,
                            double n2){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_normalAtVertex(parmesh->listgrp[0].mesh, k, n0, n1, n2));
}

int PMMG_Set_ithSols_inSolsAtVertices(PMMG_pParMesh parmesh,int i, double* s){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_ithSols_inSolsAtVertices(parmesh->listgrp[0].field,i,s));
}

int PMMG_Set_ithSol_inSolsAtVertices(PMMG_pParMesh parmesh,int i,double *s,int pos){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_ithSol_inSolsAtVertices(parmesh->listgrp[0].field,i, s, pos));
}

int PMMG_Set_scalarMet(PMMG_pParMesh parmesh, double m,int pos){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_scalarSol(parmesh->listgrp[0].met, m, pos));
}

int PMMG_Set_scalarMets(PMMG_pParMesh parmesh, double *met){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_scalarSols(parmesh->listgrp[0].met, met));
}

int PMMG_Set_vectorMet(PMMG_pParMesh parmesh, double vx,double vy, double vz,
                       int pos){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_vectorSol(parmesh->listgrp[0].met, vx, vy, vz, pos));
}

int PMMG_Set_vectorMets(PMMG_pParMesh parmesh, double *mets){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_vectorSols(parmesh->listgrp[0].met, mets));
}

int PMMG_Set_tensorMet(PMMG_pParMesh parmesh, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_tensorSol(parmesh->listgrp[0].met, m11, m12, m13, m22, m23, m33, pos));
}

int PMMG_Set_tensorMets(PMMG_pParMesh parmesh, double *mets){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_tensorSols(parmesh->listgrp[0].met, mets));
}

int PMMG_Get_meshSize(PMMG_pParMesh parmesh,int *np,int *ne,int *nprism,int *nt,
                      int *nquad,int *na){
  MMG5_pMesh mesh;
  int        ier;

  assert ( parmesh->ngrp == 1 );

  mesh = parmesh->listgrp[0].mesh;
  ier = MMG3D_Get_meshSize(mesh,np,ne,nprism,nt,nquad,na);

  return ier;
}

int PMMG_Get_solsAtVerticesSize(PMMG_pParMesh parmesh, int *nsols,int *nentities,
                                int *typSol){
  MMG5_pMesh mesh;
  MMG5_pSol  *sol;
  int        ier;

  ier = 1;

  if ( nsols )
    *nsols = 0;
  if ( nentities )
    *nentities = 0;

  if ( parmesh->listgrp && parmesh->listgrp[0].mesh && parmesh->listgrp[0].field ) {
    mesh = parmesh->listgrp[0].mesh;
    sol  = &parmesh->listgrp[0].field;

    ier  = MMG3D_Get_solsAtVerticesSize(mesh,sol,nsols,nentities,typSol);
  }

  return ier;
}

int PMMG_Get_metSize(PMMG_pParMesh parmesh,int *typEntity,int *np,int *typSol){
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  int        ier;

  ier = 1;

  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  ier  = MMG3D_Get_solSize(mesh,met,typEntity,np,typSol);

  return ier;
}

int PMMG_Get_vertex(PMMG_pParMesh parmesh, double* c0, double* c1, double* c2,
                    int* ref,int* isCorner, int* isRequired){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_vertex(parmesh->listgrp[0].mesh, c0, c1, c2, ref, isCorner, isRequired));
}

int PMMG_Get_vertices(PMMG_pParMesh parmesh, double* vertices, int* refs,
                      int* areCorners, int* areRequired){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_vertices(parmesh->listgrp[0].mesh, vertices, refs, areCorners, areRequired));
}

int PMMG_Get_tetrahedron(PMMG_pParMesh parmesh, int* v0, int* v1, int* v2,
                         int* v3,int* ref, int* isRequired){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_tetrahedron(parmesh->listgrp[0].mesh, v0, v1, v2, v3, ref, isRequired));
}

int PMMG_Get_tetrahedra(PMMG_pParMesh parmesh, int* tetra,int* refs,
                        int* areRequired){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_tetrahedra(parmesh->listgrp[0].mesh, tetra, refs, areRequired));
}

int PMMG_Get_prism(PMMG_pParMesh parmesh, int* v0, int* v1, int* v2,
                   int* v3,int* v4,int* v5,int* ref, int* isRequired){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_prism(parmesh->listgrp[0].mesh, v0, v1, v2, v3, v4, v5, ref, isRequired));
}

int PMMG_Get_prisms(PMMG_pParMesh parmesh, int* prisms,int* refs,
                    int* areRequired){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_prisms(parmesh->listgrp[0].mesh, prisms, refs, areRequired));
}

int PMMG_Get_triangle(PMMG_pParMesh parmesh, int* v0, int* v1, int* v2, int* ref,
                      int* isRequired){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_triangle(parmesh->listgrp[0].mesh, v0, v1, v2, ref, isRequired));
}

int PMMG_Get_triangles(PMMG_pParMesh parmesh, int* tria, int* refs,
                       int* areRequired){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_triangles(parmesh->listgrp[0].mesh, tria, refs, areRequired));
}

int PMMG_Get_quadrilateral(PMMG_pParMesh parmesh, int* v0, int* v1, int* v2,int* v3,
                            int* ref, int* isRequired){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_quadrilateral(parmesh->listgrp[0].mesh, v0, v1, v2, v3, ref, isRequired));
}

int PMMG_Get_quadrilaterals(PMMG_pParMesh parmesh, int* quads, int* refs,
                            int* areRequired){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_quadrilaterals(parmesh->listgrp[0].mesh, quads, refs, areRequired));
}

int PMMG_Get_edge(PMMG_pParMesh parmesh, int* e0, int* e1, int* ref,
                   int* isRidge, int* isRequired){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_edge(parmesh->listgrp[0].mesh, e0, e1, ref, isRidge, isRequired));
}
int PMMG_Get_edges(PMMG_pParMesh parmesh, int* edges,int *refs,int* areRidges,int* areRequired) {
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_edges(parmesh->listgrp[0].mesh,edges,refs,areRidges,areRequired));
}

int PMMG_Get_normalAtVertex(PMMG_pParMesh parmesh, int k, double *n0, double *n1,
                              double *n2){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_normalAtVertex(parmesh->listgrp[0].mesh, k, n0, n1, n2));
}

int PMMG_Get_ithSol_inSolsAtVertices(PMMG_pParMesh parmesh,int i,double* s,int pos){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_ithSol_inSolsAtVertices(parmesh->listgrp[0].field,i,s,pos));
}

int PMMG_Get_ithSols_inSolsAtVertices(PMMG_pParMesh parmesh,int i,double* s){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_ithSols_inSolsAtVertices(parmesh->listgrp[0].field,i,s));

}

int PMMG_Get_scalarMet(PMMG_pParMesh parmesh, double* m){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_scalarSol(parmesh->listgrp[0].met, m));
}

int PMMG_Get_scalarMets(PMMG_pParMesh parmesh, double* m){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_scalarSols(parmesh->listgrp[0].met, m));
}

int PMMG_Get_vectorMet(PMMG_pParMesh parmesh, double* vx, double* vy, double* vz){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_vectorSol(parmesh->listgrp[0].met, vx, vy, vz));
}

int PMMG_Get_vectorMets(PMMG_pParMesh parmesh, double* mets){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_vectorSols(parmesh->listgrp[0].met, mets));
}

int PMMG_Get_tensorMet(PMMG_pParMesh parmesh, double *m11,double *m12, double *m13,
                       double *m22,double *m23, double *m33){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_tensorSol(parmesh->listgrp[0].met, m11, m12, m13, m22, m23, m33));
}

int PMMG_Get_tensorMets(PMMG_pParMesh parmesh, double *mets){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_tensorSols(parmesh->listgrp[0].met, mets));
}

int PMMG_Set_numberOfNodeCommunicators(PMMG_pParMesh parmesh, int next_comm) {

  PMMG_CALLOC(parmesh,parmesh->ext_node_comm,next_comm,PMMG_Ext_comm,
              "allocate ext_comm ",return 0);
  parmesh->next_node_comm = next_comm;

  PMMG_CALLOC(parmesh,parmesh->int_node_comm,1,PMMG_Int_comm,
              "allocate int_comm ",return 0);

  return 1;
}

int PMMG_Set_numberOfFaceCommunicators(PMMG_pParMesh parmesh, int next_comm) {

  PMMG_CALLOC(parmesh,parmesh->ext_face_comm,next_comm,PMMG_Ext_comm,
              "allocate ext_comm ",return 0);
  parmesh->next_face_comm = next_comm;

  PMMG_CALLOC(parmesh,parmesh->int_face_comm,1,PMMG_Int_comm,
              "allocate int_comm ",return 0);

  return 1;
}

int PMMG_Set_ithNodeCommunicatorSize(PMMG_pParMesh parmesh, int ext_comm_index, int color_out, int nitem) {
  PMMG_pExt_comm pext_comm;

  /* Return if communicator index doesn't exist */
  if( (ext_comm_index < 0 ) || (ext_comm_index > parmesh->next_node_comm-1) ) {
    fprintf(stderr,"\n ## Error: function %s on proc %d: Communicator index should be in the range [0,ncomm-1], %d passed when %d communicators were asked. Please check your interface.\n",__func__,parmesh->myrank,ext_comm_index,parmesh->next_node_comm);
    return 0;
  }

  /* Return if communicator color doesn't exist or is myrank */
  if( (color_out < 0) || (color_out > parmesh->nprocs-1) || (color_out == parmesh->myrank) ) {
    fprintf(stderr,"\n ## Error: function %s on proc %d: Communicator color should be a proc ID in the range [0,nprocs-1] and different from my rank, %d passed while %d MPI procs are seen. Please check your interface.\n",__func__,parmesh->myrank,color_out,parmesh->nprocs);
    return 0;
  }

  /* Return if communicator size non-positive */
  if( nitem <= 0 ) {
    fprintf(stderr,"\n ## Error: function %s on proc %d: Communicator size should be strictly positive, %d passed. Please check your interface.\n",__func__,parmesh->myrank,nitem);
    return 0;
  }


  /* Get the node communicator */
  pext_comm = &parmesh->ext_node_comm[ext_comm_index];

  /* Set colors */
  pext_comm->color_in  = parmesh->myrank;
  pext_comm->color_out = color_out;

  /* Allocate communicator */
  PMMG_CALLOC(parmesh,pext_comm->int_comm_index,nitem,int,
                  "allocate int_comm_index",return 0);
  PMMG_CALLOC(parmesh,pext_comm->itosend,nitem,int,"allocate itosend",return 0);
  PMMG_CALLOC(parmesh,pext_comm->itorecv,nitem,int,"allocate itorecv",return 0);
  pext_comm->nitem          = nitem;
  pext_comm->nitem_to_share = nitem;

  return 1;
}

int PMMG_Set_ithFaceCommunicatorSize(PMMG_pParMesh parmesh, int ext_comm_index, int color_out, int nitem) {
  PMMG_pExt_comm pext_comm;

  /* Return if communicator index doesn't exist */
  if( (ext_comm_index < 0 ) || (ext_comm_index > parmesh->next_face_comm-1) ) {
    fprintf(stderr,"\n ## Error: function %s on proc %d: Communicator index should be in the range [0,ncomm-1], %d passed when %d communicators were asked. Please check your interface.\n",__func__,parmesh->myrank,ext_comm_index,parmesh->next_face_comm);
    return 0;
  }

  /* Return if communicator color doesn't exist or is myrank */
  if( (color_out < 0) || (color_out > parmesh->nprocs-1) || (color_out == parmesh->myrank) ) {
    fprintf(stderr,"\n ## Error: function %s on proc %d: Communicator color should be a proc ID in the range [0,nprocs-1] and different from my rank, %d passed while %d MPI procs are seen. Please check your interface.\n",__func__,parmesh->myrank,color_out,parmesh->nprocs);
    return 0;
  }

  /* Return if communicator size non-positive */
  if( nitem <= 0 ) {
    fprintf(stderr,"\n ## Error: function %s on proc %d: Communicator size should be strictly positive, %d passed. Please check your interface.\n",__func__,parmesh->myrank,nitem);
    return 0;
  }


  /* Get the face communicator */
  pext_comm = &parmesh->ext_face_comm[ext_comm_index];

  /* Set colors */
  pext_comm->color_in  = parmesh->myrank;
  pext_comm->color_out = color_out;

  /* Allocate communicator */
  PMMG_CALLOC(parmesh,pext_comm->int_comm_index,nitem,int,
                  "allocate int_comm_index",return 0);
  pext_comm->nitem = nitem;

  return 1;
}

int PMMG_Set_ithNodeCommunicator_nodes(PMMG_pParMesh parmesh, int ext_comm_index, int* local_index, int* global_index, int isNotOrdered) {
  PMMG_pExt_comm pext_node_comm;
  int            *oldId,nitem,i,ier;

  /* Return if communicator index doesn't exist */
  if( (ext_comm_index < 0 ) || (ext_comm_index > parmesh->next_node_comm-1) ) {
    fprintf(stderr,"\n ## Error: function %s on proc %d: Communicator index should be in the range [0,ncomm-1], %d passed when %d communicators were asked. Please check your interface.\n",__func__,parmesh->myrank,ext_comm_index,parmesh->next_node_comm);
    return 0;
  }

  /* Get the node communicator */
  pext_node_comm = &parmesh->ext_node_comm[ext_comm_index];
  nitem = pext_node_comm->nitem_to_share;

  /* Return if communicator is empty */
  if( !nitem ) {
    fprintf(stderr,"\n ## Error: function %s: 0 nodes in communicator %d on proc %d. Please check your interface.\n",__func__,ext_comm_index,parmesh->myrank);
    return 0;
  }

  /* Reorder according to global enumeration, so that the indexing matches on
   * the two sides of each pair of procs */
  ier = 1;
  if( isNotOrdered ) {
    PMMG_CALLOC(parmesh,oldId,nitem,int,"oldId",return 0);
    ier = PMMG_sort_iarray(parmesh,local_index,global_index,oldId,nitem);
    PMMG_DEL_MEM(parmesh,oldId,int,"oldId");
  }

  /* Save local and global node indices */
  if( ier ) {
    for( i = 0; i < nitem; i++ ) {
      pext_node_comm->int_comm_index[i] = local_index[i];
      pext_node_comm->itosend[i] = local_index[i];
      pext_node_comm->itorecv[i] = global_index[i];
    }
  }

  return ier;
}


int PMMG_Set_ithFaceCommunicator_faces(PMMG_pParMesh parmesh, int ext_comm_index, int* local_index, int* global_index, int isNotOrdered) {
  PMMG_pExt_comm pext_face_comm;
  int            *oldId,nitem,i,ier;

  /* Return if communicator index doesn't exist */
  if( (ext_comm_index < 0 ) || (ext_comm_index > parmesh->next_face_comm-1) ) {
    fprintf(stderr,"\n ## Error: function %s on proc %d: Communicator index should be in the range [0,ncomm-1], %d passed when %d communicators were asked. Please check your interface.\n",__func__,parmesh->myrank,ext_comm_index,parmesh->next_face_comm);
    return 0;
  }

  /* Get the face communicator */
  pext_face_comm = &parmesh->ext_face_comm[ext_comm_index];
  nitem = pext_face_comm->nitem;

  /* Return if communicator is empty */
  if( !nitem ) {
    fprintf(stderr,"\n ## Error: function %s: 0 faces in communicator %d on proc %d. Please check your interface.\n",__func__,ext_comm_index,parmesh->myrank);
    return 0;
  }

  /* Reorder according to global enumeration, so that faces indexing
   * matches on the two sides of each pair of procs */
  ier = 1;
  if( isNotOrdered ) {
    PMMG_CALLOC(parmesh,oldId,nitem,int,"oldId",return 0);
    ier = PMMG_sort_iarray(parmesh,local_index,global_index,oldId,nitem);
    PMMG_DEL_MEM(parmesh,oldId,int,"oldId");
  }

  /* Save local face index */
  for( i = 0; i < nitem; i++ )
    pext_face_comm->int_comm_index[i] = local_index[i];

  return ier;
}

int PMMG_Get_numberOfNodeCommunicators(PMMG_pParMesh parmesh, int *next_comm) {

  *next_comm = parmesh->next_node_comm;

  return 1;
}

int PMMG_Get_numberOfFaceCommunicators(PMMG_pParMesh parmesh, int *next_comm) {

  *next_comm = parmesh->next_face_comm;

  return 1;
}

int PMMG_Get_ithNodeCommunicatorSize(PMMG_pParMesh parmesh, int ext_comm_index, int *color_out, int *nitem) {
  PMMG_pExt_comm pext_node_comm;

  pext_node_comm = &parmesh->ext_node_comm[ext_comm_index];
  *color_out = pext_node_comm->color_out;
  *nitem     = pext_node_comm->nitem;

  return 1;
}

int PMMG_Get_ithFaceCommunicatorSize(PMMG_pParMesh parmesh, int ext_comm_index, int *color_out, int *nitem) {
  PMMG_pExt_comm pext_face_comm;

  pext_face_comm = &parmesh->ext_face_comm[ext_comm_index];
  *color_out = pext_face_comm->color_out;
  *nitem     = pext_face_comm->nitem;

  return 1;
}

int PMMG_Get_NodeCommunicator_nodes(PMMG_pParMesh parmesh, int** local_index) {
  PMMG_pGrp      grp;
  PMMG_pInt_comm int_node_comm;
  PMMG_pExt_comm ext_node_comm;
  MMG5_pMesh     mesh;
  int            ip,i,idx,icomm;

  /* Meshes are merged in grp 0 */
  int_node_comm = parmesh->int_node_comm;
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;


  /** 1) Store node index in intvalues */
  PMMG_CALLOC(parmesh,int_node_comm->intvalues,int_node_comm->nitem,int,"intvalues",return 0);
  for( i = 0; i < grp->nitem_int_node_comm; i++ ){
    ip   = grp->node2int_node_comm_index1[i];
    idx  = grp->node2int_node_comm_index2[i];
    parmesh->int_node_comm->intvalues[idx] = ip;
  }

  /** 2) For each external communicator, get node index from intvalues */
  for( icomm = 0; icomm < parmesh->next_node_comm; icomm++ ){
    ext_node_comm = &parmesh->ext_node_comm[icomm];
    for( i = 0; i < ext_node_comm->nitem; i++ ){
      idx = ext_node_comm->int_comm_index[i];
      local_index[icomm][i] = int_node_comm->intvalues[idx];
    }
  }

  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,int,"intvalues");
  return 1;
}

int PMMG_Get_FaceCommunicator_faces(PMMG_pParMesh parmesh, int** local_index) {
  MMG5_Hash      hash;
  PMMG_pGrp      grp;
  PMMG_pInt_comm int_face_comm;
  PMMG_pExt_comm ext_face_comm;
  MMG5_pMesh     mesh;
  MMG5_pTetra    pt;
  MMG5_pTria     ptt;
  int            kt,ie,ifac,ia,ib,ic,i,idx,icomm;

  /* Meshes are merged in grp 0 */
  int_face_comm = parmesh->int_face_comm;
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;


  /** 1) Hash triangles */
  if ( ! MMG5_hashNew(mesh,&hash,0.51*mesh->nt,1.51*mesh->nt) ) return 0;

  for (kt=1; kt<=mesh->nt; kt++) {
    ptt = &mesh->tria[kt];
    if ( !MMG5_hashFace(mesh,&hash,ptt->v[0],ptt->v[1],ptt->v[2],kt) ) {
      MMG5_DEL_MEM(mesh,hash.item);
      return 0;
    }
  }

  /** 2) Store triangle index in intvalues */
  PMMG_CALLOC(parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,"intvalues",return 0);
  for( i = 0; i < grp->nitem_int_face_comm; i++ ){
    ie   =  grp->face2int_face_comm_index1[i]/12;
    ifac = (grp->face2int_face_comm_index1[i]%12)/3;
    idx  =  grp->face2int_face_comm_index2[i];
    pt = &mesh->tetra[ie];
    ia = pt->v[MMG5_idir[ifac][0]];
    ib = pt->v[MMG5_idir[ifac][1]];
    ic = pt->v[MMG5_idir[ifac][2]];
    kt = MMG5_hashGetFace(&hash,ia,ib,ic);
    parmesh->int_face_comm->intvalues[idx] = kt;
  }

  /** 3) For each external communicator, get triangle index from intvalues */
  for( icomm = 0; icomm < parmesh->next_face_comm; icomm++ ){
    ext_face_comm = &parmesh->ext_face_comm[icomm];
    for( i = 0; i < ext_face_comm->nitem; i++ ){
      idx = ext_face_comm->int_comm_index[i];
      local_index[icomm][i] = int_face_comm->intvalues[idx];
    }
  }

  MMG5_DEL_MEM(mesh,hash.item);
  PMMG_DEL_MEM(parmesh,int_face_comm->intvalues,int,"intvalues");
  return 1;
}

int PMMG_Check_Set_NodeCommunicators(PMMG_pParMesh parmesh,int ncomm,int* nitem,
                                 int* color, int** local_index) {
  PMMG_pGrp      grp;
  PMMG_pInt_comm int_node_comm;
  PMMG_pExt_comm ext_node_comm;
  MMG5_pMesh     mesh;
  MMG5_Hash      hashPair;
  int            *values,*oldIdx,ip,idx,i,icomm,getComm;

  /* Meshes are merged in grp 0 */
  int_node_comm = parmesh->int_node_comm;
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;

  /** 1) Check number of communicators */
  if( parmesh->next_node_comm != ncomm ) {
    fprintf(stderr,"## Wrong number of node communicators on proc %d: input %d, set %d. ##\n",parmesh->myrank,ncomm,parmesh->next_node_comm);
    return 0;
  }

  /** 2) Create hash table for proc pairs */
  if( !MMG5_hashNew(mesh, &hashPair, 6*ncomm, 8*ncomm) ) return 0;
  for( icomm = 0; icomm < ncomm; icomm++ ) {
    ext_node_comm = &parmesh->ext_node_comm[icomm];
    /* Store IDs as IDs+1, so that 0 value can be used for error handling */
    if( !MMG5_hashEdge( mesh, &hashPair,
                        ext_node_comm->color_in+1, ext_node_comm->color_out+1,
                        icomm+1 ) ) {
      fprintf(stderr,"## Impossible to hash proc pair %d -- %d. ##\n",ext_node_comm->color_in,ext_node_comm->color_out);
      return 0;
    }
  }

  PMMG_CALLOC(parmesh,int_node_comm->intvalues,int_node_comm->nitem,int,"intvalues",return 0);
  PMMG_CALLOC(parmesh,values,int_node_comm->nitem,int,"values",return 0);
  PMMG_CALLOC(parmesh,oldIdx,int_node_comm->nitem,int,"oldIdx",return 0);

  /** 3) Put nodes index in intvalues */
  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    ip  = grp->node2int_node_comm_index1[i];
    idx = grp->node2int_node_comm_index2[i];
    int_node_comm->intvalues[idx] = ip;
   }


  /** 4) Check communicators size and color */
  for( icomm = 0; icomm < ncomm; icomm++ ) {
    getComm = MMG5_hashGet( &hashPair, parmesh->myrank+1, color[icomm]+1 );
    if( !getComm ) {
      fprintf(stderr,"## Interface %d --  %d not found!\n",
               parmesh->myrank,color[icomm]);
      return 0;
    }
    getComm--;
    ext_node_comm = &parmesh->ext_node_comm[getComm];
    if( color[icomm] != ext_node_comm->color_out ) {
      fprintf(stderr,"## Wrong color for node communicator %d on proc %d: input %d, set %d ##\n",icomm,parmesh->myrank,color[icomm],ext_node_comm->color_out);
      return 0;
    }
    if( nitem[icomm] != ext_node_comm->nitem ) {
      fprintf(stderr,"## Wrong size for node communicator %d on proc %d: input %d, set %d ##\n",icomm,parmesh->myrank,nitem[icomm],ext_node_comm->nitem);
      return 0;
    }
  }


  /** 5) Find input nodes in intvalues */
  for( icomm = 0; icomm < ncomm; icomm++ ) {
    getComm = MMG5_hashGet( &hashPair, parmesh->myrank+1, color[icomm]+1 );
    if( !getComm ) {
      fprintf(stderr,"## Interface %d --  %d not found!\n",
               parmesh->myrank,color[icomm]);
      return 0;
    }
    getComm--;
    ext_node_comm = &parmesh->ext_node_comm[getComm];
    
    /* Sort input data */
    PMMG_sort_iarray( parmesh,
                      values,local_index[icomm],
                      oldIdx, nitem[icomm] );
    
    /* Sort external communicator */
    for( i = 0; i < nitem[icomm]; i++ )
      values[i] = int_node_comm->intvalues[ext_node_comm->int_comm_index[i]];
    PMMG_sort_iarray( parmesh,
                      ext_node_comm->int_comm_index,values,
                      oldIdx, ext_node_comm->nitem );

    /* Check communicator against input data */
    for( i = 0; i < ext_node_comm->nitem; i++ ) {
      idx = ext_node_comm->int_comm_index[i];
      if( local_index[icomm][i] != int_node_comm->intvalues[idx] ) {
        fprintf(stderr,"## Impossible to find node %d in comm %d on proc %d. ##\n",local_index[icomm][i],icomm,parmesh->myrank);
        return 0;
      }
    }

    /* Ripristinate ext comm ordering */
    for( i = 0; i < ext_node_comm->nitem; i++ )
      values[i] = ext_node_comm->int_comm_index[i];
    for( i = 0; i < ext_node_comm->nitem; i++ )
      ext_node_comm->int_comm_index[oldIdx[i]] = values[i];
  }

  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,int,"intvalues");
  PMMG_DEL_MEM(parmesh,oldIdx,int,"oldIdx");
  PMMG_DEL_MEM(parmesh,values,int,"values");
  MMG5_DEL_MEM(mesh,hashPair.item);
  return 1;
}

int PMMG_Check_Set_FaceCommunicators(PMMG_pParMesh parmesh,int ncomm,int* nitem,
                                     int* color, int** trianodes) {
  PMMG_pGrp      grp;
  PMMG_pExt_comm ext_face_comm;
  MMG5_pMesh     mesh;
  MMG5_pTetra    pt;
  MMG5_Hash      hash;
  MMG5_Hash      hashPair;
  int            count,ie,ifac,ia,ib,ic,i,idx,icomm,getComm;

  /* Meshes are merged in grp 0 */
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;

  /** 1) Check number of communicators */
  if( parmesh->next_face_comm != ncomm ) {
    fprintf(stderr,"## Wrong number of face communicators in function %s\n",__func__);
    return 0;
  }

  /** 2) Create hash table for proc pairs */
  if( !MMG5_hashNew(mesh, &hashPair, 6*ncomm, 8*ncomm) ) {
    fprintf(stderr,"## Impossible to create hash table for process pairs in function %s\n",__func__);
    return 0;
  }
  for( icomm = 0; icomm < ncomm; icomm++ ) {
    ext_face_comm = &parmesh->ext_face_comm[icomm];
    /* Store IDs as IDs+1, so that 0 value can be used for error handling */
    if( !MMG5_hashEdge( mesh, &hashPair,
                        ext_face_comm->color_in+1, ext_face_comm->color_out+1,
                        icomm+1 ) ) {
      fprintf(stderr,"## Impossible to hash proc pair %d -- %d. ##\n",ext_face_comm->color_in,ext_face_comm->color_out);
      return 0;
    }
  }

  /** 3) Check communicators size and color */
  count = 0;
  for( icomm = 0; icomm < ncomm; icomm++ ) {
    getComm = MMG5_hashGet( &hashPair, parmesh->myrank+1, color[icomm]+1 );
    if( !getComm ) {
      fprintf(stderr,"## Interface %d --  %d not found!\n",
               parmesh->myrank,color[icomm]);
      return 0;
    }
    getComm--;
    ext_face_comm = &parmesh->ext_face_comm[getComm];
    if( color[icomm] != ext_face_comm->color_out ) {
      fprintf(stderr,"## Wrong color for face communicator %d on proc %d: input %d, set %d ##\n",icomm,parmesh->myrank,color[icomm],ext_face_comm->color_out);
      return 0;
    }
    if( nitem[icomm] != ext_face_comm->nitem ) {
      fprintf(stderr,"## Wrong size for face communicator %d on proc %d: input %d, set %d ##\n",icomm,parmesh->myrank,nitem[icomm],ext_face_comm->nitem);
      return 0;
    }
    count += nitem[icomm];
  }

  /** 4) Create boundary and hash table */
  if( MMG3D_bdryBuild(mesh) == -1) {
    fprintf(stderr,"## Impossible to build boundary in function %s\n",__func__);
    return 0;
  }
  if ( ! MMG5_hashNew(mesh,&hash,0.51*count,1.51*count) ) {
    fprintf(stderr,"## Impossible to build hash table for boundary faces in function %s\n",__func__);
    return 0;
  }

  /* Hash triangles in the internal communicator */
  for( i = 0; i < grp->nitem_int_face_comm; i++ ) {
    ie   =  grp->face2int_face_comm_index1[i]/12;
    ifac = (grp->face2int_face_comm_index1[i]%12)/3;
    idx  =  grp->face2int_face_comm_index2[i];
    pt = &mesh->tetra[ie];
    ia = pt->v[MMG5_idir[ifac][0]];
    ib = pt->v[MMG5_idir[ifac][1]];
    ic = pt->v[MMG5_idir[ifac][2]];
    /* Store ID+1 to use 0 value for error handling */
    if( !MMG5_hashFace(mesh,&hash,ia,ib,ic,idx+1) ) {
      fprintf(stderr,"## Impossible to hash face (%d,%d,%d) on proc %d. ##\n",ia,ib,ic,parmesh->myrank);
      MMG5_DEL_MEM(mesh,hash.item);
      return 0;
    }
  }

  /** 5) Find input triangles in the hash table */
  for( icomm = 0; icomm < ncomm; icomm++ ){
    for( i = 0; i < nitem[icomm]; i++ ) {
      ia = trianodes[icomm][3*i];
      ib = trianodes[icomm][3*i+1];
      ic = trianodes[icomm][3*i+2];
      if( !MMG5_hashGetFace(&hash,ia,ib,ic) ) {
        fprintf(stderr,"## Face (%d,%d,%d) not found in face communicator %d on proc %d. ##\n",ia,ib,ic,icomm,parmesh->myrank);
        MMG5_DEL_MEM(mesh,hash.item);
        return 0;
      }
    }
  }

  MMG5_DEL_MEM(mesh,hash.item);
  MMG5_DEL_MEM(mesh,hashPair.item);
  return 1;
}

int PMMG_Check_Get_NodeCommunicators(PMMG_pParMesh parmesh,
                                     int ncomm_in,int* nitem_in,
                                     int* color_in, int** local_index_in,
                                     int ncomm_out,int* nitem_out,
                                     int* color_out, int** local_index_out) {
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_Hash      hashPair;
  int            *values,*oldIdx,i,icomm,getComm,count;

  /* Meshes are merged in grp 0 */
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;

  /** 1) Check number of communicators */
  if( ncomm_in != ncomm_out) {
    fprintf(stderr,"## Wrong number of node communicators on proc %d: input %d, set %d. ##\n",parmesh->myrank,ncomm_in,ncomm_out);
    return 0;
  }

  /** 2) Create hash table for proc pairs */
  if( !MMG5_hashNew(mesh, &hashPair, 6*ncomm_in, 8*ncomm_in) ) return 0;
  for( icomm = 0; icomm < ncomm_in; icomm++ ) {
    /* Store IDs as IDs+1, so that 0 value can be used for error handling */
    if( !MMG5_hashEdge( mesh, &hashPair,
                        parmesh->myrank+1, color_out[icomm]+1,
                        icomm+1 ) ) return 0;
  }


  /** 3) Check communicators size and color */
  count = 0;
  for( icomm = 0; icomm < ncomm_in; icomm++ ) {
    getComm = MMG5_hashGet( &hashPair, parmesh->myrank+1, color_in[icomm]+1 );
    if( !getComm ) {
      fprintf(stderr,"## Interface %d --  %d not found!\n",
               parmesh->myrank,color_in[icomm]);
      return 0;
    }
    getComm--;
    if( color_in[icomm] != color_out[getComm] ) {
      fprintf(stderr,"## Wrong color for node communicator %d on proc %d: input %d, set %d ##\n",icomm,parmesh->myrank,color_in[icomm],color_out[getComm]);
      return 0;
    }
    if( nitem_in[icomm] != nitem_out[getComm] ) {
      fprintf(stderr,"## Wrong size for node communicator %d on proc %d: input %d, set %d ##\n",icomm,parmesh->myrank,nitem_in[icomm],nitem_out[getComm]);
      return 0;
    }
    if( nitem_in[icomm] > count ) count = nitem_in[icomm];
  }

  PMMG_CALLOC(parmesh,values,count,int,"values",return 0);
  PMMG_CALLOC(parmesh,oldIdx,count,int,"oldIdx",return 0);


  /** 4) Find input nodes */
  for( icomm = 0; icomm < ncomm_in; icomm++ ) {
    getComm = MMG5_hashGet( &hashPair, parmesh->myrank+1, color_in[icomm]+1 );
    if( !getComm ) {
      fprintf(stderr,"## Interface %d --  %d not found!\n",
               parmesh->myrank,color_in[icomm]);
      return 0;
    }
    getComm--;

    /* Sort input data */
    PMMG_sort_iarray( parmesh,
                      values,local_index_in[icomm],
                      oldIdx, nitem_in[icomm] );

    /* Sort external communicator */
    PMMG_sort_iarray( parmesh,
                      values, local_index_out[getComm],
                      oldIdx, nitem_out[getComm] );

    /* Check communicator against input data */
    for( i = 0; i < nitem_in[icomm]; i++ ) {
      if( local_index_in[icomm][i] != local_index_out[getComm][i] ) return 0;
    }

    /* Ripristinate ext comm ordering */
    for( i = 0; i < nitem_in[icomm]; i++ )
      values[i] = local_index_out[getComm][i];
    for( i = 0; i < nitem_in[icomm]; i++ )
      local_index_out[getComm][oldIdx[i]] = values[i];
  }

  PMMG_DEL_MEM(parmesh,oldIdx,int,"oldIdx");
  PMMG_DEL_MEM(parmesh,values,int,"values");
  MMG5_DEL_MEM(mesh,hashPair.item);
 return 1;
}

int PMMG_Check_Get_FaceCommunicators(PMMG_pParMesh parmesh,
                                     int ncomm_in,int* nitem_in,
                                     int* color_in, int** trianodes_in,
                                     int ncomm_out,int* nitem_out,
                                     int* color_out, int** trianodes_out) {
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_Hash      hash;
  MMG5_Hash      hashPair;
  int            count,ia,ib,ic,i,icomm,getComm;

  /* Meshes are merged in grp 0 */
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;

  /** 1) Check number of communicators */
  if( ncomm_in != ncomm_out ) return 0;

  /** 2) Create hash table for proc pairs */
  if( !MMG5_hashNew(mesh, &hashPair, 6*ncomm_in, 8*ncomm_in) ) return 0;
  for( icomm = 0; icomm < ncomm_in; icomm++ ) {
    /* Store IDs as IDs+1, so that 0 value can be used for error handling */
    if( !MMG5_hashEdge( mesh, &hashPair,
                        parmesh->myrank+1, color_out[icomm]+1,
                        icomm+1 ) ) {
      fprintf(stderr,"## Impossible to hash proc pair %d -- %d. ##\n",parmesh->myrank,color_out[icomm]);
      return 0;
    }
  }

  /** 3) Check communicators size and color */
  count = 0;
  for( icomm = 0; icomm < ncomm_in; icomm++ ) {
    getComm = MMG5_hashGet( &hashPair, parmesh->myrank+1, color_in[icomm]+1 );
    if( !getComm ) {
      fprintf(stderr,"## Interface %d --  %d not found!\n",
               parmesh->myrank,color_in[icomm]);
      return 0;
    }
    getComm--;
    if( color_in[icomm] != color_out[getComm] ) {
      fprintf(stderr,"## Wrong color for face communicator %d on proc %d: input %d, set %d ##\n",icomm,parmesh->myrank,color_in[icomm],color_out[getComm]);
      return 0;
    }
    if( nitem_in[icomm] != nitem_out[getComm] ) {
      fprintf(stderr,"## Wrong size for face communicator %d on proc %d: input %d, set %d ##\n",icomm,parmesh->myrank,nitem_in[icomm],nitem_out[getComm]);
      return 0;
    }
    count += nitem_in[icomm];
  }

  /** 4) Create boundary and hash table */
  if( MMG3D_bdryBuild(mesh) == -1) return 0;
  if ( ! MMG5_hashNew(mesh,&hash,0.51*count,1.51*count) ) return 0;

  /* Hash triangles in the internal communicator */
  count = 0;
  for( icomm = 0; icomm <ncomm_in; icomm++) {
    for( i = 0; i < nitem_in[icomm]; i++ ) {
      ia = trianodes_out[icomm][3*i];
      ib = trianodes_out[icomm][3*i+1];
      ic = trianodes_out[icomm][3*i+2];
      /* Store ID+1 to use 0 value for error handling */
      if( !MMG5_hashFace(mesh,&hash,ia,ib,ic,++count) ) {
        fprintf(stderr,"## Impossible to hash face (%d,%d,%d) on proc %d. ##\n",ia,ib,ic,parmesh->myrank);
         MMG5_DEL_MEM(mesh,hash.item);
        return 0;
      }
    }
  }

  /** 5) Find input triangles in the hash table */
  for( icomm = 0; icomm < ncomm_in; icomm++ ){
    for( i = 0; i < nitem_in[icomm]; i++ ) {
      ia = trianodes_in[icomm][3*i];
      ib = trianodes_in[icomm][3*i+1];
      ic = trianodes_in[icomm][3*i+2];
      if( !MMG5_hashGetFace(&hash,ia,ib,ic) ) {
        fprintf(stderr,"## Face (%d,%d,%d) not found in face communicator %d on proc %d. ##\n",ia,ib,ic,icomm,parmesh->myrank);
        MMG5_DEL_MEM(mesh,hash.item);
        return 0;
      }
    }
  }

  MMG5_DEL_MEM(mesh,hash.item);
  MMG5_DEL_MEM(mesh,hashPair.item);
  return 1;
}

/**
 * \param parmesh pointer toward parmesh structure.
 * \param idx_glob pointer to the global triangle numbering.
 * \param owner pointer to the rank of the process owning the triangle.
 * \return 1 if success, 0 if fail.
 *
 * Get global triangle numbering (starting from 1) and rank of the process owning
 * the triangle.
 * If of the triangle is simply a parallel face (but not a boundary), its owner
 * will be negative.
 */
int PMMG_Get_triangleGloNum( PMMG_pParMesh parmesh, int *idx_glob, int *owner ) {
  MMG5_pMesh mesh;
  MMG5_pTria ptr;

  if( !parmesh->info.globalNum ) {
    fprintf(stderr,"\n  ## Error: %s: Triangle global numbering has not been computed.\n",
            __func__);
    fprintf(stderr,"     Parameter PMMG_IPARAM_globalNum has to be set to 1.\n");
    fprintf(stderr,"     Please rerun ParMmg).\n");
    return 0;
  }

  assert( parmesh->ngrp == 1 );
  mesh = parmesh->listgrp[0].mesh;

  if ( mesh->nti == mesh->nt ) {
    mesh->nti = 0;
    if( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of triangles.\n",
              __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the PMMG_Get_triangleGloNum function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of triangles: %d\n ",mesh->nt);
    }
  }

  mesh->nti++;

  if ( mesh->nti > mesh->nt ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get numbering.\n",__func__);
    fprintf(stderr,"     The number of call of PMMG_Get_triangleGloNum function");
    fprintf(stderr," can not exceed the number of triangles: %d\n ",mesh->nt);
    return 0;
  }

  ptr = &mesh->tria[mesh->nti];
  *idx_glob = ptr->flag;
  *owner    = ptr->base;

  return 1;
}

/**
 * \param parmesh pointer toward parmesh structure.
 * \param idx_glob array of global triangles numbering.
 * \param owner array of ranks of processes owning each triangle.
 * \return 1 if success, 0 if fail.
 *
 * Get global triangles numbering (starting from 1) and ranks of processes owning
 * each node.
 * If of the triangle is simply a parallel face (but not a boundary), its owner
 * will be negative.
 */
int PMMG_Get_trianglesGloNum( PMMG_pParMesh parmesh, int *idx_glob, int *owner ) {
  MMG5_pMesh mesh;
  MMG5_pTria ptr;
  int        k;

  if( !parmesh->info.globalNum ) {
    fprintf(stderr,"\n  ## Error: %s: Triangles global numbering has not been computed.\n",
            __func__);
    fprintf(stderr,"     Parameter PMMG_IPARAM_globalNum has to be set to 1.\n");
    fprintf(stderr,"     Please rerun ParMmg).\n");
    return 0;
  }

  assert( parmesh->ngrp == 1 );
  mesh = parmesh->listgrp[0].mesh;

  for( k = 1; k <= mesh->nt; k++ ){
    ptr = &mesh->tria[k];
    idx_glob[k-1] = ptr->flag;
    owner[k-1]    = ptr->base;
  }

  return 1;
}

/**
 * \param parmesh pointer toward parmesh structure.
 * \param idx_glob pointer to the global node numbering.
 * \param owner pointer to the rank of the process owning the node.
 * \return 1 if success, 0 if fail.
 *
 * Get global node numbering (starting from 1) and rank of the process owning
 * the node.
 *
 */
int PMMG_Get_vertexGloNum( PMMG_pParMesh parmesh, int *idx_glob, int *owner ) {
  MMG5_pMesh  mesh;
  MMG5_pPoint ppt;

  if( !parmesh->info.globalNum ) {
    fprintf(stderr,"\n  ## Error: %s: Nodes global numbering has not been computed.\n",
            __func__);
    fprintf(stderr,"     Parameter PMMG_IPARAM_globalNum has to be set to 1.\n");
    fprintf(stderr,"     Please rerun ParMmg).\n");
    return 0;
  }

  assert( parmesh->ngrp == 1 );
  mesh = parmesh->listgrp[0].mesh;

  if ( mesh->npi == mesh->np ) {
    mesh->npi = 0;
    if( mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Warning: %s: reset the internal counter of points.\n",
              __func__);
      fprintf(stderr,"     You must pass here exactly one time (the first time ");
      fprintf(stderr,"you call the PMMG_Get_vertexGloNum function).\n");
      fprintf(stderr,"     If not, the number of call of this function");
      fprintf(stderr," exceed the number of points: %d\n ",mesh->np);
    }
  }

  mesh->npi++;

  if ( mesh->npi > mesh->np ) {
    fprintf(stderr,"\n  ## Error: %s: unable to get numbering.\n",__func__);
    fprintf(stderr,"     The number of call of PMMG_Get_vertexGloNum function");
    fprintf(stderr," can not exceed the number of points: %d\n ",mesh->np);
    return 0;
  }

  ppt = &mesh->point[mesh->npi];
  *idx_glob = ppt->tmp;
  *owner    = ppt->flag;

  return 1;
}

/**
 * \param parmesh pointer toward parmesh structure.
 * \param idx_glob array of global nodes numbering.
 * \param owner array of ranks of processes owning each node.
 * \return 1 if success, 0 if fail.
 *
 * Get global nodes numbering (starting from 1) and ranks of processes owning
 * each node.
 *
 */
int PMMG_Get_verticesGloNum( PMMG_pParMesh parmesh, int *idx_glob, int *owner ) {
  MMG5_pMesh  mesh;
  MMG5_pPoint ppt;
  int         ip;

  if( !parmesh->info.globalNum ) {
    fprintf(stderr,"\n  ## Error: %s: Nodes global numbering has not been computed.\n",
            __func__);
    fprintf(stderr,"     Parameter PMMG_IPARAM_globalNum has to be set to 1.\n");
    fprintf(stderr,"     Please rerun ParMmg).\n");
    return 0;
  }

  assert( parmesh->ngrp == 1 );
  mesh = parmesh->listgrp[0].mesh;

  for( ip = 1; ip <= mesh->np; ip++ ){
    ppt = &mesh->point[ip];
    idx_glob[ip-1] = ppt->tmp;
    owner[ip-1]    = ppt->flag;
  }

  return 1;
}

/**
 * \param parmesh pointer toward parmesh structure
 * \param owner IDs of the processes owning each interface node
 * \param idx_glob global IDs of interface nodes
 * \param nunique nb of non-redundant interface nodes on current rank
 * \param ntot totat nb of non-redundant interface nodes
 *
 * Create global IDs (starting from 1) for nodes on parallel interfaces.
 *
 */
int PMMG_Get_NodeCommunicator_owners(PMMG_pParMesh parmesh,int **owner,int **idx_glob,int *nunique, int *ntot) {
  PMMG_pInt_comm int_node_comm;
  PMMG_pExt_comm ext_node_comm;
  PMMG_pGrp      grp;
  MPI_Request    request;
  MPI_Status     status;
  int            *intvalues,*itosend,*itorecv,*iproc2comm;
  int            color,nitem;
  int            label,*nlabels,*displ,mydispl,unique;
  int            icomm,i,idx,iproc,src,dst,tag;

  /* Do this only if there is one group */
  assert( parmesh->ngrp == 1 );
  grp = &parmesh->listgrp[0];

  /* Allocate internal communicator */
  int_node_comm = parmesh->int_node_comm;
  PMMG_CALLOC(parmesh,int_node_comm->intvalues,int_node_comm->nitem,int,"intvalues",return 0);
  intvalues = int_node_comm->intvalues;

  /* Allocate label counts and offsets */
  PMMG_CALLOC(parmesh,nlabels,parmesh->nprocs,int,"nlabels",return 0);
  PMMG_CALLOC(parmesh,displ,parmesh->nprocs+1,int,"displ",return 0);

  /* Array to reorder communicators */
  PMMG_MALLOC(parmesh,iproc2comm,parmesh->nprocs,int,"iproc2comm",return 0);

  for( iproc = 0; iproc < parmesh->nprocs; iproc++ )
    iproc2comm[iproc] = PMMG_UNSET;

  for( icomm = 0; icomm < parmesh->next_node_comm; icomm++ ) {
    ext_node_comm = &parmesh->ext_node_comm[icomm];
    iproc = ext_node_comm->color_out;
    iproc2comm[iproc] = icomm;
  }

  /**
   * 1) Number and count. Analyse each communicator by external color order,
   *     fill internal communicator with (PMMG_UNSET-color) if the node is not
   *     owned by myrank, or count the node if it is owned by myrank.
   *     (With this ordering, any not-owned node has been necessarily visited
   *     by its owner color and cannot be counted)
   */
  label = 0;
  unique = 0;
  for( color = 0; color < parmesh->myrank; color++ ) {
    icomm = iproc2comm[color];

    /* Skip non-existent communicators */
    if( icomm == PMMG_UNSET ) continue;

    ext_node_comm = &parmesh->ext_node_comm[icomm];
    nitem =  ext_node_comm->nitem;

    /* Mark not-owned nodes */
    for( i = 0; i < nitem; i++ ) {
      idx = ext_node_comm->int_comm_index[i];
      /* Only the first visitor owns the ghost node */
      if( intvalues[idx] ) continue;
      intvalues[idx] = PMMG_UNSET-color;
      ++unique;
    }
  }
  for( color = parmesh->myrank+1; color < parmesh->nprocs; color++ ) {
    icomm = iproc2comm[color];

    /* Skip non-existent communicators */
    if( icomm == PMMG_UNSET ) continue;

    ext_node_comm = &parmesh->ext_node_comm[icomm];
    nitem =  ext_node_comm->nitem;

   /* Count points only on owned communicators */
    for( i = 0; i < nitem; i++ ) {
      idx = ext_node_comm->int_comm_index[i];
      /* Count point only if not already marked */
      if( intvalues[idx] ) continue;
      intvalues[idx] = ++label;
      ++unique;
    }
  }


  if ( owner ) {
    /**
     * 2)  Store owners in the output array. Not-owned nodes store a
     *      (PMMG_UNSET-color) label in the internal communicator, while nodes
     *      owned by myrank store a non-negative label.
     */
    for( icomm = 0; icomm < parmesh->next_node_comm; icomm++ ) {
      ext_node_comm = &parmesh->ext_node_comm[icomm];
      color = ext_node_comm->color_out;
      nitem = ext_node_comm->nitem;

      for( i = 0; i < nitem; i++ ) {
        idx = ext_node_comm->int_comm_index[i];
        if( intvalues[idx] < 0 )
          owner[icomm][i] = -(intvalues[idx]-PMMG_UNSET);
        else
          owner[icomm][i] = parmesh->myrank;
      }
    }
  }


  /**
   * 3) Compute a consecutive global numbering by retrieving parallel offsets
   */

  /* Get nb of labels on each proc and compute offsets */
  MPI_Allgather( &label,1,MPI_INT,
                 nlabels,1,MPI_INT,parmesh->comm );

  for( iproc = 0; iproc < parmesh->nprocs; iproc++ )
    displ[iproc+1] = displ[iproc]+nlabels[iproc];
  mydispl = displ[parmesh->myrank];

  /* Get nb of non-redundant entities on each proci and total (for output) */
  if( nunique ) *nunique = unique;
  if( ntot )    *ntot = displ[parmesh->nprocs];


  /* Add offset to the owned labels */
  for( idx = 0; idx < grp->nitem_int_node_comm; idx++ ) {
    if( intvalues[idx] <= PMMG_UNSET ) continue;
    intvalues[idx] += mydispl;
  }


  /**
   * 4) Communicate global numbering to the ghost copies.
   */
  for( icomm = 0; icomm < parmesh->next_node_comm; icomm++ ) {
    ext_node_comm = &parmesh->ext_node_comm[icomm];
    color = ext_node_comm->color_out;
    nitem = ext_node_comm->nitem;

    PMMG_CALLOC(parmesh,ext_node_comm->itosend,nitem,int,"itosend",return 0);
    PMMG_CALLOC(parmesh,ext_node_comm->itorecv,nitem,int,"itorecv",return 0);
    itosend = ext_node_comm->itosend;
    itorecv = ext_node_comm->itorecv;

    src = MG_MIN(parmesh->myrank,color);
    dst = MG_MAX(parmesh->myrank,color);
    tag = parmesh->nprocs*src+dst;

    if( parmesh->myrank == src ) {
      /* Fill send buffer from internal communicator */
      for( i = 0; i < nitem; i++ ) {
        idx = ext_node_comm->int_comm_index[i];
        itosend[i] = intvalues[idx];
      }
      MPI_CHECK( MPI_Isend(itosend,nitem,MPI_INT,dst,tag,
                           parmesh->comm,&request),return 0 );
    }
    if ( parmesh->myrank == dst ) {
      MPI_CHECK( MPI_Recv(itorecv,nitem,MPI_INT,src,tag,
                          parmesh->comm,&status),return 0 );
      /* Store recv buffer in the internal communicator */
      for( i = 0; i < nitem; i++ ) {
        idx = ext_node_comm->int_comm_index[i];
        /* Update the value only if receiving it from the owner, or if already
         *  updated by the sender */
        if( itorecv[i] > PMMG_UNSET ) intvalues[idx] = itorecv[i];
      }
    }
  }


  /**
   * 5) Store numbering results in the output array.
   */
  for( icomm = 0; icomm < parmesh->next_node_comm; icomm++ ) {
    ext_node_comm = &parmesh->ext_node_comm[icomm];
    color = ext_node_comm->color_out;
    nitem = ext_node_comm->nitem;

    for( i = 0; i < nitem; i++ ) {
      idx = ext_node_comm->int_comm_index[i];
      idx_glob[icomm][i] = intvalues[idx];
    }
  }


#ifndef NDEBUG
  /* Check global IDs */
  int *mylabels;
  PMMG_CALLOC(parmesh,mylabels,label+1,int,"mylabels",return 0);

  /* Purposely in reverse order to overwrite internal communicator */
  for( iproc = parmesh->nprocs-1; iproc >= 0; iproc-- ) {
    icomm = iproc2comm[iproc];
    if( icomm == PMMG_UNSET ) continue;

    ext_node_comm = &parmesh->ext_node_comm[icomm];
    color = ext_node_comm->color_out;
    nitem = ext_node_comm->nitem;

    itorecv = ext_node_comm->itorecv;

    src = MG_MIN(parmesh->myrank,color);
    dst = MG_MAX(parmesh->myrank,color);
    tag = parmesh->nprocs*src+dst;
    if( parmesh->myrank == src ) {
      MPI_CHECK( MPI_Isend(idx_glob[icomm],nitem,MPI_INT,dst,tag,
                            parmesh->comm,&request),return 0 );
    }
    if ( parmesh->myrank == dst ) {
      MPI_CHECK( MPI_Recv(itorecv,nitem,MPI_INT,src,tag,
                          parmesh->comm,&status),return 0 );
      for( i=0; i < nitem; i++ ) {
        idx = ext_node_comm->int_comm_index[i];
        assert( idx_glob[icomm][i] == intvalues[idx] );
        assert( idx_glob[icomm][i] == itorecv[i] );
      }
    }

    /* Mark seen labels */
    if( parmesh->myrank < color ) {
      for( i=0; i < nitem; i++ ) {
        idx = ext_node_comm->int_comm_index[i];
        if( intvalues[idx] <= mydispl ) continue;
        mylabels[intvalues[idx]-mydispl]++;
      }
    }
  }
  /* Check for holes in the seen labels */
  for( i = 1; i <= label; i++ )
    assert(mylabels[i]);

  PMMG_DEL_MEM(parmesh,mylabels,int,"mylabels");
#endif

  /* Don't free buffers before they have been received */
  MPI_CHECK( MPI_Barrier(parmesh->comm),return 0 );

  /* Free arrays */
  PMMG_DEL_MEM(parmesh,nlabels,int,"nlabels");
  PMMG_DEL_MEM(parmesh,displ,int,"displ");
  PMMG_DEL_MEM(parmesh,iproc2comm,int,"iproc2comm");

  for( icomm = 0; icomm < parmesh->next_node_comm; icomm++ ) {
    ext_node_comm = &parmesh->ext_node_comm[icomm];
    PMMG_DEL_MEM(parmesh,ext_node_comm->itosend,int,"itosend");
    PMMG_DEL_MEM(parmesh,ext_node_comm->itorecv,int,"itorecv");
  }

  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,int,"intvalues");

  return 1;
}

/**
 * \param parmesh pointer toward parmesh structure
 * \param owner IDs of processes owning each interface triangle
 * \param idx_glob global IDs of interface triangles
 * \param nunique nb of non-redundant interface triangles on current rank
 * \param ntot totat nb of non-redundant interface triangles
 *
 * Create global IDs (starting from 1) for triangles on parallel interfaces.
 *
 */
int PMMG_Get_FaceCommunicator_owners(PMMG_pParMesh parmesh,int **owner,int **idx_glob,int *nunique,int *ntot) {
  PMMG_pExt_comm ext_face_comm;
  MPI_Request    request;
  MPI_Status     status;
  int            unique;
  int            color,nitem,npairs_loc,*npairs,*displ_pair,*glob_pair_displ;
  int            src,dst,tag,sendbuffer,recvbuffer,iproc,icomm,i;

  /* Do this only if there is one group */
  assert( parmesh->ngrp == 1 );

  PMMG_CALLOC(parmesh,npairs,parmesh->nprocs,int,"npair",return 0);
  PMMG_CALLOC(parmesh,displ_pair,parmesh->nprocs+1,int,"displ_pair",return 0);


  /**
   * 1) Compute face owners and count nb of new pair faces hosted on myrank.
   */
  npairs_loc = 0;
  unique = 0;
  for( icomm = 0; icomm < parmesh->next_face_comm; icomm++ ) {
    ext_face_comm = &parmesh->ext_face_comm[icomm];
    color = ext_face_comm->color_out;
    nitem = ext_face_comm->nitem;
    unique += nitem;
    if( color > parmesh->myrank ) npairs_loc += nitem;//1;

    if ( owner ) {
      for( i = 0; i < nitem; i++ ) {
        owner[icomm][i] = MG_MIN(color,parmesh->myrank);
      }
    }
  }


  /**
   * 2) Compute global face numbering. Communicate parallel offsets on each
   *    communicator, than each process update the numbering independently.
   */

  /* Get nb of pair faces and compute pair offset */
  MPI_Allgather( &npairs_loc,1,MPI_INT,
                 npairs,1,MPI_INT,parmesh->comm );

  for( iproc = 0; iproc < parmesh->nprocs; iproc++ )
    displ_pair[iproc+1] = displ_pair[iproc]+npairs[iproc];

  PMMG_CALLOC(parmesh,glob_pair_displ,parmesh->next_face_comm+1,int,"glob_pair_displ",return 0);
  for( icomm = 0; icomm < parmesh->next_face_comm; icomm++ )
    glob_pair_displ[icomm] = displ_pair[parmesh->myrank];

  /* Store nb of non-redundant faces on each proc and in total for output */
  if( nunique ) *nunique = unique;
  if( ntot ) *ntot = displ_pair[parmesh->nprocs];


  for( icomm = 0; icomm < parmesh->next_face_comm; icomm++ ) {
    ext_face_comm = &parmesh->ext_face_comm[icomm];
    color = ext_face_comm->color_out;
    nitem = ext_face_comm->nitem;

    if( color > parmesh->myrank )
      glob_pair_displ[icomm+1] = glob_pair_displ[icomm]+nitem;//+1;
  }

  /* Compute global pair faces enumeration */
  for( icomm = 0; icomm < parmesh->next_face_comm; icomm++ ) {
    ext_face_comm = &parmesh->ext_face_comm[icomm];
    color = ext_face_comm->color_out;
    nitem = ext_face_comm->nitem;

    /* Assign global index */
    src = MG_MIN(parmesh->myrank,color);
    dst = MG_MAX(parmesh->myrank,color);
    tag = parmesh->nprocs*src+dst;
    if( parmesh->myrank == src ) {
      sendbuffer = glob_pair_displ[icomm];
      MPI_CHECK( MPI_Isend(&sendbuffer,1,MPI_INT,dst,tag,
                            parmesh->comm,&request),return 0 );
    }
    if ( parmesh->myrank == dst ) {
      MPI_CHECK( MPI_Recv(&recvbuffer,1,MPI_INT,src,tag,
                          parmesh->comm,&status),return 0 );
      glob_pair_displ[icomm] = recvbuffer;
    }
  }

  for( icomm = 0; icomm < parmesh->next_face_comm; icomm++ ) {
    ext_face_comm = &parmesh->ext_face_comm[icomm];
    color = ext_face_comm->color_out;
    nitem = ext_face_comm->nitem;

    for( i = 0; i < nitem; i++ )
      idx_glob[icomm][i] = glob_pair_displ[icomm]+i+1; /* index starts from 1 */
  }

  /* Don't free buffers before they have been received */
  MPI_CHECK( MPI_Barrier(parmesh->comm),return 0 );

  /* Free arrays */
  PMMG_DEL_MEM(parmesh,npairs,int,"npairs");
  PMMG_DEL_MEM(parmesh,displ_pair,int,"displ_pair");
  PMMG_DEL_MEM(parmesh,glob_pair_displ,int,"glob_pair_displ");


#ifndef NDEBUG
  /* Check global IDs */
  for( icomm = 0; icomm < parmesh->next_face_comm; icomm++ ) {
    ext_face_comm = &parmesh->ext_face_comm[icomm];
    color = ext_face_comm->color_out;
    nitem = ext_face_comm->nitem;
    PMMG_CALLOC(parmesh,ext_face_comm->itorecv,nitem,int,"itorecv",return 0);

    src = MG_MIN(parmesh->myrank,color);
    dst = MG_MAX(parmesh->myrank,color);
    tag = parmesh->nprocs*src+dst;
    if( parmesh->myrank == src ) {
      MPI_CHECK( MPI_Isend(idx_glob[icomm],nitem,MPI_INT,dst,tag,
                            parmesh->comm,&request),return 0 );
    }
    if ( parmesh->myrank == dst ) {
      MPI_CHECK( MPI_Recv(ext_face_comm->itorecv,nitem,MPI_INT,src,tag,
                          parmesh->comm,&status),return 0 );
      for( i = 0; i < nitem; i++ )
        assert( idx_glob[icomm][i] == ext_face_comm->itorecv[i] );
    }
  }

  for( icomm = 0; icomm < parmesh->next_face_comm; icomm++ ) {
    ext_face_comm = &parmesh->ext_face_comm[icomm];
    PMMG_DEL_MEM(parmesh,ext_face_comm->itorecv,int,"itorecv");
  }
#endif

  return 1;
}

int PMMG_Free_names(PMMG_pParMesh parmesh)
{
  PMMG_DEL_MEM ( parmesh, parmesh->meshin,char,"meshin" );
  PMMG_DEL_MEM ( parmesh, parmesh->meshout,char,"meshout" );
  PMMG_DEL_MEM ( parmesh, parmesh->metin,char,"metin" );
  PMMG_DEL_MEM ( parmesh, parmesh->metout,char,"metout" );
  PMMG_DEL_MEM ( parmesh, parmesh->lsin,char,"lsin" );
  PMMG_DEL_MEM ( parmesh, parmesh->dispin,char,"dispin" );
  PMMG_DEL_MEM ( parmesh, parmesh->fieldin,char,"fieldin" );
  PMMG_DEL_MEM ( parmesh, parmesh->fieldout,char,"fieldout" );
  return 1;
}

int PMMG_Free_all(const int starter,...)
{
  va_list argptr;
  int     ier;

  va_start(argptr, starter);

  ier = PMMG_Free_all_var(argptr);

  va_end(argptr);

  return ier;
}

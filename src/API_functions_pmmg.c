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

int PMMG_Set_inputMeshName(PMMG_pParMesh parmesh, const char* meshin) {
  MMG5_pMesh mesh;
  int        k,ier;

  ier = 1;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    ier  = MG_MIN( ier, MMG3D_Set_inputMeshName(mesh,meshin) );
  }
  return ier;
}

int PMMG_Set_inputSolsName(PMMG_pParMesh parmesh, const char* solin) {
  MMG5_pMesh mesh;
  MMG5_pSol  sol,psl;
  int        k,i,ier;

  ier = 1;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    sol  = parmesh->listgrp[k].sol;

    for ( i=0; i<mesh->nsols; ++i ) {
      psl = sol + i;
      ier  = MG_MIN( ier, MMG3D_Set_inputSolName(mesh,psl,solin) );
    }
  }
  return ier;
}

int PMMG_Set_inputMetName(PMMG_pParMesh parmesh, const char* metin) {
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  int        k,ier;

  ier = 1;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    met  = parmesh->listgrp[k].met;
    ier  = MG_MIN( ier, MMG3D_Set_inputSolName(mesh,met,metin) );
  }
  return ier;
}

int PMMG_Set_outputMeshName(PMMG_pParMesh parmesh, const char* meshout) {
  MMG5_pMesh mesh;
  int        k,ier;

  ier = 1;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    ier  = MG_MIN( ier, MMG3D_Set_outputMeshName(mesh,meshout) );
  }
  return ier;
}

int PMMG_Set_outputSolsName(PMMG_pParMesh parmesh, const char* solout) {
  MMG5_pMesh mesh;
  MMG5_pSol  sol,psl;
  int        k,i,ier;

  ier = 1;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    sol  = parmesh->listgrp[k].sol;
    for ( i=0; i<mesh->nsols; ++i ) {
      psl = sol + i;
      ier  = MG_MIN( ier, MMG3D_Set_outputSolName(mesh,psl,solout) );
    }
  }
  return ier;
}

int PMMG_Set_outputMetName(PMMG_pParMesh parmesh, const char* metout) {
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  int        k,ier;

  ier = 1;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    met  = parmesh->listgrp[k].met;
    ier  = MG_MIN( ier, MMG3D_Set_outputSolName(mesh,met,metout) );
  }
  return ier;
}

void PMMG_Init_parameters(PMMG_pParMesh parmesh,MPI_Comm comm) {
  MMG5_pMesh mesh;
  size_t     mem;
  int        k,flag;

  memset(&parmesh->info,0, sizeof(PMMG_Info));

  parmesh->info.mem    = PMMG_UNSET; /* [n/-1]   ,Set memory size to n Mbytes/keep the default value */
  parmesh->info.root   = PMMG_NUL;

  parmesh->ddebug      = PMMG_NUL;
  parmesh->niter       = PMMG_NITER;
  parmesh->info.fem    = MMG5_FEM;
  parmesh->info.repartitioning = PMMG_REDISTRIBUTION_mode;
  parmesh->info.ifc_layers = PMMG_MVIFCS_NLAYERS;
  parmesh->info.grps_ratio = PMMG_GRPS_RATIO;
  parmesh->info.loadbalancing_mode = PMMG_LOADBALANCING_metis;
  parmesh->info.contiguous_mode = PMMG_CONTIG_DEF;
  parmesh->info.target_mesh_size =  PMMG_REMESHER_TARGET_MESH_SIZE;
  parmesh->info.metis_ratio =  PMMG_RATIO_MMG_METIS;
  parmesh->info.API_mode    =  PMMG_APIDISTRIB_faces;
  parmesh->info.PROctree_mode = PMMG_USE_PROCTREE;

  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
 #warning Option -nosurf imposed by default.
    mesh->info.nosurf = 1;
  }

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

  mem = (parmesh->memGloMax-parmesh->memMax)/(MMG5_MILLION*parmesh->ngrp) - 1;

  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_mem,(int)mem);
  }

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
  sol  = &parmesh->listgrp[0].sol;

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
  int         k,npmax,xpmax,nemax,xtmax;

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

    mem = (parmesh->memGloMax-parmesh->memMax)/(MMG5_MILLION*parmesh->ngrp) - 1;

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
  case PMMG_IPARAM_metisRatio :
    parmesh->info.metis_ratio = val;
    break;
  case PMMG_IPARAM_ifcLayers :
    parmesh->info.ifc_layers = val;
    break;
  case PMMG_IPARAM_APImode :
    parmesh->info.API_mode = val;
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
//      if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_nosurf,val) ) return 0;
      fprintf(stderr,"  ## Warning: Surfacic adaptation not implemented!\n");
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
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_dparameter(mesh,NULL,MMG3D_DPARAM_hmin,val) ) {
        return 0;
      }
    }
    break;
  case PMMG_DPARAM_hmax :
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
  return(MMG3D_Set_ithSols_inSolsAtVertices(parmesh->listgrp[0].sol,i,s));
}

int PMMG_Set_ithSol_inSolsAtVertices(PMMG_pParMesh parmesh,int i,double *s,int pos){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Set_ithSol_inSolsAtVertices(parmesh->listgrp[0].sol,i, s, pos));
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

  if ( parmesh->listgrp && parmesh->listgrp[0].mesh && parmesh->listgrp[0].sol ) {
    mesh = parmesh->listgrp[0].mesh;
    sol  = &parmesh->listgrp[0].sol;

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
  return(MMG3D_Get_ithSol_inSolsAtVertices(parmesh->listgrp[0].sol,i,s,pos));
}

int PMMG_Get_ithSols_inSolsAtVertices(PMMG_pParMesh parmesh,int i,double* s){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_ithSols_inSolsAtVertices(parmesh->listgrp[0].sol,i,s));

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
 
  /* Get the node communicator */
  pext_node_comm = &parmesh->ext_node_comm[ext_comm_index];
  nitem = pext_node_comm->nitem_to_share;

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
 
  /* Get the face communicator */
  pext_face_comm = &parmesh->ext_face_comm[ext_comm_index];
  nitem = pext_face_comm->nitem;

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
  size_t         memAv,oldMemMax;

 
  /* Meshes are merged in grp 0 */
  int_node_comm = parmesh->int_node_comm;
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;


  /** 1) Store node index in intvalues */
  PMMG_TRANSFER_AVMEM_TO_PARMESH(parmesh,memAv,oldMemMax);
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
  size_t         memAv,oldMemMax;

  /* Meshes are merged in grp 0 */
  int_face_comm = parmesh->int_face_comm;
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;

  PMMG_TRANSFER_AVMEM_TO_PARMESH(parmesh,memAv,oldMemMax);
  PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh,memAv,oldMemMax);


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
  PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,mesh,memAv,oldMemMax);
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
  size_t         memAv,oldMemMax;

  /* Meshes are merged in grp 0 */
  int_node_comm = parmesh->int_node_comm;
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;

  PMMG_TRANSFER_AVMEM_TO_PARMESH(parmesh,memAv,oldMemMax);
  PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh,memAv,oldMemMax);

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

  PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,mesh,memAv,oldMemMax); 
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
  size_t         memAv,oldMemMax;

  /* Meshes are merged in grp 0 */
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;

  PMMG_TRANSFER_AVMEM_TO_PARMESH(parmesh,memAv,oldMemMax);
  PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh,memAv,oldMemMax);

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

  PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,mesh,memAv,oldMemMax);
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
  size_t         memAv,oldMemMax;

  /* Meshes are merged in grp 0 */
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;

  PMMG_TRANSFER_AVMEM_TO_PARMESH(parmesh,memAv,oldMemMax);
  PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh,memAv,oldMemMax);

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

  PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,mesh,memAv,oldMemMax); 
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
  size_t         memAv,oldMemMax;

  /* Meshes are merged in grp 0 */
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;

  PMMG_TRANSFER_AVMEM_TO_PARMESH(parmesh,memAv,oldMemMax);
  PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh,memAv,oldMemMax);

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

  PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,mesh,memAv,oldMemMax);
  MMG5_DEL_MEM(mesh,hash.item);
  MMG5_DEL_MEM(mesh,hashPair.item);
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

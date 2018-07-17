/**
 * \file API_functions_pmmg.c
 * \brief C API functions definitions for PARMMG library.
 * \copyright GNU Lesser General Public License.
 *
 * C API for PARMMG library.
 *
 */
#include "parmmg.h"

int PMMG_Init_parMesh( PMMG_pParMesh *parmesh ,MPI_Comm comm)
{
  PMMG_pGrp grp = NULL;

  /* ParMesh allocation */
  assert ( (*parmesh == NULL) && "trying to initialize non empty parmesh" );
   *parmesh = calloc( 1, sizeof(PMMG_ParMesh) );
  if ( *parmesh == NULL )
    goto fail_pmesh;

  /* Assign some values to memory related fields to begin working with */
  (*parmesh)->memGloMax = 4 * 1024L * 1024L;
  (*parmesh)->memMax = 4 * 1024L * 1024L;
  (*parmesh)->memCur = sizeof(PMMG_ParMesh);

  /** Init Group */
  (*parmesh)->ngrp = 1;
  PMMG_CALLOC(*parmesh,(*parmesh)->listgrp,1,PMMG_Grp,
              "allocating groups container", goto fail_grplst );
  grp = &(*parmesh)->listgrp[0];
  grp->mesh = NULL;
  grp->met  = NULL;
  grp->disp = NULL;
  if ( 1 != MMG3D_Init_mesh( MMG5_ARG_start,
                             MMG5_ARG_ppMesh, &grp->mesh,
                             MMG5_ARG_ppMet, &grp->met,
                             MMG5_ARG_end ) )
    goto fail_mesh;


  PMMG_Init_parameters(*parmesh,comm);

  return 1;

fail_mesh:
    PMMG_DEL_MEM(*parmesh,(*parmesh)->listgrp,1,PMMG_Grp,
                 "deallocating groups container");

fail_grplst:
  (*parmesh)->ngrp = 0;
  (*parmesh)->memMax = 0;
  (*parmesh)->memCur = 0;
  free( *parmesh );
  *parmesh = NULL;

fail_pmesh:
  return 0;
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

int PMMG_Set_inputSolName(PMMG_pParMesh parmesh, const char* solin) {
  MMG5_pMesh mesh;
  MMG5_pSol  sol;
  int        k,ier;

  ier = 1;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    sol  = parmesh->listgrp[k].sol;
    ier  = MG_MIN( ier, MMG3D_Set_inputSolName(mesh,sol,solin) );
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

int PMMG_Set_outputSolName(PMMG_pParMesh parmesh, const char* solout) {
  MMG5_pMesh mesh;
  MMG5_pSol  sol;
  int        k,ier;

  ier = 1;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    sol  = parmesh->listgrp[k].sol;
    ier  = MG_MIN( ier, MMG3D_Set_outputSolName(mesh,sol,solout) );
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
  int        k,flag;

  memset(&parmesh->info,0, sizeof(PMMG_Info));

  parmesh->info.mem    = PMMG_UNSET; /* [n/-1]   ,Set memory size to n Mbytes/keep the default value */
  parmesh->ddebug      = 0;
  parmesh->niter       = 1;

  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    /* Set Mmg verbosity to 0 */
    mesh->info.imprim = 0;
  }

  /* Init MPI data */
  parmesh->comm   = comm;

  MPI_Initialized(&flag);
  if ( flag ) {
    MPI_Comm_size( parmesh->comm, &parmesh->nprocs );
    MPI_Comm_rank( parmesh->comm, &parmesh->myrank );
  }
  else {
    parmesh->nprocs = 1;
    parmesh->myrank = 0;
  }

  /* ParMmg verbosity */
  if ( !parmesh->myrank ) {
    parmesh->info.imprim = 1;
  }
  else {
    parmesh->info.imprim = 0;
  }
  parmesh->info.imprim0  = 1;

  /* Default memory */
  PMMG_parmesh_SetMemGloMax( parmesh, 0 );
}

int PMMG_Set_iparameter(PMMG_pParMesh parmesh, int iparam,int val){
  MMG5_pMesh  mesh;
  MMG5_pSol   met;
  int         k,mem;

  switch ( iparam ) {
  case PMMG_IPARAM_verbose :
    if ( !parmesh->myrank ) {
      parmesh->info.imprim = val;
    }
    parmesh->info.imprim0 = val;

    break;
  case PMMG_IPARAM_mem :
    if ( val <= 0 ) {
      fprintf( stdout,
        "  ## Warning: maximal memory authorized must be strictly positive.\n");
      fprintf(stdout,"  Reset to default value.\n");
    } else
      parmesh->memMax = val;

    mem = (int)(val/parmesh->ngrp);

    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_mem,mem) ) return 0;
    }
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
      if ( !MMG3D_Set_iparameter(mesh,NULL,MMG3D_IPARAM_nosurf,val) ) return 0;
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
  MMG5_pSol   met;
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
      return(0);
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
  default :
    fprintf(stderr,"  ## Error: unknown type of parameter\n");
    return 0;
  }

  return 1;
}

int PMMG_Free_all( PMMG_pParMesh *parmesh )
{
  PMMG_parmesh_Free_Comm( *parmesh );

  PMMG_parmesh_Free_Listgrp( *parmesh );

  _MMG5_SAFE_FREE(*parmesh);

  return 1;
}

int PMMG_Set_grpSize(PMMG_pGrp grp, int np, int ne, int nprism, int nt,
                     int nquad, int na, int typEntity, int typSol, int typMet,
                     int *loc2glob){
  // Set mesh size
  if(!MMG3D_Set_meshSize(grp->mesh, np, ne, nprism, nt, nquad, na)){
    fprintf(stderr,"  ## Error in setting group mesh size.\n");
    return(0);
  }
  // Set physical solution size
  if(!MMG3D_Set_solSize(grp->mesh, grp->sol, typEntity, np, typSol)){
    fprintf(stderr,"  ## Error in setting group solution size.\n");
    return(0);
  }
  // Set metrics size
  if(!MMG3D_Set_solSize(grp->mesh, grp->met, typEntity, np, typMet)){
    fprintf(stderr,"  ## Error in setting group metrics size.\n");
    return(0);
  }
  return(1);
}

int PMMG_Set_vertex(PMMG_pGrp grp, double c0, double c1, double c2,
                    int ref, int pos){
  return(MMG3D_Set_vertex(grp->mesh, c0, c1, c2, ref, pos));
}

int PMMG_Set_vertices(PMMG_pGrp grp, double *vertices,int *refs){
  return(MMG3D_Set_vertices(grp->mesh, vertices, refs));
}

int PMMG_Set_tetrahedron(PMMG_pGrp grp, int v0, int v1, int v2, int v3,
                         int ref, int pos){
  return(MMG3D_Set_tetrahedron(grp->mesh, v0, v1, v2, v3, ref, pos));
}

int PMMG_Set_tetrahedra(PMMG_pGrp grp, int *tetra, int *refs){
  return(MMG3D_Set_tetrahedra(grp->mesh, tetra, refs));
}

int PMMG_Set_prism(PMMG_pGrp grp, int v0, int v1, int v2,
                   int v3, int v4, int v5, int ref, int pos){
  return(MMG3D_Set_prism(grp->mesh, v0, v1, v2, v3, v4, v5, ref, pos));
}

int PMMG_Set_prisms(PMMG_pGrp grp, int *prisms, int *refs){
  return(MMG3D_Set_prisms(grp->mesh, prisms, refs));
}

int PMMG_Set_triangle(PMMG_pGrp grp, int v0, int v1, int v2,
                      int ref,int pos){
  return(MMG3D_Set_triangle(grp->mesh, v0, v1, v2, ref, pos));
}

int PMMG_Set_triangles(PMMG_pGrp grp, int *tria, int *refs){
  return(MMG3D_Set_triangles(grp->mesh, tria, refs));
}

int PMMG_Set_quadrilateral(PMMG_pGrp grp, int v0, int v1,
                           int v2, int v3, int ref,int pos){
  return(MMG3D_Set_quadrilateral(grp->mesh, v0, v1, v2, v3, ref, pos));
}

int PMMG_Set_quadrilaterals(PMMG_pGrp grp, int *quads, int *refs){
  return(MMG3D_Set_quadrilaterals(grp->mesh, quads, refs));
}

int PMMG_Set_edge(PMMG_pGrp grp, int v0, int v1, int ref, int pos){
  return(MMG3D_Set_edge(grp->mesh, v0, v1, ref, pos));
}

int PMMG_Set_corner(PMMG_pGrp grp, int k){
  return(MMG3D_Set_corner(grp->mesh, k));
}

int PMMG_Set_requiredVertex(PMMG_pGrp grp, int k){
  return(MMG3D_Set_requiredVertex(grp->mesh, k));
}

int PMMG_Set_requiredTetrahedron(PMMG_pGrp grp, int k){
  return(MMG3D_Set_requiredTetrahedron(grp->mesh, k));
}

int PMMG_Set_requiredTetrahedra(PMMG_pGrp grp, int *reqIdx, int nreq){
  return(MMG3D_Set_requiredTetrahedra(grp->mesh, reqIdx, nreq));
}

int PMMG_Set_requiredTriangle(PMMG_pGrp grp, int k){
  return(MMG3D_Set_requiredTriangle(grp->mesh, k));
}

int PMMG_Set_requiredTriangles(PMMG_pGrp grp, int *reqIdx, int nreq){
  return(MMG3D_Set_requiredTriangles(grp->mesh, reqIdx, nreq));
}

int PMMG_Set_ridge(PMMG_pGrp grp, int k){
  return(MMG3D_Set_ridge(grp->mesh, k));
}

int PMMG_Set_requiredEdge(PMMG_pGrp grp, int k){
  return(MMG3D_Set_requiredEdge(grp->mesh, k));
}

int PMMG_Set_normalAtVertex(PMMG_pGrp grp, int k, double n0, double n1,
                            double n2){
  return(MMG3D_Set_normalAtVertex(grp->mesh, k, n0, n1, n2));
}

int PMMG_Set_scalarSol(PMMG_pGrp grp, double s,int pos){
  return(MMG3D_Set_scalarSol(grp->sol, s, pos));
}

int PMMG_Set_scalarSols(PMMG_pGrp grp, double *sol){
  return(MMG3D_Set_scalarSols(grp->sol, sol));
}

int PMMG_Set_vectorSol(PMMG_pGrp grp, double vx,double vy, double vz,
                       int pos){
  return(MMG3D_Set_vectorSol(grp->sol, vx, vy, vz, pos));
}

int PMMG_Set_vectorSols(PMMG_pGrp grp, double *sols){
  return(MMG3D_Set_vectorSols(grp->sol, sols));
}

int PMMG_Set_tensorSol(PMMG_pGrp grp, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos){
  return(MMG3D_Set_tensorSol(grp->sol, m11, m12, m13, m22, m23, m33, pos));
}

int PMMG_Set_tensorSols(PMMG_pGrp grp, double *sols){
  return(MMG3D_Set_tensorSols(grp->sol, sols));
}

int PMMG_Set_scalarMet(PMMG_pGrp grp, double m,int pos){
  return(MMG3D_Set_scalarSol(grp->met, m, pos));
}

int PMMG_Set_scalarMets(PMMG_pGrp grp, double *met){
  return(MMG3D_Set_scalarSols(grp->met, met));
}

int PMMG_Set_vectorMet(PMMG_pGrp grp, double vx,double vy, double vz,
                       int pos){
  return(MMG3D_Set_vectorSol(grp->met, vx, vy, vz, pos));
}

int PMMG_Set_vectorMets(PMMG_pGrp grp, double *mets){
  return(MMG3D_Set_vectorSols(grp->met, mets));
}

int PMMG_Set_tensorMet(PMMG_pGrp grp, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos){
  return(MMG3D_Set_tensorSol(grp->met, m11, m12, m13, m22, m23, m33, pos));
}

int PMMG_Set_tensorMets(PMMG_pGrp grp, double *mets){
  return(MMG3D_Set_tensorSols(grp->met, mets));
}

int PMMG_Get_grpSize(PMMG_pGrp grp, int* np, int* ne,int *nprism, int* nt,
                     int* nquad, int* na, int* typEntity,
                     int* typSol, int* typMet){
  if(!MMG3D_Get_meshSize(grp->mesh, np, ne, nprism, nt, nquad, na)){
    fprintf(stderr,"  ## Error in getting group mesh size.\n");
    return(0);
  }
  if(!MMG3D_Get_solSize(grp->mesh, grp->sol, typEntity, np, typSol)){
    // FIXME np is read again.
    fprintf(stderr,"  ## Error in getting group solution size.\n");
    return(0);
  }
  if(!MMG3D_Get_solSize(grp->mesh, grp->met, typEntity, np, typMet)){
    // FIXME np is read again, as typEntity.
    fprintf(stderr,"  ## Error in getting group metrics size.\n");
    return(0);
  }
  return(1);
}

int PMMG_Get_vertex(PMMG_pGrp grp, double* c0, double* c1, double* c2,
                    int* ref,int* isCorner, int* isRequired){
  return(MMG3D_Get_vertex(grp->mesh, c0, c1, c2, ref, isCorner, isRequired));
}

int PMMG_Get_vertices(PMMG_pGrp grp, double* vertices, int* refs,
                      int* areCorners, int* areRequired){
  return(MMG3D_Get_vertices(grp->mesh, vertices, refs, areCorners, areRequired));
}

int PMMG_Get_tetrahedron(PMMG_pGrp grp, int* v0, int* v1, int* v2,
                         int* v3,int* ref, int* isRequired){
  return(MMG3D_Get_tetrahedron(grp->mesh, v0, v1, v2, v3, ref, isRequired));
}

int PMMG_Get_tetrahedra(PMMG_pGrp grp, int* tetra,int* refs,
                        int* areRequired){
  return(MMG3D_Get_tetrahedra(grp->mesh, tetra, refs, areRequired));
}

int PMMG_Get_prism(PMMG_pGrp grp, int* v0, int* v1, int* v2,
                   int* v3,int* v4,int* v5,int* ref, int* isRequired){
  return(MMG3D_Get_prism(grp->mesh, v0, v1, v2, v3, v4, v5, ref, isRequired));
}

int PMMG_Get_prisms(PMMG_pGrp grp, int* prisms,int* refs,
                    int* areRequired){
  return(MMG3D_Get_prisms(grp->mesh, prisms, refs, areRequired));
}

int PMMG_Get_triangle(PMMG_pGrp grp, int* v0, int* v1, int* v2, int* ref,
                      int* isRequired){
  return(MMG3D_Get_triangle(grp->mesh, v0, v1, v2, ref, isRequired));
}

int PMMG_Get_triangles(PMMG_pGrp grp, int* tria, int* refs,
                       int* areRequired){
  return(MMG3D_Get_triangles(grp->mesh, tria, refs, areRequired));
}

int PMMG_Get_quadrilateral(PMMG_pGrp grp, int* v0, int* v1, int* v2,int* v3,
                            int* ref, int* isRequired){
  return(MMG3D_Get_quadrilateral(grp->mesh, v0, v1, v2, v3, ref, isRequired));
}

int PMMG_Get_quadrilaterals(PMMG_pGrp grp, int* quads, int* refs,
                            int* areRequired){
  return(MMG3D_Get_quadrilaterals(grp->mesh, quads, refs, areRequired));
}

int PMMG_Get_edge(PMMG_pGrp grp, int* e0, int* e1, int* ref,
                   int* isRidge, int* isRequired){
  return(MMG3D_Get_edge(grp->mesh, e0, e1, ref, isRidge, isRequired));
}

int PMMG_Get_normalAtVertex(PMMG_pGrp grp, int k, double *n0, double *n1,
                              double *n2){
  return(MMG3D_Get_normalAtVertex(grp->mesh, k, n0, n1, n2));
}

int PMMG_Get_scalarSol(PMMG_pGrp grp, double* s){
  return(MMG3D_Get_scalarSol(grp->sol, s));
}

int PMMG_Get_scalarSols(PMMG_pGrp grp, double* s){
  return(MMG3D_Get_scalarSols(grp->sol, s));
}

int PMMG_Get_vectorSol(PMMG_pGrp grp, double* vx, double* vy, double* vz){
  return(MMG3D_Get_vectorSol(grp->sol, vx, vy, vz));
}

int PMMG_Get_vectorSols(PMMG_pGrp grp, double* sols){
  return(MMG3D_Get_vectorSols(grp->sol, sols));
}

int PMMG_Get_tensorSol(PMMG_pGrp grp, double *m11,double *m12, double *m13,
                       double *m22,double *m23, double *m33){
  return(MMG3D_Get_tensorSol(grp->sol, m11, m12, m13, m22, m23, m33));
}

int PMMG_Get_tensorSols(PMMG_pGrp grp, double *sols){
  return(MMG3D_Get_tensorSols(grp->sol, sols));
}

int PMMG_Get_scalarMet(PMMG_pGrp grp, double* m){
  return(MMG3D_Get_scalarSol(grp->met, m));
}

int PMMG_Get_scalarMets(PMMG_pGrp grp, double* m){
  return(MMG3D_Get_scalarSols(grp->met, m));
}

int PMMG_Get_vectorMet(PMMG_pGrp grp, double* vx, double* vy, double* vz){
  return(MMG3D_Get_vectorSol(grp->met, vx, vy, vz));
}

int PMMG_Get_vectorMets(PMMG_pGrp grp, double* mets){
  return(MMG3D_Get_vectorSols(grp->met, mets));
}

int PMMG_Get_tensorMet(PMMG_pGrp grp, double *m11,double *m12, double *m13,
                       double *m22,double *m23, double *m33){
  return(MMG3D_Get_tensorSol(grp->met, m11, m12, m13, m22, m23, m33));
}

int PMMG_Get_tensorMets(PMMG_pGrp grp, double *mets){
  return(MMG3D_Get_tensorSols(grp->met, mets));
}

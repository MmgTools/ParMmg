/**
 * \file API_functions_pmmg.c
 * \brief C API functions definitions for PARMMG library.
 * \copyright GNU Lesser General Public License.
 *
 * C API for PARMMG library.
 *
 */
#include "parmmg.h"


/**
 * \param parmesh pointer toward a parmesh structure
 * \param comm    external communicator to be freed
 *
 * deallocate all internal communicator's fields
 */
static void PMMG_parmesh_int_comm_free( PMMG_pParMesh parmesh,
                                        PMMG_pInt_comm comm )
{
  if ( comm == NULL )
    return;

  if ( NULL != comm->intvalues ) {
    assert ( comm->nitem != 0 && "incorrect parameters in internal communicator" );
    PMMG_DEL_MEM(parmesh,comm->intvalues,comm->nitem,int,"int comm int array");
  }
  if ( NULL != comm->doublevalues ) {
    assert ( comm->nitem != 0 && "incorrect parameters in internal communicator" );
    PMMG_DEL_MEM(parmesh,
                 comm->doublevalues,comm->nitem,double,"int comm double array");
  }
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param comm    external communicator to be freed
 * \param ncomm   parameter ncomm
 *
 * deallocate all external communicators's fields
 */
static void PMMG_parmesh_ext_comm_free( PMMG_pParMesh parmesh,
                                        PMMG_pExt_comm comm, int ncomm )
{
  int i = 0;

  if ( comm == NULL )
    return;

  for( i = 0; i < ncomm; ++i ) {
    if ( NULL != comm->int_comm_index ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(parmesh,comm->int_comm_index,comm->nitem,int,"ext comm int array");
    }
    if ( NULL != comm->itosend ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(parmesh,comm->itosend,comm->nitem,int,"ext comm itosend array");
    }
    if ( NULL != comm->itorecv ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(parmesh,comm->itorecv,comm->nitem,int,"ext comm itorecv array");
    }
    if ( NULL != comm->rtosend ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(parmesh,comm->rtosend,comm->nitem,int,"ext comm rtosend array");
    }
    if ( NULL != comm->rtorecv ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(parmesh,comm->rtorecv,comm->nitem,int,"ext comm rtorecv array");
    }
  }
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param idx1    node2int_node_comm_index1 to be freed
 * \param idx2    node2int_node_comm_index2 to be freed
 * \param n       pointer to node2int_node_comm_nitem size
 *
 * Deallocate all the MMG3D meshes and their communicators and zero the size
 */
static void PMMG_parmesh_grp_comm_free( PMMG_pParMesh parmesh,
                                        int **idx1, int **idx2, int *n )
{
  PMMG_DEL_MEM(parmesh,*idx1,*n,int,"group communicator");
  PMMG_DEL_MEM(parmesh,*idx2,*n,int,"group communicator");
  *n = 0;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param listgrp group of MMG3D meshes in parmesh
 * \param ngrp    number of mmg meshes in listgrp
 *
 * Deallocate all the MMG3D meshes and their communicators
 */
void PMMG_listgrp_free( PMMG_pParMesh parmesh, PMMG_pGrp *listgrp, int ngrp )
{
  int k;

  for ( k = 0; k < ngrp; ++k )
    PMMG_grp_free( parmesh, listgrp[0] + k );

  PMMG_DEL_MEM(parmesh,*listgrp,ngrp,PMMG_Grp,"Deallocating listgrp container");
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param grp     group to free
 *
 * Deallocate all the MMG3D meshes and their communicators
 */
void PMMG_grp_free( PMMG_pParMesh parmesh, PMMG_pGrp grp )
{
  PMMG_parmesh_grp_comm_free( parmesh,
                              &grp->node2int_node_comm_index1,
                              &grp->node2int_node_comm_index2,
                              &grp->nitem_int_node_comm);
  PMMG_parmesh_grp_comm_free( parmesh,
                              &grp->edge2int_edge_comm_index1,
                              &grp->edge2int_edge_comm_index2,
                              &grp->nitem_int_edge_comm);
  PMMG_parmesh_grp_comm_free( parmesh,
                              &grp->face2int_face_comm_index1,
                              &grp->face2int_face_comm_index2,
                              &grp->nitem_int_face_comm);
  MMG3D_Free_all( MMG5_ARG_start,
                  MMG5_ARG_ppMesh, &grp->mesh,
                  MMG5_ARG_ppMet, &grp->met,
                  MMG5_ARG_end );
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * Free parmesh communicators that are allocated
 */
void PMMG_parmesh_Free_Comm( PMMG_pParMesh parmesh )
{
  PMMG_parmesh_int_comm_free( parmesh, parmesh->int_node_comm );
  PMMG_parmesh_int_comm_free( parmesh, parmesh->int_edge_comm );
  PMMG_parmesh_int_comm_free( parmesh, parmesh->int_face_comm );

  PMMG_parmesh_ext_comm_free( parmesh, parmesh->ext_node_comm, parmesh->next_node_comm );
  PMMG_DEL_MEM(parmesh, parmesh->ext_node_comm, parmesh->next_node_comm,
            PMMG_Ext_comm, "ext node comm");
  PMMG_parmesh_ext_comm_free( parmesh, parmesh->ext_edge_comm, parmesh->next_edge_comm );
  PMMG_DEL_MEM(parmesh, parmesh->ext_edge_comm, parmesh->next_edge_comm,
            PMMG_Ext_comm, "ext edge comm");
  PMMG_parmesh_ext_comm_free( parmesh, parmesh->ext_face_comm, parmesh->next_face_comm );
  PMMG_DEL_MEM(parmesh, parmesh->ext_face_comm, parmesh->next_face_comm,
            PMMG_Ext_comm, "ext face comm");
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * Free parmesh listgrp that are allocated
 */
void PMMG_parmesh_Free_Listgrp( PMMG_pParMesh parmesh )
{
  PMMG_listgrp_free( parmesh, &parmesh->listgrp, parmesh->ngrp );

  PMMG_DEL_MEM(parmesh,parmesh->listgrp,1,PMMG_Grp,"deallocating groups container");
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * Free any parmesh members that are allocated
 */
void PMMG_parmesh_Free( PMMG_pParMesh parmesh )
{
  PMMG_parmesh_Free_Comm( parmesh );

  PMMG_parmesh_Free_Listgrp( parmesh );
}


/**
 * \param parmesh pointer toward a parmesh structure
 * \param val     exit value
 *
 * Controlled parmmg termination:
 *   Deallocate parmesh struct and its allocated members
 *   If this is an unsuccessful exit call abort to cancel any remaining processes
 *   Call MPI_Finalize / exit
 */
#warning NIKOS: MPI_Finalize might not be desirable here
void PMMG_exit_and_free( PMMG_pParMesh parmesh, const int val )
{
  PMMG_parmesh_Free( parmesh );
  if ( val != PMMG_SUCCESS )
    MPI_Abort( parmesh->comm, val );
  MPI_Finalize();
  exit( val );
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param comm MPI communicator for ParMmg
 *
 * Initialization of the input parameters.
 *
 */
void PMMG_Init_parameters(PMMG_pParMesh parmesh,MPI_Comm comm) {
  MMG5_pMesh mesh;
  int        k,flag;

  memset(&parmesh->info,0, sizeof(PMMG_Info));

  parmesh->info.mem    = PMMG_UNSET; /* [n/-1]   ,Set memory size to n Mbytes/keep the default value */
  parmesh->ddebug      = 0;
  parmesh->niter       = 1;

  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[0].mesh;
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


/**
 * \param parmesh pointer toward a parmesh structure
 * \param comm MPI communicator for ParMmg
 *
 * \return 0 on error
 *         1 on success
 *
 * allocate a parmesh struct with a single mesh struct and initialize
 * some of the struct fields
 *
 */
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

/**
 * \param parmesh pointer toward a parmesh structure
 * \param iparam  parameter enumeration option
 * \param val     parameter value
 *
 * \return 0 on error
 *         1 on success
 *
 * Set integer parameters.
 *
 */
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
    return(0);
  }

  return(1);
}

/**
 * \file API_functions.c
 * \brief C API functions definitions for PARMMG library.
 * \copyright GNU Lesser General Public License.
 *
 * C API for PARMMG library.
 *
 */

#include "libparmmg.h"

/**
 * \param .
 *
 * Allocation and initialization of parmmg structures.
 *
 */
int PMMG_InitParMesh() {
  return(1);
}

/**
 * Set group size.
 */
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

/**
 * Set group mesh vertex.
 */
int PMMG_Set_vertex(PMMG_pGrp grp, double c0, double c1, double c2,
                    int ref, int pos){
  return(MMG3D_Set_vertex(grp->mesh, c0, c1, c2, ref, pos));
}

/**
 * Set group mesh vertices.
 */
int PMMG_Set_vertices(PMMG_pGrp grp, double *vertices,int *refs){
  return(MMG3D_Set_vertices(grp->mesh, vertices, refs));
}

/**
 * Set group mesh tetrahedron.
 */
int PMMG_Set_tetrahedron(PMMG_pGrp grp, int v0, int v1, int v2, int v3,
                         int ref, int pos){
  return(MMG3D_Set_tetrahedron(grp->mesh, v0, v1, v2, v3, ref, pos));
}

/**
 * Set group mesh tetrahedra.
 */
int PMMG_Set_tetrahedra(PMMG_pGrp grp, int *tetra, int *refs){
  return(MMG3D_Set_tetrahedra(grp->mesh, tetra, refs));
}

/**
 * Set group mesh prism.
 */
int PMMG_Set_prism(PMMG_pGrp grp, int v0, int v1, int v2,
                   int v3, int v4, int v5, int ref, int pos){
  return(MMG3D_Set_prism(grp->mesh, v0, v1, v2, v3, v4, v5, ref, pos));
}

/**
 * Set group mesh prisms.
 */
int PMMG_Set_prisms(PMMG_pGrp grp, int *prisms, int *refs){
  return(MMG3D_Set_prisms(grp->mesh, prisms, refs));
}

/**
 * Set group mesh triangle.
 */
int PMMG_Set_triangle(PMMG_pGrp grp, int v0, int v1, int v2,
                      int ref,int pos){
  return(MMG3D_Set_triangle(grp->mesh, v0, v1, v2, ref, pos));
}

/**
 * Set group mesh triangles.
 */
int PMMG_Set_triangles(PMMG_pGrp grp, int *tria, int *refs){
  return(MMG3D_Set_triangles(grp->mesh, tria, refs));
}

/**
 * Set group mesh quadrilateral.
 */
int PMMG_Set_quadrilateral(PMMG_pGrp grp, int v0, int v1,
                           int v2, int v3, int ref,int pos){
  return(MMG3D_Set_quadrilateral(grp->mesh, v0, v1, v2, v3, ref, pos));
}

/**
 * Set group mesh quadrilaterals.
 */
int PMMG_Set_quadrilaterals(PMMG_pGrp grp, int *quads, int *refs){
  return(MMG3D_Set_quadrilaterals(grp->mesh, quads, refs));
}

/**
 * Set group mesh edge.
 */
int PMMG_Set_edge(PMMG_pGrp grp, int v0, int v1, int ref, int pos){
  return(MMG3D_Set_edge(grp->mesh, v0, v1, ref, pos));
}

/**
 * Set group mesh corner.
 */
int PMMG_Set_corner(PMMG_pGrp grp, int k){
  return(MMG3D_Set_corner(grp->mesh, k));
}

/**
 * Set group mesh required vertex.
 */
int PMMG_Set_requiredVertex(PMMG_pGrp grp, int k){
  return(MMG3D_Set_requiredVertex(grp->mesh, k));
}

/**
 * Set group mesh required tetrahedron.
 */
int PMMG_Set_requiredTetrahedron(PMMG_pGrp grp, int k){
  return(MMG3D_Set_requiredTetrahedron(grp->mesh, k));
}

/**
 * Set group mesh required tetrahedra.
 */
int PMMG_Set_requiredTetrahedra(PMMG_pGrp grp, int *reqIdx, int nreq){
  return(MMG3D_Set_requiredTetrahedra(grp->mesh, reqIdx, nreq));
}

/**
 * Set group mesh required triangle.
 */
int PMMG_Set_requiredTriangle(PMMG_pGrp grp, int k){
  return(MMG3D_Set_requiredTriangle(grp->mesh, k));
}

/**
 * Set group mesh required triangles.
 */
int PMMG_Set_requiredTriangles(PMMG_pGrp grp, int *reqIdx, int nreq){
  return(MMG3D_Set_requiredTriangles(grp->mesh, reqIdx, nreq));
}


/**
 * Set group mesh ridge.
 */
int PMMG_Set_ridge(PMMG_pGrp grp, int k){
  return(MMG3D_Set_ridge(grp->mesh, k));
}

/**
 * Set group mesh required edge.
 */
int PMMG_Set_requiredEdge(PMMG_pGrp grp, int k){
  return(MMG3D_Set_requiredEdge(grp->mesh, k));
}

/**
 * Set group mesh normal at vertex.
 */
int PMMG_Set_normalAtVertex(PMMG_pGrp grp, int k, double n0, double n1,
                            double n2){
  return(MMG3D_Set_normalAtVertex(grp->mesh, k, n0, n1, n2));
}

/**
 * Set group scalar solution at vertex.
 */
int PMMG_Set_scalarSol(PMMG_pGrp grp, double s,int pos){
  return(MMG3D_Set_scalarSol(grp->sol, s, pos));
}

/**
 * Set group scalar solution at vertices.
 */
int PMMG_Set_scalarSols(PMMG_pGrp grp, double *sol){
  return(MMG3D_Set_scalarSols(grp->sol, sol));
}

/**
 * Set group vectorial solution at vertex.
 */
int PMMG_Set_vectorSol(PMMG_pGrp grp, double vx,double vy, double vz,
                       int pos){
  return(MMG3D_Set_vectorSol(grp->sol, vx, vy, vz, pos));
}

/**
 * Set group vectorial solution at vertices.
 */
int PMMG_Set_vectorSols(PMMG_pGrp grp, double *sols){
  return(MMG3D_Set_vectorSols(grp->sol, sols));
}

/**
 * Set group tensorial solution at vertex.
 */
int PMMG_Set_tensorSol(PMMG_pGrp grp, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos){
  return(MMG3D_Set_tensorSol(grp->sol, m11, m12, m13, m22, m23, m33, pos));
}

/**
 * Set group tensorial solution at vertices.
 */
int PMMG_Set_tensorSols(PMMG_pGrp grp, double *sols){
  return(MMG3D_Set_tensorSols(grp->sol, sols));
}

/**
 * Set group scalar metrics at vertex.
 */
int PMMG_Set_scalarMet(PMMG_pGrp grp, double m,int pos){
  return(MMG3D_Set_scalarSol(grp->met, m, pos));
}

/**
 * Set group scalar metrics at vertices.
 */
int PMMG_Set_scalarMets(PMMG_pGrp grp, double *met){
  return(MMG3D_Set_scalarSols(grp->met, met));
}

/**
 * Set group vectorial metrics at vertex.
 */
int PMMG_Set_vectorMet(PMMG_pGrp grp, double vx,double vy, double vz,
                       int pos){
  return(MMG3D_Set_vectorSol(grp->met, vx, vy, vz, pos));
}

/**
 * Set group vectorial metrics at vertices.
 */
int PMMG_Set_vectorMets(PMMG_pGrp grp, double *mets){
  return(MMG3D_Set_vectorSols(grp->met, mets));
}

/**
 * Set group tensorial metrics at vertex.
 */
int PMMG_Set_tensorMet(PMMG_pGrp grp, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos){
  return(MMG3D_Set_tensorSol(grp->met, m11, m12, m13, m22, m23, m33, pos));
}

/**
 * Set group tensorial metrics at vertices.
 */
int PMMG_Set_tensorMets(PMMG_pGrp grp, double *mets){
  return(MMG3D_Set_tensorSols(grp->met, mets));
}

/**
 * Get group sizes.
 */
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

/**
 * Get group mesh vertex.
 */
int PMMG_Get_vertex(PMMG_pGrp grp, double* c0, double* c1, double* c2,
                    int* ref,int* isCorner, int* isRequired){
  return(MMG3D_Get_vertex(grp->mesh, c0, c1, c2, ref, isCorner, isRequired));
}

/**
 * Get group mesh vertices.
 */
int PMMG_Get_vertices(PMMG_pGrp grp, double* vertices, int* refs,
                      int* areCorners, int* areRequired){
  return(MMG3D_Get_vertices(grp->mesh, vertices, refs, areCorners, areRequired));
}

/**
 * Get group mesh tetrahedron.
 */
int PMMG_Get_tetrahedron(PMMG_pGrp grp, int* v0, int* v1, int* v2,
                         int* v3,int* ref, int* isRequired){
  return(MMG3D_Get_tetrahedron(grp->mesh, v0, v1, v2, v3, ref, isRequired));
}

/**
 * Get group mesh tetrahedra.
 */
int PMMG_Get_tetrahedra(PMMG_pGrp grp, int* tetra,int* refs,
                        int* areRequired){
  return(MMG3D_Get_tetrahedra(grp->mesh, tetra, refs, areRequired));
}

/**
 * Get group mesh prism.
 */
int PMMG_Get_prism(PMMG_pGrp grp, int* v0, int* v1, int* v2,
                   int* v3,int* v4,int* v5,int* ref, int* isRequired){
  return(MMG3D_Get_prism(grp->mesh, v0, v1, v2, v3, v4, v5, ref, isRequired));
}

/**
 * Get group mesh prisms.
 */
int PMMG_Get_prisms(PMMG_pGrp grp, int* prisms,int* refs,
                    int* areRequired){
  return(MMG3D_Get_prisms(grp->mesh, prisms, refs, areRequired));
}

/**
 * Get group mesh triangle.
 */
int PMMG_Get_triangle(PMMG_pGrp grp, int* v0, int* v1, int* v2, int* ref,
                      int* isRequired){
  return(MMG3D_Get_triangle(grp->mesh, v0, v1, v2, ref, isRequired));
}

/**
 * Get group mesh triangles.
 */
int PMMG_Get_triangles(PMMG_pGrp grp, int* tria, int* refs,
                       int* areRequired){
  return(MMG3D_Get_triangles(grp->mesh, tria, refs, areRequired));
}

/**
 * Get group mesh quadrilateral.
 */
int PMMG_Get_quadrilateral(PMMG_pGrp grp, int* v0, int* v1, int* v2,int* v3,
                            int* ref, int* isRequired){
  return(MMG3D_Get_quadrilateral(grp->mesh, v0, v1, v2, v3, ref, isRequired));
}

/**
 * Get group mesh quadrilaterals.
 */
int PMMG_Get_quadrilaterals(PMMG_pGrp grp, int* quads, int* refs,
                            int* areRequired){
  return(MMG3D_Get_quadrilaterals(grp->mesh, quads, refs, areRequired));
}

/**
 * Get group mesh edge.
 */
int PMMG_Get_edge(PMMG_pGrp grp, int* e0, int* e1, int* ref,
                   int* isRidge, int* isRequired){
  return(MMG3D_Get_edge(grp->mesh, e0, e1, ref, isRidge, isRequired));
}

/**
 * Get group mesh normal at vertex.
 */
int PMMG_Get_normalAtVertex(PMMG_pGrp grp, int k, double *n0, double *n1,
                              double *n2){
  return(MMG3D_Get_normalAtVertex(grp->mesh, k, n0, n1, n2));
}


/**
 * Get group scalar solution at vertex.
 */
int PMMG_Get_scalarSol(PMMG_pGrp grp, double* s){
  return(MMG3D_Get_scalarSol(grp->sol, s));
}

/**
 * Get group scalar solution at vertices.
 */
int PMMG_Get_scalarSols(PMMG_pGrp grp, double* s){
  return(MMG3D_Get_scalarSols(grp->sol, s));
}

/**
 * Get group vectorial solution at vertex.
 */
int PMMG_Get_vectorSol(PMMG_pGrp grp, double* vx, double* vy, double* vz){
  return(MMG3D_Get_vectorSol(grp->sol, vx, vy, vz));
}

/**
 * Get group vectorial solution at vertices.
 */
int PMMG_Get_vectorSols(PMMG_pGrp grp, double* sols){
  return(MMG3D_Get_vectorSols(grp->sol, sols));
}

/**
 * Get group tensorial solution at vertex.
 */
int PMMG_Get_tensorSol(PMMG_pGrp grp, double *m11,double *m12, double *m13,
                       double *m22,double *m23, double *m33){
  return(MMG3D_Get_tensorSol(grp->sol, m11, m12, m13, m22, m23, m33));
}

/**
 * Get group tensorial solution at vertices.
 */
int PMMG_Get_tensorSols(PMMG_pGrp grp, double *sols){
  return(MMG3D_Get_tensorSols(grp->sol, sols));
}


/**
 * Get group scalar metrics at vertex.
 */
int PMMG_Get_scalarMet(PMMG_pGrp grp, double* m){
  return(MMG3D_Get_scalarSol(grp->met, m));
}

/**
 * Get group scalar metrics at vertices.
 */
int PMMG_Get_scalarMets(PMMG_pGrp grp, double* m){
  return(MMG3D_Get_scalarSols(grp->met, m));
}

/**
 * Get group vectorial metrics at vertex.
 */
int PMMG_Get_vectorMet(PMMG_pGrp grp, double* vx, double* vy, double* vz){
  return(MMG3D_Get_vectorSol(grp->met, vx, vy, vz));
}

/**
 * Get group vectorial metrics at vertices.
 */
int PMMG_Get_vectorMets(PMMG_pGrp grp, double* mets){
  return(MMG3D_Get_vectorSols(grp->met, mets));
}

/**
 * Get group tensorial metrics at vertex.
 */
int PMMG_Get_tensorMet(PMMG_pGrp grp, double *m11,double *m12, double *m13,
                       double *m22,double *m23, double *m33){
  return(MMG3D_Get_tensorSol(grp->met, m11, m12, m13, m22, m23, m33));
}

/**
 * Get group tensorial metrics at vertices.
 */
int PMMG_Get_tensorMets(PMMG_pGrp grp, double *mets){
  return(MMG3D_Get_tensorSols(grp->met, mets));
}

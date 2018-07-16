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
  PMMG_grp_comm_free( parmesh,&grp->node2int_node_comm_index1,
                              &grp->node2int_node_comm_index2,
                              &grp->nitem_int_node_comm);
  PMMG_grp_comm_free( parmesh,&grp->edge2int_edge_comm_index1,
                              &grp->edge2int_edge_comm_index2,
                              &grp->nitem_int_edge_comm);
  PMMG_grp_comm_free( parmesh,&grp->face2int_face_comm_index1,
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
  PMMG_int_comm_free( parmesh, parmesh->int_node_comm );
  PMMG_int_comm_free( parmesh, parmesh->int_edge_comm );
  PMMG_int_comm_free( parmesh, parmesh->int_face_comm );

  PMMG_ext_comm_free( parmesh, parmesh->ext_node_comm, parmesh->next_node_comm );
  PMMG_DEL_MEM(parmesh, parmesh->ext_node_comm, parmesh->next_node_comm,
            PMMG_Ext_comm, "ext node comm");
  PMMG_ext_comm_free( parmesh, parmesh->ext_edge_comm, parmesh->next_edge_comm );
  PMMG_DEL_MEM(parmesh, parmesh->ext_edge_comm, parmesh->next_edge_comm,
            PMMG_Ext_comm, "ext edge comm");
  PMMG_ext_comm_free( parmesh, parmesh->ext_face_comm, parmesh->next_face_comm );
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

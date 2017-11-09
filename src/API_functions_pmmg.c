/**
 * \file API_functions_pmmg.c
 * \brief C API functions definitions for PARMMG library.
 * \copyright GNU Lesser General Public License.
 *
 * C API for PARMMG library.
 *
 */
#include "parmmg.h"


static void pmesh_int_comm_free( PMMG_pParMesh pmesh, PMMG_pint_comm comm )
{
  if ( comm == NULL )
    return;

  if ( NULL != comm->intvalues ) {
    assert ( comm->nitem != 0 && "incorrect parameters in internal communicator" );
    PMMG_DEL_MEM(pmesh,comm->intvalues,comm->nitem,int,"int comm int array");
  }
  if ( NULL != comm->doublevalues ) {
    assert ( comm->nitem != 0 && "incorrect parameters in internal communicator" );
    PMMG_DEL_MEM(pmesh,comm->doublevalues,comm->nitem,double,"int comm double array");
  }
}

static void pmesh_ext_comm_free( PMMG_pParMesh pmesh, PMMG_pext_comm comm, int ncomm )
{
  int i = 0;

  if ( comm == NULL )
    return;

  for( i = 0; i < ncomm; ++i ) {
    if ( NULL != comm->int_comm_index ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(pmesh,comm->int_comm_index,comm->nitem,int,"ext comm int array");
    }
    if ( NULL != comm->itosend ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(pmesh,comm->itosend,comm->nitem,int,"ext comm itosend array");
    }
    if ( NULL != comm->itorecv ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(pmesh,comm->itorecv,comm->nitem,int,"ext comm itorecv array");
    }
    if ( NULL != comm->rtosend ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(pmesh,comm->rtosend,comm->nitem,int,"ext comm rtosend array");
    }
    if ( NULL != comm->rtorecv ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(pmesh,comm->rtorecv,comm->nitem,int,"ext comm rtorecv array");
    }
  }
}

static void pmesh_grp_comm_free( PMMG_pParMesh parmesh, int *idx1, int *idx2, int *n )
{
  PMMG_DEL_MEM(parmesh,idx1,*n,int,"group communicator");
  PMMG_DEL_MEM(parmesh,idx2,*n,int,"group communicator");
  *n = 0;
}

void PMMG_grp_free( PMMG_pParMesh parmesh, PMMG_pGrp *listgrp, int ngrp )
{
  int k = 0;
  for ( k = 0; k < ngrp; ++k ) {
    pmesh_grp_comm_free( parmesh,
        (*listgrp)[k].node2int_node_comm_index1,
        (*listgrp)[k].node2int_node_comm_index2,
        &(*listgrp)[k].nitem_int_node_comm);
    pmesh_grp_comm_free( parmesh,
        (*listgrp)[k].edge2int_edge_comm_index1,
        (*listgrp)[k].edge2int_edge_comm_index2,
        &(*listgrp)[k].nitem_int_edge_comm);
    pmesh_grp_comm_free( parmesh,
        (*listgrp)[k].face2int_face_comm_index1,
        (*listgrp)[k].face2int_face_comm_index2,
        &(*listgrp)[k].nitem_int_face_comm);
    MMG3D_Free_all( MMG5_ARG_start,
                    MMG5_ARG_ppMesh, &(*listgrp)[k].mesh,
                    MMG5_ARG_ppMet, &(*listgrp)[k].met,
                    MMG5_ARG_end );
  }
  PMMG_DEL_MEM(parmesh,*listgrp,ngrp,PMMG_Grp,"Deallocating listgrp container");
}

void PMMG_PMesh_Free( PMMG_pParMesh parmesh )
{
  PMMG_grp_free( parmesh, &parmesh->listgrp, parmesh->ngrp );

  pmesh_int_comm_free( parmesh, parmesh->int_node_comm );
  pmesh_int_comm_free( parmesh, parmesh->int_edge_comm );
  pmesh_int_comm_free( parmesh, parmesh->int_face_comm );

  pmesh_ext_comm_free( parmesh, parmesh->ext_node_comm, parmesh->next_node_comm );
  PMMG_DEL_MEM(parmesh, parmesh->ext_node_comm, parmesh->next_node_comm,
            PMMG_ext_comm, "ext node comm");
  pmesh_ext_comm_free( parmesh, parmesh->ext_edge_comm, parmesh->next_edge_comm );
  PMMG_DEL_MEM(parmesh, parmesh->ext_edge_comm, parmesh->next_edge_comm,
            PMMG_ext_comm, "ext edge comm");
  pmesh_ext_comm_free( parmesh, parmesh->ext_face_comm, parmesh->next_face_comm );
  PMMG_DEL_MEM(parmesh, parmesh->ext_face_comm, parmesh->next_face_comm,
            PMMG_ext_comm, "ext face comm");
  PMMG_DEL_MEM(parmesh,parmesh->listgrp,1,PMMG_Grp,"deallocating groups container");
}

void PMMG_exit_and_free( PMMG_pParMesh parmesh, const int val )
{
  PMMG_PMesh_Free( parmesh );
  MPI_Abort( MPI_COMM_WORLD, val );
  MPI_Finalize();
  exit( val );
}

/** return:
 *   PMMG_SUCCESS
 *   PMMG_FAILURE Failed to allocate pmmg struct or initialize mesh/sol
 *
 *   Allocate the main pmmg struct with a listgrp of one mesh and initialize
 *   maximum available memory
 */
int PMMG_Init_parMesh( PMMG_pParMesh *parmesh )
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

  PMMG_PMesh_SetMemGloMax( *parmesh, 0 );

  return PMMG_SUCCESS;

fail_mesh:
    PMMG_DEL_MEM(*parmesh,(*parmesh)->listgrp,1,PMMG_Grp,"deallocating groups container");
fail_grplst:
  (*parmesh)->ngrp = 0;
  (*parmesh)->memMax = 0;
  (*parmesh)->memCur = 0;
   free( *parmesh );
   *parmesh = NULL;
fail_pmesh:
  return PMMG_FAILURE;
}


/** return:
 *
 *  Sets the maximum amount of memory that is available to a parmmg process to
 *  memReq Mb.
 *  If memReq is zero then it is set to half of the available memory physically
 *  available on the machine. On multicore machines the available memory is
 *  shared equally to pmmg processes. If memReq is negative or more than the
 *  detected available memory, then again it is set to the detected available
 *  memory
 */
void PMMG_PMesh_SetMemGloMax( PMMG_pParMesh parmesh, long long int memReq )
{
  long long int maxAvail = 0;
  MPI_Comm comm_shm = 0;
  int size_shm = 1;
  const int million = 1024 * 1024;

  assert ( (parmesh != NULL) && "trying to set glo max mem in empty parmesh" );
  MPI_Comm_split_type( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                       &comm_shm );
  MPI_Comm_size( comm_shm, &size_shm );

  maxAvail = _MMG5_memSize();
  // if detection failed => default value of _MMG5_MEMMAX Mo
  if ( maxAvail == 0 )
    maxAvail = _MMG5_MEMMAX << 20;

  // Multiple MPI processes may be running on the same node => distribute equally
  if ( (memReq > 0) && ((memReq * million) < maxAvail) )
    parmesh->memGloMax = (memReq * million) / size_shm;
  else
    parmesh->memGloMax = (maxAvail * 50) / (size_shm * 100);

  fprintf ( stdout,
            "Requested %lld Mb max memory usage. Max memory limit set to %lld Mb\n",
            memReq, parmesh->memGloMax/million );
}


/**
 * \param parmesh parmesh structure to adjust
 * \param percent integer value bewtween 0 and 100
 *
 * \return 1 if success, 0 if fail
 *
 * Set the maximum memory parmesh and the meshes in listgrp can use.
 * The total memory available is split between the parmesh structure and the
 * listgrp structures according to the percentage specified by the percent
 * input variable:
 *   percent % of the available mem is assigned to pmesh.memMax
 *   (100-percent)/100 are assigned to the mesh[i].memMax
 */
int PMMG_PMesh_SetMemMax( PMMG_pParMesh parmesh, int percent )
{
  MMG5_pMesh mesh;
  int        remaining_ngrps;
  long long  available;
  int        i = 0;

  assert ( (0 < percent) && (100 > percent) && "percent has to be >0 and <100" );

  parmesh->memMax = parmesh->memGloMax * percent / 100;
  available       = parmesh->memGloMax - parmesh->memMax;
  remaining_ngrps = parmesh->ngrp;
  for ( i = 0; i < parmesh->ngrp; ++i ) {
    mesh = parmesh->listgrp[i].mesh;
    mesh->memMax = available/remaining_ngrps;
    /* Not enough memory: set the minimal memory to be able to continue */
    if ( mesh->memMax < mesh->memCur ) {
      mesh->memMax = mesh->memCur;
    }
    available -= mesh->memMax;
    --remaining_ngrps;
    if ( available < 0 ) {
      fprintf(stderr,"\n  ## Error: %s: not enough memory\n",__func__);
      return 0;
    }
  }
  return 1;
}


int PMMG_Set_iparameter(PMMG_pParMesh parmesh, int iparam,int val){
  MMG5_pMesh  mesh;
  MMG5_pSol   met;
  int         k,mem;

  switch ( iparam ) {
    /* Integer parameters */
  case PMMG_IPARAM_verbose :
    for ( k=0; k<parmesh->ngrp; ++k ) {
      mesh = parmesh->listgrp[k].mesh;
      if ( !MMG3D_Set_iparameter( mesh, NULL, MMG3D_IPARAM_verbose, val ) )
        return 0;
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
  default :
    fprintf(stderr,"  ## Error: unknown type of parameter\n");
    return(0);
  }

  return(1);
}

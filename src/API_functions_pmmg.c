/**
 * \file API_functions_pmmg.c
 * \brief C API functions definitions for PARMMG library.
 * \copyright GNU Lesser General Public License.
 *
 * C API for PARMMG library.
 *
 */
#include "parmmg.h"


#warning NIKOS: on NULL input should it fail or simply return?
static void pmesh_int_comm_free( PMMG_pParMesh pmesh, PMMG_pint_comm comm )
{
  if ( comm == NULL )
    return;

  if ( NULL != comm->intvalues ) {
    assert ( comm->nitem != 0 && "incorrect parameters in internal communicator" );
    PMMG_FREE(pmesh,comm->intvalues,comm->nitem,int,"int comm int array");
  }
  if ( NULL != comm->doublevalues ) {
    assert ( comm->nitem != 0 && "incorrect parameters in internal communicator" );
    PMMG_FREE(pmesh,comm->doublevalues,comm->nitem,double,"int comm double array");
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
      PMMG_FREE(pmesh,comm->int_comm_index,comm->nitem,int,"ext comm int array");
    }
    if ( NULL != comm->itosend ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_FREE(pmesh,comm->itosend,comm->nitem,int,"ext comm itosend array");
    }
    if ( NULL != comm->itorecv ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_FREE(pmesh,comm->itorecv,comm->nitem,int,"ext comm itorecv array");
    }
    if ( NULL != comm->rtosend ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_FREE(pmesh,comm->rtosend,comm->nitem,int,"ext comm rtosend array");
    }
    if ( NULL != comm->rtorecv ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_FREE(pmesh,comm->rtorecv,comm->nitem,int,"ext comm rtorecv array");
    }
  }
}

void PMMG_PMesh_Free( PMMG_pParMesh parmesh )
{
  int k = 0;

#warning NIKOS: does this deallocate mesh->point,mesh->tetra,mesh->tria,esh->edge?
  for ( k = 0; k < parmesh->ngrp; ++k ) {
    MMG3D_Free_all( MMG5_ARG_start,
                    MMG5_ARG_ppMesh, &parmesh->listgrp[k].mesh,
                    MMG5_ARG_ppMet, &parmesh->listgrp[k].met,
                    MMG5_ARG_end );
  }

  pmesh_int_comm_free( parmesh, parmesh->int_node_comm );
  pmesh_int_comm_free( parmesh, parmesh->int_edge_comm );
  pmesh_int_comm_free( parmesh, parmesh->int_face_comm );

  pmesh_ext_comm_free( parmesh, parmesh->ext_node_comm, parmesh->next_node_comm );
  PMMG_FREE(parmesh, parmesh->ext_node_comm, parmesh->next_node_comm,
            PMMG_ext_comm, "ext node comm");
  pmesh_ext_comm_free( parmesh, parmesh->ext_edge_comm, parmesh->next_edge_comm );
  PMMG_FREE(parmesh, parmesh->ext_edge_comm, parmesh->next_edge_comm,
            PMMG_ext_comm, "ext edge comm");
  pmesh_ext_comm_free( parmesh, parmesh->ext_face_comm, parmesh->next_face_comm );
  PMMG_FREE(parmesh, parmesh->ext_face_comm, parmesh->next_face_comm,
            PMMG_ext_comm, "ext face comm");
  PMMG_FREE(parmesh,parmesh->listgrp,1,PMMG_Grp,"deallocating groups container");
}

void PMMG_exit_and_free( PMMG_pParMesh parmesh, const int val )
{
  PMMG_PMesh_Free( parmesh );
  MPI_Finalize();
  exit( val );
}

// return:
//   PMMG_SUCCESS
//   PMMG_FAILURE Failed to allocate structs
int PMMG_Init_parMesh( PMMG_pParMesh *parmesh )
{
  PMMG_pGrp grp = NULL;

  /* ParMesh allocation */
  assert ( (*parmesh == NULL) && "trying to initialize non empty parmesh" );
  *parmesh = calloc( 1, sizeof(PMMG_ParMesh) );
  if ( *parmesh == NULL )
    goto fail_pmesh;
  (*parmesh)->memGloMax = 4 * 1024L * 1024L; // Assign a number to begin with
  (*parmesh)->memMax = 4 * 1024L * 1024L; // Assign a number to begin with
  (*parmesh)->memCur = sizeof(PMMG_ParMesh); // Update memory usage

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

  return PMMG_SUCCESS;

fail_mesh:
    PMMG_FREE(*parmesh,(*parmesh)->listgrp,1,PMMG_Grp,"deallocating groups container");
fail_grplst:
  (*parmesh)->ngrp = 0;
  (*parmesh)->memMax = 0;
  (*parmesh)->memCur = 0;
   free( *parmesh );
   *parmesh = NULL;
fail_pmesh:
  return PMMG_FAILURE;
}


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

  // If total physical memory is known, use 50% of it.
  // Otherwise try use a default value of _MMG5_MEMMAX Mo
  // Since multiple MPI jobs may be running on the same node, adjust by size_shm
  maxAvail = _MMG5_memSize() / size_shm;
  maxAvail = maxAvail * 50 / 100;
  if ( maxAvail == 0 )
    maxAvail = ( _MMG5_MEMMAX << 20 ) / size_shm;

  if ( memReq > 0 && memReq < maxAvail )
    parmesh->memGloMax = memReq / size_shm;
  else
    parmesh->memGloMax = maxAvail;

  fprintf ( stdout,
            "Requested %lld Mb max memory usage. Max memory limit set to %lld \n",
            memReq, parmesh->memGloMax );
}


/** Calculate the maximum available memory a parmesh or a mesh can use.
 *  This is equal to the global max memory minus the sum of parmesh/mesh's memCur
 *  except thy memCur of the particular parmesh or mesh struct of which memMax
 *  we are updating
 */
long long int PMMG_PMesh_SetMemMax( PMMG_pParMesh parmesh, long long int memCur )
{
  int i = 0;
  long long int total = 0;

  // Total mem usage of all meshes in listgrp
  for ( i = 0; i < parmesh->ngrp; ++i )
    total += parmesh->listgrp[i].mesh->memCur;

  total += parmesh->memCur;

  // Minus the current struct's current mem usage
  total -= memCur;

  return parmesh->memGloMax - total;
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

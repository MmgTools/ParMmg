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
  MPI_Comm comm;

  comm = parmesh->comm;

  if ( !PMMG_Free_all( &parmesh ) ) {
    fprintf(stderr,"  ## Warning: unable to clean the parmmg memory.\n"
            " Possible memory leak.\n");
  }

  if ( val != PMMG_SUCCESS )
    MPI_Abort( parmesh->comm, val );

  MPI_Finalize();
  exit( val );
}

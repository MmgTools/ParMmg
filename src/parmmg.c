/**
 * \file parmmg.c
 * \brief main file for the parmmg application
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "parmmg.h"

mytime         PMMG_ctim[TIMEMAX];

/**
 *
 * Print elapsed time at end of process.
 *
 */
static void PMMG_endcod() {
  char   stim[32];

  chrono(OFF,&PMMG_ctim[0]);
  printim(PMMG_ctim[0].gdif,stim);
  fprintf(stdout,"\n   ELAPSED TIME  %s\n",stim);
}

/**
 * \param  parmesh pointer to parmesh structure
 * \param  rank MPI process rank
 *
 * \return PMMG_LOWFAILURE    if fail before mesh scaling
 *         PMMG_STRONGFAILURE if fail after mesh scaling,
 *         PMMG_SUCCESS
 *
 * Mesh preprocessing: set function pointers, scale mesh, perform mesh
 * analysis and display length and quality histos.
 */
static int PMMG_preprocessMesh( PMMG_pParMesh parmesh, int rank )
{
  MMG5_pMesh mesh;
  MMG5_pSol  met;

  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  assert ( ( mesh != NULL ) && ( met != NULL ) && "Preprocessing empty args");

  if ( parmesh->info.imprim )
    fprintf(stdout,"\n   -- PHASE 2 : ANALYSIS\n");

  /** Function setters (must be assigned before quality computation) */
  _MMG3D_Set_commonFunc();

  /** Mesh scaling and quality histogram */
  if ( !_MMG5_scaleMesh(mesh,met) )
    return PMMG_LOWFAILURE;

  /** specific meshing */
  if ( mesh->info.optim && !met->np ) {
    if ( !MMG3D_doSol(mesh,met) )
      return PMMG_STRONGFAILURE;

    _MMG3D_solTruncatureForOptim(mesh,met);
  }
  if ( mesh->info.hsiz > 0. ) {
    if ( !MMG3D_Set_constantSize(mesh,met) ) {
     if ( !_MMG5_unscaleMesh(mesh,met) )
       return PMMG_STRONGFAILURE;
     return PMMG_STRONGFAILURE;
    }
  }

  MMG3D_setfunc(mesh,met);

  if ( !_MMG3D_tetraQual( mesh, met, 0 ) )
    return PMMG_STRONGFAILURE;

  if ( abs(parmesh->info.imprim) > 0 )
    if ( !_MMG3D_inqua(mesh,met) )
      return PMMG_STRONGFAILURE;

  /** Mesh analysis */
  if ( !_MMG3D_analys(mesh) )
    return PMMG_STRONGFAILURE;

  if ( !rank )
    if ( parmesh->info.imprim > 1 && met->m )
      _MMG3D_prilen(mesh,met,0);

  /** Mesh unscaling */
  if ( !_MMG5_unscaleMesh(mesh,met) )
    return PMMG_STRONGFAILURE;

  if ( parmesh->info.imprim )
    fprintf(stdout,"   -- PHASE 2 COMPLETED.\n");

  return PMMG_SUCCESS;
}

/**
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \return \ref PMMG_SUCCESS if success.
 *         \ref PMMG_LOWFAILURE if failed but a conform mesh is saved.
 *         \ref PMMG_STRONGFAILURE if failed and we can't save the mesh.
 *
 * Main program for PARMMG executable: perform parallel mesh adaptation.
 *
 */
int main( int argc, char *argv[] )
{
  PMMG_pParMesh parmesh = NULL;
  MMG5_pMesh    mesh = NULL;
  MMG5_pSol     met = NULL;
  int           i,rank = 0;
  int           ier = 0;
  int           iresult;
  long int      tmpmem;
  char          stim[32];

  // Shared memory communicator: processes that are on the same node, sharing
  //    local memory and can potentially communicate without using the network
  MPI_Comm comm_shm = 0;
  int      rank_shm = 0;
  int      size_shm = 1;

  /** Init MPI */
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if ( !rank ) {
    fprintf(stdout,"  -- PARMMG3d, Release %s (%s) \n",PMMG_VER,PMMG_REL);
    fprintf(stdout,"     %s\n",PMMG_CPY);
    fprintf(stdout,"     %s %s\n\n",__DATE__,__TIME__);

    fprintf(stdout,"  -- MMG3d,    Release %s (%s) \n",MG_VER,MG_REL);
    fprintf(stdout,"     %s\n",MG_CPY);
  }

  if ( !rank ) {
    atexit(PMMG_endcod);
  }

  tminit(PMMG_ctim,TIMEMAX);
  chrono(ON,&PMMG_ctim[0]);

  /* Allocate the main pmmg struct and assign default values */
  if ( 1 != PMMG_Init_parMesh( &parmesh,MPI_COMM_WORLD ) ) {
    MPI_Abort( MPI_COMM_WORLD, PMMG_STRONGFAILURE );
    MPI_Finalize();
    return PMMG_FAILURE;
  }


  /* reset default values for file names */
  if ( 1 != MMG3D_Free_names(MMG5_ARG_start,
                             MMG5_ARG_ppMesh, &parmesh->listgrp[0].mesh,
                             MMG5_ARG_ppMet,  &parmesh->listgrp[0].met,
                             MMG5_ARG_end) )
    PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );

  /* Init memMax sizes. Only one mesh for now => pmmg structs do not need much */
  if ( !PMMG_parmesh_SetMemMax(parmesh, 20) )
    PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );

  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  if ( 1 != PMMG_parsar( argc, argv, parmesh ) )
    PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );

  if ( parmesh->ddebug ) {
    MPI_Comm_split_type( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                         &comm_shm );
    MPI_Comm_rank( comm_shm, &rank_shm );
    MPI_Comm_size( comm_shm, &size_shm );

    if ( !rank_shm )
      printf("\n     %d MPI PROCESSES (%d ON LOCAL NODE):\n",parmesh->nprocs,
             size_shm);

    printf("         MPI RANK %d (LOCAL RANK %d)\n", parmesh->myrank,rank_shm );
  }

  /* load data */
  chrono(ON,&PMMG_ctim[1]);
  if ( !rank ) {
    fprintf(stdout,"\n  -- INPUT DATA: LOADING MESH ON RANK 0\n");
  }

  iresult = 1;
  ier = 1;

  if ( !parmesh->myrank ) {
    if ( MMG3D_loadMesh( mesh, mesh->namein ) != 1 ) {
      ier = 0;
      goto check_mesh_loading;
    }
    // SetMemMax is used again here because loadMesh overwrites memMax
    if ( 1 != PMMG_parmesh_SetMemMax(parmesh, 20) ) {
      ier = 0;
      goto check_mesh_loading;
    }

    if ( MMG3D_loadSol( mesh, met, met->namein ) == -1 ) {
      if ( parmesh->info.imprim ) {
        fprintf(stderr,"\n  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
      }
      ier = 0;
      goto check_mesh_loading;
    }
  }

  chrono(OFF,&PMMG_ctim[1]);
  if ( !rank ) {
    printim(PMMG_ctim[1].gdif,stim);
    fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);
  }

check_mesh_loading:
  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( iresult != 1 )
    PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );

  /** Check input data */
  if ( PMMG_SUCCESS != PMMG_check_inputData( parmesh ) )
    PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );

  chrono(ON,&PMMG_ctim[2]);
  if ( parmesh->info.imprim ) {
    fprintf(stdout,"\n  %s\n   MODULE PARMMG3D: IMB-LJLL : %s (%s)\n  %s\n",
            PMMG_STR,PMMG_VER,PMMG_REL,PMMG_STR);

    fprintf(stdout,"\n   -- PHASE 1 : DISTRIBUTE MESH AMONG PROCESSES\n");
  }

  /** Send mesh to other procs */
  if ( 1 != PMMG_bcast_mesh( parmesh ) )
    PMMG_exit_and_free( parmesh,PMMG_STRONGFAILURE );

  /** Mesh preprocessing: set function pointers, scale mesh, perform mesh
   * analysis and display length and quality histos. */
  ier = PMMG_preprocessMesh( parmesh, parmesh->myrank );
  if ( ier != PMMG_SUCCESS ) {
    if ( ( ier == PMMG_LOWFAILURE ) || !(_MMG5_unscaleMesh( mesh, met )) )
      PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );
    return PMMG_LOWFAILURE;
  }

  /** Send mesh partionning to other procs */
  if ( !PMMG_distribute_mesh( parmesh ) ) {
    if ( 1 != _MMG5_unscaleMesh( mesh, met ) )
      PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );
    PMMG_exit_and_free( parmesh, PMMG_LOWFAILURE );
  }

  chrono(OFF,&PMMG_ctim[2]);
  if ( parmesh->info.imprim ) {
    printim(PMMG_ctim[2].gdif,stim);
    fprintf(stdout,"   -- PHASE 1 COMPLETED.     %s\n",stim);
  }

  /** Remeshing */
  chrono(ON,&PMMG_ctim[3]);
  if ( parmesh->info.imprim ) {
    fprintf( stdout,"\n  -- PHASE 3 : %s MESHING\n",
             met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC" );
  }

  ier = PMMG_parmmglib1(parmesh);

  chrono(OFF,&PMMG_ctim[3]);
  if ( ier == PMMG_STRONGFAILURE ) {
    PMMG_exit_and_free( parmesh, ier );
  } else if ( ier == PMMG_LOWFAILURE ) {
#warning SAVE THE MESH?
    PMMG_exit_and_free( parmesh, PMMG_SUCCESS );
  } else {
    if ( parmesh->info.imprim ) {
      printim(PMMG_ctim[3].gdif,stim);
      fprintf(stdout,"  -- PHASE 3 COMPLETED.     %s\n",stim);
      fprintf(stdout,"\n  %s\n   END OF MODULE PARMMG3d: IMB-LJLL \n  %s\n",
            PMMG_STR,PMMG_STR);
    }
  }

  /** Unscaling */
  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;
  if ( 1 != _MMG5_unscaleMesh( mesh, met ) )
    PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );

  /** Merge all the meshes on the proc 0 */
  chrono(ON,&PMMG_ctim[4]);
  if ( parmesh->info.imprim ) {
    fprintf( stdout,"\n   -- PHASE 4 : MERGE MESHES OVER PROCESSORS\n" );
  }

  if ( !PMMG_merge_parmesh( parmesh ) )
    PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );

  chrono(OFF,&PMMG_ctim[4]);
  if ( parmesh->info.imprim ) {
    printim(PMMG_ctim[4].gdif,stim);
    fprintf( stdout,"   -- PHASE 4 COMPLETED.     %s\n",stim );
  }

  if ( !parmesh->myrank ) {
    chrono(ON,&PMMG_ctim[5]);
    if (  parmesh->info.imprim ) {
      fprintf( stdout,"\n   -- PHASE 5 : MESH PACKED UP\n" );
    }

    /** All the memory is devoted to the mesh **/
    PMMG_parmesh_Free_Comm(parmesh);
    tmpmem = parmesh->memMax - parmesh->memCur;
    parmesh->memMax = parmesh->memCur;
    parmesh->listgrp[0].mesh->memMax += tmpmem;

    if ( MMG3D_hashTetra( mesh, 0 ) ) {
      if ( -1 == MMG3D_bdryBuild( mesh ) ) {
        PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );
      }
    } else {
      /** Impossible to rebuild the triangle **/
      fprintf(stdout,"\n\n\n  -- IMPOSSIBLE TO SAVE THE BOUNDARY TRIANGLES\n\n\n");
    }

    chrono(OFF,&PMMG_ctim[5]);
    if (  parmesh->info.imprim ) {
      printim(PMMG_ctim[5].gdif,stim);
      fprintf( stdout,"   -- PHASE 5 COMPLETED.     %s\n",stim );
    }

    /* Write mesh */
    chrono(ON,&PMMG_ctim[6]);
    if ( parmesh->info.imprim )
      fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",mesh->nameout);

    if ( 1 != MMG3D_saveMesh( mesh, mesh->nameout ) )
      PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );
    if ( met->m && 1 != MMG3D_saveSol( mesh, met, met->nameout ) )
      PMMG_exit_and_free( parmesh, PMMG_LOWFAILURE );

   chrono(OFF,&PMMG_ctim[6]);
    if ( parmesh->info.imprim )
      fprintf(stdout,"  -- WRITING COMPLETED\n");

  }

  PMMG_exit_and_free( parmesh, ier );
}

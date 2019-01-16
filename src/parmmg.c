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
  PMMG_pGrp     grp;
  int           rank;
  int           ier,iresult,ierSave,msh;
  char          stim[32];

  // Shared memory communicator: processes that are on the same node, sharing
  //    local memory and can potentially communicate without using the network
  MPI_Comm comm_shm = 0;
  int      rank_shm = 0;
  int      size_shm = 1;

  /** Initializations: MPI, mesh, and memory */
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
  if ( 1 != PMMG_Init_parMesh( PMMG_ARG_start,
                               PMMG_ARG_ppParMesh,&parmesh,
                               PMMG_ARG_dim,3,
                               PMMG_ARG_MPIComm,MPI_COMM_WORLD,
                               PMMG_ARG_end) ) {
    MPI_Abort( MPI_COMM_WORLD, PMMG_STRONGFAILURE );
    MPI_Finalize();
    return PMMG_FAILURE;
  }


  /* reset default values for file names */
  if ( 1 != MMG3D_Free_names(MMG5_ARG_start,
                             MMG5_ARG_ppMesh, &parmesh->listgrp[0].mesh,
                             MMG5_ARG_ppMet,  &parmesh->listgrp[0].met,
                             MMG5_ARG_end) )
    PMMG_RETURN_AND_FREE( parmesh, PMMG_STRONGFAILURE );

  /* Init memMax sizes. Only one mesh for now => pmmg structs do not need much */
  if ( !PMMG_parmesh_SetMemMax(parmesh, 20) )
    PMMG_RETURN_AND_FREE( parmesh, PMMG_STRONGFAILURE );

  if ( 1 != PMMG_parsar( argc, argv, parmesh ) )
    PMMG_RETURN_AND_FREE( parmesh, PMMG_STRONGFAILURE );

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

  /** load data */
  chrono(ON,&PMMG_ctim[1]);
  if ( rank==parmesh->info.root ) {
    fprintf(stdout,"\n  -- INPUT DATA: LOADING MESH ON RANK %d\n",
            parmesh->info.root);
  }

  grp = &parmesh->listgrp[0];
  iresult = 1;
  ier     = 1;

  if ( 1 != PMMG_loadMesh_centralized(parmesh,grp->mesh->namein) ) {
    ier = 0;
    goto check_mesh_loading;
  }

  if ( -1 == PMMG_loadMet_centralized( parmesh, grp->met->namein ) ) {
    if ( parmesh->info.imprim ) {
      fprintf(stderr,"\n  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
    }
    ier = 0;
    goto check_mesh_loading;
  }

  if ( !PMMG_parsop(parmesh) )
    PMMG_RETURN_AND_FREE( parmesh, PMMG_STRONGFAILURE );

  chrono(OFF,&PMMG_ctim[1]);
  if ( rank==parmesh->info.root ) {
    printim(PMMG_ctim[1].gdif,stim);
    fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);
  }

check_mesh_loading:
  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( iresult != 1 )
    PMMG_RETURN_AND_FREE( parmesh, PMMG_LOWFAILURE );

  /** Main call */
  if ( parmesh->listgrp[0].mesh->mark ) {
    /* Save a local parameters file containing the default parameters */
    printf("  ## Error: default parameter file saving not yet implemented.\n");
    ier = 2;//PMMG_defaultOption(grp->mesh,grp->met);
    PMMG_RETURN_AND_FREE(parmesh,ier);
  }
  else {
    /* Parallel remeshing */
    ier = PMMG_parmmglib_centralized(parmesh);
  }

  /** Check result and save output files */
  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MAX, parmesh->comm );
  if ( iresult == PMMG_STRONGFAILURE ) { PMMG_RETURN_AND_FREE( parmesh, ier ); }

  chrono(ON,&PMMG_ctim[1]);

  grp = &parmesh->listgrp[0];
  if ( parmesh->info.imprim ) {
    fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",grp->mesh->nameout);
  }

  if ( parmesh->listgrp && parmesh->listgrp[0].mesh ) {
    grp = &parmesh->listgrp[0];
    MMG5_chooseOutputFormat(grp->mesh,&msh);

    if ( !msh ) {
      ierSave = PMMG_saveMesh_centralized(parmesh,grp->mesh->nameout);
    }
    else {
      printf("  ## Error: Gmsh output format not yet implemented.\n");
      ierSave = 2;//PMMG_saveMshMesh_centralized(grp->mesh,grp->met,mesh->nameout);
    }

    if ( !ierSave ) {
      PMMG_RETURN_AND_FREE(parmesh,PMMG_STRONGFAILURE);
    }

    if ( !msh ) {
      if ( !PMMG_saveMet_centralized(parmesh,grp->mesh->nameout) ) {
        PMMG_RETURN_AND_FREE(parmesh,PMMG_STRONGFAILURE);
      }

      if ( grp->sol && !PMMG_saveAllSols_centralized(parmesh,grp->sol->nameout) ) {
        PMMG_RETURN_AND_FREE(parmesh,PMMG_STRONGFAILURE);
      }
    }
  }

  chrono(OFF,&PMMG_ctim[1]);
  if ( parmesh->info.imprim )
    fprintf(stdout,"  -- WRITING COMPLETED\n");

  PMMG_RETURN_AND_FREE( parmesh, iresult );
}

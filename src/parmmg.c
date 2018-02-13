#include "parmmg.h"

/**
 * \param  mesh pointer to mesh structure
 * \param  met
 * \param  rank MPI process rank
 *
 * \return PMMG_LOWFAILURE    if fail before mesh scaling
 *         PMMG_STRONGFAILURE if fail after mesh scaling,
 *         PMMG_SUCCESS
 *
 * Mesh preprocessing: set function pointers, scale mesh, perform mesh
 * analysis and display length and quality histos.
 */
static int PMMG_preprocessMesh( MMG5_pMesh mesh, MMG5_pSol met, int rank )
{
  assert ( ( mesh != NULL ) && ( met != NULL ) && "Preprocessing empty args");

  if ( !rank && mesh->info.imprim )
    fprintf(stdout,"\n   -- PHASE 2 : ANALYSIS\n");

  /** Function setters (must be assigned before quality computation) */
  _MMG3D_Set_commonFunc();

  /** Mesh scaling and quality histogram */
#warning Do we need to scale here (for the analysis step) ?
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

  if ( (!rank) && abs(mesh->info.imprim) > 0 )
    if ( !_MMG3D_inqua(mesh,met) )
      return PMMG_STRONGFAILURE;

  /** Mesh analysis */
  if ( !_MMG3D_analys(mesh) )
    return PMMG_STRONGFAILURE;

  if ( !rank )
    if ( mesh->info.imprim > 1 && met->m )
      _MMG3D_prilen(mesh,met,0);

  if ( !rank && mesh->info.imprim )
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
  int           rank = 0;
  int           ier = 0;
  long int      tmpmem;

  // Shared memory communicator: processes that are on the same node, sharing
  //    local memory and can potentially communicate without using the network
  MPI_Comm comm_shm = 0;
  int rank_shm = 0;
  int size_shm = 1;

  /** Init MPI */
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if ( !rank ) {
    fprintf(stdout,"  -- PARMMG3d, Release %s (%s) \n",PMMG_VER,PMMG_REL);
    fprintf(stdout,"  -- MMG3d,    Release %s (%s) \n",MG_VER,MG_REL);
    fprintf(stdout,"     %s\n",PMMG_CPY);
    fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);
  }

  /* Allocate the main pmmg struct and assign default values */
  if ( PMMG_SUCCESS != PMMG_Init_parMesh( &parmesh ) ) {
    MPI_Abort( MPI_COMM_WORLD, PMMG_FAILURE );
    MPI_Finalize();
    return PMMG_FAILURE;
  }
  parmesh->niter = 1;
  parmesh->comm = MPI_COMM_WORLD;
  parmesh->myrank = rank;
  MPI_Comm_size( parmesh->comm, &parmesh->nprocs );

  /* reset default values for file names */
  if ( 1 != MMG3D_Free_names(MMG5_ARG_start,
                             MMG5_ARG_ppMesh, &parmesh->listgrp[0].mesh,
                             MMG5_ARG_ppMet,  &parmesh->listgrp[0].met,
                             MMG5_ARG_end) )
    PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );

  // Init memMax sizes. Only one mesh for now => pmmg structs do not need much
  if ( !PMMG_parmesh_SetMemMax(parmesh, 20) )
    PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );

  MPI_Comm_split_type( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
    &comm_shm );
  MPI_Comm_rank( comm_shm, &rank_shm );
  MPI_Comm_size( comm_shm, &size_shm );
  printf(" ++++NIKOS: %d/%d shared memory rank %d of %d \n\n",
    parmesh->nprocs, parmesh->myrank, rank_shm, size_shm );

  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  if ( PMMG_SUCCESS != PMMG_parsar( argc, argv, parmesh ) )
    PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );

  if ( !parmesh->myrank && mesh->info.imprim )
    fprintf(stdout,"\n   -- PHASE 0 : LOADING MESH ON rank 0\n");

  if ( !parmesh->myrank ) {
    if ( 1 != MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_verbose,5) )
      PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );

    if ( MMG3D_loadMesh( mesh, mesh->namein ) != 1 )
      PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );
    // SetMemMax is used again here because loadMesh overwrites memMax
    PMMG_parmesh_SetMemMax(parmesh, 20);

    if ( MMG3D_loadSol( mesh, met, met->namein ) == -1 )
      PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );
  } else {
    if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_verbose,0) )
      PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );
  }

  /** Check input data */
  if ( PMMG_SUCCESS != PMMG_check_inputData( parmesh ) )
    PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );

  if ( !parmesh->myrank && mesh->info.imprim )
    fprintf(stdout,"\n   -- PHASE 1 : DISTRIBUTE MESH AMONG PROCESSES\n");

  /** Send mesh to other procs */
  if ( PMMG_SUCCESS != PMMG_bcast_mesh( parmesh ) )
    PMMG_exit_and_free( parmesh,PMMG_STRONGFAILURE );

  /** Mesh preprocessing: set function pointers, scale mesh, perform mesh
   * analysis and display length and quality histos. */
#warning NIKOS: I think I have messed up this error handling
  ier = PMMG_preprocessMesh( mesh, met, parmesh->myrank );
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
  if ( (!parmesh->myrank) && parmesh->listgrp[0].mesh->info.imprim )
    fprintf(stdout,"   -- PHASE 1 COMPLETED.\n");

  /** Remeshing */
  if ( (!parmesh->myrank) && parmesh->listgrp[0].mesh->info.imprim )
    fprintf( stdout,
             "\n  -- PHASE 3 : %s MESHING\n",
             met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC" );

  ier = PMMG_parmmglib1(parmesh);
  if ( ier == PMMG_STRONGFAILURE )
    PMMG_exit_and_free( parmesh, ier );
  else if ( (!parmesh->myrank) && parmesh->listgrp[0].mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 3 COMPLETED.\n");

  /** Unscaling */
  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;
  if ( 1 != _MMG5_unscaleMesh( mesh, met ) )
    PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );

  /** Merge all the meshes on the proc 0 */
  if ( !parmesh->myrank && mesh->info.imprim )
    fprintf( stdout,"\n   -- PHASE 4 : MERGE MESHES OVER PROCESSORS\n" );

#warning NIKOS: this is the only function that hasnt been revised in regards to error handling/returned values.Lots of unaccounted allocations
  if ( !PMMG_merge_parmesh( parmesh ) )
    PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );

  if ( !parmesh->myrank && mesh->info.imprim )
    fprintf( stdout,"   -- PHASE 4 COMPLETED.\n" );

  if ( !parmesh->myrank ) {
    if (  mesh->info.imprim )
      fprintf( stdout,"\n   -- PHASE 5 : MESH PACKED UP\n" );

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
      fprintf(stdout,"\n\n\n  -- IMPOSSIBLE TO SAVE THE BDRY TRIANGLE\n\n\n");
    }
    
    if (  mesh->info.imprim )
      fprintf( stdout,"   -- PHASE 5 COMPLETED.\n" );

    /* Write mesh */
    if ( 1 != MMG3D_saveMesh( mesh, mesh->nameout ) )
      PMMG_exit_and_free( parmesh, PMMG_STRONGFAILURE );
    if ( 1 != MMG3D_saveSol( mesh, met, met->nameout ) )
      PMMG_exit_and_free( parmesh, PMMG_LOWFAILURE );
  }

  PMMG_exit_and_free( parmesh, ier );
}

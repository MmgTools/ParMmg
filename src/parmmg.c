#include "parmmg.h"

/**
 * \param parmesh pointer toward the mesh structure.
 * \return -1 if fail before mesh scaling, 0 if fail after mesh scaling,
 * 1 if success.
 *
 * Mesh preprocessing: set function pointers, scale mesh, perform mesh
 * analysis and display length and quality histos.
 *
 */
static inline
int PMMG_preprocessMesh(PMMG_pParMesh parmesh) {
  MMG5_pMesh      mesh;
  MMG5_pSol       met;

  mesh   = parmesh->listgrp[0].mesh;
  met    = parmesh->listgrp[0].met;

  if ( !parmesh->myrank && mesh->info.imprim )
    fprintf(stdout,"\n   -- PHASE 2 : ANALYSIS\n");

  /** Function setters (must be assigned before quality computation) */
  _MMG3D_Set_commonFunc();
  MMG3D_setfunc(mesh,met);

  /** Mesh scaling and quality histogram */
#warning Do we need to scale here (for the analysis step) ?
  if ( !_MMG5_scaleMesh(mesh,met) ) return -1;

  if ( !_MMG3D_tetraQual( mesh, met, 0 ) ) return 0;

  if ( (!parmesh->myrank) && abs(mesh->info.imprim) > 0 ) {
    if ( !_MMG3D_inqua(mesh,met) ) return 0;
  }

  /** specific meshing */
  if ( mesh->info.optim && !met->np ) {
    if ( !MMG3D_doSol(mesh,met) ) return 0;

    _MMG3D_scalarSolTruncature(mesh,met);
  }

  /** Mesh analysis */
  mesh->info.hausd = 0.002;
  if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_nosurf,0 ) )
    return PMMG_STRONGFAILURE;

  if ( !_MMG3D_analys(mesh) ) return 0;

  if ( !parmesh->myrank )
    if ( mesh->info.imprim > 1 && met->m ) _MMG3D_prilen(mesh,met,0);

  if ( !parmesh->myrank && mesh->info.imprim )
    fprintf(stdout,"   -- PHASE 2 COMPLETED.\n");

  return 1;
}

/**
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \return \ref PMMG_SUCCESS if success.
 * \return \ref PMMG_LOWFAILURE if failed but a conform mesh is saved.
 * \return \ref PMMG_STRONGFAILURE if failed and we can't save the mesh.
 *
 * Main program for PARMMG executable: perform parallel mesh adaptation.
 *
 */
int main(int argc,char *argv[]) {
  PMMG_pParMesh    parmesh=NULL;
  PMMG_pGrp        grp;
  MMG5_pMesh       mesh;
  MMG5_pSol        met;
  int              ier,rank;

  /** Init MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if ( !rank ) {
    fprintf(stdout,"  -- PARMMG3d, Release %s (%s) \n",PMMG_VER,PMMG_REL);
    fprintf(stdout,"  -- MMG3d,    Release %s (%s) \n",MG_VER,MG_REL);
    fprintf(stdout,"     %s\n",PMMG_CPY);
    fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);
  }

  /** Assign default values */
  if ( !PMMG_Init_parMesh(PMMG_ARG_start,
                          PMMG_ARG_ppParMesh,&parmesh,
                          PMMG_ARG_end) )
    _PMMG_RETURN_AND_FREE(parmesh,PMMG_STRONGFAILURE);

  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  if ( !parmesh->myrank ) {
    if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_verbose,5) )
      _PMMG_RETURN_AND_FREE(parmesh, PMMG_STRONGFAILURE);

    /** Read sequential mesh */
    char fmesh[ 20 + 6 ] = "m";
    char fmet[ 20 + 5 ] = "m";
    if ( argc > 1 ) {
      strncpy ( fmesh, argv[1], 20 );
      strncpy ( fmet, argv[1], 20 );
    }
    strncat ( fmesh, ".mesh", 6 );
    strncat ( fmet, ".sol", 5 );
    if ( PMMG_loadMesh( parmesh, fmesh ) < 1 )
      _PMMG_RETURN_AND_FREE( parmesh, PMMG_STRONGFAILURE );
    if ( PMMG_loadSol( parmesh, fmet ) < 0 )
      _PMMG_RETURN_AND_FREE( parmesh, PMMG_STRONGFAILURE );
  } else {
    if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_verbose,0) )
      _PMMG_RETURN_AND_FREE(parmesh, PMMG_STRONGFAILURE);
  }

  /** Check input data */
  if ( !_PMMG_check_inputData(parmesh) )
    _PMMG_RETURN_AND_FREE(parmesh,PMMG_STRONGFAILURE);

  if ( !parmesh->myrank && mesh->info.imprim )
    fprintf(stdout,"\n   -- PHASE 1 : MESH DISTRIBUTION\n");

  /** Send mesh to other procs */
  if ( !PMMG_bcastMesh(parmesh) )
    _PMMG_RETURN_AND_FREE(parmesh,PMMG_STRONGFAILURE);

  /** Mesh preprocessing: set function pointers, scale mesh, perform mesh
   * analysis and display length and quality histos. */
  ier = PMMG_preprocessMesh(parmesh) ;
  if ( ier <= 0 ) {
    if ( ( ier == -1 ) || !(_MMG5_unscaleMesh( mesh, met )) )
      _PMMG_RETURN_AND_FREE( parmesh, PMMG_STRONGFAILURE );
    return PMMG_LOWFAILURE;
  }

  /** Send mesh partionning to other procs */
  if ( !PMMG_distributeMesh(parmesh) ) {
    if ( !_MMG5_unscaleMesh(mesh,met) )
      _PMMG_RETURN_AND_FREE(parmesh,PMMG_STRONGFAILURE);
    _PMMG_RETURN_AND_FREE(parmesh,PMMG_LOWFAILURE);
  }
  if ( !parmesh->myrank && mesh->info.imprim )
    fprintf(stdout,"   -- PHASE 1 COMPLETED.\n");

  /** Remeshing */
  if ( !parmesh->myrank && mesh->info.imprim )
    fprintf(stdout,"\n  -- PHASE 3 : %s MESHING\n",met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC");

  ier = _PMMG_parmmglib1(parmesh);

  if ( !parmesh->myrank && mesh->info.imprim )
    fprintf(stdout,"  -- PHASE 3 COMPLETED.\n");

  if ( ier!= PMMG_STRONGFAILURE ) {
    /** Unscaling */
    if ( !_MMG5_unscaleMesh(mesh,met) )
      _PMMG_RETURN_AND_FREE(parmesh,PMMG_STRONGFAILURE);

    /** Merge all the meshes on the proc 0 */
    if ( !parmesh->myrank && mesh->info.imprim )
      fprintf(stdout,"\n   -- PHASE 4 : MERGE MESH\n");

    if ( !PMMG_mergeParMesh(parmesh,0) )
      _PMMG_RETURN_AND_FREE(parmesh,PMMG_STRONGFAILURE);
    if ( !parmesh->myrank && mesh->info.imprim )
      fprintf(stdout,"   -- PHASE 4 COMPLETED.\n");

    if ( !parmesh->myrank ) {
      if (  mesh->info.imprim )
        fprintf(stdout,"\n   -- PHASE 5 : MESH PACKED UP\n");

      if ( !MMG3D_hashTetra(parmesh->listgrp[0].mesh,0) )
            _PMMG_RETURN_AND_FREE(parmesh,PMMG_STRONGFAILURE);

      if ( _MMG3D_bdryBuild(parmesh->listgrp[0].mesh) < 0 )
        _PMMG_RETURN_AND_FREE(parmesh,PMMG_STRONGFAILURE);

      if (  mesh->info.imprim )
        fprintf(stdout,"   -- PHASE 5 COMPLETED.\n");

      /* Write mesh */
#warning : for the moment, we only write a mesh named out.mesh
      if ( !PMMG_saveMesh(parmesh,"out.mesh") )
            _PMMG_RETURN_AND_FREE(parmesh,PMMG_STRONGFAILURE);
      if ( !PMMG_saveSol (parmesh,"out.sol" ) )
        _PMMG_RETURN_AND_FREE(parmesh,PMMG_LOWFAILURE);
    }
  }

  _PMMG_RETURN_AND_FREE(parmesh,ier);
}

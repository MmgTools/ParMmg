/**
 * \file libparmmg.c
 * \brief Wrapper for the parallel remeshing library.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (InriaSoft)
 * \author Nikos Pattakos (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Mmain library functions (parallel remeshing starting from centralized or
 * distributed data.
 *
 */

#include "parmmg.h"


/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Check the validity of the input mesh data (tetra orientation, solution
 * compatibility with respect to the provided mesh, Mmg options).
 *
 */
int PMMG_check_inputData(PMMG_pParMesh parmesh)
{
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  int        k;

  if ( parmesh->info.imprim > 0 )
    fprintf(stdout,"\n  -- PMMG: CHECK INPUT DATA\n");

  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    met  = parmesh->listgrp[k].met;

    /* Check options */
    if ( mesh->info.lag > -1 ) {
      fprintf(stderr,
              "  ## Error: lagrangian mode unavailable (MMG3D_IPARAM_lag):\n");
      return 0;
    } else if ( mesh->info.iso ) {
      fprintf(stderr,"  ## Error: level-set discretisation unavailable"
              " (MMG3D_IPARAM_iso):\n");
      return 0;
    } else if ( mesh->info.optimLES && met->size==6 ) {
      fprintf(stdout,"  ## Error: strong mesh optimization for LES methods"
              " unavailable (MMG3D_IPARAM_optimLES) with an anisotropic metric.\n");
      return 0;
    }
    /* specific meshing */
    if ( met->np ) {
      if ( mesh->info.optim ) {
        printf("\n  ## ERROR: MISMATCH OPTIONS: OPTIM OPTION CAN NOT BE USED"
               " WITH AN INPUT METRIC.\n");
        return 0;
      }

      if ( mesh->info.hsiz>0. ) {
        printf("\n  ## ERROR: MISMATCH OPTIONS: HSIZ OPTION CAN NOT BE USED"
               " WITH AN INPUT METRIC.\n");
        return 0;
      }
    }

    if ( mesh->info.optim &&  mesh->info.hsiz>0. ) {
      printf("\n  ## ERROR: MISMATCH OPTIONS: HSIZ AND OPTIM OPTIONS CAN NOT BE USED"
             " TOGETHER.\n");
      return 0;
    }

    /* load data */
    MMG5_warnOrientation(mesh);

    if ( met->np && (met->np != mesh->np) ) {
      fprintf(stdout,"  ## WARNING: WRONG METRIC NUMBER. IGNORED\n");
      MMG5_DEL_MEM(mesh,met->m);
      met->np = 0;
    } else if ( met->size!=1 && met->size!=6 ) {
      fprintf(stderr,"  ## ERROR: WRONG DATA TYPE.\n");
      return 0;
    }
  }

  return 1;
}

/**
 * \param  parmesh pointer to parmesh structure
 *
 * \return PMMG_SUCCESS if success, PMMG_LOWFAILURE if fail and return an
 * unscaled mesh, PMMG_STRONGFAILURE if fail and return a scaled mesh.
 *
 * Mesh preprocessing: set function pointers, scale mesh, perform mesh
 * analysis and display length and quality histos.
 */
static int PMMG_preprocessMesh( PMMG_pParMesh parmesh )
{
  MMG5_pMesh mesh;
  MMG5_pSol  met;

  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  assert ( ( mesh != NULL ) && ( met != NULL ) && "Preprocessing empty args");

  /** Function setters (must be assigned before quality computation) */
  MMG3D_Set_commonFunc();

  /** Mesh scaling and quality histogram */
  if ( !MMG5_scaleMesh(mesh,met) ) {
    return PMMG_LOWFAILURE;
  }

  /** specific meshing */
  if ( mesh->info.optim && !met->np ) {
    if ( !MMG3D_doSol(mesh,met) ) {
      return PMMG_STRONGFAILURE;
    }
    MMG3D_solTruncatureForOptim(mesh,met);
  }

  if ( mesh->info.hsiz > 0. ) {
    if ( !MMG3D_Set_constantSize(mesh,met) ) {
      return PMMG_STRONGFAILURE;
    }
  }

  MMG3D_setfunc(mesh,met);

  if ( !MMG3D_tetraQual( mesh, met, 0 ) ) {
    return PMMG_STRONGFAILURE;
  }

  if ( parmesh->info.imprim > 0  ||  parmesh->info.imprim < -1 ) {
    if ( !PMMG_qualhisto(parmesh,PMMG_INQUA) ) {
        return PMMG_STRONGFAILURE;
    }
  }

  /** Mesh analysis */
  if ( !MMG3D_analys(mesh) ) {
    return PMMG_STRONGFAILURE;
  }

  if ( parmesh->info.imprim > 4 && (!mesh->info.iso) && met->m ) {
    MMG3D_prilen(mesh,met,0);
  }

  /** Mesh unscaling */
  if ( !MMG5_unscaleMesh(mesh,met) ) {
    return PMMG_STRONGFAILURE;
  }

  return PMMG_SUCCESS;
}



int PMMG_parmmglib_centralized(PMMG_pParMesh parmesh) {
  MMG5_pMesh    mesh;
  MMG5_pSol     met;
  int           ier;
  int           iresult,ierlib;
  long int      tmpmem;
  mytime        ctim[TIMEMAX];
  char          stim[32];

 if ( parmesh->info.imprim >= 0 ) {
    fprintf(stdout,"\n  %s\n   MODULE PARMMGLIB_CENTRALIZED: IMB-LJLL : "
            "%s (%s)\n  %s\n",PMMG_STR,PMMG_VER,PMMG_REL,PMMG_STR);
  }

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /** Check input data */
  chrono(ON,&(ctim[1]));

  ier = PMMG_check_inputData( parmesh );
  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !iresult ) return PMMG_LOWFAILURE;

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( parmesh->info.imprim > 0 ) {
    fprintf(stdout,"  -- CHECK INPUT DATA COMPLETED.     %s\n",stim);
  }


  chrono(ON,&(ctim[2]));
  if ( parmesh->info.imprim > 0 ) {
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS AND MESH DISTRIBUTION\n");
  }

  /** Send mesh to other procs */
  if ( parmesh->info.imprim > 2 ) {
    chrono(ON,&(ctim[6]));
    fprintf(stdout,"\n  -- BCAST" );
  }
  ier = PMMG_bcast_mesh( parmesh );
  if ( ier!=1 ) return PMMG_LOWFAILURE;

  if ( parmesh->info.imprim > 2 ) {
    chrono(OFF,&(ctim[6]));
    printim(ctim[6].gdif,stim);
    fprintf(stdout,"\n  -- BCAST COMPLETED    %s\n",stim );
  }

  /** Mesh preprocessing: set function pointers, scale mesh, perform mesh
   * analysis and display length and quality histos. */
  if ( parmesh->info.imprim > 2 ) {
    chrono(ON,&(ctim[7]));
    fprintf(stdout,"\n  -- ANALYSIS" );
  }
  ier = PMMG_preprocessMesh( parmesh );
  if ( parmesh->info.imprim > 2 ) {
    chrono(OFF,&(ctim[7]));
    printim(ctim[7].gdif,stim);
    fprintf(stdout,"\n  -- ANALYSIS COMPLETED    %s\n",stim );
  }

  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;
  if ( (ier==PMMG_STRONGFAILURE) && MMG5_unscaleMesh( mesh, met ) ) {
    ier = PMMG_LOWFAILURE;
  }
  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MAX, parmesh->comm );
  if ( iresult!=PMMG_SUCCESS ) {
    return iresult;
  }

  /** Send mesh partionning to other procs */
  if ( parmesh->info.imprim > 2 ) {
    chrono(ON,&(ctim[8]));
    fprintf(stdout,"\n  -- PARTITIONING" );
  }
  if ( !PMMG_distribute_mesh( parmesh ) ) {
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_LOWFAILURE);
  }
  if ( parmesh->info.imprim > 2 ) {
    chrono(OFF,&(ctim[8]));
    printim(ctim[8].gdif,stim);
    fprintf(stdout,"\n  -- PARTITIONING COMPLETED    %s\n",stim );
  }

  chrono(OFF,&(ctim[2]));
  if ( parmesh->info.imprim > 0 ) {
    printim(ctim[2].gdif,stim);
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);
  }

  /** Remeshing */
  chrono(ON,&(ctim[3]));
  if ( parmesh->info.imprim > 0 ) {
    fprintf( stdout,"\n  -- PHASE 2 : %s MESHING\n",
             met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC" );
  }

  ier = PMMG_parmmglib1(parmesh);
  MPI_Allreduce( &ier, &ierlib, 1, MPI_INT, MPI_MAX, parmesh->comm );

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( parmesh->info.imprim > 0) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
  }
  if ( ierlib == PMMG_STRONGFAILURE ) {
    return ierlib;
  }

  /** Merge all the meshes on the proc 0 */
  chrono(ON,&(ctim[4]));
  if ( parmesh->info.imprim > 0 ) {
    fprintf( stdout,"\n   -- PHASE 3 : MERGE MESHES OVER PROCESSORS\n" );
  }

  iresult = PMMG_merge_parmesh( parmesh );
  if ( !iresult ) {
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_STRONGFAILURE);
  }

  chrono(OFF,&(ctim[4]));
  if ( parmesh->info.imprim > 0 ) {
    printim(ctim[4].gdif,stim);
    fprintf( stdout,"   -- PHASE 3 COMPLETED.     %s\n",stim );
  }

  if ( !parmesh->myrank ) {
    /** Boundaries reconstruction */
    chrono(ON,&(ctim[5]));
    if (  parmesh->info.imprim > 0 ) {
      fprintf( stdout,"\n   -- PHASE 4 : MESH PACKED UP\n" );
    }

    tmpmem = parmesh->memMax - parmesh->memCur;
    parmesh->memMax = parmesh->memCur;
    parmesh->listgrp[0].mesh->memMax += tmpmem;

    mesh = parmesh->listgrp[0].mesh;
    if ( (!MMG3D_hashTetra( mesh, 0 )) || (-1 == MMG3D_bdryBuild( mesh )) ) {
      /** Impossible to rebuild the triangle */
      fprintf(stdout,"\n\n\n  -- IMPOSSIBLE TO BUILD THE BOUNDARY MESH\n\n\n");
      PMMG_CLEAN_AND_RETURN(parmesh,PMMG_LOWFAILURE);
    }

    chrono(OFF,&(ctim[5]));
    if (  parmesh->info.imprim > 0 ) {
      printim(ctim[5].gdif,stim);
      fprintf( stdout,"   -- PHASE 4 COMPLETED.     %s\n",stim );
    }
  }

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( parmesh->info.imprim >= 0 ) {
    fprintf(stdout,"\n   PARMMGLIB_CENTRALIZED: ELAPSED TIME  %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE PARMMGLIB_CENTRALIZED: IMB-LJLL \n  %s\n",
            PMMG_STR,PMMG_STR);
  }

  PMMG_CLEAN_AND_RETURN(parmesh,ierlib);
}

int PMMG_parmmglib_distributed(PMMG_pParMesh parmesh) {
  MMG5_pMesh       mesh;
  MMG5_pSol        met;
  int              k,ier,iresult,ierlib;
  long int         tmpmem;
  mytime           ctim[TIMEMAX];
  char             stim[32];


  if ( parmesh->info.imprim >= 0 ) {
    fprintf(stdout,"\n  %s\n   MODULE PARMMGLIB_DISTRIBUTED: IMB-LJLL : "
            "%s (%s)\n  %s\n",PMMG_STR,PMMG_VER,PMMG_REL,PMMG_STR);
  }

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /** Check input data */
  chrono(ON,&(ctim[1]));

  ier = PMMG_check_inputData( parmesh );
  MPI_CHECK( MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm ),
             return PMMG_LOWFAILURE);
  if ( !iresult ) return PMMG_LOWFAILURE;

  chrono(OFF,&(ctim[1]));
  printim(ctim[1].gdif,stim);
  if ( parmesh->info.imprim > 0 ) {
    fprintf(stdout,"  -- CHECK INPUT DATA COMPLETED.     %s\n",stim);
  }


  chrono(ON,&(ctim[2]));
  if ( parmesh->info.imprim > 0 ) {
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  }

  assert ( parmesh->ngrp < 2 );
  if ( parmesh->ngrp ) {
    /** Mesh preprocessing: set function pointers, scale mesh, perform mesh
     * analysis and display length and quality histos. */
    ier  = PMMG_preprocessMesh( parmesh );
    mesh = parmesh->listgrp[k].mesh;
    met  = parmesh->listgrp[k].met;
    if ( (ier==PMMG_STRONGFAILURE) && MMG5_unscaleMesh( mesh, met ) ) {
      ier = PMMG_LOWFAILURE;
    }
  }
  else { ier = PMMG_SUCCESS; }

  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MAX, parmesh->comm );
  if ( iresult!=PMMG_SUCCESS ) {
    return iresult;
  }

  chrono(OFF,&(ctim[2]));
  if ( parmesh->info.imprim > 0 ) {
    printim(ctim[1].gdif,stim);
    fprintf(stdout,"   -- PHASE 1 COMPLETED.     %s\n",stim);
  }

  /** Remeshing */
  chrono(ON,&(ctim[3]));
  if ( parmesh->info.imprim > 0 ) {
    fprintf( stdout,"\n  -- PHASE 2 : %s MESHING\n",
             met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC" );
  }

  ier = PMMG_parmmglib1(parmesh);
  MPI_Allreduce( &ier, &ierlib, 1, MPI_INT, MPI_MAX, parmesh->comm );

  chrono(OFF,&(ctim[3]));
  printim(ctim[3].gdif,stim);
  if ( parmesh->info.imprim > 0 ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
  }
  if ( ierlib == PMMG_STRONGFAILURE ) {
    return ierlib;
  }

  /** Boundaries reconstruction */
  chrono(ON,&(ctim[4]));
  if ( parmesh->info.imprim > 0 )
    fprintf(stdout,"\n   -- PHASE 3 : MESH PACKED UP\n");

  /** All the memory is devoted to the mesh **/
  PMMG_parmesh_Free_Comm(parmesh);
  tmpmem = parmesh->memMax - parmesh->memCur;
  parmesh->memMax = parmesh->memCur;
  parmesh->listgrp[0].mesh->memMax += tmpmem;

  mesh = parmesh->listgrp[0].mesh;
  if ( (!MMG3D_hashTetra( mesh, 0 )) || ( -1 == MMG3D_bdryBuild(parmesh->listgrp[0].mesh) ) ) {
    /** Impossible to rebuild the triangle **/
    fprintf(stdout,"\n\n\n  -- IMPOSSIBLE TO BUILD THE BOUNDARY MESH\n\n\n");
    return PMMG_LOWFAILURE;
  }

  chrono(OFF,&(ctim[4]));
  if ( parmesh->info.imprim > 0 ) {
    printim(ctim[4].gdif,stim);
    fprintf(stdout,"\n   -- PHASE 3 COMPLETED.     %s\n",stim);
  }

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( parmesh->info.imprim >= 0 ) {
    fprintf(stdout,"\n   PARMMGLIB_DISTRIBUTED: ELAPSED TIME  %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE PARMMGLIB_DISTRIBUTED: IMB-LJLL \n  %s\n",
            PMMG_STR,PMMG_STR);
  }

  return ierlib;
}

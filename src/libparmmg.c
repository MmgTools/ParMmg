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

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
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
 * \param parmesh pointer toward the parmesh structure
 * \param mesh pointer toward the mesh structure
 *
 * \remark Modeled after the MMG3D_analys function, with additional node/face
 *         communicators indices construction. Also build face comms from node
 *         ones (before tags are passed to edges and nodes).
 */
int PMMG_analys_buildComm(PMMG_pParMesh parmesh,MMG5_pMesh mesh) {
  PMMG_pExt_comm ext_node_comm;
  MMG5_Hash      hash;
  int            i,iext_comm,iext;

  /**--- stage 1: data structures for surface */
  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"\n  ** SURFACE ANALYSIS\n");

  /* create tetra adjacency */
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Hashing problem (1). Exit program.\n");
    return 0;
  }

  /* create prism adjacency */
  if ( !MMG3D_hashPrism(mesh) ) {
    fprintf(stderr,"\n  ## Prism hashing problem. Exit program.\n");
    return 0;
  }
  /* compatibility triangle orientation w/r tetras */
  if ( !MMG5_bdryPerm(mesh) ) {
    fprintf(stderr,"\n  ## Boundary orientation problem. Exit program.\n");
    return 0;
  }

  /* identify surface mesh */
  if ( !MMG5_chkBdryTria(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    return 0;
  }
  MMG5_freeXTets(mesh);
  MMG5_freeXPrisms(mesh);

  /* Set surface triangles to required in nosurf mode or for parallel boundaries */
  MMG3D_set_reqBoundaries(mesh);


  /* create surface adjacency */
  if ( !MMG3D_hashTria(mesh,&hash) ) {
    MMG5_DEL_MEM(mesh,hash.item);
    fprintf(stderr,"\n  ## Hashing problem (2). Exit program.\n");
    return 0;
  }

  /* build hash table for geometric edges */
  if ( !MMG5_hGeom(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem (0). Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    MMG5_DEL_MEM(mesh,mesh->htab.geom);
    return 0;
  }

  /**--- stage 2: surface analysis */
  if ( abs(mesh->info.imprim) > 5  || mesh->info.ddebug )
    fprintf(stdout,"  ** SETTING TOPOLOGY\n");

  /* identify connexity */
  if ( !MMG5_setadj(mesh) ) {
    fprintf(stderr,"\n  ## Topology problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /* check for ridges */
  if ( mesh->info.dhd > MMG5_ANGLIM && !MMG5_setdhd(mesh) ) {
    fprintf(stderr,"\n  ## Geometry problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /* identify singularities */
  if ( !MMG5_singul(mesh) ) {
    fprintf(stderr,"\n  ## MMG5_Singularity problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug )
    fprintf(stdout,"  ** DEFINING GEOMETRY\n");

  /* define (and regularize) normals */
  if ( !MMG5_norver(mesh) ) {
    fprintf(stderr,"\n  ## Normal problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /* set bdry entities to tetra */
  if ( !MMG5_bdrySet(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    return 0;
  }

  /* For both API modes, build communicators indices and set xtetra as PARBDY */
  switch( parmesh->info.API_mode ) {
    case PMMG_APIDISTRIB_faces :
      /* Set face communicators indexing */
      if( !PMMG_build_faceCommIndex( parmesh ) ) return 0;
      /* Convert tria index into iel face index (it needs a valid cc field in
       * each tria), and tag xtetra face as PARBDY before the tag is transmitted
       * to edges and nodes */
      PMMG_tria2elmFace_coords( parmesh );
      break;
    case PMMG_APIDISTRIB_nodes :
      /* Set node communicators indexing */
      if( !PMMG_build_nodeCommIndex( parmesh ) ) return 0;
      /* Build face comms from node ones and set xtetra tags */
      if ( !PMMG_build_faceCommFromNodes(parmesh) ) return PMMG_STRONGFAILURE;
      break;
  }

  /* set non-manifold edges sharing non-intersecting multidomains as required */
  if ( abs(mesh->info.imprim) > 5  || mesh->info.ddebug )
    fprintf(stdout,"  ** UPDATING TOPOLOGY AT NON-MANIFOLD POINTS\n");

  if ( !MMG5_setNmTag(mesh,&hash) ) {
    fprintf(stderr,"\n  ## Non-manifold topology problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    return 0;
  }

  /* check subdomains connected by a vertex and mark these vertex as corner and required */
  MMG5_chkVertexConnectedDomains(mesh);

  /* build hash table for geometric edges */
  if ( !mesh->na && !MMG5_hGeom(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem (0). Exit program.\n");
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    MMG5_DEL_MEM(mesh,mesh->htab.geom);
    return 0;
  }

  /* Update edges tags and references for xtetras */
  if ( !MMG5_bdryUpdate(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    return 0;
  }

  /* define geometry for non manifold points */
  if ( !MMG3D_nmgeom(mesh) ) return 0;

  /* release memory */
  MMG5_DEL_MEM(mesh,mesh->htab.geom);
  MMG5_DEL_MEM(mesh,mesh->adjt);
  MMG5_DEL_MEM(mesh,mesh->tria);
  mesh->nt = 0;

  if ( mesh->nprism ) MMG5_DEL_MEM(mesh,mesh->adjapr);

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
int PMMG_preprocessMesh( PMMG_pParMesh parmesh )
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

  if ( !PMMG_qualhisto(parmesh,PMMG_INQUA) ) {
    return PMMG_STRONGFAILURE;
  }

  /** Mesh analysis */
  if ( !MMG3D_analys(mesh) ) {
    return PMMG_STRONGFAILURE;
  }

  if ( parmesh->info.imprim0 > PMMG_VERB_ITWAVES && (!mesh->info.iso) && met->m ) {
    PMMG_prilen(parmesh,0);
  }

  /** Mesh unscaling */
  if ( !MMG5_unscaleMesh(mesh,met) ) {
    return PMMG_STRONGFAILURE;
  }

  return PMMG_SUCCESS;
}

/**
 * \param  parmesh pointer to parmesh structure
 *
 * \return PMMG_SUCCESS if success, PMMG_LOWFAILURE if fail and return an
 * unscaled mesh, PMMG_STRONGFAILURE if fail and return a scaled mesh.
 *
 * Mesh preprocessing (distributed input mesh): set function pointers,
 * scale mesh, perform mesh analysis, display length and quality histos,
 * and build communicators.
 */
int PMMG_preprocessMesh_distributed( PMMG_pParMesh parmesh )
{
  MMG5_pMesh mesh;
  MMG5_pSol  met;

  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  assert ( ( mesh != NULL ) && ( met != NULL ) && "Preprocessing empty args");

  /** Set distributed API mode. Interface faces OR nodes need to be set by the
   * user through the API interface at this point, meening that the
   * corresponding external comm is set to the correct size, and filled with
   * local entity indices (for node comms, also itosend and itorecv arrays are
   * filled with local/global node IDs).
   * .
   * - By default, if interface faces are provided by the user, the input nodes
   *   are ignored and node communicators are built from the face ones.
   * - Otherwise, interfaces nodes serve to build node communicators, then face
   *   comms are built from node ones.
   */
  if( parmesh->next_face_comm ) {
    parmesh->info.API_mode = PMMG_APIDISTRIB_faces;
  } else if( parmesh->next_node_comm ) {
    parmesh->info.API_mode = PMMG_APIDISTRIB_nodes;
  } else {
    fprintf(stderr," ## Error: %s: parallel faces or nodes must be set through the API interface\n",__func__);
    return 0;
  }

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

  if ( !PMMG_qualhisto(parmesh,PMMG_INQUA) ) {
    return PMMG_STRONGFAILURE;
  }

  /** Mesh analysis, face/node communicators indices construction (depending
   * from the API mode), build face comms from node ones */
  if ( !PMMG_analys_buildComm(parmesh,mesh) ) {
    return PMMG_STRONGFAILURE;
  }

  if ( parmesh->info.imprim > PMMG_VERB_ITWAVES && (!mesh->info.iso) && met->m ) {
    MMG3D_prilen(mesh,met,0);
  }

  /** Mesh unscaling */
  if ( !MMG5_unscaleMesh(mesh,met) ) {
    return PMMG_STRONGFAILURE;
  }

  /** Build node communicators from face ones (here because the (mesh needs to
   * be unscaled) */
  if( parmesh->info.API_mode == PMMG_APIDISTRIB_faces ) {
    PMMG_parmesh_ext_comm_free( parmesh,parmesh->ext_node_comm,parmesh->next_node_comm);
    PMMG_DEL_MEM(parmesh, parmesh->ext_node_comm,PMMG_Ext_comm,"ext node comm");
    parmesh->next_node_comm = 0;
    PMMG_DEL_MEM(parmesh, parmesh->ext_node_comm,PMMG_Int_comm,"int node comm");
    PMMG_CALLOC(parmesh,parmesh->int_node_comm,1,PMMG_Int_comm,"int node comm",return 0);
    if ( !PMMG_build_nodeCommFromFaces(parmesh) ) {
      return PMMG_STRONGFAILURE;
    }
  }

  assert ( PMMG_check_extFaceComm ( parmesh ) );
  assert ( PMMG_check_intFaceComm ( parmesh ) );
  assert ( PMMG_check_extNodeComm ( parmesh ) );
  assert ( PMMG_check_intNodeComm ( parmesh ) );

  return PMMG_SUCCESS;
}

int PMMG_parmmglib_centralized(PMMG_pParMesh parmesh) {
  MMG5_pMesh    mesh;
  MMG5_pSol     met;
  int           ier;
  int           iresult,ierlib;
  long int      tmpmem;
  mytime        ctim[TIMEMAX];
  int8_t        tim;
  char          stim[32];

 if ( parmesh->info.imprim > PMMG_VERB_NO ) {
    fprintf(stdout,"\n  %s\n   MODULE PARMMGLIB_CENTRALIZED: IMB-LJLL : "
            "%s (%s)\n  %s\n",PMMG_STR,PMMG_VER,PMMG_REL,PMMG_STR);
  }

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /** Check input data */
  tim = 1;
  chrono(ON,&(ctim[tim]));

  ier = PMMG_check_inputData( parmesh );
  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !iresult ) return PMMG_LOWFAILURE;

  chrono(OFF,&(ctim[tim]));
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf(stdout,"  -- CHECK INPUT DATA COMPLETED.     %s\n",stim);
  }


  chrono(ON,&(ctim[2]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS AND MESH DISTRIBUTION\n");
  }

  /** Send mesh to other procs */
  tim = 6;
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    chrono(ON,&(ctim[tim]));
    fprintf(stdout,"\n  -- BCAST" );
  }
  ier = PMMG_bcast_mesh( parmesh );
  if ( ier!=1 ) return PMMG_LOWFAILURE;

  if ( parmesh->info.imprim >= PMMG_VERB_STEPS  ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"\n  -- BCAST COMPLETED    %s\n",stim );
  }

  /** Mesh preprocessing: set function pointers, scale mesh, perform mesh
   * analysis and display length and quality histos. */
  tim = 7;
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    chrono(ON,&(ctim[tim]));
    fprintf(stdout,"\n  -- ANALYSIS" );
  }
  ier = PMMG_preprocessMesh( parmesh );
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
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
  tim = 8;
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    chrono(ON,&(ctim[tim]));
    fprintf(stdout,"\n  -- PARTITIONING" );
  }
  if ( !PMMG_distribute_mesh( parmesh ) ) {
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_LOWFAILURE);
  }
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"\n  -- PARTITIONING COMPLETED    %s\n",stim );
  }

  chrono(OFF,&(ctim[2]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    printim(ctim[2].gdif,stim);
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);
  }

  /** Remeshing */
  tim = 3;
  chrono(ON,&(ctim[tim]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf( stdout,"\n  -- PHASE 2 : %s MESHING\n",
             met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC" );
  }

  ier = PMMG_parmmglib1(parmesh);
  MPI_Allreduce( &ier, &ierlib, 1, MPI_INT, MPI_MAX, parmesh->comm );

  chrono(OFF,&(ctim[tim]));
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
  }
  if ( ierlib == PMMG_STRONGFAILURE ) {
    return ierlib;
  }

  /** Merge all the meshes on the proc 0 */
  tim = 4;
  chrono(ON,&(ctim[tim]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf( stdout,"\n   -- PHASE 3 : MERGE MESHES OVER PROCESSORS\n" );
  }

  iresult = PMMG_merge_parmesh( parmesh );
  if ( !iresult ) {
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_STRONGFAILURE);
  }

  chrono(OFF,&(ctim[tim]));
  if ( parmesh->info.imprim >  PMMG_VERB_VERSION  ) {
    printim(ctim[tim].gdif,stim);
    fprintf( stdout,"   -- PHASE 3 COMPLETED.     %s\n",stim );
  }

  if ( !parmesh->myrank ) {
    /** Boundaries reconstruction */
    tim = 5;
    chrono(ON,&(ctim[tim]));
    if (  parmesh->info.imprim > PMMG_VERB_VERSION ) {
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

    chrono(OFF,&(ctim[tim]));
    if (  parmesh->info.imprim >  PMMG_VERB_VERSION ) {
      printim(ctim[tim].gdif,stim);
      fprintf( stdout,"   -- PHASE 4 COMPLETED.     %s\n",stim );
    }
  }

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( parmesh->info.imprim >= PMMG_VERB_VERSION ) {
    fprintf(stdout,"\n   PARMMGLIB_CENTRALIZED: ELAPSED TIME  %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE PARMMGLIB_CENTRALIZED: IMB-LJLL \n  %s\n",
            PMMG_STR,PMMG_STR);
  }

  PMMG_CLEAN_AND_RETURN(parmesh,ierlib);
}

int PMMG_parmmglib_distributed(PMMG_pParMesh parmesh) {
  MMG5_pMesh       mesh;
  MMG5_pSol        met;
  int              ier,iresult,ierlib;
  long int         tmpmem;
  mytime           ctim[TIMEMAX];
  int8_t           tim;
  char             stim[32];


  if ( parmesh->info.imprim >= PMMG_VERB_VERSION ) {
    fprintf(stdout,"\n  %s\n   MODULE PARMMGLIB_DISTRIBUTED: IMB-LJLL : "
            "%s (%s)\n  %s\n",PMMG_STR,PMMG_VER,PMMG_REL,PMMG_STR);
  }

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /** Check input data */
  tim = 1;
  chrono(ON,&(ctim[tim]));

  ier = PMMG_check_inputData( parmesh );
  MPI_CHECK( MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm ),
             return PMMG_LOWFAILURE);
  if ( !iresult ) return PMMG_LOWFAILURE;

  chrono(OFF,&(ctim[tim]));
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf(stdout,"  -- CHECK INPUT DATA COMPLETED.     %s\n",stim);
  }

  tim = 2;
  chrono(ON,&(ctim[tim]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  }

  assert ( parmesh->ngrp < 2 );
  if ( parmesh->ngrp ) {
    /** Mesh preprocessing: set function pointers, scale mesh, perform mesh
     * analysis and display length and quality histos. */
    ier  = PMMG_preprocessMesh_distributed( parmesh );
    mesh = parmesh->listgrp[0].mesh;
    met  = parmesh->listgrp[0].met;
    if ( (ier==PMMG_STRONGFAILURE) && MMG5_unscaleMesh( mesh, met ) ) {
      ier = PMMG_LOWFAILURE;
    }
  }
  else { ier = PMMG_SUCCESS; }

  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MAX, parmesh->comm );
  if ( iresult!=PMMG_SUCCESS ) {
    return iresult;
  }

  chrono(OFF,&(ctim[tim]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"   -- PHASE 1 COMPLETED.     %s\n",stim);
  }

  /** Remeshing */
  tim = 3;
  chrono(ON,&(ctim[tim]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf( stdout,"\n  -- PHASE 2 : %s MESHING\n",
             met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC" );
  }

  ier = PMMG_parmmglib1(parmesh);
  MPI_Allreduce( &ier, &ierlib, 1, MPI_INT, MPI_MAX, parmesh->comm );

  chrono(OFF,&(ctim[tim]));
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
  }
  if ( ierlib == PMMG_STRONGFAILURE ) {
    return ierlib;
  }

  /** Boundaries reconstruction */
  tim = 4;
  chrono(ON,&(ctim[tim]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n   -- PHASE 3 : MESH PACKED UP\n");

  /** All the memory is devoted to the mesh **/
  tmpmem = parmesh->memMax - parmesh->memCur;
  parmesh->memMax = parmesh->memCur;
  parmesh->listgrp[0].mesh->memMax += tmpmem;

  mesh = parmesh->listgrp[0].mesh;
  if ( (!MMG3D_hashTetra( mesh, 0 )) || ( -1 == MMG3D_bdryBuild(parmesh->listgrp[0].mesh) ) ) {
    /** Impossible to rebuild the triangle **/
    fprintf(stdout,"\n\n\n  -- IMPOSSIBLE TO BUILD THE BOUNDARY MESH\n\n\n");
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_LOWFAILURE);
  }

  chrono(OFF,&(ctim[tim]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"\n   -- PHASE 3 COMPLETED.     %s\n",stim);
  }

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( parmesh->info.imprim >= PMMG_VERB_VERSION ) {
    fprintf(stdout,"\n   PARMMGLIB_DISTRIBUTED: ELAPSED TIME  %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE PARMMGLIB_DISTRIBUTED: IMB-LJLL \n  %s\n",
            PMMG_STR,PMMG_STR);
  }

  PMMG_CLEAN_AND_RETURN(parmesh,ierlib);
}

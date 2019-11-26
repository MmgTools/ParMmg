/* =============================================================================
**  This file is part of the parmmg software package for parallel tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux, 2017-
**
**  parmmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  parmmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with parmmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the parmmg distribution only if you accept them.
** =============================================================================
*/

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
#include "git_log_pmmg.h"


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
  MMG5_Hash      hash;

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
      PMMG_parmesh_ext_comm_free( parmesh,parmesh->ext_face_comm,parmesh->next_face_comm);
      PMMG_DEL_MEM(parmesh, parmesh->ext_face_comm,PMMG_Ext_comm,"ext face comm");
      parmesh->next_face_comm = 0;
      PMMG_DEL_MEM(parmesh, parmesh->int_face_comm,PMMG_Int_comm,"int face comm");
      if ( !PMMG_build_faceCommFromNodes(parmesh) ) return PMMG_STRONGFAILURE;
      break;
  }

  /* Tag parallel faces on material interfaces as boundary */
  if( !PMMG_parbdySet( parmesh ) ) {
    fprintf(stderr,"\n  ## Unable to recognize parallel faces on material interfaces. Exit program.\n");
    return 0;
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
  if ( !MMG5_scaleMesh(mesh,met,NULL) ) {
    return PMMG_LOWFAILURE;
  }

  /** specific meshing */
  if ( mesh->info.optim && !met->np ) {
    if ( !MMG3D_doSol(mesh,met) ) {
      return PMMG_STRONGFAILURE;
    }
    MMG5_solTruncatureForOptim(mesh,met);
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

  if ( !PMMG_qualhisto(parmesh,PMMG_INQUA,1) ) {
    return PMMG_STRONGFAILURE;
  }

  /** Mesh analysis */
  if ( !MMG3D_analys(mesh) ) {
    return PMMG_STRONGFAILURE;
  }

//  if ( parmesh->info.imprim0 > PMMG_VERB_ITWAVES && (!mesh->info.iso) && met->m ) {
//    PMMG_prilen(parmesh,0);
//  }

  /** Mesh unscaling */
  if ( !MMG5_unscaleMesh(mesh,met,NULL) ) {
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

  /** Check distributed API mode. Interface faces OR nodes need to be set by the
   * user through the API interface at this point, meening that the
   * corresponding external comm is set to the correct size, and filled with
   * local entity indices (for node comms, also itosend and itorecv arrays are
   * filled with local/global node IDs).
  */
  if( parmesh->info.API_mode == PMMG_APIDISTRIB_faces && !parmesh->next_face_comm ) {
    fprintf(stderr," ## Error: %s: parallel interface faces must be set through the API interface\n",__func__);
    return 0;
  } else if( parmesh->info.API_mode == PMMG_APIDISTRIB_nodes && !parmesh->next_node_comm ) {
    fprintf(stderr," ## Error: %s: parallel interface nodes must be set through the API interface\n",__func__);
    return 0;
  }

  /** Function setters (must be assigned before quality computation) */
  MMG3D_Set_commonFunc();

  /** Mesh scaling and quality histogram */
  if ( !MMG5_scaleMesh(mesh,met,NULL) ) {
    return PMMG_LOWFAILURE;
  }

  /** specific meshing */
  if ( mesh->info.optim && !met->np ) {
    if ( !MMG3D_doSol(mesh,met) ) {
      return PMMG_STRONGFAILURE;
    }
    MMG5_solTruncatureForOptim(mesh,met);
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

  /** Mesh analysis, face/node communicators indices construction (depending
   * from the API mode), build face comms from node ones */
  if ( !PMMG_analys_buildComm(parmesh,mesh) ) {
    return PMMG_STRONGFAILURE;
  }

  if ( parmesh->info.imprim > PMMG_VERB_ITWAVES && (!mesh->info.iso) && met->m ) {
    MMG3D_prilen(mesh,met,0);
  }

  /** Mesh unscaling */
  if ( !MMG5_unscaleMesh(mesh,met,NULL) ) {
    return PMMG_STRONGFAILURE;
  }

  /** Build node communicators from face ones (here because the (mesh needs to
   * be unscaled) */
  if( parmesh->info.API_mode == PMMG_APIDISTRIB_faces ) {
    PMMG_parmesh_ext_comm_free( parmesh,parmesh->ext_node_comm,parmesh->next_node_comm);
    PMMG_DEL_MEM(parmesh, parmesh->ext_node_comm,PMMG_Ext_comm,"ext node comm");
    parmesh->next_node_comm = 0;
    PMMG_DEL_MEM(parmesh, parmesh->int_node_comm,PMMG_Int_comm,"int node comm");
    PMMG_CALLOC(parmesh,parmesh->int_node_comm,1,PMMG_Int_comm,"int node comm",return 0);
    if ( !PMMG_build_nodeCommFromFaces(parmesh) ) {
      return PMMG_STRONGFAILURE;
    }
  }

  if ( !PMMG_qualhisto(parmesh,PMMG_INQUA,0) ) {
    return PMMG_STRONGFAILURE;
  }

  assert ( PMMG_check_extFaceComm ( parmesh ) );
  assert ( PMMG_check_intFaceComm ( parmesh ) );
  assert ( PMMG_check_extNodeComm ( parmesh ) );
  assert ( PMMG_check_intNodeComm ( parmesh ) );

  return PMMG_SUCCESS;
}

int PMMG_distributeMesh_centralized_timers( PMMG_pParMesh parmesh,mytime *ctim ) {
  MMG5_pMesh    mesh;
  MMG5_pSol     met;
  int           ier,iresult;
  int8_t        tim;
  char          stim[32];

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

//  /** Send mesh to other procs */
//  tim = 6;
//  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
//    chrono(ON,&(ctim[tim]));
//    fprintf(stdout,"\n  -- BCAST" );
//  }
//  ier = PMMG_bcast_mesh( parmesh );
//  if ( ier!=1 ) return PMMG_LOWFAILURE;
//
//  if ( parmesh->info.imprim >= PMMG_VERB_STEPS  ) {
//    chrono(OFF,&(ctim[tim]));
//    printim(ctim[tim].gdif,stim);
//    fprintf(stdout,"\n  -- BCAST COMPLETED    %s\n",stim );
//  }

  /** Mesh preprocessing: set function pointers, scale mesh, perform mesh
   * analysis and display length and quality histos. */
  if( parmesh->myrank == parmesh->info.root ) {
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
    if ( (ier==PMMG_STRONGFAILURE) && MMG5_unscaleMesh( mesh, met, NULL ) ) {
      ier = PMMG_LOWFAILURE;
    }

    /* Memory repartition */
    if ( !PMMG_parmesh_updateMemMax( parmesh,50,1 ) ) ier = 3;

  }
//  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MAX, parmesh->comm );
//  if ( iresult!=PMMG_SUCCESS ) {
//    return iresult;
//  }

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

  /** Function setters (must be assigned before quality computation) */
  if( parmesh->myrank != parmesh->info.root ) {
    mesh = parmesh->listgrp[0].mesh;
    met  = parmesh->listgrp[0].met;
    MMG3D_Set_commonFunc();
    MMG3D_setfunc(mesh,met);
  }

  chrono(OFF,&(ctim[2]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    printim(ctim[2].gdif,stim);
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);
  }

  iresult = PMMG_SUCCESS;
  return iresult;
}

int PMMG_parmmglib_centralized(PMMG_pParMesh parmesh) {
  PMMG_pGrp     grp;
  MMG5_pMesh    mesh;
  MMG5_pSol     met;
  int           ier;
  int           iresult,ierlib,npmax,xpmax,nemax,xtmax;
  long int      tmpmem;
  mytime        ctim[TIMEMAX];
  int8_t        tim;
  char          stim[32];

 if ( parmesh->info.imprim > PMMG_VERB_NO ) {
    fprintf(stdout,"\n  %s\n   MODULE PARMMGLIB_CENTRALIZED: IMB-LJLL : "
            "%s (%s)\n  %s\n",PMMG_STR,PMMG_VER,PMMG_REL,PMMG_STR);
    fprintf(stdout,"     git branch: %s\n",PMMG_GIT_BRANCH);
    fprintf(stdout,"     git commit: %s\n",PMMG_GIT_COMMIT);
    fprintf(stdout,"     git date:   %s\n\n",PMMG_GIT_DATE);
  }

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /* Distribute the mesh */
  ier = PMMG_distributeMesh_centralized_timers( parmesh, ctim );
  if( ier != PMMG_SUCCESS ) return ier;

  grp    = &parmesh->listgrp[0];
  mesh   = grp->mesh;
  met    = grp->met;

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

#warning remove the lib_centralized and lib_distributed library to have modular centralized input + annalysis or parallel input + analysis , libparmmg1 call, then centralized or distributed output
  switch ( parmesh->info.fmtout ) {
  case ( MMG5_FMT_VtkPvtu ):
    // Distributed Output
#warning boundaries arent rebuilded

    break;
  default:
    // Centralized Output
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
      mesh  = parmesh->listgrp[0].mesh;
      npmax = mesh->npmax;
      nemax = mesh->nemax;
      xpmax = mesh->xpmax;
      xtmax = mesh->xtmax;
      mesh->npmax = mesh->np;
      mesh->nemax = mesh->ne;
      mesh->xpmax = mesh->xp;
      mesh->xtmax = mesh->xt;

      if ( !PMMG_setMemMax_realloc( mesh, npmax, xpmax, nemax, xtmax ) ) {
        fprintf(stdout,"\n\n\n  -- LACK OF MEMORY\n\n\n");
        PMMG_CLEAN_AND_RETURN(parmesh,PMMG_LOWFAILURE);
      }

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
  int              ier,iresult,ierlib,nepax,npmax,nemax,xpmax,xtmax;
  long int         tmpmem;
  mytime           ctim[TIMEMAX];
  int8_t           tim;
  char             stim[32];


  if ( parmesh->info.imprim >= PMMG_VERB_VERSION ) {
    fprintf(stdout,"\n  %s\n   MODULE PARMMGLIB_DISTRIBUTED: IMB-LJLL : "
            "%s (%s)\n  %s\n",PMMG_STR,PMMG_VER,PMMG_REL,PMMG_STR);
    fprintf(stdout,"     git branch: %s\n",PMMG_GIT_BRANCH);
    fprintf(stdout,"     git commit: %s\n",PMMG_GIT_COMMIT);
    fprintf(stdout,"     git date:   %s\n\n",PMMG_GIT_DATE);
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
    if ( (ier==PMMG_STRONGFAILURE) && MMG5_unscaleMesh( mesh, met, NULL ) ) {
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

  mesh  = parmesh->listgrp[0].mesh;
  npmax = mesh->npmax;
  nemax = mesh->nemax;
  xpmax = mesh->xpmax;
  xtmax = mesh->xtmax;
  mesh->npmax = mesh->np;
  mesh->nemax = mesh->ne;
  mesh->xpmax = mesh->xp;
  mesh->xtmax = mesh->xt;
  if ( !PMMG_setMemMax_realloc( mesh, npmax, xpmax, nemax, xtmax ) ) {
    fprintf(stdout,"\n\n\n  -- LACK OF MEMORY\n\n\n");
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_LOWFAILURE);
  }

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

int PMMG_distributeMesh_centralized( PMMG_pParMesh parmesh ) {
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  int ier,iresult;

  ier = PMMG_check_inputData( parmesh );
  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !iresult ) return PMMG_LOWFAILURE;

  /** Send mesh to other procs */
  ier = PMMG_bcast_mesh( parmesh );
  if ( ier!=1 ) return PMMG_LOWFAILURE;

  /** Mesh preprocessing: set function pointers, scale mesh, perform mesh
   * analysis and display length and quality histos. */
  ier = PMMG_preprocessMesh( parmesh );

  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;
  if ( (ier==PMMG_STRONGFAILURE) && MMG5_unscaleMesh( mesh, met, NULL ) ) {
    ier = PMMG_LOWFAILURE;
  }
  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MAX, parmesh->comm );
  if ( iresult!=PMMG_SUCCESS ) {
    return iresult;
  }

  /** Send mesh partionning to other procs */
  if ( !PMMG_partBcast_mesh( parmesh ) ) {
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_LOWFAILURE);
  }

  return PMMG_SUCCESS;
}

/**
 * \param parmesh pointer toward parmesh structure
 * \param color_out array of interface colors
 * \param ifc_node_loc local IDs of interface nodes
 * \param ifc_node_glob global IDs of interface nodes
 * \param next_node_comm number of node interfaces
 * \param nitem_node_comm number of nodes on each interface
 *
 * Create global IDs for nodes on parallel interfaces.
 *
 */
int PMMG_color_intfcNode(PMMG_pParMesh parmesh,int *color_out,
                         int **ifc_node_loc,int **ifc_node_glob,
                         int next_node_comm,int *nitem_node_comm) {
  MMG5_pMesh     mesh;
  MMG5_pPoint    ppt;
  MPI_Request    request;
  MPI_Status     status;
  int            npairs_loc,*npairs,*displ_pair,*glob_pair_displ;
  int            src,dst,tag,sendbuffer,recvbuffer,iproc,icomm,i,idx;

  mesh = parmesh->listgrp[0].mesh;

  PMMG_CALLOC(parmesh,npairs,parmesh->nprocs,int,"npair",return 0);
  PMMG_CALLOC(parmesh,displ_pair,parmesh->nprocs+1,int,"displ_pair",return 0);

  /* Use points flag to mark interface points:
   * - interface points are initialised as PMMG_UNSET;
   * - once a point is given a global ID, it is flagged with it so that other
   *   interfaces can see it.
   */
  for( icomm = 0; icomm < next_node_comm; icomm++ )
    for( i=0; i < nitem_node_comm[icomm]; i++ ) {
      ppt = &mesh->point[ifc_node_loc[icomm][i]];
      ppt->flag = PMMG_UNSET;
    }

  /* Count nb of new pair nodes hosted on proc */
  npairs_loc = 0;
  for( icomm = 0; icomm < next_node_comm; icomm++ )
    if( color_out[icomm] > parmesh->myrank ) npairs_loc += nitem_node_comm[icomm];//1;

  /* Get nb of pair nodes and compute pair offset */
  MPI_Allgather( &npairs_loc,1,MPI_INT,
                 npairs,1,MPI_INT,parmesh->comm );

  for( iproc = 0; iproc < parmesh->nprocs; iproc++ )
    displ_pair[iproc+1] = displ_pair[iproc]+npairs[iproc];

 
  PMMG_CALLOC(parmesh,glob_pair_displ,next_node_comm+1,int,"glob_pair_displ",return 0); 
  for( icomm = 0; icomm < next_node_comm; icomm++ )
    glob_pair_displ[icomm] = displ_pair[parmesh->myrank];
  for( icomm = 0; icomm < next_node_comm; icomm++ ) {
    if( color_out[icomm] > parmesh->myrank )
      glob_pair_displ[icomm+1] = glob_pair_displ[icomm]+nitem_node_comm[icomm];//+1;
  }

  /* Compute global pair nodes enumeration (injective, non-surjective map) */
  for( icomm = 0; icomm < next_node_comm; icomm++ ) {
    
    /* Assign global index */
    src = fmin(parmesh->myrank,color_out[icomm]);
    dst = fmax(parmesh->myrank,color_out[icomm]);
    tag = parmesh->nprocs*src+dst;
    if( parmesh->myrank == src ) {
      sendbuffer = glob_pair_displ[icomm];
      MPI_CHECK( MPI_Isend(&sendbuffer,1,MPI_INT,dst,tag,
                            parmesh->comm,&request),return 0 );
    }
    if ( parmesh->myrank == dst ) {
      MPI_CHECK( MPI_Recv(&recvbuffer,1,MPI_INT,src,tag,
                          parmesh->comm,&status),return 0 );
      glob_pair_displ[icomm] = recvbuffer;
    }
  }


  /* Each proc buils global IDs if color_in < color_out, then sends IDs to
   * color_out;
   *  */
  for( icomm = 0; icomm < next_node_comm; icomm++ ) {
    src = fmin(parmesh->myrank,color_out[icomm]);
    dst = fmax(parmesh->myrank,color_out[icomm]);
    tag = parmesh->nprocs*src+dst;
    /* Recv IDs from previous proc */
    if ( parmesh->myrank == dst ) {
      MPI_CHECK( MPI_Recv(ifc_node_glob[icomm],nitem_node_comm[icomm],MPI_INT,src,tag,
                          parmesh->comm,&status),return 0 );
      /* Update flag so that you can use it to build your own IDs */
      for( i=0; i < nitem_node_comm[icomm]; i++ ) {
        ppt = &mesh->point[ifc_node_loc[icomm][i]];
        ppt->flag = ifc_node_glob[icomm][i];
      }
    }
    /* Build your own IDs and send them to next proc */
    if( parmesh->myrank == src ) {
      idx = 1; /* index starts from 1 */
      for( i=0; i < nitem_node_comm[icomm]; i++ ) {
        ppt = &mesh->point[ifc_node_loc[icomm][i]];
        if( ppt->flag == PMMG_UNSET ) {
          ppt->flag = glob_pair_displ[icomm]+idx++;
        }
        ifc_node_glob[icomm][i] = ppt->flag;
      }
      MPI_CHECK( MPI_Isend(ifc_node_glob[icomm],nitem_node_comm[icomm],MPI_INT,dst,tag,
                            parmesh->comm,&request),return 0 );
    }
  }

  /* Free arrays */
  PMMG_DEL_MEM(parmesh,npairs,int,"npairs");
  PMMG_DEL_MEM(parmesh,displ_pair,int,"displ_pair");
  PMMG_DEL_MEM(parmesh,glob_pair_displ,int,"glob_pair_displ");


  /* Check global IDs */
  int **itorecv;
  PMMG_CALLOC(parmesh,itorecv,next_node_comm,int*,"itorecv pointer",return 0); 
  for( icomm = 0; icomm < next_node_comm; icomm++ ) {
    PMMG_CALLOC(parmesh,itorecv[icomm],nitem_node_comm[icomm],int,"itorecv",return 0); 
    
    src = fmin(parmesh->myrank,color_out[icomm]);
    dst = fmax(parmesh->myrank,color_out[icomm]);
    tag = parmesh->nprocs*src+dst;
    if( parmesh->myrank == src ) {
      MPI_CHECK( MPI_Isend(ifc_node_glob[icomm],nitem_node_comm[icomm],MPI_INT,dst,tag,
                            parmesh->comm,&request),return 0 );
    }
    if ( parmesh->myrank == dst ) {
      MPI_CHECK( MPI_Recv(itorecv[icomm],nitem_node_comm[icomm],MPI_INT,src,tag,
                          parmesh->comm,&status),return 0 );
      for( i=0; i < nitem_node_comm[icomm]; i++ )
        assert( ifc_node_glob[icomm][i] == itorecv[icomm][i] );
    }
  }

  for( icomm = 0; icomm < next_node_comm; icomm++ )
    PMMG_DEL_MEM(parmesh,itorecv[icomm],int,"itorecv"); 
  PMMG_DEL_MEM(parmesh,itorecv,int*,"itorecv pointer"); 
  
  return 1;
}

/**
 * \param parmesh pointer toward parmesh structure
 * \param color_out array of interface colors
 * \param ifc_tria_loc local IDs of interface triangles
 * \param ifc_tria_glob global IDs of interface triangles
 * \param next_face_comm number of triangle interfaces
 * \param nitem_face_comm number of triangles on each interface
 *
 * Create global IDs for triangles on parallel interfaces.
 *
 */
int PMMG_color_intfcTria(PMMG_pParMesh parmesh,int *color_out,
                         int **ifc_tria_loc,int **ifc_tria_glob,
                         int next_face_comm,int *nitem_face_comm) {
  MPI_Request    request;
  MPI_Status     status;
  int            npairs_loc,*npairs,*displ_pair,*glob_pair_displ;
  int            src,dst,tag,sendbuffer,recvbuffer,iproc,icomm,i;

  PMMG_CALLOC(parmesh,npairs,parmesh->nprocs,int,"npair",return 0);
  PMMG_CALLOC(parmesh,displ_pair,parmesh->nprocs+1,int,"displ_pair",return 0);

  /* Count nb of new pair faces hosted on proc */
  npairs_loc = 0;
  for( icomm = 0; icomm < next_face_comm; icomm++ )
    if( color_out[icomm] > parmesh->myrank ) npairs_loc += nitem_face_comm[icomm];//1;

  /* Get nb of pair faces and compute pair offset */
  MPI_Allgather( &npairs_loc,1,MPI_INT,
                 npairs,1,MPI_INT,parmesh->comm );

  for( iproc = 0; iproc < parmesh->nprocs; iproc++ )
    displ_pair[iproc+1] = displ_pair[iproc]+npairs[iproc];

  
  PMMG_CALLOC(parmesh,glob_pair_displ,next_face_comm+1,int,"glob_pair_displ",return 0); 
  for( icomm = 0; icomm < next_face_comm; icomm++ )
    glob_pair_displ[icomm] = displ_pair[parmesh->myrank];
  for( icomm = 0; icomm < next_face_comm; icomm++ ) {
    if( color_out[icomm] > parmesh->myrank )
      glob_pair_displ[icomm+1] = glob_pair_displ[icomm]+nitem_face_comm[icomm];//+1;
  }

  /* Compute global pair faces enumeration (injective, non-surjective map) */
  for( icomm = 0; icomm < next_face_comm; icomm++ ) {
    
    /* Assign global index */
    src = fmin(parmesh->myrank,color_out[icomm]);
    dst = fmax(parmesh->myrank,color_out[icomm]);
    tag = parmesh->nprocs*src+dst;
    if( parmesh->myrank == src ) {
      sendbuffer = glob_pair_displ[icomm];
      MPI_CHECK( MPI_Isend(&sendbuffer,1,MPI_INT,dst,tag,
                            parmesh->comm,&request),return 0 );
    }
    if ( parmesh->myrank == dst ) {
      MPI_CHECK( MPI_Recv(&recvbuffer,1,MPI_INT,src,tag,
                          parmesh->comm,&status),return 0 );
      glob_pair_displ[icomm] = recvbuffer;
    }
  }

  for( icomm = 0; icomm < next_face_comm; icomm++ )
    for( i=0; i < nitem_face_comm[icomm]; i++ )
      ifc_tria_glob[icomm][i] = glob_pair_displ[icomm]+i+1; /* index starts from 1 */

  /* Free arrays */
  PMMG_DEL_MEM(parmesh,npairs,int,"npairs");
  PMMG_DEL_MEM(parmesh,displ_pair,int,"displ_pair");
  PMMG_DEL_MEM(parmesh,glob_pair_displ,int,"glob_pair_displ");


  /* Check global IDs */
  int **itorecv;
  PMMG_CALLOC(parmesh,itorecv,next_face_comm,int*,"itorecv pointer",return 0); 
  for( icomm = 0; icomm < next_face_comm; icomm++ ) {
    PMMG_CALLOC(parmesh,itorecv[icomm],nitem_face_comm[icomm],int,"itorecv",return 0); 
    
    src = fmin(parmesh->myrank,color_out[icomm]);
    dst = fmax(parmesh->myrank,color_out[icomm]);
    tag = parmesh->nprocs*src+dst;
    if( parmesh->myrank == src ) {
      MPI_CHECK( MPI_Isend(ifc_tria_glob[icomm],nitem_face_comm[icomm],MPI_INT,dst,tag,
                            parmesh->comm,&request),return 0 );
    }
    if ( parmesh->myrank == dst ) {
      MPI_CHECK( MPI_Recv(itorecv[icomm],nitem_face_comm[icomm],MPI_INT,src,tag,
                          parmesh->comm,&status),return 0 );
      for( i=0; i < nitem_face_comm[icomm]; i++ )
        assert( ifc_tria_glob[icomm][i] == itorecv[icomm][i] );
    }
  }

  for( icomm = 0; icomm < next_face_comm; icomm++ )
    PMMG_DEL_MEM(parmesh,itorecv[icomm],int,"itorecv");
  PMMG_DEL_MEM(parmesh,itorecv,int*,"itorecv pointer"); 
 
  return 1;
}

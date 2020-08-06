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

/* Declared in the header, but defined at compile time */
int (*PMMG_interp4bar)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTetra pt,int,PMMG_barycoord *barycoord);
int (*PMMG_interp3bar)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTria ptr,int,PMMG_barycoord *barycoord);
int (*PMMG_interp2bar)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTria ptr,int ip,int l,PMMG_barycoord *barycoord);

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
  /* Don't reset the hmin value computed when unscaling the mesh */
  if ( !parmesh->info.sethmin ) {
    mesh->info.sethmin = 1;
  }
  /* Don't reset the hmax value computed when unscaling the mesh */
  if ( !parmesh->info.sethmax ) {
    mesh->info.sethmax = 1;
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
  PMMG_setfunc(parmesh);

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

  if ( parmesh->info.imprim0 > PMMG_VERB_ITWAVES && (!mesh->info.iso) && met->m ) {
    PMMG_prilen(parmesh,0,1);
  }

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
  if( parmesh->nprocs >1 ) {
    if( parmesh->info.API_mode == PMMG_APIDISTRIB_faces && !parmesh->next_face_comm ) {
      fprintf(stderr," ## Error: %s: parallel interface faces must be set through the API interface\n",__func__);
      return PMMG_STRONGFAILURE;
    } else if( parmesh->info.API_mode == PMMG_APIDISTRIB_nodes && !parmesh->next_node_comm ) {
      fprintf(stderr," ## Error: %s: parallel interface nodes must be set through the API interface\n",__func__);
      return PMMG_STRONGFAILURE;
    }
  }

  /** Function setters (must be assigned before quality computation) */
  MMG3D_Set_commonFunc();

  /** Mesh scaling and quality histogram */
  if ( !MMG5_scaleMesh(mesh,met,NULL) ) {
    return PMMG_LOWFAILURE;
  }
#warning hmin/hmax computed on each proc while we want a global value from the global bounding box and/or the global metric field...
  /* Don't reset the hmin value computed when unscaling the mesh */
  if ( !parmesh->info.sethmin ) {
    mesh->info.sethmin = 1;
  }
  /* Don't reset the hmax value computed when unscaling the mesh */
  if ( !parmesh->info.sethmax ) {
    mesh->info.sethmax = 1;
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
  PMMG_setfunc(parmesh);

  if ( !MMG3D_tetraQual( mesh, met, 0 ) ) {
    return PMMG_STRONGFAILURE;
  }

  if ( parmesh->info.imprim > PMMG_VERB_ITWAVES && (!mesh->info.iso) && met->m ) {
#warning: Luca: check this function
    MMG3D_prilen(mesh,met,0);
  }

  /** Mesh unscaling */
  if ( !MMG5_unscaleMesh(mesh,met,NULL) ) {
    return PMMG_STRONGFAILURE;
  }

  /** Mesh analysis: check triangles, create xtetras */
  if ( !PMMG_analys_tria(parmesh,mesh) ) {
    return PMMG_STRONGFAILURE;
  }
  /* For both API modes, build communicators indices */
  switch( parmesh->info.API_mode ) {
    case PMMG_APIDISTRIB_faces :
      /* Set face communicators indexing */
      if( !PMMG_build_faceCommIndex( parmesh ) ) return 0;
      break;
    case PMMG_APIDISTRIB_nodes :
      /* Set node communicators indexing */
      if( !PMMG_build_nodeCommIndex( parmesh ) ) return 0;
      break;
  }


  if( parmesh->info.API_mode == PMMG_APIDISTRIB_faces ) {
    /* Convert tria index into iel face index (it needs a valid cc field in
     * each tria), and tag xtetra face as PARBDY before the tag is transmitted
     * to edges and nodes */
    PMMG_tria2elmFace_coords( parmesh );
  }

  /* For both API modes, build communicators indices and set xtetra as PARBDY */
  switch( parmesh->info.API_mode ) {
    case PMMG_APIDISTRIB_faces :
      /* Build node communicators from face ones (here because the (mesh needs
       * to be unscaled) */
      PMMG_parmesh_ext_comm_free( parmesh,parmesh->ext_node_comm,parmesh->next_node_comm);
      PMMG_DEL_MEM(parmesh, parmesh->ext_node_comm,PMMG_Ext_comm,"ext node comm");
      parmesh->next_node_comm = 0;
      PMMG_DEL_MEM(parmesh, parmesh->int_node_comm,PMMG_Int_comm,"int node comm");
      PMMG_CALLOC(parmesh,parmesh->int_node_comm,1,PMMG_Int_comm,"int node comm",return 0);
      if ( !PMMG_build_nodeCommFromFaces(parmesh) ) return PMMG_STRONGFAILURE;
      break;
    case PMMG_APIDISTRIB_nodes :
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

  /** Geometrical mesh analysis */
  if ( !PMMG_analys(parmesh,mesh) ) {
    return PMMG_STRONGFAILURE;
  }


  if ( !PMMG_qualhisto(parmesh,PMMG_INQUA,0) ) {
    return PMMG_STRONGFAILURE;
  }

  MMG5_DEL_MEM(mesh,mesh->tria);
  mesh->nt = 0;

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

  } else {
    ier = PMMG_SUCCESS;
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

  /** Function setters (must be assigned before quality computation) */
  if( parmesh->myrank != parmesh->info.root ) {
    mesh = parmesh->listgrp[0].mesh;
    met  = parmesh->listgrp[0].met;
    MMG3D_Set_commonFunc();
    MMG3D_setfunc(mesh,met);
    PMMG_setfunc(parmesh);
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

  /* I/O check: if an input field name is provided but the output one is not,
   compute automatically an output solution field name. */
  if ( parmesh->fieldin &&  *parmesh->fieldin ) {
    ier = PMMG_Set_outputSolsName(parmesh,NULL);
    if ( !ier ) {
      fprintf(stdout,"  ## Warning: %s: rank %d: an input field name is"
              " provided without an output one.\n"
              "            : the saving process may fail.\n",
              __func__,parmesh->myrank);
    }
  }


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
  case ( PMMG_UNSET ):
    /* No output */
    break;
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
      // Add saving at parallel format
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
  int              ier,iresult,ierlib,npmax,nemax,xpmax,xtmax;
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

  /* I/O check: if an input field name is provided but the output one is not,
   compute automatically an output solution field name. */
  if ( parmesh->fieldin &&  *parmesh->fieldin ) {
    ier = PMMG_Set_outputSolsName(parmesh,NULL);
    if ( !ier ) {
      fprintf(stdout,"  ## Warning: %s: rank %d: an input field name is"
              " provided without an output one.\n"
              "            : the saving process may fail.\n",
              __func__,parmesh->myrank);
    }
  }

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

  /** Check input data */
  ier = PMMG_check_inputData( parmesh );
  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !iresult ) return PMMG_LOWFAILURE;

  /** Mesh preprocessing: set function pointers, scale mesh, perform mesh
   * analysis and display length and quality histos. */
  if( parmesh->myrank == parmesh->info.root ) {
    if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
      fprintf(stdout,"\n  -- ANALYSIS" );
    }
    ier = PMMG_preprocessMesh( parmesh );
    if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
      fprintf(stdout,"\n  -- ANALYSIS COMPLETED\n");
    }

    mesh = parmesh->listgrp[0].mesh;
    met  = parmesh->listgrp[0].met;
    if ( (ier==PMMG_STRONGFAILURE) && MMG5_unscaleMesh( mesh, met, NULL ) ) {
      ier = PMMG_LOWFAILURE;
    }

    /* Memory repartition */
    if ( !PMMG_parmesh_updateMemMax( parmesh,50,1 ) ) ier = 3;

  } else {
    ier = PMMG_SUCCESS;
  }

  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MAX, parmesh->comm );
  if ( iresult!=PMMG_SUCCESS ) {
    return iresult;
  }

  /** Send mesh partionning to other procs */
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    fprintf(stdout,"\n  -- PARTITIONING" );
  }
  if ( !PMMG_distribute_mesh( parmesh ) ) {
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_LOWFAILURE);
  }
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    fprintf(stdout,"\n  -- PARTITIONING COMPLETED\n");
  }

  /** Function setters (must be assigned before quality computation) */
  if( parmesh->myrank != parmesh->info.root ) {
    mesh = parmesh->listgrp[0].mesh;
    met  = parmesh->listgrp[0].met;
    MMG3D_Set_commonFunc();
    MMG3D_setfunc(mesh,met);
  }

  /* Memory repartition */
  if ( !PMMG_parmesh_updateMemMax( parmesh,50,1 ) ) return 3;

  iresult = PMMG_SUCCESS;
  return iresult;
}

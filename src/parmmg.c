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
 * \file parmmg.c
 * \brief main file for the parmmg application
 * \author Cecile Dobrzynski (Bx INP/Inria)
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
 * \todo refactoring to improve readibility
 */
int main( int argc, char *argv[] )
{
  PMMG_pParMesh parmesh = NULL;
  PMMG_pGrp     grp;
  int           rank;
  int           ier,iermesh,iresult,ierSave,fmtin,fmtout;
  int8_t        tim,distributedInput;
  char          stim[32],*ptr;

  // Shared memory communicator: processes that are on the same node, sharing
  //    local memory and can potentially communicate without using the network
  MPI_Comm comm_shm = 0;
  int      rank_shm = 0;

  /** Initializations: MPI, mesh, and memory */
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if ( !rank ) {
    fprintf(stdout,"  -- PARMMG, Release %s (%s) \n",PMMG_VERSION_RELEASE,PMMG_RELEASE_DATE);
    fprintf(stdout,"     %s\n",PMMG_COPYRIGHT);
    fprintf(stdout,"     %s %s\n\n",__DATE__,__TIME__);

    fprintf(stdout,"  -- MMG3D,    Release %s (%s) \n",MMG_VERSION_RELEASE,MMG_RELEASE_DATE);
    fprintf(stdout,"     %s\n",MMG_COPYRIGHT);
  }

  if ( !rank ) {
    atexit(PMMG_endcod);
  }

  tminit(PMMG_ctim,TIMEMAX);
  chrono(ON,&PMMG_ctim[0]);

  /* Allocate the main pmmg struct and assign default values */
  if ( 1 != PMMG_Init_parMesh( PMMG_ARG_start,
                               PMMG_ARG_ppParMesh,&parmesh,
                               PMMG_ARG_pLs,
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
                             MMG5_ARG_ppLs,   &parmesh->listgrp[0].ls,
                             MMG5_ARG_end) )
    PMMG_RETURN_AND_FREE( parmesh, PMMG_STRONGFAILURE );

  /* Init memMax sizes. Only one mesh for now => pmmg structs do not need much */
  if ( !PMMG_parmesh_SetMemMax(parmesh) )
    PMMG_RETURN_AND_FREE( parmesh, PMMG_STRONGFAILURE );

  /* Read command line */
  if ( 1 != PMMG_parsar( argc, argv, parmesh ) )
    PMMG_RETURN_AND_FREE( parmesh, PMMG_STRONGFAILURE );

  if ( parmesh->ddebug ) {
    MPI_Comm_split_type( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                         &comm_shm );
    MPI_Comm_rank( comm_shm, &rank_shm );
    MPI_Comm_size( comm_shm, &parmesh->size_shm );

    if ( !rank_shm )
      printf("\n     %d MPI PROCESSES (%d ON LOCAL NODE):\n",parmesh->nprocs,
             parmesh->size_shm);

    printf("         MPI RANK %d (LOCAL RANK %d)\n", parmesh->myrank,rank_shm );
  }

  /** load data */
  tim = 1;
  chrono(ON,&PMMG_ctim[tim]);
  if ( rank==parmesh->info.root && parmesh->info.imprim > PMMG_VERB_NO ) {
    fprintf(stdout,"\n  -- INPUT DATA: LOADING MESH ON RANK %d\n",
            parmesh->info.root);
  }

  grp = &parmesh->listgrp[0];
  iresult = 1;
  ier     = 1;

  ptr   = MMG5_Get_filenameExt(parmesh->meshin);

  fmtin = PMMG_Get_format(ptr,MMG5_FMT_MeditASCII);

  /* Compute default output format */
  ptr = MMG5_Get_filenameExt(parmesh->meshout);

  /* Format from output mesh name */
  fmtout = PMMG_Get_format(ptr,fmtin);

  distributedInput = 0;

  switch ( fmtin ) {
  case ( MMG5_FMT_MeditASCII ): case ( MMG5_FMT_MeditBinary ):

    // Algiane: Dirty (to be discussed, I don't have a clean solution)
    iermesh = PMMG_loadMesh_centralized(parmesh,parmesh->meshin);
    MPI_Bcast( &iermesh,     1, MPI_INT, parmesh->info.root, parmesh->comm );

    if ( 1 != iermesh ) {
      /* try to load distributed mesh */
      int ier_loc = PMMG_loadMesh_distributed(parmesh,parmesh->meshin);
      MPI_Allreduce( &ier_loc, &iermesh, 1, MPI_INT, MPI_MIN, parmesh->comm);
      if ( iermesh == 1 ) {
        distributedInput = 1;
        if ( parmesh->info.fmtout == PMMG_FMT_Centralized ) {
          /* Force centralization for output medit format: fmtout contains
           * already suitable info */
          parmesh->info.fmtout = fmtout;
        }
        else if ( fmtout == MMG5_FMT_MeditASCII ) {
          /* Distributed ASCII Medit */
          parmesh->info.fmtout = PMMG_FMT_DistributedMeditASCII;
        }
        else if ( fmtout == MMG5_FMT_MeditBinary ) {
          /* Distributed Binary Medit */
          parmesh->info.fmtout = PMMG_FMT_DistributedMeditBinary;
        }
        else if ( parmesh->info.fmtout != PMMG_UNSET ) {
          /* other format output: store format into parmesh->info.fmtout */
          parmesh->info.fmtout = fmtout;
        }
      }
    }
    else {
      /* Centralized mesh */
      if ( parmesh->info.fmtout == PMMG_FMT_Distributed ) {
        /* Force distribution */
        if ( fmtout == MMG5_FMT_MeditASCII ) {
          /* Distributed ASCII Medit */
          parmesh->info.fmtout = PMMG_FMT_DistributedMeditASCII;
        }
        else if ( fmtout == MMG5_FMT_MeditBinary ) {
          /* Distributed Binary Medit */
          parmesh->info.fmtout = PMMG_FMT_DistributedMeditBinary;
        }
        else if ( parmesh->info.fmtout != PMMG_UNSET ) {
          /* format is already good: store it in parmesh->info.fmtout */
          parmesh->info.fmtout = fmtout;
        }
      }
      else if ( parmesh->info.fmtout != PMMG_UNSET ){
        /* format is already good: store it in parmesh->info.fmtout */
        parmesh->info.fmtout = fmtout;
      }
    }

    if ( 1 != iermesh && rank == parmesh->info.root ) {
      if ( iermesh==0 ) {
        fprintf(stderr,"  ** %s  NOT FOUND.\n",parmesh->meshin);
        fprintf(stderr,"  ** UNABLE TO OPEN INPUT FILE.\n");
      }
      ier = 0;
      goto check_mesh_loading;
    }

    if ( grp->mesh->info.lag >= 0 ) {

      if ( rank == parmesh->info.root ) {
        fprintf(stderr,"\n  ## ERROR: LAGRANGIAN MODE I/O NOT YET IMPLEMENTED\n");
      }
      ier = 0;
      goto check_mesh_loading;

      // Check what happens: where is stored the input displacement name?
      // Mmg code
      // /* In Lagrangian mode, the name of the displacement file has been parsed in ls */
      //if ( !MMG3D_Set_inputSolName(grp->mesh,grp->disp,grp->ls->namein) ) {
      //  MMG5_RETURN_AND_FREE(mesh,met,ls,disp,MMG5_STRONGFAILURE);
      //}
      //MMG5_DEL_MEM(mesh,ls->namein);
    }

    if ( grp->mesh->info.lag >= 0 || grp->mesh->info.iso ) {
      /* displacement or isovalue are mandatory */
      if( !distributedInput ) {
        iermesh = ( PMMG_loadSol_centralized( parmesh, parmesh->lsin ) );
      }
      else {
        int ier_loc = PMMG_loadSol_distributed( parmesh, parmesh->lsin );
        MPI_Allreduce( &ier_loc, &iermesh, 1, MPI_INT, MPI_MIN, parmesh->comm);
      }
      if ( iermesh < 1 ) {
        if ( rank == parmesh->info.root ) {
          fprintf(stderr,"\n  ## ERROR: UNABLE TO LOAD SOLUTION FILE.\n");
        }
        ier = 0;
        goto check_mesh_loading;
      }
    }
    else {
      /* Facultative metric */
      if ( !distributedInput ) {
        iermesh = PMMG_loadMet_centralized( parmesh, parmesh->metin );
      }
      else {
        int ier_loc = PMMG_loadMet_distributed( parmesh, parmesh->metin );
        MPI_Allreduce( &ier_loc, &iermesh, 1, MPI_INT, MPI_MIN, parmesh->comm);
      }
      if ( -1 == iermesh  ) {
        if ( rank == parmesh->info.root ) {
          fprintf(stderr,"\n  ## ERROR: WRONG DATA TYPE OR WRONG SOLUTION NUMBER.\n");
        }
        ier = 0;
        goto check_mesh_loading;
      }
    }
    /* In iso mode: read metric if any */
    if ( grp->mesh->info.iso) {
      if ( parmesh->metin ) {
        if ( !distributedInput ) {
          iermesh = PMMG_loadMet_centralized( parmesh, parmesh->metin );
        }
        else {
          int ier_loc = PMMG_loadMet_distributed( parmesh, parmesh->metin );
          MPI_Allreduce( &ier_loc, &iermesh, 1, MPI_INT, MPI_MIN, parmesh->comm);
        }
        if ( -1 == iermesh ) {
          if ( rank == parmesh->info.root ) {
            fprintf(stderr,"\n  ## ERROR: UNABLE TO LOAD METRIC.\n");
          }
          ier = 0;
          goto check_mesh_loading;
        }
      }
      else {
        /* Give a name to the metric if not provided for distributed metric output */
        if ( !MMG5_Set_inputSolName(grp->mesh,grp->met,"") ) {
          fprintf(stdout,"  ## WARNING: Unable to give a name to the metric.\n");
        }
        else {
          ier = PMMG_Set_name(parmesh,&parmesh->metin,grp->met->namein,"mesh.sol");
          if (!ier) {
            fprintf(stdout,"  ## ERROR: Unable to give a name to the metric.\n");
            PMMG_RETURN_AND_FREE( parmesh, PMMG_LOWFAILURE );
          }
        }
      }
    }

    /* Read solutions field if provided */
    if ( parmesh->fieldin && *parmesh->fieldin ) {
      if ( !distributedInput ) {
        iermesh = PMMG_loadAllSols_centralized(parmesh,parmesh->fieldin);
      }
      else {
        int ier_loc = PMMG_loadAllSols_distributed(parmesh,parmesh->fieldin);
        MPI_Allreduce( &ier_loc, &iermesh, 1, MPI_INT, MPI_MIN, parmesh->comm);
      }
      if ( iermesh < 1 ) {
        if ( rank == parmesh->info.root ) {
          fprintf(stderr,"\n  ## ERROR: UNABLE TO LOAD FIELDS.\n");
        }
        ier = 0;
        goto check_mesh_loading;
      }
    }
    break;

  case PMMG_FMT_HDF5:
    ier = PMMG_loadMesh_hdf5( parmesh, parmesh->meshin );
    parmesh->info.fmtout = fmtout;
    distributedInput = 1;

    break;

  default:
    if ( rank == parmesh->info.root ) {
      fprintf(stderr,"  ** I/O AT FORMAT %s NOT IMPLEMENTED.\n",MMG5_Get_formatName(fmtin) );
    }
    ier = 0;
    goto check_mesh_loading;
  }

  if ( !PMMG_parsop(parmesh) ) {
    ier = 0;
    goto check_mesh_loading;
  }

  chrono(OFF,&PMMG_ctim[tim]);
  if ( parmesh->info.imprim > PMMG_VERB_NO ) {
    printim(PMMG_ctim[tim].gdif,stim);
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
  else if ( !distributedInput ) {
    /* Parallel remeshing starting from a centralized mesh */
    ier = PMMG_parmmg_centralized(parmesh);
  }
  else {
    /* Parallel remeshing starting from a distributed mesh */
    ier = PMMG_parmmg_distributed(parmesh);
  }

  /** Check result and save output files */
  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MAX, parmesh->comm );
  if ( iresult == PMMG_STRONGFAILURE ) { PMMG_RETURN_AND_FREE( parmesh, ier ); }

  tim = 2;
  chrono(ON,&PMMG_ctim[tim]);

  grp = &parmesh->listgrp[0];
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    if ( parmesh->info.fmtout == PMMG_FMT_DistributedMeditASCII ) {
      char *basename = MMG5_Remove_ext ( parmesh->meshout,".mesh" );
      fprintf(stdout,"\n  -- WRITING DATA FILE %s.<rankid>.mesh\n",basename);
      MMG5_SAFE_FREE ( basename );
    }
    else if ( parmesh->info.fmtout == PMMG_FMT_DistributedMeditBinary ) {
      char *basename = MMG5_Remove_ext ( parmesh->meshout,".meshb" );
      fprintf(stdout,"\n  -- WRITING DATA FILE %s.<rankid>.meshb\n",basename);
      MMG5_SAFE_FREE ( basename );
    }
    else if ( parmesh->info.fmtout == MMG5_FMT_VtkPvtu ) {
      char *basename = MMG5_Remove_ext ( parmesh->meshout,".pvtu" );
      int i, rename=0;
      for(i=0;basename[i]!='\0';i++) {
        if(basename[i]=='.') {
          basename[i] = '-';
          rename = 1;
        }
      }
      fprintf(stdout,"\n  -- WRITING DATA FILES %s.pvtu\n",basename);
      if (rename) fprintf(stdout,"       ## WARNING: Filename has been changed: "
              "%s => %s.pvtu\n",parmesh->meshout,basename);
      MMG5_SAFE_FREE ( basename );
    }
    else {
      fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",parmesh->meshout);
    }
    if (grp->field) {
      fprintf(stdout,"       Writing mesh, metric and fields.\n");
    }
    else {
      fprintf(stdout,"       Writing mesh and metric.\n");
    }
  }

  /* Initialize ierSave to 1 because for centralized output it will be assigned
   * only on root rank. */
  ierSave = 1;

  if ( parmesh->listgrp && parmesh->listgrp[0].mesh ) {
    grp = &parmesh->listgrp[0];

    switch ( parmesh->info.fmtout ) {
    case ( PMMG_UNSET ):
      /* No output */
      if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
        printf("     ... SKIPPING!\n");
      }
      break;

    case ( MMG5_FMT_VtkPvtu ):
      if (grp->field) {
        ier = PMMG_savePvtuMesh_and_allData(parmesh,parmesh->meshout);
      }
      else{
        ier = PMMG_savePvtuMesh(parmesh,parmesh->meshout);
      }
      MPI_Allreduce( &ier, &ierSave, 1, MPI_INT, MPI_MIN, parmesh->comm );
      break;

    case ( MMG5_FMT_GmshASCII ): case ( MMG5_FMT_GmshBinary ):
    case ( MMG5_FMT_VtkVtu ):
    case ( MMG5_FMT_VtkVtk ):
      printf("  ## Error: Output format not yet implemented.\n");
      PMMG_RETURN_AND_FREE(parmesh,PMMG_STRONGFAILURE);
      break;

    case ( PMMG_FMT_Distributed ):
    case ( PMMG_FMT_DistributedMeditASCII):
      ptr = MMG5_Get_filenameExt( parmesh->meshout );
      if ( (!ptr) || strcmp(ptr,".mesh") ) {
        /* Add .mesh extension to avoid saving at binary format */
        PMMG_REALLOC(parmesh,parmesh->meshout,strlen(parmesh->meshout)+6,
                     strlen(parmesh->meshout)+1,char,"",);
        strcat(parmesh->meshout,".mesh");
      }

    case ( PMMG_FMT_DistributedMeditBinary):
      assert ( parmesh->meshout );

      ier = PMMG_saveMesh_distributed(parmesh,parmesh->meshout);
      if ( ier ) {
        if ( parmesh->listgrp[0].met && parmesh->listgrp[0].met->m ) {
          ier = PMMG_saveMet_distributed(parmesh,parmesh->metout);
        }
      }
      if ( ier &&  grp->field ) {
        ier = PMMG_saveAllSols_distributed(parmesh,parmesh->fieldout);
      }
      MPI_Allreduce( &ier, &ierSave, 1, MPI_INT, MPI_MIN, parmesh->comm );

      break;

    case ( PMMG_FMT_HDF5 ):
      ier = PMMG_saveMesh_hdf5(parmesh,parmesh->meshout);
      MPI_Allreduce( &ier, &ierSave, 1, MPI_INT, MPI_MIN, parmesh->comm );

      break;

    default:
      ierSave = PMMG_saveMesh_centralized(parmesh,parmesh->meshout);

      if ( ierSave && parmesh->listgrp[0].met && parmesh->listgrp[0].met->m ) {
        ierSave = PMMG_saveMet_centralized(parmesh,parmesh->metout);
      }

      if ( ierSave && grp->field ) {
        ierSave = PMMG_saveAllSols_centralized(parmesh,parmesh->fieldout);
      }
      break;
    }
  }

  /* Check output success */
  if ( ierSave<1 ) {
    PMMG_RETURN_AND_FREE(parmesh,PMMG_STRONGFAILURE);
  }

  chrono(OFF,&PMMG_ctim[tim]);
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"  -- WRITING COMPLETED\n");

  PMMG_RETURN_AND_FREE( parmesh, iresult );
}

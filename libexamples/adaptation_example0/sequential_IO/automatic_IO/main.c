/**
 * Example of use of the parmmg library (basic use of mesh adaptation)
 *
 * \author Algiane Froehly (InriaSoft)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>

/** Include the parmmg library hader file */
// if the header file is in the "include" directory
// #include "libparmmg.h"
// if the header file is in "include/parmmg"
#include "parmmg/libparmmg.h"

int main(int argc,char *argv[]) {
  PMMG_pParMesh   parmesh;
  int             ier,rank,i,nsols;
  char            *filename,*metname,*solname,*fileout,*tmp;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if ( !rank ) fprintf(stdout,"  -- TEST PARMMGLIB \n");

  solname = NULL;
  metname = NULL;
  tmp     = NULL;

  if ( (argc<3) && !rank ) {
    printf(" Usage: %s filein fileout [[-sol metfile]/[-met metfile]] [-field solfile] \n",argv[0]);
    return 1;
  }

  /* Name and path of the mesh file */
  filename = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
  if ( filename == NULL ) {
    perror("  ## Memory problem: calloc");
    MPI_Abort(MPI_COMM_WORLD,2);
    exit(EXIT_FAILURE);
  }
  strcpy(filename,argv[1]);

  fileout = (char *) calloc(strlen(argv[2]) + 1, sizeof(char));
  if ( fileout == NULL ) {
    perror("  ## Memory problem: calloc");
    MPI_Abort(MPI_COMM_WORLD,2);
    exit(EXIT_FAILURE);
  }
  strcpy(fileout,argv[2]);

  i = 2;
  while ( ++i<argc ) {
    tmp = (char *) calloc(strlen(argv[i]) + 1, sizeof(char));
    strcpy(tmp,argv[i]);

    if ( (!strcmp(tmp, "-met"))||(!strcmp(tmp, "-sol")) ) {
      metname = (char *) calloc(strlen(argv[i+1]) + 1, sizeof(char));
      strcpy(metname,argv[i+1]);
      ++i;
    }

    else if ( !strcmp(tmp,"-field") ) {
      solname = (char *) calloc(strlen(argv[i+1]) + 1, sizeof(char));
      strcpy(solname,argv[i+1]);
      ++i;
    }

    else {
      printf("Unexpected argument: %s \n",tmp);
      MPI_Abort(MPI_COMM_WORLD,2);
      return 1;
    }
  }

  /** ------------------------------ STEP   I -------------------------- */
  /** 1) Initialisation of th parmesh structures */
  /* args of InitMesh:
   * PMMG_ARG_start: we start to give the args of a variadic func
   * PMMG_ARG_ppParMesh: next arg will be a pointer over a PMMG_pParMesh
   * &parmesh: pointer toward your PMMG_pParMesh
   * MMG5_ARG_pMesh: initialization of a mesh inside the parmesh.
   * MMG5_ARG_pMet: init a metric inside the parmesh
   * PMMG_ARG_dim: next arg will be the mesh dimension
   * 3: mesh dimension
   * PMMG_MPIComm: next arg will be the MPI COmmunicator
   * MPI_COMM_WORLD: MPI communicator
   *
   */
  parmesh = NULL;

  PMMG_Init_parMesh(PMMG_ARG_start,
                    PMMG_ARG_ppParMesh,&parmesh,
                    PMMG_ARG_pMesh,PMMG_ARG_pMet,
                    PMMG_ARG_dim,3,PMMG_ARG_MPIComm,MPI_COMM_WORLD,
                    PMMG_ARG_end);

  /** 2) Build mesh in MMG5 format */
  /** Two solutions: just use the PMMG_loadMesh_centralized function that will
      read a .mesh(b) file formatted or manually set your mesh using the
      PMMG_Set* functions */

  /** with PMMG_loadMesh_centralized function */
  if ( PMMG_loadMesh_centralized(parmesh,filename) != 1 ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  /** 3) Try to load a metric in PMMG format */
  /** Two solutions: just use the PMMG_loadMet_centralized function that will
      read a .sol(b) file formatted or manually set your metric using the PMMG_Set*
      functions */

  /** With PMMG_loadMet_centralized function */
  if ( metname ) {
    if ( PMMG_loadMet_centralized(parmesh,metname) != 1 ) {
      printf("Unable to load metric file.\n");
      MPI_Abort(MPI_COMM_WORLD,2);
      exit(EXIT_FAILURE);
    }
  }

  /** 4) Build solutions in PMMG format */
  /** Two solutions: just use the PMMG_loadAllSols_centralized function that
      will read a .sol(b) file formatted or manually set your solutions using
      the PMMG_Set* functions */

  /** With PMMG_loadAllSols_centralized function */

  if ( solname ) {
    if ( PMMG_loadAllSols_centralized(parmesh,solname) != 1 ) {
      printf("Unable to load solutions file.\n");
      MPI_Abort(MPI_COMM_WORLD,3);
      exit(EXIT_FAILURE);
    }
    /* Compute automatically output solution name from the output mesh path
       and the input field name */
    PMMG_Set_outputMeshName(parmesh,fileout);
    PMMG_Set_inputSolsName(parmesh,solname);
    PMMG_Set_outputSolsName(parmesh,NULL);
  }
  /** ------------------------------ STEP  II -------------------------- */
  /** remesh function */
  ier = PMMG_parmmglib_centralized(parmesh);

  if ( ier != PMMG_STRONGFAILURE ) {
    /** ------------------------------ STEP III -------------------------- */
    /** get results */
    /** Two solutions: just use the PMMG_saveMesh/PMMG_saveSol functions
        that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
        using the PMMG_getMesh/PMMG_getSol functions */

    /** 1) Automatically save the mesh */
    if ( PMMG_saveMesh_centralized(parmesh,fileout) != 1 ) {
      fprintf(stdout,"UNABLE TO SAVE MESH\n");
      ier = PMMG_STRONGFAILURE;
    }

    /** 2) Automatically save the metric */
    if ( PMMG_saveMet_centralized(parmesh,fileout) != 1 ) {
      fprintf(stdout,"UNABLE TO SAVE METRIC\n");
      ier = PMMG_LOWFAILURE;
    }

    /** 3) Automatically save the solutions if needed */
    PMMG_Get_solsAtVerticesSize(parmesh,&nsols,NULL,NULL);
    if ( nsols ) {
      /* Use default output solution name */
      if ( PMMG_saveAllSols_centralized(parmesh,NULL) != 1 ) {
        fprintf(stdout,"UNABLE TO SAVE SOLUTIONS\n");
        ier = PMMG_LOWFAILURE;
      }
    }
  }
  else {
    fprintf(stdout,"BAD ENDING OF PARMMGLIB: UNABLE TO SAVE MESH\n");
  }

  /** 4) Free the PMMG5 structures */
  PMMG_Free_all(PMMG_ARG_start,
                PMMG_ARG_ppParMesh,&parmesh,
                PMMG_ARG_end);

  free(filename);
  filename = NULL;

  free(fileout);
  fileout = NULL;

  MPI_Finalize();

  return ier;
}

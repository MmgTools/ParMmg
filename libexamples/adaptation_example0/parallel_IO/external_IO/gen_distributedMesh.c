/**
 * Example of use of the parmmg library with a distributed input mesh (basic
 * use of mesh adaptation)
 *
 * This example show how to set and get parallel mesh interfaces using the
 * ParMmg setters and getters, starting from a global mesh which is
 * automatically partitioned and distributed among the processes.
 * Depending on the command line option "niter", the programs performs a dry run
 * of ParMMG without remeshing steps, to the purpose of checking parallel
 * interface consistency, or a true run of ParMMG with parallel remeshing.
 * Depending on the command line option "API_mode", either face or node
 * interfaces are set.

 * \author Luca Cirrottola (Inria)
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
#include "libparmmg.h"


int main(int argc,char *argv[]) {
  PMMG_pParMesh   parmesh;
  int             rank,nprocs;
  int             opt,API_mode;
  char            *filename,*fileout;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

  if ( !rank ) fprintf(stdout,"  -- TEST PARMMGLIB \n");

  if ( (argc!=4) && !rank ) {
    printf(" Usage: %s filein fileout io_option\n",argv[0]);
    printf("     API_mode = 0   to Get/Set the parallel interfaces through triangles\n");
    printf("     API_mode = 1   to Get/Set the parallel interfaces through nodes\n");
    return 1;
  }

  /* Name and path of the mesh file */
  filename = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
  if ( filename == NULL ) {
    perror("  ## Memory problem: calloc");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  strcpy(filename,argv[1]);

  fileout = (char *) calloc(strlen(argv[2]) + 9 + 4, sizeof(char));
  if ( fileout == NULL ) {
    perror("  ## Memory problem: calloc");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  strcpy(fileout,argv[2]);

  /* Option to Set mesh entities vertex by vertex */
  opt      = 1;

  /* Get API mode (face or node interfaces) */
  API_mode = atoi(argv[3]);

  /** ------------------------------ STEP   I -------------------------- */
  /** Each process loads a global mesh.
   */

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

  /** ------------------------------ STEP  II -------------------------- */
  /** Preprocess and partition the mesh.
   */
  // Perform 0 iteration => analysis only
  PMMG_Set_iparameter(parmesh,PMMG_IPARAM_niter,0);

  // Do not centralize meshes at the end of the library call
  PMMG_Set_iparameter(parmesh,PMMG_IPARAM_distributedOutput,1);

  // Set API mode
  PMMG_Set_iparameter(parmesh,PMMG_IPARAM_APImode,API_mode);

  PMMG_parmmglib_centralized(parmesh);

  /** ------------------------------ STEP  III ------------------------- */
  /** Get parallel interfaces and swap meshes, so that you can use meshIN to
   * initialize a new mesh in parmesh */

  /** 1) Recover parallel interfaces */
  /** save mesh and interfaces **/
  PMMG_saveMesh_distributed(parmesh,fileout);

  /** 5) Free the PMMG5 structures */
  PMMG_Free_all(PMMG_ARG_start,
                PMMG_ARG_ppParMesh,&parmesh,
                PMMG_ARG_end);

  free(filename);
  filename = NULL;

  free(fileout);
  fileout = NULL;

  MPI_Finalize();

  return PMMG_SUCCESS;
}

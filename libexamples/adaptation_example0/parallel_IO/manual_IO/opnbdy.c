/**
 * Example of use of the parmmg library (basic use of mesh adaptation).
 *
 * This example show how to set and get parallel mesh interfaces using the
 * ParMmg setters and getters, starting from a manual partitioning of a global
 * mesh.
 * Depending on the command line option "niter", the programs performs a dry run
 * of ParMMG without remeshing steps, to the purpose of checking parallel
 * interface consistency, or a true run of ParMMG with parallel remeshing.
 * Depending on the command line option "API_mode", either face or node
 * interfaces are set.
 *
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
// if the header file is in the "include" directory
#include "libparmmg.h"
// if the header file is in "include/parmmg"
//#include "parmmg/libparmmg.h"

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

void get_local_mesh(int np, int ne, int nt, int *pmask, int *inv_pmask,
                    int *emask, int *tmask, int *inv_tmask,
                    double *pcoor, double *pcoor_all, int *pref, int *pref_all,
                    int *evert, int *evert_all, int *eref, int *eref_all,
                    int *tvert, int *tvert_all, int *tref, int *tref_all,
                    double *met, double *met_all,int ncomm,
                    int *ntifc, int **ifc_tria_loc, int **ifc_tria_glob,
                    int *npifc, int **ifc_nodes_loc, int **ifc_nodes_glob);

/* Main program */
int main(int argc,char *argv[]) {
  PMMG_pParMesh   parmesh;
  int             ier,ierlib,rank,nprocs;
  char            *fileout;


  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

  if ( !rank ) {
    fprintf(stdout,"  -- TEST PARMMGLIB: manual mesh partition along a direction"
            " normal to the physical interface and opnbdy mode.\n");
    fprintf(stdout,"  -- The input matching centralized mesh (cube.mesh)"
            " is provided in the opnbdy directory for visualization "
            "(even if not used).\n"
            " -- Part of meshes local to each procs are also provided "
            " (cube-P0.mesh and cube-P1.mesh)\n" );
  }

  if ( (nprocs!=2) && !rank ) {
    printf(" ## Error: nprocs != 2.\n"
           "    This example is manually distributed on 2 procs so it must be"
           "    runned on the same number of procs.\n");
    return 0;
  }

  if ( (argc!=2) && !rank ) {
    printf(" Usage: %s fileout\n",argv[0]);
    return 0;
  }

  fileout = (char *) calloc(strlen(argv[1]) + 6 + 4, sizeof(char));
  if ( fileout == NULL ) {
    perror("  ## Memory problem: calloc");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  strcpy(fileout,argv[1]);

  /** ------------------------------ STEP   I -------------------------- */
  /** Each process initialize the parmesh structure, then stores its part of the
   * mesh
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

  /** 2) On each proc, give the local mesh in MMG5 format */
  int nVertices       = 8;
  int nTetrahedra     = 6;
  int nTriangles      = 14;
  int nPrisms         = 0;
  int nQuadrilaterals = 0;
  int nEdges          = 0;

  /* Set local mesh size (the same on the two procs in this particular case) */
  if ( PMMG_Set_meshSize(parmesh,nVertices,nTetrahedra,nPrisms,nTriangles,
                         nQuadrilaterals,nEdges) != 1 ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  if ( rank==0 ) {
    /* ------- Give local vertices */
    /* local point number 1: (1., 0., 0.) - ref 0 */
    if ( PMMG_Set_vertex(parmesh,1.,0.,0., 0, 1) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local point number 2: (0.5, 0., 0.) - ref 0 */
    if ( PMMG_Set_vertex(parmesh,0.5,0.,0., 0, 2) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local point number 3: (0.5, 0., 1.) - ref 0 */
    if ( PMMG_Set_vertex(parmesh,0.5,0.,1., 0, 3) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local point number 4: (0.5, 1., 1.) - ref 0 */
    if ( PMMG_Set_vertex(parmesh,0.5,1.,1., 0, 4) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local point number 5: ( 1., 0., 1.) - ref 0 */
    if ( PMMG_Set_vertex(parmesh,1.,0.,1., 0, 5) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local point number 6: ( 1., 1., 1.) - ref 0 */
    if ( PMMG_Set_vertex(parmesh,1.,1.,1., 0, 6) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local point number 7: (0.5, 1., 0.) - ref 0 */
    if ( PMMG_Set_vertex(parmesh,0.5,1.,0., 0, 7) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local point number 8: ( 1., 1., 0.) - ref 0 */
    if ( PMMG_Set_vertex(parmesh,1.,1.,0., 0, 8) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    /* ------- Give local tetra */
    /* local tetra number 1: 1, 2, 3, 4 - ref 0 */
    if ( PMMG_Set_tetrahedron(parmesh,1,2,3,4,0,1) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tetra number 2: 4, 5, 1, 6 - ref 0 */
    if ( PMMG_Set_tetrahedron(parmesh,4,5,1,6,0,2) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tetra number 3: 7, 1, 8, 4 - ref 0 */
    if ( PMMG_Set_tetrahedron(parmesh,7,1,8,4,0,3) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tetra number 4: 7, 4, 2, 1 - ref 0 */
    if ( PMMG_Set_tetrahedron(parmesh,7,4,2,1,0,4) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tetra number 5: 6, 1, 4, 8 - ref 0 */
    if ( PMMG_Set_tetrahedron(parmesh,6,1,4,8,0,5) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tetra number 6: 1, 3, 5, 4 - ref 0 */
    if ( PMMG_Set_tetrahedron(parmesh,1,3,5,4,0,6) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    /* ------- Give local tria: open boudary is provided with ref 10,
       interface triangles are provided with ref 0 */
    /* local tria number 1: 2, 3, 4 - ref 0 */
    if ( PMMG_Set_triangle(parmesh,2,3,4,0,1) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tria number 2: 2, 7, 4 - ref 0 */
    if ( PMMG_Set_triangle(parmesh,2,7,4,0,2) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tria number 3: 1, 3, 2 - ref 4 */
    if ( PMMG_Set_triangle(parmesh,1,3,2,4,3) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tria number 4: 5, 1, 6 - ref 4 */
    if ( PMMG_Set_triangle(parmesh,5,1,6,4,4) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tria number 5: 4, 5, 6 - ref 4 */
    if ( PMMG_Set_triangle(parmesh,4,5,6,4,5) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tria number 6: 7, 4, 8 - ref 4 */
    if ( PMMG_Set_triangle(parmesh,7,4,8,4,6) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tria number 7: 7, 8, 1 - ref 4 */
    if ( PMMG_Set_triangle(parmesh,7,8,1,4,7) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tria number 8: 7, 1, 2 - ref 4 */
    if ( PMMG_Set_triangle(parmesh,7,1,2,4,8) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tria number 9: 6, 8, 4 - ref 4 */
    if ( PMMG_Set_triangle(parmesh,6,8,4,4,9) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tria number 10: 6, 1, 8 - ref 4 */
    if ( PMMG_Set_triangle(parmesh,6,1,8,4,10) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tria number 11: 3, 5, 4 - ref 4 */
    if ( PMMG_Set_triangle(parmesh,3,5,4,4,11) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tria number 12: 1, 5, 3 - ref 4 */
    if ( PMMG_Set_triangle(parmesh,1,5,3,4,12) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tria number 13: 1, 4, 2 - ref 10 */
    if ( PMMG_Set_triangle(parmesh,1,4,2,10,13) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local tria number 14: 6, 4, 1 - ref 10 */
    if ( PMMG_Set_triangle(parmesh,6,4,1,10,14) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    //  MMG3D_saveMesh(parmesh->listgrp[0].mesh,"cube-P0.mesh");
  }
  else { /* rank == 1 */
    /* ------- Give local vertices */
    /* local point number 1 : ( 0  , 0.,   0.) - ref 0 */
    if ( PMMG_Set_vertex(parmesh,0  , 0., 0., 0,1 ) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local point number 2 : ( 0.5, 0.,   0.) - ref 0 */
    if ( PMMG_Set_vertex(parmesh,0.5, 0., 0., 0,2 ) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local point number 3 : ( 0.5, 0.,   1.) - ref 0 */
    if ( PMMG_Set_vertex(parmesh,0.5, 0., 1., 0,3 ) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local point number 4 : ( 0  , 0.,   1.) - ref 0 */
    if ( PMMG_Set_vertex(parmesh,0  , 0., 1., 0,4 ) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local point number 5 : ( 0  , 1.,   0.) - ref 0 */
    if ( PMMG_Set_vertex(parmesh,0  , 1., 0., 0,5 ) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local point number 6 : ( 0.5, 1.,   0.) - ref 0 */
    if ( PMMG_Set_vertex(parmesh,0.5, 1., 0., 0,6 ) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local point number 7 : ( 0.5, 1.,   1.) - ref 0 */
    if ( PMMG_Set_vertex(parmesh,0.5, 1., 1., 0,7 ) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    /* local point number 8 : ( 0  , 1.,   1.) - ref 0 */
    if ( PMMG_Set_vertex(parmesh,0  , 1., 1., 0,8 ) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    /* ------- Give local tetra */
    if ( PMMG_Set_tetrahedron(parmesh,1, 4, 2, 8, 0, 1) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_tetrahedron(parmesh,8, 3, 2, 7, 0, 2) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_tetrahedron(parmesh,5, 2, 6, 8, 0, 3) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_tetrahedron(parmesh,5, 8, 1, 2, 0, 4) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_tetrahedron(parmesh,7, 2, 8, 6, 0, 5) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_tetrahedron(parmesh,2, 4, 3, 8, 0, 6) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    /* ------- Give local tria */
    if ( PMMG_Set_triangle(parmesh,1, 4, 8, 3, 1) !=1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_triangle(parmesh,1, 2, 4, 3, 2) !=1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_triangle(parmesh,8, 3, 7, 3, 3) !=1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_triangle(parmesh,5, 8, 6, 3, 4) !=1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_triangle(parmesh,5, 6, 2, 3, 5) !=1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_triangle(parmesh,5, 2, 1, 3, 6) !=1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_triangle(parmesh,5, 1, 8, 3, 7) !=1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_triangle(parmesh,7, 6, 8, 3, 8) !=1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_triangle(parmesh,4, 3, 8, 3, 9) !=1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_triangle(parmesh,2, 3, 4, 3, 10) !=1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_triangle(parmesh,2, 3, 4, 3, 10) !=1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_triangle(parmesh,2, 6, 7, 0, 11) !=1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_triangle(parmesh,2, 3, 7, 0, 12) !=1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_triangle(parmesh,8, 7, 2,10,13) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( PMMG_Set_triangle(parmesh,8, 1, 2,10,14) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    // MMG3D_saveMesh(parmesh->listgrp[0].mesh,"cube-P1.mesh");
  } /* end local mesh setting */

  /** 3) Initialization of interface communicators in ParMMG.
   *     The user can choose between providing triangles (faces) interface
   *     information (through the PMMG_APIDISTRIB_faces parameter), or nodes
   *     interface information (through the PMMG_APIDISTRIB_nodes parameter).
   */

  /* Set API mode */
  if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_APImode, PMMG_APIDISTRIB_faces ) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  if( !rank ) printf("\n--- API mode: Setting face communicators\n");

  /** a) Set the number of pair of processors that communicate through faces */
  // We have only the pair P0-P1 that communicate (as we are on 2 procs)
  ier = PMMG_Set_numberOfFaceCommunicators(parmesh,1);

  /** b) For each interface (proc pair) seen by the current rank: */
  /* b.1) set nb. of triangles on the interface and the rank of the outward
   * proc */
  // On proc 0: we have the communicator number 0 to fill: we share 2 triangles with proc 1
  // On proc 1: we have the communicator number 0 to fill: we share 2 triangles with proc 0
  ier = PMMG_Set_ithFaceCommunicatorSize(parmesh, 0,(rank+1)%nprocs,2);

  /* b.2) Set local and global index for each entity on the interface */
  // The tria number 1 of proc 0 matches with the tria number 12 of proc 1 and
  // its global index is 2. The tria number 2 of proc 0 matches with the tria
  // number 11 of proc 1 and their global index is 1;

  if ( !rank ) {
    /* rank 0 */
    int local_index[2] = {1,2};
    int global_index[2] = {2,1};

    PMMG_Set_ithFaceCommunicator_faces(parmesh, 0,local_index,global_index, 1 );
  }
  else {
    /* rank 1 */
    int local_index[2] = {11,12};
    int global_index[2]= {1,2};

    PMMG_Set_ithFaceCommunicator_faces(parmesh, 0,local_index,global_index, 1 );
  }

  /** ------------------------------ STEP II -------------------------- */
  /** remesh step */

  /* Force centralized output */
  if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_distributedOutput,0) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  /* Set opnbdy preservation */
  if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_opnbdy, 1 ) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  /* Switch off surface adaptation until ready on opnbdy */
  if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_nosurf, 1 ) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  /* try to generate mesh of edge sizes of size 0.1 */
  if( !PMMG_Set_dparameter( parmesh, PMMG_DPARAM_hsiz, 0.1 ) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  /* remeshing function */
  ierlib = PMMG_parmmglib_distributed( parmesh );

  if ( ierlib == PMMG_SUCCESS ) {

    /** ------------------------------ STEP V  ---------------------------- */
    /** get results */
    /** Two solutions: just use the PMMG_saveMesh/PMMG_saveSol functions
        that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
        using the PMMG_getMesh/PMMG_getSol functions */

    /** 1) Save automatically the mesh (see the manual/main.c example for an
     * example of getting manually the mesh) */
    PMMG_saveMesh_centralized(parmesh,fileout);

  }
  else if ( ierlib == PMMG_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF PARMMGLIB: UNABLE TO SAVE MESH\n");
  }


  /** 5) Free the PMMG5 structures */
  PMMG_Free_all(PMMG_ARG_start,
                PMMG_ARG_ppParMesh,&parmesh,
                PMMG_ARG_end);

  free(fileout); fileout = NULL;

  MPI_Finalize();

  return ierlib;
}

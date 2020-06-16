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
#include "mmg3d.h" // for developpers only: to use MMG3D_bdryBuild

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

int main(int argc,char *argv[]) {
  PMMG_pParMesh   parmesh;
  MMG5_pMesh      mesh;
  int             ier,rank,nprocs,i;
  int             opt,API_mode;
  char            *filename,*metname,*solname,*fileout,*metout,*tmp;
  int             nVertices,nTetrahedra,nTriangles,nEdges;


  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

  if ( !rank ) fprintf(stdout,"  -- TEST PARMMGLIB \n");

  solname = NULL;
  metname = NULL;
  tmp     = NULL;

  if ( (argc!=4) && !rank ) {
    printf(" Usage: %s fileout io_option\n",argv[0]);
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
  sprintf(fileout, "%s-P%02d", fileout, rank );
  strcat(fileout,".mesh");


  metout = (char *) calloc(strlen(argv[2]) + 9 + 4, sizeof(char));
  if ( metout == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(metout,argv[2]);
  sprintf(metout, "%s-P%02d", metout, rank );
  strcat(metout,"-met.sol");

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

  /** 3) Try to load a metric in PMMG format */
  /** Two solutions: just use the PMMG_loadMet_centralized function that will
      read a .sol(b) file formatted or manually set your metric using the PMMG_Set*
      functions */

  /** With PMMG_loadMet_centralized function */
  if ( metname ) {
    if ( PMMG_loadMet_centralized(parmesh,filename) !=1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }

  /** 4) Build solutions in PMMG format */
  /** Two solutions: just use the PMMG_loadAllSols_centralized function that
      will read a .sol(b) file formatted or manually set your solutions using
      the PMMG_Set* functions */

  /** With PMMG_loadAllSols_centralized function */

  if ( solname ) {
    if ( PMMG_loadAllSols_centralized(parmesh,filename) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }

  /** ------------------------------ STEP  II -------------------------- */
  /** Preprocess and partition the mesh.
   */

  if( PMMG_distributeMesh_centralized(parmesh) != PMMG_SUCCESS ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }


  /** ------------------------------ STEP  III ------------------------- */
  /** Get parallel interfaces and swap meshes, so that you can use meshIN to
   * initialize a new mesh in parmesh */

  /* Create boundary entities */
  mesh = parmesh->listgrp[0].mesh;
  if( MMG3D_bdryBuild(mesh) == -1) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }


  /** 1) Recover parallel interfaces */

  int n_node_comm,n_face_comm,*nitem_node_comm,*nitem_face_comm;
  int *color_node, *color_face,**face_owner,nunique_face,ntot_face;
  int **idx_node_loc,**idx_node_glob,**node_owner,nunique_node,ntot_node;
  int **idx_face_loc,**idx_face_glob;
  int icomm;

  /* Get number of node interfaces */
  ier = PMMG_Get_numberOfNodeCommunicators(parmesh,&n_node_comm);
  if ( ier!=1 ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  /* Get outward proc rank and number of nodes on each interface */
  color_node      = (int *) malloc(n_node_comm*sizeof(int));
  nitem_node_comm = (int *) malloc(n_node_comm*sizeof(int));
  for( icomm = 0; icomm < n_node_comm; icomm++ ) {
    ier = PMMG_Get_ithNodeCommunicatorSize(parmesh, icomm,
                                           &color_node[icomm],
                                           &nitem_node_comm[icomm]);
    if ( ier !=1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }

  /* Get IDs of nodes on each interface */
  idx_node_loc  = (int **) malloc(n_node_comm*sizeof(int *));
  idx_node_glob = (int **) malloc(n_node_comm*sizeof(int *));
  node_owner    = (int **) malloc(n_node_comm*sizeof(int *));
  for( icomm = 0; icomm < n_node_comm; icomm++ ) {
    idx_node_loc[icomm]  = (int *) malloc(nitem_node_comm[icomm]*sizeof(int));
    idx_node_glob[icomm] = (int *) malloc(nitem_node_comm[icomm]*sizeof(int));
    node_owner[icomm]    = (int *) malloc(nitem_node_comm[icomm]*sizeof(int));
  }
  ier = PMMG_Get_NodeCommunicator_nodes(parmesh, idx_node_loc);
  if ( ier!=1 ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  /* Get number of face interfaces */
  ier = PMMG_Get_numberOfFaceCommunicators(parmesh,&n_face_comm);
  if ( ier!=1 ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  /* Get outward proc rank and number of faces on each interface */
  color_face      = (int *) malloc(n_face_comm*sizeof(int));
  nitem_face_comm = (int *) malloc(n_face_comm*sizeof(int));
  for( icomm = 0; icomm < n_face_comm; icomm++ ) {
    ier = PMMG_Get_ithFaceCommunicatorSize(parmesh, icomm,
                                           &color_face[icomm],
                                           &nitem_face_comm[icomm]);
    if ( ier!=1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }

  /* Get IDs of triangles on each interface */
  idx_face_loc  = (int **) malloc(n_face_comm*sizeof(int *));
  idx_face_glob = (int **) malloc(n_face_comm*sizeof(int *));
  face_owner    = (int **) malloc(n_face_comm*sizeof(int *));
  for( icomm = 0; icomm < n_face_comm; icomm++ ) {
    idx_face_loc[icomm]  = (int *) malloc(nitem_face_comm[icomm]*sizeof(int));
    idx_face_glob[icomm] = (int *) malloc(nitem_face_comm[icomm]*sizeof(int));
    face_owner[icomm]    = (int *) malloc(nitem_face_comm[icomm]*sizeof(int));
  }
  ier = PMMG_Get_FaceCommunicator_faces(parmesh, idx_face_loc);
  if ( ier!=1 ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  /* Get triangle nodes */
  nVertices   = 0;
  nTetrahedra = 0;
  nTriangles  = 0;
  nEdges      = 0;
  if ( PMMG_Get_meshSize(parmesh,&nVertices,&nTetrahedra,NULL,&nTriangles,NULL,
                         &nEdges) !=1 ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  int *ref       = (int*)calloc(nTriangles,sizeof(int));
  int *required  = (int*)calloc(nTriangles,sizeof(int));
  int *triaNodes = (int*)calloc(3*nTriangles,sizeof(int));

  if ( PMMG_Get_triangles(parmesh,triaNodes,ref,required) != 1 ) {
    fprintf(stderr,"Unable to get mesh triangles\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  /* Color interface triangles with a unique global enumeration that encompasses
   * all interface triangles currently present in the global mesh, and assign
   * a owner partition to each of them.
   */
  if( !PMMG_Get_FaceCommunicator_owners(parmesh,face_owner,idx_face_glob,&nunique_face,&ntot_face) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  /* Color interface nodes with a unique global enumeration that encompasses
   * interface nodes currently present in the global mesh, and assign
   * a owner partition to each of them.
   */
  if( !PMMG_Get_NodeCommunicator_owners(parmesh,node_owner,idx_node_glob,&nunique_node,&ntot_node) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }


  /** save mesh and interfaces **/
  char filemesh[256];
  sprintf(filemesh,"%s_out.%d.mesh",filename,parmesh->myrank);
  MMG3D_saveMesh(parmesh->listgrp[0].mesh,filemesh);

  FILE *fid;
  sprintf(filemesh,"%s_out.%d.mesh_parFaces",filename,parmesh->myrank); 
  fid = fopen(filemesh,"w");
  fprintf(fid,"ParallelTriangleCommunicators\n%d\n",n_face_comm);
  for( icomm = 0; icomm < n_face_comm; icomm++ )
    fprintf(fid,"%d %d\n",color_face[icomm],nitem_face_comm[icomm]);
  fprintf(fid,"\n"); 
  fprintf(fid,"ParallelTriangles\n");
  for( icomm = 0; icomm < n_face_comm; icomm++ )
    for( i = 0; i < nitem_face_comm[icomm]; i++ )
      fprintf(fid,"%d %d %d\n",idx_face_loc[icomm][i],idx_face_glob[icomm][i],color_face[icomm]);
  fclose(fid);

  sprintf(filemesh,"%s_out.%d.mesh_parNodes",filename,parmesh->myrank);
  fid = fopen(filemesh,"w");
  fprintf(fid,"ParallelVertexCommunicators\n%d\n",n_node_comm);
  for( icomm = 0; icomm < n_node_comm; icomm++ )
    fprintf(fid,"%d %d\n",color_node[icomm],nitem_node_comm[icomm]);
  fprintf(fid,"\n");
  fprintf(fid,"ParallelVertices\n");
  for( icomm = 0; icomm < n_node_comm; icomm++ )
    for( i = 0; i < nitem_node_comm[icomm]; i++ )
      fprintf(fid,"%d %d %d\n",idx_node_loc[icomm][i],idx_node_glob[icomm][i],color_node[icomm]);
  fclose(fid);

  /* Free arrays */
  for( icomm = 0; icomm < n_node_comm; icomm++ ) {
    free(idx_node_loc[icomm]);
    free(idx_node_glob[icomm]);
    free(node_owner[icomm]);
  }
  free(idx_node_loc);
  free(idx_node_glob);
  free(node_owner);

  for( icomm = 0; icomm < n_face_comm; icomm++ ) {
    free(idx_face_loc[icomm]);
    free(idx_face_glob[icomm]);
    free(face_owner[icomm]);
  }
  free(idx_face_loc);
  free(idx_face_glob);
  free(face_owner);

  /** 5) Free the PMMG5 structures */
  PMMG_Free_all(PMMG_ARG_start,
                PMMG_ARG_ppParMesh,&parmesh,
                PMMG_ARG_end);

  free(filename);
  filename = NULL;

  free(fileout);
  fileout = NULL;

  free(metout);
  metout = NULL;


  MPI_Finalize();

  return PMMG_SUCCESS;
}

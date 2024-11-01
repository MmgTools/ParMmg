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

/** Include the parmmg and libmmg3d library hader file */
#include "libparmmg.h"
#include "libmmg3d.h" // for developpers only: to use MMG3D_bdryBuild
#include "libmmg3d_private.h" // for developpers only: to use MMG3D_bdryBuild

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

int main(int argc,char *argv[]) {
  PMMG_pParMesh   parmesh;
  MMG5_pMesh      mesh,meshIN;
  MMG5_pSol       solIN;
  MMG5_pPoint     ppt;
  MMG5_pTria      ptt;
  MMG5_pTetra     pt;
  int             ip,ie,ier,ierlib,rank,nprocs,i,k;
  int             opt,API_mode,niter;
  char            *filename,*metname,*solname,*fileout,*metout,*tmp;
  FILE            *inm;
  int             pos,nreq,nc,nr;
  int             nVertices,nTetrahedra,nTriangles,nEdges;
  int             nodeGloNumber,nodeOwner;


  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );


  if ( !rank ) fprintf(stdout,"  -- TEST PARMMGLIB \n");

  solname = NULL;
  metname = NULL;
  tmp     = NULL;

  if ( (argc!=5) && !rank ) {
    printf(" Usage: %s fileout io_option\n",argv[0]);
    printf("     niter    = 0   to perform a dry run of Parmmg and check paralle interfaces construction\n");
    printf("     niter    = [n] to perform [n] iterations of remeshing\n");
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

  /* Get number of remeshing iterations */
  niter    = atoi(argv[3]);

  /* Get API mode (face or node interfaces) */
  API_mode = atoi(argv[4]);

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
  if ( metname )
    PMMG_loadMet_centralized(parmesh,filename);

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
  int **faceNodes;
  int icomm;

  /* Get number of node interfaces */
  ier = PMMG_Get_numberOfNodeCommunicators(parmesh,&n_node_comm);

  /* Get outward proc rank and number of nodes on each interface */
  color_node      = (int *) malloc(n_node_comm*sizeof(int));
  nitem_node_comm = (int *) malloc(n_node_comm*sizeof(int));
  for( icomm = 0; icomm < n_node_comm; icomm++ )
    ier = PMMG_Get_ithNodeCommunicatorSize(parmesh, icomm,
                                           &color_node[icomm],
                                           &nitem_node_comm[icomm]);


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

  /* Get number of face interfaces */
  ier = PMMG_Get_numberOfFaceCommunicators(parmesh,&n_face_comm);

  /* Get outward proc rank and number of faces on each interface */
  color_face      = (int *) malloc(n_face_comm*sizeof(int));
  nitem_face_comm = (int *) malloc(n_face_comm*sizeof(int));
  for( icomm = 0; icomm < n_face_comm; icomm++ )
    ier = PMMG_Get_ithFaceCommunicatorSize(parmesh, icomm,
                                           &color_face[icomm],
                                           &nitem_face_comm[icomm]);

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

  /* Get triangle nodes */
  nVertices   = 0;
  nTetrahedra = 0;
  nTriangles  = 0;
  nEdges      = 0;
  if ( PMMG_Get_meshSize(parmesh,&nVertices,&nTetrahedra,NULL,&nTriangles,NULL,
                         &nEdges) !=1 ) {
    ier = PMMG_STRONGFAILURE;
  }

  int *ref       = (int*)calloc(nTriangles,sizeof(int));
  int *required  = (int*)calloc(nTriangles,sizeof(int));
  int *triaNodes = (int*)calloc(3*nTriangles,sizeof(int));

  if ( PMMG_Get_triangles(parmesh,triaNodes,ref,required) != 1 ) {
    fprintf(stderr,"Unable to get mesh triangles\n");
    ier = PMMG_STRONGFAILURE;
  }

  /* Color interface triangles with a unique global enumeration that encompasses
   * all interface triangles currently present in the global mesh, and assign a
   * owner partition to each of them.
   */
  if( !PMMG_Get_FaceCommunicator_owners(parmesh,face_owner,idx_face_glob,&nunique_face,&ntot_face) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  /* Color interface nodes with a unique global enumeration that encompasses
   * all interface nodes currently present in the global mesh, and assign a
   * owner partition to each of them.
   */
  if( !PMMG_Get_NodeCommunicator_owners(parmesh,node_owner,idx_node_glob,&nunique_node,&ntot_node) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

/*
  printf("Rank %d, my nunique_face %d, ntot %d\n",rank,nunique_face,ntot_face);
  printf("Rank %d, my nunique_node %d, ntot %d\n",rank,nunique_node,ntot_node);

  for( icomm = 0; icomm < n_face_comm; icomm++ )
    for( i = 0; i < nitem_face_comm[icomm]; i++ )
      printf("IN rank %d comm %d color %d tria loc %d glob %d owner %d\n",rank,icomm,color_face[icomm],idx_face_loc[icomm][i],idx_face_glob[icomm][i],face_owner[icomm][i]);


  for( icomm = 0; icomm < n_node_comm; icomm++ )
    for( i = 0; i < nitem_node_comm[icomm]; i++ )
      printf("IN rank %d comm %d color %d node loc %d glob %d owner %d\n",rank,icomm,color_node[icomm],idx_node_loc[icomm][i],idx_node_glob[icomm][i],node_owner[icomm][i]);
*/


  /** 2) Create input mesh "meshIN" and swap it with parmmg mesh, so that
   *     ParMMG structures can be freed. */
  meshIN = NULL;
  solIN = NULL;
  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&meshIN,MMG5_ARG_ppMet,&solIN,
                  MMG5_ARG_end);

  if ( MMG3D_Set_meshSize(meshIN,mesh->np,mesh->ne,mesh->nprism,mesh->nt,
                          mesh->nquad,mesh->na) != 1 ) exit(EXIT_FAILURE);

  /* Swap meshes in order to free and recreate parmesh structures*/
  mesh = meshIN;
  meshIN = parmesh->listgrp[0].mesh;
  parmesh->listgrp[0].mesh = mesh;
  mesh = parmesh->listgrp[0].mesh;


  /* Free parmesh structures */
  PMMG_Free_all(PMMG_ARG_start,
                PMMG_ARG_ppParMesh,&parmesh,
                PMMG_ARG_end);
  parmesh = NULL;


  /** ------------------------------ STEP  IV -------------------------- */
  /** Recreate parmesh structures and initialize the distributed mesh using
   *  meshIN as input */

  /* Recreate parmesh structures */
  PMMG_Init_parMesh(PMMG_ARG_start,
                    PMMG_ARG_ppParMesh,&parmesh,
                    PMMG_ARG_pMesh,PMMG_ARG_pMet,
                    PMMG_ARG_dim,3,PMMG_ARG_MPIComm,MPI_COMM_WORLD,
                    PMMG_ARG_end);


  if ( PMMG_Set_meshSize(parmesh,meshIN->np,meshIN->ne,meshIN->nprism,meshIN->nt,
                               meshIN->nquad,meshIN->na) != 1 ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  /* Set points, vertex by vertex */
  for( ip = 1; ip <= meshIN->np; ip++ ) {
    ppt = &meshIN->point[ip];
    if ( PMMG_Set_vertex(parmesh,ppt->c[0],ppt->c[1],ppt->c[2],
                         ppt->ref, ip) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }

  /* Set elements, tetra by tetra */
  for( ie = 1; ie <= meshIN->ne; ie++ ) {
    pt = &meshIN->tetra[ie];
    if ( PMMG_Set_tetrahedron(parmesh,pt->v[0],pt->v[1],pt->v[2],pt->v[3],
                              pt->ref, ie) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }

  /* Vertex by vertex: for each triangle, give the vertices index, the
   * reference and the position of the triangle */
  for ( ie = 0; ie <= meshIN->nt; ie++ ) {
    ptt = &meshIN->tria[ie];
    if ( PMMG_Set_triangle(parmesh,
                           ptt->v[0],ptt->v[1],ptt->v[2],
                           ptt->ref,ie) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }

  /**    Initialization of interface communicators in ParMMG.
   *     The user can choose between providing triangles (faces) interface
   *     information (through the PMMG_APIDISTRIB_faces parameter), or nodes
   *     interface information (through the PMMG_APIDISTRIB_nodes parameter).
   */

  /* Set API mode */
  if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_APImode, API_mode ) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  };

  /* Set triangles or nodes interfaces depending on API mode */
  switch( API_mode ) {

    case PMMG_APIDISTRIB_faces :
      if( !rank ) printf("\n--- API mode: Setting face communicators\n");

      /* Set the number of interfaces */
      ier = PMMG_Set_numberOfFaceCommunicators(parmesh, n_face_comm);

      /* Loop on each interface (proc pair) seen by the current rank) */
      for( icomm = 0; icomm < n_face_comm; icomm++ ) {

        /* Set nb. of entities on interface and rank of the outward proc */
        ier = PMMG_Set_ithFaceCommunicatorSize(parmesh, icomm,
                                               color_face[icomm],
                                               nitem_face_comm[icomm]);

        /* Set local and global index for each entity on the interface */
        ier = PMMG_Set_ithFaceCommunicator_faces(parmesh, icomm,
                                                 idx_face_loc[icomm],
                                                 idx_face_glob[icomm], 1 );
      }
      break;

    case PMMG_APIDISTRIB_nodes :
      if( !rank ) printf("\n--- API mode: Setting node communicators\n");

      /* Set the number of interfaces */
      ier = PMMG_Set_numberOfNodeCommunicators(parmesh, n_node_comm);

      /* Loop on each interface (proc pair) seen by the current rank) */
      for( icomm = 0; icomm < n_node_comm; icomm++ ) {

        /* Set nb. of entities on interface and rank of the outward proc */
        ier = PMMG_Set_ithNodeCommunicatorSize(parmesh, icomm,
                                               color_node[icomm],
                                               nitem_node_comm[icomm]);

        /* Set local and global index for each entity on the interface */
        ier = PMMG_Set_ithNodeCommunicator_nodes(parmesh, icomm,
                                                 idx_node_loc[icomm],
                                                 idx_node_glob[icomm], 1 );
      }
      break;
  }

//  /** save mesh and interfaces **/
//  char filemesh[48];
//  sprintf(filemesh,"mesh_in.%d.mesh",parmesh->myrank);
//  MMG3D_saveMesh(parmesh->listgrp[0].mesh,filemesh);
//
//  FILE *fid;
//  sprintf(filemesh,"parFaces_in.%d",parmesh->myrank);
//  fid = fopen(filemesh,"w");
//  fprintf(fid,"%d\n",n_face_comm);
//  for( icomm = 0; icomm < n_face_comm; icomm++ ) {
//    fprintf(fid,"\n%d\n%d\n",color_face[icomm],nitem_face_comm[icomm]);
//    for( i = 0; i < nitem_face_comm[icomm]; i++ )
//      fprintf(fid,"%d %d\n",idx_face_loc[icomm][i],idx_face_glob[icomm][i]);
//  }
//  fclose(fid);
//
//  sprintf(filemesh,"parNodes_in.%d",parmesh->myrank);
//  fid = fopen(filemesh,"w");
//  fprintf(fid,"%d\n",n_node_comm);
//  for( icomm = 0; icomm < n_node_comm; icomm++ ) {
//    fprintf(fid,"\n%d\n%d\n",color_node[icomm],nitem_node_comm[icomm]);
//    for( i = 0; i < nitem_node_comm[icomm]; i++ )
//      fprintf(fid,"%d %d\n",idx_node_loc[icomm][i],idx_node_glob[icomm][i]);
//  }
//  fclose(fid);

  /** ------------------------------ STEP V ---------------------------- */
  /** remesh step */

  /* Set number of iterations */
  if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_niter, niter ) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  };

  /* Don't remesh the surface */
  if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_nosurf, 1 ) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  };

  /* Compute output nodes and triangles global numbering */
  if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_globalNum, 1 ) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  };



  /* remeshing function */
  ierlib = PMMG_parmmglib_distributed( parmesh );

  if ( ierlib == PMMG_SUCCESS ) {

    /** If no remeshing is performed (zero remeshing iterations), check set
     * parallel interfaces against input data. */
    if( !niter ) {

      /* Check matching of input interface nodes with the set ones */
      if( !PMMG_Check_Set_NodeCommunicators(parmesh,n_node_comm,nitem_node_comm,
                                         color_node,idx_node_loc) ) {
        printf("### Wrong set node communicators!\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }

      /* Get input triangle nodes */
      faceNodes     = (int **) malloc(n_face_comm*sizeof(int *));
      for( icomm = 0; icomm < n_face_comm; icomm++ ) {
        faceNodes[icomm]     = (int *) malloc(3*nitem_face_comm[icomm]*sizeof(int));
        for( i = 0; i < nitem_face_comm[icomm]; i++ ) {
          pos = idx_face_loc[icomm][i];
          faceNodes[icomm][3*i]   = triaNodes[3*(pos-1)];
          faceNodes[icomm][3*i+1] = triaNodes[3*(pos-1)+1];
          faceNodes[icomm][3*i+2] = triaNodes[3*(pos-1)+2];
        }
      }

      /* Check matching of input interface triangles with the set ones */
      if( !PMMG_Check_Set_FaceCommunicators(parmesh,n_face_comm,nitem_face_comm,
                                            color_face,faceNodes) ) {
        printf("### Wrong set face communicators!\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    }

    /** ------------------------------ STEP  VI -------------------------- */
    /** recover parallel interfaces */

    int n_node_comm_out,n_face_comm_out;
    int *nitem_node_comm_out,*nitem_face_comm_out;
    int *color_node_out, *color_face_out;
    int **idx_face_loc_out;
    int **idx_node_loc_out;

    /* Get number of node interfaces */
    ier = PMMG_Get_numberOfNodeCommunicators(parmesh,&n_node_comm_out);

    /* Get outward proc rank and number of nodes on each interface */
    color_node_out      = (int *) malloc(n_node_comm_out*sizeof(int));
    nitem_node_comm_out = (int *) malloc(n_node_comm_out*sizeof(int));
    for( icomm = 0; icomm < n_node_comm_out; icomm++ )
      ier = PMMG_Get_ithNodeCommunicatorSize(parmesh, icomm,
                                             &color_node_out[icomm],
                                             &nitem_node_comm_out[icomm]);

    /* Get IDs of nodes on each interface */
    idx_node_loc_out = (int **) malloc(n_node_comm_out*sizeof(int *));
    for( icomm = 0; icomm < n_node_comm_out; icomm++ )
      idx_node_loc_out[icomm] = (int *) malloc(nitem_node_comm_out[icomm]*sizeof(int));
    ier = PMMG_Get_NodeCommunicator_nodes(parmesh, idx_node_loc_out);

    /* Get number of face interfaces */
    ier = PMMG_Get_numberOfFaceCommunicators(parmesh,&n_face_comm_out);

    /* Get outward proc rank and number of faces on each interface */
    color_face_out      = (int *) malloc(n_face_comm_out*sizeof(int));
    nitem_face_comm_out = (int *) malloc(n_face_comm_out*sizeof(int));
    for( icomm = 0; icomm < n_face_comm_out; icomm++ )
      ier = PMMG_Get_ithFaceCommunicatorSize(parmesh, icomm,
                                             &color_face_out[icomm],
                                             &nitem_face_comm_out[icomm]);

    /* Get IDs of triangles on each interface */
    idx_face_loc_out = (int **) malloc(n_face_comm_out*sizeof(int *));
    for( icomm = 0; icomm < n_face_comm_out; icomm++ )
      idx_face_loc_out[icomm] = (int *) malloc(nitem_face_comm_out[icomm]*sizeof(int));
    ier = PMMG_Get_FaceCommunicator_faces(parmesh, idx_face_loc_out);

/*
    for( icomm = 0; icomm < n_node_comm_out; icomm++ )
      for( i = 0; i < nitem_node_comm_out[icomm]; i++ )
        printf("OUT rank %d comm %d color %d node %d\n",parmesh->myrank,icomm,color_node_out[icomm],idx_node_loc_out[icomm][i]);

    for( icomm = 0; icomm < n_face_comm_out; icomm++ )
      for( i = 0; i < nitem_face_comm_out[icomm]; i++ )
        printf("OUT rank %d comm %d color %d tria %d\n",parmesh->myrank,icomm,color_face_out[icomm],idx_face_loc_out[icomm][i]);
*/

    /** If no remeshing is performed (zero remeshing iterations), check
     *  retrieved parallel interfaces against input data. */
    if( !niter ) {

      /* Check matching of input interface nodes with the output ones */
      if( !PMMG_Check_Get_NodeCommunicators(parmesh,
                                            n_node_comm,nitem_node_comm,
                                            color_node,idx_node_loc,
                                            n_node_comm_out,nitem_node_comm_out,
                                            color_node_out,idx_node_loc_out) ) {
        printf("### Wrong retrieved node communicators!\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }

      /* Get output triangles (as the boundary is re-generated, triangle IDs
       * have changed) */
      int** faceNodes_out = (int **) malloc(n_face_comm_out*sizeof(int *));
      for( icomm = 0; icomm < n_face_comm_out; icomm++ ) {
        faceNodes_out[icomm] = (int *) malloc(3*nitem_face_comm_out[icomm]*sizeof(int));
        for( i = 0; i < nitem_face_comm_out[icomm]; i++ ) {
          pos = idx_face_loc_out[icomm][i];
          faceNodes_out[icomm][3*i]   = triaNodes[3*(pos-1)];
          faceNodes_out[icomm][3*i+1] = triaNodes[3*(pos-1)+1];
          faceNodes_out[icomm][3*i+2] = triaNodes[3*(pos-1)+2];
        }
      }

      /* Check matching of input interface triangles with the output ones */
      if( !PMMG_Check_Get_FaceCommunicators(parmesh,
                                            n_face_comm,nitem_face_comm,
                                            color_face,faceNodes,
                                            n_face_comm_out,nitem_face_comm_out,
                                            color_face_out,faceNodes_out) ) {
        printf("### Wrong retrieved face communicators!\n");
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }

      for( icomm = 0; icomm < n_face_comm; icomm++ )
        free(faceNodes[icomm]);
      free(faceNodes);
      for( icomm = 0; icomm < n_face_comm_out; icomm++ )
        free(faceNodes_out[icomm]);
      free(faceNodes_out);
    }

    free(ref);
    free(required);
    free(triaNodes);
    free(nitem_node_comm);
    free(nitem_face_comm);
    free(color_node);
    free(color_face);
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

    free(nitem_node_comm_out);
    free(nitem_face_comm_out);
    free(color_node_out);
    free(color_face_out);
    for( icomm = 0; icomm < n_node_comm_out; icomm++ )
      free(idx_node_loc_out[icomm]);
    free(idx_node_loc_out);
    for( icomm = 0; icomm < n_face_comm_out; icomm++ )
      free(idx_face_loc_out[icomm]);
    free(idx_face_loc_out);


    /** ------------------------------ STEP VII --------------------------- */
    /** get results */
    /** Two solutions: just use the PMMG_saveMesh/PMMG_saveSol functions
        that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
        using the PMMG_getMesh/PMMG_getSol functions */

    /** 1) Get the mesh with ParMmg getters and save it at the Medit file format */
    if( !(inm = fopen(fileout,"w")) ) {
      fprintf(stderr,"  ** UNABLE TO OPEN OUTPUT MESH FILE.\n");
      exit(EXIT_FAILURE);
    }
    fprintf(inm,"MeshVersionFormatted 2\n");
    fprintf(inm,"\nDimension 3\n");

    /** a) get the size of the mesh: vertices, tetra, triangles, edges and
     * allocate the arrays to receive data */
    nVertices   = 0;
    nTetrahedra = 0;
    nTriangles  = 0;
    nEdges      = 0;
    if ( PMMG_Get_meshSize(parmesh,&nVertices,&nTetrahedra,NULL,&nTriangles,NULL,
                           &nEdges) !=1 ) {
      ier = PMMG_STRONGFAILURE;
    }

    /* Table to store the vertices */
    double *vert = (double*)calloc((nVertices)*3,sizeof(double));
    if ( !vert ) {
      perror("  ## Memory problem: point calloc");
      nVertices = 0;
      ier = PMMG_STRONGFAILURE;
    }

    /* Table to store the tetra */
    int *tetra = (int*)calloc((nTetrahedra)*4,sizeof(int));
    if ( !tetra ) {
      perror("  ## Memory problem: tetra calloc");
      nTetrahedra = 0;
      ier = PMMG_STRONGFAILURE;
    }

    /* Table to store the tria */
    int *tria = (int*)calloc((nTriangles)*3,sizeof(int));
    if ( !tria ) {
      perror("  ## Memory problem: tria calloc");
      nTriangles = 0;
      ier = PMMG_STRONGFAILURE;
    }

    /* Table to store the edges */
    int *edge = (int*)calloc((nEdges)*2,sizeof(int));
    if ( !edge ) {
      perror("  ## Memory problem: edge calloc");
      nEdges = 0;
      ier = PMMG_STRONGFAILURE;
    }

    /* Table to store the vertices/tetra/triangles/edges references */
    int *ref = (int*)calloc(MAX4(nVertices,nTetrahedra,nTriangles,nEdges),sizeof(int));
    if ( !ref ) {
      perror("  ## Memory problem: ref calloc");
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    /* Table to know if a vertex is corner */
    int *corner = (int*)calloc(nVertices,sizeof(int));
    if ( !corner ) {
      perror("  ## Memory problem: corner calloc");
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    /* Table to know if a vertex/tetra/tria/edge is required */
    int *required = (int*)calloc(MAX4(nVertices,nTetrahedra,nTriangles,nEdges),sizeof(int));
    if ( !required ) {
      perror("  ## Memory problem: required calloc");
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    /* Table to know if an edge delimits a sharp angle */
    int *ridge = (int*)calloc(nEdges ,sizeof(int));
    if ( !ridge ) {
      perror("  ## Memory problem: ridge calloc");
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    /* Table for global node numbering */
    int *nodeGloNum  = (int *) malloc((2*nVertices)*sizeof(int));
    if ( !nodeGloNum ) {
      perror("  ## Memory problem: global nodes numbering calloc");
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }


    /** a) Nodes global numbering recovering */
    if( !opt ) {
      /* By array */
      if( !PMMG_Get_verticesGloNum( parmesh,nodeGloNum,&nodeGloNum[nVertices]) ) {
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    } else {
      /* Vertex by vertex */
      for( k = 1; k <= nVertices; k++ ) {
        if( !PMMG_Get_vertexGloNum( parmesh, &nodeGloNumber, &nodeOwner ) ) {
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
      }
    }


    /** b) Vertex recovering */
    nreq = nc = 0;
    fprintf(inm,"\nVertices\n%d\n",nVertices);

    if ( !opt ) {
      /* By array */
      if ( PMMG_Get_vertices(parmesh,vert,ref,corner,required) != 1 ) {
        fprintf(inm,"Unable to get mesh vertices \n");
        ier = PMMG_STRONGFAILURE;
      }
      for ( k=0; k<nVertices; k++ ) {
        if ( corner && corner[k] )  nc++;
        if ( required && required[k] )  nreq++;
      }
    }
    else {
      /* Vertex by vertex */
      for ( k=0; k<nVertices; k++ ) {
        pos = 3*k;
        if ( PMMG_Get_vertex(parmesh,&(vert[pos]),&(vert[pos+1]),&(vert[pos+2]),
                             &(ref[k]),&(corner[k]),&(required[k])) != 1 ) {
          fprintf(inm,"Unable to get mesh vertex %d \n",k);
          ier = PMMG_STRONGFAILURE;
        }
        if ( corner && corner[k] )  nc++;
        if ( required && required[k] )  nreq++;
      }
    }
    for ( k=0; k<nVertices; k++ ) {
      pos = 3*k;
      fprintf(inm,"%.15lg %.15lg %.15lg %d \n",vert[pos],vert[pos+1],vert[pos+2],ref[k]);
    }

    fprintf(inm,"\nCorners\n%d\n",nc);
    for ( k=0; k<nVertices; k++ ) {
      if ( corner && corner[k] )  fprintf(inm,"%d \n",k);
    }
    fprintf(inm,"\nRequiredVertices\n%d\n",nreq);
    for ( k=0; k<nVertices; k++ ) {
      if ( required && required[k] )  fprintf(inm,"%d \n",k);
    }
    free(corner);
    corner = NULL;

    /** d) Triangles recovering */
    nreq = 0;
    fprintf(inm,"\nTriangles\n%d\n",nTriangles);

    if ( !opt ) {
      /* By array */
      if ( PMMG_Get_triangles(parmesh,tria,ref,required) != 1 ) {
        fprintf(inm,"Unable to get mesh triangles\n");
        ier = PMMG_STRONGFAILURE;
      }
      for ( k=0; k<nTriangles; k++ ) {
        if ( required && required[k] )  nreq++;
      }
    }
    else {
      /* Triangle by triangle */
      for ( k=0; k<nTriangles; k++ ) {
        pos = 3*k;
        if ( PMMG_Get_triangle(parmesh,&(tria[pos]),&(tria[pos+1]),&(tria[pos+2]),
                               &(ref[k]),&(required[k])) != 1 ) {
          fprintf(inm,"Unable to get mesh triangle %d \n",k);
          ier = PMMG_STRONGFAILURE;
        }
        if ( required && required[k] )  nreq++;
      }
    }
    for ( k=0; k<nTriangles; k++ ) {
      pos = 3*k;
      fprintf(inm,"%d %d %d %d \n",tria[pos],tria[pos+1],tria[pos+2],ref[k]);
    }


    fprintf(inm,"\nRequiredTriangles\n%d\n",nreq);
    for ( k=0; k<nTriangles; k++ ) {
      if ( required && required[k] )  fprintf(inm,"%d \n",k);
    }

    /** e) Edges recovering */
    nreq = 0;nr = 0;
    fprintf(inm,"\nEdges\n%d\n",nEdges);

    if ( !opt ) {
      /* By array */
      if ( PMMG_Get_edges(parmesh,edge,ref,ridge,required) != 1 ) {
        fprintf(inm,"Unable to get mesh edges\n");
        ier = PMMG_STRONGFAILURE;
      }
      for ( k=0; k<nEdges; k++ ) {
        if ( ridge && ridge[k] )  nr++;
        if ( required && required[k] )  nreq++;
      }
    }
    else {
      /* Edge by edge */
      for ( k=0; k<nEdges; k++ ) {
        pos = 2*k;
        if ( PMMG_Get_edge(parmesh,&(edge[pos]),&(edge[pos+1]),
                           &(ref[k]),&(ridge[k]),&(required[k])) != 1 ) {
          fprintf(inm,"Unable to get mesh edge %d \n",k);
          ier = PMMG_STRONGFAILURE;
        }
        if ( ridge && ridge[k] )  nr++;
        if ( required && required[k] )  nreq++;
      }
    }
    for ( k=0; k<nEdges; k++ ) {
      pos = 2*k;
      fprintf(inm,"%d %d %d \n",edge[pos],edge[pos+1],ref[k]);
    }

    fprintf(inm,"\nRequiredEdges\n%d\n",nreq);
    for ( k=0; k<nEdges; k++ ) {
      if ( required && required[k] )  fprintf(inm,"%d \n",k);
    }
    fprintf(inm,"\nRidges\n%d\n",nr);
    for ( k=0; k<nEdges; k++ ) {
      if ( ridge && ridge[k] )  fprintf(inm,"%d \n",k);
    }

    /** c) Tetra recovering */
    nreq = 0;
    fprintf(inm,"\nTetrahedra\n%d\n",nTetrahedra);

    if ( !opt ) {
      /* By array */
      if ( PMMG_Get_tetrahedra(parmesh,tetra,ref,required) != 1 ) {
        fprintf(inm,"Unable to get mesh tetra\n");
        ier = PMMG_STRONGFAILURE;
      }
      for ( k=0; k<nTetrahedra; k++ ) {
        if ( required && required[k] )  nreq++;
      }
    }
    else {
      /* Tetra by tetra */
      for ( k=0; k<nTetrahedra; k++ ) {
        pos = 4*k;
        if ( PMMG_Get_tetrahedron(parmesh,
                                  &(tetra[pos  ]),&(tetra[pos+1]),
                                  &(tetra[pos+2]),&(tetra[pos+3]),
                                  &(ref[k]),&(required[k])) != 1 ) {
          fprintf(inm,"Unable to get mesh tetra %d \n",k);
          ier = PMMG_STRONGFAILURE;
        }
        if ( required && required[k] )  nreq++;
      }
    }
    for ( k=0; k<nTetrahedra; k++ ) {
      pos = 4*k;
      fprintf(inm,"%d %d %d %d %d \n",
              tetra[pos],tetra[pos+1],tetra[pos+2],tetra[pos+3],ref[k]);
    }

    fprintf(inm,"\nRequiredTetrahedra\n%d\n",nreq);
    for ( k=0; k<nTetrahedra; k++ ) {
      if ( required && required[k] )  fprintf(inm,"%d \n",k);
    }

    fprintf(inm,"\nEnd\n");
    fclose(inm);

    free(vert)       ; vert     = NULL;
    free(tetra)      ; tetra    = NULL;
    free(tria)       ; tria     = NULL;
    free(edge)       ; edge     = NULL;
    free(ref)        ; ref      = NULL;
    free(required)   ; required = NULL;
    free(ridge)      ; ridge    = NULL;
    free(nodeGloNum) ; nodeGloNum = NULL;

    /** 3) Get the metric with ParMmg getters */
    if ( !(inm = fopen(metout,"w")) ) {
      fprintf(stderr,"  ** UNABLE TO OPEN OUTPUT METRIC FILE.\n");
      exit(EXIT_FAILURE);
    }
    fprintf(inm,"MeshVersionFormatted 2\n");
    fprintf(inm,"\nDimension 3\n");

    /** a) get the size of the metric: type of entity to which apply the
        metric(SolAtVertices,...), number of entities to which apply the metric,
        type of solution (scalar, tensor...) */
    nVertices = 0;
    int typEntity,typSol;
    if ( PMMG_Get_metSize(parmesh,&typEntity,&nVertices,&typSol) != 1 ) {
      printf("Unagle to get metric size\n");
      nVertices = 0;
      ier = PMMG_LOWFAILURE;
    }

    /* We set a scalar metric so the output metric must be scalar */
    if ( ( typEntity != MMG5_Vertex )  || ( typSol != MMG5_Scalar ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    /** b) Vertex recovering */
    double *sol = (double*)calloc(nVertices+1 ,sizeof(double));

    fprintf(inm,"\nSolAtVertices\n%d\n",nVertices);
    fprintf(inm,"1 1 \n\n");
    if ( !opt ) {
      /* by array */
      if ( PMMG_Get_scalarMets(parmesh,sol) != 1 ) {
        fprintf(inm,"Unable to get metrics\n");
        ier = PMMG_LOWFAILURE;
      }
    }
    else {
      for ( k=0; k<nVertices; k++ ) {
        /* Vertex by vertex */
        if ( PMMG_Get_scalarMet(parmesh,&sol[k]) != 1 ) {
          fprintf(inm,"Unable to get metrics %d \n",k);
          ier = PMMG_LOWFAILURE;
        }
      }
    }
    for ( k=0; k<nVertices; k++ ) {
      fprintf(inm,"%.15lg \n",sol[k]);
    }

    fprintf(inm,"\nEnd\n");
    fclose(inm);

    free(sol);

    /** 4) Get the solutions with ParMmg getters */
    // To implement when ParMmg will be ready
  }
  else if ( ierlib == PMMG_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF PARMMGLIB: UNABLE TO SAVE MESH\n");
  }


  /* Free auxiliary mesh structures */
  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&meshIN,
                 MMG5_ARG_ppMet,&solIN,
                 MMG5_ARG_end);


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

  return ierlib;
}

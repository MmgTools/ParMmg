/**
 * Example of use of the parmmg library with a distributed input mesh (basic
 * use of mesh adaptation)
 *
 * This example show how to set and get parallel mesh interfaces using the
 * ParMmg setters and getters, starting from a global mesh which is
 * automatically partitioned and distributed among the processes.
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

/** Include the parmmg and mmg3d library header file */
#include "libparmmg.h"
#include "libmmg3d.h"

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

int main(int argc,char *argv[]) {
  PMMG_pParMesh   parmesh;
  int             ier,ierlib,rank,k;
  int             opt,API_mode,niter;
  char            *filename,*metname,*solname,*fileout,*metout,*tmp;
  FILE            *inm;
  int             pos,nreq,nc,nr;
  int             nVertices,nTetrahedra,nTriangles,nEdges;


  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if ( !rank ) fprintf(stdout,"  -- TEST PARMMGLIB: show manual recovering of parallel mesh\n");

  solname = NULL;
  metname = NULL;
  tmp     = NULL;

  if ( (argc!=4) && !rank ) {
    printf(" Usage: %s filein fileout io_option\n",argv[0]);
    printf("     API_mode = 0   to Get/Set the parallel interfaces through triangles\n");
    printf("     API_mode = 1   to Get/Set the parallel interfaces through nodes\n");
    return 1;
  }

  /* Get API mode (face or node interfaces) */
  API_mode = atoi(argv[3]);


  /* Name and path of the mesh file */
  filename = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
  if ( filename == NULL ) {
    perror("  ## Memory problem: calloc");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  strcpy(filename,argv[1]);

  fileout = (char *) calloc(strlen(argv[2]) + 8, sizeof(char));
  if ( fileout == NULL ) {
    perror("  ## Memory problem: calloc");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  strcpy(fileout,argv[2]);


  metout = (char *) calloc(strlen(argv[2]) + 1, sizeof(char));
  if ( metout == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(metout,argv[2]);

  /* Option to Set mesh entities vertex by vertex */
  opt      = 1;

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

  /* Load mesh and communicators */
  if ( !PMMG_loadMesh_distributed(parmesh,filename) ) {
    fprintf ( stderr, "Error: Unable to load %s distributed mesh.\n",filename);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }


  /** ------------------------------ STEP II ---------------------------- */
  /** remesh step */

  /* Set the number of remeshing iterations */
  niter = 1;
  if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_niter, niter ) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  };

  /* Remesh the surface */
  if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_nosurf, 0 ) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  };

  if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_angle, 45 ) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  };

  if( !PMMG_Set_dparameter( parmesh, PMMG_DPARAM_hsiz, 1.0 ) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  };

  if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_verbose, 6 ) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  };


  /* remeshing function */
  ierlib = PMMG_parmmglib_distributed( parmesh );

  if ( ierlib == PMMG_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF PARMMGLIB: UNABLE TO SAVE MESH\n");

    /* Free the PMMG5 structures */
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
  }


  /** ------------------------------ STEP III --------------------------- */
  /** get results */
  /** Two solutions: just use the
      PMMG_saveMesh_distributed/PMMG_saveMet_distributed functions that will
      write .mesh(b)/.sol formatted files or manually get your mesh/sol using
      the PMMG_getMesh/PMMG_getSol functions */

  /** 1) Get the mesh with ParMmg getters and save it at the Medit file format */
  sprintf(fileout,"%s_%d.mesh",fileout,rank);
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

  /** b) get local mesh */

  /* allocate the arrays to receive data */
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

  /** Vertex recovering */
  nreq = nc = 0;
  fprintf(inm,"\nVertices\n%d\n",nVertices);

  /* By array */
  if ( PMMG_Get_vertices(parmesh,vert,ref,corner,required) != 1 ) {
    fprintf(inm,"Unable to get mesh vertices \n");
    ier = PMMG_STRONGFAILURE;
  }
  for ( k=0; k<nVertices; k++ ) {
    if ( corner && corner[k] )  nc++;
    if ( required && required[k] )  nreq++;
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

  /** Triangles recovering */
  nreq = 0;
  fprintf(inm,"\nTriangles\n%d\n",nTriangles);

  /* By array */
  if ( PMMG_Get_triangles(parmesh,tria,ref,required) != 1 ) {
    fprintf(inm,"Unable to get mesh triangles\n");
    ier = PMMG_STRONGFAILURE;
  }
  for ( k=0; k<nTriangles; k++ ) {
    if ( required && required[k] )  nreq++;
  }

  for ( k=0; k<nTriangles; k++ ) {
    pos = 3*k;
    fprintf(inm,"%d %d %d %d \n",tria[pos],tria[pos+1],tria[pos+2],ref[k]);
  }


  fprintf(inm,"\nRequiredTriangles\n%d\n",nreq);
  for ( k=0; k<nTriangles; k++ ) {
    if ( required && required[k] )  fprintf(inm,"%d \n",k);
  }

  /** Edges recovering */
  nreq = 0;nr = 0;
  fprintf(inm,"\nEdges\n%d\n",nEdges);

  /* By array */
  if ( PMMG_Get_edges(parmesh,edge,ref,ridge,required) != 1 ) {
    fprintf(inm,"Unable to get mesh edges\n");
    ier = PMMG_STRONGFAILURE;
  }
  for ( k=0; k<nEdges; k++ ) {
    if ( ridge && ridge[k] )  nr++;
    if ( required && required[k] )  nreq++;
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

  /** Tetra recovering */
  nreq = 0;
  fprintf(inm,"\nTetrahedra\n%d\n",nTetrahedra);

  /* By array */
  if ( PMMG_Get_tetrahedra(parmesh,tetra,ref,required) != 1 ) {
    fprintf(inm,"Unable to get mesh tetra\n");
    ier = PMMG_STRONGFAILURE;
  }
  for ( k=0; k<nTetrahedra; k++ ) {
    if ( required && required[k] )  nreq++;
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
  free(vert)    ; vert     = NULL;
  free(tetra)   ; tetra    = NULL;
  free(tria)    ; tria     = NULL;
  free(edge)    ; edge     = NULL;
  free(ref)     ; ref      = NULL;
  free(required); required = NULL;
  free(ridge)   ; ridge    = NULL;

  /** c) parallel information recovering (communicators) */
  /** recover parallel interfaces (depending on choosed mode (API_mode)) */
  if ( API_mode == 0 ) {
    int icomm,i;
    int n_face_comm_out;
    int *nitem_face_comm_out;
    int *color_face_out;
    int **idx_face_loc_out,**idx_face_glob_out;

    /* Get number of face interfaces */
    ier = PMMG_Get_numberOfFaceCommunicators(parmesh,&n_face_comm_out);
    fprintf(inm,"\nParallelTriangleCommunicators\n%d\n",n_face_comm_out);

    /* Get outward proc rank and number of faces on each interface */
    color_face_out      = (int *) malloc(n_face_comm_out*sizeof(int));
    nitem_face_comm_out = (int *) malloc(n_face_comm_out*sizeof(int));
    for( icomm = 0; icomm < n_face_comm_out; icomm++ )
      ier = PMMG_Get_ithFaceCommunicatorSize(parmesh, icomm,
                                             &color_face_out[icomm],
                                             &nitem_face_comm_out[icomm]);

    for( icomm = 0; icomm < n_face_comm_out; icomm++ ) {
      fprintf(inm,"%d %d\n",color_face_out[icomm],nitem_face_comm_out[icomm]);
    }

    /* Get  of triangles on each interface */
    idx_face_loc_out  = (int **) malloc(n_face_comm_out*sizeof(int *));
    idx_face_glob_out = (int **) malloc(n_face_comm_out*sizeof(int *));

    for( icomm = 0; icomm < n_face_comm_out; icomm++ ) {
      idx_face_loc_out[icomm]  = (int *) malloc(nitem_face_comm_out[icomm]*sizeof(int));
      idx_face_glob_out[icomm] = (int *) malloc(nitem_face_comm_out[icomm]*sizeof(int));
    }

    /* Get local index of interface nodes */
    ier = PMMG_Get_FaceCommunicator_faces(parmesh, idx_face_loc_out);

    /* If needed, get: the entity owner (process id), the entity global index,
     * the number of interface tria on rank, the total number of interface
     * triangles. Here we only get the global index of each interface node.
     */
    ier = PMMG_Get_FaceCommunicator_owners(parmesh,NULL,idx_face_glob_out,NULL,NULL);

    fprintf(inm,"\nParallelCommunicatorFaces\n");
    for( icomm = 0; icomm < n_face_comm_out; icomm++ ) {
      for ( i=0; i<nitem_face_comm_out[icomm];++i ) {
        fprintf(inm,"%d %d %d\n",idx_face_loc_out[icomm][i],idx_face_glob_out[icomm][i],icomm);
      }
    }

    free(nitem_face_comm_out);
    free(color_face_out);
    for( icomm = 0; icomm < n_face_comm_out; icomm++ ) {
      free(idx_face_loc_out[icomm]);
      free(idx_face_glob_out[icomm]);
    }

    free(idx_face_loc_out);
    free(idx_face_glob_out);
  }
  else {
    int icomm,i;
    int n_node_comm_out;
    int *nitem_node_comm_out;
    int *color_node_out;
    int **idx_node_loc_out,**idx_node_glob_out;

    /* Get number of node interfaces */
    ier = PMMG_Get_numberOfNodeCommunicators(parmesh,&n_node_comm_out);
    fprintf(inm,"\nParallelVertexCommunicators\n%d\n",n_node_comm_out);

    /* Get outward proc rank and number of nodes on each interface */
    color_node_out      = (int *) malloc(n_node_comm_out*sizeof(int));
    nitem_node_comm_out = (int *) malloc(n_node_comm_out*sizeof(int));
    for( icomm = 0; icomm < n_node_comm_out; icomm++ )
      ier = PMMG_Get_ithNodeCommunicatorSize(parmesh, icomm,
                                             &color_node_out[icomm],
                                             &nitem_node_comm_out[icomm]);

    for( icomm = 0; icomm < n_node_comm_out; icomm++ ) {
      fprintf(inm,"%d %d\n",color_node_out[icomm],nitem_node_comm_out[icomm]);
    }

    /* Get IDs of nodes on each interface */
    idx_node_loc_out  = (int **) malloc(n_node_comm_out*sizeof(int *));
    idx_node_glob_out = (int **) malloc(n_node_comm_out*sizeof(int *));

    for( icomm = 0; icomm < n_node_comm_out; icomm++ ) {
      idx_node_loc_out[icomm]  = (int *) malloc(nitem_node_comm_out[icomm]*sizeof(int));
      idx_node_glob_out[icomm] = (int *) malloc(nitem_node_comm_out[icomm]*sizeof(int));
    }

    /* Get local index of interface nodes */
    ier = PMMG_Get_NodeCommunicator_nodes(parmesh, idx_node_loc_out);

    /* If needed, get: the entity owner (process id), the entity global index,
     * the number of interface nodes on rank, the total number of interface
     * nodes. Here we only get the global index of each interface node.
     */
    ier = PMMG_Get_NodeCommunicator_owners(parmesh,NULL,idx_node_glob_out,NULL,NULL);

    fprintf(inm,"\nParallelCommunicatorVertices\n");
    for( icomm = 0; icomm < n_node_comm_out; icomm++ ) {
      for ( i=0; i<nitem_node_comm_out[icomm];++i ) {
        fprintf(inm,"%d %d %d\n",idx_node_loc_out[icomm][i],idx_node_glob_out[icomm][i],icomm);
      }
    }

    free(nitem_node_comm_out);
    free(color_node_out);
    for( icomm = 0; icomm < n_node_comm_out; icomm++ ) {
      free(idx_node_glob_out[icomm]);
      free(idx_node_loc_out[icomm]);
    }
    free(idx_node_glob_out);
    free(idx_node_loc_out);
  }

  fprintf(inm,"\nEnd\n");
  fclose(inm);


  /** d) Get the metric with ParMmg getters */
  if ( !(inm = fopen(metout,"w")) ) {
    fprintf(stderr,"  ** UNABLE TO OPEN OUTPUT METRIC FILE.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(inm,"MeshVersionFormatted 2\n");
  fprintf(inm,"\nDimension 3\n");

  /** get the size of the metric: type of entity to which apply the
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

  /** Solution recovering */
  double *sol = (double*)calloc(nVertices+1 ,sizeof(double));

  fprintf(inm,"\nSolAtVertices\n%d\n",nVertices);
  fprintf(inm,"1 1 \n\n");

  /* by array */
  if ( PMMG_Get_scalarMets(parmesh,sol) != 1 ) {
    fprintf(inm,"Unable to get metrics\n");
    ier = PMMG_LOWFAILURE;
  }
  for ( k=0; k<nVertices; k++ ) {
    fprintf(inm,"%.15lg \n",sol[k]);
  }
  fprintf(inm,"\nEnd\n");
  fclose(inm);

  free(sol);

  /** ------------------------------ STEP  VI -------------------------- */

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

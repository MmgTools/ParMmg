/**
 * Example of use of the parmmg library (basic use of mesh adaptation).
 *
 * This example show how to set and get the mesh, metric and solutions using the
 * ParMmg setters and getters. Depending of the command line option "opt", the
 * values are setted entities by entities or the entire value array is directly
 * setted.
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

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))


int main(int argc,char *argv[]) {
  PMMG_pParMesh   parmesh;
  int             ier,rank,i,nsols,k,opt;
  char            *filename,*metname,*solname,*fileout,*metout,*solout,*tmp;
  FILE            *inm;
  int             pos,nreq,nc,nr;


  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if ( !rank ) fprintf(stdout,"  -- TEST PARMMGLIB \n");

  solname = NULL;
  metname = NULL;
  solout  = NULL;
  metout  = NULL;
  tmp     = NULL;

  if ( (argc!=3) && !rank ) {
    printf(" Usage: %s fileout io_option\n",argv[0]);
    printf("     io_option = 0 to Get/Set the mesh/metric/solution by array\n");
    printf("     io_option = 1 to Get/Set the mesh/metric/solution vertex by vertex\n");
    return 1;
  }

  fileout = (char *) calloc(strlen(argv[1]) + 6, sizeof(char));
  if ( fileout == NULL ) {
    perror("  ## Memory problem: calloc");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  strcpy(fileout,argv[1]);
  strcat(fileout,".mesh");

  metout = (char *) calloc(strlen(argv[1]) + 9, sizeof(char));
  if ( metout == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(metout,argv[1]);
  strcat(metout,"-met.sol");

  solout = (char *) calloc(strlen(argv[1]) + 13, sizeof(char));
  if ( solout == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(solout,argv[1]);
  strcat(solout,"-solphys.sol");

  opt = atoi(argv[2]);

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

  /** with ParMmg setters */
  /** a) give the size of the mesh */
  int nVertices       = 12;
  int nTetrahedra     = 12;
  int nPrisms         = 0;
  int nTriangles      = 20;
  int nQuadrilaterals = 0;
  int nEdges          = 0;

  if(parmesh->myrank == parmesh->info.root){
    if ( PMMG_Set_meshSize(parmesh,nVertices,nTetrahedra,nPrisms,nTriangles,
                              nQuadrilaterals,nEdges) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }

  /** b) give the vertices (12 vertices with 3 coor = array of size 36) */
  double vert_coor[36] = { 0.0, 0.0, 0.0,
                           0.5, 0.0, 0.0,
                           0.5, 0.0, 1.0,
                           0.0, 0.0, 1.0,
                           0.0, 1.0, 0.0,
                           0.5, 1.0, 0.0,
                           0.5, 1.0, 1.0,
                           0.0, 1.0, 1.0,
                           1.0, 0.0, 0.0,
                           1.0, 1.0, 0.0,
                           1.0, 0.0, 1.0,
                           1.0, 1.0, 1.0  };

  int  vert_ref[12] = {0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  };

  if(parmesh->myrank == parmesh->info.root){
    if ( !opt ) {
      /* By array: give the array of the vertices coordinates such as the
       * coordinates of the k^th point are stored in vert_coor[3*(k-1)],
       * vert_coor[3*(k-1)+1] and vert_coor[3*(k-1)+2] */
      if ( PMMG_Set_vertices(parmesh,vert_coor,vert_ref) != 1 ) {
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    }
    else {
      /* Vertex by vertex: for each vertex, give the coordinates, the reference
         and the position in mesh of the vertex */
      for ( k=0; k<nVertices; ++k ) {
        pos = 3*k;
        if ( PMMG_Set_vertex(parmesh,vert_coor[pos],vert_coor[pos+1],vert_coor[pos+2],
                             vert_ref[k], k+1) != 1 ) {
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  /** c) give the tetrahedras (12 tetra with 4 vertices = array of size 48) */
  int tetra_vert[48] = { 1,  4,  2,  8,
                         8,  3,  2,  7,
                         5,  2,  6,  8,
                         5,  8,  1,  2,
                         7,  2,  8,  6,
                         2,  4,  3,  8,
                         9,  2,  3,  7,
                         7, 11,  9, 12,
                         6,  9, 10,  7,
                         6,  7,  2,  9,
                         12, 9,  7, 10,
                         9,  3, 11,  7  };

  int tetra_ref[12] = {1  ,1  ,1  ,1  ,1  ,1  ,2  ,2  ,2  ,2  ,2  ,2  };

  if(parmesh->myrank == parmesh->info.root){
    if ( !opt ) {
      /* By array: give the array of the tetra vertices and the array of the tetra
       * references. The array of the tetra vertices is such as the four
       * vertices of the k^th tetra are respectively stored in
       * tetra_vert[4*(k-1)],tetra_vert[4*(k-1)+1],tetra_vert[4*(k-1)+2] and
       * tetra_vert[4*(k-1)+3]. */
      if ( PMMG_Set_tetrahedra(parmesh,tetra_vert,tetra_ref) != 1 ) {
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    }
    else {
      /* Vertex by vertex: for each tetrahedra,
        give the vertices index, the reference and the position of the tetra */
      for ( k=0; k<nTetrahedra; ++k ) {
        pos = 4*k;
  
        if ( PMMG_Set_tetrahedron(parmesh,tetra_vert[pos],tetra_vert[pos+1],
                                  tetra_vert[pos+2],tetra_vert[pos+3],tetra_ref[k],k+1) != 1 ) {
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  /** d) give the triangles (20 tria with 3 vertices = array of size 60  */
  int tria_vert[60] = { 1,  4,  8,
                        1,  2,  4,
                        8,  3,  7,
                        5,  8,  6,
                        5,  6,  2,
                        5,  2,  1,
                        5,  1,  8,
                        7,  6,  8,
                        4,  3,  8,
                        2,  3,  4,
                        9,  3,  2,
                        11, 9, 12,
                        7, 11, 12,
                        6,  7, 10,
                        6, 10,  9,
                        6,  9,  2,
                        12,10,  7,
                        12, 9, 10,
                        3, 11,  7,
                        9, 11,  3  };

  int tria_ref[20] = { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                       4, 4, 4, 4, 4, 4, 4, 4, 4, 4  };

  if(parmesh->myrank == parmesh->info.root){
    if ( !opt ) {
      /* By array: give the array of the tria vertices and the array of the tria
       * references. The array of the tria vertices is such as the three
       * vertices of the k^th tria are stored in
       * tria_vert[3*(k-1)], tria_vert[3*(k-1)+1] and tria_vert[4*(k-1)+2] */
      if ( PMMG_Set_triangles(parmesh,tria_vert,tria_ref) != 1 ) {
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    }
    else {
      /* Vertex by vertex: for each triangle, give the vertices index, the
       * reference and the position of the triangle */
      for ( k=0; k<nTriangles; ++k ) {
        pos = 3*k;
        if ( PMMG_Set_triangle(parmesh,
                               tria_vert[pos],tria_vert[pos+1],tria_vert[pos+2],
                               tria_ref[k],k+1) != 1 ) {
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  /** 3) Build metric in ParMmg format */
  /** Two solutions: just use the PMMG_loadMet_centralized function that will
      read a .sol(b) file formatted or manually set your sol using the
      PMMG_Set* functions */

  /** Manually set of the metric */
  /** a) give info for the metric structure: metric applied on vertex entities,
      number of vertices, the metric is scalar*/
  if(parmesh->myrank == parmesh->info.root){
    if ( PMMG_Set_metSize(parmesh,MMG5_Vertex,nVertices,MMG5_Scalar) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }

  /** b) give solutions values and positions */
  double met[12] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

  if(parmesh->myrank == parmesh->info.root){
    if ( !opt ) {
      /* by array */
      if ( PMMG_Set_scalarMets(parmesh,met) != 1 ) {
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    }
    else {
      /* vertex by vertex */
      for ( k=0; k<nVertices ; k++ ) {
        if ( PMMG_Set_scalarMet(parmesh,met[k],k+1) != 1 ) {
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  /** 4) Build solutions in PMMG format */
  /** Two solutions: just use the PMMG_loadAllSols_centralized function that
      will read a .sol(b) file formatted or manually set your solutions using
      the PMMG_Set* functions */

  /** With parmmg setters: 3 sols per vertex, the first is scalar,
      the second vectorial, the third tensorial  */
  const int nSolsAtVertices = 3; // 3 solutions per vertex
  int solType[3] = {MMG5_Scalar,MMG5_Vector,MMG5_Tensor};

  if(parmesh->myrank == parmesh->info.root){
    if ( PMMG_Set_solsAtVerticesSize(parmesh,nSolsAtVertices,nVertices,solType) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }

  /** b) give solutions values and positions:
   - First solution (scalar) is equal to x^2 + y^2 + z^2
   - Second (vector) is (x,y,z)
   - Third (Tensor) is (100,0,0,100/(z+1),0,100/(z*z+1))
  */
  double scalar_sol[12],vector_sol[36],tensor_sol[72];
  for ( k=0; k<nVertices; k++ ) {
    pos = 3*k;

    /* First solution */
    scalar_sol[k] = vert_coor[pos]*vert_coor[pos]
      + vert_coor[pos+1]*vert_coor[pos+1]
      + vert_coor[pos+2]*vert_coor[pos+2];

    /* Second */
    vector_sol[3*k]   = vert_coor[pos];
    vector_sol[3*k+1] = vert_coor[pos+1];
    vector_sol[3*k+2] = vert_coor[pos+2];

    /* Third */
    tensor_sol[6*k]   = 100.;
    tensor_sol[6*k+1] = 0.;
    tensor_sol[6*k+2] = 0.;
    tensor_sol[6*k+3] = 100./(vert_coor[pos+2]+1.);
    tensor_sol[6*k+4] = 0.;
    tensor_sol[6*k+5] = 100./(vert_coor[pos+2]*vert_coor[pos+2]+1.);
  }

  if(parmesh->myrank == parmesh->info.root){
    if ( !opt ) {
      /* Give the solution by array */
      /* First solution */
      if ( PMMG_Set_ithSols_inSolsAtVertices(parmesh,1,scalar_sol) != 1 ) {
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
      /* Second */
      if ( PMMG_Set_ithSols_inSolsAtVertices(parmesh,2,vector_sol) != 1 ) {
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
      /* Third */
      if ( PMMG_Set_ithSols_inSolsAtVertices(parmesh,3,tensor_sol) != 1 ) {
        MPI_Finalize();
        exit(EXIT_FAILURE);
      }
    }
    else {
      /* Vertex by vertex */
      for ( k=0; k<nVertices; k++ ) {
        /* First solution */
        if ( PMMG_Set_ithSol_inSolsAtVertices(parmesh,1,&(scalar_sol[k]),k+1) != 1 ) {
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
        /* Second */
        pos = 3*k;
        if ( PMMG_Set_ithSol_inSolsAtVertices(parmesh,2,&(vector_sol[pos]),k+1) != 1 ) {
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
        /* Third */
        pos = 6*(k-1);
        if ( PMMG_Set_ithSol_inSolsAtVertices(parmesh,3,&(tensor_sol[pos]),k+1) != 1 ) {
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
      }
    }
  }

  //PMMG_saveAllSols_centralized(parmesh,"init-solphys.sol");

  /** ------------------------------ STEP  II -------------------------- */
  /** remesh function */
  ier = PMMG_parmmglib_centralized(parmesh);

  if ( ier != PMMG_STRONGFAILURE ) {

    /** ------------------------------ STEP III -------------------------- */
    /** get results */
    /** Two solutions: just use the PMMG_saveMesh/PMMG_saveSol functions
        that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
        using the PMMG_getMesh/PMMG_getSol functions */

    /** 1) Get the mesh with ParMmg getters and save it at the Medti file format */
    if ( !rank ) {
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

      free(tetra)   ; tetra    = NULL;
      free(ref)     ; ref      = NULL;
      free(required); required = NULL;
      free(ridge)   ; ridge    = NULL;

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

      /** 4) Get the solutions with ParMmg getters */
      // To implement when ParMmg will be ready
    }
  }
  else if ( ier == PMMG_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF PARMMGLIB: UNABLE TO SAVE MESH\n");
  }

  /** 5) Free the PMMG5 structures */
  PMMG_Free_all(PMMG_ARG_start,
                PMMG_ARG_ppParMesh,&parmesh,
                PMMG_ARG_end);

  free(fileout); fileout = NULL;
  free(solout) ; solout  = NULL;
  free(metout) ; metout  = NULL;

  MPI_Finalize();

  return ier;
}

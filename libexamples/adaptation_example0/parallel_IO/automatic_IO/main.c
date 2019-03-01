/**
 * Example of use of the parmmg library with a distributed input mesh (basic
 * use of mesh adaptation)
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
#include "parmmg.h"
// if the header file is in "include/parmmg"
//#include "parmmg/libparmmg.h"

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

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
int color_intfcNode(PMMG_pParMesh parmesh,int *color_out,
                    int **ifc_node_loc,int **ifc_node_glob,
                    int next_node_comm,int *nitem_node_comm) {
  MMG5_pMesh     mesh;
  MMG5_pPoint    ppt;
  MPI_Request    request;
  MPI_Status     status;
  int            npairs_loc,*npairs,*displ_pair,*glob_pair_displ;
  int            src,dst,tag,sendbuffer,recvbuffer,iproc,icomm,iloc,i,idx;

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

  /* Count nb of new pairs hosted on proc */
  npairs_loc = 0;
  for( icomm = 0; icomm < next_node_comm; icomm++ )
    if( color_out[icomm] > parmesh->myrank ) npairs_loc += nitem_node_comm[icomm];//1;

  /* Get nb of pairs and compute pair offset */
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

  /* Compute global pair enumeration (injective, non-surjective map) */
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
int color_intfcTria(PMMG_pParMesh parmesh,int *color_out,
                    int **ifc_tria_loc,int **ifc_tria_glob,
                    int next_face_comm,int *nitem_face_comm) {
  MPI_Request    request;
  MPI_Status     status;
  int            npairs_loc,*npairs,*displ_pair,*glob_pair_displ;
  int            src,dst,tag,sendbuffer,recvbuffer,iproc,icomm,iloc,i;

  PMMG_CALLOC(parmesh,npairs,parmesh->nprocs,int,"npair",return 0);
  PMMG_CALLOC(parmesh,displ_pair,parmesh->nprocs+1,int,"displ_pair",return 0);

  /* Count nb of new pairs hosted on proc */
  npairs_loc = 0;
  for( icomm = 0; icomm < next_face_comm; icomm++ )
    if( color_out[icomm] > parmesh->myrank ) npairs_loc += nitem_face_comm[icomm];//1;

  /* Get nb of pairs and compute pair offset */
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

  /* Compute global pair enumeration (injective, non-surjective map) */
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


int main(int argc,char *argv[]) {
  PMMG_pParMesh   parmesh;
  MMG5_pMesh      mesh,meshIN;
  MMG5_pSol       met,solIN;
  MMG5_pPoint     ppt;
  MMG5_pTria      ptt;
  MMG5_pTetra     pt;
  int             ip,ie,ier,ierlib,iresult,rank,i,k,nsols;
  int             opt,API_mode;
  char            *filename,*metname,*solname,*fileout,*metout,*tmp;
  FILE            *inm;
  int             pos,nreq,nc,nr;
  int             nVertices,nTetrahedra,nPrisms,nTriangles,nQuadrilaterals,nEdges;


  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if ( !rank ) fprintf(stdout,"  -- TEST PARMMGLIB \n");

  solname = NULL;
  metname = NULL;
  tmp     = NULL;

  if ( (argc!=4) && !rank ) {
    printf(" Usage: %s fileout io_option\n",argv[0]);
    printf("     API_mode = 0 to Get/Set the parallel interfaces through triangles\n");
    printf("     API_mode = 1 to Get/Set the parallel interfaces through nodes\n");
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

  opt      = 1; /* Set mesh entities vertex by vertex */
  API_mode = atoi(argv[3]);

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
  /** distribute mesh on procs */
 
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
  if ( (ier==PMMG_STRONGFAILURE) && MMG5_unscaleMesh( mesh, met ) ) {
    ier = PMMG_LOWFAILURE;
  }
  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MAX, parmesh->comm );
  if ( iresult!=PMMG_SUCCESS ) {
    return iresult;
  }

  /** Send mesh partionning to other procs */
  if ( !PMMG_distribute_mesh( parmesh ) ) {
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_LOWFAILURE);
  }


  /** ------------------------------ STEP  III ------------------------- */
  /** Swap meshes, so that you can use meshIN to initialize mesh */

  /* Create boundary entities */
  mesh = parmesh->listgrp[0].mesh;
  if( MMG3D_bdryBuild(mesh) == -1) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }


  /** recover parallel interfaces */

  int next_node_comm,next_face_comm,*nitem_node_comm,*nitem_face_comm;
  int *color_node_out, *color_face_out;
  int icomm,ncomm;

  ier = PMMG_Get_numberOfNodeCommunicators(parmesh,&next_node_comm);
  ier = PMMG_Get_numberOfFaceCommunicators(parmesh,&next_face_comm);

  color_node_out = (int *) malloc(next_node_comm*sizeof(int));
  nitem_node_comm = (int *) malloc(next_node_comm*sizeof(int));
  for( icomm=0; icomm<next_node_comm; icomm++ )
    ier = PMMG_Get_ithNodeCommunicatorSize(parmesh, icomm,
                                           &color_node_out[icomm],
                                           &nitem_node_comm[icomm]);

  color_face_out = (int *) malloc(next_face_comm*sizeof(int));
  nitem_face_comm = (int *) malloc(next_face_comm*sizeof(int));
  for( icomm=0; icomm<next_face_comm; icomm++ )
    ier = PMMG_Get_ithFaceCommunicatorSize(parmesh, icomm,
                                           &color_face_out[icomm],
                                           &nitem_face_comm[icomm]);

  int **out_tria_loc,**out_tria_glob,**out_node_loc,**out_node_glob;
 
  out_node_loc  = (int **) malloc(next_node_comm*sizeof(int *));
  out_node_glob = (int **) malloc(next_node_comm*sizeof(int *));
  for( icomm=0; icomm<next_node_comm; icomm++ ) {
    out_node_loc[icomm]  = (int *) malloc(nitem_node_comm[icomm]*sizeof(int));
    out_node_glob[icomm] = (int *) malloc(nitem_node_comm[icomm]*sizeof(int));
  }
  ier = PMMG_Get_NodeCommunicator_nodes(parmesh, out_node_loc);

  out_tria_loc  = (int **) malloc(next_face_comm*sizeof(int *));
  out_tria_glob = (int **) malloc(next_face_comm*sizeof(int *));
  for( icomm=0; icomm<next_face_comm; icomm++ ) {
    out_tria_loc[icomm]  = (int *) malloc(nitem_face_comm[icomm]*sizeof(int));
    out_tria_glob[icomm] = (int *) malloc(nitem_face_comm[icomm]*sizeof(int));
  }
  ier = PMMG_Get_FaceCommunicator_faces(parmesh, out_tria_loc);

  /* Color interface triangles with global enumeration */
  if( !color_intfcTria(parmesh,color_face_out,out_tria_loc,out_tria_glob,
                       next_face_comm,nitem_face_comm) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  /* Color interface nodes with global enumeration */
  if( !color_intfcNode(parmesh,color_node_out,out_node_loc,out_node_glob,
                       next_node_comm,nitem_node_comm) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

/*
  for( icomm=0; icomm<next_face_comm; icomm++ )
    for( i=0; i < nitem_face_comm[icomm]; i++ )
      printf("rank %d comm %d tria loc %d glob %d\n",parmesh->myrank,icomm,out_tria_loc[icomm][i],out_tria_glob[icomm][i]);


  for( icomm=0; icomm<next_node_comm; icomm++ )
    for( i=0; i < nitem_node_comm[icomm]; i++ )
      printf("rank %d comm %d node loc %d glob %d\n",parmesh->myrank,icomm,out_node_loc[icomm][i],out_node_glob[icomm][i]);
*/

  /* Create input mesh: meshIN */
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


  /* Free and recreate parmesh structures */
  PMMG_Free_all(PMMG_ARG_start,
                PMMG_ARG_ppParMesh,&parmesh,
                PMMG_ARG_end);
  parmesh = NULL;

  PMMG_Init_parMesh(PMMG_ARG_start,
                    PMMG_ARG_ppParMesh,&parmesh,
                    PMMG_ARG_pMesh,PMMG_ARG_pMet,
                    PMMG_ARG_dim,3,PMMG_ARG_MPIComm,MPI_COMM_WORLD,
                    PMMG_ARG_end);



  /** ------------------------------ STEP  IV -------------------------- */
  /** Initialize the distributed mesh using meshIN as input */

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


  /* Set triangles or nodes interfaces depending on API mode */
  switch( API_mode ) {
    
    case PMMG_APIDISTRIB_faces :
      if( !rank ) printf("\n--- API mode: Setting face communicators\n");
      
      ier = PMMG_Set_numberOfFaceCommunicators(parmesh, next_face_comm);
      for( icomm=0; icomm<next_face_comm; icomm++ ) {
        ier = PMMG_Set_ithFaceCommunicatorSize(parmesh, icomm,
                                               color_face_out[icomm],
                                               nitem_face_comm[icomm]);
        ier = PMMG_Set_ithFaceCommunicator_faces(parmesh, icomm,
                                                 out_tria_loc[icomm],
                                                 out_tria_glob[icomm], 1 );
      }
      break;
    
    case PMMG_APIDISTRIB_nodes :
      if( !rank ) printf("\n--- API mode: Setting node communicators\n");
      
      ier = PMMG_Set_numberOfNodeCommunicators(parmesh, next_node_comm);
      for( icomm=0; icomm<next_node_comm; icomm++ ) {
        ier = PMMG_Set_ithNodeCommunicatorSize(parmesh, icomm,
                                               color_node_out[icomm],
                                               nitem_node_comm[icomm]);
        ier = PMMG_Set_ithNodeCommunicator_nodes(parmesh, icomm,
                                                 out_node_loc[icomm],
                                                 out_node_glob[icomm], 1 );
      }
      break;
  }
 
  for( icomm = 0; icomm <next_face_comm; icomm++ ) {
    free(out_tria_loc[icomm]);
    free(out_tria_glob[icomm]);
  }
  free(nitem_face_comm);
  free(color_face_out);
  for( icomm = 0; icomm <next_node_comm; icomm++ ) {
    free(out_node_loc[icomm]);
    free(out_node_glob[icomm]);
  }
  free(nitem_node_comm);
  free(color_node_out);


  /** ------------------------------ STEP V ---------------------------- */
  /** remesh function */
  ierlib = PMMG_parmmglib_distributed( parmesh );

  if ( ierlib != PMMG_STRONGFAILURE ) {

    /** ------------------------------ STEP  VI -------------------------- */
    /** recover parallel interfaces */
  
    ier = PMMG_Get_numberOfNodeCommunicators(parmesh,&next_node_comm);
    ier = PMMG_Get_numberOfFaceCommunicators(parmesh,&next_face_comm);
  
    color_node_out = (int *) malloc(next_node_comm*sizeof(int));
    nitem_node_comm = (int *) malloc(next_node_comm*sizeof(int));
    for( icomm=0; icomm<next_node_comm; icomm++ )
      ier = PMMG_Get_ithNodeCommunicatorSize(parmesh, icomm,
                                             &color_node_out[icomm],
                                             &nitem_node_comm[icomm]);
 
    color_face_out = (int *) malloc(next_face_comm*sizeof(int));
    nitem_face_comm = (int *) malloc(next_face_comm*sizeof(int));
    for( icomm=0; icomm<next_face_comm; icomm++ )
      ier = PMMG_Get_ithFaceCommunicatorSize(parmesh, icomm,
                                             &color_face_out[icomm],
                                             &nitem_face_comm[icomm]);
  
   
    out_node_loc = (int **) malloc(next_node_comm*sizeof(int *));
    for( icomm=0; icomm<next_node_comm; icomm++ )
      out_node_loc[icomm] = (int *) malloc(nitem_node_comm[icomm]*sizeof(int));
    ier = PMMG_Get_NodeCommunicator_nodes(parmesh, out_node_loc);
 
    out_tria_loc = (int **) malloc(next_face_comm*sizeof(int *));
    for( icomm=0; icomm<next_face_comm; icomm++ )
      out_tria_loc[icomm] = (int *) malloc(nitem_face_comm[icomm]*sizeof(int));
    ier = PMMG_Get_FaceCommunicator_faces(parmesh, out_tria_loc);


    for( icomm=0; icomm<next_node_comm; icomm++ )
      for( i=0; i < nitem_node_comm[icomm]; i++ )
        printf("rank %d comm %d node %d\n",parmesh->myrank,icomm,out_node_loc[icomm][i]);
 
    for( icomm=0; icomm<next_face_comm; icomm++ )
      for( i=0; i < nitem_face_comm[icomm]; i++ )
        printf("rank %d comm %d tria %d\n",parmesh->myrank,icomm,out_tria_loc[icomm][i]);


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
  else if ( ierlib == PMMG_STRONGFAILURE ) {
    fprintf(stdout,"BAD ENDING OF PARMMGLIB: UNABLE TO SAVE MESH\n");
  }




  /** 5) Free the PMMG5 structures */
  PMMG_Free_all(PMMG_ARG_start,
                PMMG_ARG_ppParMesh,&parmesh,
                PMMG_ARG_end);

  free(filename);
  filename = NULL;

  free(fileout);
  fileout = NULL;

  MPI_Finalize();

  return ierlib;
}

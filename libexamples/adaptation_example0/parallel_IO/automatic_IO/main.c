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
// #include "libparmmg.h"
// if the header file is in "include/parmmg"
#include "parmmg.h"
#include "linkedlist_pmmg.h"

int color_intfcNode(PMMG_pParMesh parmesh,int *color_out,
                    int **ifc_node_loc,int **ifc_node_glob,
                    int next_node_comm,int *nitem_node_comm) {
  PMMG_lnkdList  **proclist;
  MPI_Request    request;
  MPI_Status     status;
  int            *oldIdx;
  int            npairs_loc,*npairs,*displ_pair,*glob_pair_displ;
  int            src,dst,tag,sendbuffer,recvbuffer,iproc,icomm,iloc,i;

  PMMG_CALLOC(parmesh,npairs,parmesh->nprocs,int,"npair",return 0);
  PMMG_CALLOC(parmesh,displ_pair,parmesh->nprocs+1,int,"displ_pair",return 0);

  /* Count nb of new pairs hosted on proc */
  npairs_loc = 0;
  for( icomm = 0; icomm < next_node_comm; icomm++ )
    if( color_out[icomm] > parmesh->myrank ) npairs_loc += nitem_node_comm[icomm];//1;

  /* Get nb of pairs and compute pair offset */
  MPI_Allgather( &npairs_loc,1,MPI_INT,
                 npairs,1,MPI_INT,parmesh->comm );

  for( iproc = 0; iproc < parmesh->nprocs; iproc++ )
    displ_pair[iproc+1] = displ_pair[iproc]+npairs[iproc];
  printf("rank %d npairs_loc %d displ %d\n",parmesh->myrank,npairs_loc,displ_pair[parmesh->myrank]);

  
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

  for( icomm = 0; icomm < next_node_comm; icomm++ ) {
    /* Get face comm in permuted ordered */
    printf("color_in %d color_out %d id %d\n",parmesh->myrank,color_out[icomm],glob_pair_displ[icomm]);
  }
 
  for( icomm = 0; icomm < next_node_comm; icomm++ )
    for( i=0; i < nitem_node_comm[icomm]; i++ )
      ifc_node_glob[icomm][i] = glob_pair_displ[icomm]+i;

  return 1;
}


int color_intfcTria(PMMG_pParMesh parmesh,int *color_out,
                    int **ifc_tria_loc,int **ifc_tria_glob,
                    int next_face_comm,int *nitem_face_comm) {
  PMMG_lnkdList  **proclist;
  MPI_Request    request;
  MPI_Status     status;
  int            *oldIdx;
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
  printf("rank %d npairs_loc %d displ %d\n",parmesh->myrank,npairs_loc,displ_pair[parmesh->myrank]);

  
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

  for( icomm = 0; icomm < next_face_comm; icomm++ ) {
    /* Get face comm in permuted ordered */
    printf("color_in %d color_out %d id %d\n",parmesh->myrank,color_out[icomm],glob_pair_displ[icomm]);
  }
 
  for( icomm = 0; icomm < next_face_comm; icomm++ )
    for( i=0; i < nitem_face_comm[icomm]; i++ )
      ifc_tria_glob[icomm][i] = glob_pair_displ[icomm]+i;

  return 1;
}


int main(int argc,char *argv[]) {
  PMMG_pParMesh   parmesh;
  MMG5_pMesh      mesh,meshIN;
  MMG5_pSol       met,solIN;
  MMG5_pPoint     ppt;
  MMG5_pTetra     pt;
  int             ip,ie,ier,ierlib,iresult,rank,i,nsols;
  int             opt,API_mode;
  char            *filename,*metname,*solname,*fileout,*tmp;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  if ( !rank ) fprintf(stdout,"  -- TEST PARMMGLIB \n");

  solname = NULL;
  metname = NULL;
  tmp     = NULL;

  if ( (argc!=5) && !rank ) {
    printf(" Usage: %s fileout io_option\n",argv[0]);
    printf("     io_option = 0 to Get/Set the mesh/metric/solution by array\n");
    printf("     io_option = 1 to Get/Set the mesh/metric/solution vertex by vertex\n");
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

  fileout = (char *) calloc(strlen(argv[2]) + 1, sizeof(char));
  if ( fileout == NULL ) {
    perror("  ## Memory problem: calloc");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  strcpy(fileout,argv[2]);

  opt      = atoi(argv[3]);
  API_mode = atoi(argv[4]);

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

  /* Create meshIN */
  meshIN = NULL;
  solIN = NULL;
  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&meshIN,MMG5_ARG_ppMet,&solIN,
                  MMG5_ARG_end);
 
  if ( MMG3D_Set_meshSize(meshIN,mesh->np,mesh->ne,mesh->nprism,mesh->nt,
                          mesh->nquad,mesh->na) != 1 ) exit(EXIT_FAILURE);

  /* Swap meshes */
  mesh = meshIN;
  meshIN = parmesh->listgrp[0].mesh;
  parmesh->listgrp[0].mesh = mesh;
  mesh = parmesh->listgrp[0].mesh;


  int *local_face_index,*global_face_index;
  PMMG_CALLOC(parmesh,local_face_index,meshIN->nt,int,"local tria indices",return 0);
  PMMG_CALLOC(parmesh,global_face_index,meshIN->nt,int,"global tria indices",return 0);

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

  for( icomm=0; icomm<next_node_comm; icomm++ )
    for( i=0; i < nitem_node_comm[icomm]; i++ )
      printf("rank %d comm %d node %d\n",parmesh->myrank,icomm,out_node_loc[icomm][i]);

  out_tria_loc  = (int **) malloc(next_face_comm*sizeof(int *));
  out_tria_glob = (int **) malloc(next_face_comm*sizeof(int *));
  for( icomm=0; icomm<next_face_comm; icomm++ ) {
    out_tria_loc[icomm]  = (int *) malloc(nitem_face_comm[icomm]*sizeof(int));
    out_tria_glob[icomm] = (int *) malloc(nitem_face_comm[icomm]*sizeof(int));
  }
  ier = PMMG_Get_FaceCommunicator_faces(parmesh, out_tria_loc);

  for( icomm=0; icomm<next_face_comm; icomm++ )
    for( i=0; i < nitem_face_comm[icomm]; i++ )
      printf("rank %d comm %d tria %d\n",parmesh->myrank,icomm,out_tria_loc[icomm][i]);




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

  ierlib = PMMG_SUCCESS;
//  switch( API_mode ) {
//    
//    case PMMG_APIDISTRIB_faces :
//      if( !rank ) printf("\n--- API mode: Setting face communicators\n");
//      
//      ier = PMMG_Set_numberOfFaceCommunicators(parmesh, ncomm);
//      for( icomm=0; icomm<ncomm; icomm++ ) {
//        ier = PMMG_Set_ithFaceCommunicatorSize(parmesh, icomm,
//                                               color_face_out[icomm],
//                                               nitem_face_comm[icomm]);
//        ier = PMMG_Set_ithFaceCommunicator_faces(parmesh, icomm,
//                                                 out_tria_loc[icomm],
//                                                 out_tria_glob[icomm], 1 );
//      }
//      break;
//    
//    case PMMG_APIDISTRIB_nodes :
//      if( !rank ) printf("\n--- API mode: Setting node communicators\n");
//      
//      ier = PMMG_Set_numberOfNodeCommunicators(parmesh, ncomm);
//      for( icomm=0; icomm<ncomm; icomm++ ) {
//        ier = PMMG_Set_ithNodeCommunicatorSize(parmesh, icomm,
//                                               color_node_out[icomm],
//                                               nitem_node_comm[icomm]);
//        ier = PMMG_Set_ithNodeCommunicator_nodes(parmesh, icomm,
//                                                 out_node_loc[icomm],
//                                                 out_node_glob[icomm], 1 );
//      }
//      break;
//  }
 

  
  /* Color interface triangles with global enumeration */
  if( !color_intfcTria(parmesh,color_face_out,out_tria_loc,out_tria_glob,
                       next_face_comm,nitem_face_comm) ) {
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

//  if ( ier != PMMG_STRONGFAILURE ) {
//    /** ------------------------------ STEP III -------------------------- */
//    /** get results */
//    /** Two solutions: just use the PMMG_saveMesh/PMMG_saveSol functions
//        that will write .mesh(b)/.sol formatted files or manually get your mesh/sol
//        using the PMMG_getMesh/PMMG_getSol functions */
//
//    /** 1) Automatically save the mesh */
//    if ( PMMG_saveMesh_centralized(parmesh,fileout) != 1 ) {
//      fprintf(stdout,"UNABLE TO SAVE MESH\n");
//      ier = PMMG_STRONGFAILURE;
//    }
//
//    /** 2) Automatically save the metric */
//    if ( PMMG_saveMet_centralized(parmesh,fileout) != 1 ) {
//      fprintf(stdout,"UNABLE TO SAVE METRIC\n");
//      ier = PMMG_LOWFAILURE;
//    }
//
//    /** 3) Automatically save the solutions if needed */
//    PMMG_Get_solsAtVerticesSize(parmesh,&nsols,NULL,NULL);
//    if ( nsols ) {
//      if ( PMMG_saveAllSols_centralized(parmesh,fileout) != 1 ) {
//        fprintf(stdout,"UNABLE TO SAVE SOLUTIONS\n");
//        ier = PMMG_LOWFAILURE;
//      }
//    }
//  }
//  else {
//    fprintf(stdout,"BAD ENDING OF PARMMGLIB: UNABLE TO SAVE MESH\n");
//  }

  /** 4) Free the PMMG5 structures */
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

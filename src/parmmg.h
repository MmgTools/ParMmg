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
 * \file parmmg.h
 * \brief internal functions headers for parmmg
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#ifndef _PARMMG_H
#define _PARMMG_H

#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <mpi_pmmg.h>
#include "hdf_pmmg.h"

#include "metis.h"
#include "libparmmg.h"
#include "interpmesh_pmmg.h"
#include "libmmg3d.h"
#include "libmmg3d_private.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \def PMMG_NUL
 *
 * Null value
 *
 */
#define PMMG_NUL     0

/**
 * \def PMMG_NITER
 *
 * Default number of iterations
 *
 */
#define PMMG_NITER   3

/**
 * \def PMMG_IMPRIM
 *
 * Default verbosity
 *
 */
#define PMMG_IMPRIM   1

/**
 * \def PMMG_MMG_IMPRIM
 *
 * Default verbosity for Mmg
 *
 */
#define PMMG_MMG_IMPRIM   -1

/**
 *
 * Size of quality histogram arrays
 *
 */
#define PMMG_QUAL_HISSIZE 5

/**
 *
 * Size of mpi datatype for quality histo computation
 *
 */
#define PMMG_QUAL_MPISIZE 4

/**
 *
 * Size of mpi datatype for length histo computation
 *
 */
#define PMMG_LENSTATS_MPISIZE 12

/**
 *
 * print input quality histogram
 *
 */
#define PMMG_INQUA 1

/**
 *
 * print output quality histogram
 *
 */
#define PMMG_OUTQUA 2

/**
 *
 * no verbosity for pmmg library
 *
 */
#define PMMG_VERB_NO -1

/**
 *
 * minimal verbosity for pmmg library: print library version and duration
 *
 */
#define PMMG_VERB_VERSION 0

/**
 *
 * low verbosity for pmmg library: print main steps + quality histo overview
 *
 */
#define PMMG_VERB_QUAL 1

/**
 *
 * average verbosity for pmmg library: add parmmg steps information
 *
 */
#define PMMG_VERB_STEPS 3

/**
 *
 * average verbosity for pmmg library: add waves information
 *
 */
#define PMMG_VERB_ITWAVES 4

/**
 *
 * detailed verbosity for pmmg library: add detailed quality histo
 *
 */
#define PMMG_VERB_DETQUAL 5

/**
 *
 * Split groups for redistribution (split_grps)
 *
 */
#define PMMG_GRPSPL_DISTR_TARGET 1

/**
 *
 * Split groups for mmg (split_grps)
 *
 */
#define PMMG_GRPSPL_MMG_TARGET 2

/**
 *
 * Use custom partitioning saved in the reference field (1=yes, 0=no)
 *
 */
#define PMMG_PREDEF_PART 0

/**
 * \enum PMMG_Format
 * \brief Type of supported file format
 */
enum PMMG_Format {
  PMMG_FMT_MeditASCII  = MMG5_FMT_MeditASCII, /*!< ASCII Medit (.mesh) */
  PMMG_FMT_MeditBinary = MMG5_FMT_MeditBinary,/*!< Binary Medit (.meshb) */
  PMMG_FMT_GmshASCII   = MMG5_FMT_GmshASCII,  /*!< ASCII Gmsh */
  PMMG_FMT_GmshBinary  = MMG5_FMT_GmshBinary, /*!< Binary Gmsh */
  PMMG_FMT_VtkPvtp     = MMG5_FMT_VtkPvtp,    /*!< VTK pvtp */
  PMMG_FMT_VtkPvtu     = MMG5_FMT_VtkPvtu,    /*!< VTK pvtu */
  PMMG_FMT_VtkVtu      = MMG5_FMT_VtkVtu,     /*!< VTK vtu */
  PMMG_FMT_VtkVtp      = MMG5_FMT_VtkVtp,     /*!< VTK vtp */
  PMMG_FMT_VtkVtk      = MMG5_FMT_VtkVtk,     /*!< VTK vtk */
  PMMG_FMT_Tetgen      = MMG5_FMT_Tetgen,     /*!< Tetgen or Triangle */
  PMMG_FMT_Centralized,                       /*!< Centralized Setters/Getters */
  PMMG_FMT_Distributed,                       /*!< Distributed Setters/Getters */
  PMMG_FMT_DistributedMeditASCII,             /*!< Distributed ASCII Medit (.mesh) */
  PMMG_FMT_DistributedMeditBinary,            /*!< Distributed Binary Medit (.meshb) */
  PMMG_FMT_HDF5,                              /*!< HDF5 format */
  PMMG_FMT_Unknown,                           /*!< Unrecognized */
};

/**< Subgroups target size for a fast remeshing step */
static const int PMMG_REMESHER_TARGET_MESH_SIZE = -30000000;

/**< Subgroups target size for a fast remeshing step */
static const int PMMG_REMESHER_NGRPS_MAX = 100;

/**< Number of metis node per mmg mesh... to test*/
static const int PMMG_RATIO_MMG_METIS = -100;

/**< Subgroups target size for a fast remeshing step */
static const int PMMG_REDISTR_NGRPS_MAX = 1000;

/**< Subgroups minimum size to try to avoid empty partitions */
static const int PMMG_REDISTR_NELEM_MIN = 6;

/**< Allowed imbalance ratio between current and demanded groups size */
static const double PMMG_GRPS_RATIO = 2.0;

/**< Number of elements layers for interface displacement */
static const int PMMG_MVIFCS_NLAYERS = 2;

/**
 * \param parmesh pointer toward a parmesh structure
 * \param val     exit value
 *
 * Controlled parmmg termination:
 *   Deallocate parmesh struct and its allocated members
 *   If this is an unsuccessful exit call abort to cancel any remaining processes
 *   Call MPI_Finalize / exit
 */

#define PMMG_RETURN_AND_FREE(parmesh,val) do                            \
  {                                                                     \
                                                                        \
    if ( !PMMG_Free_all( PMMG_ARG_start,                                \
                         PMMG_ARG_ppParMesh,&parmesh,                   \
                         PMMG_ARG_end) ) {                              \
      fprintf(stderr,"  ## Warning: unable to clean the parmmg memory.\n" \
              " Possible memory leak.\n");                              \
    }                                                                   \
                                                                        \
    MPI_Finalize();                                                     \
    return(val);                                                        \
                                                                        \
  } while(0)

/**
 * Clean the mesh, the metric and the solutions and return \a val.
 */
#define PMMG_CLEAN_AND_RETURN(parmesh,val)do                            \
  {                                                                     \
    int kgrp, ksol;                                                     \
                                                                        \
    for ( kgrp=0; kgrp<parmesh->ngrp; ++kgrp ) {                        \
      if ( parmesh->listgrp[kgrp].mesh ) {                              \
        parmesh->listgrp[kgrp].mesh->npi = parmesh->listgrp[kgrp].mesh->np; \
        parmesh->listgrp[kgrp].mesh->nti = parmesh->listgrp[kgrp].mesh->nt; \
        parmesh->listgrp[kgrp].mesh->nai = parmesh->listgrp[kgrp].mesh->na; \
        parmesh->listgrp[kgrp].mesh->nei = parmesh->listgrp[kgrp].mesh->ne; \
      }                                                                 \
                                                                        \
      if ( parmesh->listgrp[kgrp].met )                                 \
        parmesh->listgrp[kgrp].met->npi  = parmesh->listgrp[kgrp].met->np; \
                                                                        \
      if ( parmesh->listgrp[kgrp].ls )                                 \
        parmesh->listgrp[kgrp].ls->npi  = parmesh->listgrp[kgrp].ls->np; \
                                                                        \
      if ( parmesh->listgrp[kgrp].mesh ) {                              \
        for ( ksol=0; ksol<parmesh->listgrp[kgrp].mesh->nsols; ++ksol ) { \
          parmesh->listgrp[kgrp].field[ksol].npi  = parmesh->listgrp[kgrp].field[ksol].np; \
        }                                                               \
      }                                                                 \
    }                                                                   \
                                                                        \
    return val;                                                         \
                                                                        \
  }while(0)


#define ERROR_AT(msg1,msg2)                                          \
  fprintf( stderr, "%s %s function: %s, file: %s, line: %d \n", \
           msg1, msg2, __func__, __FILE__, __LINE__ )

#define MEM_CHK_AVAIL(mesh,bytes,msg) do {                            \
  if ( (mesh)->memCur + (bytes) > (mesh)->memMax ) {                  \
    char diag[1024];                                                  \
    snprintf(diag, 1024, " Allocation of %ld bytes exceeds max %ld: ", bytes, (mesh)->memMax); \
    ERROR_AT(msg, diag);                                              \
    stat = PMMG_FAILURE;                                              \
  } else if ( (mesh)->memCur + (bytes) < 0  ) {                       \
    ERROR_AT(msg," Tried to free more mem than allocated: " );        \
    stat = PMMG_SUCCESS;                                              \
  }                                                                   \
  else {                                                              \
    stat = PMMG_SUCCESS;                                              \
  } } while(0)

#define PMMG_DEL_MEM(mesh,ptr,type,msg) do {                \
    size_t size_to_free;                                    \
                                                            \
    if ( ptr ) {                                            \
      size_to_free = myfree( ptr );                         \
      assert ( (mesh)->memCur >= size_to_free );            \
      (mesh)->memCur -= size_to_free;                       \
      (ptr) = NULL;                                         \
    }                                                       \
  } while(0)

#define PMMG_MALLOC(mesh,ptr,size,type,msg,on_failure) do { \
  int    stat = PMMG_SUCCESS;                               \
  size_t size_to_allocate;                                  \
                                                            \
  (ptr) = NULL;                                             \
  if ( (size) != 0 ) {                                      \
    size_to_allocate = (size)*sizeof(type);                 \
    MEM_CHK_AVAIL(mesh,size_to_allocate,msg );              \
    if ( stat == PMMG_SUCCESS ) {                           \
      (ptr) = (type*)mymalloc( size_to_allocate );          \
      if ( (ptr) == NULL ) {                                \
        ERROR_AT( msg, " malloc failed: " );                \
        on_failure;                                         \
      } else {                                              \
        (mesh)->memCur += size_to_allocate;                 \
        stat = PMMG_SUCCESS;                                \
      }                                                     \
    } else {                                                \
      on_failure;                                           \
    }                                                       \
  } } while(0)

#define PMMG_CALLOC(mesh,ptr,size,type,msg,on_failure) do { \
  int    stat = PMMG_SUCCESS;                               \
  size_t size_to_allocate;                                  \
                                                            \
  (ptr) = NULL;                                             \
  if ( (size) != 0 ) {                                      \
    size_to_allocate = (size)*sizeof(type);                 \
    MEM_CHK_AVAIL(mesh,size_to_allocate,msg);               \
    if ( stat == PMMG_SUCCESS ) {                           \
      (ptr) = (type*)mycalloc( (size), sizeof(type) );      \
      if ( (ptr) == NULL ) {                                \
        ERROR_AT(msg," calloc failed: ");                   \
        on_failure;                                         \
      } else {                                              \
        (mesh)->memCur += size_to_allocate;                 \
      }                                                     \
    } else {                                                \
      on_failure;                                           \
    }                                                       \
  } } while(0)

#define PMMG_REALLOC(mesh,ptr,newsize,oldsize,type,msg,on_failure) do { \
  int    stat = PMMG_SUCCESS;                                           \
  size_t size_to_allocate,size_to_add,size_to_increase;                 \
  type*  tmp;                                                           \
                                                                        \
  if ( (ptr) == NULL ) {                                                \
    assert(((oldsize)==0) && "NULL pointer pointing to non 0 sized memory?"); \
    PMMG_MALLOC(mesh,ptr,(newsize),type,msg,on_failure);                \
  } else if ((newsize)==0) {                                            \
    PMMG_DEL_MEM(mesh,ptr,type,msg);                                    \
  } else if ((newsize) < (oldsize)) {                                   \
    size_to_allocate = (newsize)*sizeof(type);                          \
    tmp = (type *)myrealloc((ptr),size_to_allocate,                     \
                            (oldsize)*sizeof(type));                    \
    if ( tmp == NULL ) {                                                \
      ERROR_AT(msg," Realloc failed: ");                                \
      PMMG_DEL_MEM(mesh,ptr,type,msg);                                  \
      on_failure;                                                       \
    } else {                                                            \
      (ptr) = tmp;                                                      \
      (mesh)->memCur -= (((oldsize)*sizeof(type))-size_to_allocate);    \
    }                                                                   \
  } else if ((newsize) > (oldsize)) {                                   \
    size_to_add = ((newsize)-(oldsize))*sizeof(type);                   \
    size_to_allocate = (newsize)*sizeof(type);                          \
    size_to_increase = (oldsize)*sizeof(type);                          \
                                                                        \
    MEM_CHK_AVAIL(mesh,size_to_add,msg);                                \
    if ( stat == PMMG_SUCCESS ) {                                       \
      tmp = (type *)myrealloc((ptr),size_to_allocate,size_to_increase); \
      if ( tmp == NULL ) {                                              \
        ERROR_AT(msg, " Realloc failed: " );                            \
        PMMG_DEL_MEM(mesh,ptr,type,msg);                                \
        on_failure;                                                     \
      } else {                                                          \
        (ptr) = tmp;                                                    \
        (mesh)->memCur += ( size_to_add );                              \
      }                                                                 \
    }                                                                   \
    else {                                                              \
      on_failure;                                                       \
    }                                                                   \
  }                                                                     \
  } while(0)

#define PMMG_RECALLOC(mesh,ptr,newsize,oldsize,type,msg,on_failure) do {                \
    int my_stat = PMMG_SUCCESS;                                                         \
                                                                                        \
    PMMG_REALLOC(mesh,ptr,newsize,oldsize,type,msg,my_stat=PMMG_FAILURE;on_failure;);   \
    if ( (my_stat == PMMG_SUCCESS ) && ((newsize) > (oldsize)) ) {                      \
      memset( (ptr) + oldsize, 0, ((size_t)((newsize)-(oldsize)))*sizeof(type));        \
    }                                                                                   \
  } while(0)


/**
 * \param parmesh pointer toward a parmesh structure
 * \param mesh pointer toward a mesh structure
 * \param on_failure instruction to execute if fail
 *
 * Check the allowed memory. */
#define PMMG_MEM_CHECK(parmesh,mesh,on_failure) do {             \
    size_t memGloMax,memMax;                                      \
    memGloMax = parmesh->memGloMax;                               \
    memMax = mesh->memMax;                                        \
    if( memMax != memGloMax ) {                                   \
      fprintf(stderr,"\n  ## Error: %s: allowed memory mismatch." \
                     " Maximal: %zu -- global: %zu\n",            \
              __func__,memMax,memGloMax);                         \
      on_failure;                                                 \
    }                                                             \
  } while(0)


/* Input */
int PMMG_Set_name(PMMG_pParMesh,char **,const char* name,const char* defname);
int PMMG_check_inputData ( PMMG_pParMesh parmesh );
int PMMG_preprocessMesh( PMMG_pParMesh parmesh );
int PMMG_preprocessMesh_distributed( PMMG_pParMesh parmesh );
int PMMG_parsar( int argc, char *argv[], PMMG_pParMesh parmesh );
void PMMG_setfunc( PMMG_pParMesh parmesh );

/* Mesh analysis */
void PMMG_Analys_Init_SurfNormIndex( MMG5_pTetra pt );
int PMMG_Analys_Get_SurfNormalIndex( MMG5_pTetra pt,int ifac,int i );
int PMMG_boulernm(PMMG_pParMesh parmesh,MMG5_pMesh mesh,MMG5_Hash *hash,int start,int ip,int *ng,int *nr);
int PMMG_boulen(PMMG_pParMesh parmesh,MMG5_pMesh mesh,int start,int ip,int iface,double t[3]);
int PMMG_analys_tria(PMMG_pParMesh parmesh,MMG5_pMesh mesh);
int PMMG_analys(PMMG_pParMesh parmesh,MMG5_pMesh mesh,MPI_Comm comm);
int PMMG_update_analys(PMMG_pParMesh parmesh);
int PMMG_hashPar( MMG5_pMesh mesh,MMG5_HGeom *pHash );
int PMMG_hashPar_pmmg( PMMG_pParMesh parmesh,MMG5_HGeom *pHash );
int PMMG_hashOldPar_pmmg( PMMG_pParMesh parmesh,MMG5_pMesh mesh,MMG5_Hash *hash );

/* Isovalue discretization functions */
int  PMMG_ls(PMMG_pParMesh parmesh);
int  PMMG_cuttet_ls(PMMG_pParMesh parmesh);
int  PMMG_resetRef_ls(PMMG_pParMesh parmesh,MMG5_pMesh mesh);
int  PMMG_setref_ls(PMMG_pParMesh parmesh,MMG5_pMesh mesh, MMG5_pSol sol);
int  PMMG_snpval_ls(PMMG_pParMesh parmesh,MMG5_pMesh mesh,MMG5_pSol sol);

int      PMMG_hashUpdate_all(MMG5_Hash *hash,MMG5_int a,MMG5_int b,MMG5_int k,MMG5_int s);
MMG5_int PMMG_hashGet_all(MMG5_Hash *hash,MMG5_int a,MMG5_int b,MMG5_int *k,MMG5_int *s);

void PMMG_nosplit_sort(MMG5_pMesh mesh,MMG5_int k,int ifac,MMG5_int *tetra_sorted,MMG5_int *node_sorted);
void PMMG_split1_sort(MMG5_pMesh mesh,MMG5_int k,int ifac,uint8_t tau[4],MMG5_int ne_tmp,MMG5_int *tetra_sorted,MMG5_int *node_sorted);
void PMMG_split2sf_sort(MMG5_pMesh mesh,MMG5_int k,int ifac,uint8_t tau[4],int imin,MMG5_int ne_tmp,MMG5_int *tetra_sorted,MMG5_int *node_sorted);
void PMMG_split3cone_sort(MMG5_pMesh mesh,MMG5_int k,int ifac,uint8_t tau[4],int ia,int ib,MMG5_int ne_tmp,MMG5_int *tetra_sorted,MMG5_int *node_sorted);
void PMMG_split4op_sort(MMG5_pMesh mesh,MMG5_int k,int ifac,uint8_t tau[4],int imin01,int imin23,MMG5_int ne_tmp,MMG5_int *tetra_sorted,MMG5_int *node_sorted);
int  PMMG_sort_vertices(MMG5_pMesh mesh,MMG5_int k,MMG5_int *v_t,int ifac);
void PMMG_sort_tetra(MMG5_int *tetra,MMG5_int *node,MMG5_int *v_t0,MMG5_int *v_t1,MMG5_int *v_t2);
void PMMG_swap_vertices(MMG5_int *a);
void PMMG_swap_ints(int *a, int *b);
void PMMG_swap_3int_arrays(int *a, int *b);
int  PMMG_compare_3ints_array(int *a, int *b);

/* Internal library */
void PMMG_setfunc( PMMG_pParMesh parmesh );
int PMMG_parmmglib1 ( PMMG_pParMesh parmesh );

/* Mesh distrib */
int PMMG_bdryUpdate( MMG5_pMesh mesh );
int PMMG_bcast_mesh ( PMMG_pParMesh parmesh );
int PMMG_partBcast_mesh( PMMG_pParMesh parmesh );
int PMMG_splitPart_grps( PMMG_pParMesh,int,int,int );
int PMMG_split_grps( PMMG_pParMesh parmesh,int grpIdOld,int ngrp,idx_t *part,int fitMesh );

/* Load Balancing */
int PMMG_interactionMap(PMMG_pParMesh parmesh,int **interactions,int **interaction_map);
int PMMG_transfer_all_grps(PMMG_pParMesh parmesh,idx_t *part,int);
int PMMG_distribute_grps( PMMG_pParMesh parmesh,int partitioning_mode );
int PMMG_loadBalancing( PMMG_pParMesh parmesh,int partitioning_mode );
int PMMG_split_n2mGrps( PMMG_pParMesh,int,int,int );
double PMMG_computeWgt( MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra pt,int ifac );
void PMMG_computeWgt_mesh( MMG5_pMesh mesh,MMG5_pSol met,int tag );

/* Mesh interpolation */
int PMMG_oldGrps_newGroup( PMMG_pParMesh parmesh,int igrp );
int PMMG_oldGrps_fillGroup( PMMG_pParMesh parmesh,int igrp );
int PMMG_update_oldGrps( PMMG_pParMesh parmesh );
int PMMG_interpMetricsAndFields( PMMG_pParMesh parmesh,int* );
int PMMG_copyMetricsAndFields_point( MMG5_pMesh mesh, MMG5_pMesh oldMesh, MMG5_pSol met, MMG5_pSol oldMet, MMG5_pSol,MMG5_pSol, int* permNodGlob,uint8_t);

/* Communicators building and unallocation */
void PMMG_parmesh_int_comm_free( PMMG_pParMesh,PMMG_pInt_comm);
void PMMG_parmesh_ext_comm_free( PMMG_pParMesh,PMMG_pExt_comm,int);
void PMMG_grp_comm_free( PMMG_pParMesh ,int**,int**,int*);
void PMMG_node_comm_free( PMMG_pParMesh );
void PMMG_edge_comm_free( PMMG_pParMesh );
int PMMG_Compute_verticesGloNum( PMMG_pParMesh parmesh,MPI_Comm comm );
int PMMG_Compute_trianglesGloNum( PMMG_pParMesh parmesh,MPI_Comm comm );
int PMMG_color_commNodes( PMMG_pParMesh parmesh,MPI_Comm comm );
void PMMG_tria2elmFace_flags( PMMG_pParMesh parmesh );
void PMMG_tria2elmFace_coords( PMMG_pParMesh parmesh );
int PMMG_tria_highestcoord( MMG5_pMesh mesh, MMG5_int *v_t);
int PMMG_build_nodeCommIndex( PMMG_pParMesh parmesh );
int PMMG_build_faceCommIndex( PMMG_pParMesh parmesh );
int PMMG_build_nodeCommFromFaces( PMMG_pParMesh parmesh, MPI_Comm comm );
int PMMG_build_faceCommFromNodes( PMMG_pParMesh parmesh, MPI_Comm comm );
int PMMG_build_simpleExtNodeComm( PMMG_pParMesh parmesh );
int PMMG_build_intNodeComm( PMMG_pParMesh parmesh );
int PMMG_build_completeExtNodeComm( PMMG_pParMesh parmesh, MPI_Comm comm );
int PMMG_build_edgeComm( PMMG_pParMesh,MMG5_pMesh,MMG5_HGeom *hpar,MPI_Comm);
int PMMG_build_completeExtEdgeComm( PMMG_pParMesh parmesh, MPI_Comm comm  );

int PMMG_pack_faceCommunicators(PMMG_pParMesh parmesh);
int PMMG_pack_nodeCommunicators(PMMG_pParMesh parmesh);

/* Communicators checks */
int PMMG_check_intFaceComm( PMMG_pParMesh parmesh );
int PMMG_check_extFaceComm( PMMG_pParMesh parmesh, MPI_Comm comm );
int PMMG_check_intNodeComm( PMMG_pParMesh parmesh);
int PMMG_check_extNodeComm( PMMG_pParMesh parmesh, MPI_Comm comm );
int PMMG_check_extEdgeComm( PMMG_pParMesh parmesh, MPI_Comm comm );

/* Tags */
void PMMG_tag_par_node(MMG5_pPoint ppt);
void PMMG_tag_par_edge(MMG5_pxTetra pxt,int j);
void PMMG_tag_par_face(MMG5_pxTetra pxt,int j);
void PMMG_tag_par_tria(MMG5_pTria ptt);
void PMMG_untag_par_node(MMG5_pPoint ppt);
void PMMG_untag_par_edge(MMG5_pxTetra pxt,int j);
void PMMG_untag_par_face(MMG5_pxTetra pxt,int j);
int  PMMG_resetOldTag(PMMG_pParMesh parmesh);
int  PMMG_updateTag(PMMG_pParMesh parmesh);
void PMMG_updateTagRef_node(PMMG_pParMesh parmesh,MMG5_pMesh mesh);
int  PMMG_parbdySet( PMMG_pParMesh parmesh );
int  PMMG_parbdyTria( PMMG_pParMesh parmesh );

/* Mesh merge */
int PMMG_mergeGrpJinI_interfacePoints_addGrpJ( PMMG_pParMesh,PMMG_pGrp,PMMG_pGrp);
int PMMG_mergeGrps_interfacePoints( PMMG_pParMesh parmesh );
int PMMG_mergeGrpJinI_internalPoints( PMMG_pGrp,PMMG_pGrp grpJ );
int PMMG_mergeGrpJinI_interfaceTetra( PMMG_pParMesh,PMMG_pGrp,PMMG_pGrp );
int PMMG_mergeGrpJinI_internalTetra( PMMG_pGrp,PMMG_pGrp );
int PMMG_merge_grps ( PMMG_pParMesh parmesh,int );

/* Move interfaces */
int PMMG_part_getInterfaces( PMMG_pParMesh parmesh,int *part,int *ngrps,int target );
int PMMG_part_getProcs( PMMG_pParMesh parmesh,int *part );
int PMMG_fix_contiguity( PMMG_pParMesh parmesh,int *counter );
int PMMG_fix_contiguity_centralized( PMMG_pParMesh parmesh,idx_t *part );
int PMMG_fix_contiguity_split( PMMG_pParMesh parmesh,idx_t ngrp,idx_t *part );
int PMMG_part_moveInterfaces( PMMG_pParMesh parmesh,int *vtxdist,int *map,int *base_front );
int PMMG_mark_interfacePoints( PMMG_pParMesh parmesh,MMG5_pMesh mesh,int* vtxdist,int* priorityMap );
int PMMG_init_ifcDirection( PMMG_pParMesh parmesh,int **vtxdist,int **map );
int PMMG_set_ifcDirection( PMMG_pParMesh parmesh,int **vtxdist,int **map );
int PMMG_get_ifcDirection( PMMG_pParMesh parmesh,int *vtxdist,int *map,int color0,int color1 );

/* Packing */
int PMMG_update_node2intPackedTetra( PMMG_pGrp grp );
int PMMG_mark_packedTetra(MMG5_pMesh mesh,int *ne);
int PMMG_update_node2intPackedVertices( PMMG_pGrp grp );
int PMMG_packTetra ( PMMG_pParMesh parmesh, int igrp );

/* Memory */
int  PMMG_link_mesh( MMG5_pMesh mesh );
void PMMG_listgrp_free( PMMG_pParMesh parmesh, PMMG_pGrp *listgrp, int ngrp );
void PMMG_grp_free( PMMG_pParMesh parmesh, PMMG_pGrp grp );
int  PMMG_parmesh_SetMemMax( PMMG_pParMesh parmesh);
int  PMMG_setMeshSize( MMG5_pMesh,int,int,int,int,int );
int  PMMG_setMeshSize_alloc( MMG5_pMesh );
int  PMMG_setMeshSize_realloc( MMG5_pMesh,int,int,int,int);
int  PMMG_fitMeshSize( PMMG_pParMesh parmesh, PMMG_pGrp );
int  PMMG_updateMeshSize( PMMG_pParMesh parmesh,int fitMesh);
void PMMG_parmesh_SetMemGloMax( PMMG_pParMesh parmesh );
void PMMG_parmesh_Free_Comm( PMMG_pParMesh parmesh );
void PMMG_parmesh_Free_Listgrp( PMMG_pParMesh parmesh );
void PMMG_destroy_int( PMMG_pParMesh,void **ptr[],size_t,char*);
int  PMMG_clean_emptyMesh( PMMG_pParMesh parmesh, PMMG_pGrp listgrp, int ngrp );
int  PMMG_resize_extComm ( PMMG_pParMesh,PMMG_pExt_comm,int,int* );
int  PMMG_resize_extCommArray ( PMMG_pParMesh,PMMG_pExt_comm*,int,int*);

/* Tools */
int PMMG_copy_mmgInfo ( MMG5_Info *info, MMG5_Info *info_cpy );

/* Quality */
int PMMG_qualhisto( PMMG_pParMesh parmesh,int,int,MPI_Comm comm );
int PMMG_prilen( PMMG_pParMesh parmesh,int8_t,int,MPI_Comm comm );
int PMMG_tetraQual( PMMG_pParMesh parmesh,int8_t metRidTyp );

/* Variadic_pmmg.c */
int PMMG_Init_parMesh_var_internal(va_list argptr,int callFromC);
int PMMG_Free_all_var(va_list argptr);

const char* PMMG_Get_pmmgArgName(int typArg);

/* Private I/Os and APIs*/
int PMMG_loadMesh_hdf5_i(PMMG_pParMesh parmesh, int *load_entities, const char *filename);
int PMMG_saveMesh_hdf5_i(PMMG_pParMesh parmesh, int *save_entities, const char *filename);
int PMMG_Set_defaultIOEntities_i(int io_entities[PMMG_IO_ENTITIES_size] );
int PMMG_Set_IOEntities_i(int io_entities[PMMG_IO_ENTITIES_size], int target, int val);
int PMMG_Get_format( char *ptr, int fmt );

#ifdef __cplusplus
}
#endif

#endif

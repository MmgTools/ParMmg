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

#include "libparmmg.h"

#ifdef __cplusplus
extern "C" {
#endif


#define ERROR_AT(msg1,msg2)                                          \
  fprintf( stderr, msg1 msg2 " function: %s, file: %s, line: %d \n", \
           __func__, __FILE__, __LINE__ )

#define MEM_CHK_AVAIL(mesh,bytes,msg) do {                            \
  if ( (mesh)->memCur + (bytes) > (mesh)->memMax ) {                  \
    ERROR_AT(msg," Exceeded max memory allowed: " );                  \
    stat = PMMG_FAILURE;                                              \
  } else if ( (mesh)->memCur + (bytes) < 0  ) {                       \
    ERROR_AT(msg," Tried to free more mem than allocated: " );        \
    stat = PMMG_FAILURE;                                              \
  }                                                                   \
  else {                                                              \
    stat = PMMG_SUCCESS;                                              \
  } } while(0)

#define PMMG_DEL_MEM(mesh,ptr,size,type,msg) do {           \
    int stat = PMMG_SUCCESS;                                \
    if ( (size) != 0 ) {                                    \
      long long size_to_free = -(size)*sizeof(type);        \
      MEM_CHK_AVAIL(mesh,size_to_free,msg);                 \
      if ( stat == PMMG_SUCCESS )                           \
        (mesh)->memCur += size_to_free;                     \
      free( ptr );                                          \
      (ptr) = NULL;                                         \
    }                                                       \
  } while(0)

#define PMMG_MALLOC(mesh,ptr,size,type,msg,on_failure) do { \
  int stat = PMMG_SUCCESS;                                  \
  if ( (size) != 0 ) {                                      \
    long long size_to_allocate = (size)*sizeof(type);       \
    MEM_CHK_AVAIL(mesh,size_to_allocate,msg );              \
    if ( stat == PMMG_SUCCESS ) {                           \
      (ptr) = malloc( size_to_allocate );                   \
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
  int stat = PMMG_SUCCESS;                                  \
  if ( (size) != 0 ) {                                      \
    long long size_to_allocate = (size)*sizeof(type);       \
    MEM_CHK_AVAIL(mesh,size_to_allocate,msg);               \
    if ( stat == PMMG_SUCCESS ) {                           \
      (ptr) = calloc( (size), sizeof(type) );               \
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
    int stat = PMMG_SUCCESS;                                            \
    if ( (ptr) == NULL ) {                                              \
      assert(((oldsize)==0) && "NULL pointer pointing to non 0 sized memory?");\
      PMMG_MALLOC(mesh,ptr,(newsize),type,msg,on_failure);              \
    } else if ((newsize)==0) {                                          \
      PMMG_DEL_MEM(mesh,ptr,(oldsize),type,msg);                        \
    } else if ((newsize) < (oldsize)) {                                 \
      long long size_to_allocate = (newsize)*sizeof(type);              \
      (ptr) = realloc( (ptr), size_to_allocate);                        \
      if ( ptr == NULL ) {                                              \
        ERROR_AT(msg," Realloc failed: ");                              \
        on_failure;                                                     \
      } else {                                                          \
        (mesh)->memCur -= ((long long)((oldsize)*sizeof(type))-size_to_allocate);\
      }                                                                 \
    } else if ((newsize) > (oldsize)) {                                 \
long long size_to_allocate = ((long long)(newsize)-(long long)(oldsize))*sizeof(type);\
      MEM_CHK_AVAIL(mesh,size_to_allocate,msg);                         \
      if ( stat == PMMG_SUCCESS ) {                                     \
        (ptr) = realloc((ptr), (long long)(newsize) * sizeof(type));    \
        if ( (ptr) == NULL ) {                                          \
          ERROR_AT(msg, " Realloc failed: " );                          \
          on_failure;                                                   \
        } else {                                                        \
          (mesh)->memCur += ( size_to_allocate );      \
        }                                                               \
      }                                                                 \
    }                                                                   \
  } while(0)

#define PMMG_RECALLOC(mesh,ptr,newsize,oldsize,type,msg,on_failure) do { \
  int my_stat = PMMG_SUCCESS;                                            \
  PMMG_REALLOC(mesh,ptr,newsize,oldsize,type,msg,on_failure;my_stat=PMMG_FAILURE);\
  if ( (my_stat == PMMG_SUCCESS ) && ((newsize) > (oldsize)) )           \
    memset( (ptr) + oldsize, 0, ((size_t)((newsize)-(oldsize)))*sizeof(type));     \
  } while(0)


/* Input */
int PMMG_check_inputData ( PMMG_pParMesh parmesh );
int PMMG_parsar( int argc, char *argv[], PMMG_pParMesh parmesh );

/* Internal library */
int PMMG_parmmglib1 ( PMMG_pParMesh parmesh );

/* Mesh distrib */
int PMMG_bdryUpdate( MMG5_pMesh mesh );
int PMMG_bcast_mesh ( PMMG_pParMesh parmesh );
int PMMG_grpSplit_setMeshSize( MMG5_pMesh,int,int,int,int,int );
int PMMG_split_grps( PMMG_pParMesh,int,int );

/* Load Balancing */
int PMMG_distribute_grps( PMMG_pParMesh parmesh );
int PMMG_loadBalancing( PMMG_pParMesh parmesh );
int PMMG_split_n2mGrps( PMMG_pParMesh,int,int );

/* Communicators building and unallocation */
void PMMG_int_comm_free( PMMG_pParMesh,PMMG_pint_comm);
void PMMG_ext_comm_free( PMMG_pParMesh,PMMG_pext_comm,int);
void PMMG_grp_comm_free( PMMG_pParMesh ,int**,int**,int*);
void PMMG_node_comm_free( PMMG_pParMesh );

int PMMG_build_nodeCommFromFaces( PMMG_pParMesh parmesh );
int PMMG_build_simpleExtNodeComm( PMMG_pParMesh parmesh );
int PMMG_build_intNodeComm( PMMG_pParMesh parmesh );
int PMMG_build_completeExtNodeComm( PMMG_pParMesh parmesh );

int PMMG_pack_faceCommunicators(PMMG_pParMesh parmesh);
int PMMG_pack_nodeCommunicators(PMMG_pParMesh parmesh);

/* Communicators checks */
int PMMG_check_intFaceComm( PMMG_pParMesh parmesh );
int PMMG_check_extFaceComm( PMMG_pParMesh parmesh );
int PMMG_check_intNodeComm( PMMG_pParMesh parmesh );
int PMMG_check_extNodeComm( PMMG_pParMesh parmesh );

/* Mesh merge */
int PMMG_mergeGrpJinI_interfacePoints_addGrpJ( PMMG_pParMesh,PMMG_pGrp,PMMG_pGrp);
int PMMG_mergeGrps_interfacePoints( PMMG_pParMesh parmesh );
int PMMG_mergeGrpJinI_internalPoints( PMMG_pGrp,PMMG_pGrp grpJ );
int PMMG_mergeGrpJinI_interfaceTetra( PMMG_pParMesh,PMMG_pGrp,PMMG_pGrp );
int PMMG_mergeGrpJinI_internalTetra( PMMG_pGrp,PMMG_pGrp );
int PMMG_merge_grps ( PMMG_pParMesh parmesh );

/* Packing */
int PMMG_update_node2intPackedTetra( PMMG_pGrp grp );
int PMMG_mark_packedTetra(MMG5_pMesh mesh,int *ne);
int PMMG_update_node2intPackedVertices( PMMG_pGrp grp );
int PMMG_packTetra ( PMMG_pParMesh parmesh, int igrp );

/* Memory */
int  PMMG_link_mesh( MMG5_pMesh mesh );
void PMMG_listgrp_free( PMMG_pParMesh parmesh, PMMG_pGrp *listgrp, int ngrp );
void PMMG_grp_free( PMMG_pParMesh parmesh, PMMG_pGrp grp );
int  PMMG_parmesh_SetMemMax( PMMG_pParMesh parmesh, int percent);
int  PMMG_parmesh_updateMemMax( PMMG_pParMesh parmesh,int percent,int fitMesh);
void PMMG_parmesh_SetMemGloMax( PMMG_pParMesh parmesh, long long int memReq );
void PMMG_parmesh_Free_Comm( PMMG_pParMesh parmesh );
void PMMG_parmesh_Free_Listgrp( PMMG_pParMesh parmesh );

void PMMG_exit_and_free( PMMG_pParMesh parmesh, const int val );
#ifdef __cplusplus
}
#endif

#endif

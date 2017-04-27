/**
 * \file parmmg.h
 * \brief internal functions headers for parmmg
 * \author
 * \version
 * \date 11 2016
 * \copyright
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

#include "libparmmg.h"

#ifdef __cplusplus
extern "C" {
#endif


/** Free allocated pointers of parmesh, call MPI_Finalize and return value val */
#define PMMG_RETURN_AND_FREE( parmesh, val ) do { int k;       \
  for ( k = 0; k < parmesh->ngrp; ++k ) {                      \
    grp  = &parmesh->listgrp[k];                               \
    mesh = grp->mesh;                                          \
    met  = grp->met;                                           \
    MMG3D_Free_all( MMG5_ARG_start,                            \
                    MMG5_ARG_ppMesh, &parmesh->listgrp[k].mesh,\
                    MMG5_ARG_ppMet, &parmesh->listgrp[k].met,  \
                    MMG5_ARG_end );                            \
  }                                                            \
  _MMG5_SAFE_FREE(parmesh);                                    \
  MPI_Finalize();                                              \
  return(val);                                                 \
  } while(0)

#define PMMG_ERRORMSG(message1, message2) fprintf( stderr, message1 message2   \
  "function: %s, file: %s, line: %d \n", __func__, __FILE__, __LINE__)

#define PMMG_MEM_CHK_AVAIL(parmesh,bytes,message,stat) do {                                           \
  if (   ( (parmesh)->memCur + (bytes) > (parmesh)->memMax )                                          \
      || ( (parmesh)->memCur + (bytes) < 0 ) ) {                                                      \
    PMMG_ERRORMSG(message, " Exceeded max memory allowed or tried to free more mem than allocated: ");\
    stat = 0;                                                                                         \
  } } while ( 0 )

#define PMMG_FREE(parmesh,ptr,bytes,message) do {   \
 int stat = 1;                                      \
 PMMG_MEM_CHK_AVAIL(parmesh,-(bytes),message,stat );\
 if ( stat ) {                                      \
   (parmesh)->memCur -= (bytes);                    \
   free( ptr );                                     \
   ptr = NULL;                                      \
 } } while( 0 )

#define PMMG_MALLOC(parmesh,ptr,bytes,message) do { \
 int stat = 1;                                      \
 PMMG_MEM_CHK_AVAIL(parmesh,(bytes),message,stat);  \
 if ( stat ) {                                      \
   ptr = malloc( (bytes) );                         \
   if ( ptr == NULL ) {                             \
     PMMG_ERRORMSG(message, " Malloc failed: ");    \
     return ( 0 );                                  \
   } else {                                         \
    (parmesh)->memCur += (bytes);                   \
 } } } while( 0 )

#define PMMG_CALLOC(parmesh,ptr,size,type,message) do {         \
  int stat = 1;                                                 \
  PMMG_MEM_CHK_AVAIL(parmesh,(size)*sizeof(type),message,stat); \
  if ( stat ) {                                                 \
   ptr = calloc( (size), sizeof(type) );                        \
   if ( ptr == NULL ) {                                         \
     PMMG_ERRORMSG(message, " calloc failed: ");                \
     return ( 0 );                                              \
   } else {                                                     \
    (parmesh)->memCur += (size) * sizeof(type);                 \
   } } } while( 0 )

#define PMMG_REALLOC(mesh,ptr,newsize,oldsize,type,message) do { \
    puts("ADD ME"; exit(EXIT_FAILURE);                           \
  } while( 0 )

#define PMMG_RECALLOC(mesh,ptr,newsize,oldsize,type,message) do { \
    puts("ADD ME"; exit(EXIT_FAILURE);                            \
  } while( 0 )


void PMMG_Init_parameters( PMMG_pParMesh mesh );

int PMMG_Init_parMesh_var( va_list argptr );
int PMMG_check_inputData ( PMMG_pParMesh parmesh );

int PMMG_parmmglib1 ( PMMG_pParMesh parmesh );

int  PMMG_mergeGrps ( PMMG_pParMesh parmesh );
int  PMMG_bcastMesh ( PMMG_pParMesh parmesh );
int  PMMG_bdryStore ( MMG5_pMesh mesh );
int  PMMG_bdryUpdate( MMG5_pMesh mesh );

/**
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 1.
 *
 * Store command line arguments.
 *
 * \remark no matching fortran function.
 *
 */
int PMMG_parsar( int argc, char *argv[], PMMG_pParMesh parmesh );

#ifdef __cplusplus
}
#endif

#endif

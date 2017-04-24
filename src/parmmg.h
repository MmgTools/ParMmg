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
#define PMMG_FREE_AND_RETURN( parmesh, val ) do { int k;       \
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

#define PMMG_MEM_ERRORMSG  do {                                                \
    fprintf( stderr,                                                           \
          " Exceeded max memory allowed or tried to free more than allocated. "\
          "function: %s, file: %s, line: %d \n", __func__, __FILE__, __LINE__);\
    return ( PMMG_STRONGFAILURE );                                             \
  } while ( 0 )

#define PMMG_FREE ( mesh, ptr, size, type, message ) do {                    \
    int status = 1;                                                          \
    (mesh)->memCur -= (long long)( size * sizeof(type) );                    \
    _MMG5_CHK_MEM ( mesh, size, message, status = 0 );                       \
    if ( status ) {                                                          \
      _MMG5_SAFE_FREE( ptr );                                                \
    } else {                                                                 \
      (mesh)->memCur += (long long)( size );                                 \
      MMG_MEM_ERRORMSG;                                                      \
    }                                                                        \
  } while( 0 )

#define PMMG_MALLOC ( mesh, ptr, size, type, message ) do {                  \
    int status = 1;                                                          \
    _MMG5_ADD_MEM ( mesh, size * sizeof(type), message, status = 0 )         \
    if ( status ) {                                                          \
      _MMG5_SAFE_MALLOC( ptr, size, type );                                  \
    } else {                                                                 \
      MMG_MEM_ERRORMSG;                                                      \
    }                                                                        \
  } while( 0 )

#define PMMG_CALLOC ( mesh, ptr, size, type, message ) do {                  \
    int status = 1;                                                          \
    _MMG5_ADD_MEM( mesh, size * sizeof(type), message, status = 0 );         \
    if ( status ) {                                                          \
      _MMG5_SAFE_CALLOC( ptr, size, type );                                  \
    } else {                                                                 \
      MMG_MEM_ERRORMSG;                                                      \
    }                                                                        \
  } while( 0 )

#define PMMG_REALLOC ( mesh, ptr, newsize, oldsize, type, message ) do {     \
    int status = 1;                                                          \
    _MMG5_ADD_MEM( mesh, newsize * sizeof(type), message, status = 0 );      \
    if ( status ) {                                                          \
      _MMG5_SAFE_REALLOC( ptr, newsize, type, message );                     \
      (mesh)->memCur -= (long long)( oldsize * sizeof(type) );               \
    } else {                                                                 \
      MMG_MEM_ERRORMSG;                                                      \
    }                                                                        \
  } while( 0 )

#define PMMG_RECALLOC ( mesh, ptr, newsize, oldsize, type, message ) do {    \
    int status = 1;                                                          \
    _MMG5_ADD_MEM( mesh, newsize * sizeof(type), message, status = 0 );      \
    if ( status ) {                                                          \
      PMMG_REALLOC( ptr, newsize, type, message );                           \
      (mesh)->memCur -= (long long)( oldsize * sizeof(type) );               \
      if ( oldsize < newsize )                                               \
        memset( ptr + oldsize, 0, (newSize - oldsize) * sizeof(type) );      \
    } else {                                                                 \
      MMG_MEM_ERRORMSG;                                                      \
    }                                                                        \
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

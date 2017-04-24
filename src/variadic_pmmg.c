/**
 * \file variadic_pmmg.c
 * \brief C API functions definitions for PARMMG library.
 * \copyright GNU Lesser General Public License.
 *
 * Variadic functions for PARMMG library.
 *
 */

#include "parmmg.h"

void PMMG_Init_parameters( PMMG_pParMesh parmesh )
{
  parmesh->memMax = _MMG5_memSize();
  /* If total physical memory is known, use 50% of it.
   * Otherwise try to use a default value of 800 Mo */
  if ( parmesh->memMax ) {
    parmesh->memMax = parmesh->memMax*50/100;
  } else {
    printf("  Maximum memory arbitrarily set to default value: %d Mo.\n",_MMG5_MEMMAX);
    parmesh->memMax = _MMG5_MEMMAX << 20;
  }
}

// return:
//   0 on error: Incorrect arguments passed
//   1 on success
int PMMG_Init_parMesh_var( va_list argptr )
{
  PMMG_pParMesh *parmesh = NULL;
  PMMG_pGrp     grp = NULL;
  int           typArg;
  int           parmeshCount = 0;
  MPI_Comm      comm = MPI_COMM_WORLD;

  while ( ( typArg = va_arg( argptr, int ) ) != PMMG_ARG_end ) {
    switch ( typArg ) {
    case( PMMG_ARG_ppParMesh ):
      parmesh = va_arg( argptr, PMMG_pParMesh* );
      ++parmeshCount;
      break;
//NIKOS TODO: if this option is available we have to also handle communicators creation/initialization
//    case(PMMG_ARG_ngroups):
//      ngrps = va_arg(argptr,int);
//      break;
    case( PMMG_ARG_comm ):
      comm = va_arg( argptr, MPI_Comm );
      break;
    default:
      fprintf(stderr, "  ## Error: PMMG_Init_parMesh:\n unexpected argument type: %d\n", typArg );
      fprintf(stderr, " Argument type must be one of the PMMG_ARG* preprocessor variable:"
                      " PMMG_ARG_ppParMesh, PMMG5_ARG_comm\n" );
      return ( 0 );
    }
  }

  if ( parmeshCount != 1 ) {
    fprintf( stderr, "  ## Error: PMMG_Init_parMesh:\n you have to provide exactly one"
                     " parmesh structure to be initialized and hold the mesh.\n" );
    return ( 0 );
  }

  /* ParMesh allocation */
  if ( *parmesh )
    _MMG5_SAFE_FREE( *parmesh );
  _MMG5_SAFE_CALLOC( *parmesh, 1, PMMG_ParMesh );
  PMMG_Init_parameters( *parmesh );

  _MMG5_ADD_MEM( parmesh, sizeof( PMMG_ParMesh ), "This will never be printed...", return( 0 ) );
#warning remove when the program will run
  (*parmesh)->ddebug = 1;

  /* MPI Data */
  (*parmesh)->comm = comm;
  MPI_Comm_rank( comm, &(*parmesh)->myrank );
  MPI_Comm_size( comm, &(*parmesh)->nprocs );

  /** Init Group */
  (*parmesh)->ngrp = 1;
  PMMG_CALLOC ( (*parmesh)->listgrp, (*parmesh)->ngrp, PMMG_Grp, parmesh, return( 0 ) );

  grp = &(*parmesh)->listgrp[0];
  grp->mesh = NULL;
  grp->met  = NULL;
  grp->disp = NULL;
  MMG3D_Init_mesh( MMG5_ARG_start,
                   MMG5_ARG_ppMesh, &grp->mesh,
                   MMG5_ARG_ppMet, &grp->met,
                   MMG5_ARG_end );
  return ( 1 );
}

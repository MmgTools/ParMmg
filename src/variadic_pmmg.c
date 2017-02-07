/**
 * \file variadic_pmmg.c
 * \brief C API functions definitions for PARMMG library.
 * \copyright GNU Lesser General Public License.
 *
 * Variadic functions for PARMMG library.
 *
 */

#include "libparmmg.h"

int _PMMG_Init_parMesh_var( va_list argptr ) {
  PMMG_pParMesh  *parmesh;
  PMMG_pGrp      grp;
  int            *argc;
  int            typArg;
  int            parmeshCount,rank,size;
  int            ngrps,k;
  char           ***argv;
  MPI_Comm       comm;

  parmesh      = NULL;
  comm         = MPI_COMM_WORLD;
  ngrps        = 1;
  parmeshCount = 0;

  while ( (typArg = va_arg(argptr,int)) != PMMG_ARG_end )
  {
    switch ( typArg )
    {
    case(PMMG_ARG_ppParMesh):
      parmesh = va_arg(argptr,PMMG_pParMesh*);
      ++parmeshCount;
      break;
    case(PMMG_ARG_ngroups):
      ngrps = va_arg(argptr,int);
      break;
    case(PMMG_ARG_comm):
      comm = va_arg(argptr,MPI_Comm);
      break;
    default:
      fprintf(stderr,"  ## Error: PMMG_Init_parMesh:\n"
              " unexpected argument type: %d\n",typArg);
      fprintf(stderr," Argument type must be one"
              " of the PMMG_ARG* preprocessor variable:"
              " PMMG_ARG_ppParMesh, MMG5_ARG_ngroups,"
              "  MMG5_ARG_comm\n");
      return(0);
    }
  }

  if ( parmeshCount !=1 ) {
    fprintf(stderr,"  ## Error: PMMG_Init_parMesh:\n"
            " you need to initialize the parmesh structure that"
            " will contain your mesh.\n");
    return(0);
  }

  /* ParMesh allocation */
  if ( *parmesh )  _MMG5_SAFE_FREE(*parmesh);
  _MMG5_SAFE_CALLOC(*parmesh,1,PMMG_ParMesh);

#warning to remove when the program will run
  (*parmesh)->ddebug = 1;

  /* MPI Data */
  (*parmesh)->comm   = comm;
  MPI_Comm_rank(comm,(*parmesh)->myrank);
  MPI_Comm_size((comm,(*parmesh)->nprocs);

  /** Init Group */
  (*parmesh)->ngrp = ngrps;
  _MMG5_SAFE_CALLOC((*parmesh)->listgrp,(*parmesh)->ngrp,PMMG_Grp);

  for ( k=0; k<(*parmesh)->ngrp; ++k ) {
    grp = &(*parmesh)->listgrp[k];

    grp->mesh = NULL;
    grp->sol  = NULL;

    MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&grp->mesh,
                    MMG5_ARG_ppMet,&grp->sol,
                    MMG5_ARG_end);

  }

  return 1;
}

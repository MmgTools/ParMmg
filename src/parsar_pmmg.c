#include "parmmg.h"

// Helper macro. Used only in this file. On success:
//   copies the contents of fromV[fromC] (which are argv[argc]) to toV[toC]
//   updates toC
#define ARGV_APPEND(parmesh,fromV,toV,fromC,toC,msg,on_failure)   do {       \
  PMMG_MALLOC(parmesh, toV[ toC ], strlen( fromV[ fromC ] ) + 1, char, msg,  \
              on_failure);                                                   \
  strncpy( toV[ toC ], fromV[ fromC ], strlen( fromV[ fromC ] ) + 1 );       \
  toV[ toC ][strlen( fromV[ fromC ] )] = '\0';                               \
  ++toC;                                                                     \
}while(0)

// Free custom argv allocations
static void
argv_cleanup( PMMG_pParMesh parmesh, char **mmgArgv, int mmgArgc, int argc )
{
  int i;
  for ( i = 0; i < mmgArgc; ++i )
    PMMG_DEL_MEM(parmesh, mmgArgv[i], strlen( mmgArgv[i] ), char, "Deallocating mmgargv[i]: " );
  PMMG_DEL_MEM(parmesh, mmgArgv, argc, char*, "Deallocating mmgargv: " );
}

static void
defaultValues( PMMG_pParMesh parmesh, const int rank )
{
  if ( rank == 0 ) {
    fprintf( stdout, "\n\n\tParMMG\nDefault parameter values:\n\n");
    fprintf( stdout,"# of remeshing iterations (-niter)  : 1\n");
    fprintf( stdout, "\n\n\tMMG");
    _MMG5_mmgDefaultValues( parmesh->listgrp[0].mesh );
  }
  PMMG_exit_and_free( parmesh, PMMG_SUCCESS );
}

static void
usage( PMMG_pParMesh parmesh, char * const progname )
{
  if ( parmesh->myrank == 0 ) {
    fprintf( stdout, "\n\n\tParMMG\nDefault parameter values:\n\n");
    fprintf( stdout, "-niter n  Number of remeshing iterations\n");
    fprintf( stdout, "\n\n\tMMG");
    MMG3D_usage( progname );
  }
  PMMG_exit_and_free( parmesh, PMMG_SUCCESS );
}

/**
* \param argc the argument count parameter from main
* \param argv the argument values parameter from main
* \param parmesh pointer to active parmesh
*
* \return PMMG_SUCCESS
*         PMMG_FAILURE
* PARSe ARguments from the command line function
* it works on top of the MMG3D_parsar function, ie:
*   arguments exclusively for PMMG_parsar are handled here and not added in the
*     mmgArgc/mmgArgv list.
*   arguments common with MMG3D_parsar are handled here and also added in the
*     mmgArgc/mmgArgv list and passed to MMG3D_parsar.
*   arguments that only affect MMG3D_parsar are simply added in the
*     mmgArgc/mmgArgv list and passed to MMG3D_parsar.
*/
int PMMG_parsar( int argc, char *argv[], PMMG_pParMesh parmesh )
{
  int i = 0;
  const long int million = 1024L * 1024L;
  int ret_val = PMMG_SUCCESS;

  int mmgArgc = 0;
  char** mmgArgv = NULL;

  for ( i = 1; i < argc; ++i )
    if ( !strcmp( argv[ i ],"-val" ) )
      defaultValues( parmesh, parmesh->myrank );
    else if ( ( !strcmp( argv[ i ],"-?" ) ) || ( !strcmp( argv[ i ],"-h" ) ) )
      usage( parmesh, argv[0] );

  // Create a new set of argc/argv variables adding only the the cl options that
  // mmg has to process
  // Overallocating as they are at most argc. Trying to avoid the overallocation
  // is not worth any effort, these are ~kb
  PMMG_MALLOC(parmesh, mmgArgv, argc, char*, " copy of argv for mmg: ",
              ret_val = PMMG_FAILURE; goto fail_mmgargv);

  // First argument is always argv[0] ie prog name
  i = 0;
  ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc, " mmgArgv[0] for mmg: ",
              ret_val = PMMG_FAILURE; goto fail_proc);

  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch( argv[i][1] ) {
      case 'm':  /* memory */
        if ( ++i < argc && isdigit( argv[i][0] ) ) {
          if ( ( atoi(argv[ i ]) > _MMG5_memSize() ) || ( atoi(argv[ i ]) < 0 ) )
            fprintf( stderr,
                     "Erroneous mem size requested, using default: %lld\n",
                     parmesh->memGloMax );
          else
            PMMG_PMesh_SetMemGloMax( parmesh, atoi( argv[i] ) * million);
          if ( 1 == MMG3D_Set_iparameter( parmesh->listgrp[0].mesh,
                                          parmesh->listgrp[0].met,
                                          MMG3D_IPARAM_mem,
                                          atoi( argv[i] ) ) ) {
            ret_val = PMMG_FAILURE;
            goto fail_proc;
          }
        } else {
          fprintf( stderr, "Missing argument option %c\n", argv[i-1][1] );
          usage( parmesh, argv[0] );
        }
      break;

      case 'n':  // number of adaptation iterations
        if ( ( 0 == strncmp( argv[i], "-niter", 5 ) ) && ( ( i + 1 ) < argc ) ) {
          ++i;
          if ( isdigit( argv[i][0] ) && ( atoi( argv[i] ) > 0 ) ) {
            parmesh->niter = atoi( argv[i] );
          } else {
            parmesh->niter = 1;
            fprintf( stderr,
                     "Erroneous adaptation iterations requested, using default: %d\n",
                     parmesh->niter );
          }
        } else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = PMMG_FAILURE; goto fail_proc );
        }
      break;
      default:
        ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                    " adding to mmgArgv for mmg: ",
                    ret_val = PMMG_FAILURE; goto fail_proc);

      break;
      }
    } else {
      ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                  " adding to mmgArgv for mmg: ",
                  ret_val = PMMG_FAILURE; goto fail_proc);
    }
    ++i;
  }

  // parmmg finished parsing arguments, the rest will e handled by mmg3d
  if ( 1 != MMG3D_parsar( mmgArgc, mmgArgv,
                          parmesh->listgrp[0].mesh,
                          parmesh->listgrp[0].met ) ) {
    ret_val = PMMG_FAILURE;
    goto fail_proc;
  }

fail_proc:
  argv_cleanup( parmesh, mmgArgv, mmgArgc, argc );
fail_mmgargv:
  return ret_val;
}

#undef ARGV_APPEND

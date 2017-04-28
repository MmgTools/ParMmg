#include "parmmg.h"

// Helper macro used only in this file:
//   copies the contents of fromV[fromC] (which are argv[argc]) to
//   toV[toC] and update toC
#define APPEND_ARGV(parmesh,fromV,toV,fromC,toC,message) do {                  \
    PMMG_MALLOC( parmesh, toV[(toC)],                                          \
                 (strlen( fromV[(fromC)] ) + 1 ) * sizeof(char), message );    \
    strncpy( toV[toC], fromV[fromC], strlen( fromV[fromC] ) + 1 );             \
    toV[toC][strlen( fromV[fromC] )]='\0';                                     \
    ++toC;                                                                     \
  }while(0)


static void defaultValues( const PMMG_pParMesh parmesh )
{
  fprintf( stdout, "\n\n\tParMMG\nDefault parameter values:\n\n");

  fprintf( stdout, "\n\n\tMMG");
  _MMG5_mmgDefaultValues( parmesh->listgrp[0].mesh );

  exit( EXIT_SUCCESS );
}

static void usage( char *progname )
{
  fprintf( stdout, "\n\n\tParMMG\nDefault parameter values:\n\n");

  fprintf( stdout, "\n\n\tMMG");
  MMG3D_usage( progname );
}

/**
* \param argc the argument count parameter from main
* \param argv the argument values parameter from main
* \param parmesh pointer to active parmesh
*
* \return 0 on success
*         1 on failure
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

  int mmgArgc = 0;
  char** mmgArgv = NULL;

  for ( i = 1; i < argc; ++i )
    if ( !strcmp( argv[ i ],"-val" ) )
      defaultValues( parmesh );
    else if ( ( !strcmp( argv[ i ],"-?" ) ) || ( !strcmp( argv[ i ],"-h" ) ) )
      usage( argv[0] );

  // create a new set of argc/argv variables adding only the the cl options
  // intended for mmg
  // overallocating as there will never be more than argc.
  // not worth the trouble to try and decrease it as it is less than a kb anyway
  PMMG_MALLOC(parmesh,mmgArgv,argc * sizeof(char*), " copy of argv for mmg: ");

  // First argument is always argv[0] ie prog name
  i = 0;
  APPEND_ARGV(parmesh,argv,mmgArgv,i,mmgArgc," mmgArgv[0] for mmg: ");

  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
      case 'm':  /* memory */
        if ( ++i < argc && isdigit( argv[i][0] ) ) {
          if ( (atoi(argv[ i ]) > _MMG5_memSize()) || (atoi(argv[ i ]) < 0) )
            fprintf( stderr,
                     "Erroneous mem size requested, using default: %lld\n",
                     parmesh->memMax );
          else
            parmesh->memMax = atoi(argv[i]);
        } else {
          fprintf( stderr, "Missing argument option %c\n", argv[i-1][1] );
          usage( argv[0] );
        }
      break;

      default:
        APPEND_ARGV(parmesh,argv,mmgArgv,i,mmgArgc,
                    " adding to mmgArgv for mmg: ");
      break;
      }
    } else {
      APPEND_ARGV(parmesh,argv,mmgArgv,i,mmgArgc,
                  " adding to mmgArgv for mmg: ");
    }
    ++i;
  }

  // parmmg finished parsing arguments, the rest will e handled by mmg3d
  MMG3D_parsar( mmgArgc, mmgArgv, parmesh->listgrp[0].mesh,
                parmesh->listgrp[0].met );

  // Free allocations and return
  for ( i = 0; i < mmgArgc; ++i )
    PMMG_FREE(parmesh,mmgArgv[i],strlen(mmgArgv[i]),"Deallocating mmgargv[i]: ");
  PMMG_FREE(parmesh,mmgArgv,argc * sizeof(char*),"Deallocating mmgargv: ");

  return 0;
}

#undef APPEND_ARGV

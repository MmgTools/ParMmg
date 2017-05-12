#include "parmmg.h"

// Helper macro used only in this file:
//   copies the contents of fromV[fromC] (which are argv[argc]) to
//   toV[toC] and update toC
//+++++NIKOS TODO: I DO NOT LIKE THIS APPROACH USING THIS MACRO, PLEASE REVISE
//+++++NIKOS TODO: returning without freeing is a memory leak here
#define APPEND_ARGV(parmesh,fromV,toV,fromC,toC,message) do {                  \
    PMMG_MALLOC( parmesh, toV[(toC)],                                          \
                 (strlen( fromV[(fromC)] ) + 1 ) * sizeof(char), message );    \
    strncpy( toV[toC], fromV[fromC], strlen( fromV[fromC] ) + 1 );             \
    toV[toC][strlen( fromV[fromC] )]='\0';                                     \
    ++toC;                                                                     \
  }while(0)


static void defaultValues( const MMG5_pMesh mesh, const int rank )
{
  if ( rank == 0 ) {
    fprintf( stdout, "\n\n\tParMMG\nDefault parameter values:\n\n");
    fprintf( stdout,"# of remeshing iterations (-niter)  : 1\n");
    fprintf( stdout, "\n\n\tMMG");
    _MMG5_mmgDefaultValues( mesh );
  }
  MPI_Finalize();
  exit( EXIT_SUCCESS );
}

static void usage( char * const progname, const int rank )
{
  if ( rank == 0 ) {
    fprintf( stdout, "\n\n\tParMMG\nDefault parameter values:\n\n");
    fprintf( stdout, "-niter n  Number of remeshing iterations\n");
    fprintf( stdout, "\n\n\tMMG");
    MMG3D_usage( progname );
  }
  MPI_Finalize();
  exit( EXIT_SUCCESS );
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
  const long int million = 1048576L;

  int mmgArgc = 0;
  char** mmgArgv = NULL;

  for ( i = 1; i < argc; ++i )
    if ( !strcmp( argv[ i ],"-val" ) )
      defaultValues( parmesh->listgrp[0].mesh, parmesh->myrank );
    else if ( ( !strcmp( argv[ i ],"-?" ) ) || ( !strcmp( argv[ i ],"-h" ) ) )
      usage( argv[0], parmesh->myrank );

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
            parmesh->memMax = atoi(argv[i]) * million;
          if ( !MMG3D_Set_iparameter( parmesh->listgrp[0].mesh,
            parmesh->listgrp[0].met, MMG3D_IPARAM_mem,atoi(argv[i])) )
            return 1; //!!!!!NIKOS TODO
        } else {
          fprintf( stderr, "Missing argument option %c\n", argv[i-1][1] );
          usage( argv[0], parmesh->myrank );
        }
      break;

      case 'n':  // number of adaptation iterations
        if ( (0 == strncmp( argv[i], "-niter", 5 )) && ((i+1) < argc) ) {
          ++i;
          if ( isdigit( argv[i][0] ) && (atoi( argv[i] ) > 0) ) {
            parmesh->niter = atoi( argv[i] );
            printf( "\t\t\t\t +++++++++++++ASD : adding: %d \n\n\n\n",  atoi( argv[i] ) );
          } else {
            parmesh->niter = 1;
            fprintf( stderr,
                     "Erroneous adaptation iterations requested, using default: %d\n",
                     parmesh->niter );
          }
        } else {
          printf( "\t\t\t\t +++++++++++++ASD : den etairiaksene: %s \n\n\n\n",  argv[i] );
          APPEND_ARGV(parmesh,argv,mmgArgv,i,mmgArgc,
                      " adding to mmgArgv for mmg: ");
        }
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

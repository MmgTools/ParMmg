#include "parmmg.h"

static void defaultValues( const PMMG_pParMesh parmesh )
{
  fprintf( stdout, "\n\n\tParMMG\nDefault parameter values:\n\n");
  fprintf( stdout, "Number of groups = 1\n");

  fprintf( stdout, "\n\n\tMMG");
  _MMG5_mmgDefaultValues( parmesh->listgrp[0].mesh );

  exit( EXIT_SUCCESS );
}

int PMMG_parsar( int argc, char *argv[], PMMG_pParMesh parmesh )
{
  int i = 0;

  int mmgAlloced = argc;
  int mmgArgc = 0;
  char** mmgArgv = NULL;

  for ( i = 1; i < argc; ++i )
    if ( !strcmp( argv[ i ],"-val" ) )
      defaultValues( parmesh );

//NIKOS DOING:
//     handle the -h option by printing the pmmg specifics and then the mmg -h option
//     handle parmmg specific options and do not add them in the mmgargc/argv list
//     handle common parmmg/mmg options and add them to the mmgArgc/argv list
//  _MMG5_SAFE_CALLOC( mmgArgv, allocated, char* );

//  if ( mmgArgc >= allocated ) {
//    allocated *= 2;
//    _MMG5_SAFE_REALLOC( mmgArgv, allocated, char*, "increasing argv array" );
//  }
//  strncpy( mmgArgv[mmgArgc], argv[mmgArgc], strlen( argv[mmgArgc] ) + 1 );
//  ++mmgArgc;
//     pass the created list to mmg3d_parsar

//  _MMG5_SAFE_CALLOC( mmgArgv[0], strlen( argv[0] ) + 1, char );
//  strncpy( mmgArgv[0], argv[0], strlen( argv[0] ) + 1 );
//  mmgArgv[0][strlen( argv[0] )]='\0';

  //MMG3D_parsar( mmgArgc, mmgArgv, parmesh->listgrp[0].mesh, parmesh->listgrp[0].met );
  MMG3D_parsar( argc, argv, parmesh->listgrp[0].mesh, parmesh->listgrp[0].met );

//  _MMG5_SAFE_FREE( mmgArgv );
  return EXIT_SUCCESS;
}

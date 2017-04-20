#include "parmmg.h"

void PMMG_defaultValues( PMMG_pParMesh parmesh )
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
  int count = 0;
  int mmgArgc = 1;
  char ** mmgArgv = NULL;
  
  for ( i = 1; i < argc; ++i )
    if ( !strcmp( argv[ i ],"-val" ) )
      PMMG_defaultValues( parmesh );

//NIKOS DOING: Now, handle the -h option by printing the pmmg specifics and then the mmg as usual 

  _MMG5_SAFE_CALLOC( mmgArgv, argc, char* );
  _MMG5_SAFE_CALLOC( mmgArgv[0], strlen( argv[0] ) + 1, char );
  strncpy( mmgArgv[0], argv[0], strlen( argv[0] ) + 1 );
  mmgArgv[0][strlen( argv[0] )]='\0';
//NIKOS DOING: Parse the arguments: 
//               if matching argument is found, do what you have to do
//               else add it to mmgArgv/mmgArgc

  MMG3D_parsar( mmgArgc, mmgArgv, parmesh->listgrp[0].mesh, parmesh->listgrp[0].sol );

  return EXIT_SUCCESS;
}

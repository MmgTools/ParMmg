/* =============================================================================
**  This file is part of the parmmg software package for parallel tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux, 2017-
**
**  parmmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  parmmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with parmmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the parmmg distribution only if you accept them.
** =============================================================================
*/

/**
 * \file libparmmg_tools.c
 * \brief C API functions definitions for PARMMG library.
 * \copyright GNU Lesser General Public License.
 *
 * C API for PARMMG library.
 *
 */
#include "parmmg.h"

/*! Helper macro used only in this file: copies the contents of fromV[fromC]
 *  to toV[toC] updates toC */
#define ARGV_APPEND(parmesh,fromV,toV,fromC,toC,msg,on_failure)   do {  \
    PMMG_MALLOC(parmesh, toV[ toC ], strlen( fromV[ fromC ] ) + 1, char, msg, \
                on_failure);                                            \
    strncpy( toV[ toC ], fromV[ fromC ], strlen( fromV[ fromC ] ) + 1 ); \
    toV[ toC ][strlen( fromV[ fromC ] )] = '\0';                        \
    ++toC;                                                              \
  }while(0)

/**
 * \param parmesh pointer to pmmg structure
 * \param mmgArgv pointer to argv like buffer
 * \param mmgArgc pointer to argc like buffer
 * \param argc    actual argc value
 *
 * Free the allocations of the custom created argc/argv wrapper that is passed
 * to mmg to parse the command line options
 */
static void
PMMG_argv_cleanup( PMMG_pParMesh parmesh, char **mmgArgv, int mmgArgc, int argc )
{
  int i;
  for ( i = 0; i < mmgArgc; ++i )
    PMMG_DEL_MEM(parmesh, mmgArgv[i],char, "Deallocating mmgargv[i]: " );
  PMMG_DEL_MEM(parmesh, mmgArgv,char*, "Deallocating mmgargv: " );
}

int PMMG_defaultValues( PMMG_pParMesh parmesh )
{
  int ier = 1;

  if ( !parmesh->myrank ) {
    fprintf(stdout,"\n\n");
    fprintf(stdout,"  --- ParMMG ---\n");
    fprintf(stdout,"default parameter values:\n");
    fprintf(stdout,"\n** Generic options\n");
    fprintf(stdout,"verbosity                 (-v)      : %d\n",
            parmesh->info.imprim);

    fprintf(stdout,"maximal memory size       (-m)      : %zu MB\n",
            parmesh->memGloMax/MMG5_MILLION);
    fprintf(stdout,"\n** Parameters\n");
    fprintf( stdout,"# of remeshing iterations (-niter)        : %d\n",parmesh->niter);
    fprintf( stdout,"loadbalancing_mode (not yet customizable) : PMMG_LOADBALANCING_metis\n");
    fprintf( stdout,"target mesh size for Mmg (-mesh-size) : %d\n",abs(PMMG_REMESHER_TARGET_MESH_SIZE));
    fprintf( stdout,"ratio: # meshes / # metis super nodes (-metis-ratio) : %d\n",
             abs(PMMG_RATIO_MMG_METIS) );
    fprintf( stdout,"# of layers for interface displacement (-nlayers) : %d\n",PMMG_MVIFCS_NLAYERS);

#ifdef USE_SCOTCH
    fprintf(stdout,"SCOTCH renumbering                  : enabled\n");
#else
    fprintf(stdout,"SCOTCH renumbering                  : disabled\n");
#endif

    if ( parmesh->listgrp[0].mesh ) {
      fprintf(stdout,"\n  --- MMG ---");
      if ( !MMG3D_defaultValues( parmesh->listgrp[0].mesh ) ) {
        ier = 0;
      }
    }
  }

  return ier;
}

int PMMG_usage( PMMG_pParMesh parmesh, char * const prog )
{
  if ( !parmesh->myrank ) {
    fprintf(stdout,"\nUsage: %s [-v [n]] [opts..] filein [fileout]\n",
            prog);

    fprintf(stdout,"\n** Generic options :\n");
    fprintf(stdout,"-h         Print this message\n");
    fprintf(stdout,"-v [n]     Tune ParMmg level of verbosity, [-10..10]\n");
    fprintf(stdout,"-mmg-v [n] Tune Mmg level of verbosity, [-10..10]\n");
    fprintf(stdout,"-m [n]     Set maximal memory size to n Mbytes\n");
    fprintf(stdout,"-d         Turn on debug mode for ParMmg\n");
    fprintf(stdout,"-mmg-d     Turn on debug mode for Mmg\n");
    fprintf(stdout,"-val       Print the default parameters values\n");
    //fprintf(stdout,"-default  Save a local parameters file for default parameters"
    //        " values\n");

    fprintf(stdout,"\n**  File specifications\n");
    fprintf(stdout,"-in  file  input triangulation\n");
    fprintf(stdout,"-out file  output triangulation\n");
    fprintf(stdout,"-sol file  load solution or metric file\n");

    fprintf(stdout,"\n**  Parameters\n");
    fprintf(stdout,"-niter        val  number of remeshing iterations\n");
    fprintf(stdout,"-mesh-size    val  target mesh size for the remesher\n");
    fprintf(stdout,"-metis-ratio  val  number of metis super nodes per mesh\n");
    fprintf(stdout,"-nlayers      val  number of layers for interface displacement\n");

    //fprintf(stdout,"-ar     val  angle detection\n");
    //fprintf(stdout,"-nr          no angle detection\n");
    fprintf(stdout,"-hmin   val  minimal mesh size\n");
    fprintf(stdout,"-hmax   val  maximal mesh size\n");
    fprintf(stdout,"-hsiz   val  constant mesh size\n");
    // fprintf(stdout,"-hausd  val  control Hausdorff distance\n");
    fprintf(stdout,"-hgrad  val  control gradation\n");
    // fprintf(stdout,"-ls     val  create mesh of isovalue val (0 if no argument provided)\n");
    fprintf(stdout,"-A           enable anisotropy (without metric file).\n");
    // fprintf(stdout,"-opnbdy      preserve input triangles at the interface of"
    //        " two domains of the same reference.\n");

#ifdef USE_ELAS
    // fprintf(stdout,"-lag [0/1/2] Lagrangian mesh displacement according to mode 0/1/2\n");
#endif
#ifndef PATTERN
    fprintf(stdout,"-octree val  Specify the max number of points per octree cell \n");
#endif
#ifdef USE_SCOTCH
    fprintf(stdout,"-rn [n]      Turn on or off the renumbering using SCOTCH [1/0] \n");
#endif
    fprintf(stdout,"\n");

    fprintf(stdout,"-nofem       do not force Mmg to create a finite element mesh \n");
    fprintf(stdout,"-optim       mesh optimization\n");
    fprintf(stdout,"-optimLES    strong mesh optimization for LES computations\n");
    fprintf(stdout,"-noinsert    no point insertion/deletion \n");
    fprintf(stdout,"-noswap      no edge or face flipping\n");
    fprintf(stdout,"-nomove      no point relocation\n");
    fprintf(stdout,"-nosurf      no surface modifications\n");
    fprintf(stdout,"\n\n");

  }

  return 1;
}

int PMMG_parsar( int argc, char *argv[], PMMG_pParMesh parmesh )
{
  int        val,i  = 0;
  int        ret_val = 1;
  int        mmgArgc = 0;
  char**     mmgArgv = NULL;

  assert ( parmesh->ngrp == 1 && "distributed input not yet implemented" );

  /** Parse arguments specific to parMmg then add to mmgArgv the mmg arguments
   * and call the mmg3d parser. */
  for ( i = 1; i < argc; ++i ) {
    if ( !strcmp( argv[ i ],"-val" ) ) {
      RUN_ON_ROOT_AND_BCAST( PMMG_defaultValues(parmesh),0,
                             parmesh->myrank,ret_val=0; goto fail_mmgargv);
      ret_val = 0;
      goto fail_mmgargv;
    }
    else if ( ( !strcmp( argv[ i ],"-?" ) ) || ( !strcmp( argv[ i ],"-h" ) ) ) {
      RUN_ON_ROOT_AND_BCAST( PMMG_usage(parmesh, argv[0]),0,
                             parmesh->myrank,ret_val=0; goto fail_mmgargv);
      ret_val = 0;
      goto fail_mmgargv;
    }
  }

  /* Create a new set of argc/argv variables adding only the the cl options that
     mmg has to process
     Overallocating as they are at most argc. Trying to avoid the overallocation
     is not worth any effort, these are ~kb */
  PMMG_MALLOC(parmesh, mmgArgv, argc, char*, " copy of argv for mmg: ",
              ret_val = 0; goto fail_mmgargv);

  /* First argument is always argv[0] ie prog name */
  i = 0;
  ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc, " mmgArgv[0] for mmg: ",
              ret_val = 0; goto fail_proc);

  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch( argv[i][1] ) {
      case 'm':
        if ( !strcmp(argv[i],"-mmg-v") ) {

          /* Mmg verbosity */
          if ( ++i < argc ) {
            if ( isdigit(argv[i][0]) ||
                 (argv[i][0]=='-' && isdigit(argv[i][1])) ) {
              val = atoi(argv[i]);

              if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_mmgVerbose,val) ) {
                ret_val = 0;
                goto fail_proc;
              }
            }
            else {
              i--;
            }
          }
          else {
            fprintf( stderr, "\nMissing argument option %c\n", argv[i-1][1] );
             ret_val = 0;
            goto fail_proc;
          }
        }
        else if ( !strcmp(argv[i],"-mesh-size") ) {

          /* Remesher target mesh size */
          if ( ++i < argc ) {
            if ( isdigit(argv[i][0]) ) {
              val = atoi(argv[i]);

              if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_meshSize,val) ) {
                ret_val = 0;
                goto fail_proc;
              }
            }
            else {
              fprintf( stderr, "\nMissing argument option %c\n", argv[i-1][1] );
              ret_val = 0;
              goto fail_proc;
            }
          }
          else {
            fprintf( stderr, "\nMissing argument option %c\n", argv[i-1][1] );
            ret_val = 0;
            goto fail_proc;
          }
        }
        else if ( !strcmp(argv[i],"-metis-ratio") ) {

          /* Number of metis super nodes per mesh */
          if ( ++i < argc ) {
            if ( isdigit(argv[i][0]) ) {
              val = atoi(argv[i]);

              if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_metisRatio,val) ) {
                ret_val = 0;
                goto fail_proc;
              }
            }
            else {
              fprintf( stderr, "\nMissing argument option %c\n", argv[i-1][1] );
              ret_val = 0;
              goto fail_proc;
            }
          }
          else {
            fprintf( stderr, "\nMissing argument option %c\n", argv[i-1][1] );
            ret_val = 0;
            goto fail_proc;
          }
        }
        else if ( !strcmp(argv[i],"-mmg-d") ) {
          if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_mmgDebug,val) ) {
            ret_val = 0;
            goto fail_proc;
          }
        }
        else {
          /* memory */
          if ( ++i < argc && isdigit( argv[i][0] ) ) {
            if ( ( atoi(argv[ i ]) > MMG5_memSize() ) || ( atoi(argv[ i ]) < 0 ) ) {
              fprintf( stderr,
                       "\nErroneous mem size requested (%s)\n",argv[i] );
              ret_val = 0;
              goto fail_proc;
            }
            else {
              parmesh->info.mem = atoi( argv[i] );
              PMMG_parmesh_SetMemGloMax( parmesh );
            }
            PMMG_parmesh_SetMemMax( parmesh, 20 );
          } else {
            fprintf( stderr, "\nMissing argument option %c\n", argv[i-1][1] );
            ret_val = 0;
            goto fail_proc;
          }
        }
        break;

      case 'n':  /* number of adaptation iterations */
        if ( ( 0 == strncmp( argv[i], "-niter", 5 ) ) && ( ( i + 1 ) < argc ) ) {
          ++i;
          if ( isdigit( argv[i][0] ) && ( atoi( argv[i] ) > 0 ) ) {
            parmesh->niter = atoi( argv[i] );
          } else {
            parmesh->niter = PMMG_NITER;
            fprintf( stderr,
                     "\nWrong number of adaptation iterations (%s).\n",argv[i]);

            ret_val = 0;
            goto fail_proc;
          }
        } else if ( ( 0 == strncmp( argv[i], "-nlayers", 5 ) ) && ( ( i + 1 ) < argc ) ) {
          ++i;
          if ( isdigit( argv[i][0] ) && ( atoi( argv[i] ) > 0 ) ) {
            parmesh->info.ifc_layers = atoi( argv[i] );
          } else {
            parmesh->info.ifc_layers = PMMG_MVIFCS_NLAYERS;
            fprintf( stderr,
                     "\nWrong number of layers for interface displacement (%s).\n",argv[i]);

            ret_val = 0;
            goto fail_proc;
          }
        }else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;

      case 'd':  /* debug */
        if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_debug,1) )  {
          ret_val = 0;
          goto fail_proc;
        }
        break;
#ifdef USE_SCOTCH
      case 'r':
        if ( !strcmp(argv[i],"-rn") ) {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;
#endif

      case 'v':  /* verbosity */
        if ( ++i < argc ) {
          if ( isdigit(argv[i][0]) ||
               (argv[i][0]=='-' && isdigit(argv[i][1])) ) {
            if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_verbose,atoi(argv[i])) ) {
              ret_val = 0;
              goto fail_proc;
            }
          }
          else
            i--;
        }
        else {
          fprintf(stderr,"\nMissing argument option %c\n",argv[i-1][1]);
           ret_val = 0;
          goto fail_proc;
        }
        break;

      default:
        ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                    " adding to mmgArgv for mmg: ",
                    ret_val = 0; goto fail_proc);

        break;
      }
    } else {
      ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                  " adding to mmgArgv for mmg: ",
                  ret_val = 0; goto fail_proc);
    }
    ++i;
  }

  // parmmg finished parsing arguments, the rest will be handled by mmg3d
  if ( 1 != MMG3D_parsar( mmgArgc, mmgArgv,
                          parmesh->listgrp[0].mesh,
                          parmesh->listgrp[0].met,
                          NULL ) ) {
    ret_val = 0;
    goto fail_proc;
  }
  parmesh->info.fem = parmesh->listgrp[0].mesh->info.fem;

fail_proc:
  PMMG_argv_cleanup( parmesh, mmgArgv, mmgArgc, argc );
fail_mmgargv:
  return ret_val;
}

#undef ARGV_APPEND


int PMMG_parsop ( PMMG_pParMesh parmesh )
{
  MMG5_pMesh mesh;
  int        ier;

  assert ( parmesh->ngrp == 1 && "distributed input not yet implemented" );
  mesh = parmesh->listgrp[0].mesh;

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  mesh->info.imprim = parmesh->info.imprim;

  ier = MMG3D_parsop(mesh,parmesh->listgrp[0].met);

  /* Restore the mesh verbosity */
  mesh->info.imprim = parmesh->info.mmg_imprim;

  return ier;
}

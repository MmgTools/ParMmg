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
    memcpy( toV[ toC ], fromV[ fromC ], (strlen( fromV[ fromC ] ) + 1)*sizeof(char) ); \
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
    fprintf( stdout,"repartitioning mode                       : PMMG_REDISTRIBUTION_ifc_displacement\n");
//    fprintf( stdout,"loadbalancing_mode (not yet customizable) : PMMG_LOADBALANCING_metis\n");
//    fprintf( stdout,"target mesh size for Mmg (-mesh-size) : %d\n",abs(PMMG_REMESHER_TARGET_MESH_SIZE));
//    fprintf( stdout,"ratio: # meshes / # metis super nodes (-metis-ratio) : %d\n",abs(PMMG_RATIO_MMG_METIS) );
    fprintf( stdout,"# of layers for interface displacement (-nlayers) : %d\n",PMMG_MVIFCS_NLAYERS);
    fprintf( stdout,"allowed imbalance between current and desired groups size (-groups-ratio) : %f\n",PMMG_GRPS_RATIO);

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
    fprintf(stdout,"-in    file  input triangulation\n");
    fprintf(stdout,"-out   file  output triangulation\n");
    fprintf(stdout,"-sol   file  load level-set, displacement or metric file\n");
    fprintf(stdout,"-field file  load sol field to interpolate from init onto final mesh\n");
    fprintf(stdout,"-noout       do not write output triangulation\n");

    fprintf(stdout,"\n**  Parameters\n");
    fprintf(stdout,"-niter        val  number of remeshing iterations\n");
    fprintf(stdout,"-mesh-size    val  target mesh size for the remesher\n");
    fprintf(stdout,"-metis-ratio  val  number of metis super nodes per mesh\n");
    fprintf(stdout,"-nlayers      val  number of layers for interface displacement\n");
    fprintf(stdout,"-groups-ratio val  allowed imbalance between current and desired groups size\n");
    fprintf(stdout,"-nobalance         switch off load balancing of the output mesh\n");

    //fprintf(stdout,"-ar     val  angle detection\n");
    //fprintf(stdout,"-nr          no angle detection\n");
    fprintf(stdout,"-hmin         val  minimal mesh size\n");
    fprintf(stdout,"-hmax         val  maximal mesh size\n");
    fprintf(stdout,"-hsiz         val  constant mesh size\n");
    // fprintf(stdout,"-hausd  val  control Hausdorff distance\n");
    fprintf(stdout,"-hgrad        val  control gradation\n");
    fprintf(stdout,"-hgradreq     val  control gradation from required entities\n");
    // fprintf(stdout,"-ls     val  create mesh of isovalue val (0 if no argument provided)\n");
    fprintf(stdout,"-A                 enable anisotropy (without metric file).\n");
    // fprintf(stdout,"-opnbdy      preserve input triangles at the interface of"
    //        " two domains of the same reference.\n");

#ifdef USE_ELAS
    // fprintf(stdout,"-lag [0/1/2] Lagrangian mesh displacement according to mode 0/1/2\n");
#endif
#ifndef PATTERN
    fprintf(stdout,"-octree       val  Specify the max number of points per octree cell \n");
#endif
#ifdef USE_SCOTCH
    fprintf(stdout,"-rn [n]            Turn on or off the renumbering using SCOTCH [1/0] \n");
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

  assert ( parmesh->ngrp == 1 && "Not available for more than 1 group per proc.\n");

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
      case 'c':
        if ( !strcmp(argv[i],"-centralized-output") ) {
          /* force centralized output: only relevant using medit distributed
           * input or library call */
          if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_distributedOutput,0) )  {
            ret_val = 0;
            goto fail_proc;
          }
        }
        else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;
      case 'f':
        if ( !strcmp(argv[i],"-field") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( ! PMMG_Set_inputSolsName(parmesh,argv[i]) ) {
              RUN_ON_ROOT_AND_BCAST( PMMG_usage(parmesh, argv[0]),0,
                                     parmesh->myrank,ret_val=0; goto fail_mmgargv);
              ret_val = 0;
              goto fail_mmgargv;
            }
          }
          else {
            RUN_ON_ROOT_AND_BCAST( PMMG_usage(parmesh, argv[0]),0,
                                   parmesh->myrank,ret_val=0; goto fail_mmgargv);
            ret_val = 0;
            goto fail_mmgargv;
          }
        }
        else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;
      case 'h':
        if ( !strcmp(argv[i],"-hmin") && ++i < argc ) {
          if ( !PMMG_Set_dparameter(parmesh,PMMG_DPARAM_hmin,atof(argv[i])) ) {
            ret_val = 0;
            goto fail_proc;
          }
        }
        else if ( !strcmp(argv[i],"-hmax") && ++i < argc ) {
          if ( !PMMG_Set_dparameter(parmesh,PMMG_DPARAM_hmax,atof(argv[i])) ) {
            ret_val = 0;
            goto fail_proc;
          }
        }
        else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;
      case 'g':
        if ( !strcmp(argv[i],"-groups-ratio") ) {

          if ( ++i < argc ) {
            if ( isdigit(argv[i][0]) ) {

              if ( !PMMG_Set_dparameter(parmesh,PMMG_DPARAM_groupsRatio,atof(argv[i])) ) {
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
        else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;

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
        else if ( !strcmp(argv[i],"-m") ) {
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
        else {
          /* else : what happens with -met option... to treat */
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;

      case 'n':  /* number of adaptation iterations */
        if ( ( 0 == strncmp( argv[i], "-niter", 5 ) ) && ( ( i + 1 ) < argc ) ) {
          ++i;
          if ( isdigit( argv[i][0] ) && ( atoi( argv[i] ) >= 0 ) ) {
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
        } else if ( 0 == strncmp( argv[i], "-nobalance", 9 ) ) {
          parmesh->info.nobalancing = MMG5_ON;
        } else if ( 0 == strncmp( argv[i], "-nofem", 5 ) ) {
          if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_nofem,1) )  {
            ret_val = 0;
            goto fail_proc;
          }
        } else if ( 0 == strncmp( argv[i], "-noout", 5 ) ) {
          parmesh->info.fmtout = PMMG_UNSET;
        } else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;

      case 'd':
        if ( !strcmp(argv[i],"-distributed-output") ) {
          /* force distributed output: only relevant using medit centralized
           * input or library call */
          if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_distributedOutput,1) )  {
            ret_val = 0;
            goto fail_proc;
          }
        }
        else if ( !strcmp(argv[i],"-d") ) {
          /* debug */
          if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_debug,1) )  {
            ret_val = 0;
            goto fail_proc;
          }
        }
        else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
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

      case 's':
        if ( 0 == strncmp( argv[i], "-surf", 4 ) ) {
          parmesh->listgrp[0].mesh->info.nosurf = 0;
        }
        else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;
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

  /* Store mesh names into the parmesh if needed */
  if ( !parmesh->meshin ) {
    assert ( parmesh->listgrp[0].mesh->namein );
    PMMG_Set_name(parmesh,&parmesh->meshin,
                  parmesh->listgrp[0].mesh->namein,"mesh.mesh");
  }
  if ( !parmesh->meshout ) {
    assert ( parmesh->listgrp[0].mesh->nameout );
    PMMG_Set_name(parmesh,&parmesh->meshout,
                  parmesh->listgrp[0].mesh->nameout,"mesh.o.mesh");
  }
  if ( (!parmesh->metin) && parmesh->listgrp[0].met && parmesh->listgrp[0].met->namein ) {
    PMMG_Set_name(parmesh,&parmesh->metin,
                  parmesh->listgrp[0].met->namein,"mesh.sol");
  }
  if ( (!parmesh->metout) && parmesh->listgrp[0].met && parmesh->listgrp[0].met->nameout ) {
    PMMG_Set_name(parmesh,&parmesh->metout,
                  parmesh->listgrp[0].met->nameout,"mesh.o.sol");
  }
  if ( (!parmesh->lsin) && parmesh->listgrp[0].ls && parmesh->listgrp[0].ls->namein ) {
    PMMG_Set_name(parmesh,&parmesh->lsin,
                  parmesh->listgrp[0].ls->namein,"mesh.sol");
  }
  if ( (!parmesh->dispin) && parmesh->listgrp[0].disp && parmesh->listgrp[0].disp->namein ) {
    PMMG_Set_name(parmesh,&parmesh->dispin,
                  parmesh->listgrp[0].disp->namein,"mesh.sol");
  }

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

void PMMG_setfunc( PMMG_pParMesh parmesh ) {
  MMG5_pSol met = parmesh->listgrp[0].met;

  if( met && met->size == 6 ) {

    PMMG_interp4bar = PMMG_interp4bar_ani;
    PMMG_interp3bar = PMMG_interp3bar_ani;
    PMMG_interp2bar = PMMG_interp2bar_ani;

  } else {

    PMMG_interp4bar = PMMG_interp4bar_iso;
    PMMG_interp3bar = PMMG_interp3bar_iso;
    PMMG_interp2bar = PMMG_interp2bar_iso;

  }

}

/**
 * \param fp file pointer
 * \param bin 1 if file is at Medit binary format
 * \param linesep line separator
 *
 * \return 1 if success, 0 if fail
 *
 * Search the End keyword in an Medit file and position the file pointer before
 * this keyword.
 *
 */
static inline
int PMMG_search_filePosition ( FILE *fp,int bin,char linesep ) {
  int  ier;
  long pos;
  char c;

  /* Return if file pointer is not valid */
  if ( !fp ) {
    return 1;
  }

  if ( !bin ) {
    /* Search last char of file before EOF */
    ier =  fseek(fp, -1, SEEK_END);

    if ( ier ) {
      /* Empty file (fseek has failed) */
      return 1;
    }

    /* File exists and is not empty: search the last char that is not a line
     * separator */
    c = fgetc(fp);
    while ( c == linesep ) {
      ier = fseek(fp, -2, SEEK_CUR);
      if ( ier ) {
        /* The file contains only one line break */
        break;
      }
      c = fgetc(fp);
    }

    /* Search the end of the previous line and compute the length of the last one */
    int len = 1;
    while ( c != linesep ) {
      ier = fseek(fp, -2, SEEK_CUR);
      if ( ier ) {
        /* The file contains only one line break */
        break;
      }
      ++len;
      c = fgetc(fp);
    }

    if ( ier ) {
      /* fseek has failed because we reach the file beginning: go back to the very
       * beginning of file. */
      rewind(fp);
    }
    else {
      char chaine[MMG5_FILESTR_LGTH];
      pos = ftell(fp);

      ier = fscanf(fp,"%127s",&chaine[0]);

      if ( !strncmp(chaine,"End",strlen("End") ) ) {

        rewind(fp);
        fseek(fp,pos,SEEK_SET);
      }
    }
  }
  else {
    /* Append to file (because we are not able to compute the new position of
     * the end key so it brokes medit to remove it). */
    ier = fseek(fp,0, SEEK_END);
    if ( ier ) {
      /* No EOF... why? */
      fprintf(stderr,"  ## Error: %s: Unexpected error\n.",__func__);
      return 0;
    }
  }
  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename file name (if null, print on stdout).
 *
 * \return 0 if fail, 1 otherwise
 *
 * Print parallel communicator in ASCII format: communicators are printed in
 * stdout if no provided filename.
 * Otherwise:
 *   - if the file doesn't exists with Medit extension, it is created;
 *   - if it exists, communicators are appened at the end of the existing Medit
 *      file, removing a possibly existing "End"" keyword
 *
 */
int PMMG_printCommunicator( PMMG_pParMesh parmesh,const char* filename ) {
  PMMG_pExt_comm ext_comm;
  int   **idx_loc,**idx_glob;
  int   bin,ncomm,color,nitem;
  int   icomm,i,ier,ier_glob;
  FILE *fid;

  /** Step 1: find where to write communicators */
  bin = 0;
  ier = 1;
  if ( filename ) {
    ier = MMG3D_openMesh(PMMG_VERB_NO,filename,&fid,&bin,"rb+","rb+");
    if ( ier == -1 ) {
      /* Memory issue: nothing to do */
      fprintf(stderr," ** lack of memory\n");
    }
    else if ( !ier ) {
      /* File creation */
      ier = MMG3D_openMesh(parmesh->info.imprim,filename,&fid,&bin,"w+","wb+");
    }
    else {
      /* File exists: search for the End keyword and position the file pointer
       * before */
      ier = PMMG_search_filePosition(fid,bin,'\n');
      if ( !ier ) {
        fprintf(stderr,"  ## Error: %s: Unable to position file pointer."
                " Exit.\n",__func__);
      }
    }
  }
  else {
    fid = stdout;
  }

  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);
  if ( ier_glob != 1 ) {
    return ier_glob;
  }

  /** Step 2: compute communicators for output */
  if( parmesh->info.API_mode == PMMG_APIDISTRIB_faces ) {
    MMG5_SAFE_MALLOC ( idx_loc,  parmesh->next_face_comm, int *,
                       fprintf(stderr," ** lack of memory\n");ier=0);
    MMG5_SAFE_MALLOC ( idx_glob, parmesh->next_face_comm, int *,
                       fprintf(stderr," ** lack of memory\n");ier=0);

    for( icomm = 0; icomm < parmesh->next_face_comm; icomm++ ) {
      ext_comm = &parmesh->ext_face_comm[icomm];
      nitem    = ext_comm->nitem;
      MMG5_SAFE_MALLOC ( idx_loc[icomm],nitem, int,
                         fprintf(stderr," ** lack of memory\n");ier=0);
      MMG5_SAFE_MALLOC ( idx_glob[icomm],nitem, int,
                         fprintf(stderr," ** lack of memory\n");ier=0);
    }

    ier = PMMG_Get_FaceCommunicator_faces(parmesh, idx_loc);

    MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);
    if ( ier_glob != 1 ) {
      fprintf(stderr,"\n  ## Error: %s: unable to compute face communicators.\n",
              __func__);
      return ier_glob;
    }
    ier = PMMG_Get_FaceCommunicator_owners(parmesh,NULL,idx_glob,NULL,NULL);

  }
  else {
    MMG5_SAFE_MALLOC ( idx_loc,  parmesh->next_node_comm, int *,
                       fprintf(stderr," ** lack of memory\n");ier=0);
    MMG5_SAFE_MALLOC ( idx_glob, parmesh->next_node_comm, int *,
                       fprintf(stderr," ** lack of memory\n");ier=0);

    for( icomm = 0; icomm < parmesh->next_node_comm; icomm++ ) {
      ext_comm = &parmesh->ext_node_comm[icomm];
      nitem    = ext_comm->nitem;
      MMG5_SAFE_MALLOC ( idx_loc[icomm],nitem, int,
                         fprintf(stderr," ** lack of memory\n");ier=0);
      MMG5_SAFE_MALLOC ( idx_glob[icomm],nitem, int,
                         fprintf(stderr," ** lack of memory\n");ier=0);
    }

    ier = PMMG_Get_NodeCommunicator_nodes(parmesh, idx_loc);

    MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);
    if ( ier_glob != 1 ) {
      fprintf(stderr,"\n  ## Error: %s: unable to compute node communicators.\n",
              __func__);

      return ier_glob;
    }

    ier = PMMG_Get_NodeCommunicator_owners(parmesh,NULL,idx_glob,NULL,NULL);
  }

  /** Step 3: file saving */
  if ( !bin ) {
    if( parmesh->info.API_mode == PMMG_APIDISTRIB_faces ) {
      ncomm = parmesh->next_face_comm;
      fprintf(fid,"\nParallelTriangleCommunicators\n%d\n",ncomm);
      for( icomm = 0; icomm < ncomm; icomm++ ) {
        ext_comm = &parmesh->ext_face_comm[icomm];
        color = ext_comm->color_out;
        nitem = ext_comm->nitem;
        fprintf(fid,"%d %d\n",color,nitem);
      }
      fprintf(fid,"\nParallelCommunicatorTriangles\n");
      for( icomm = 0; icomm < ncomm; icomm++ ) {
        ext_comm = &parmesh->ext_face_comm[icomm];
        color = ext_comm->color_out;
        nitem = ext_comm->nitem;
        if( idx_glob ) {
          for( i = 0; i < nitem; i++ ) {
            fprintf(fid,"%d %d %d\n",idx_loc[icomm][i],idx_glob[icomm][i],icomm);
          }
        }
        else {
          for( i = 0; i < nitem; i++ ) {
            fprintf(fid,"%d -1 %d\n",idx_loc[icomm][i],icomm);
          }
        }
      }

      /* Free mem */
      for( icomm = 0; icomm < parmesh->next_face_comm; icomm++ ) {
        MMG5_SAFE_FREE ( idx_loc[icomm] );
        MMG5_SAFE_FREE ( idx_glob[icomm] );
      }
      MMG5_SAFE_FREE ( idx_loc );
      MMG5_SAFE_FREE ( idx_glob );

    } else if( parmesh->info.API_mode == PMMG_APIDISTRIB_nodes ) {
      ncomm = parmesh->next_node_comm;
      fprintf(fid,"ParallelVertexCommunicators\n%d\n",ncomm);
      for( icomm = 0; icomm < ncomm; icomm++ ) {
        ext_comm = &parmesh->ext_node_comm[icomm];
        color = ext_comm->color_out;
        nitem = ext_comm->nitem;
        fprintf(fid,"%d %d\n",color,nitem);
      }
      fprintf(fid,"\nParallelCommunicatorVertices\n");
      for( icomm = 0; icomm < ncomm; icomm++ ) {
        ext_comm = &parmesh->ext_node_comm[icomm];
        color = ext_comm->color_out;
        nitem = ext_comm->nitem;
        if( idx_glob ) {
          for( i = 0; i < nitem; i++ ) {
            fprintf(fid,"%d %d %d\n",idx_loc[icomm][i],idx_glob[icomm][i],icomm);
          }
        }
        else {
          for( i = 0; i < nitem; i++ ) {
            fprintf(fid,"%d -1 %d\n",idx_loc[icomm][i],icomm);
          }
        }
      }

      /* Free mem */
      for( icomm = 0; icomm < parmesh->next_node_comm; icomm++ ) {
        MMG5_SAFE_FREE ( idx_loc[icomm] );
        MMG5_SAFE_FREE ( idx_glob[icomm] );
      }
      MMG5_SAFE_FREE ( idx_loc );
      MMG5_SAFE_FREE ( idx_glob );
    }
    fprintf(fid,"\n\nEnd\n");
  }
  else {
    fprintf(stderr,"  ## Error: %s: Binary file format not yet implemented"
            " for communicators. Exit.\n",__func__);
    return 0;
  }

  if( filename ) {
    fclose(fid);
  }

  return 1;
}

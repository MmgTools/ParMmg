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

#define PMMG_UNRECOGNIZED_ARG                                           \
  do {                                                                  \
    PMMG_ERROR_ARG("\nUnrecognized option %s\n",pmmgArgv,i);            \
  } while(0)

#define PMMG_ERROR_ARG(mess,argv_s,i)                                   \
  do {                                                                  \
    RUN_ON_ROOT_AND_BCAST(                                              \
      fprintf(stderr,mess,argv_s[i]) &&                                 \
      fprintf(stderr,"Please, run %s -h command to get help.\n",argv_s[0]) && \
      0 ,parmesh->info.root,parmesh->myrank,                            \
      ret_val = 0;goto clean );                                         \
  } while(0)

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
    fprintf(stdout,"-met   file  load metric file\n");
    fprintf(stdout,"-field file  load sol field to interpolate from init onto final mesh\n");
    fprintf(stdout,"-noout       do not write output triangulation\n");
    fprintf(stdout,"-centralized-output centralized output (Medit format only)\n");
    fprintf(stdout,"-distributed-output distributed output (Medit format only)\n");

    fprintf(stdout,"\n**  Mode specifications (mesh adaptation by default)\n");
    fprintf(stdout,"-ls     val create mesh of isovalue val (0 if no argument provided)\n");

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
  int        mmgArgc = 0, pmmgArgc = 0;
  char       **mmgArgv = NULL,**pmmgArgv = NULL;

  assert ( parmesh->ngrp == 1 && "Not available for more than 1 group per proc.\n");

  /** First step: search if user want to see the default parameters values or is
   * asking for help */
  for ( i = 1; i < argc; ++i ) {
    if ( !strcmp( argv[ i ],"-val" ) ) {
      RUN_ON_ROOT_AND_BCAST( (PMMG_defaultValues(parmesh) && 0),0,
                             parmesh->myrank,ret_val=0; goto clean);
    }
    else if ( ( !strcmp( argv[ i ],"-?" ) ) || ( !strcmp( argv[ i ],"-h" ) ) ) {
      RUN_ON_ROOT_AND_BCAST( (PMMG_usage(parmesh, argv[0]) && 0),0,
                             parmesh->myrank,ret_val=0; goto clean);
    }
  }


  /** Second step: intercept ParMmg args that exists in Mmg but asks for a
   * specific treatment ( m, v, d) */

  /* Create a new set of argc/argv variables adding only the the cl options that
     mmg has to process
     Overallocating as they are at most argc. Trying to avoid the overallocation
     is not worth any effort, these are ~kb */
  MMG5_SAFE_MALLOC( mmgArgv, argc, char*,ret_val = 0; goto clean);

  /* First argument is always argv[0] ie prog name */
  i = 0;
  MMG_ARGV_APPEND(argv, mmgArgv, i, mmgArgc,ret_val = 0; goto clean);

  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch( argv[i][1] ) {
      case 'm':
        if ( !strcmp(argv[i],"-m") ) {
          /* memory */
          if ( ++i < argc && isdigit( argv[i][0] ) ) {
            if ( ( atoi(argv[ i ]) > MMG5_memSize() ) || ( atoi(argv[ i ]) < 0 ) ) {
              fprintf( stderr,
                       "\nErroneous mem size requested (%s)\n",argv[i] );
              ret_val = 0;
              goto clean;
            }
            else {
              parmesh->info.mem = atoi( argv[i] );
              PMMG_parmesh_SetMemGloMax( parmesh );
            }
            PMMG_parmesh_SetMemMax( parmesh );
          } else {
            PMMG_ERROR_ARG("\nMissing argument option %s\n",argv,i-1);
          }
        }
        else {
          /*  Arg starts by '-m' but doesn't have to be intercepted: Append to
           *  list of args to send to Mmg */
          MMG_ARGV_APPEND(argv, mmgArgv, i, mmgArgc,ret_val = 0; goto clean);
        }
        break;

      case 'd':
        if ( !strcmp(argv[i],"-d") ) {
          /* debug */
          if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_debug,1) )  {
            ret_val = 0;
            goto clean;
          }
        }
        else {
          /*  Arg starts by '-d' but doesn't have to be intercepted: Append to
           *  list of args to send to Mmg */
          MMG_ARGV_APPEND(argv, mmgArgv, i, mmgArgc,ret_val = 0; goto clean);
        }
        break;
      case 'v':  /* verbosity */
        if ( !strcmp(argv[i],"-v") ) {
          if ( ++i < argc && ( isdigit(argv[i][0]) ||
               (argv[i][0]=='-' && isdigit(argv[i][1])) ) ) {
            if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_verbose,atoi(argv[i])) ) {
              ret_val = 0;
              goto clean;
            }
          }
          else {
            i--;
            PMMG_ERROR_ARG("\nMissing argument option for %s\n",mmgArgv,i);
          }
        }
        else {
          /*  Arg starts by '-v' but doesn't have to be intercepted: Append to
           *  list of args to send to Mmg */
          MMG_ARGV_APPEND(argv, mmgArgv, i, mmgArgc,ret_val = 0; goto clean);
        }
        break;
      default:
        /* Arg starts by '-' but doesn't have to be intercepted: Append to list
         * of args to send to Mmg */
        MMG_ARGV_APPEND(argv, mmgArgv, i, mmgArgc,ret_val = 0; goto clean);
        break;
      }
    }
    else {
      /* Arg doesn't start with '-': Append to list of args to send to Mmg */
      MMG_ARGV_APPEND(argv, mmgArgv, i, mmgArgc,ret_val = 0; goto clean);
    }
    ++i;
  }

  /** Third step: Let Mmg parse args it knows among remaining and append in
     pmmgArgv structure unknown ones */
  MMG5_SAFE_MALLOC( pmmgArgv, mmgArgc, char*,ret_val = 0; goto clean);

  i = 0;

  MMG_ARGV_APPEND(argv, pmmgArgv, i, pmmgArgc,ret_val = 0; goto clean);
  MMG5_pMesh mesh = parmesh->listgrp[0].mesh;
  MMG5_pSol  met  = parmesh->listgrp[0].met;
  MMG5_pSol  sol  = parmesh->listgrp[0].ls; // Ok for now as // disp is not planned
  MMG3D_storeknownar(mmgArgc,mmgArgv,mesh,met,sol,&pmmgArgc,pmmgArgv);

  /** Fourth step: parse remaining args with parmmg */
  i = 1;
  while ( i < pmmgArgc ) {
    if ( *pmmgArgv[i] == '-' ) {
      switch( pmmgArgv[i][1] ) {
      case 'c':
        if ( !strcmp(pmmgArgv[i],"-centralized-output") ) {
          /* force centralized output: only relevant using medit distributed
           * input or library call */
          if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_distributedOutput,0) )  {
            ret_val = 0;
            goto clean;
          }
        }
        else {
          PMMG_UNRECOGNIZED_ARG;
        }
        break;
      case 'f':
        if ( !strcmp(pmmgArgv[i],"-field") ) {
          if ( ++i < pmmgArgc && isascii(pmmgArgv[i][0]) && pmmgArgv[i][0]!='-' ) {
            if ( ! PMMG_Set_inputSolsName(parmesh,pmmgArgv[i]) ) {
              fprintf(stderr,"\nUnable to set filname for %s\n",pmmgArgv[i-1]);
              ret_val = 0;
              goto clean;
            }
          }
          else {
            PMMG_ERROR_ARG("\nMissing filname for %s\n",pmmgArgv,i-1);
          }
        }
        else {
          PMMG_UNRECOGNIZED_ARG;
        }
        break;
      case 'g':
        if ( !strcmp(pmmgArgv[i],"-groups-ratio") ) {

          if ( ++i < pmmgArgc ) {
            if ( isdigit(pmmgArgv[i][0]) ) {

              if ( !PMMG_Set_dparameter(parmesh,PMMG_DPARAM_groupsRatio,atof(pmmgArgv[i])) ) {
                ret_val = 0;
                goto clean;
              }
            }
            else {
              i--;
            }
          }
          else {
            PMMG_ERROR_ARG("\nMissing argument option %s\n",pmmgArgv,i-1);
          }
        }
        else {
          PMMG_UNRECOGNIZED_ARG;
        }
        break;

      case 'm':
        if ( !strcmp(pmmgArgv[i],"-mmg-v") ) {
          /* Mmg verbosity */
          if ( ++i < pmmgArgc ) {
            if ( isdigit(pmmgArgv[i][0]) ||
                 (pmmgArgv[i][0]=='-' && isdigit(pmmgArgv[i][1])) ) {
              val = atoi(pmmgArgv[i]);

              if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_mmgVerbose,val) ) {
                ret_val = 0;
                goto clean;
              }
            }
            else {
              i--;
            }
          }
          else {
            PMMG_ERROR_ARG("\nMissing argument option %s\n",pmmgArgv,i-1);
          }
        }
        else if ( !strcmp(pmmgArgv[i],"-mesh-size") ) {

          /* Remesher target mesh size */
          if ( ++i < pmmgArgc && isdigit(pmmgArgv[i][0]) ) {
            val = atoi(pmmgArgv[i]);
            if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_meshSize,val) ) {
              ret_val = 0;
              goto clean;
            }
          }
          else {
            PMMG_ERROR_ARG("\nMissing argument option %s\n",pmmgArgv,i-1);
          }
        }
        else if ( !strcmp(pmmgArgv[i],"-metis-ratio") ) {

          /* Number of metis super nodes per mesh */
          if ( ++i < pmmgArgc && isdigit(pmmgArgv[i][0]) ) {
            val = atoi(pmmgArgv[i]);

            if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_metisRatio,val) ) {
              ret_val = 0;
              goto clean;
            }
          }
          else {
            PMMG_ERROR_ARG("\nMissing argument option %s\n",pmmgArgv,i-1);
          }
        }
        else if ( !strcmp(pmmgArgv[i],"-mmg-d") ) {
          if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_mmgDebug,val) ) {
            ret_val = 0;
            goto clean;
          }
        }
        else {
          PMMG_UNRECOGNIZED_ARG;
        }
        break;

      case 'n':  /* number of adaptation iterations */
        if ( ( 0 == strncmp( pmmgArgv[i], "-niter", 5 ) ) && ( ( i + 1 ) < pmmgArgc ) ) {
          ++i;
          if ( isdigit( pmmgArgv[i][0] ) && ( atoi( pmmgArgv[i] ) >= 0 ) ) {
            parmesh->niter = atoi( pmmgArgv[i] );
          } else {
            parmesh->niter = PMMG_NITER;
            fprintf( stderr,
                     "\nWrong number of adaptation iterations (%s).\n",pmmgArgv[i]);

            ret_val = 0;
            goto clean;
          }
        } else if ( ( 0 == strncmp( pmmgArgv[i], "-nlayers", 5 ) ) && ( ( i + 1 ) < pmmgArgc ) ) {
          ++i;
          if ( isdigit( pmmgArgv[i][0] ) && ( atoi( pmmgArgv[i] ) > 0 ) ) {
            parmesh->info.ifc_layers = atoi( pmmgArgv[i] );
          } else {
            parmesh->info.ifc_layers = PMMG_MVIFCS_NLAYERS;
            fprintf( stderr,
                     "\nWrong number of layers for interface displacement (%s).\n",pmmgArgv[i]);

            ret_val = 0;
            goto clean;
          }
        } else if ( 0 == strncmp( pmmgArgv[i], "-nobalance", 9 ) ) {
          parmesh->info.nobalancing = MMG5_ON;
        } else if ( 0 == strncmp( pmmgArgv[i], "-noout", 5 ) ) {
          parmesh->info.fmtout = PMMG_UNSET;
        }
        else {
          PMMG_UNRECOGNIZED_ARG;
        }
        break;

      case 'd':
        if ( !strcmp(pmmgArgv[i],"-distributed-output") ) {
          /* force distributed output: only relevant using medit centralized
           * input or library call */
          if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_distributedOutput,1) )  {
            ret_val = 0;
            goto clean;
          }
        }
        else {
          PMMG_UNRECOGNIZED_ARG;
        }
        break;
      case 'v':  /* verbosity */
        if ( ++i < pmmgArgc ) {
          if ( isdigit(pmmgArgv[i][0]) ||
               (pmmgArgv[i][0]=='-' && isdigit(pmmgArgv[i][1])) ) {
            if ( !PMMG_Set_iparameter(parmesh,PMMG_IPARAM_verbose,atoi(pmmgArgv[i])) ) {
              ret_val = 0;
              goto clean;
            }
          }
          else
            i--;
        }
        else {
          fprintf(stderr,"\nMissing argument option %s\n",pmmgArgv[i-1]);
          ret_val = 0;
          goto clean;
        }
        break;

      default:
        PMMG_UNRECOGNIZED_ARG;
        break;
      }
    } else {
      if ( parmesh->meshin == NULL && mesh->namein == NULL ) {
        if ( !PMMG_Set_inputMeshName(parmesh,pmmgArgv[i]) ) {
          ret_val = 0;
          goto clean;
        }
      }
      else if ( parmesh->meshout == NULL && mesh->nameout == NULL ) {
        if ( !PMMG_Set_outputMeshName(parmesh,pmmgArgv[i]) ) {
           ret_val = 0;
           goto clean;
        }
      }
      else {
        PMMG_ERROR_ARG("\nArgument %s ignored\n",pmmgArgv,i);
      }
    }
    ++i;
  }

  /** Step 5: Transfer options parsed by Mmg toward ParMmg (if needed) and raise
   * errors for unsupported options */
  parmesh->info.iso = parmesh->listgrp[0].mesh->info.iso;
  parmesh->info.fem = parmesh->listgrp[0].mesh->info.fem;
  parmesh->info.sethmin = parmesh->listgrp[0].mesh->info.sethmin;
  parmesh->info.sethmax = parmesh->listgrp[0].mesh->info.sethmax;

  if ( parmesh->listgrp[0].mesh->info.isosurf ) {

    if ( parmesh->myrank == parmesh->info.root ) {
      fprintf(stderr," ## Error: Splitting boundaries on isovalue not yet"
              " implemented.");
    }
    ret_val = 0;
    goto clean;
  }

  if ( parmesh->listgrp[0].mesh->info.lag >=0 ) {

    if ( parmesh->myrank == parmesh->info.root ) {
      fprintf(stderr," ## Error: Lagrangian motion not yet implemented.");
    }
    ret_val = 0;
    goto clean;
  }

  if( parmesh->listgrp[0].mesh->info.opnbdy ) {
    if ( parmesh->info.root == parmesh->myrank ) {
      fprintf(stderr," ## Warning: Surface adaptation not supported with opnbdy."
              "\nSetting nosurf on.\n");
    }
    if ( !MMG3D_Set_iparameter(parmesh->listgrp[0].mesh,NULL,MMG3D_IPARAM_nosurf,1) ) {
      ret_val = 0;
      goto clean;
    }
  }

  /** Step 6: Sychronize parmesh and mesh file names */
  if ( parmesh->meshin ) {
    /* Input mesh name provided without -in command line arg */
    assert ( mesh->namein );
  }
  else {
    if ( mesh->namein ) {
      /* Input mesh name provided with -in command line arg */
      PMMG_Set_name(parmesh,&parmesh->meshin,mesh->namein,"mesh.mesh");
    }
    else {
      /* Input mesh name not provided */
      if ( parmesh->myrank==parmesh->info.root ) {
        fprintf(stderr,"\nMissing input mesh name.\n");
        fprintf(stderr,"Please, run %s -h command to get help.\n",argv[0]);
      }
      ret_val = 0;
      goto clean;
    }
  }

  if ( parmesh->meshout ) {
    /* Output mesh name provided without -out command line arg */
    assert ( mesh->nameout );
  }
  else {
    if ( mesh->nameout ) {
      /* Output mesh name provided with -out command line arg */
      PMMG_Set_name(parmesh,&parmesh->meshout,mesh->nameout,"mesh.o.mesh");
    }
    else {
      /* Output mesh name not provided */
      char *data;
      MMG5_SAFE_CALLOC(data,strlen(parmesh->meshin)+3,char,return 0);
      strncpy(data,parmesh->meshin,strlen(parmesh->meshin)+3);

      char *ext = MMG5_Get_filenameExt(data);
      if ( ext && !strncmp ( ext,".h5",strlen(".h5") ) ) {
        /* .h5 extension is unknown by Mmg: fix this */
        *ext = '\0';
        strcat(data,".o.h5");
        MMG5_Set_outputMeshName( mesh,data );
      }
      else {
        /* Let Mmg deal automatically with all other file formats */
        MMG5_Set_outputMeshName( mesh,"" );
      }
      MMG5_SAFE_FREE(data);

      assert ( mesh->nameout );
      PMMG_Set_name(parmesh,&parmesh->meshout,
                    parmesh->listgrp[0].mesh->nameout,"mesh.o.mesh");
    }
  }

  /* Metric and solution names are always directly parsed inside met, ls and
   * disp field: in adaptation mode, if the metric name is provided using the
   * -sol arg and if a ls/disp structure is allocated inside the parmesh (for
   * now we deal only with the ls case), the metric name has been stored into
   * the ls/disp name and we have to transfer it in the metric structure. */
  if ( met->namein==NULL &&
       !(mesh->info.iso || mesh->info.isosurf || mesh->info.lag>=0) ) {

    if ( sol->namein ) {
      /* A solution name has been provided using -sol option (facultative) */
      if ( !MMG3D_Set_inputSolName(mesh,met,sol->namein) ) {
        RUN_ON_ROOT_AND_BCAST( (PMMG_usage(parmesh, argv[0]) && 0),0,
                               parmesh->myrank,ret_val=0; goto clean);
      }
      MMG5_DEL_MEM(mesh,sol->namein);
    }
  }

  /* If no input solution name has been parse, assign default name to the
   * suitable data structure (metin in adp mode, lsin in ls mode, dispin in lag
   * mode) */
  MMG5_pSol tmp = NULL;
  if ( mesh->info.iso || mesh->info.isosurf ) {
    tmp = parmesh->listgrp[0].ls;
  }
  else if ( mesh->info.lag >=0 ) {
    tmp = parmesh->listgrp[0].disp;
  }
  else {
    tmp = parmesh->listgrp[0].met;
  }
  assert ( tmp );

  if ( tmp->namein == NULL ) {
    if ( !MMG3D_Set_inputSolName(mesh,tmp,"") ) {
      ret_val = 0;
      goto clean;
    }
  }

  /* Assign default output metric name */
  if ( met->nameout == NULL ) {
    if ( !MMG3D_Set_outputSolName(mesh,met,"") )
      return 0;
  }

  /* Transfer solution names into the parmesh */
  assert ( !parmesh->metin );
  assert ( !parmesh->metout );
  assert ( !parmesh->lsin );
  assert ( !parmesh->dispin );

  if ( met && met->namein ) {
    PMMG_Set_name(parmesh,&parmesh->metin,met->namein,"mesh.sol");
  }
  if ( parmesh->listgrp[0].ls && parmesh->listgrp[0].ls->namein ) {
    PMMG_Set_name(parmesh,&parmesh->lsin,
                  parmesh->listgrp[0].ls->namein,"mesh.sol");
  }
  if ( parmesh->listgrp[0].disp && parmesh->listgrp[0].disp->namein ) {
    PMMG_Set_name(parmesh,&parmesh->dispin,
                  parmesh->listgrp[0].disp->namein,"mesh.sol");
  }

  if ( met && met->nameout ) {
    PMMG_Set_name(parmesh,&parmesh->metout,met->nameout,"mesh.o.sol");
  }

clean:
  MMG5_argv_cleanup( mmgArgv, mmgArgc );
  MMG5_argv_cleanup( pmmgArgv, pmmgArgc );

  return ret_val;
}

#undef ARGV_APPEND


int PMMG_parsop ( PMMG_pParMesh parmesh )
{
  MMG5_pMesh mesh;
  int        ier;

  /* We may have ngrp=0 if distributed inputs have been provided on a different
   * number of processes than the ones used for computation */
  assert ( parmesh->ngrp <= 1 && "more than one group per rank not implemented");
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

  /* Don't print communicators  outside the adaptation loop */
  if( parmesh->iter == PMMG_UNSET ) return 1;

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
    fprintf(fid,"\nNumberOfPartitions\n%d\n",parmesh->nprocs);

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

int PMMG_Get_tetFromTria(PMMG_pParMesh parmesh, int ktri, int* ktet, int* iface ){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_tetFromTria(parmesh->listgrp[0].mesh, ktri, ktet, iface));
}

int PMMG_Get_tetsFromTria(PMMG_pParMesh parmesh, int ktri, int ktet[2], int iface[2] ){
  assert ( parmesh->ngrp == 1 );
  return(MMG3D_Get_tetsFromTria(parmesh->listgrp[0].mesh, ktri, ktet, iface));
}

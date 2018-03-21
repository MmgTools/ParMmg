#include "parmmg.h"

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
    PMMG_DEL_MEM(parmesh, mmgArgv[i], strlen( mmgArgv[i] ), char, "Deallocating mmgargv[i]: " );
  PMMG_DEL_MEM(parmesh, mmgArgv, argc, char*, "Deallocating mmgargv: " );
}

/**
 * \param parmesh pointer to pmmg structure
 * \param rank    process's MPI rank
 *
 * set the mmg default values in the first mesh in the listgrp
 * of the pmmg struct of rank 0
 */
static void
PMMG_defaultValues( PMMG_pParMesh parmesh, const int rank )
{
  if ( rank == 0 ) {
    fprintf( stdout, "\n\n\tParMMG\nDefault parameter values:\n\n");
    fprintf( stdout,"# of remeshing iterations (-niter)  : 1\n");
    fprintf( stdout, "\n\n\tMMG");
    _MMG5_mmgDefaultValues( parmesh->listgrp[0].mesh );
  }
  PMMG_exit_and_free( parmesh, PMMG_SUCCESS );
}

/**
 * \param parmesh  pointer to pmmg structure
 * \param progname program name string
 *
 * print the command line usage of the parmmg tool
 */
static void
PMMG_usage( PMMG_pParMesh parmesh, char * const progname )
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
 * \param parmesh pointer to pmmg structure
 * \param memReq  size of memory in Mb. If memReq is zero then it is
 *                automatically set to half of the machine's available memory.
 *                On machines with multicore processors (ie most of today's cpus)
 *                the total available memory is shared equally to pmmg processes
 *                running on the same machine.
 *                If memReq is negative or more than the detected available
 *                memory, then the requested value is discarded and the maximum
 *                allowed memory is set to half the detected available memory.
 *
 *  Sets the maximum amount of memory that a parmmg process is allowed to use.
 *  This includes both the memory used for the parmmg struct and the mmg structs
 *  in listgrp
 */
void PMMG_parmesh_SetMemGloMax( PMMG_pParMesh parmesh, long long int memReq )
{
  long long int maxAvail = 0;
  MPI_Comm comm_shm = 0;
  int size_shm = 1;
  int flag;
  const int million = 1024 * 1024;

  assert ( (parmesh != NULL) && "trying to set glo max mem in empty parmesh" );

  MPI_Initialized( &flag );

  if ( flag ) {
#warning here we must use parmesh->comm and not MPI_COMM_WORLD
    MPI_Comm_split_type( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                         &comm_shm );
    MPI_Comm_size( comm_shm, &size_shm );
  }
  else {
    size_shm = 1;
  }

  maxAvail = _MMG5_memSize();
  // if detection failed => default value of _MMG5_MEMMAX Mo
  if ( maxAvail == 0 )
    maxAvail = _MMG5_MEMMAX << 20;

  // Multiple MPI processes may be running on the same node => distribute equally
  if ( (memReq > 0) && ((memReq * million) < maxAvail) )
    parmesh->memGloMax = (memReq * million) / size_shm;
  else
    parmesh->memGloMax = (maxAvail * 50) / (size_shm * 100);

  fprintf ( stdout,
            "Requested %lld Mb max memory usage. Max memory limit set to %lld Mb\n",
            memReq, parmesh->memGloMax/million );
}


/**
 * \param mesh pointer toward the mesh structure
 * \param met pointer toward the metric structure
 *
 * \return 0 if fail, 1 otherwise
 *
 * memory repartition between the point, xpoint, tetra, xtetra arrays for the
 * memMax amout of memory available.
 *
 */
static inline
int PMMG_memOption_memRepartition(MMG5_pMesh mesh,MMG5_pSol met) {
  long long  million = 1048576L;
  long long  usedMem,avMem,reservedMem;
  long       castedVal;
  int        ctri,npadd,bytes;

  /* init allocation need 38 octets */
  reservedMem = 38 +  (long long)
    (mesh->nprism*sizeof(MMG5_Prism) + mesh->xpr*sizeof(MMG5_Prism) );

  /* Compute the needed initial memory */
  usedMem = reservedMem + (mesh->np+1)*sizeof(MMG5_Point)
    + (mesh->xp+1)*sizeof(MMG5_xPoint) + (mesh->ne+1)*sizeof(MMG5_Tetra)
    + (mesh->xt+1)*sizeof(MMG5_xTetra);

  if ( mesh->adja )
    usedMem += (4*mesh->ne+1)*sizeof(int);

  if ( met->m )
    usedMem += met->size*(mesh->np+1)*sizeof(double);

  if ( usedMem > mesh->memMax  ) {
    fprintf(stderr,"\n  ## Error: %s: %lld Mo of memory ",__func__,mesh->memMax/million);
    castedVal =  _MMG5_SAFELL2LCAST(usedMem/million+1);
    fprintf(stderr,"is not enough to load mesh. You need to ask %ld Mo minimum\n",
            castedVal);
    return 0;
  }

  ctri = 2;
  /* Euler-poincare: ne = 6*np; nt = 2*np; na = np/5 *
   * point+tria+tets+adja+adjt+sol+item */
  bytes = sizeof(MMG5_Point) + sizeof(MMG5_xPoint) +
    6*sizeof(MMG5_Tetra) + ctri*sizeof(MMG5_xTetra);

  if ( mesh->adja )
    bytes += 4*6*sizeof(int);

  if ( met->m )
    bytes += met->size*sizeof(double);

  avMem = mesh->memMax-usedMem;

  npadd = (int) ( (double)avMem/bytes );
  mesh->npmax = MG_MIN(mesh->npmax,mesh->np+npadd);
  mesh->xpmax = MG_MIN(mesh->xpmax,mesh->xp+npadd);
  mesh->nemax = MG_MIN(mesh->nemax,6*npadd+mesh->ne);
  mesh->xtmax = MG_MIN(mesh->xtmax,ctri*npadd+mesh->xt);

  met->npmax  = mesh->npmax;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    castedVal = _MMG5_SAFELL2LCAST(mesh->memMax/million);
    fprintf(stdout,"  MAXIMUM MEMORY AUTHORIZED (Mo)    %ld\n",
            castedVal);

    fprintf(stdout,"  _MMG3D_NPMAX    %d\n",mesh->npmax);
    fprintf(stdout,"  _MMG3D_XPMAX    %d\n",mesh->xpmax);
    fprintf(stdout,"  _MMG3D_NEMAX    %d\n",mesh->nemax);
    fprintf(stdout,"  _MMG3D_XTMAX    %d\n",mesh->xtmax);
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Keep track of empty links for the tetra and points array.
 *
 */
int PMMG_link_mesh( MMG5_pMesh mesh ) {
  MMG5_pTetra pt;
  int k,iadr;

  /* keep track of empty links */
  if ( mesh->npmax > mesh->np ) {
    mesh->npnil = mesh->np + 1;
    for (k=mesh->npnil; k<=mesh->npmax; k++) {
      /* Set tangent field of point to 0 */
      mesh->point[k].n[0] = 0;
      mesh->point[k].n[1] = 0;
      mesh->point[k].n[2] = 0;
      /* link */
      if(k<mesh->npmax-1) mesh->point[k].tmp  = k+1;
    }
  }
  else {
    assert ( mesh->np == mesh->npmax );
    mesh->npnil = 0;
  }

  if ( mesh->nemax > mesh->ne ) {
    mesh->nenil = mesh->ne + 1;
    for (k=mesh->nenil; k<=mesh->nemax; k++) {
      pt = &mesh->tetra[k];
      memset(pt,0,sizeof(MMG5_Tetra));
      iadr = 4*(k-1) + 1;
      if ( mesh->adja )
        memset(&mesh->adja[iadr],0,4*sizeof(int));
      
      if(k<mesh->nemax-1) pt->v[3] = k+1;
    }
  }
  else {
    assert ( mesh->ne == mesh->nemax );
    mesh->nenil = 0;
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param npmax_old old maximum number of points.
 * \param xpmax_old old maximum number of boundary points.
 * \param nemax_old old maximum number of tetra.
 * \param xtmax_old old maximum number of boundary tetra.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Reallocation of the array fields of the mesh for the given xpmax,
 * npmax,xtmax, nemax.
 *
 */
int PMMG_setMemMax_realloc( MMG5_pMesh mesh,int npmax_old,int xpmax_old,
                             int nemax_old,int xtmax_old ) {

  PMMG_RECALLOC(mesh,mesh->point,mesh->npmax+1,npmax_old+1,MMG5_Point,
                "vertices array", return 0);

  PMMG_RECALLOC(mesh,mesh->xpoint,mesh->xpmax+1,xpmax_old+1,MMG5_xPoint,
              "boundary vertices array", return 0);

  PMMG_RECALLOC(mesh,mesh->tetra,mesh->nemax+1,nemax_old+1,MMG5_Tetra,
              "tetra array", return 0);

  if ( mesh->adja ) {
    PMMG_RECALLOC(mesh,mesh->adja,4*mesh->nemax+5,4*nemax_old+5,int,
                  "adja array", return 0);
  }

  PMMG_RECALLOC(mesh,mesh->xtetra,mesh->xtmax+1,xtmax_old+1,MMG5_xTetra,
              "boundary tetra array", return 0);

  return ( PMMG_link_mesh( mesh ) );
}

/**
 * \param parmesh parmesh structure to adjust
 * \param percent integer value bewtween 0 and 100
 * \param fitMesh if 1, set maximum mesh size at its exact size.
 *
 * \return 1 if success, 0 if fail
 *
 * Set the maximum memory parmesh and the meshes in listgrp can use.
 * The total memory available is split between the parmesh structure and the
 * listgrp structures according to the percentage specified by the percent
 * input variable:
 *   percent % of the available mem is assigned to pmesh.memMax
 *   (100-percent)/100 are assigned to the mesh[i].memMax
 */
int PMMG_parmesh_SetMemMax( PMMG_pParMesh parmesh, int percent )
{
  MMG5_pMesh mesh;
  long long  available;
  long long  million = 1048576L;
  int        remaining_ngrps;
  int        i = 0;

  assert ( (0 < percent) && (100 > percent) && "percent has to be >0 and <100" );

  parmesh->memMax = parmesh->memGloMax * percent / 100;
  available       = parmesh->memGloMax - parmesh->memMax;

  if ( available < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: all the memory is used for the communicators\n",
            __func__);
    return 0;
  }

  remaining_ngrps = parmesh->ngrp;
  for ( i = 0; i < parmesh->ngrp; ++i ) {
    mesh = parmesh->listgrp[i].mesh;
    mesh->memMax = available/remaining_ngrps;

    /* Not enough memory: set the minimal memory to be able to continue */
    if ( mesh->memMax < mesh->memCur ) {
      mesh->memMax = mesh->memCur;
    }
    /* Force the mmg3d zaldy function to find the wanted memMax value (in MMG3D_loadMesh) */
    mesh->info.mem = mesh->memMax/million;

    /* Count the remaining available memory */
    available -= mesh->memMax;
    --remaining_ngrps;
    if ( available < 0 ) {
      fprintf(stderr,"\n  ## Error: %s: not enough memory\n",__func__);
      return 0;
    }
  }
  return 1;
}

/**
 * \param parmesh parmesh structure to adjust
 * \param percent increasing ratio of the parmesh memory (>100).
 * \param fitMesh if 1, set maximum mesh size at its exact size.
 *
 * \return 1 if success, 0 if fail
 *
 * Update the memory repartition between the communicators and the groups.
 *
 */
int PMMG_parmesh_updateMemMax( PMMG_pParMesh parmesh, int percent, int fitMesh )
{
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  long long  available;
  long long  million = 1048576L;
  int        remaining_ngrps,npmax_old,xpmax_old,nemax_old,xtmax_old;
  int        i = 0;

  if ( parmesh->memCur >= parmesh->memMax )
    parmesh->memMax = (parmesh->memCur * percent)/100;

  available       = parmesh->memGloMax - parmesh->memMax;

  if ( available < 0 ) {
    fprintf(stderr,"\n  ## Error: %s: all the memory is used for the communicators\n",
            __func__);
    return 0;
  }

  remaining_ngrps = parmesh->ngrp;
  for ( i = 0; i < parmesh->ngrp; ++i ) {
    mesh = parmesh->listgrp[i].mesh;
    mesh->memMax = available/remaining_ngrps;

    /* Not enough memory: set the minimal memory to be able to continue */
    if ( mesh->memMax < mesh->memCur ) {
      mesh->memMax = mesh->memCur;
    }
    /* Force the mmg3d zaldy function to find the wanted memMax value (in MMG3D_loadMesh) */
    mesh->info.mem = mesh->memMax/million;

    /* Memory repartition for the MMG meshes arrays */
    npmax_old = mesh->npmax;
    xpmax_old = mesh->xpmax;
    nemax_old = mesh->nemax;
    xtmax_old = mesh->xtmax;
    if ( fitMesh ) {
      mesh->npmax = mesh->np;
      mesh->xpmax = mesh->xp;
      mesh->nemax = mesh->ne;
      mesh->xtmax = mesh->xt;
    }
    else {
      mesh->npmax = MG_MAX(1.5*mesh->np,_MMG3D_NPMAX);
      mesh->xpmax = 1.5*mesh->xp;
      mesh->nemax = MG_MAX(1.5*mesh->ne,_MMG3D_NEMAX);
      mesh->xtmax = MG_MAX(1.5*mesh->xt,_MMG3D_NTMAX);
    }

    met = parmesh->listgrp[i].met;
    if ( !PMMG_memOption_memRepartition(mesh,met) ) return 0;

    if ( !PMMG_setMemMax_realloc(mesh,npmax_old,xpmax_old,nemax_old,xtmax_old) )
      return 0;

    met->np     = mesh->np;
    met->npmax  = mesh->npmax;
    if ( met->m )
      PMMG_REALLOC(mesh,met->m,met->size*(met->npmax+1),met->size*(npmax_old+1),
                   double,"metric array",return 0);

    /* Count the remaining available memory */
    available -= mesh->memMax;
    --remaining_ngrps;
    if ( available < 0 ) {
      fprintf(stderr,"\n  ## Error: %s: not enough memory\n",__func__);
      return 0;
    }
  }
  return 1;
}

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

/**
* \param argc the argument count parameter from main
* \param argv the argument values parameter from main
* \param parmesh pointer to active parmesh
*
* \return 1 on success
*         0 on failure
*
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
  int ret_val = 1;
  int mmgArgc = 0;
  char** mmgArgv = NULL;

  for ( i = 1; i < argc; ++i )
    if ( !strcmp( argv[ i ],"-val" ) )
      PMMG_defaultValues( parmesh, parmesh->myrank );
    else if ( ( !strcmp( argv[ i ],"-?" ) ) || ( !strcmp( argv[ i ],"-h" ) ) )
      PMMG_usage( parmesh, argv[0] );

  // Create a new set of argc/argv variables adding only the the cl options that
  // mmg has to process
  // Overallocating as they are at most argc. Trying to avoid the overallocation
  // is not worth any effort, these are ~kb
  PMMG_MALLOC(parmesh, mmgArgv, argc, char*, " copy of argv for mmg: ",
              ret_val = 0; goto fail_mmgargv);

  // First argument is always argv[0] ie prog name
  i = 0;
  ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc, " mmgArgv[0] for mmg: ",
              ret_val = 0; goto fail_proc);

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
            PMMG_parmesh_SetMemGloMax( parmesh, atoi( argv[i] ) );
          PMMG_parmesh_SetMemMax( parmesh, 20 );
        } else {
          fprintf( stderr, "Missing argument option %c\n", argv[i-1][1] );
          PMMG_usage( parmesh, argv[0] );
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
                      ret_val = 0; goto fail_proc );
        }
      break;

      case 'd':  // number of adaptation iterations
        parmesh->ddebug = 1;
        ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                    " adding to mmgArgv for mmg: ",
                    ret_val = 0; goto fail_proc );
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

  // parmmg finished parsing arguments, the rest will e handled by mmg3d
  if ( 1 != MMG3D_parsar( mmgArgc, mmgArgv,
                          parmesh->listgrp[0].mesh,
                          parmesh->listgrp[0].met ) ) {
    ret_val = 0;
    goto fail_proc;
  }

fail_proc:
  PMMG_argv_cleanup( parmesh, mmgArgv, mmgArgc, argc );
fail_mmgargv:
  return ret_val;
}

#undef ARGV_APPEND

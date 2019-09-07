/**
 * \file zaldy_pmmg.c
 * \brief Memory management
 * \copyright GNU Lesser General Public License.
 *
 * Memory management for PARMMG library.
 *
 */

#include "parmmg.h"

/**
 * \param parmesh pointer to pmmg structure
 * \param memReq  size of memory in Mb. If memReq is zero then it is
 *                automatically set to half of the machine's available memory.
 *                On machines with multicore processors
 *                the total available memory is shared equally to pmmg processes
 *                running on the same machine.
 *                If memReq is negative or more than the detected available
 *                memory, then the requested value is discarded and the maximum
 *                allowed memory is set to half the detected available memory.
 *
 *  Sets the maximum amount of memory that a parmmg process is allowed
 *  to use depending on the user specifications and the available
 *  memory. On machines with multicore processors the memory is shared
 *  equally to pmmg processes running on the same machine.  This
 *  includes both the memory used for the parmmg struct and the mmg
 *  structs in listgrp.
 */
void PMMG_parmesh_SetMemGloMax( PMMG_pParMesh parmesh )
{
  size_t   maxAvail = 0;
  MPI_Comm comm_shm = 0;
  int      size_shm = 1;
  int      flag;

  assert ( (parmesh != NULL) && "trying to set glo max mem in empty parmesh" );

  /** Step 1: Get the numper of processes per node */
  MPI_Initialized( &flag );

  if ( flag ) {
    MPI_Comm_split_type( parmesh->comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                         &comm_shm );
    MPI_Comm_size( comm_shm, &size_shm );
  }
  else {
    size_shm = 1;
  }

  /** Step 2: Set maximal memory per process depending on the -m option setting */
  maxAvail = MMG5_memSize()/size_shm;

  if ( parmesh->info.mem <= 0 ) {
    /* Nos users specifications */
    if ( !maxAvail ) {
      /* default value when not able to compute the available memory = 800 MB */
      printf("  Maximum memory per process set to default value: %d MB.\n",MMG5_MEMMAX);
      parmesh->memGloMax = MMG5_MEMMAX << 20;
    }
    else {
      /* maximal memory = 50% of total physical memory */
      parmesh->memGloMax = maxAvail * MMG5_MEMPERCENT;
    }
  }
  else {
    /* memory asked by user if possible, otherwise total physical memory */
    if ( maxAvail && (size_t)parmesh->info.mem*MMG5_MILLION > maxAvail ) {
      fprintf(stderr,"\n  ## Warning: %s: asking for %d MB of memory per process ",
              __func__,parmesh->info.mem);
      fprintf(stderr,"when only %zu available.\n",maxAvail/MMG5_MILLION);
    }
    else {
      parmesh->memGloMax= (size_t)parmesh->info.mem*MMG5_MILLION;
    }
  }

  if ( abs(parmesh->info.imprim) > 4 || parmesh->ddebug ) {
    fprintf(stdout,"  MAXIMUM MEMORY AUTHORIZED PER PROCESS (MB)    %zu\n",
            parmesh->memGloMax/MMG5_MILLION);
  }
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
  size_t     usedMem,avMem,reservedMem;
  int        ctri,npadd,bytes;

  /* init allocation need 38 octets */
  reservedMem = 38 +  (size_t)
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
    fprintf(stderr,"\n  ## Error: %s: %zu Mo of memory ",__func__,mesh->memMax/MMG5_MILLION);
    fprintf(stderr,"is not enough to load mesh. You need to ask %zu Mo per process minimum.\n",
            usedMem/MMG5_MILLION+1);
    fprintf(stderr,"\nTry to use the -m option to impose the maximal memory per process.\n");
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
    fprintf(stdout,"  MAXIMUM MEMORY PER PROCESS AUTHORIZED (Mo)    %zu\n",
            mesh->memMax/MMG5_MILLION);

    fprintf(stdout,"  MMG3D_NPMAX    %d\n",mesh->npmax);
    fprintf(stdout,"  MMG3D_XPMAX    %d\n",mesh->xpmax);
    fprintf(stdout,"  MMG3D_NEMAX    %d\n",mesh->nemax);
    fprintf(stdout,"  MMG3D_XTMAX    %d\n",mesh->xtmax);
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
    for (k=mesh->npnil; k<mesh->npmax-1; k++)
      mesh->point[k].tmp  = k+1;
  }
  else {
    assert ( mesh->np == mesh->npmax );
    mesh->npnil = 0;
  }

  if ( mesh->nemax > mesh->ne ) {
    mesh->nenil = mesh->ne + 1;
    for (k=mesh->nenil; k<mesh->nemax-1; k++)
      mesh->tetra[k].v[3] = k+1;
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
 *
 * \return 1 if success, 0 if fail
 *
 * Set the maximum memory that parmesh and the meshes in listgrp can use.
 * The total memory available is split between the parmesh structure and the
 * listgrp structures according to the percentage specified by the percent
 * input variable:
 *   percent % of the available mem is assigned to pmesh.memMax
 *   (100-percent)/100 are assigned to the mesh[i].memMax
 */
int PMMG_parmesh_SetMemMax( PMMG_pParMesh parmesh, int percent )
{
  MMG5_pMesh mesh;
  size_t     available;
  int        remaining_ngrps;
  int        i = 0;

  assert ( (0 < percent) && (100 > percent) && "percent has to be >0 and <100" );

  parmesh->memMax = parmesh->memGloMax * percent / 100;

  if ( parmesh->memGloMax <= parmesh->memMax ) {
    fprintf(stderr,"\n  ## Error: %s: all the memory is used for the communicators\n",
            __func__);
    return 0;
  }
  available       = parmesh->memGloMax - parmesh->memMax;

  remaining_ngrps = parmesh->ngrp;
  for ( i = 0; i < parmesh->ngrp; ++i ) {
    mesh = parmesh->listgrp[i].mesh;
    mesh->memMax = available/remaining_ngrps;

    /* Not enough memory: set the minimal memory to be able to continue */
    if ( mesh->memMax < mesh->memCur ) {
      mesh->memMax = mesh->memCur;
    }
    /* Force the mmg3d zaldy function to find the wanted memMax value (in MMG3D_loadMesh) */
    mesh->info.mem = mesh->memMax/MMG5_MILLION;

    /* Count the remaining available memory */
    if ( available < mesh->memMax ) {
      fprintf(stderr,"\n  ## Error: %s: not enough memory\n",__func__);
      return 0;
    }

    available -= mesh->memMax;
    --remaining_ngrps;
  }
  return 1;
}

/**
 * \param parmesh pointer toward a parmesh.
 * \param mesh pointer toward the mesh that we want to fit.
 * \param met pointer toward the metric that we want to fit.
 *
 * \return 1 if success, 0 if fail.
 *
 * Set the memMax field of a mesh to memCur and realloc the mesh in order to use
 * the less possible memory.
 *
 */
int PMMG_parmesh_fitMesh( PMMG_pParMesh parmesh, MMG5_pMesh mesh, MMG5_pSol met ) {
  int npmax_old,xpmax_old,nemax_old,xtmax_old;
  int ier = 1;

  npmax_old = mesh->npmax;
  xpmax_old = mesh->xpmax;
  nemax_old = mesh->nemax;
  xtmax_old = mesh->xtmax;

  mesh->npmax = mesh->np;
  mesh->xpmax = mesh->xp;
  mesh->nemax = mesh->ne;
  mesh->xtmax = mesh->xt;

  met->npmax = mesh->npmax;

  if ( !PMMG_setMemMax_realloc(mesh,npmax_old,xpmax_old,
                               nemax_old,xtmax_old) ) ier = 0;

  if ( met->m ) {
    PMMG_REALLOC(parmesh,met->m,met->size*(met->npmax+1),
                 met->size*(npmax_old+1),double,"metric_array",
                 assert(0);ier = 0;);
  }
  mesh->memMax = mesh->memCur;

  return ier;
}

/**
 * \param parmesh parmesh structure to adjust
 * \param percent ratio of the available memory to give to the parmesh
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
  size_t     available;
  int        remaining_ngrps,npmax_old,xpmax_old,nemax_old,xtmax_old;
  int        k,i = 0;

  parmesh->memMax = parmesh->memCur;

  if (  parmesh->memGloMax <=  parmesh->memMax ) {
    fprintf(stderr,"\n  ## Error: %s: all the memory is used for the communicators\n",
            __func__);
    return 0;
  }

  available       = parmesh->memGloMax - parmesh->memMax;

  for ( k=0; k<parmesh->ngrp; ++k ) {
    parmesh->listgrp[k].mesh->memMax = parmesh->listgrp[k].mesh->memCur;

    if ( available < parmesh->listgrp[k].mesh->memMax ) {
      fprintf(stderr,"\n  ## Error: %s: all the memory is used for the communicators\n",
              __func__);
      return 0;
    }
    available -= parmesh->listgrp[k].mesh->memMax;
 }

  parmesh->memMax += percent * available/100;

  remaining_ngrps = parmesh->ngrp;
  for ( i = 0; i < parmesh->ngrp; ++i ) {
    mesh = parmesh->listgrp[i].mesh;
    mesh->memMax = available/remaining_ngrps;

    /* Not enough memory: set the minimal memory to be able to continue */
    if ( mesh->memMax < mesh->memCur ) {
      mesh->memMax = mesh->memCur;
    }
    /* Force the mmg3d zaldy function to find the wanted memMax value (in MMG3D_loadMesh) */
    mesh->info.mem = mesh->memMax/MMG5_MILLION;

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
      mesh->npmax = MG_MAX(1.5*mesh->np,MMG3D_NPMAX);
      mesh->xpmax = 1.5*mesh->xp;
      mesh->nemax = MG_MAX(1.5*mesh->ne,MMG3D_NEMAX);
      mesh->xtmax = MG_MAX(1.5*mesh->xt,MMG3D_NTMAX);
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
    if ( available < mesh->memMax ) {
      fprintf(stderr,"\n  ## Error: %s: not enough memory\n",__func__);
      return 0;
    }

    available -= mesh->memMax;
    --remaining_ngrps;
  }
  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param ext_comm pointer toward the external communicator to resize
 * \param newSize new size of the external comm
 * \param oldSize old size of the external comm
 *
 * \return 0 if fail, 1 otherwise
 *
 * Resize the int_comm_index array of an external communicator.
 *
 */
int PMMG_resize_extComm ( PMMG_pParMesh parmesh,PMMG_pExt_comm ext_comm,
                          int newSize,int *oldSize ) {

  if ( newSize == *oldSize ) return 1;

  PMMG_REALLOC(parmesh,ext_comm->int_comm_index,newSize,*oldSize,int,
               "int_comm_index",return 0);

  *oldSize = newSize;

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param ext_comm pointer toward the array of external communicators to resize
 * \param newSize new size of the array of external comm
 * \param oldSize old size of the array of external comm
 *
 * \return 0 if fail, 1 otherwise
 *
 * Resize an array of external communicators.
 *
 */
int PMMG_resize_extCommArray ( PMMG_pParMesh parmesh,PMMG_pExt_comm *ext_comm,
                               int newSize,int *oldSize ) {
  int            k;

  if ( newSize == *oldSize ) return 1;

  /* Free the external communicators that will be deleted */
  for ( k=newSize; k<*oldSize; ++k ) {

    if ( !PMMG_resize_extComm ( parmesh,(*ext_comm)+k,0,&(*ext_comm+k)->nitem ) )
      return 0;

    if ( (*ext_comm+k)->itosend )
      PMMG_DEL_MEM ( parmesh,(*ext_comm+k)->itosend,int,"itosend" );

    if ( (*ext_comm+k)->itorecv )
      PMMG_DEL_MEM ( parmesh,(*ext_comm+k)->itorecv,int,"itorecv" );

    if ( (*ext_comm+k)->rtosend )
      PMMG_DEL_MEM ( parmesh,(*ext_comm+k)->rtosend, double,"rtosend" );

    if ( (*ext_comm+k)->rtorecv )
      PMMG_DEL_MEM ( parmesh,(*ext_comm+k)->rtorecv,int,"rtorecv" );

  }

  PMMG_REALLOC(parmesh,*ext_comm,newSize,*oldSize,PMMG_Ext_comm,"ext_comm",return 0);

  if ( newSize > *oldSize )
    memset( *ext_comm + *oldSize, 0x0,
            ((size_t)((newSize)-(*oldSize)))*sizeof(PMMG_Ext_comm));

  *oldSize = newSize;

  return 1;
}

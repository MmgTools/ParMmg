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
  int      flag;

  assert ( (parmesh != NULL) && "trying to set glo max mem in empty parmesh" );

  /** Step 1: Get the number of processes per node */
  MPI_Initialized( &flag );

  if ( flag ) {
    MPI_Comm_split_type( parmesh->comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                         &comm_shm );
    MPI_Comm_size( comm_shm, &parmesh->size_shm );
    MPI_Comm_free( &comm_shm );
  }
  else {
    parmesh->size_shm = 1;
  }

  /** Step 2: Set maximal memory per process depending on the -m option setting:
  - if the user doesn't provides a memory value or provides an invalid value: we equirepartite the memory over the MPI processes of the node. Functions that consumes different amounts of memory depending on the process have to manage internally the memory repartition (for example the \a PMMG_loadMesh_centralized function).
  - if the user provides a valid memory value (under or equal to the physical memory), it is used as is guessing that the user know what he is asking (it may be useful during the parallel calls of Mmg to not have a memory equirepartition as some process may use a smaller amount of memory than others but we are not able to predict it.
 */
  maxAvail = MMG5_memSize();

  if ( parmesh->info.mem <= 0 ) {
    /* Nos users specifications: equirepartition of */
    if ( !maxAvail ) {
      /* default value when not able to compute the available memory = 800 MB */
      printf("  Maximum memory per process set to default value: %d MB.\n",MMG5_MEMMAX);
      parmesh->memGloMax = (MMG5_MEMMAX/parmesh->size_shm) << 20;
    }
    else {
      /* maximal memory = equirepartition of total physical memory over the MPI processes of the node. */
      parmesh->memGloMax = maxAvail/parmesh->size_shm;
    }
  }
  else {
    int memOverflow = 0;
    /* Memory asked by user if possible (authorized to ask the entire memory nod per process, independently of the number of process per node). */
    if ( maxAvail && (size_t)parmesh->info.mem*MMG5_MILLION > maxAvail ) {
      /* User asks for more than the memory of the node */
      fprintf(stdout,"\n  ## Warning: %s: asking for %d MB of memory per process ",
              __func__,parmesh->info.mem);
      fprintf(stdout,"when only %zu available on the node.\n",maxAvail/MMG5_MILLION);
      memOverflow = 1;
    }
    else {
      if ( (size_t)parmesh->info.mem*MMG5_MILLION > maxAvail/parmesh->size_shm ) {
      /* User asks for more than the equirepartition of the node memory across the MPI processes */
        fprintf(stdout,"\n  ## Warning: %s: asking for %d MB per MPI process with %d process per node and %zu MB available on the node.\n",
                __func__,parmesh->info.mem,parmesh->size_shm,maxAvail/MMG5_MILLION);
        memOverflow = 1;
      }
    }

    /* In all cases, impose what the user ask */
    parmesh->memGloMax= (size_t)parmesh->info.mem*MMG5_MILLION;

    if ( memOverflow ) {
      fprintf(stdout,"              The program may run out of memory and be killed (Signal 9 or SIGKILL error).\n\n");
    }
  }

  if ( abs(parmesh->info.imprim) > 4 || parmesh->ddebug ) {
    fprintf(stdout,"  MAXIMUM MEMORY AUTHORIZED PER PROCESS (MB)    %zu\n",
            parmesh->memGloMax/MMG5_MILLION);
  }
}

/**
 * \param parmesh parmesh structure
 *
 * \return 1 if success, 0 if fail
 *
 * Set the maximum memory that parmesh and the meshes in listgrp can use.
 */
int PMMG_parmesh_SetMemMax( PMMG_pParMesh parmesh ) {
  MMG5_pMesh mesh;
  size_t     memMax;
  int        i;

  memMax = parmesh->memGloMax;

  parmesh->memMax = memMax;
  for( i = 0; i < parmesh->ngrp; ++i ) {
    mesh = parmesh->listgrp[i].mesh;
    mesh->memMax = parmesh->memGloMax;

    /* Hack to not let Mmg recomputes the available memory by itself (it has no
     * knowledge that it is called in parallel) */
    mesh->info.mem = mesh->memMax/MMG5_MILLION;
  }

  return 1;
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
  size_t     npadd;
  int        ctri,bytes;

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


  /** Try to estimate the memody usage of adding a point (with the related
   *  tetra, tria, so:ution...) in order to reduce the maximum size when
   *  possible. */

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

  /* The number of points that can be added is approximately given by the
   * ratio between the available memory and the memory usage of a
   * point+related structures */
  npadd = (size_t) ( (double)avMem/bytes );

  /* Shrink size if too big */
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
  int k;

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
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Allocation of the array fields of the mesh for the given npmax, xpmax, nemax,
 * xtmax.
 *
 */
int PMMG_setMeshSize_alloc( MMG5_pMesh mesh ) {

  PMMG_CALLOC(mesh,mesh->point,mesh->npmax+1,MMG5_Point,
              "vertices array", return 0);

  PMMG_CALLOC(mesh,mesh->xpoint,mesh->xpmax+1,MMG5_xPoint,
              "boundary vertices array", return 0);

  PMMG_CALLOC(mesh,mesh->tetra,mesh->nemax+1,MMG5_Tetra,
              "tetra array", return 0);

  PMMG_CALLOC(mesh,mesh->xtetra,mesh->xtmax+1,MMG5_xTetra,
              "boundary tetra array", return 0);

  if ( mesh->nt ) {
    PMMG_CALLOC(mesh,mesh->tria,mesh->nt+1,MMG5_Tria,
                "triangles array", return 0);
  }

  return ( PMMG_link_mesh( mesh ) );
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
 * Reallocation of the array fields of the mesh for the given npmax,
 * xpmax, nemax, xtmax.
 *
 */
int PMMG_setMeshSize_realloc( MMG5_pMesh mesh,int npmax_old,int xpmax_old,
                              int nemax_old,int xtmax_old ) {

  if ( !npmax_old )
    PMMG_CALLOC(mesh, mesh->point, mesh->npmax+1, MMG5_Point,
                "vertices array", return 0);
  else
    PMMG_RECALLOC(mesh,mesh->point,mesh->npmax+1,npmax_old+1,MMG5_Point,
                  "vertices array", return 0);
  if ( !xpmax_old )
    PMMG_CALLOC(mesh, mesh->xpoint, mesh->xpmax+1, MMG5_xPoint,
                "boundary vertices array", return 0);
  else
    PMMG_RECALLOC(mesh,mesh->xpoint,mesh->xpmax+1,xpmax_old+1,MMG5_xPoint,
                  "boundary vertices array", return 0);

  if ( !nemax_old )
    PMMG_CALLOC(mesh, mesh->tetra, mesh->nemax+1, MMG5_Tetra,
                "tetra array", return 0);
  else
    PMMG_RECALLOC(mesh,mesh->tetra,mesh->nemax+1,nemax_old+1,MMG5_Tetra,
                  "tetra array", return 0);

  if ( mesh->adja ) {
    PMMG_RECALLOC(mesh,mesh->adja,4*mesh->nemax+5,4*nemax_old+5,int,
                  "adja array", return 0);
  }

  if ( !xtmax_old )
    PMMG_CALLOC(mesh, mesh->xtetra, mesh->xtmax+1, MMG5_xTetra,
                "boundary tetra array", return 0);
  else
    PMMG_RECALLOC(mesh,mesh->xtetra,mesh->xtmax+1,xtmax_old+1,MMG5_xTetra,
                  "boundary tetra array", return 0);

  return ( PMMG_link_mesh( mesh ) );
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param np number of vertices.
 * \param ne number of tetrahedra.
 * \param nt number of triangles.
 * \param xp number of boundary point
 * \param xt number of boundary tetra
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Check the input mesh size and assign their values to the mesh.
 *
 */
int PMMG_setMeshSize_initData(MMG5_pMesh mesh, int np, int ne, int nt,
                              int xp, int xt ) {

  if ( ( (mesh->info.imprim > PMMG_VERB_DETQUAL) || mesh->info.ddebug ) &&
       ( mesh->point || mesh->xpoint || mesh->tetra || mesh->xtetra) )
    fprintf(stderr,"\n  ## Warning: %s: old mesh deletion.\n",__func__);

  if ( !np ) {
    fprintf(stderr,"  ** MISSING DATA:\n");
    fprintf(stderr,"     Your mesh must contains at least points.\n");
    return(0);
  }
  if ( !ne && (mesh->info.imprim > PMMG_VERB_DETQUAL || mesh->info.ddebug) ) {
    fprintf(stderr,"  ** WARNING:\n");
    fprintf(stderr,"     Your mesh don't contains tetrahedra.\n");
  }

  if ( mesh->point )
    MMG5_DEL_MEM(mesh,mesh->point);
  if ( mesh->tetra )
    MMG5_DEL_MEM(mesh,mesh->tetra);
  if ( mesh->prism )
    MMG5_DEL_MEM(mesh,mesh->prism);
  if ( mesh->tria )
    MMG5_DEL_MEM(mesh,mesh->tria);
  if ( mesh->quadra )
    MMG5_DEL_MEM(mesh,mesh->quadra);
  if ( mesh->edge )
    MMG5_DEL_MEM(mesh,mesh->edge);

  mesh->np  = np;
  mesh->ne  = ne;
  mesh->nt  = nt;
  mesh->xp  = xp;
  mesh->xt  = xt;

  mesh->npi = mesh->np;
  mesh->nei = mesh->ne;
  mesh->nti = mesh->nt;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param np number of vertices.
 * \param ne number of tetrahedra.
 * \param nt number of triangles.
 * \param xp number of boundary points.
 * \param xt number of boundary tetra.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Check the input mesh size and assign their values to the mesh.
 *
 */
int PMMG_setMeshSize(MMG5_pMesh mesh,int np,int ne,int nt,int xp,int xt ) {

  /* Check input data and set mesh->ne/na/np/nt to the suitable values */
  if ( !PMMG_setMeshSize_initData(mesh,np,ne,nt,xp,xt) )
    return 0;

  mesh->npmax  = mesh->np;
  mesh->nemax  = mesh->ne;
  mesh->ntmax  = mesh->nt;
  mesh->xpmax  = mesh->xp;
  mesh->xtmax  = mesh->xt;

  /* Mesh allocation and linkage */
  if ( !PMMG_setMeshSize_alloc( mesh ) ) return 0;

  return(1);

}

/**
 * \param parmesh pointer toward a parmesh.
 * \param grp pointer toward the grp that we want to fit.
 *
 * \return 1 if success, 0 if fail.
 *
 * Set the memMax field of a mesh to memCur and realloc the mesh in order to use
 * the less possible memory.
 *
 */
int PMMG_fitMeshSize( PMMG_pParMesh parmesh, PMMG_pGrp grp ) {
  const MMG5_pMesh mesh  = grp->mesh;
  const MMG5_pSol  met   = grp->met;
  const MMG5_pSol  disp  = grp->disp;
  const MMG5_pSol  ls    = grp->ls;
  const MMG5_pSol  field = grp->field;
  MMG5_pSol        psl;

  int npmax_old,xpmax_old,nemax_old,xtmax_old,is;
  int ier = 1;

  npmax_old = mesh->npmax;
  xpmax_old = mesh->xpmax;
  nemax_old = mesh->nemax;
  xtmax_old = mesh->xtmax;

  mesh->npmax = mesh->np;
  mesh->xpmax = mesh->xp;
  mesh->nemax = mesh->ne;
  mesh->xtmax = mesh->xt;

  if ( met ) {
    met->npmax = mesh->npmax;
  }
  if ( ls ) {
    ls->npmax   = mesh->npmax;
  }
  if ( disp ) {
    disp->npmax = mesh->npmax;
  }

  if ( mesh->nsols ) {
    assert ( field );
    for ( is=0; is<mesh->nsols; ++is ) {
      psl = &field[is];
      psl->npmax = mesh->npmax;
    }
  }

  if ( !PMMG_setMeshSize_realloc(mesh,npmax_old,xpmax_old,
                                 nemax_old,xtmax_old) ) ier = 0;

  if ( met && met->m ) {
    PMMG_REALLOC(mesh,met->m,met->size*(met->npmax+1),
                 met->size*(npmax_old+1),double,"metric_array",
                 ier = 0;);
  }
  if ( ls && ls->m ) {
    PMMG_REALLOC(mesh,ls->m,ls->size*(ls->npmax+1),
                 ls->size*(npmax_old+1),double,"ls_array",
                 ier = 0;);
  }
  if ( disp && disp->m ) {
    PMMG_REALLOC(mesh,disp->m,disp->size*(disp->npmax+1),
                 disp->size*(npmax_old+1),double,"disp_array",
                 ier = 0;);
  }
  if ( mesh->nsols ) {
    for ( is=0; is<mesh->nsols; ++is ) {
      psl = &field[is];
      if ( psl && psl->m ) {
        PMMG_REALLOC(mesh,psl->m,psl->size*(psl->npmax+1),
                     psl->size*(npmax_old+1),double,"field_array",
                     ier = 0;);
      }
    }
  }

  return ier;
}

/**
 * \param parmesh parmesh structure to adjust
 * \param fitMesh if 1, set maximum mesh size at its exact size.
 *
 * \return 1 if success, 0 if fail
 *
 * Update the size of the group meshes.
 *
 */
int PMMG_updateMeshSize( PMMG_pParMesh parmesh, int fitMesh )
{
  MMG5_pMesh mesh;
  MMG5_pSol  met,ls,disp,field,psl;
  size_t     available,used,delta;
  int        remaining_ngrps,npmax_old,xpmax_old,nemax_old,xtmax_old;
  int        i,is;

  for ( i = 0; i < parmesh->ngrp; ++i ) {
    mesh = parmesh->listgrp[i].mesh;

    /* Force the MMG5_memOption_memSet function to find the wanted memMax value
     * in MMG3D_Set_meshSize, MMG3D_Set_iparameter, MMG3D_zaldy (for i/o) */
    mesh->info.mem = parmesh->memGloMax/MMG5_MILLION;

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
      mesh->npmax = 1.5*mesh->np;
      mesh->xpmax = 1.5*mesh->xp;
      mesh->nemax = 1.5*mesh->ne;
      mesh->xtmax = 1.5*mesh->xt;
    }

    met = parmesh->listgrp[i].met;
    if ( !PMMG_memOption_memRepartition(mesh,met) ) return 0;

    if ( !PMMG_setMeshSize_realloc(mesh,npmax_old,xpmax_old,nemax_old,xtmax_old) )
      return 0;

    if ( met ) {
      met->np    = mesh->np;
      met->npmax = mesh->npmax;
    }

    ls = parmesh->listgrp[i].ls;
    if ( ls ) {
      ls->np      = mesh->np;
      ls->npmax   = mesh->npmax;
    }

    disp = parmesh->listgrp[i].disp;
    if ( disp ) {
      disp->np    = mesh->np;
      disp->npmax = mesh->npmax;
    }

    field = parmesh->listgrp[i].field;
    if ( mesh->nsols ) {
      assert ( field );
      for ( is=0; is<mesh->nsols; ++is ) {
        psl = &field[is];
        psl->np    = mesh->np;
        psl->npmax = mesh->npmax;
      }
    }

    if ( met && met->m )
      PMMG_REALLOC(mesh,met->m,met->size*(met->npmax+1),met->size*(npmax_old+1),
                   double,"metric array",return 0);

    if ( ls && ls->m ) {
      PMMG_REALLOC(mesh,ls->m,ls->size*(ls->npmax+1),
                   ls->size*(npmax_old+1),double,"ls_array",
                   return 0);
    }

    if ( disp && disp->m ) {
      PMMG_REALLOC(mesh,disp->m,disp->size*(disp->npmax+1),
                   disp->size*(npmax_old+1),double,"disp_array",
                   return 0);
    }

    if ( mesh->nsols ) {
      for ( is=0; is<mesh->nsols; ++is ) {
        psl = field + is;
        if ( psl && psl->m ) {
          PMMG_REALLOC(mesh,psl->m,psl->size*(psl->npmax+1),
                       psl->size*(npmax_old+1),double,"field_array",
                       return 0);
        }
      }
    }

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
      PMMG_DEL_MEM ( parmesh,(*ext_comm+k)->rtorecv,double,"rtorecv" );

  }

  PMMG_REALLOC(parmesh,*ext_comm,newSize,*oldSize,PMMG_Ext_comm,"ext_comm",return 0);

  if ( newSize > *oldSize )
    memset( *ext_comm + *oldSize, 0x0,
            ((size_t)((newSize)-(*oldSize)))*sizeof(PMMG_Ext_comm));

  *oldSize = newSize;

  return 1;
}

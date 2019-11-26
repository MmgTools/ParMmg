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
 * \file distributemesh.c
 * \brief Distribute the mesh on the processors.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */
#include "parmmg.h"
#include "mpitypes_pmmg.h"
#include "metis_pmmg.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param start index of the starting tetrahedra.
 * \param ip index of the point in the mesh.
 * \param iploc local index of the point in the tetrahedra \a start.
 * \param visited integer array, equal to ip for each tetrahedron already
 * visited when looking for node ip.
 * \param list pointer toward the list of the tetra in the volumic ball of
 * \a ip.
 * \return 0 if fail and the number of the tetra in the ball otherwise.
 *
 * Fill the volumic ball (i.e. filled with tetrahedra) of point \a ip in tetra
 * \a start. Results are stored under the form \f$4*kel + jel\f$, kel = number
 * of the tetra, jel = local index of p within kel.
 *
 */
int PMMG_boulevolp (MMG5_pMesh mesh, int start, int ip, int iploc, int *visited,
                    int *list){
  MMG5_pTetra  pt1;
  int    *adja,ilist,cur,k,k1;
  int8_t  i,j,l;

  /* Store initial tetrahedron */
  visited[start] = ip;
  list[0] = 4*start + iploc;
  ilist=1;

  /* Explore list and travel by adjacency through elements sharing p */
  cur = 0;
  while ( cur < ilist ) {
    k = list[cur] / 4;
    i = list[cur] % 4; // index of point p in tetra k
    adja = &mesh->adja[4*(k-1)+1];

    for (l=0; l<3; l++) {
      i  = MMG5_inxt3[i];
      k1 = adja[i];
      if ( !k1 )  continue;
      k1 /= 4;
      pt1 = &mesh->tetra[k1];
      if ( visited[k1] == ip )  continue;
      visited[k1] = ip;
      for (j=0; j<4; j++)
        if ( pt1->v[j] == ip )  break;
      assert(j<4);
      /* overflow */
      if ( ilist > MMG3D_LMAX-3 )  return 0;
      list[ilist] = 4*k1+j;
      ilist++;
    }
    cur++;
  }
  return ilist;
}

/**
 * \param xtetra pointer toward a table containing the xtetra structures.
 * \param *perm pointer toward the permutation table (to perform in place
 * permutations).
 * \param ind1 index of the first xtetra to swap.
 * \param ind2 index of the second xtetra to swap.
 *
 * Swap two xtetra in the table of xtetrahedras.
 */
static void PMMG_swapxTetra( MMG5_pxTetra xtetra, int* perm, int ind1, int ind2 )
{
  MMG5_xTetra pxttmp;
  int         tmp;

  /** 1- swap the xtetra */
  memcpy(&pxttmp      ,&xtetra[ind2],sizeof(MMG5_xTetra));
  memcpy(&xtetra[ind2],&xtetra[ind1],sizeof(MMG5_xTetra));
  memcpy(&xtetra[ind1],&pxttmp      ,sizeof(MMG5_xTetra));

  /** 2- swap the permutation table */
  tmp        = perm[ind2];
  perm[ind2] = perm[ind1];
  perm[ind1] = tmp;
}

/**
 * \param xpoint pointer toward a table containing the xpoint structures.
 * \param *perm pointer toward the permutation table (to perform in place
 * permutations).
 * \param ind1 index of the first xpoint to swap.
 * \param ind2 index of the second xpoint to swap.
 *
 * Swap two xpoint in the table of xpoints.
 */
static void PMMG_swapxPoint( MMG5_pxPoint xpoint, int* perm, int ind1, int ind2 )
{
  MMG5_xPoint pxptmp;
  int         tmp;

  /** 1- swap the xpoint */
  memcpy(&pxptmp      ,&xpoint[ind2],sizeof(MMG5_xPoint));
  memcpy(&xpoint[ind2],&xpoint[ind1],sizeof(MMG5_xPoint));
  memcpy(&xpoint[ind1],&pxptmp      ,sizeof(MMG5_xPoint));

  /** 2- swap the permutation table */
  tmp        = perm[ind2];
  perm[ind2] = perm[ind1];
  perm[ind1] = tmp;
}

/**
 * \param point pointer toward a table containing the point structures.
 * \param met pointer toward a table containing the metric structure.
 * \param *perm pointer toward the permutation table (to perform in place
 * permutations).
 * \param ind1 index of the first xpoint to swap.
 * \param ind2 index of the second xpoint to swap.
 * \param metsize size of the metric (1=iso,6=aniso).
 *
 * Swap two points in the table of points.
 */
static void PMMG_swapPoint( MMG5_pPoint point, double* met,int* perm,
                            int ind1, int ind2, int metsiz )
{
  MMG5_Point ppttmp;
  MMG5_Sol   mettmp;
  int        tmp,addr2,addr1;

  /** 1- swap the xpoint */
  memcpy(&ppttmp      ,&point[ind2], sizeof(MMG5_Point));
  memcpy(&point[ind2] ,&point[ind1], sizeof(MMG5_Point));
  memcpy(&point[ind1] ,&ppttmp      ,sizeof(MMG5_Point));

  /** 2- swap the mets */
  if ( met ) {
    addr1 = ind1*metsiz;
    addr2 = ind2*metsiz;
    memcpy(&mettmp    ,&met[addr2],metsiz*sizeof(double));
    memcpy(&met[addr2],&met[addr1],metsiz*sizeof(double));
    memcpy(&met[addr1],&mettmp    ,metsiz*sizeof(double));
  }

  /** 3- swap the permutation table */
  tmp        = perm[ind2];
  perm[ind2] = perm[ind1];
  perm[ind1] = tmp;
}

/**
 * \param parmesh pointer toward a parmesh structure.
 *
 * \return 0 on all procs if fail
 *         1 on all procs if success
 *
 * Send the initial mesh (with tria and edges) from proc 0 toward the other
 * procs.
 *
 */
int PMMG_bcast_mesh( PMMG_pParMesh parmesh )
{
  PMMG_pGrp    grp;
  MMG5_pMesh   mesh;
  MMG5_pSol    met;
  MPI_Datatype mpi_light_point, mpi_light_tetra, mpi_tria,mpi_edge;
  int          k,rank,root,ier,ieresult,isMet;

  /** Proc 0 send the mesh to the other procs */
  grp    = &parmesh->listgrp[0];
  mesh   = grp->mesh;
  met    = grp->met;
  rank   = parmesh->myrank;
  root   = parmesh->info.root;
  isMet  = met->m ? 1 : 0;

  /* Minimize the memory used by the mesh */
  ier = 1;
  if ( rank == root ) {
    PMMG_RECALLOC(mesh,mesh->point,mesh->np+1,mesh->npmax+1,MMG5_Point,
                  "vertices array", ier = 6);
    if ( ier ) {
      mesh->npmax = mesh->np;
    }
    if ( isMet ) {
      PMMG_RECALLOC(mesh,met->m,met->size*(met->np+1),met->size*(met->npmax+1),double,
                    "metric array", ier = 6);
    }
    if ( ier ) {
      met->npmax = met->np;
    }

    PMMG_RECALLOC(mesh,mesh->tetra,mesh->ne+1,mesh->nemax+1,MMG5_Tetra,
                  "tetra array", ier = 6);
    if ( ier ) {
      mesh->nemax = mesh->ne;
    }
  }

  /* Mesh */
  MPI_CHECK( MPI_Bcast( &mesh->np,     1, MPI_INT,       root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->ne,     1, MPI_INT,       root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->nt,     1, MPI_INT,       root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->na,     1, MPI_INT,       root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->ntmax,  1, MPI_INT,       root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->memMax, 1, MPI_LONG_LONG, root, parmesh->comm ), ier=6);
  /* Metric */
  MPI_CHECK( MPI_Bcast( &met->size,  1, MPI_INT, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &met->type,  1, MPI_INT, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &met->npmax, 1, MPI_INT, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &met->np,    1, MPI_INT, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &isMet,      1, MPI_INT, root, parmesh->comm ), ier=6);
  /* Info */
  MPI_CHECK( MPI_Bcast( &mesh->info.dhd,       1, MPI_DOUBLE, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.hmin,      1, MPI_DOUBLE, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.hmax,      1, MPI_DOUBLE, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.hsiz,      1, MPI_DOUBLE, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.hgrad,     1, MPI_DOUBLE, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.hgradreq,  1, MPI_DOUBLE, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.hausd,     1, MPI_DOUBLE, root, parmesh->comm ), ier=6);

  MPI_CHECK( MPI_Bcast( &mesh->info.delta,     1, MPI_DOUBLE, root, parmesh->comm ), ier=6); 
  MPI_CHECK( MPI_Bcast( &mesh->info.min,       3, MPI_DOUBLE, root, parmesh->comm ), ier=6);

  MPI_CHECK( MPI_Bcast( &mesh->info.ls,        1, MPI_DOUBLE, root, parmesh->comm ), ier=6);

  MPI_CHECK( MPI_Bcast( &mesh->info.npar,      1, MPI_INT, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.opnbdy,    1, MPI_INT, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.renum,     1, MPI_INT, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.PROctree,  1, MPI_INT, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.nmat,      1, MPI_INT, root, parmesh->comm ), ier=6);

  MPI_CHECK( MPI_Bcast( &mesh->info.nreg,      1, MPI_CHAR, root, parmesh->comm ), ier=6); 
  MPI_CHECK( MPI_Bcast( &mesh->info.imprim,    1, MPI_CHAR, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.ddebug,    1, MPI_CHAR, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.iso,       1, MPI_CHAR, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.lag,       1, MPI_CHAR, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.parTyp,    1, MPI_CHAR, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.optim,     1, MPI_CHAR, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.optimLES,  1, MPI_CHAR, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.noinsert,  1, MPI_CHAR, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.noswap,    1, MPI_CHAR, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.nomove,    1, MPI_CHAR, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.nosurf,    1, MPI_CHAR, root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.inputMet,  1, MPI_CHAR, root, parmesh->comm ), ier=6);

  /* affectation of old refs in ls-mode */
  if ( mesh->info.nmat ) {
    
    if( rank != root )
      MMG5_SAFE_CALLOC(mesh->info.mat,mesh->info.nmat,MMG5_Mat, ier = 0);

    if ( ier ) {
      for ( k=0; k<mesh->info.nmat; ++k ) {
        MPI_CHECK( MPI_Bcast( &mesh->info.mat[k].dospl, 1, MPI_CHAR, root, parmesh->comm ), ier=6);
        MPI_CHECK( MPI_Bcast( &mesh->info.mat[k].ref,   1, MPI_INT, root, parmesh->comm ), ier=6);
        MPI_CHECK( MPI_Bcast( &mesh->info.mat[k].rin,   1, MPI_INT, root, parmesh->comm ), ier=6);
        MPI_CHECK( MPI_Bcast( &mesh->info.mat[k].rex,   1, MPI_INT, root, parmesh->comm ), ier=6);
      }
    }
  }

  /* local parameters */
  if ( mesh->info.npar ) {

    if( rank != root )
      MMG5_SAFE_CALLOC(mesh->info.par,mesh->info.npar,MMG5_Par, ier = 0);

    if ( ier ) {
      for ( k=0; k<mesh->info.npar; ++k ) {
        MPI_CHECK( MPI_Bcast( &mesh->info.par[k].hmin,    1, MPI_DOUBLE, root, parmesh->comm ), ier=6);
        MPI_CHECK( MPI_Bcast( &mesh->info.par[k].hmax,    1, MPI_DOUBLE, root, parmesh->comm ), ier=6);
        MPI_CHECK( MPI_Bcast( &mesh->info.par[k].hausd,   1, MPI_DOUBLE, root, parmesh->comm ), ier=6);
        MPI_CHECK( MPI_Bcast( &mesh->info.par[k].ref,     1, MPI_INT, root, parmesh->comm ), ier=6);
        MPI_CHECK( MPI_Bcast( &mesh->info.par[k].elt,     1, MPI_CHAR, root, parmesh->comm ), ier=6);
      }
    }
  }

  mesh->npmax = mesh->npi = mesh->np;
  mesh->npnil = 0;
  mesh->nemax = mesh->nei = mesh->ne;
  mesh->nenil = 0;
  mesh->ntmax = mesh->nti = mesh->nt;
  mesh->xtmax = mesh->ntmax;

  met->npmax = met->npi = met->np;
  met->ver = mesh->ver;
  met->dim = mesh->dim;

  if ( rank != root ) {
    PMMG_CALLOC(mesh,mesh->point,mesh->npmax+1,MMG5_Point,"initial vertices"  , ier=6);
    PMMG_CALLOC(mesh,mesh->tetra,mesh->nemax+1,MMG5_Tetra,"initial tetrahedra", ier=6);

    if ( mesh->nt )
      PMMG_CALLOC(mesh,mesh->tria,mesh->nt+1,MMG5_Tria,"initial triangles", ier=6);

    if ( mesh->na )
      PMMG_CALLOC(mesh,mesh->edge,mesh->na+1,MMG5_Edge,"initial edges", ier=6);

    if ( isMet )
      PMMG_CALLOC(mesh,met->m,met->size*(met->npmax+1),double,"initial metric", ier=6);
  }

  if ( ier<6 && !PMMG_create_MPI_lightPoint( &mpi_light_point ) ) { ier=6; }
  if ( ier<6 && !PMMG_create_MPI_lightTetra( &mpi_light_tetra ) ) { ier=5; }
  if ( ier<5 && mesh->nt && !PMMG_create_MPI_Tria( &mpi_tria ) )   { ier=4; }
  if ( ier<4 && mesh->na && !PMMG_create_MPI_Edge( &mpi_edge ) )   { ier=3; }

  MPI_CHECK ( MPI_Allreduce( &ier,&ieresult,1,MPI_INT,MPI_MAX,parmesh->comm ),
              ier = ieresult = 2 );

  if ( ieresult<2 ) {

    MPI_CHECK( MPI_Bcast(mesh->point,mesh->np+1,mpi_light_point,root,parmesh->comm), ier=2);
    MPI_CHECK( MPI_Bcast(mesh->tetra,mesh->ne+1,mpi_light_tetra,root,parmesh->comm), ier=2);
    if ( mesh->nt )
      MPI_CHECK( MPI_Bcast( mesh->tria,mesh->nt+1,mpi_tria,root,parmesh->comm ), ier=2);
    if ( mesh->na )
      MPI_CHECK( MPI_Bcast( mesh->edge,mesh->na+1,mpi_edge,root,parmesh->comm ), ier=2);
    if ( met->m )
      MPI_CHECK( MPI_Bcast(met->m,met->size*(met->npmax+1),MPI_DOUBLE,root,
                           parmesh->comm ), ier=2);
  }

  /* Deallocations */
  if ( ier < 6 ) {
    MPI_Type_free( &mpi_light_point );
    if ( ier < 5 ) {
      MPI_Type_free( &mpi_light_tetra );
      if ( ier < 4 ) {
        if ( mesh->nt )
          MPI_Type_free( &mpi_tria );
        if ( ier < 3 ) {
          if ( mesh->na )
            MPI_Type_free( &mpi_edge );
        }
      }
    }
  }

  MPI_CHECK ( MPI_Allreduce( &ier,&ieresult,1,MPI_INT,MPI_MAX,parmesh->comm ),
              ieresult = 2 );

  return ieresult==1;

}

/**
 * \param mesh pointer toward a MMG5 mesh structure
 * \param rank MPI rank
 *
 * \return 1 if success, 0 if fail
 *
 * Pack the tetrahedra and remove those ones that are not on the processor (the
 * proc index is stored in pt->mark).
 *
 */
static inline
int PMMG_packTetraOnProc(MMG5_pMesh mesh, int rank) {
  MMG5_pTetra pt,ptnew;
  int         ne,nbl,k,iadr;

  ne  = 0;
  nbl = 1;
  for ( k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];

    if ( (!MG_EOK(pt)) || (pt->mark != rank) ) continue;

    pt->v[0] = mesh->point[pt->v[0]].tmp;
    pt->v[1] = mesh->point[pt->v[1]].tmp;
    pt->v[2] = mesh->point[pt->v[2]].tmp;
    pt->v[3] = mesh->point[pt->v[3]].tmp;

    ++ne;

    if ( k!=nbl ) {
      assert( nbl==pt->flag && "non-matching index" );
      ptnew = &mesh->tetra[nbl];
      memcpy(ptnew,pt,sizeof(MMG5_Tetra));
    }
    ++nbl;
  }
  mesh->ne = mesh->nei = ne;

  /* Treat empty tetra (prepare to link mesh) */
  for( k=mesh->ne+1; k<= mesh->nemax; k++ ) {
    pt = &mesh->tetra[k];
    memset(pt,0,sizeof(MMG5_Tetra));
    iadr = 4*(k-1) + 1;
    if ( mesh->adja )
      memset(&mesh->adja[iadr],0,4*sizeof(int));
  }

  return 1;
}

/**
 * \param mesh pointer toward a MMG5 mesh structure
 * \param met pointer toward a MMG5 solution structure
 * \param pointPerm array of new point positions
 * \param xPointPerm array of new xPoint positions
 * \param xTetraPerm array of new xTetra positions
 *
 * \return 1 if success, 0 if fail
 *
 * Permute the point and boundary entities with respect to the provided
 * permutation arrays
 *
 */
static inline
int PMMG_permuteMesh(MMG5_pMesh mesh,MMG5_pSol met,
                     int *pointPerm,int *xPointPerm,int *xTetraPerm,
                     int np,int nxp,int nxt) {
  int k;

  /** Compact xtetra on the proc: in place permutations */
  for ( k=1; k<=mesh->xt; ++k )
    while ( xTetraPerm[k] != k && xTetraPerm[k] )
      PMMG_swapxTetra(mesh->xtetra,xTetraPerm,k,xTetraPerm[k]);

  /** Compact vertices on the proc: in place permutations */
  for ( k=1; k<=mesh->np; ++k )
    while ( pointPerm[k] != k && pointPerm[k] )
      PMMG_swapPoint(mesh->point,met->m,pointPerm,k,pointPerm[k],met->size);

  /** Compact xpoint on the proc: in place permutations */
  for ( k=1; k<=mesh->xp; ++k )
    while ( xPointPerm[k] != k && xPointPerm[k] )
      PMMG_swapxPoint(mesh->xpoint,xPointPerm,k,xPointPerm[k]);

  mesh->np = mesh->npi = np;
  met->np  = met->npi  = np;
  mesh->xp = nxp;
  mesh->xt = nxt;

  /* Treat empty points (prepare to link mesh) */
  for( k=mesh->np+1;k<=mesh->npmax;k++ ) {
    /* Set tangent field of point to 0 */
    mesh->point[k].n[0] = 0;
    mesh->point[k].n[1] = 0;
    mesh->point[k].n[2] = 0;
  }
  mesh->point[mesh->npmax-1].tmp = 0;
  mesh->point[mesh->npmax].tmp   = 0;

  return 1;
}

/**
 * \param parmesh pointer toward a PMMG parmesh structure.
 * \param part pointer toward the metis array containing the partitions.
 * \param mesh pointer a MMG5 mesh structure.
 * \param np  number of point inside the local mesh.
 * \param nxp number of boundary points inside the local mesh.
 * \param nxt number of boundary tetra inside the local mesh.
 * \param pointPerm array of new point positions.
 * \param xPointPerm array of new xPoint positions.
 * \param xTetraPerm array of new xTetra positions.
 * \param shared_pt  array of the number of points shared with each other procs.
 * \param shared_face array of the number of faces shared with other procs.
 * \param seen_shared_pt binary array that contains 1 if a point is shared with
 * a given proc
 *
 * \return 0 if fail, 1 if success
 *
 * Mark the mesh entities that will stay on the processor, numbered it and
 * create the permutation arrays to obtain the final mesh. Count and mark the
 * nodes and faces shared with other processors.
 *
 */
static inline
int PMMG_mark_localMesh(PMMG_pParMesh parmesh,idx_t *part,MMG5_pMesh mesh,
                        int *np,int *nxp,int *nxt,
                        int **pointPerm,int **xPointPerm,int **xTetraPerm,
                        int **shared_pt,int **shared_face,
                        int8_t **seen_shared_pt) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_pPoint  ppt;
  int          nprocs,rank,rankVois,ret_val,k,kvois,j,ip,iploc,l,ne,newsize,
               list[MMG3D_LMAX+2],ilist,k1,rankVois1,cur,*flag;
  int8_t       ifac;

  ret_val = 1;

  rank   = parmesh->myrank;
  nprocs = parmesh->nprocs;

  /** Remove the part of the mesh that are not on the proc rank */
  (*shared_pt) = (*shared_face) = NULL;
  PMMG_CALLOC(parmesh,(*shared_pt),nprocs,int,"shared_pt array",
              ret_val = 0;goto fail_alloc1);
  PMMG_CALLOC(parmesh,(*seen_shared_pt),nprocs*mesh->np,int8_t,"seen_shared_pt array",
              ret_val = 0; goto fail_alloc2);
  PMMG_CALLOC(parmesh,(*shared_face),nprocs,int,"shared_face array",
              ret_val = PMMG_FAILURE;goto fail_alloc3);
  PMMG_CALLOC(parmesh,(*pointPerm),mesh->np+1,int,"pointPerm",
              ret_val = 0;goto fail_alloc4);
  PMMG_CALLOC(parmesh,(*xTetraPerm),mesh->xtmax+1,int,"xTetraPerm",
              ret_val = 0; goto fail_alloc5);
  PMMG_CALLOC(parmesh,(*xPointPerm),mesh->xpmax+1,int,"xPointPerm",
              ret_val = 0;goto fail_alloc6);

  (*nxp) = 0;
  (*nxt) = 0;
  (*np)  = 0;

  /** Visited tetra flag array */
  PMMG_CALLOC(parmesh,flag,mesh->ne+1,int,"Visited tetra flag array",return 0;);
  for ( k=1; k<=mesh->ne; k++ )
    flag[k] = 0;

  /* Reset the tmp field of points (it will be used to store the local index of
   * the point on each proc) */
  for ( k=1; k<=mesh->np; k++ )
    mesh->point[k].tmp = 0;

  /* Reset the flag field of tetras (it will be used to store the local index of
   * the tetra on each proc) */
  ne = 0;
  for ( k=1; k<=mesh->ne; k++ )
    mesh->tetra[k].flag = 0;

  /** Mark mesh entities that will stay on the proc and count the number of
   * point that must be communicated to the other procs  */
  for ( k=1; k<=mesh->ne; k++ ) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    pt->mark = part[k-1];
    if ( pt->mark != rank ) continue;

    pt->flag = ++ne;

    for ( ifac=0; ifac<4; ifac++ ) {
      kvois = mesh->adja[4*k-3+ifac]/4;
      if ( kvois )
        rankVois = part[kvois-1];
      else
        rankVois = rank;

      /* Mark the interfaces between two procs */
      if ( rank != rankVois ) {
        ++(*shared_face)[rankVois];
        if ( !pt->xt ) {
          if ( (mesh->xt + 1) > mesh->xtmax ) {
            /* realloc of xtetras table */
            newsize = MG_MAX((1.+mesh->gap)*mesh->xtmax,mesh->xtmax+1);
            PMMG_RECALLOC(mesh,mesh->xtetra,newsize+1,
                          mesh->xtmax+1,MMG5_xTetra,"larger xtetra ",
                          ret_val = 0;goto fail_alloc7);
            PMMG_RECALLOC(parmesh,(*xTetraPerm),newsize+1,
                          mesh->xtmax+1,int,"larger xtetra permutation table ",
                          ret_val = 0; goto fail_alloc7);
            mesh->xtmax = newsize;
          }
          ++mesh->xt;
          pt->xt = mesh->xt;
        }
        pxt = &mesh->xtetra[pt->xt];
        /* Parallel face (if already boundary, make it recognizable as a true
         * boundary) */
        if ( pxt->ftag[ifac] & MG_BDY ) pxt->ftag[ifac] |= MG_PARBDYBDY;
        pxt->ftag[ifac] |= (MG_PARBDY + MG_BDY + MG_REQ + MG_NOSURF);

        /* Parallel edges */
        for ( j=0; j<3; ++j )
          pxt->tag[MMG5_iarf[ifac][j]] |= (MG_PARBDY + MG_BDY + MG_REQ + MG_NOSURF);
      }

      for ( j=0; j<3; ++j ) {
        iploc = MMG5_idir[ifac][j];
        ip    = pt->v[iploc];
        ppt   = &mesh->point[ip];

        if ( rankVois != rank ) {
          /* Check if the node is also shared among procs which are not
           * communicating through an element face ==> Build the ball of the
           * node and scan the elements partition number, start from kvois. */
          flag[k] = ip;
          for (l=0; l<4; l++)
            if ( mesh->tetra[kvois].v[l] == ip )  break;
          ilist = PMMG_boulevolp(mesh,kvois,ip,l,flag,list);

          for (cur = 0; cur<ilist; cur++ ) {
            k1 = list[cur] / 4;
            rankVois1 = part[k1-1];
            if ( rankVois1 == rank ) continue; // Skip if not to be communicated

            /* Count (and mark) each point that will be shared between me and the
             * proc rankVois */
            if ( !(*seen_shared_pt)[nprocs*(ip-1)+rankVois1] ) {
              (*seen_shared_pt)[nprocs*(ip-1)+rankVois1] = 1;
              ++(*shared_pt)[rankVois1];

              /* Mark parallel vertex */
              ppt->tag |= (MG_PARBDY + MG_BDY + MG_REQ + MG_NOSURF);

// TO REMOVE WHEN MMG WILL BE READY
              if ( !ppt->xp ) {
                if ( (mesh->xp+1) > mesh->xpmax ) {
                  /* realloc of xtetras table */
                  newsize = MG_MAX((1.+mesh->gap)*mesh->xpmax,mesh->xpmax+1);
                  PMMG_RECALLOC(mesh,mesh->xpoint,newsize+1,
                                mesh->xpmax+1,MMG5_xPoint,"larger xpoint ",
                                ret_val = 0;goto fail_alloc7);
                  PMMG_RECALLOC(parmesh,(*xPointPerm),newsize+1,
                                mesh->xpmax+1,int,"larger xpoint permutation table ",
                                ret_val = 0; goto fail_alloc7);
                  mesh->xpmax = newsize;
                }
                ++mesh->xp;
                ppt->xp = mesh->xp;
                if ( ppt->tmp ) {
                  (*xPointPerm)[ppt->xp] = ++(*nxp);
                  ppt->xp                = (*nxp);
                }
              }
// END TO REMOVE WHEN MMG WILL BE READY
            }
          }
        }

        if ( !ppt->tmp ) {
          /* Mark the new point index and update the table of permutation */
          ppt->tmp         = ++(*np);
          (*pointPerm)[ip] = (*np);

          /* update the table of permutation for xPoint if needed */
          if ( ppt->xp ) {
            (*xPointPerm)[ppt->xp] = ++(*nxp);
            ppt->xp                = (*nxp);
          }
        }
      }
    }

    /* update the table of permutation for xTetra if needed */
    if ( !pt->xt ) continue;

    (*xTetraPerm)[pt->xt] = ++(*nxt);
    pt->xt = (*nxt);
  }

  PMMG_DEL_MEM(parmesh,flag,int,"deallocate visited tetra flag array");

fail_alloc1:
  return ret_val;

fail_alloc7:
  PMMG_DEL_MEM(parmesh,*xPointPerm,int,"deallocate xPointPerm");
fail_alloc6:
  PMMG_DEL_MEM(parmesh,*xTetraPerm,int,"deallocate xTetraPerm");
fail_alloc5:
  PMMG_DEL_MEM(parmesh,*pointPerm,int,"deallocate pointPerm");
fail_alloc4:
  PMMG_DEL_MEM(parmesh,*shared_face,int,"deallocate shared_face");
fail_alloc3:
  PMMG_DEL_MEM(parmesh,*seen_shared_pt,int8_t,"deallocate seen_shared_pt");
fail_alloc2:
  PMMG_DEL_MEM(parmesh,*shared_pt,int,"deallocate shared_pt");
  return ret_val;
}

static inline
int PMMG_create_empty_communicators( PMMG_pParMesh parmesh ) {
  PMMG_pGrp       grp;

  grp    = &parmesh->listgrp[0];

  /** Internal communicators allocation */
  parmesh->next_node_comm = 0;
  parmesh->next_face_comm = 0;
  grp->nitem_int_node_comm = 0;
  grp->nitem_int_face_comm = 0;

  PMMG_CALLOC(parmesh,parmesh->int_node_comm,1,PMMG_Int_comm,
              "allocating int_node_comm",return 0);
  PMMG_CALLOC(parmesh,parmesh->int_face_comm,1,PMMG_Int_comm,
              "allocating int_face_comm",return 0);
  parmesh->int_node_comm->nitem = 0;
  parmesh->int_face_comm->nitem = 0;

  return 1;
}

/**
 * \param parmesh pointer toward a PMMG parmesh structure.
 * \param part pointer toward the metis array containing the partitions.
 * \param shared_pt pointer toward the array of the number of points shared
 * with each other procs.
 * \param shared_face pointer toward the array of the number of faces shared
 * \param seen_shared_pt pointer toward a binary array that contains 1 if a
 * point is shared with a given proc
 *
 * \return 0 if fail, 1 if success
 *
 * Create the face and node communicators when distributing a mesh over multiple
 * processors
 *
 * \warning the parmesh structure must contains only 1 group.
 *
 */
static inline
int PMMG_create_communicators(PMMG_pParMesh parmesh,idx_t *part,int *shared_pt,
                              int *shared_face,int8_t *seen_shared_pt) {
  PMMG_pGrp       grp;
  MMG5_pMesh      mesh;
  MMG5_pTetra     pt,ptvois;
  PMMG_pExt_comm  pext_node_comm,pext_face_comm;
  int             rank,nprocs,rankCur,rankVois;
  int             next_node_comm,next_face_comm;
  int             nitem_int_node_comm,nitem_int_face_comm;
  int            *node2int_node_comm_index1,*node2int_node_comm_index2;
  int            *face2int_face_comm_index1,*face2int_face_comm_index2;
  int             inIntComm;
  int            *idx,kvois,k,i,j,ip,iploc, iplocvois;
  int8_t          ifac,ifacvois;

  rank   = parmesh->myrank;
  nprocs = parmesh->nprocs;
  grp    = &parmesh->listgrp[0];
  mesh   = grp->mesh;

  /** Count the number of external node/face communicators and initialize it */
  next_node_comm = next_face_comm = 0;
  for ( k=0; k<nprocs; ++k ) {
    if ( shared_pt[k] )   ++next_node_comm;
    if ( shared_face[k] ) ++next_face_comm;
  }

  PMMG_CALLOC(parmesh,parmesh->ext_node_comm,next_node_comm,PMMG_Ext_comm,
              "allocate ext_node_comm ",return 0);
  parmesh->next_node_comm = next_node_comm;
  PMMG_CALLOC(parmesh,parmesh->ext_face_comm,next_face_comm,PMMG_Ext_comm,
              "allocate ext_face_comm ",return 0);
  parmesh->next_face_comm = next_face_comm;

  /** Count internal communicators items before erasing the shared_face array */
  /* Internal node comm */
  nitem_int_node_comm = 0;
  for ( k=1; k<=mesh->np; k++ ) {
    if ( !mesh->point[k].tmp )  continue;

    for ( j=0; j<nprocs; ++j ) {
      if ( seen_shared_pt[nprocs*(k-1)+j] ) {
        ++nitem_int_node_comm;
        break;
      }
    }
  }
  /* Internal face comm */
  nitem_int_face_comm = 0;
  for ( k=0; k<parmesh->nprocs; k++ ) {
    nitem_int_face_comm += shared_face[k];
  }

  /** External communicators allocation */
  next_node_comm = next_face_comm = 0;
  for ( k=0; k<nprocs; ++k ) {
    /* Node ext comm */
    if ( shared_pt[k] ) {
      pext_node_comm = &parmesh->ext_node_comm[next_node_comm];
      pext_node_comm->color_in  = rank;
      pext_node_comm->color_out = k;
      PMMG_CALLOC(parmesh,pext_node_comm->int_comm_index,shared_pt[k],int,
                  "allocate comm idx",return 0);
      pext_node_comm->nitem     = shared_pt[k];
      /* Use shared_pt to store the idx of the external communicator me->k */
      shared_pt[k] = next_node_comm++;
    }

    /* Face ext comm */
    if ( shared_face[k] ) {
      pext_face_comm            = &parmesh->ext_face_comm[next_face_comm];
      pext_face_comm->color_in  = rank;
      pext_face_comm->color_out = k;
      PMMG_CALLOC(parmesh,pext_face_comm->int_comm_index,shared_face[k],int,
                  "allocate comm idx",return 0);
      pext_face_comm->nitem     = shared_face[k];
      /* Use shared_pt to store the idx of the external communicator me->k */
      shared_face[k] = next_face_comm++;
    }
  }

  /** Internal communicators allocation */
  /* Internal node comm */
  assert ( !grp->nitem_int_node_comm );
  PMMG_CALLOC(parmesh,grp->node2int_node_comm_index1,nitem_int_node_comm,int,
              "node2int_node_comm_index1 ",return 0);
  grp->nitem_int_node_comm = nitem_int_node_comm;

  PMMG_CALLOC(parmesh,grp->node2int_node_comm_index2, nitem_int_node_comm,int,
              "alloc node2int_node_comm_index2 ",
              PMMG_DEL_MEM(parmesh,grp->node2int_node_comm_index1,int,
                           "free node2int_node_comm_index1 ");
              grp->nitem_int_node_comm=0; return 0);

  /* Internal face comm */
  assert ( !grp->nitem_int_face_comm );
  grp->nitem_int_face_comm = nitem_int_face_comm;
  PMMG_CALLOC(parmesh,grp->face2int_face_comm_index1,nitem_int_face_comm,int,
              "alloc face2int_face_comm_index1 ",return 0);

  PMMG_CALLOC(parmesh,grp->face2int_face_comm_index2,nitem_int_face_comm,int,
              "alloc face2int_face_comm_index2 ",
              PMMG_DEL_MEM(parmesh,grp->face2int_face_comm_index1,int,
                           "free face2int_face_comm_index1 ");
              grp->nitem_int_face_comm = 0;return 0);

  /** Travel through the mesh and fill the communicators */
  node2int_node_comm_index1 = grp->node2int_node_comm_index1;
  node2int_node_comm_index2 = grp->node2int_node_comm_index2;
  face2int_face_comm_index1 = grp->face2int_face_comm_index1;
  face2int_face_comm_index2 = grp->face2int_face_comm_index2;

  /* Node Communicators */
  /* Idx is used to store the external communicator cursor */
  PMMG_CALLOC(parmesh,idx,parmesh->next_node_comm,int,"allocating idx",return 0);

  i = 0;
  for ( k=1; k<=mesh->np; k++ ) {
    if ( !mesh->point[k].tmp ) continue;

    inIntComm = 0;
    for ( j=0; j<nprocs; ++j ) {
      pext_node_comm = &parmesh->ext_node_comm[shared_pt[j]];

      if ( seen_shared_pt[nprocs*(k-1)+j] == 1 ) {
        /* Add point in external communicator */
        pext_node_comm->int_comm_index[idx[shared_pt[j]]++] = i;

        if ( !inIntComm ) {
          /* Add point in internal communicator */
          inIntComm = 1;
          node2int_node_comm_index1[i] = mesh->point[k].tmp;
          node2int_node_comm_index2[i] = i;
        }
      }
    }
    /* Increment internal comm cursor */
    if ( inIntComm )
      ++i;
  }
  PMMG_DEL_MEM(parmesh,idx,int,"deallocating idx");

  /* Face Communicators */
  /* Idx is used to store the external communicator cursor */
  PMMG_CALLOC(parmesh,idx,parmesh->next_face_comm,int,"allocating idx",return 0);

  i = 0;
  for ( k=1; k<=mesh->ne; k++ ) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    rankCur = pt->mark;

    for ( ifac=0; ifac<4; ifac++ ) {
      kvois = mesh->adja[4*k-3+ifac]/4;
      if ( (!kvois) || k>kvois ) continue;

      rankVois = part[kvois-1];
      if ( rankCur == rankVois ) continue;

      ptvois = &mesh->tetra[kvois];

      /* Find a common starting point inside the face for both tetra */
      iploc = 0; // We impose the starting point in k
      ip    = pt->v[MMG5_idir[ifac][iploc]];

      ifacvois = mesh->adja[4*k-3+ifac]%4;
      for ( iplocvois=0; iplocvois < 3; ++iplocvois )
        if ( ptvois->v[MMG5_idir[ifacvois][iplocvois]] == ip ) break;
      assert ( iplocvois < 3 );

      if ( rankCur == rank ) {
        /* Add the elt k to communicators */
        pext_face_comm = &parmesh->ext_face_comm[shared_face[rankVois]];
        pext_face_comm->int_comm_index[idx[shared_face[rankVois]]++] = i;

        face2int_face_comm_index1[i] = 12*pt->flag + 3*ifac + iploc ;
        face2int_face_comm_index2[i] = i;
        ++i;
      }
      else if ( rankVois == rank ) {
        /* Add the elt kvois to communicators */
        pext_face_comm = &parmesh->ext_face_comm[shared_face[rankCur]];
        pext_face_comm->int_comm_index[idx[shared_face[rankCur]]++] = i;

        face2int_face_comm_index1[i] = 12*mesh->tetra[kvois].flag +
          3*ifacvois + iplocvois;
        face2int_face_comm_index2[i] = i;
        ++i;
      }
    }
  }
  PMMG_DEL_MEM(mesh,mesh->adja,int,"dealloc mesh adja");
  PMMG_DEL_MEM(parmesh,idx,int,"deallocating idx");

  PMMG_CALLOC(parmesh,parmesh->int_node_comm,1,PMMG_Int_comm,
              "allocating int_node_comm",return 0);
  PMMG_CALLOC(parmesh,parmesh->int_face_comm,1,PMMG_Int_comm,
              "allocating int_face_comm",return 0);

  /* We have 1 Grp per proc, thus : int_*_comm->nitem : nitem_int_*_comm */
  parmesh->int_node_comm->nitem = nitem_int_node_comm;
  parmesh->int_face_comm->nitem = nitem_int_face_comm;

  return 1;
}

/**
 * \param mesh pointer toward a MMG5 mesh structure
 * \param met pointer toward a MMG5 solution structure
 * \param rank rank of the MPI process
 * \param np number of points of the local mesh
 * \param nxp number of xpoints in the local mesh
 * \param nxt number of local xtetra in the local mesh
 * \param pointPerm array of new point positions
 * \param xPointPerm array of new xPoint positions
 * \param xTetraPerm array of new xTetra positions
 *
 * \return 1 if success, 0 if fail
 *
 * Create the mesh that will stay on the processor rank.
 *
 */
static inline
int PMMG_create_localMesh(MMG5_pMesh mesh,MMG5_pSol met,int rank,int np,int nxp,
                          int nxt,int *pointPerm,int *xPointPerm,int*xTetraPerm) {

  /** Compact tetrahedra on the proc */
  if ( !PMMG_packTetraOnProc(mesh,rank) ) return 0;

  /** Mesh permutations */
  if ( !PMMG_permuteMesh(mesh,met,pointPerm,xPointPerm,xTetraPerm,np,nxp,nxt) )
    return 0;

  if ( !PMMG_link_mesh( mesh ) ) return 0;

  /** Update xtetra edge tags */
  if ( PMMG_SUCCESS != PMMG_bdryUpdate( mesh ) ) return 0;

  /** Adjacency reconstruction */
  if ( !MMG3D_hashTetra( mesh, 0 ) ) return 0;

//  if ( parmesh->ddebug ) {
//    grplst_meshes_to_saveMesh( parmesh->listgrp, 1, parmesh->myrank, "End_distributeMesh_proc");
//    if ( met )
//      PMMG_saveSol( parmesh, filename );
//  }
  return 1;
}

/**
 * \param parmesh pointer toward a PMMG parmesh structure.
 * \param part pointer toward the metis array containing the partitions.
 *
 * \return 0 (on all procs) if fail, 1 otherwise
 *
 * Delete parts of the mesh not on the processor.
 */
int PMMG_partBcast_mesh( PMMG_pParMesh parmesh )
{
  PMMG_pGrp      grp = NULL;
  MMG5_pMesh     mesh = NULL;
  MMG5_pSol      met = NULL;
  idx_t          *part = NULL;
  int            nprocs,rank,ne,old_np,np,nxt,nxp;
  int            *shared_pt,*shared_face;
  int            *pointPerm = NULL, *xTetraPerm = NULL, *xPointPerm = NULL;
  int8_t         *seen_shared_pt = NULL;
  int            ier,ieresult;
  MPI_Datatype   metis_dt;

  ier = 1;

  /** Proc 0 send the mesh to the other procs */
  nprocs = parmesh->nprocs;
  rank   = parmesh->myrank;

  assert(parmesh->ngrp == 1);
  grp    = &parmesh->listgrp[0];
  mesh   = grp->mesh;
  met    = grp->met;

  /** Call metis for partionning */
  ne = mesh->ne;
  PMMG_CALLOC ( parmesh,part,ne,idx_t,"allocate metis buffer", ier=5 );

  if ( (!parmesh->myrank) && nprocs > 1 ) {
    if ( !PMMG_part_meshElts2metis( parmesh, part, parmesh->nprocs ) ) {
      ier = 5;
    }
    if( !PMMG_fix_contiguity_centralized( parmesh,part ) ) ier = 5;
  }

  /** Send the partition data to the other procs */
  if ( IDXTYPEWIDTH == 32 )
    metis_dt = MPI_INT32_T;
  else if ( IDXTYPEWIDTH == 64 )
    metis_dt = MPI_INT64_T;
  else {
    fprintf(stderr,"  ## Error: %s: unable to detect the metis integer width (%d).\n",
            __func__,IDXTYPEWIDTH);
    ier = 4;
  }

  MPI_CHECK( MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MAX, parmesh->comm ),
             ier = ieresult = 4 );
  if ( ieresult>1 ) goto unalloc;

  MPI_CHECK( MPI_Bcast( &part[0], ne, metis_dt, 0, parmesh->comm ),
             ier = 3 );

  /* Memory repartition */
  if ( !PMMG_parmesh_updateMemMax( parmesh,50,1 ) ) ier = 3;

  /** Mark the mesh to detect entities that will stay on the proc as well as
   * shared entites with the other procs */
  old_np = mesh->np;
  if ( !PMMG_mark_localMesh(parmesh,part,mesh,&np,&nxp,&nxt,&pointPerm,
                            &xPointPerm,&xTetraPerm,&shared_pt,&shared_face,
                            &seen_shared_pt) ) ier = 2;

  /** Communicators creation */
  if ( !PMMG_create_communicators(parmesh,part,shared_pt,shared_face,seen_shared_pt) )
    ier = 2;

  PMMG_DEL_MEM(parmesh,parmesh->int_node_comm->intvalues,int,"intvalues");

  /** Local mesh creation */
  if ( !PMMG_create_localMesh(mesh,met,rank,np,nxp,nxt,pointPerm,xPointPerm,xTetraPerm) )
    ier = 2;

  /** Check grps contiguity */
  ier = PMMG_checkAndReset_grps_contiguity( parmesh );

unalloc:
  if ( ier < 5 ) {
    PMMG_DEL_MEM(parmesh,part,idx_t,"deallocate metis buffer");

    MPI_CHECK( MPI_Allreduce( &ier,&ieresult,1,MPI_INT,MPI_MAX,parmesh->comm ),
               ieresult = 4);

    if ( ier < 3 ) {
      PMMG_DEL_MEM(parmesh,xPointPerm,int,"deallocate metis buffer");
      PMMG_DEL_MEM(parmesh,xTetraPerm,int,"deallocate xTetraPerm");
      PMMG_DEL_MEM(parmesh,pointPerm,int,"deallocate pointPerm");
      PMMG_DEL_MEM(parmesh,shared_face,int,"deallocate shared_face");
      PMMG_DEL_MEM(parmesh,seen_shared_pt,int8_t,"deallocate seen_shared_pt");
      PMMG_DEL_MEM(parmesh,shared_pt,int,"deallocate shared_pt");
    }
  }

#ifdef DEBUG
  if ( ieresult<2 && parmesh->nprocs > 1 ) {
    assert ( PMMG_check_extFaceComm ( parmesh ) );
    assert ( PMMG_check_intNodeComm ( parmesh ) );
    assert ( PMMG_check_extNodeComm ( parmesh ) );
  }
#endif

  return ieresult==1;
}

/**
 * \param parmesh pointer toward a PMMG parmesh structure.
 * \param part pointer toward the metis array containing the partitions.
 *
 * \return 0 (on all procs) if fail, 1 otherwise
 *
 * Delete parts of the mesh not on the processor.
 */
int PMMG_distribute_mesh( PMMG_pParMesh parmesh )
{
  PMMG_pGrp  grp;
  MMG5_pMesh mesh;
  idx_t      *part;
  int        igrp,ier,ieresult;
  size_t     available,oldMemMax;

  ier = 1;

  if( !PMMG_create_empty_communicators( parmesh ) ) return 0;
  if( parmesh->nprocs == 1 ) return 1;

  /** Proc 0 send the mesh to the other procs */
  if( parmesh->myrank == parmesh->info.root ) {

    grp    = &parmesh->listgrp[0];
    mesh   = grp->mesh;

    /** Call metis for partionning */
    PMMG_CALLOC ( parmesh,part,mesh->ne,idx_t,"allocate metis buffer", ier=5 );

    if ( !PMMG_part_meshElts2metis( parmesh, part, parmesh->nprocs ) ) {
      ier = 5;
    }
    if( !PMMG_fix_contiguity_centralized( parmesh,part ) ) ier = 5;

    /* Split grp 0 into (nprocs) groups */
    ier = PMMG_split_grps( parmesh,0,parmesh->nprocs,part,1 );

    /** Check grps contiguity */
    ier = PMMG_checkAndReset_grps_contiguity( parmesh );

    PMMG_DEL_MEM(parmesh,part,idx_t,"deallocate metis buffer");
  }

  {
    int igrp;
    for( igrp = 0; igrp < parmesh->ngrp; igrp++ ) {
      grp    = &parmesh->listgrp[igrp];
      int i,iel,ifac;
      MMG5_pTetra pt;
      MMG5_pxTetra pxt;
      for( i = 0; i < grp->nitem_int_face_comm; i++ ) {
        iel  =  grp->face2int_face_comm_index1[i] / 12;
        ifac = (grp->face2int_face_comm_index1[i] % 12) / 3;
        pt = &grp->mesh->tetra[iel];
        assert( pt->xt );
        pxt = &grp->mesh->xtetra[pt->xt];
        assert( pxt->ftag[ifac] & MG_PARBDY );
      }
    }
  }

  /** Distribute the groups over the processors */
  if( parmesh->myrank != parmesh->info.root ) parmesh->ngrp = 0;
 
  PMMG_CALLOC ( parmesh,part,parmesh->nprocs+1,idx_t,"allocate metis buffer", ier=5 );
  for( igrp = 0; igrp <= parmesh->nprocs; igrp++ ) part[igrp] = igrp;
  ier = PMMG_transfer_all_grps(parmesh,part);
  if ( ier <= 0 ) {
    fprintf(stderr,"\n  ## Group distribution problem.\n");
  }

  assert( parmesh->ngrp = 1);
  grp = &parmesh->listgrp[0];
  mesh = grp->mesh;
 
  PMMG_TRANSFER_AVMEM_TO_PARMESH(parmesh,available,oldMemMax);
  PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh,available,oldMemMax);
  if ( (!mesh->adja) && !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Error: %s: tetra hashing problem. Exit program.\n",
            __func__);
    return 0;
  }
  PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,mesh,available,oldMemMax);

  {
    grp = &parmesh->listgrp[0];
    int i,iel,ifac;
    MMG5_pTetra pt;
    MMG5_pxTetra pxt;
    for( i = 0; i < grp->nitem_int_face_comm; i++ ) {
      iel  =  grp->face2int_face_comm_index1[i] / 12;
      ifac = (grp->face2int_face_comm_index1[i] % 12) / 3;
      pt = &grp->mesh->tetra[iel];
      assert( pt->xt );
      pxt = &grp->mesh->xtetra[pt->xt];
      assert( pxt->ftag[ifac] & MG_PARBDY );
    }
  }

//  char basename[48];
//  sprintf(basename,"my_transfer_");
//  PMMG_listgrp_to_saveMesh( parmesh,basename ); 

  /* Check the communicators */
  assert ( PMMG_check_intNodeComm(parmesh) && "Wrong internal node comm" );
  assert ( PMMG_check_intFaceComm(parmesh) && "Wrong internal face comm" );
  assert ( PMMG_check_extNodeComm(parmesh) && "Wrong external node comm" );
  assert ( PMMG_check_extFaceComm(parmesh) && "Wrong external face comm" );

  /* The part array is deallocated when groups to be sent are merged */

ieresult = 1;
  return ieresult==1;
}

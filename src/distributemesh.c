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
#include "mpitypes.h"

/**
 * \param xtetra pointer toward a table containing the xtetra structures.
 * \param *perm pointer toward the permutation table (to perform in place
 * permutations).
 * \param ind1 index of the first xtetra to swap.
 * \param ind2 index of the second xtetra to swap.
 *
 * Swap two xtetra in the table of xtetrahedras.
 *
 */
static inline
void PMMG_swapxTetra(MMG5_pxTetra xtetra, int* perm, int ind1, int ind2) {
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
 *
 */
static inline
void PMMG_swapxPoint(MMG5_pxPoint xpoint, int* perm, int ind1, int ind2) {
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
 * \param sol pointer toward a table containing the solution structures.
 * \param *perm pointer toward the permutation table (to perform in place
 * permutations).
 * \param ind1 index of the first xpoint to swap.
 * \param ind2 index of the second xpoint to swap.
 * \param solsize size of the solution.
 *
 * Swap two points in the table of points.
 *
 */
static inline
void PMMG_swapPoint(MMG5_pPoint point,double* sol,int* perm,
                    int ind1,int ind2,int solsiz) {
  MMG5_Point ppttmp;
  MMG5_Sol   soltmp;
  int        tmp,addr2,addr1;

  /** 1- swap the xpoint */
  memcpy(&ppttmp      ,&point[ind2], sizeof(MMG5_Point));
  memcpy(&point[ind2] ,&point[ind1], sizeof(MMG5_Point));
  memcpy(&point[ind1] ,&ppttmp      ,sizeof(MMG5_Point));

  /** 2- swap the sols */
  if ( sol ) {
    addr1 = ind1*solsiz;
    addr2 = ind2*solsiz;
    memcpy(&soltmp    ,&sol[addr2],solsiz*sizeof(double));
    memcpy(&sol[addr2],&sol[addr1],solsiz*sizeof(double));
    memcpy(&sol[addr1],&soltmp    ,solsiz*sizeof(double));
  }

  /** 3- swap the permutation table */
  tmp        = perm[ind2];
  perm[ind2] = perm[ind1];
  perm[ind1] = tmp;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Send the initial mesh from proc 0 toward the other procs.
 *
 */
int PMMG_bcastMesh(PMMG_pParMesh parmesh) {
  PMMG_pGrp       grp;
  MMG5_pMesh      mesh;
  MMG5_pSol       sol;
  MPI_Datatype    mpi_light_point,mpi_light_tetra,mpi_tria,mpi_edge;
  int             rank;

  /** Proc 0 send the mesh to the other procs */
  grp    = &parmesh->listgrp[0];
  mesh   = grp->mesh;
  sol    = grp->sol;
  rank   = parmesh->myrank;

  /* Mesh */
  MPI_Bcast( &mesh->np,     1, MPI_INT,       0, parmesh->comm);
  MPI_Bcast( &mesh->ne,     1, MPI_INT,       0, parmesh->comm);
  MPI_Bcast( &mesh->nt,     1, MPI_INT,       0, parmesh->comm);
  MPI_Bcast( &mesh->na,     1, MPI_INT,       0, parmesh->comm);
  MPI_Bcast( &mesh->ntmax,  1, MPI_INT,       0, parmesh->comm);
  MPI_Bcast( &mesh->memMax, 1, MPI_LONG_LONG, 0, parmesh->comm);

  mesh->nemax = mesh->nei = mesh->ne;
  mesh->nenil = 0;
  mesh->npmax = mesh->npi = mesh->np;
  mesh->npnil = 0;
  mesh->nti   = mesh->nt;
  mesh->xtmax = mesh->ntmax;

  /* Solution */
  MPI_Bcast( &sol->size,    1, MPI_INT,       0, parmesh->comm);
  MPI_Bcast( &sol->npmax,   1, MPI_INT,       0, parmesh->comm);
  MPI_Bcast( &sol->np,      1, MPI_INT,       0, parmesh->comm);

  sol->npi = sol->np;
  sol->ver = mesh->ver;
  sol->dim = mesh->dim;

  if ( rank ) {
    _MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(MMG5_Point),"initial vertices",
                  fprintf(stderr,"  Exit program.\n");
                  return(0));
    _MMG5_SAFE_CALLOC(mesh->point,mesh->npmax+1,MMG5_Point);

    _MMG5_ADD_MEM(mesh,(mesh->nemax+1)*sizeof(MMG5_Tetra),"initial tetrahedra",
                fprintf(stderr,"  Exit program.\n");
                return(0));
    _MMG5_SAFE_CALLOC(mesh->tetra,mesh->nemax+1,MMG5_Tetra);

    if ( mesh->nt ) {
      _MMG5_ADD_MEM(mesh,(mesh->nt+1)*sizeof(MMG5_Tria),"initial triangles",
                    return(0));
      _MMG5_SAFE_CALLOC(mesh->tria,mesh->nt+1,MMG5_Tria);
    }

    if ( mesh->na ) {
      _MMG5_ADD_MEM(mesh,(mesh->na+1)*sizeof(MMG5_Edge),"initial edges",
                    return(0));
      _MMG5_SAFE_CALLOC(mesh->edge,mesh->na+1,MMG5_Edge);
    }

    if ( sol->npmax ) {
      _MMG5_ADD_MEM(mesh,(sol->size*(sol->npmax+1))*sizeof(double),"initial solution",
                    fprintf(stderr,"  Exit program.\n");
                    return(0));
      _MMG5_SAFE_CALLOC(sol->m,(sol->size*(sol->npmax+1)),double);
    }
  }

  if ( !PMMG_create_MPI_lightPoint (mesh->point,  &mpi_light_point ) ) return(0);
  if ( !PMMG_create_MPI_lightTetra (mesh->tetra,  &mpi_light_tetra ) ) return(0);
  if ( mesh->nt && !PMMG_create_MPI_Tria(mesh->tria,   &mpi_tria   ) ) return(0);
  if ( mesh->na && !PMMG_create_MPI_Edge(mesh->edge,   &mpi_edge   ) ) return(0);

  MPI_Bcast( mesh->point,  mesh->np+1,    mpi_light_point,  0, parmesh->comm);
  MPI_Bcast( mesh->tetra,  mesh->ne+1,    mpi_light_tetra,  0, parmesh->comm);
  if ( mesh->nt ) MPI_Bcast( mesh->tria, mesh->nt+1, mpi_tria, 0, parmesh->comm);
  if ( mesh->na ) MPI_Bcast( mesh->edge, mesh->na+1, mpi_edge, 0, parmesh->comm);
  if ( sol->m )
    MPI_Bcast( sol->m,sol->size*(sol->npmax+1),MPI_DOUBLE, 0, parmesh->comm);

  MPI_Type_free(&mpi_light_point);
  MPI_Type_free(&mpi_light_tetra);
  if ( mesh->nt ) MPI_Type_free(&mpi_tria);
  if ( mesh->na ) MPI_Type_free(&mpi_edge);

  return 1;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param part pointer toward an array of int containing the partitions.
 * \return 0 if fail, 1 otherwise.
 *
 * Delete parts of the mesh not on the processor.
 *
 */
int PMMG_distributeMesh(PMMG_pParMesh parmesh) {
  PMMG_pGrp       grp;
  MMG5_pMesh      mesh;
  MMG5_pSol       sol;
  MMG5_pTetra     pt,ptnew;
  MMG5_pxTetra    pxt;
  MMG5_pPoint     ppt;
  PMMG_pext_comm  pext_comm;
  idx_t           *part;
  int             nprocs,rank,np,ne,nxt,nxp;
  int             *pointPerm,*xPointPerm,*xTetraPerm;
  int             ip,iploc,ifac,i,j,k,*idx,kvois,rankVois;
  int             *node2int_node_comm_index1,*node2int_node_comm_index2;
  int             nitem_int_node_comm,next_node_comm,*seenRanks;
  int             inIntComm,nbl;
  int8_t          *pointRanks;


  /** Proc 0 send the mesh to the other procs */
  nprocs = parmesh->nprocs;
  grp    = &parmesh->listgrp[0];
  mesh   = grp->mesh;
  sol    = grp->sol;
  rank   = parmesh->myrank;

  /** Call metis for partionning*/
  _MMG5_SAFE_CALLOC(part,(parmesh->listgrp[0].mesh)->ne,idx_t);

  if ( nprocs > 1 && !PMMG_metispartitioning(parmesh,part) ) return 0;

  /** Remove the part of the mesh that are not on the proc rank */
  _MMG5_SAFE_CALLOC(seenRanks,nprocs,int);
  _MMG5_SAFE_CALLOC(pointRanks,nprocs*mesh->np,int8_t);

  _MMG5_SAFE_CALLOC(pointPerm,mesh->np+1,int);
  _MMG5_SAFE_CALLOC(xTetraPerm,mesh->xtmax+1,int);
  _MMG5_SAFE_CALLOC(xPointPerm,mesh->xp+1,int);

  nxp = 0;
  nxt = 0;
  np  = 0;
  ne  = 0;
  nitem_int_node_comm = 0;

  /* Reset the tmp field of points */
  for ( k=1; k<=mesh->np; k++ )
    mesh->point[k].tmp = 0;

  /** Mark mesh entities that will stay on the proc and count the number of
   * point that must be communicate with the other procs  */
  for ( k=1; k<=mesh->ne; k++ ) {
    pt = &mesh->tetra[k];
    pt->mark = part[k-1];

    if ( pt->mark != rank )
      continue;

    for ( ifac=0; ifac<4; ifac++ ) {
      kvois    = mesh->adja[4*k-3+ifac]/4;

      if ( kvois )
        rankVois = part[kvois-1];
      else
        rankVois = rank;

      /* Mark the interfaces between two procs */
      if ( rank != rankVois ) {
        if ( !pt->xt ) {
          ++mesh->xt;
          if ( mesh->xt > mesh->xtmax ) {
            /* realloc of xtetras table */
            _MMG5_SAFE_RECALLOC(xTetraPerm,mesh->xtmax,(int)(0.2*mesh->xtmax)+1,
                                int,"larger tetra permutation table");

            _MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                               "larger xtetra table",
                               mesh->xt--;
                               fprintf(stderr,"  Exit program.\n");
                               return(0));

          }
          pt->xt = mesh->xt;
        }
        pxt = &mesh->xtetra[pt->xt];
        /* Parallel face */
        pxt->ftag[ifac] |= (MG_PARBDY + MG_BDY + MG_REQ);

        /* Parallel edges */
#warning using the MG_REQ tag, we will loose the "true" required tags
        for ( j=0; j<3; ++j ) {
          pxt->tag[_MMG5_iarf[ifac][j]] |= (MG_PARBDY + MG_BDY + MG_REQ);
        }

      }

      for ( j=0; j<3; ++j ) {
        iploc = _MMG5_idir[ifac][j];
        ip    = pt->v[iploc];
        ppt   = &mesh->point[ip];

        /* Count (and mark) each point that will be shared between me and the
         * proc rankVois */
        if ( rankVois != rank && !pointRanks[nprocs*(ip-1)+rankVois] ) {
          pointRanks[nprocs*(ip-1)+rankVois] = 1;
          ++seenRanks[rankVois];

          /* Mark parallel vertex */
          ppt->tag |= (MG_PARBDY + MG_BDY + MG_REQ);
        }

        if ( !ppt->tmp ) {
          /* Mark the new point index and update the table of permutation */
          ppt->tmp = ++np;
          pointPerm[ip] = np;

          /* update the table of permutation for xPoint if needed */
          if ( ppt->xp ) {
            xPointPerm[ppt->xp] = ++nxp;
            ppt->xp             = nxp;
          }
        }
      }
    }

    /* Update the table of permutation of the Tetra and the tetra vertices indices */
    for ( j=0; j<4; ++j )
      pt->v[j] = mesh->point[pt->v[j]].tmp;

    /* update the table of permutation for xTetra if needed */
    if ( !pt->xt ) continue;

    xTetraPerm[pt->xt] = ++nxt;
    pt->xt = nxt;
  }
  _MMG5_SAFE_FREE(mesh->adja);

  /** Count the number of external node communicators and initialize it */
  next_node_comm = 0;
  for ( k=0; k<nprocs; ++k ) {
    if ( seenRanks[k] ) ++next_node_comm;
  }

  parmesh->next_node_comm = next_node_comm;
  _MMG5_SAFE_CALLOC(parmesh->ext_node_comm,next_node_comm,PMMG_ext_comm);

  next_node_comm = 0;
  for ( k=0; k<nprocs; ++k ) {
    if ( seenRanks[k] ) {
      pext_comm = &parmesh->ext_node_comm[next_node_comm];
      pext_comm->color_in  = rank;
      pext_comm->color_out = k;
      pext_comm->nitem     = seenRanks[k];
      _MMG5_SAFE_CALLOC(pext_comm->int_comm_index,pext_comm->nitem,int);
      /* Use seenRanks to store the idx of the external communicator me->k */
      seenRanks[k] = next_node_comm++;
    }
  }

  /** Initialize the internal node communicator */
  nitem_int_node_comm = 0;
  for ( k=1; k<=mesh->np; k++ ) {

    if ( !mesh->point[k].tmp )  continue;

    for ( j=0; j<nprocs; ++j ) {
      if ( pointRanks[nprocs*(k-1)+j] ) {
        ++nitem_int_node_comm;
        break;
      }
    }
  }

  grp->nitem_int_node_comm = nitem_int_node_comm;
  _MMG5_SAFE_CALLOC(grp->node2int_node_comm_index1,nitem_int_node_comm,int);
  _MMG5_SAFE_CALLOC(grp->node2int_node_comm_index2,nitem_int_node_comm,int);

  /** Travel through the mesh and fill the communicators */
  i = 0;
  node2int_node_comm_index1 = grp->node2int_node_comm_index1;
  node2int_node_comm_index2 = grp->node2int_node_comm_index2;

  /* Idx is used to store the external communicator cursor */
  _MMG5_SAFE_CALLOC(idx,parmesh->next_node_comm,int);

  for ( k=1; k<=mesh->np; k++ ) {
    if ( !mesh->point[k].tmp )  continue;

    inIntComm = 0;
    for ( j=0; j<nprocs; ++j ) {
      pext_comm = &parmesh->ext_node_comm[seenRanks[j]];

      if ( pointRanks[nprocs*(k-1)+j] == 1 ) {
        /* Add point in external communicator */
        pext_comm->int_comm_index[idx[seenRanks[j]]++] = i;

        if ( !inIntComm ) {
          /* Add point in internal communicator */
          inIntComm = 1;
          node2int_node_comm_index1[i] = mesh->point[k].tmp;
          node2int_node_comm_index2[i] = i;
        }
      }
    }
    /* Increment internal comm cursor */
    if ( inIntComm )  ++i;
  }
  _MMG5_SAFE_CALLOC(parmesh->int_node_comm,1,PMMG_int_comm);
  /* We have 1 Grp per proc, thus : int_node_comm.nitem : nitem_int_node_comm */
  parmesh->int_node_comm->nitem = nitem_int_node_comm;

  /** Compact tetrahedra on the proc */
  ne  = 0;
  nbl = 1;
  for ( k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];

    if ( pt->mark!=rank ) continue;
    ++ne;

    if ( k!=nbl ) {
      ptnew = &mesh->tetra[nbl];
      memcpy(ptnew,pt,sizeof(MMG5_Tetra));
    }
    ++nbl;
  }
  mesh->ne = ne;

  /** Compact xtetra on the proc */
  for ( k=1; k<=mesh->xt; ++k ) {
    while ( xTetraPerm[k] != k && xTetraPerm[k] )
      PMMG_swapxTetra(mesh->xtetra,xTetraPerm,k,xTetraPerm[k]);
  }
  mesh->xt = nxt;

  /** Compact vertices on the proc: in place permutations */
  for ( k=1; k<=mesh->np; ++k ) {
    while ( pointPerm[k] != k && pointPerm[k] )
      PMMG_swapPoint(mesh->point,sol->m,pointPerm,k,pointPerm[k],sol->size);
  }
  mesh->np = np;

  /** Compact xpoint on the proc: in place permutations */
  for ( k=1; k<=mesh->xp; ++k ) {
    while ( xPointPerm[k] != k && xPointPerm[k] )
      PMMG_swapxPoint(mesh->xpoint,xPointPerm,k,xPointPerm[k]);
  }
  mesh->xp = nxp;

  _MMG5_SAFE_FREE(part);

  _MMG5_SAFE_FREE(seenRanks);
  _MMG5_SAFE_FREE(pointRanks);
  _MMG5_SAFE_FREE(idx);

  _MMG5_SAFE_FREE(pointPerm);
  _MMG5_SAFE_FREE(xPointPerm);
  _MMG5_SAFE_FREE(xTetraPerm);

  /** Update xtetra edge tags */
  if ( !PMMG_bdryUpdate(mesh) ) return 0;

  /** Adjacency reconstruction */
  if ( !MMG3D_hashTetra(parmesh->listgrp[0].mesh,0) ) return(0);

  if ( parmesh->ddebug ) {
    /* sprintf(filename,"End_distributeMesh_proc%d.mesh",rank); */
    /* _MMG3D_bdryBuild(parmesh->listgrp[0].mesh); */
    /* PMMG_saveMesh(parmesh,filename); */
    /* if ( sol ) PMMG_saveSol(parmesh,filename); */
  }

  return(1);
}

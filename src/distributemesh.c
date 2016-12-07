/**
 * \file distributemesh.c
 * \brief Distribute the mesh on the processors.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */


#include "libparmmg.h"

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
 * \param *perm pointer toward the permutation table (to perform in place
 * permutations).
 * \param ind1 index of the first xpoint to swap.
 * \param ind2 index of the second xpoint to swap.
 *
 * Swap two points in the table of points.
 *
 */
static inline
void PMMG_swapPoint(MMG5_pPoint point, int* perm, int ind1, int ind2) {
  MMG5_Point ppttmp;
  int         tmp;

  /** 1- swap the xpoint */
  memcpy(&ppttmp      ,&point[ind2], sizeof(MMG5_Point));
  memcpy(&point[ind2] ,&point[ind1], sizeof(MMG5_Point));
  memcpy(&point[ind1] ,&ppttmp      ,sizeof(MMG5_Point));

  /** 2- swap the permutation table */
  tmp        = perm[ind2];
  perm[ind2] = perm[ind1];
  perm[ind1] = tmp;
}


/**
 * \param parmesh pointer toward the mesh structure.
 * \param part pointer toward an array of int containing the partitions.
 * \return 0 if fail, 1 otherwise.
 *
 * Delete parts of the mesh not on the processor.
 *
 */
int PMMG_distributeMesh(PMMG_pParMesh parmesh,int *part) {
  PMMG_pGrp       grp;
  MMG5_pMesh      mesh;
#warning add solution/metric communication
  MMG5_pSol       sol;
  MMG5_pTetra     pt,ptnew;
  MMG5_pxTetra    pxt;
  MMG5_pPoint     ppt, pptnew;
  MMG5_pxPoint    pxp;
  PMMG_pext_comm  pext_comm;
  int             nprocs,rank,np,ne,nbl,nxt,nxp;
  int             *pointPerm,*xPointPerm,*xTetraPerm;
  int             ip,iploc,ifac,i,j,k,*idx,kvois,rankVois;
  int             *node2int_node_comm_index1,*node2int_node_comm_index2;
  int             nitem_int_node_comm,next_node_comm,*seenRanks;
  int             inIntComm;
  char            filename[11];
  int8_t          *pointRanks;

  nprocs = parmesh->nprocs;
  grp    = parmesh->listgrp;
  mesh   = grp[0].mesh;
  rank   = parmesh->myrank;

#warning to trash
  printf( " je suis le proc %d\n", parmesh->myrank);

  _MMG5_SAFE_CALLOC(seenRanks,nprocs,int);
  _MMG5_SAFE_CALLOC(pointRanks,nprocs*mesh->np,int8_t);

  _MMG5_SAFE_CALLOC(pointPerm,mesh->np+1,int);
  _MMG5_SAFE_CALLOC(xTetraPerm,mesh->xt+1,int);
  _MMG5_SAFE_CALLOC(xPointPerm,mesh->xp+1,int);

  nxp = 0;
  nxt = 0;
  np  = 0;
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

      for ( j=0; j<3; ++j ) {
        iploc = _MMG5_idir[ifac][j];
        ip    = pt->v[iploc];
        ppt   = &mesh->point[ip];

        /* Count (and mark) each point that will be shared between me and the
         * proc rankVois */

        if ( rankVois != rank && !pointRanks[nprocs*(ip-1)+rankVois] ) {
          pointRanks[nprocs*(ip-1)+rankVois] = 1;
          ++seenRanks[rankVois];
        }

        if ( !ppt->tmp ) {
          /* Mark the new point index and update the table of permutation */
          ppt->tmp = ++np;
          pointPerm[ip] = np;

          /* update the table of permutation for xPoint if needed */
          if ( ppt->xp ) {
            pxp                 = &mesh->xpoint[ppt->xp];
            xPointPerm[ppt->xp] = ++nxp;
            ppt->xp             = nxp;
          }
        }
      }
    }
    for ( j=0; j<4; ++j )
      pt->v[j] = mesh->point[pt->v[j]].tmp;

    /* update the table of permutation for xTetra if needed */
    if ( !pt->xt ) continue;

    pxt    = &mesh->xtetra[pt->xt];
    xTetraPerm[pt->xt] = ++nxt;
    pt->xt = nxt;
  }

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

  /** Compact xtetra on the proc: in place permutations */
  for ( k=1; k<=mesh->xt; ++k ) {
    while ( xTetraPerm[k] != k && xTetraPerm[k] )
      PMMG_swapxTetra(mesh->xtetra,xTetraPerm,k,xTetraPerm[k]);
  }
  mesh->xt = nxt;

  /** Compact vertices on the proc: in place permutations */
  for ( k=1; k<=mesh->np; ++k ) {
    while ( pointPerm[k] != k && pointPerm[k] )
      PMMG_swapPoint(mesh->point,pointPerm,k,pointPerm[k]);
  }
  mesh->np = np;

  /** Compact xpoint on the proc: in place permutations */
  for ( k=1; k<=mesh->xp; ++k ) {
    while ( xPointPerm[k] != k && xPointPerm[k] )
      PMMG_swapxPoint(mesh->xpoint,xPointPerm,k,xPointPerm[k]);
  }
  mesh->xp = nxp;

#warning Try to remove the adjacency reconstruction and packing (we need to adapt Mmg to allow to provide xTetra/points instead of Triangles)

  /** Tetra adjacency reconstruction */
  _MMG5_SAFE_FREE(parmesh->listgrp[0].mesh->adja);
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"  ## PMMG Hashing problem (1). Exit program.\n");
    return(0);
  }

  /* Pack the mesh */
  if ( !_MMG3D_packMesh(mesh,NULL,NULL ) ) {
    fprintf(stderr,"  ## PMMG Packing problem (1). Exit program.\n");
    return(0);
  }

  sprintf(filename,"proc%d.mesh",rank);
  MMG3D_saveMesh(mesh,filename);

  _MMG5_SAFE_FREE(seenRanks);
  _MMG5_SAFE_FREE(pointRanks);
  _MMG5_SAFE_FREE(idx);

  _MMG5_SAFE_FREE(pointPerm);
  _MMG5_SAFE_FREE(xPointPerm);
  _MMG5_SAFE_FREE(xTetraPerm);

  return(1);
}

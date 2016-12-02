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
  PMMG_pGrp    grp;
  MMG5_pMesh   mesh;
#warning add solution/metric communication
  MMG5_pSol    sol;
  MMG5_pTetra  pt,ptnew;
  MMG5_pxTetra pxt;
  MMG5_pPoint  ppt, pptnew;
  MMG5_pxPoint pxp;
  int          rank,np,ne,nbl,nxt,nxp;
  int          *pointPerm,*xPointPerm,*xTetraPerm;
  int          ip,iploc,ifac,j,k,kvois,rankVois;
  char         filename[11];

  grp  = parmesh->listgrp;
  mesh = grp[0].mesh;
  rank = parmesh->myrank;

#warning to trash
  printf( " je suis le proc %d\n", parmesh->myrank);

  _MMG5_SAFE_CALLOC(pointPerm,mesh->np,int);
  _MMG5_SAFE_CALLOC(xTetraPerm,mesh->xt,int);
  _MMG5_SAFE_CALLOC(xPointPerm,mesh->xp,int);

  nxp = 0;
  nxt = 0;
  np  = 0;


  /** Reset the tmp field of points */
  for ( k=1; k<=mesh->np; k++ )
    mesh->point[k].tmp = 0;

  /** Count and mark mesh entities over the proc */
  for ( k=1; k<=mesh->ne; k++ ) {
    pt = &mesh->tetra[k];
    pt->mark = part[k-1];

    if ( pt->mark != rank )
      continue;

    for ( ifac=0; ifac<4; ifac++ ) {
      kvois    = mesh->adja[4*k-3+ifac]/4;
      rankVois = part[kvois-1];

      for ( j=0; j<3; ++j ) {
        iploc = _MMG5_idir[ifac][j];
        ip    = pt->v[iploc];
        ppt   = &mesh->point[ip];

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
        pt->v[iploc] = ppt->tmp;
      }
    }
    /* update the table of permutation for xTetra if needed */
    if ( !pt->xt ) continue;

    pxt    = &mesh->xtetra[pt->xt];
    xTetraPerm[pt->xt] = ++nxt;
    pt->xt = nxt;
  }

  /** Compact tetrahedra */
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

  /** Compact xtetra: in place permutations */
  for ( k=1; k<=mesh->xt; ++k ) {
    while ( xTetraPerm[k] != k && xTetraPerm[k] )
      PMMG_swapxTetra(mesh->xtetra,xTetraPerm,k,xTetraPerm[k]);
  }
  mesh->xt = nxt;

  /** Compact vertices: in place permutations */
  for ( k=1; k<=mesh->np; ++k ) {
    while ( pointPerm[k] != k && pointPerm[k] )
      PMMG_swapPoint(mesh->point,pointPerm,k,pointPerm[k]);
  }
  mesh->np = np;

  /** Compact xpoint: in place permutations */
  for ( k=1; k<=mesh->xp; ++k ) {
    while ( xPointPerm[k] != k && xPointPerm[k] )
      PMMG_swapxPoint(mesh->xpoint,xPointPerm,k,xPointPerm[k]);
  }
  mesh->xp = nxp;

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

  _MMG5_SAFE_FREE(pointPerm);
  _MMG5_SAFE_FREE(xPointPerm);
  _MMG5_SAFE_FREE(xTetraPerm);

  return(1);
}

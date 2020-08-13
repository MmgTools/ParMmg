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
 * \file hash_pmmg.c
 * \brief Hash function.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \author Nikolas Pattakos (Inria)
 * \author Luca Cirrottola (Inria)
 * \version
 * \copyright
 *
 */
#include "mmg3d.h"
#include "parmmg.h"

/**
 * \param mesh pointer toward a MMG5 mesh structure.
 * \param pHash pointer to the edge hash table.
 * \return PMMG_FAILURE
 *         PMMG_SUCCESS
 *
 * Hash the edges. Use the assumption that all paralle edges are seen by a
 * MG_PARBDY face on an xtetra.
 *
 */
int PMMG_hashPar( MMG5_pMesh mesh,MMG5_HGeom *pHash ) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  int          k,na;
  int          ifac,j,ia,i1,i2;

  /** Allocation of hash table to store parallel edges */
  na = (int)(mesh->np*0.2); // Euler-Poincare

  if ( 1 != MMG5_hNew( mesh, pHash, na, 3 * na ) ) return PMMG_FAILURE;

  /** Store parallel edges */
  for (k=1; k<=mesh->ne; ++k) {
    pt = &mesh->tetra[k];
    if ( !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];
    for( ifac = 0; ifac < 4; ifac++ ) {
      if ( !(pxt->ftag[ifac] & MG_PARBDY) ) continue;
      for ( j=0; j<3; j++ ) {
        ia = MMG5_iarf[ifac][j];
        i1 = MMG5_iare[ia][0];
        i2 = MMG5_iare[ia][1];
        MMG5_hEdge( mesh,pHash,pt->v[i1],pt->v[i2],pxt->edg[ia],pxt->tag[ia] );
      }
    }
  }

  return PMMG_SUCCESS;
}

/**
 * \param mesh pointer toward a MMG5 mesh structure.
 * \return PMMG_FAILURE
 *         PMMG_SUCCESS
 *
 * Hash the MG_PARBDY edges and store their tags (and refs). Then loop over the
 * tetrahedra to update the edge tags of the xtetra.
 *
 */
int PMMG_bdryUpdate( MMG5_pMesh mesh )
{
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_HGeom   hash;
  int          k,edg;
  int16_t      tag;
  int8_t       i,i1,i2;


  assert ( !mesh->htab.geom );

  /* Hash the MG_PARBDY edges */
  if( PMMG_hashPar(mesh,&hash) != PMMG_SUCCESS ) return PMMG_FAILURE;

  /** Update xtetra edge tag if needed */
  for (k=1; k<=mesh->ne; ++k) {

    pt = &mesh->tetra[k];
    if ( !pt->xt )
      continue;

    pxt = &mesh->xtetra[pt->xt];
    for ( i=0; i<6; ++i ) {
      if ( pxt->tag[i] & MG_PARBDY )
        continue;

      i1 = MMG5_iare[i][0];
      i2 = MMG5_iare[i][1];
      if ( !MMG5_hGet( &hash, pt->v[i1], pt->v[i2], &edg, &tag ) )
        continue;

      pxt->tag[i] |= tag;
      pxt->edg[i]  = edg;
    }
  }

  MMG5_DEL_MEM(mesh,hash.geom);

  if ( mesh->edge )
    PMMG_DEL_MEM(mesh,mesh->edge,MMG5_Edge,"deallocating edges");

  return PMMG_SUCCESS;
}

/**
 * \file hash_pmmg.c
 * \brief Hash function.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \author Nikolas Pattakos (Inria)
 * \version
 * \copyright
 *
 */
#include "mmg3d.h"
#include "parmmg.h"

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
  int          k,edg,na;
  int16_t      tag;
  int8_t       i,i1,i2;

  /** Allocation of hash table to store parallel edges */
  assert ( !mesh->htab.geom );

  na = (int)(mesh->np*0.2); // Euler-Poincare

  if ( 1 != _MMG5_hNew( mesh, &hash, na, 3 * na ) )
    return PMMG_FAILURE;

  /** Store parallel edges */
  for (k=1; k<=mesh->ne; ++k) {

    pt = &mesh->tetra[k];
    if ( !pt->xt )
      continue;

    pxt = &mesh->xtetra[pt->xt];
    for ( i=0; i<6; ++i ) {
      if ( !(pxt->tag[i] & MG_PARBDY) )
        continue;

      i1 = _MMG5_iare[i][0];
      i2 = _MMG5_iare[i][1];
      _MMG5_hEdge( mesh, &hash, pt->v[i1], pt->v[i2], pxt->edg[i], pxt->tag[i] );
    }
  }

  /** Update xtetra edge tag if needed */
  for (k=1; k<=mesh->ne; ++k) {

    pt = &mesh->tetra[k];
    if ( !pt->xt )
      continue;

    pxt = &mesh->xtetra[pt->xt];
    for ( i=0; i<6; ++i ) {
      if ( pxt->tag[i] & MG_PARBDY )
        continue;

      i1 = _MMG5_iare[i][0];
      i2 = _MMG5_iare[i][1];
      if ( !_MMG5_hGet( &hash, pt->v[i1], pt->v[i2], &edg, &tag ) )
        continue;

      pxt->tag[i] |= tag;
      pxt->edg[i]  = edg;
    }
  }

  PMMG_DEL_MEM(mesh,hash.geom,(hash.max+1),MMG5_hgeom,"Edge hash table");

  if ( mesh->edge )
    PMMG_DEL_MEM(mesh,mesh->edge,mesh->na+1,MMG5_Edge,"deallocating edges");

  return PMMG_SUCCESS;
}

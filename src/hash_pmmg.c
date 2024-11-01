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
 * \author Laetitia Mottet (UBordeaux)
 * \version
 * \copyright
 *
 */
#include "parmmg.h"

int PMMG_hashOldPar_pmmg( PMMG_pParMesh parmesh,MMG5_pMesh mesh,MMG5_Hash *hash ) {
  MMG5_pTetra    pt;
  MMG5_pxTetra   pxt;
  MMG5_pPoint    ppt;
  MMG5_hedge    *ph;
  int            k,kk,hmax,ia,ib;
  int8_t         i,i1,i2,j,l;
  unsigned int   key;

  /* adjust hash table params */
  hmax =(int)(3.71*6*mesh->xt);
  hash->siz  = 6*mesh->xt;
  hash->max  = hmax + 1;
  hash->nxt  = hash->siz;
  MMG5_ADD_MEM(mesh,(hash->max+1)*sizeof(MMG5_hedge),"hash table",return 0);
  MMG5_SAFE_CALLOC(hash->item,hash->max+1,MMG5_hedge,return 0);

  for (k=hash->siz; k<hash->max; k++)
    hash->item[k].nxt = k+1;

  /* loop on tetrahedra */
  for( k = 1; k <= mesh->ne; k++ ) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    if( !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];

    /* loop on vertices */
    for( i = 0; i < 4; i++ ) {
      ppt = &mesh->point[pt->v[i]];
      if( !(ppt->tag & MG_BDY) || !(ppt->tag & MG_OLDPARBDY) ) continue;

      /* loop on edges touching the vertex */
      for( j = 0; j < 3; j++ ) {
        l = MMG5_arpt[i][j];
        if( !(pxt->tag[l] & MG_NOM ) ) continue;

        i1 = MMG5_iare[l][0];
        i2 = MMG5_iare[l][1];

        /* compute key */
        ia  = MG_MIN(pt->v[i1],pt->v[i2]);
        ib  = MG_MAX(pt->v[i1],pt->v[i2]);
        key = (MMG5_KA*ia + MMG5_KB*ib) % hash->siz;
        ph  = &hash->item[key];

        /* store edge */
        if ( ph->a == 0 ) {
          ph->a = ia;
          ph->b = ib;
          ph->k = 6*k + l;
          ph->nxt = 0;
          ++ph->s;
          continue;
        }
        while ( ph->a ) {
          if ( ph->a == ia && ph->b == ib ) {
            ++ph->s;
            break;
          }
          else if ( !ph->nxt ) {
            ph->nxt = hash->nxt;
            ph = &hash->item[ph->nxt];
            assert(ph);

            if ( hash->nxt >= hash->max-1 ) {
              if ( mesh->info.ddebug ) {
                fprintf(stderr,"\n  ## Warning: %s: memory alloc problem (edge):"
                        " %d\n",__func__,hash->max);
              }
              MMG5_TAB_RECALLOC(mesh,hash->item,hash->max,MMG5_GAP,MMG5_hedge,
                                 "MMG5_edge",
                                 MMG5_DEL_MEM(mesh,hash->item);
                                 return 0);

              ph = &hash->item[hash->nxt];

              for (kk=ph->nxt; kk<hash->max; kk++)
                hash->item[kk].nxt = kk+1;
            }

            hash->nxt = ph->nxt;
            ph->a = ia;
            ph->b = ib;
            ph->k = 6*k + l;
            ph->nxt = 0;
            ++ph->s;
            break;
          }
          else
            ph = &hash->item[ph->nxt];
        }
      }
    }
  }

  return 1;
}

/**
 * \param mesh pointer toward a MMG5 mesh structure.
 * \param pHash pointer to the edge hash table.
 * \return PMMG_FAILURE
 *         PMMG_SUCCESS
 *
 * Hash the parallel edges. Only use face communicators to this purpose.
 *
 */
int PMMG_hashPar_fromFaceComm( PMMG_pParMesh parmesh,MMG5_HGeom *pHash ) {
  PMMG_pGrp    grp = &parmesh->listgrp[0];
  MMG5_pMesh   mesh = grp->mesh;
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  PMMG_pInt_comm int_face_comm;
  MMG5_int       na,ie;
  int          i,ifac,j,ia,i1,i2;

  assert( parmesh->ngrp == 1 );

  /** Allocation of hash table to store parallel edges */
  na = (MMG5_int)(mesh->np*0.2); // Euler-Poincare

  if ( 1 != MMG5_hNew( mesh, pHash, na, 3 * na ) ) return PMMG_FAILURE;

  /** Store parallel edges */
  for( i = 0; i < grp->nitem_int_face_comm; i++ ) {
    ie   =  grp->face2int_face_comm_index1[i]/12;
    ifac = (grp->face2int_face_comm_index1[i]%12)/3;
    pt = &mesh->tetra[ie];
    for ( j=0; j<3; j++ ) {
      ia = MMG5_iarf[ifac][j];
      /* Get edge vertices and hash it */
      i1 = MMG5_iare[ia][0];
      i2 = MMG5_iare[ia][1];
      MMG5_hEdge( mesh,pHash,pt->v[i1],pt->v[i2],0,MG_PARBDY);
    }
  }

  return PMMG_SUCCESS;
}


/**
 * \param mesh pointer toward a MMG5 mesh structure.
 * \param pHash pointer to the edge hash table.
 *
 * \return PMMG_FAILURE
 *         PMMG_SUCCESS
 *
 * Hash the edges belonging to parallel faces and store their tags with the
 * additionnal MG_PARBDY tag.
 *
 * \remark Use the assumption that all paralle edges are
 * seen by a MG_PARBDY face on an xtetra.
 *
 */
int PMMG_hashParTag_fromXtet( MMG5_pMesh mesh,MMG5_HGeom *pHash ) {
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
        /* Get edge vertices and hash it */
        i1 = MMG5_iare[ia][0];
        i2 = MMG5_iare[ia][1];
        MMG5_hEdge( mesh,pHash,pt->v[i1],pt->v[i2],0,MG_NOTAG);
        /* Tag edge and nodes as parallel */
        pxt->tag[ia] |= MG_PARBDY;
        mesh->point[pt->v[i1]].tag |= MG_PARBDY;
        mesh->point[pt->v[i2]].tag |= MG_PARBDY;
        MMG5_hTag( pHash,pt->v[i1],pt->v[i2],pxt->edg[ia],pxt->tag[ia] );
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
  uint16_t     tag;
  int8_t       i,i1,i2;


  assert ( !mesh->htab.geom );

  /* Hash the MG_PARBDY edges */
  if( PMMG_hashParTag_fromXtet(mesh,&hash) != PMMG_SUCCESS ) return PMMG_FAILURE;

  /** Update xtetra edge tag if needed */
  for (k=1; k<=mesh->ne; ++k) {

    pt = &mesh->tetra[k];
    if ( !pt->xt )
      continue;

    pxt = &mesh->xtetra[pt->xt];
    for ( i=0; i<6; ++i ) {
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

/**
 * \param mesh pointer to the mesh structure.
 * \return 1 if success, 0 otherwise.
 *
 * - Remove double triangles from tria array.
 *
 * - Remove triangles that do not belong to a boundary (non opnbdy mode) from
 *   tria array.
 *
 * - Check the matching between actual and given number of faces in the mesh:
 * Count the number of faces in mesh and compare this number to the number of
 *   given triangles.
 *
 * - If the founded number exceed the given one, add the missing
 *   boundary triangles (call to MMG5_bdryTria). Do nothing otherwise.
 *
 * - Fill the adjacency relationship between prisms and tetra (fill adjapr with
 *   a negative value to mark this special faces).
 *
 * - Set to required the triangles at interface betwen prisms and tet.
 *
 */
int PMMG_chkBdryTria(MMG5_pMesh mesh, MMG5_int* permtria) {
  MMG5_int       ntmesh,ntpres;
  int            ier;
  MMG5_Hash      hashElt;

  /** Step 1: scan the mesh and count the boundaries */
  ier = MMG5_chkBdryTria_countBoundaries(mesh,&ntmesh,&ntpres);

  /** Step 2: detect the extra boundaries (that will be ignored) provided by the
   * user */
  if ( mesh->nt ) {
    ier = MMG5_chkBdryTria_hashBoundaries(mesh,ntmesh,&hashElt);
    // Travel through the tria, flag those that are not in the hash tab or
    // that are stored more that once.
    ier = MMG5_chkBdryTria_flagExtraTriangles(mesh,&ntpres,&hashElt);
    // Delete flagged triangles
    ier = MMG5_chkBdryTria_deleteExtraTriangles(mesh, permtria);
  }
  ntmesh +=ntpres;

  /** Step 3: add the missing boundary triangles or, if the mesh contains
   * prisms, set to required the triangles at interface betwen prisms and tet */
  ier = MMG5_chkBdryTria_addMissingTriangles(mesh,ntmesh,ntpres);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param hash pointer toward the hash table of edges.
 * \param a index of the first extremity of the edge.
 * \param b index of the second extremity of the edge.
 *
 * \param s If ls mode: 1 for a parallel edge that belongs
 * to at least one element whose reference has to be splitted (either because we
 * are not in multi-mat mode or because the reference is split in multi-mat
 * mode). To avoid useless checks, some non parallel edges may be marked.
 * If the edge belongs only to non-split references, s has to be 0.
 *
 * \return PMMG_SUCCESS if success, PMMG_FAILURE if fail (edge is not found).
 *
 * Update the value of the s field stored along the edge \f$[a;b]\f$
 *
 */
int PMMG_hashUpdate_s(MMG5_Hash *hash, MMG5_int a,MMG5_int b,MMG5_int s) {
  MMG5_hedge  *ph;
  MMG5_int     key;
  MMG5_int    ia,ib;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (MMG5_KA*(int64_t)ia + MMG5_KB*(int64_t)ib) % hash->siz;
  ph  = &hash->item[key];

  while ( ph->a ) {
    if ( ph->a == ia && ph->b == ib ) {
      ph->s = s;
      return PMMG_SUCCESS;
    }

    if ( !ph->nxt ) return PMMG_FAILURE;

    ph = &hash->item[ph->nxt];

  }

  return PMMG_FAILURE;
}

/**
 * \param hash pointer toward the hash table of edges.
 * \param a index of the first extremity of the edge.
 * \param b index of the second extremity of the edge.
 * \param k index of new point along the edge [a,b].
 * \param s If ls mode in ParMmg: index of new point in internal edge communicator;
 *          otherwise, the value stored in variable s.
 * \return PMMG_SUCCESS if success, PMMG_FAILURE if fail (edge is not found).
 *
 * Find the index of the new point stored along the edge \f$[a;b]\f$ (similar to MMG5_hashGet in mmg).
 * If ls mode in ParMmg: find the index of the new point in internal edge communicator;
 * otherwise, find the value stored in variable s.
 *
 */
MMG5_int PMMG_hashGet_all(MMG5_Hash *hash,MMG5_int a,MMG5_int b,MMG5_int *k,MMG5_int *s) {
  MMG5_hedge  *ph;
  MMG5_int    key;
  MMG5_int    ia,ib;

  if ( !hash->item ) return 0;

  ia  = MG_MIN(a,b);
  ib  = MG_MAX(a,b);
  key = (MMG5_KA*(int64_t)ia + MMG5_KB*(int64_t)ib) % hash->siz;
  ph  = &hash->item[key];

  if ( !ph->a )  return PMMG_FAILURE;
  if ( ph->a == ia && ph->b == ib )  {
    *k = ph->k;
    *s = ph->s;
    return PMMG_SUCCESS;
  }
  while ( ph->nxt ) {
    ph = &hash->item[ph->nxt];
    if ( ph->a == ia && ph->b == ib ) {
      *k = ph->k;
      *s = ph->s;
      return PMMG_SUCCESS;
    }
  }
  return PMMG_FAILURE;
}

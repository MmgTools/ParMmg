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
 * \file boulep_pmmg.c
 * \brief Functions to travel volume/surface point ball for parallel analysis.
 * \author Luca Cirrottola (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Ball travel functions on parallel points.
 *
 */

#include "parmmg.h"

/**
 * \param pt pointer toward a tetrahedron. structure.
 *
 * Function to init the index of the surface ( 0 or 1 ) that is connected to
 * a special (manifold) edge.
 *
 * \warning It uses the tetra->mark field.
 */
void PMMG_Analys_Init_SurfNormIndex( MMG5_pTetra pt ) {
  pt->mark = 0;
}

/**
 * \param pt pointer toward a tetrahedron structure.
 * \param ifac index of the local tetrahedron face the point belongs to.
 * \param iloc local index of the point on the tetrahedron face.
 *
 * Function to switch the index of the surface ( from 0 to 1 ) that is connected
 * to special (manifold) edge through point (iloc) on face (ifac) on a given
 * tetrahedron.
 * A bitwise flag is used for each local point on a local tetrahedron face,
 * because the same point can belong to two different surfaces (separated by
 * a special, manifold edge) on the same tetrahedron.
 *
 * \warning It uses the tetra->mark field.
 */
void PMMG_Analys_Set_SurfNormIndex( MMG5_pTetra pt,int ifac,int iloc ) {
  pt->mark |= (1 << (3*ifac+iloc));
}

/**
 * \param pt pointer toward a tetrahedron structure.
 * \param ifac index of the local tetrahedron face the point belongs to.
 * \param iloc local index of the point on the tetrahedron face.
 *
 * Function to get the index of the surface ( 0 or 1 ) that is connected
 * to special (manifold) edge through point (iloc) on face (ifac) on a given
 * tetrahedron.
 *
 * \warning It uses the tetra->mark field.
 */
int PMMG_Analys_Get_SurfNormalIndex( MMG5_pTetra pt,int ifac,int iloc ) {
  return pt->mark & (1 << (3*ifac+iloc));
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param hash pointer toward an allocated hash table.
 * \param start index of the starting tetrahedra.
 * \param ip local index of the point in the tetrahedra \a start.
 * \param ng pointer toward the number of ridges.
 * \param nr pointer toward the number of reference edges.
 * \return ns the number of special edges passing through ip, -1 if fail.
 *
 * Count the number of ridges and reference edges incident to
 * the vertex \a ip when ip is non-manifold.
 * \remark Same as MMG5_boulernm(), but skip edges whose extremity is flagged with
 * a rank lower than myrank.
 *
 */
int PMMG_boulernm(PMMG_pParMesh parmesh,MMG5_pMesh mesh,MMG5_Hash *hash,int start,int ip,int *ng,int *nr){
  MMG5_pTetra    pt,pt1;
  MMG5_pxTetra   pxt;
  MMG5_hedge    *ph;
  int            *adja,nump,ilist,base,cur,k,k1,ns;
  int            list[MMG3D_LMAX+2];
  int            key,ia,ib,jj,a,b;
  int8_t         j,l,i,i1;
  uint8_t        ie;

  /* reset the hash table */
  for ( k=0;  k<=hash->max; ++k ) {
    hash->item[k].a = 0;
    hash->item[k].b = 0;
  }

  for ( k=0;  k<=hash->siz; ++k ) {
    hash->item[k].nxt = 0;
  }
  for (k=hash->siz; k<hash->max; k++) {
    hash->item[k].nxt = k+1;
  }

  base = ++mesh->base;
  pt   = &mesh->tetra[start];
  nump = pt->v[ip];

  /* Store initial tetrahedron */
  pt->flag = base;
  list[0] = 4*start + ip;
  ilist = 1;

  *ng = *nr = ns = 0;

  /* Explore list and travel by adjacency through elements sharing p */
  cur = 0;
  while ( cur < ilist ) {
    k = list[cur] / 4;
    i = list[cur] % 4; // index of point p in tetra k
    pt = &mesh->tetra[k];

    /* Count the number of ridge of ref edges passing through ip. */
    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];
      for (l=0; l<3; ++l) {
        ie = MMG5_arpt[i][l];
        if( MMG5_iare[ie][0] == i )
          i1 = MMG5_iare[ie][1];
        else
          i1 = MMG5_iare[ie][0];

        /* Skip parallel boundaries that will be analyzed by another process. No
         * need to skip simple parallel edges, as there is no adjacent through
         * them. */
        if( (pxt->tag[ie] & MG_PARBDYBDY) &&
            (mesh->point[pt->v[i]].flag < parmesh->myrank) ) {
           /* do nothing */
        } else if ( MG_EDG(pxt->tag[ie]) ) {
          /* Seek if we have already seen the edge. If not, hash it and
           * increment ng or nr.*/
          a = pt->v[MMG5_iare[ie][0]];
          b = pt->v[MMG5_iare[ie][1]];
          ia  = MG_MIN(a,b);
          ib  = MG_MAX(a,b);
          key = (MMG5_KA*ia + MMG5_KB*ib) % hash->siz;
          ph  = &hash->item[key];

          if ( ph->a == ia && ph->b == ib )
            continue;
          else if ( ph->a ) {
            while ( ph->nxt && ph->nxt < hash->max ) {
              ph = &hash->item[ph->nxt];
              if ( ph->a == ia && ph->b == ib )  continue;
            }
            ph->nxt   = hash->nxt;
            ph        = &hash->item[hash->nxt];

            if ( hash->nxt >= hash->max-1 ) {
              if ( mesh->info.ddebug )
                fprintf(stderr,"\n  ## Warning: %s: memory alloc problem (edge):"
                        " %d\n",__func__,hash->max);
              MMG5_TAB_RECALLOC(mesh,hash->item,hash->max,MMG5_GAP,MMG5_hedge,
                                 "MMG5_edge",return -1);
              /* ph pointer may be false after realloc */
              ph        = &hash->item[hash->nxt];

              for (jj=ph->nxt; jj<hash->max; jj++)  hash->item[jj].nxt = jj+1;
            }
            hash->nxt = ph->nxt;
          }

          /* insert new edge */
          ph->a = ia;
          ph->b = ib;
          ph->nxt = 0;

          if ( pxt->tag[ie] & MG_GEO )
            ++(*ng);
          else if ( pxt->tag[ie] & MG_REF )
            ++(*nr);
          ++ns;
        }
      }
    }

    /* Continue to travel */
    adja = &mesh->adja[4*(k-1)+1];

    for (l=0; l<3; l++) {
      i  = MMG5_inxt3[i];
      k1 = adja[i];
      if ( !k1 )  continue;
      k1 /= 4;
      pt1 = &mesh->tetra[k1];
      if ( pt1->flag == base )  continue;
      pt1->flag = base;
      for (j=0; j<4; j++)
        if ( pt1->v[j] == nump )  break;
      assert(j<4);
      /* overflow */
      if ( ilist > MMG3D_LMAX-3 )  return 0;
      list[ilist] = 4*k1+j;
      ilist++;
    }
    cur++;
  }

  return ns;
}

/**
 * \param mesh pointer toward the mesh  structure.
 * \param start tetra index.
 * \param ip point index.
 * \param iface face index.
 * \param t computed tangent vector.
 * \return 0 if point is singular, 1 otherwise.
 *
 * Mark tetra according to the color of their vertices near a special edge.
 * \remark Modeled after MMG5_boulenm.
 *
 */
int PMMG_boulen(PMMG_pParMesh parmesh,MMG5_pMesh mesh,int start,int ip,int iface,double t[3]) {
  MMG5_pTetra   pt;
  MMG5_pPoint   p0,p1,ppt;
  double   dd,l0,l1;
  int      base,nump,nr,nnm,k,piv,na,nb,adj,nvstart,fstart,aux,ip0,ip1;
  int     *adja,color;
  int16_t  tag;
  int8_t   iopp,ipiv,indb,inda,i,isface;
  int8_t   indedg[4][4] = { {-1,0,1,2}, {0,-1,3,4}, {1,3,-1,5}, {2,4,5,-1} };

  base = ++mesh->base;
  nr  = nnm = 0;
  ip0 = ip1 = 0;

  memset(t,0x00,3*sizeof(double));

  pt   = &mesh->tetra[start];
  nump = pt->v[ip];
  k    = start;

  na   = pt->v[ip];
  nb   = pt->v[MMG5_idir[iface][MMG5_inxt2[MMG5_idirinv[iface][ip]]]];
  piv  = pt->v[MMG5_idir[iface][MMG5_iprv2[MMG5_idirinv[iface][ip]]]];

  iopp   = iface;
  fstart = 4*k+iopp;
  color = 0;
  do {

    if ( pt->xt ) {
      for ( inda=0; inda<4; inda++ ){
        if ( pt->v[inda]==na ) break;
      }
      for ( indb=0; indb<4; indb++ ){
        if ( pt->v[indb]==nb ) break;
      }
      assert( (inda < 4) && (indb < 4));
      tag = mesh->xtetra[pt->xt].tag[indedg[inda][indb]];

    }

    else  tag = 0;

    /* count special edges and switch surface color if a surface has been hit */
    if ( MG_EDG(tag) && !(tag & MG_NOM) ) {
      nr++;
      color++;
      if ( !ip0 )
        ip0 = nb;
      else
        ip1 = nb;
    } else if ( tag & MG_NOM ) {
      nnm++;
    }

    /* assign color to the surface if MG_EDG has been crossed */
    if( color % 2 ) {
      PMMG_Analys_Set_SurfNormIndex( pt,iopp,MMG5_idirinv[iopp][inda] );
    }

    /* A boundary face has been hit : change travel edge */
    aux     = nb;
    nb      = piv;
    piv     = aux;
    nvstart = k;
    adj     = k;

    /* Now unfold shell of edge (na,nb) starting from k (included) */
    do {
      k = adj;
      pt = &mesh->tetra[k];
      adja = &mesh->adja[4*(k-1)+1];
      if ( pt->flag != base ) {
        for (i=0; i<4; i++)
          if ( pt->v[i] == nump )  break;
        assert(i<4);
        pt->flag = base;
      }

      /* identification of edge number in tetra k */
      if ( !MMG3D_findEdge(mesh,pt,k,na,nb,1,NULL,&i) ) {
        fprintf(stderr,"  ## Error: rank %d, function %s: edge not found.\n",parmesh->myrank,__func__);
        return -1;
      }

      /* set sense of travel */
      if ( pt->v[ MMG5_ifar[i][0] ] == piv ) {
        adj = adja[ MMG5_ifar[i][0] ] / 4;
        ipiv = MMG5_ifar[i][1];
        iopp = MMG5_ifar[i][0];
        piv = pt->v[ipiv];
      }
      else {
        adj = adja[ MMG5_ifar[i][1] ] / 4;
        ipiv = MMG5_ifar[i][0];
        iopp = MMG5_ifar[i][1];
        piv = pt->v[ipiv];
      }
      isface = ((adja[iopp] == 0) || (mesh->tetra[adj].ref != pt->ref));
    }
    while ( adj && (adj != nvstart) && !isface );
  }
  while ( 4*k+iopp != fstart );

  assert( color == 2 );
  if ( nnm > 0 || color > 2 ) {
    fprintf(stderr,"  ## Error: rank %d, function %s: found non-manifold point.\n",parmesh->myrank,__func__);
    return 0;
  }

  assert( ip0 && ip1 );
  if ( ip0 == ip1 ) {
    fprintf(stderr,"  ## Error: rank %d, function %s: only one special edge found.\n",parmesh->myrank,__func__);
    return 0;
  }

  p0 = &mesh->point[ip0];
  p1 = &mesh->point[ip1];
  ppt = &mesh->point[nump];

  l0 = (ppt->c[0] - p0->c[0])*(ppt->c[0] - p0->c[0]) \
    + (ppt->c[1] - p0->c[1])*(ppt->c[1] - p0->c[1]) + (ppt->c[2] - p0->c[2])*(ppt->c[2] - p0->c[2]);
  l1 = (ppt->c[0] - p1->c[0])*(ppt->c[0] - p1->c[0]) \
    + (ppt->c[1] - p1->c[1])*(ppt->c[1] - p1->c[1]) + (ppt->c[2] - p1->c[2])*(ppt->c[2] - p1->c[2]);
  l0 = sqrt(l0);
  l1 = sqrt(l1);

  if ( (l0 < MMG5_EPSD2) || (l1 < MMG5_EPSD2) ) {
    t[0] = p1->c[0] - p0->c[0];
    t[1] = p1->c[1] - p0->c[1];
    t[2] = p1->c[2] - p0->c[2];
  }
  else if ( l0 < l1 ) {
    dd = l0 / l1;
    t[0] = dd*(p1->c[0] - ppt->c[0]) + ppt->c[0] - p0->c[0];
    t[1] = dd*(p1->c[1] - ppt->c[1]) + ppt->c[1] - p0->c[1];
    t[2] = dd*(p1->c[2] - ppt->c[2]) + ppt->c[2] - p0->c[2];
  }
  else {
    dd = l1 / l0;
    t[0] = dd*(p0->c[0] - ppt->c[0]) + ppt->c[0] - p1->c[0];
    t[1] = dd*(p0->c[1] - ppt->c[1]) + ppt->c[1] - p1->c[1];
    t[2] = dd*(p0->c[2] - ppt->c[2]) + ppt->c[2] - p1->c[2];
  }

  dd = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
  if ( dd > MMG5_EPSD2 ) {
    dd = 1.0 / sqrt(dd);
    t[0] *= dd;
    t[1] *= dd;
    t[2] *= dd;
  }

  return 1;

}

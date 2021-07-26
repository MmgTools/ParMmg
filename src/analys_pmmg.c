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
 * \file analys_pmmg.c
 * \brief Functions for parallel mesh analysis.
 * \author Luca Cirrottola (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Parallel mesh analysis.
 *
 */

#include "parmmg.h"

static inline
void PMMG_hashNorver_loopVariable_init( int *ie,int *ifac,int *iloc ) {
  *ie   = 1;
  *ifac = 0;
  *iloc = 0;
  return;
}

static inline
void PMMG_hashNorver_loopVariable_next( int *ie,int *ifac,int *iloc ) {
  *iloc += 1;
  if( *iloc == 3 ) {
    *iloc  = 0;
    *ifac += 1;
    if( *ifac == 4 ) {
      *ifac = 0;
      *ie  += 1;
    }
  }
  return;
}

/*
 *                ip       ip       ip      ip          (...)
 *               *        *        *       * < . . . . *
 *              /        / ^        \       .         ^
 *             /     p  /   \ d      \       . (adj) .
 *            / =>  u  /     \ o   => \  =>   .     .
 *           /        / (ie)  \ w      \       .   .
 *          v        v         \ n      v       v .
 *         *        * --------> *        *       *
 *       np+ip1   ip1            ip2    np+ip2   ip2
 */
static inline
int PMMG_hashNorver_loop( MMG5_pMesh mesh,
                          MMG5_pTetra *pt,MMG5_pxTetra *pxt,MMG5_pPoint *ppt,
                          int *ie,int *ifac,int *iloc,int *ip,int *ip1,int *ip2 ) {
  int16_t tag;

  /*
   * Triple loop on:
   * - boundary tetra,
   *   - well-oriented true boundary face,
   *     - manifold parallel face vertices.
   *
   * Exit as soon as a valid entity is found.
   */
  int8_t found = 0;
  for( int ie_ = (*ie); ie_ <= mesh->ne && !found; ie_++ ) {
    *pt = &mesh->tetra[ie_];
    if( !MG_EOK(*pt) ) continue;

    /* Stay on boundary tetra */
    if( !(*pt)->xt ) continue;
    *pxt = &mesh->xtetra[(*pt)->xt];

    /* Loop on faces */
    for( int ifac_ = (*ifac); ifac_ < 4 && !found; ifac_++ ) {
      /* Get face tag */
      tag = (*pxt)->ftag[ifac_];

      /* Skip face with reversed orientation (it will be analyzed by another
       * tetra, on this or on another process) */
      if( !MG_GET((*pxt)->ori,ifac_) ) continue;

      /* Skip internal faces */
      if( !(tag & MG_BDY) ||
          ((tag & MG_PARBDY) && !(tag & MG_PARBDYBDY)) )
        continue;

      /* Loop on face vertices */
      for( int iloc_ = (*iloc); iloc_ < 3 && !found; iloc_++ ) {
        *ip = (*pt)->v[MMG5_idir[ifac_][iloc_]];
        *ppt = &mesh->point[*ip];
        if( ((*ppt)->tag & MG_PARBDY) && !((*ppt)->tag & MG_NOM) ) {
          /* Get extremities of the upstream and downstream edges of the point */
          *ip1 = (*pt)->v[MMG5_idir[ifac_][MMG5_inxt2[iloc_]]];
          *ip2 = (*pt)->v[MMG5_idir[ifac_][MMG5_iprv2[iloc_]]];
          found = 1;
        }
      }
    }
  }

  return found;
}

static inline
int PMMG_hashNorver_sweep( PMMG_pParMesh parmesh,MMG5_pMesh mesh,
                           MMG5_HGeom *hash,MMG5_HGeom *hpar,
                           MMG5_pTetra pt,MMG5_pxTetra pxt,
                           int ifac,int iloc,int ip,int ip1,int ip2,
                           int8_t *updated,int8_t *updated_par ) {
  MMG5_pEdge pa;
  int        edg,j;
  int16_t    color_old,color_new;

  /* Get old triangle color */
  color_old = pt->mark & (1 << 3*ifac+iloc);

  /* Get upstream edge color */
  if( !MMG5_hGet( hash,ip,mesh->np+ip1,&edg,&color_new ) ) return 0;

  /* If new color differs from the old one, color tria of point ip and
   * downstream edge, and flag update if the edge color has been changed */
  if( color_new && !color_old ) {
    /* Set tria color */
    pt->mark |= 1 << 3*ifac+iloc;
    /* Set downstream edge color, without crossing ridge */
    if( !pxt->tag[MMG5_iarf[ifac][MMG5_inxt2[iloc]]] & MG_GEO ) {
      if( !MMG5_hTag( hash,ip,mesh->np+ip2,edg,color_new ) ) return 0;
      /* Mark update */
      *updated = 1;
      /* Check if the edge is parallel */
      if( MMG5_hGet(hpar,ip,ip2,&edg,&color_old) ) {
        assert( !color_old );
        /* Get node position on the edge */
        pa = &mesh->edge[edg];
        if( pa->a == ip )
          j = 0;
        else {
          assert( pa->b == ip );
          j = 1;
        }
        /* Set edge color in the internal communicator */
        parmesh->int_edge_comm->intvalues[2*(edg-1)+j] = (int)color_new;
        /* Mark update */
        *updated_par = 1;
      }
    }
  }

  return 1;
}

static inline
int PMMG_hashNorver_paredge2edge( MMG5_pMesh mesh,MMG5_HGeom *hash,
                                  int *intvalues,int idx ) {
  MMG5_pEdge pa;
  int edg,j,i[2],ip,ip1;
  int16_t color_old,color_new;

  /* Get edge and node position on edge */
  pa = &mesh->edge[idx+1];
  i[0] = pa->a;
  i[1] = pa->b;

  /* Loop on the two oriented edges */
  for( j = 0; j < 2; j++ ) {
    /* Get new color from internal communicator */
    color_new = (int8_t)intvalues[2*idx+j];
    ip  = i[j];
    ip1 = i[(j+1)%2];

    /* Get old color (if edge exists locally) */
    if( !MMG5_hGet( hash,ip,mesh->np+ip1,&edg,&color_old ) ) return 1;

    /* Update local upstream color */
    if( color_new && !color_old) {
      if( !MMG5_hTag( hash,ip,mesh->np+ip1,edg,color_new ) ) return 0;
    }
  }

  return 1;
}

int PMMG_hashNorver_locIter( PMMG_pParMesh parmesh,MMG5_pMesh mesh,
                             MMG5_HGeom *hash,MMG5_HGeom *hpar,
                             int8_t *updated_par ){
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_pPoint  ppt;
  int          ie,ifac,i,ip,ip1,ip2;
  int8_t       updated;

  /* First sweep loop */
  PMMG_hashNorver_loopVariable_init( &ie, &ifac, &i );
  while( PMMG_hashNorver_loop( mesh, &pt, &pxt, &ppt,
                               &ie, &ifac, &i, &ip, &ip1, &ip2 ) ) {

    /* Swipe upstream edge -> triangle -> downstream edge */
    if( !PMMG_hashNorver_sweep( parmesh,mesh,hash,hpar,pt,pxt,ifac,i,ip,ip1,ip2,
                                &updated,updated_par ) ) return 0;

    /* Get next counter iteration */
    PMMG_hashNorver_loopVariable_next( &ie, &ifac, &i );
  }

  /* Iterate while an edge color has been updated */
  while( updated ) {
    /* Sweep loop */
    PMMG_hashNorver_loopVariable_init( &ie, &ifac, &i );
    while( PMMG_hashNorver_loop( mesh, &pt, &pxt, &ppt,
                                 &ie, &ifac, &i, &ip, &ip1, &ip2 ) ) {

      /* Swipe upstream edge -> triangle -> downstream edge */
      if( !PMMG_hashNorver_sweep( parmesh,mesh,hash,hpar,pt,pxt,ifac,i,ip,ip1,ip2,
                                  &updated,updated_par ) ) return 0;

      /* Get next counter iteration */
      PMMG_hashNorver_loopVariable_next( &ie, &ifac, &i );
    }
  }

  return 1;
}

int PMMG_hashNorver_compExt( const void *a,const void *b ) {
  return ( *(int*)a - *(int*)b );
}

int PMMG_hashNorver_communication_ext( PMMG_pParMesh parmesh,MMG5_pMesh mesh ) {
  PMMG_pGrp      grp;
  PMMG_pExt_comm ext_node_comm;
  double         *rtosend,*rtorecv,*doublevalues;
  int            *itosend,*itorecv,*intvalues,*npshift,myshift;
  int            k,nitem,color,i,idx,j,pos,d,ip;
  MPI_Comm       comm;
  MPI_Status     status;

  assert( parmesh->ngrp == 1 );
  grp = &parmesh->listgrp[0];
  assert( grp->mesh == mesh );

  comm = parmesh->comm;
  intvalues = parmesh->int_node_comm->intvalues;
  doublevalues = parmesh->int_node_comm->doublevalues;

  /* Compute nodes global offset */
  PMMG_CALLOC(parmesh,npshift,parmesh->nprocs+1,int,"npshift", return 0);
  MPI_CHECK( MPI_Allgather(&mesh->np,1,MPI_INT,&npshift[1],1,MPI_INT,comm),
             return 0 );

  for( k = 1; k <= parmesh->nprocs; ++k )
    npshift[k] += npshift[k-1];
  myshift = npshift[parmesh->myrank];
  PMMG_DEL_MEM(parmesh,npshift,int,"npshift");


  /** Exchange values on the interfaces among procs */
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    nitem         = ext_node_comm->nitem;
    color         = ext_node_comm->color_out;

    PMMG_CALLOC(parmesh,ext_node_comm->itosend,2*nitem,int,"itosend array",
                return 0);
    PMMG_CALLOC(parmesh,ext_node_comm->itorecv,2*nitem,int,"itorecv array",
                return 0);
    PMMG_MALLOC(parmesh,ext_node_comm->rtosend,6*nitem,double,"rtosend array",
                return 0);
    PMMG_MALLOC(parmesh,ext_node_comm->rtorecv,6*nitem,double,"rtorecv array",
                return 0);
    itosend = ext_node_comm->itosend;
    itorecv = ext_node_comm->itorecv;
    rtosend = ext_node_comm->rtosend;
    rtorecv = ext_node_comm->rtorecv;

    /* Fill buffers */
    for( i = 0; i < nitem; ++i ) {
      idx  = ext_node_comm->int_comm_index[i];
      /* send all entries (their meaningfullness will be checked at receiving) */
      for( j = 0; j < 2; j++ ) {
        itosend[2*i+j] = intvalues[2*idx+j]+myshift;
        for( d = 0; d < 3; d++ )
          rtosend[6*i+3*j+d] = doublevalues[6*idx+3*j+d];
      }
    }

    /* Communication */
    MPI_CHECK(
      MPI_Sendrecv(itosend,2*nitem,MPI_INT,color,MPI_ANALYS_TAG+2,
                   itorecv,2*nitem,MPI_INT,color,MPI_ANALYS_TAG+2,
                   comm,&status),return 0 );
    MPI_CHECK(
      MPI_Sendrecv(rtosend,6*nitem,MPI_DOUBLE,color,MPI_ANALYS_TAG+3,
                   rtorecv,6*nitem,MPI_DOUBLE,color,MPI_ANALYS_TAG+3,
                   comm,&status),return 0 );
  }

  /* Fill internal communicator */
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    color = ext_node_comm->color_out;

    itorecv = ext_node_comm->itorecv;
    rtorecv = ext_node_comm->rtorecv;

    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx  = ext_node_comm->int_comm_index[i];
      /* find first free position */
      pos = 0;
      while( intvalues[2*idx+pos] && pos < 2 ) {
        pos++;
      }
      /* copy valid entries into free positions */
      j = 0;
      while( itorecv[2*i+j] && j < 2 ) {
        assert( pos < 3 );
        intvalues[2*idx+pos] = itorecv[2*i+j]-myshift;
        for( d = 0; d < 3; d++ ) {
          doublevalues[6*idx+3*pos+d] = rtorecv[6*i+3*j+d];
        }
        pos++;
        j++;
      }
    }
  }

  /* At this point all parallel manifold points should have both ridge
   * extremities */
  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    ip  = grp->node2int_node_comm_index1[i];
    idx = grp->node2int_node_comm_index2[i];

    qsort( &intvalues[2*idx], 2, sizeof(int),PMMG_hashNorver_compExt );
    mesh->point[ip].tmp = intvalues[2*idx];
  }


  /* Free memory */
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    PMMG_DEL_MEM(parmesh,ext_node_comm->itosend,int,"itosend array");
    PMMG_DEL_MEM(parmesh,ext_node_comm->itorecv,int,"itorecv array");
    PMMG_DEL_MEM(parmesh,ext_node_comm->rtosend,double,"rtosend array");
    PMMG_DEL_MEM(parmesh,ext_node_comm->rtorecv,double,"rtorecv array");
  }

  return 1;

}

int PMMG_hashNorver_communication_init( PMMG_pParMesh parmesh ) {
  PMMG_pExt_comm ext_edge_comm;
  int k,nitem,color;

  for ( k = 0; k < parmesh->next_edge_comm; ++k ) {
    ext_edge_comm = &parmesh->ext_edge_comm[k];
    nitem         = ext_edge_comm->nitem;
    color         = ext_edge_comm->color_out;

    PMMG_CALLOC(parmesh,ext_edge_comm->itosend,2*nitem,int,"itosend array",
                return 0);
    PMMG_CALLOC(parmesh,ext_edge_comm->itorecv,2*nitem,int,"itorecv array",
                return 0);
  }

  return 1;
}

int PMMG_hashNorver_communication_free( PMMG_pParMesh parmesh ) {
  PMMG_pExt_comm ext_edge_comm;
  int k;

  for ( k = 0; k < parmesh->next_edge_comm; ++k ) {
    ext_edge_comm = &parmesh->ext_edge_comm[k];
    PMMG_DEL_MEM(parmesh,ext_edge_comm->itosend,int,"itosend array");
    PMMG_DEL_MEM(parmesh,ext_edge_comm->itorecv,int,"itorecv array");
  }

  return 1;
}

int PMMG_hashNorver_communication( PMMG_pParMesh parmesh ){
  PMMG_pExt_comm ext_edge_comm;
  int            *itosend,*itorecv,*intvalues;
  int            k,nitem,color,i,idx,j;
  MPI_Comm       comm;
  MPI_Status     status;

  comm = parmesh->comm;
  intvalues = parmesh->int_edge_comm->intvalues;

  /** Exchange values on the interfaces among procs */
  for ( k = 0; k < parmesh->next_edge_comm; ++k ) {
    ext_edge_comm = &parmesh->ext_edge_comm[k];
    nitem         = ext_edge_comm->nitem;
    color         = ext_edge_comm->color_out;

    itosend = ext_edge_comm->itosend;
    itorecv = ext_edge_comm->itorecv;

    /* Fill buffers */
    for ( i=0; i<nitem; ++i ) {
      idx  = ext_edge_comm->int_comm_index[i];
      for( j = 0; j < 2; j++ ) {
        itosend[2*i+j] = intvalues[2*idx+j];
       }
    }

    /* Communication */
    MPI_CHECK(
      MPI_Sendrecv(itosend,2*nitem,MPI_INT,color,MPI_ANALYS_TAG+2,
                   itorecv,2*nitem,MPI_INT,color,MPI_ANALYS_TAG+2,
                   comm,&status),return 0 );
  }

  /* Fill internal communicator */
  for ( k = 0; k < parmesh->next_edge_comm; ++k ) {
    ext_edge_comm = &parmesh->ext_edge_comm[k];

    itorecv = ext_edge_comm->itorecv;

    for ( i=0; i<ext_edge_comm->nitem; ++i ) {
      idx  = ext_edge_comm->int_comm_index[i];

      for( j = 0; j < 2; j++ )
        intvalues[2*idx+j] |= itorecv[2*i+j];
    }
  }

  return 1;
}

int PMMG_hashNorver( PMMG_pParMesh parmesh,MMG5_pMesh mesh,MMG5_HGeom *hpar ){
  PMMG_pGrp      grp;
  PMMG_pInt_comm int_node_comm,int_edge_comm,int_face_comm;
  MMG5_pTetra    pt;
  MMG5_pxTetra   pxt;
  MMG5_pPoint    ppt;
  MMG5_HGeom     hash;
  int            ie,ifac,ia,i,ip,ip1,ip2,iter,idx;
  int            *intvalues;
  int16_t        tag;

  assert( parmesh->ngrp == 1 );
  grp = &parmesh->listgrp[0];
  assert( mesh = grp->mesh );

  int_node_comm = parmesh->int_node_comm;
  int_edge_comm = parmesh->int_edge_comm;
  int_face_comm = parmesh->int_face_comm;

  PMMG_CALLOC(parmesh,int_node_comm->doublevalues,6*int_node_comm->nitem,double,"node doublevalues",return 0);
  PMMG_CALLOC(parmesh,int_node_comm->intvalues,2*int_node_comm->nitem,int,"node intvalues",return 0);
  PMMG_CALLOC(parmesh,int_edge_comm->intvalues,2*int_edge_comm->nitem,int,"edge intvalues",return 0);


  /* Store internal communicator index on the point itself */
  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    ip  = grp->node2int_node_comm_index1[i];
    idx = grp->node2int_node_comm_index2[i];
    ppt = &mesh->point[ip];
    ppt->tmp = idx;
  }


  /* Hash table used to store edges touching a parallel point.
   * Assume that in the worst case each parallel faces has the three edges in
   * the table, plus two other internal edges. */
  if ( !MMG5_hNew(mesh,&hash,3*int_face_comm->nitem,5*int_face_comm->nitem) )
    return 0;

  /** 1) Loop on edges touching an old parallel point and insert them in the
   *     hash table. Find local ridge extremities. */
  intvalues = int_node_comm->intvalues;
  PMMG_hashNorver_loopVariable_init( &ie, &ifac, &i );
  while( PMMG_hashNorver_loop( mesh, &pt, &pxt, &ppt,
                               &ie, &ifac, &i, &ip, &ip1, &ip2 ) ) {

    /* Store extremity of the upstream edge ip->ip1 if it is a ridge */
    if( pxt->tag[MMG5_iarf[ifac][MMG5_iprv2[i]]] & MG_GEO ) {
      /* store extremity in the first free position in the internal node
       * communicator */
      if( !intvalues[2*ppt->tmp] )
        intvalues[2*ppt->tmp] = ip1;
      else {
        intvalues[2*ppt->tmp+1] = ip1;
      }
    }

    /* Store upstream edge ip0->ip1 */
    if( !MMG5_hEdge( mesh,&hash,ip,mesh->np+ip1,0,0 ) ) return 0;

    /* Store downstream edge ip0->ip2 */
    if( !MMG5_hEdge( mesh,&hash,ip,mesh->np+ip2,0,0 ) ) return 0;


    /* Get next counter iteration */
    PMMG_hashNorver_loopVariable_next( &ie, &ifac, &i );
  }

  /** 2) Parallel exchange of ridge extremities, and update color on second
   *     extremity. */
  if( !PMMG_hashNorver_communication_ext( parmesh,mesh ) ) return 0;

  PMMG_hashNorver_loopVariable_init( &ie, &ifac, &i );
  while( PMMG_hashNorver_loop( mesh, &pt, &pxt, &ppt,
                               &ie, &ifac, &i, &ip, &ip1, &ip2 ) ) {

    /* Switch edge color if its extremity is found */
    if( ip1 == ppt->tmp )
      if( !MMG5_hTag( &hash,ip,mesh->np+ip1,0,1 ) ) return 0;

    /* Get next counter iteration */
    PMMG_hashNorver_loopVariable_next( &ie, &ifac, &i );
  }


  /** 3) Propagate surface colors:
   *     - Do local iterations to converge on local colors.
   *     - Do parallel iterations to:
   *       - communicate colors,
   *       - get colors from parallel edges, and
   *       - propagate colors through local iterations.
   */
  int8_t updated_par;

  /* Reset tetra mark */
  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    if( !MG_EOK(pt) ) continue;
    pt->mark = 0;
  }

  /* 3.1) Local update iterations */
  if( !PMMG_hashNorver_locIter( parmesh,mesh,&hash,hpar,&updated_par ) ) return 0;

  /* 3.2) Parallel update iterations */
  if( !PMMG_hashNorver_communication_init( parmesh ) ) return 0;
  while( updated_par ) {

    /* 3.2.1) Parallel communication */
    if( !PMMG_hashNorver_communication( parmesh ) ) return 0;


    /* 3.2.2) Get color from parallel edges */
    for( i = 0; i < grp->nitem_int_edge_comm; i++ ){
      idx = grp->edge2int_edge_comm_index2[i];
      if( !PMMG_hashNorver_paredge2edge( mesh,&hash,int_edge_comm->intvalues,
            idx ) ) return 0;
    }

    /* 3.2.3) Local update iterations */
    if( !PMMG_hashNorver_locIter( parmesh,mesh,&hash,hpar,&updated_par ) ) return 0;
  }
  if( !PMMG_hashNorver_communication_free( parmesh ) ) return 0;


  /* Free memory */
  MMG5_DEL_MEM(mesh,hash.geom);
  PMMG_DEL_MEM(parmesh,int_edge_comm->intvalues,int,"edge intvalues");
  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,int,"node intvalues");
  PMMG_DEL_MEM(parmesh,int_node_comm->doublevalues,double,"node doublevalues");

  return 1;
}

/**
 * \param parmesh pointer to the parmesh structure
 * \param mesh pointer to the mesh structure
 *
 * \return 1 if success, 0 if failure.
 *
 * Compute continuous geometric support (normal and tangent vectors) on
 * non-manifold MG_OLDPARBDY points.
 *
 * \remark Analogous to the MMG3D_nmgeom function, but it only travels on
 * old parallel points.
 * \remark Normal and tangent vectors on these points are overwritten.
 *
 */
int PMMG_update_nmgeom(PMMG_pParMesh parmesh,MMG5_pMesh mesh){
  MMG5_pTetra     pt;
  MMG5_pPoint     p0;
  MMG5_pxPoint    pxp;
  int             k,base;
  int             *adja;
  double          n[3],t[3];
  int             ip;
  int8_t          i,j,ier;

  for( ip = 1; ip <= mesh->np; ip++ ) {
    mesh->point[ip].flag = mesh->base;
  }

  base = ++mesh->base;
  for (k=1; k<=mesh->ne; k++) {
    pt   = &mesh->tetra[k];
    if( !MG_EOK(pt) ) continue;
    adja = &mesh->adja[4*(k-1)+1];
    for (i=0; i<4; i++) {
      if ( adja[i] ) continue;
      for (j=0; j<3; j++) {
        ip = MMG5_idir[i][j];
        p0 = &mesh->point[pt->v[ip]];
        if ( p0->flag == base )  continue;
        else if ( !(p0->tag & MG_OLDPARBDY) ) continue;
        else if ( !(p0->tag & MG_NOM) )  continue;

        p0->flag = base;
        ier = MMG5_boulenm(mesh,k,ip,i,n,t);

        if ( ier < 0 )
          return 0;
        else if ( !ier ) {
          p0->tag |= MG_REQ;
          p0->tag &= ~MG_NOSURF;
        }
        else {
          if ( !p0->xp ) {
            ++mesh->xp;
            if(mesh->xp > mesh->xpmax){
              MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,MMG5_GAP,MMG5_xPoint,
                                 "larger xpoint table",
                                 mesh->xp--;
                                 fprintf(stderr,"  Exit program.\n");return 0;);
            }
            p0->xp = mesh->xp;
          }
          pxp = &mesh->xpoint[p0->xp];
          memcpy(pxp->n1,n,3*sizeof(double));
          memcpy(p0->n,t,3*sizeof(double));
        }
      }
    }
  }
  /* Deal with the non-manifold points that do not belong to a surface
   * tetra (a tetra that has a face without adjacent)*/
  for (k=1; k<=mesh->ne; k++) {
    pt   = &mesh->tetra[k];
    if( !MG_EOK(pt) ) continue;
    
    for (i=0; i<4; i++) {
      p0 = &mesh->point[pt->v[i]];
      if ( !(p0->tag & MG_OLDPARBDY) ) continue;
      else if ( p0->tag & MG_PARBDY || p0->tag & MG_REQ || !(p0->tag & MG_NOM) || p0->xp ) continue;
      ier = MMG5_boulenmInt(mesh,k,i,t);
      if ( ier ) {
        ++mesh->xp;
        if(mesh->xp > mesh->xpmax){
          MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,MMG5_GAP,MMG5_xPoint,
                            "larger xpoint table",
                            mesh->xp--;
                            fprintf(stderr,"  Exit program.\n");return 0;);
        }
        p0->xp = mesh->xp;
        pxp = &mesh->xpoint[p0->xp];
        memcpy(p0->n,t,3*sizeof(double));
      }
      else {
        p0->tag |= MG_REQ;
        p0->tag &= ~MG_NOSURF;
      }
    }
  }

  /*for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    if ( !(p0->tag & MG_NOM) || p0->xp ) continue;
    p0->tag |= MG_REQ;
    p0->tag &= ~MG_NOSURF;
  }*/
  
  return 1;
}

/**
 * \param mesh pointer toward the mesh  structure.
 * \return 0 if point is singular, 1 otherwise.
 *
 * Repristinate singularity tags on MG_OLDPARBDY points, based on the tags
 * of incident edges.
 *
 */
static inline
int PMMG_update_singul(PMMG_pParMesh parmesh,MMG5_pMesh mesh) {
  MMG5_pTetra         ptet;
  MMG5_pPoint         ppt;
  MMG5_Hash           hash;
  int                 k,i;
  int                 nc, nre, ng, nrp,ier;

  /* Second: seek the non-required non-manifold points and try to analyse
   * whether they are corner or required. */

  /* Hash table used by boulernm to store the special edges passing through
   * a given point */
  if ( ! MMG5_hashNew(mesh,&hash,mesh->np,(int)(3.71*mesh->np)) ) return 0;

  nc = nre = 0;
  ++mesh->base;
  for (k=1; k<=mesh->ne; ++k) {
    ptet = &mesh->tetra[k];
    if ( !MG_EOK(ptet) ) continue;

    for ( i=0; i<4; ++i ) {
      ppt = &mesh->point[ptet->v[i]];

      /* Skip non-previously-parallel points */
      if ( !(ppt->tag & MG_OLDPARBDY) ) continue;

      if ( (!MG_VOK(ppt)) || (ppt->flag==mesh->base)  ) continue;
      ppt->flag = mesh->base;

      if ( (!MG_EDG(ppt->tag)) || MG_SIN(ppt->tag) ) continue;

      ier = MMG5_boulernm(mesh,&hash, k, i, &ng, &nrp);
      if ( ier < 0 ) return 0;
      else if ( !ier ) continue;

      if ( (ng+nrp) > 2 ) {
        ppt->tag |= MG_CRN + MG_REQ;
        ppt->tag &= ~MG_NOSURF;
        nre++;
        nc++;
      }
      else if ( (ng == 1) && (nrp == 1) ) {
        ppt->tag |= MG_REQ;
        ppt->tag &= ~MG_NOSURF;
        nre++;
      }
      else if ( ng == 1 && !nrp ){
        ppt->tag |= MG_CRN + MG_REQ;
        ppt->tag &= ~MG_NOSURF;
        nre++;
        nc++;
      }
      else if ( ng == 1 && !nrp ){
        ppt->tag |= MG_CRN + MG_REQ;
        ppt->tag &= ~MG_NOSURF;
        nre++;
        nc++;
      }
    }
  }

  /* Free the edge hash table */
  MMG5_DEL_MEM(mesh,hash.item);

  if ( mesh->info.ddebug || abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"     %d corner and %d required vertices added\n",nc,nre);

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * Update continuous geometry information (normal and tangent vectors) on
 * previouly parallel points, as this infomation has not been built during
 * initial analysis.
 *
 * The algorithm analyzes manifold surfaces and gives colors for each of the
 * vertices seen by a triangles, depending whether a special edge has been
 * crossed when travelling the surface ball of the vertex (color 1)  or not
 * (color 0).
 * This color is used to know whether the triangle normal vector should
 * contribute to the first (0) average normal vector of the point, or the second
 * (1) if the point is a special edge.
 *
 * \remark Normal and tangent vectors are overzritten on these points.
 *
 */
int PMMG_update_norver( PMMG_pParMesh parmesh,MMG5_pMesh mesh ) {
  MMG5_pTetra    pt;
  MMG5_Tria      tt;
  MMG5_pPoint    ppt;
  MMG5_pxTetra   pxt;
  MMG5_pxPoint   pxp;
  double         n[3],dd;
  int            *adja,ip,ie,ifac,i,iloc,d,base,k1;

  base = ++mesh->base;

  /* Reset points flag and source element fields */
  for( ip = 1; ip <= mesh->np; ip++ ) {
    ppt = &mesh->point[ip];
    ppt->flag = 0;
    ppt->s = 0;
  }

  /* Reset tetra mark: it will be used by PMMG_boulen to bitwise flag the
   * element vertices depending on the portion of surface near a MG_EDG edge.
   * Initialize point source element field with tetra, face, and local index
   * on tetra */
  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    if( !MG_EOK(pt) || !pt->xt ) continue;
    PMMG_Analys_Init_SurfNormIndex( pt );

    pxt = &mesh->xtetra[pt->xt];

    for( ifac = 0; ifac < 4; ifac++ ) {

      /* skip faces not on boundary or already seen depending on orientation */
      if( !(pxt->ftag[ifac] & MG_BDY) || !MG_GET(pxt->ori,ifac) ) continue;

      for( i = 0; i < 3; i++ ) {
        iloc = MMG5_idir[ifac][i];
        ip = pt->v[iloc];
        ppt = &mesh->point[ip];
        ppt->s = 16*ie+4*ifac+iloc;
      }
    }
  }

  /* Initialize flags on non-parallel, non-singular, manifold, old interface
   * points.
   * Reset xpoint or reference/allocate it if needed.
   * Flag tetra vertices based on the portion of surface they see near a
   * NG_EDG edge. */
  for( ip = 1; ip <= mesh->np; ip++ ) {
    ppt = &mesh->point[ip];

    if( !MG_VOK(ppt) || !(ppt->tag & MG_BDY) ) continue;

    /* Flag points to analyze */
    if( !(ppt->tag & MG_OLDPARBDY) || ppt->tag & MG_PARBDY ||
        MG_SIN(ppt->tag) || ppt->tag & MG_NOM ) continue;

    ppt->flag = base;

    /* Reset normal vector or allocate/reference new xpoint */
    if( ppt->xp ) {
      pxp = &mesh->xpoint[ppt->xp];
      memset(pxp->n1,0x00,3*sizeof(double));
      memset(pxp->n2,0x00,3*sizeof(double));
    } else {
      ++mesh->xp;
      if(mesh->xp > mesh->xpmax){
        MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,MMG5_GAP,MMG5_xPoint,
                          "larger xpoint table",
                          mesh->xp--;return 0;);
      }
      ppt->xp = mesh->xp;
    }

    /* color portion of surface touching the special point */
    if( MG_EDG(ppt->tag) ) {
      ie   =  ppt->s / 16;
      ifac = (ppt->s % 16) / 4;
      iloc = (ppt->s % 16) % 4;
      assert( ie );
      if( !PMMG_boulen(parmesh,mesh,ie,iloc,ifac,ppt->n) ) {
        fprintf(stderr,"  ## Error: rank %d, function %s: failed boulen.\n",parmesh->myrank,__func__);
        return 0;
      }
    }
  }


  /* compute and attribute normal vectors */
  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];

    if( !MG_EOK(pt) || !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];

    for( ifac = 0; ifac < 4; ifac++ ) {

      /* skip faces not on boundary or already seen depending on orientation */
      if( !(pxt->ftag[ifac] & MG_BDY) || !MG_GET(pxt->ori,ifac) ) continue;


      /* get virtual triangle and compute normal */
      MMG5_tet2tri(mesh,ie,ifac,&tt);
      MMG5_nortri(mesh,&tt,n);
      assert( n[0]*n[0]+n[1]*n[1]+n[2]*n[2] > MMG5_EPSD2 );

      /* loop on face vertices and set contribution if the vertex is flagged */
      for( i = 0; i < 3; i++ ) {
        iloc = MMG5_idir[ifac][i];
        ip = pt->v[iloc];
        ppt = &mesh->point[ip];

        if( ppt->flag != base ) continue;
        assert( !(pxt->ftag[ifac] & MG_PARBDY ) );

        assert( ppt->xp );
        pxp = &mesh->xpoint[ppt->xp];

        if( PMMG_Analys_Get_SurfNormalIndex( pt,ifac,i ) ) {
          for( d = 0; d < 3; d++ )
            pxp->n2[d] += n[d];
        } else {
          for( d = 0; d < 3; d++ )
            pxp->n1[d] += n[d];
        }
      }
    }
  }

  /* normalize normals */
  for( ip = 1; ip <= mesh->np; ip++ ) {
    ppt = &mesh->point[ip];
    if( ppt->flag != base ) continue;

    pxp = &mesh->xpoint[ppt->xp];

    /* normalize first normal vector */
    dd = 0.0;
    for( d = 0; d < 3; d++ )
      dd += pxp->n1[d]*pxp->n1[d];
    if ( dd <= MMG5_EPSD2 ) {
      fprintf(stderr,"  ## Error: rank %d, function %s: computed null normal vector for first surface.\n",parmesh->myrank,__func__);
      return 0;
    }
    dd = 1.0 / sqrt(dd);
    for( d = 0; d < 3; d++ )
      pxp->n1[d] *= dd;

    /* normalize second normal vector if present */
    if( MG_EDG(ppt->tag) ) {
      dd = 0.0;
      for( d = 0; d < 3; d++ )
        dd += pxp->n2[d]*pxp->n2[d];
      if ( dd <= MMG5_EPSD2 ) {
        fprintf(stderr,"  ## Error: rank %d, function %s: computed null normal vector for second surface.\n",parmesh->myrank,__func__);
        return 0;
      }
      dd = 1.0 / sqrt(dd);
      for( d = 0; d < 3; d++ )
        pxp->n2[d] *= dd;

      if( ppt->tag & MG_GEO ) {
        /* update tangent as cross product of normal vectors if geometric edge
         * (this cannot be done on a reference edge, as it can lie on a C1
         * surface) */
        ppt->n[0] = pxp->n1[1]*pxp->n2[2] - pxp->n1[2]*pxp->n2[1];
        ppt->n[1] = pxp->n1[2]*pxp->n2[0] - pxp->n1[0]*pxp->n2[2];
        ppt->n[2] = pxp->n1[0]*pxp->n2[1] - pxp->n1[1]*pxp->n2[0];
        dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
        if ( dd > MMG5_EPSD2 ) {
          dd = 1.0 / sqrt(dd);
          ppt->n[0] *= dd;
          ppt->n[1] *= dd;
          ppt->n[2] *= dd;
        }
      }
    }
  }

#ifndef NDEBUG
  for( ip = 1; ip <= mesh->np; ip++ ) {
    ppt = &mesh->point[ip];
    if( !MG_VOK(ppt) || !(ppt->tag & MG_BDY) ) continue;
    if( ppt->tag & MG_OLDPARBDY && !(ppt->tag & MG_PARBDY) &&
        !MG_SIN(ppt->tag) && !(ppt->tag & MG_NOM) ) {
      assert(ppt->xp);
      pxp = &mesh->xpoint[ppt->xp];
      assert( pxp->n1[0]*pxp->n1[0]+pxp->n1[1]*pxp->n1[1]+pxp->n1[2]*pxp->n1[2] > 0.0 );
      if( MG_EDG(ppt->tag) )
        assert( pxp->n2[0]*pxp->n2[0]+pxp->n2[1]*pxp->n2[1]+pxp->n2[2]*pxp->n2[2] > 0.0 );
    }
  }
#endif

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * Update continuous geometry information (normal and tangent vectors) on
 * previouly parallel points, as this infomation has not been built during
 * initial analysis.
 * Also update singular tags as they were erased on parallel points.
 *
 */
int PMMG_update_analys(PMMG_pParMesh parmesh) {
  PMMG_pGrp  grp;
  MMG5_pMesh mesh;
  int        igrp,ip;

  /* Loop on groups */
  for( igrp = 0; igrp < parmesh->ngrp; igrp++ ) {
    grp = &parmesh->listgrp[igrp];
    mesh = grp->mesh;

    for( ip = 1; ip <= mesh->np; ip++ )
      mesh->point[ip].flag = mesh->base;

    /* create tetra adjacency */
    if ( !MMG3D_hashTetra(mesh,0) ) {
      fprintf(stderr,"\n  ## Hashing problem (1). Exit program.\n");
      return 0;
    }

    /* build normal and tangent vectors on previously parallel points */
    if( !PMMG_update_norver(parmesh,mesh) ) {
      fprintf(stderr,"  ## Error: rank %d, function %s: cannot update normal vectors on surface.\n",
              parmesh->myrank,__func__);
      return 0;
    }

    /* repristinate singularity tags on points that were previously parallel
     * (so, their required tags have been erased) */
    if( !PMMG_update_singul(parmesh,mesh) ) {
      fprintf(stderr,"  ## Error: rank %d, function %s: cannot update singularities on surface.\n",
              parmesh->myrank,__func__);
      return 0;
    }

    /* define geometry for non manifold points */
    if( !PMMG_update_nmgeom(parmesh,mesh) ) {
      fprintf(stderr,"  ## Error: rank %d, function %s: cannot update geometric support on non-manifold edges.\n",
              parmesh->myrank,__func__);
      return 0;
    }
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * \return 1 if success, 0 if fail.
 * Seek the non-required non-manifold points and try to analyse whether they are
 * corner or required.
 *
 * \remark We don't know how to travel through the shell of a non-manifold point
 * by triangle adjacency. Thus the work done here can't be performed in the \ref
 * MMG5_singul function.
 * \remark Modeled after the sequential MMG5_setVertexNmTag function.
 *
 */
static inline
int PMMG_setVertexNmTag(PMMG_pParMesh parmesh,MMG5_pMesh mesh) {
  PMMG_pGrp      grp;
  PMMG_pInt_comm int_node_comm;
  PMMG_pExt_comm ext_node_comm;
  MPI_Comm       comm;
  MPI_Status     status;
  MMG5_pTetra    pt;
  MMG5_pPoint    ppt;
  MMG5_Hash      hash;
  int            nc,xp,nr,ns0,ns1,nre;
  int            ip,idx,iproc,k,i,j,d;
  int            nitem,color,tag;
  int            *intvalues,*itosend,*itorecv,*iproc2comm;
  double         *doublevalues,*rtosend,*rtorecv;

  comm   = parmesh->comm;
  assert( parmesh->ngrp == 1 );
  grp = &parmesh->listgrp[0];
  int_node_comm = parmesh->int_node_comm;

  /* Allocate intvalues to accomodate xp and nr for each point, doublevalues to
   * accomodate two edge vectors at most. */
  PMMG_CALLOC(parmesh,int_node_comm->intvalues,2*int_node_comm->nitem,int,"intvalues",return 0);
  PMMG_CALLOC(parmesh,int_node_comm->doublevalues,6*int_node_comm->nitem,double,"doublevalues",return 0);
  intvalues    = int_node_comm->intvalues;
  doublevalues = int_node_comm->doublevalues;

  /* Store source tetra for every boundary point */
  for( ip = 1; ip <= mesh->np; ip++ ) {
    ppt = &mesh->point[ip];
    if( !MG_VOK(ppt) ) continue;
    ppt->flag = 0;
  }
  for( k = 1; k <= mesh->ne; k++ ) {
    pt = &mesh->tetra[k];
    for( i = 0; i < 4; i++ ) {
      ppt = &mesh->point[pt->v[i]];
      if( ppt->flag ) continue;
      ppt->s = 4*k+i;
      ppt->flag++;
    }
  }

  /* Array to reorder communicators */
  PMMG_MALLOC(parmesh,iproc2comm,parmesh->nprocs,int,"iproc2comm",return 0);

  for( iproc = 0; iproc < parmesh->nprocs; iproc++ )
    iproc2comm[iproc] = PMMG_UNSET;

  for( k = 0; k < parmesh->next_node_comm; k++ ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    iproc = ext_node_comm->color_out;
    iproc2comm[iproc] = k;
  }

  /* Flag parallel points with the lowest rank they see in order to analyse
   * them only once. */
  for( iproc = parmesh->nprocs-1; iproc >= 0; iproc-- ) {
    k = iproc2comm[iproc];
    if( k == PMMG_UNSET ) continue;
    ext_node_comm = &parmesh->ext_node_comm[k];
    for( i = 0; i < ext_node_comm->nitem; i++ ) {
      idx = ext_node_comm->int_comm_index[i];
      intvalues[idx] = iproc;
    }
  }

  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    ip  = grp->node2int_node_comm_index1[i];
    idx = grp->node2int_node_comm_index2[i];
    ppt = &mesh->point[ip];
    if( !MG_VOK(ppt) ) continue;
    ppt->flag = intvalues[idx];
    /* Reset intvalues in order to reuse it to count special edges */
    intvalues[idx] = 0;
  }

  /* Hash table used by boulernm to store the special edges passing through
   * a given point */
  if ( ! MMG5_hashNew(mesh,&hash,mesh->np,(int)(3.71*mesh->np)) ) return 0;


  /** Local singularity analysis */
  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    ip  = grp->node2int_node_comm_index1[i];
    idx = grp->node2int_node_comm_index2[i];
    ppt = &mesh->point[ip];
    if ( !MG_VOK(ppt) || ( ppt->tag & MG_CRN ) )
      continue;
    else if ( MG_EDG(ppt->tag) ) {
      /* Count the number of ridges passing through the point (xp) and the
       * number of ref edges (nr).
       * Edges on communicators only once as they are flagged with their
       * lowest seen rank. */
      ns0 = PMMG_boulernm(parmesh, mesh, &hash, ppt->s/4, ppt->s%4, &xp, &nr);
      assert( ns0 == xp+nr );

      /* Add nb of ridges/refs to intvalues */
      intvalues[2*idx]   = xp;
      intvalues[2*idx+1] = nr;

    }
  }

  /** Exchange values on the interfaces among procs */
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    nitem         = ext_node_comm->nitem;
    color         = ext_node_comm->color_out;

    PMMG_CALLOC(parmesh,ext_node_comm->itosend,2*nitem,int,"itosend array",
                return 0);
    PMMG_CALLOC(parmesh,ext_node_comm->itorecv,2*nitem,int,"itorecv array",
                return 0);
    itosend = ext_node_comm->itosend;
    itorecv = ext_node_comm->itorecv;

    /* Fill buffers */
    for ( i=0; i<nitem; ++i ) {
      idx  = ext_node_comm->int_comm_index[i];
      for( j = 0; j < 2; j++ ) {
        itosend[2*i+j] = intvalues[2*idx+j];
      }
    }

    MPI_CHECK(
      MPI_Sendrecv(itosend,2*nitem,MPI_INT,color,MPI_ANALYS_TAG,
                   itorecv,2*nitem,MPI_INT,color,MPI_ANALYS_TAG,
                   comm,&status),return 0 );

  }

  /** First pass: Sum nb. of singularities, Store received edge vectors in
   *  doublevalues if there is room for them.
   */
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    itorecv = ext_node_comm->itorecv;
    rtorecv = ext_node_comm->rtorecv;

    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx  = ext_node_comm->int_comm_index[i];

      /* Current nb. of singularities == current nb. of stored edge vectors */
      ns0 = intvalues[2*idx]+intvalues[2*idx+1];

      /* Increment nb. of singularities */
      intvalues[2*idx]   += itorecv[2*i];
      intvalues[2*idx+1] += itorecv[2*i+1];

    }
  }

  /** Second pass: Analysis.
   *  Flag intvalues[2*idx]   as PMMG_UNSET if the point is required.
   *  Flag intvalues[2*idx+1] as PMMG_UNSET if the point is corner.
   */
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    itosend = ext_node_comm->itosend;
    itorecv = ext_node_comm->itorecv;
    rtosend = ext_node_comm->rtosend;
    rtorecv = ext_node_comm->rtorecv;

    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx  = ext_node_comm->int_comm_index[i];

      xp = intvalues[2*idx];
      nr = intvalues[2*idx+1];
      if ( (xp == PMMG_UNSET) || (nr == PMMG_UNSET) )  continue;

      if ( !(xp+nr) )  continue;
      if ( (xp+nr) > 2 ) {
        intvalues[2*idx]   = PMMG_UNSET;
        intvalues[2*idx+1] = PMMG_UNSET;
        nre++;
        nc++;
      }
      else if ( (xp == 1) && (nr == 1) ) {
        intvalues[2*idx]   = PMMG_UNSET;
        nre++;
      }
      else if ( xp == 1 && !nr ){
        intvalues[2*idx]   = PMMG_UNSET;
        intvalues[2*idx+1] = PMMG_UNSET;
        nre++;
        nc++;
      }
      else if ( nr == 1 && !xp ){
        intvalues[2*idx]   = PMMG_UNSET;
        intvalues[2*idx+1] = PMMG_UNSET;
        nre++;
        nc++;
      }
    }
  }

  /** Third pass: Tag points */
  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    ip  = grp->node2int_node_comm_index1[i];
    idx = grp->node2int_node_comm_index2[i];
    ppt = &mesh->point[ip];
    if( intvalues[2*idx]   == PMMG_UNSET ) { /* is required */
      ppt->tag |= MG_REQ;
      ppt->tag &= ~MG_NOSURF;
    }
    if( intvalues[2*idx+1] == PMMG_UNSET ) { /* is corner */
      ppt->tag |= MG_CRN;
    }
  }

  /* Free the edge hash table */
  MMG5_DEL_MEM(mesh,hash.item);

  if ( mesh->info.ddebug || abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"     %d corner and %d required vertices added\n",nc,nre);

  PMMG_DEL_MEM(parmesh,iproc2comm,int,"iproc2comm");
  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,int,"intvalues");
  PMMG_DEL_MEM(parmesh,int_node_comm->doublevalues,double,"doublevalues");
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    PMMG_DEL_MEM(parmesh,ext_node_comm->itosend,int,"itosend array");
    PMMG_DEL_MEM(parmesh,ext_node_comm->itorecv,int,"itorecv array");
    PMMG_DEL_MEM(parmesh,ext_node_comm->rtosend,double,"rtosend array");
    PMMG_DEL_MEM(parmesh,ext_node_comm->rtorecv,double,"rtorecv array");
  }

  return 1;
}

/**
 * \param mesh pointer towar the mesh structure.
 * \param hash edges hash table.
 * \return 1 if success, 0 if failed.
 *
 * Set tags to non-manifold edges and vertices. Not done before because we need
 * the \ref MMG5_xTetra table.
 *
 * \warning if fail, the edge hash table \a hash is not freed.
 *
 */
int PMMG_setNmTag(PMMG_pParMesh parmesh, MMG5_pMesh mesh, MMG5_Hash *hash) {

//  /* First: seek edges at the interface of two distinct domains and mark it as
//   * required */
//  if ( !MMG5_setEdgeNmTag(mesh,hash) ) return 0;

  /* Second: seek the non-required non-manifold points and try to analyse
   * whether they are corner or required. */
  if ( !PMMG_setVertexNmTag(parmesh,mesh) ) return 0;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param mesh pointer toward the mesh structure.
 * \param adjt pointer toward the table of triangle adjacency.
 * \param start index of triangle where we start to work.
 * \param ip index of vertex on which we work.
 * \param list pointer toward the computed list of GEO vertices incident to \a ip.
 * \param listref pointer toward the corresponding edge references
 * \param ng pointer toward the number of ridges.
 * \param nr pointer toward the number of reference edges.
 * \param lmax maxmum size for the ball of the point \a ip.
 * \return The number of edges incident to the vertex \a ip.
 *
 * Store edges and count the number of ridges and reference edges incident to
 * the vertex \a ip.
 * \remark Same as MMG5_bouler(), but skip edges whose extremity is flagged with
 * a rank lower than myrank.
 *
 */
int PMMG_bouler(PMMG_pParMesh parmesh,MMG5_pMesh mesh,int *adjt,int start,int ip,
                 int *list,int *listref,int *ng,int *nr,int lmax) {
  MMG5_pTria    pt;
  int           *adja,k,ns;
  int           i,i1,i2;

  pt  = &mesh->tria[start];
  if ( !MG_EOK(pt) )  return 0;

  /* check other triangle vertices */
  k  = start;
  i  = ip;
  *ng = *nr = ns = 0;

  do {
    i1 = MMG5_inxt2[i];
    i2 = MMG5_iprv2[i];
    /* Skip parallel boundaries that will be analyzed by another process. No
     * need to skip simple parallel edges, as there is no adjacent through
     * them. */
    if( (pt->tag[i1] & MG_PARBDYBDY || pt->tag[i1] & MG_BDY) &&
        (mesh->point[pt->v[i2]].flag < parmesh->myrank) ) {
      /* do nothing */
    } else {
      if ( MG_EDG(pt->tag[i1]) ) {
        if ( pt->tag[i1] & MG_GEO )
          *ng = *ng + 1;
        else if ( pt->tag[i1] & MG_REF )
          *nr = *nr + 1;
        ns++;
        list[ns] = pt->v[i2];
        listref[ns] = pt->edg[i1];
        if ( ns > lmax-2 )  return -ns;
      }
    }
    adja = &adjt[3*(k-1)+1];
    k  = adja[i1] / 3;
    i  = adja[i1] % 3;
    i  = MMG5_inxt2[i];
    pt = &mesh->tria[k];
  }
  while ( k && k != start );

  /* reverse loop */
  if ( k != start ) {
    k = start;
    i = ip;
    do {
      pt = &mesh->tria[k];
      i2 = MMG5_iprv2[i];
      /* Skip parallel boundaries that will be analyzed by another process. No
       * need to skip simple parallel edges, as there is no adjacent through
       * them. */
      if( (pt->tag[i2] & MG_PARBDYBDY || pt->tag[i2] & MG_BDY) &&
          (mesh->point[pt->v[i2]].flag < parmesh->myrank) ) {
        /* do nothing */
      } else {
        if ( MG_EDG(pt->tag[i2]) ) {
          i1 = MMG5_inxt2[i];
          if ( pt->tag[i2] & MG_GEO )
            *ng = *ng + 1;
          else if ( pt->tag[i2] & MG_REF )
            *nr = *nr + 1;
          ns++;
          list[ns] = pt->v[i1];
          listref[ns] = pt->edg[i2];
          if ( ns > lmax-2 )  return -ns;
        }
      }
      adja = &adjt[3*(k-1)+1];
      k = adja[i2] / 3;
      i = adja[i2] % 3;
      i = MMG5_iprv2[i];
    }
    while ( k && k != start );
  }

  return ns;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param mesh pointer toward the mesh structure
 *
 * Check for singularities.
 * \remark Modeled after the MMG5_singul function.
 */
int PMMG_singul(PMMG_pParMesh parmesh,MMG5_pMesh mesh) {
  PMMG_pGrp      grp;
  PMMG_pInt_comm int_node_comm;
  PMMG_pExt_comm ext_node_comm;
  MPI_Comm       comm;
  MPI_Status     status;
  MMG5_pTria     pt;
  MMG5_pPoint    ppt,p1;
  double         ux,uy,uz,vx,vy,vz,dd;
  int            list[MMG3D_LMAX+2],listref[MMG3D_LMAX+2],nc,xp,nr,ns0,ns1,nre;
  int            ip,idx,iproc,k,i,j,d;
  int            nitem,color,tag;
  int            *intvalues,*itosend,*itorecv,*iproc2comm;
  double         *doublevalues,*rtosend,*rtorecv;

  comm   = parmesh->comm;
  assert( parmesh->ngrp == 1 );
  grp = &parmesh->listgrp[0];
  int_node_comm = parmesh->int_node_comm;

  /* Allocate intvalues to accomodate xp and nr for each point, doublevalues to
   * accomodate two edge vectors at most. */
  PMMG_CALLOC(parmesh,int_node_comm->intvalues,2*int_node_comm->nitem,int,"intvalues",return 0);
  PMMG_CALLOC(parmesh,int_node_comm->doublevalues,6*int_node_comm->nitem,double,"doublevalues",return 0);
  intvalues    = int_node_comm->intvalues;
  doublevalues = int_node_comm->doublevalues;

  /* Store source triangle for every boundary point */
  for( ip = 1; ip <= mesh->np; ip++ ) {
    ppt = &mesh->point[ip];
    if( !MG_VOK(ppt) || !(ppt->tag & MG_PARBDY ) ) continue;
    ppt->flag = 0;
    ppt->s = 0;
  }
  for( k = 1; k <= mesh->nt; k++ ) {
    pt = &mesh->tria[k];
    /* give a valid source triangle (not a PARBDY where no adjacency is
     * provided) */
    tag = pt->tag[0] & pt->tag[1] & pt->tag[2];
    if ( !MG_EOK(pt) || ((tag & MG_PARBDY) && !(tag & MG_PARBDYBDY)) )  continue;
    for( i = 0; i < 3; i++ ) {
      ppt = &mesh->point[pt->v[i]];
      if( ppt->flag ) continue;
      ppt->s = 3*k+i;
      ppt->flag++;
    }
  }

  /* Array to reorder communicators */
  PMMG_MALLOC(parmesh,iproc2comm,parmesh->nprocs,int,"iproc2comm",return 0);

  for( iproc = 0; iproc < parmesh->nprocs; iproc++ )
    iproc2comm[iproc] = PMMG_UNSET;

  for( k = 0; k < parmesh->next_node_comm; k++ ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    iproc = ext_node_comm->color_out;
    iproc2comm[iproc] = k;
  }

  /* Flag parallel points with the lowest rank they see in order to analyse
   * them only once. */
  for( iproc = parmesh->nprocs-1; iproc >= 0; iproc-- ) {
    k = iproc2comm[iproc];
    if( k == PMMG_UNSET ) continue;
    ext_node_comm = &parmesh->ext_node_comm[k];
    for( i = 0; i < ext_node_comm->nitem; i++ ) {
      idx = ext_node_comm->int_comm_index[i];
      intvalues[idx] = iproc;
    }
  }

  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    ip  = grp->node2int_node_comm_index1[i];
    idx = grp->node2int_node_comm_index2[i];
    ppt = &mesh->point[ip];
    if( !MG_VOK(ppt) ) continue;
    ppt->flag = intvalues[idx];
    /* Reset intvalues in order to reuse it to count special edges */
    intvalues[idx] = 0;
  }


  /** Exchange point tags that could have been changed by PMMG_setdhd
   *  (expecially for PARBDYBDY points that don't belong to any PARBDYBDY
   *  triangle on the current proc) */

  /* Fill communicator */
  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    ip  = grp->node2int_node_comm_index1[i];
    idx = grp->node2int_node_comm_index2[i];
    ppt = &mesh->point[ip];
    intvalues[idx] = ppt->tag;
  }

  /* Allocate buffers with the size needed by the singularity analysis, but only
   * use part of theme here to exchange tags */
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    nitem         = ext_node_comm->nitem;
    color         = ext_node_comm->color_out;

    PMMG_CALLOC(parmesh,ext_node_comm->itosend,2*nitem,int,"itosend array",
                return 0);
    PMMG_CALLOC(parmesh,ext_node_comm->itorecv,2*nitem,int,"itorecv array",
                return 0);
    PMMG_CALLOC(parmesh,ext_node_comm->rtosend,6*nitem,double,"rtosend array",
                return 0);
    PMMG_CALLOC(parmesh,ext_node_comm->rtorecv,6*nitem,double,"rtorecv array",
                return 0);
    itosend = ext_node_comm->itosend;
    itorecv = ext_node_comm->itorecv;
    rtosend = ext_node_comm->rtosend;
    rtorecv = ext_node_comm->rtorecv;

    /* Fill buffers */
    for ( i=0; i<nitem; ++i ) {
      idx  = ext_node_comm->int_comm_index[i];
      itosend[i] = intvalues[idx];
    }

    MPI_CHECK(
      MPI_Sendrecv(itosend,nitem,MPI_INT,color,MPI_ANALYS_TAG,
                   itorecv,nitem,MPI_INT,color,MPI_ANALYS_TAG,
                   comm,&status),return 0 );
  }

  /* Get tags and reset buffers and communicator */
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    nitem         = ext_node_comm->nitem;
    color         = ext_node_comm->color_out;

    itosend = ext_node_comm->itosend;
    itorecv = ext_node_comm->itorecv;
    rtosend = ext_node_comm->rtosend;
    rtorecv = ext_node_comm->rtorecv;

    /* Fill buffers */
    for ( i=0; i<nitem; ++i ) {
      idx  = ext_node_comm->int_comm_index[i];
      intvalues[idx] |= itorecv[i];
      itosend[i] = 0;
      itorecv[i] = 0;
    }
  }

  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    ip  = grp->node2int_node_comm_index1[i];
    idx = grp->node2int_node_comm_index2[i];
    ppt = &mesh->point[ip];
    ppt->tag |= intvalues[idx];
    intvalues[idx] = 0;
  }


  /** Local singularity analysis */
  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    ip  = grp->node2int_node_comm_index1[i];
    idx = grp->node2int_node_comm_index2[i];
    ppt = &mesh->point[ip];

    /* skip points that do not belong to a true boundary triangle on this proc */
    if( !ppt->s ) continue;

    if ( !MG_VOK(ppt) || ( ppt->tag & MG_CRN ) || ( ppt->tag & MG_NOM ) )
      continue;
    else if ( MG_EDG(ppt->tag) ) {
      /* Count the number of ridges passing through the point (xp) and the
       * number of ref edges (nr).
       * Edges on communicators only once as they are flagged with their
       * lowest seen rank. */
      ns0 = PMMG_bouler(parmesh,mesh,mesh->adjt,ppt->s/3,ppt->s%3,list,listref,&xp,&nr,MMG3D_LMAX);
      assert( ns0 == xp+nr );

      /* Add nb of ridges/refs to intvalues */
      intvalues[2*idx]   = xp;
      intvalues[2*idx+1] = nr;

      /* Go to next point if too many singularities */
      if ( ns0 > 2 ) continue;

      /* Add edge vectors to doublevalues */
      for( j = 0; j < ns0; j++ ) {
        p1 = &mesh->point[list[j+1]];
        for( d = 0; d < 3; d++ )
          doublevalues[6*idx+3*j+d] = p1->c[d]-ppt->c[d];
      }
    }
  }

  /** Exchange values on the interfaces among procs */
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    nitem         = ext_node_comm->nitem;
    color         = ext_node_comm->color_out;

    itosend = ext_node_comm->itosend;
    itorecv = ext_node_comm->itorecv;
    rtosend = ext_node_comm->rtosend;
    rtorecv = ext_node_comm->rtorecv;

    /* Fill buffers */
    for ( i=0; i<nitem; ++i ) {
      idx  = ext_node_comm->int_comm_index[i];
      for( j = 0; j < 2; j++ ) {
        itosend[2*i+j] = intvalues[2*idx+j];
        /* Take every value, its meaning will be evaluated at recv */
        for( d = 0; d < 3; d++ )
          rtosend[6*i+3*j+d] = doublevalues[6*idx+3*j+d];
      }
    }

    MPI_CHECK(
      MPI_Sendrecv(itosend,2*nitem,MPI_INT,color,MPI_ANALYS_TAG,
                   itorecv,2*nitem,MPI_INT,color,MPI_ANALYS_TAG,
                   comm,&status),return 0 );

    MPI_CHECK(
      MPI_Sendrecv(rtosend,6*nitem,MPI_DOUBLE,color,MPI_ANALYS_TAG+1,
                   rtorecv,6*nitem,MPI_DOUBLE,color,MPI_ANALYS_TAG+1,
                   comm,&status),return 0 );
  }

  /** First pass: Sum nb. of singularities, Store received edge vectors in
   *  doublevalues if there is room for them.
   */
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    itorecv = ext_node_comm->itorecv;
    rtorecv = ext_node_comm->rtorecv;

    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx  = ext_node_comm->int_comm_index[i];

      /* Current nb. of singularities == current nb. of stored edge vectors */
      ns0 = intvalues[2*idx]+intvalues[2*idx+1];

      /* Increment nb. of singularities */
      intvalues[2*idx]   += itorecv[2*i];
      intvalues[2*idx+1] += itorecv[2*i+1];

      /* Go to next if too many singularities */
      ns1 = intvalues[2*idx]+intvalues[2*idx+1];
      if( ns1 > 2 ) continue;

      if( ns0 < 2 ) { /* If there is room for another edge vector */
        for( j = 0; j < ns1-ns0; j++ ) { /* Loop on newly received vectors */
          for( d = 0; d < 3; d++ ) { /* Store them in doublevalues */
            doublevalues[6*idx+3*(ns0+j)+d] = rtorecv[6*i+3*j+d];
          }
        }
      }
    }
  }

  /** Second pass: Analysis.
   *  Flag intvalues[2*idx]   as PMMG_UNSET if the point is required.
   *  Flag intvalues[2*idx+1] as PMMG_UNSET if the point is corner.
   */
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    itosend = ext_node_comm->itosend;
    itorecv = ext_node_comm->itorecv;
    rtosend = ext_node_comm->rtosend;
    rtorecv = ext_node_comm->rtorecv;

    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx  = ext_node_comm->int_comm_index[i];

      xp = intvalues[2*idx];
      nr = intvalues[2*idx+1];
      if ( (xp == PMMG_UNSET) || (nr == PMMG_UNSET) )  continue;

      if ( !(xp+nr) )  continue;
      if ( (xp+nr) > 2 ) {
        intvalues[2*idx]   = PMMG_UNSET;
        intvalues[2*idx+1] = PMMG_UNSET;
        nre++;
        nc++;
      }
      else if ( (xp == 1) && (nr == 1) ) {
        intvalues[2*idx]   = PMMG_UNSET;
        nre++;
      }
      else if ( xp == 1 && !nr ){
        intvalues[2*idx]   = PMMG_UNSET;
        intvalues[2*idx+1] = PMMG_UNSET;
        nre++;
        nc++;
      }
      else if ( nr == 1 && !xp ){
        intvalues[2*idx]   = PMMG_UNSET;
        intvalues[2*idx+1] = PMMG_UNSET;
        nre++;
        nc++;
      }
      /* check ridge angle */
      else {
        ux = doublevalues[6*idx];
        uy = doublevalues[6*idx+1];
        uz = doublevalues[6*idx+2];
        vx = doublevalues[6*idx+3];
        vy = doublevalues[6*idx+4];
        vz = doublevalues[6*idx+5];
        dd = (ux*ux + uy*uy + uz*uz) * (vx*vx + vy*vy + vz*vz);
        if ( fabs(dd) > MMG5_EPSD ) {
          dd = (ux*vx + uy*vy + uz*vz) / sqrt(dd);
          if ( dd > -mesh->info.dhd ) {
            intvalues[2*idx+1] = PMMG_UNSET;
            nc++;
          }
        }
      }
    }
  }

  /** Third pass: Tag points */
  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    ip  = grp->node2int_node_comm_index1[i];
    idx = grp->node2int_node_comm_index2[i];
    ppt = &mesh->point[ip];
    if( intvalues[2*idx]   == PMMG_UNSET ) { /* is required */
      ppt->tag |= MG_REQ;
      ppt->tag &= ~MG_NOSURF;
    }
    if( intvalues[2*idx+1] == PMMG_UNSET ) { /* is corner */
      ppt->tag |= MG_CRN;
    }
  }

  if ( abs(mesh->info.imprim) > 3 && nre > 0 )
    fprintf(stdout,"     %d corners, %d singular points detected\n",nc,nre);


  PMMG_DEL_MEM(parmesh,iproc2comm,int,"iproc2comm");
  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,int,"intvalues");
  PMMG_DEL_MEM(parmesh,int_node_comm->doublevalues,double,"doublevalues");
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    PMMG_DEL_MEM(parmesh,ext_node_comm->itosend,int,"itosend array");
    PMMG_DEL_MEM(parmesh,ext_node_comm->itorecv,int,"itorecv array");
    PMMG_DEL_MEM(parmesh,ext_node_comm->rtosend,double,"rtosend array");
    PMMG_DEL_MEM(parmesh,ext_node_comm->rtorecv,double,"rtorecv array");
  }

  return 1;
}

/**
 * \param parmesh pointer to the parmesh structure
 * \param mesh pointer to the mesh structure
 * \param pHash pointer to the parallel edges hash table
 *
 * \return 1 if success, 0 if failure.
 *
 * Check dihedral angle to detect ridges on parallel edges.
 *
 * The integer Communicator is dimensioned to store the number of triangles seen
 * by a parallel edge on each partition, and a "flag" to check the references of
 * the seen triangles. This "flag" is initialized with the first seen triangle,
 * and it is switched to PMMG_UNSET if the second reference differs.
 * The double communicator is dimensioned to store at most two normal vectors
 * per edge. If the dihedral angle check detects a ridge, the edge "flag" is
 * switched to 2*PMMG_UNSET if the edge is not a reference edge, or to
 * 3*PMMG_UNSET if the edge is also a reference edge.
 */
int PMMG_setdhd(PMMG_pParMesh parmesh,MMG5_pMesh mesh,MMG5_HGeom *pHash ) {
  PMMG_pGrp      grp;
  PMMG_pInt_comm int_edge_comm;
  PMMG_pExt_comm ext_edge_comm;
  MMG5_pTetra    pt;
  MMG5_pTria     ptr;
  int            *intvalues,*itorecv,*itosend;
  double         *doublevalues,*rtorecv,*rtosend;
  int            nitem,color,nt0,nt1;
  double         n1[3],n2[3],dhd;
  int            *adja,k,kk,ne,nr,nm,j;
  int            i,ii,i1,i2;
  int            idx,edg,d;
  int16_t        tag;
  MPI_Comm       comm;
  MPI_Status     status;

  assert( parmesh->ngrp == 1 );
  grp = &parmesh->listgrp[0];
  assert( mesh == grp->mesh );

  comm = parmesh->comm;
  int_edge_comm = parmesh->int_edge_comm;

  /* Allocate edge intvalues to tag non-manifold and reference edges */
  PMMG_CALLOC(parmesh,int_edge_comm->intvalues,2*int_edge_comm->nitem,int,
              "intvalues",return 0);
  intvalues = int_edge_comm->intvalues;

  /* Allocate edge doublevalues to store triangles normals */
  PMMG_CALLOC(parmesh,int_edge_comm->doublevalues,6*int_edge_comm->nitem,double,
              "doublevalues",return 0);
  doublevalues = int_edge_comm->doublevalues;


  /** Loop on true boundary triangles and store the normal in the edge internal
   *  communicator where the triangle touches a parallel edge. */
  for( k = 1; k <= mesh->nt; k++ ) {
    ptr = &mesh->tria[k];
    if( !MG_EOK(ptr) )  continue;

    /* Skip faces that are just parallel or analyzed by another process */
    tag = ptr->tag[0] & ptr->tag[1] & ptr->tag[2];
    if( (tag & MG_PARBDY) && !(tag & MG_PARBDYBDY) ) continue;

    /* triangle normal */
    MMG5_nortri(mesh,ptr,n1);

    /* Get parallel edge touched by a boundary face and store normal vectors */
    for (i=0; i<3; i++) {
      /* Skip non-manifold edges */
      if ( (ptr->tag[i] & MG_NOM) ) continue;

      i1 = MMG5_inxt2[i];
      i2 = MMG5_inxt2[i1];
      if ( !MMG5_hGet( pHash, ptr->v[i1], ptr->v[i2], &edg, &tag ) ) continue;
      idx = edg-1;

      /* Count how many times the edge is seen locally */
      intvalues[2*idx]++;
      /* Do not store anything else for non-manifold */
      if( intvalues[2*idx] > 2 ) continue;
      /* Store reference if the edge is visited for the first time */
      if( intvalues[2*idx] == 1 )
        intvalues[2*idx+1] = ptr->ref;
      /* Check for ref edge if the edge is visited for the second time, and
       * unset the reference if a ref edge is found */
      if( intvalues[2*idx] == 2 )
        if( intvalues[2*idx+1] != ptr->ref )
          intvalues[2*idx+1] = PMMG_UNSET;
      /* Store the normal in the free position */
      for( d = 0; d < 3; d++ )
        doublevalues[6*idx+3*(intvalues[2*idx]-1)+d] = n1[d];

    }
  }

  /** Exchange values on the interfaces among procs */
  for ( k = 0; k < parmesh->next_edge_comm; ++k ) {
    ext_edge_comm = &parmesh->ext_edge_comm[k];
    nitem         = ext_edge_comm->nitem;
    color         = ext_edge_comm->color_out;

    PMMG_CALLOC(parmesh,ext_edge_comm->itosend,2*nitem,int,"itosend array",
                return 0);
    PMMG_CALLOC(parmesh,ext_edge_comm->itorecv,2*nitem,int,"itorecv array",
                return 0);
    PMMG_CALLOC(parmesh,ext_edge_comm->rtosend,6*nitem,double,"rtosend array",
                return 0);
    PMMG_CALLOC(parmesh,ext_edge_comm->rtorecv,6*nitem,double,"rtorecv array",
                return 0);
    itosend = ext_edge_comm->itosend;
    itorecv = ext_edge_comm->itorecv;
    rtosend = ext_edge_comm->rtosend;
    rtorecv = ext_edge_comm->rtorecv;

    /* Fill buffers */
    for ( i=0; i<nitem; ++i ) {
      idx  = ext_edge_comm->int_comm_index[i];
      for( j = 0; j < 2; j++ ) {
        itosend[2*i+j] = intvalues[2*idx+j];
        /* Take every value, its meaning will be evaluated at recv */
        for( d = 0; d < 3; d++ )
          rtosend[6*i+3*j+d] = doublevalues[6*idx+3*j+d];
      }
    }

    MPI_CHECK(
      MPI_Sendrecv(itosend,2*nitem,MPI_INT,color,MPI_ANALYS_TAG+2,
                   itorecv,2*nitem,MPI_INT,color,MPI_ANALYS_TAG+2,
                   comm,&status),return 0 );

    MPI_CHECK(
      MPI_Sendrecv(rtosend,6*nitem,MPI_DOUBLE,color,MPI_ANALYS_TAG+3,
                   rtorecv,6*nitem,MPI_DOUBLE,color,MPI_ANALYS_TAG+3,
                   comm,&status),return 0 );
  }

  /** First pass: Increment the number of seen triangles, check for reference
   *  edges and mark them with PMMG_UNSET, and store new triangles normals if
   *  there is room for them. */
  for ( k = 0; k < parmesh->next_edge_comm; ++k ) {
    ext_edge_comm = &parmesh->ext_edge_comm[k];

    itorecv = ext_edge_comm->itorecv;
    rtorecv = ext_edge_comm->rtorecv;

    for ( i=0; i<ext_edge_comm->nitem; ++i ) {
      idx  = ext_edge_comm->int_comm_index[i];

      /* Current nb. of triangles */
      nt0 = intvalues[2*idx];

      /* Accumulate the number of triangles touching the edge in intvalues */
      intvalues[2*idx] += itorecv[2*i];

      /* Go to the next if non-manifold */
      nt1 = intvalues[2*idx];
      if( nt1 > 2 ) continue;

      /* Check new ref */
      if( (nt1 == 2) && (nt0 == 1) ) {
        if( intvalues[2*idx+1] != itorecv[2*i+1] )
          intvalues[2*idx+1] = PMMG_UNSET;
      }

      if( nt0 < 2 ) { /* If there is room for another normal vector */
        for( j = 0; j < nt1-nt0; j++ ) { /* Loop on newly received vectors */
          for( d = 0; d < 3; d++ ) { /* Store them in doublevalues */
            doublevalues[6*idx+3*(nt0+j)+d] = rtorecv[6*i+3*j+d];
          }
        }
      }
    }
  }

  /** Second pass: Check dihedral angle and mark geometric edge with
   *  2*PMMG_UNSET (3x if it is already a reference edge) */
  for ( k = 0; k < parmesh->next_edge_comm; ++k ) {
    ext_edge_comm = &parmesh->ext_edge_comm[k];

    for ( i=0; i<ext_edge_comm->nitem; ++i ) {
      idx  = ext_edge_comm->int_comm_index[i];

      nt1 = intvalues[2*idx];

      if( nt1 == 2 ) {
        for( d = 0; d < 3; d++ ) {
          n1[d] = doublevalues[6*idx+d];
          n2[d] = doublevalues[6*idx+3+d];
        }
        dhd = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2];
        if ( dhd <= mesh->info.dhd ) {
          if( intvalues[2*idx+1] != PMMG_UNSET )
            intvalues[2*idx+1] = 2*PMMG_UNSET;
          else
            intvalues[2*idx+1] = 3*PMMG_UNSET;
        }
      }
    }
  }

  /** Third pass: Loop on triangles to tag edges and points */
  ne = nr = nm = 0;
  for( k = 1; k <= mesh->nt; k++ ) {
    ptr = &mesh->tria[k];
    if( !MG_EOK(ptr) )  continue;

    /* Skip faces that are just parallel */
    tag = ptr->tag[0] & ptr->tag[1] & ptr->tag[2];
    if( (tag & MG_PARBDY) && !(tag & MG_PARBDYBDY || tag & MG_BDY) )
      continue;

    /* Get parallel edge touched by a MG_BDY face and store normal vectors */
    for (i=0; i<3; i++) {
      /* Skip non-manifold edges */
      if ( (ptr->tag[i] & MG_NOM) ) continue;

      i1 = MMG5_inxt2[i];
      i2 = MMG5_inxt2[i1];
      if ( !MMG5_hGet( pHash, ptr->v[i1], ptr->v[i2], &edg, &tag ) ) continue;
      idx = edg-1;

      if( intvalues[2*idx] == 1 ) { /* no adjacent */
        ptr->tag[i] |= MG_GEO;
        i1 = MMG5_inxt2[i];
        i2 = MMG5_inxt2[i1];
        mesh->point[ptr->v[i1]].tag |= MG_GEO;
        mesh->point[ptr->v[i2]].tag |= MG_GEO;
        nr++;
      } else {
        if( (intvalues[2*idx+1] ==   PMMG_UNSET) ||
            (intvalues[2*idx+1] == 3*PMMG_UNSET) ) { /* reference edge */
          ptr->tag[i]   |= MG_REF;
          i1 = MMG5_inxt2[i];
          i2 = MMG5_inxt2[i1];
          mesh->point[ptr->v[i1]].tag |= MG_REF;
          mesh->point[ptr->v[i2]].tag |= MG_REF;
          ne++;
        }
        if( (intvalues[2*idx+1] == 2*PMMG_UNSET) ||
            (intvalues[2*idx+1] == 3*PMMG_UNSET) ) { /* geometric edge */
          ptr->tag[i]   |= MG_GEO;
          i1 = MMG5_inxt2[i];
          i2 = MMG5_inxt2[i1];
          mesh->point[ptr->v[i1]].tag |= MG_GEO;
          mesh->point[ptr->v[i2]].tag |= MG_GEO;
          nr++;
        }
        if( intvalues[2*idx] > 2 ) { /* non-manifold edge */
          ptr->tag[i] |= MG_GEO + MG_NOM;
          i1 = MMG5_inxt2[i];
          i2 = MMG5_inxt2[i1];
          mesh->point[ptr->v[i1]].tag |= MG_GEO + MG_NOM;
          mesh->point[ptr->v[i2]].tag |= MG_GEO + MG_NOM;
          nm++;
        }
      }
    }
  }


  if ( abs(mesh->info.imprim) > 3 && nr > 0 )
    fprintf(stdout,"     %d ridges, %d edges updated\n",nr,ne);

  PMMG_DEL_MEM(parmesh,int_edge_comm->intvalues,int,"intvalues");
  PMMG_DEL_MEM(parmesh,int_edge_comm->doublevalues,double,"doublevalues");
  for ( k = 0; k < parmesh->next_edge_comm; ++k ) {
    ext_edge_comm = &parmesh->ext_edge_comm[k];
    PMMG_DEL_MEM(parmesh,ext_edge_comm->itosend,int,"itosend array");
    PMMG_DEL_MEM(parmesh,ext_edge_comm->itorecv,int,"itorecv array");
    PMMG_DEL_MEM(parmesh,ext_edge_comm->rtosend,double,"rtosend array");
    PMMG_DEL_MEM(parmesh,ext_edge_comm->rtorecv,double,"rtorecv array");
  }

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param mesh pointer toward the mesh structure
 * \return 0 if fail, 1 if success.
 *
 * Check all boundary triangles.
 */
int PMMG_analys_tria(PMMG_pParMesh parmesh,MMG5_pMesh mesh) {

  /**--- stage 1: data structures for surface */
  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"\n  ** SURFACE ANALYSIS\n");

  /* create tetra adjacency */
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Hashing problem (1). Exit program.\n");
    return 0;
  }

  /* create prism adjacency */
  if ( !MMG3D_hashPrism(mesh) ) {
    fprintf(stderr,"\n  ## Prism hashing problem. Exit program.\n");
    return 0;
  }
  /* compatibility triangle orientation w/r tetras */
  if ( !MMG5_bdryPerm(mesh) ) {
    fprintf(stderr,"\n  ## Boundary orientation problem. Exit program.\n");
    return 0;
  }

  /* identify surface mesh */
  if ( !MMG5_chkBdryTria(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    return 0;
  }
  MMG5_freeXTets(mesh);
  MMG5_freeXPrisms(mesh);

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param mesh pointer toward the mesh structure
 *
 * \remark Modeled after the MMG3D_analys function, it doesn't deallocate the
 * tria structure in order to be able to build communicators.
 */
int PMMG_analys(PMMG_pParMesh parmesh,MMG5_pMesh mesh) {
  MMG5_Hash      hash;
  MMG5_HGeom     hpar;
  size_t         myavailable,oldMemMax;

  /* Tag parallel triangles on material interfaces as boundary */
  if( !PMMG_parbdyTria( parmesh ) ) {
    fprintf(stderr,"\n  ## Unable to recognize parallel triangles on material interfaces. Exit program.\n");
    return 0;
  }


  /* Set surface triangles to required in nosurf mode or for parallel boundaries */
  MMG3D_set_reqBoundaries(mesh);


  /* create surface adjacency */
  if ( !MMG3D_hashTria(mesh,&hash) ) {
    MMG5_DEL_MEM(mesh,hash.item);
    fprintf(stderr,"\n  ## Hashing problem (2). Exit program.\n");
    return 0;
  }

  /* build hash table for geometric edges */
  if ( !MMG5_hGeom(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem (0). Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    MMG5_DEL_MEM(mesh,mesh->htab.geom);
    return 0;
  }

  /**--- stage 2: surface analysis */
  if ( abs(mesh->info.imprim) > 5  || mesh->info.ddebug )
    fprintf(stdout,"  ** SETTING TOPOLOGY\n");

  /* identify connexity */
  if ( !MMG5_setadj(mesh) ) {
    fprintf(stderr,"\n  ## Topology problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /* Hash parallel edges */
  if( PMMG_hashPar_pmmg( parmesh,&hpar ) != PMMG_SUCCESS ) return 0;

  /* Build edge communicator */
  if( !PMMG_build_edgeComm( parmesh,mesh,&hpar ) ) return 0;

  /* check for ridges: check dihedral angle using adjacent triangle normals */
  if ( mesh->info.dhd > MMG5_ANGLIM && !MMG5_setdhd(mesh) ) {
    fprintf(stderr,"\n  ## Geometry problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  if ( mesh->info.dhd > MMG5_ANGLIM && !PMMG_setdhd( parmesh,mesh,&hpar ) ) {
    fprintf(stderr,"\n  ## Geometry problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /* identify singularities on interior points */
  if ( !MMG5_singul(mesh) ) {
    fprintf(stderr,"\n  ## MMG5_singul problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }


  /* identify singularities on parallel points */
  if ( !PMMG_singul(parmesh,mesh) ) {
    fprintf(stderr,"\n  ## PMMG_singul problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }


  if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug )
    fprintf(stdout,"  ** DEFINING GEOMETRY\n");

  /* define (and regularize) normals: create xpoints */
  if ( !MMG5_norver( mesh ) ) {
    fprintf(stderr,"\n  ## Normal problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /* set bdry entities to tetra: create xtetra and set references */
  if ( !MMG5_bdrySet(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    return 0;
  }

  /* set non-manifold edges sharing non-intersecting multidomains as required */
  if ( abs(mesh->info.imprim) > 5  || mesh->info.ddebug )
    fprintf(stdout,"  ** UPDATING TOPOLOGY AT NON-MANIFOLD POINTS\n");

  if ( !MMG5_setNmTag(mesh,&hash) ) {
    fprintf(stderr,"\n  ## Non-manifold topology problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    return 0;
  }

  if ( !PMMG_setNmTag(parmesh,mesh,&hash) ) {
    fprintf(stderr,"\n  ## Non-manifold topology problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    return 0;
  }


  /* check subdomains connected by a vertex and mark these vertex as corner and required */
#warning Luca: check that parbdy are skipped
  MMG5_chkVertexConnectedDomains(mesh);

  /* build hash table for geometric edges */
  if ( !mesh->na && !MMG5_hGeom(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem (0). Exit program.\n");
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    MMG5_DEL_MEM(mesh,mesh->htab.geom);
    return 0;
  }

  /* Update edges tags and references for xtetras */
  if ( !MMG5_bdryUpdate(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    return 0;
  }

  /* define geometry for non manifold points */
  if ( !MMG3D_nmgeom(mesh) ) return 0;


#ifdef USE_POINTMAP
  /* Initialize source point with input index */
  int ip;
  for( ip = 1; ip <= mesh->np; ip++ )
    mesh->point[ip].src = ip;
#endif

  /* release memory */
  PMMG_edge_comm_free( parmesh );
  MMG5_DEL_MEM(mesh,hpar.geom);
  MMG5_DEL_MEM(mesh,mesh->htab.geom);
  MMG5_DEL_MEM(mesh,mesh->adjt);
  MMG5_DEL_MEM(mesh,mesh->edge);
  mesh->na = 0;

  if ( mesh->nprism ) MMG5_DEL_MEM(mesh,mesh->adjapr);

  return 1;
}

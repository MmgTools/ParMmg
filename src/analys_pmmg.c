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

int PMMG_update_nmgeom(PMMG_pParMesh parmesh,MMG5_pMesh mesh){
  MMG5_pTetra     pt;
  MMG5_pPoint     p0;
  MMG5_pxPoint    pxp;
  int             k,base;
  int             *adja;
  double          n[3],t[3];
  int8_t          i,j,ip,ier;

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
        pxp->nnor = 1;
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

/** compute normals at C1 vertices, for C0: tangents */
int PMMG_update_norver( PMMG_pParMesh parmesh,MMG5_pMesh mesh ) {
  MMG5_pTria     pt;
  MMG5_pPoint    ppt;
  MMG5_xPoint    *pxp;
  double         n[3],dd;
  int            *adja,k,kk,ng,nn,nt,nf,nnr;
  int            i,ii,i1;

  /* compute normals + tangents */
  nn = ng = nt = nf = 0;
  ++mesh->base;
  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    adja = &mesh->adjt[3*(k-1)+1];
    for (i=0; i<3; i++) {
      ppt = &mesh->point[pt->v[i]];
      if ( !(ppt->tag & MG_OLDPARBDY) || ppt->tag & MG_PARBDY || ppt->tag & MG_CRN || ppt->tag & MG_NOM || ppt->flag == mesh->base )  continue;

      /* C1 point */
      if ( !MG_EDG(ppt->tag) ) {

        if ( (!mesh->nc1) ||
             ppt->n[0]*ppt->n[0]+ppt->n[1]*ppt->n[1]+ppt->n[2]*ppt->n[2]<=MMG5_EPSD2 ) {
          if ( !MMG5_boulen(mesh,mesh->adjt,k,i,ppt->n) ) {
            ++nf;
            continue;
          }
          else ++nn;
        }

        ++mesh->xp;
        if(mesh->xp > mesh->xpmax){
          MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,MMG5_GAP,MMG5_xPoint,
                             "larger xpoint table",
                             mesh->xp--;return 0;);
        }
        ppt->xp = mesh->xp;
        pxp = &mesh->xpoint[ppt->xp];
        memcpy(pxp->n1,ppt->n,3*sizeof(double));
        ppt->n[0] = ppt->n[1] = ppt->n[2] = 0.;
        ppt->flag = mesh->base;

      }

      /* along ridge-curve */
      i1  = MMG5_inxt2[i];
      if ( !MG_EDG(pt->tag[i1]) )  continue;
      else if ( !MMG5_boulen(mesh,mesh->adjt,k,i,n) ) {
        ++nf;
        continue;
      }
      ++mesh->xp;
      if(mesh->xp > mesh->xpmax){
        MMG5_TAB_RECALLOC(mesh,mesh->xpoint,mesh->xpmax,MMG5_GAP,MMG5_xPoint,
                           "larger xpoint table",
                           mesh->xp--;return 0;);
      }
      ppt->xp = mesh->xp;
      pxp = &mesh->xpoint[ppt->xp];
      memcpy(pxp->n1,n,3*sizeof(double));

      if ( pt->tag[i1] & MG_GEO && adja[i1] > 0 ) {
        kk = adja[i1] / 3;
        ii = adja[i1] % 3;
        ii = MMG5_inxt2[ii];
        if ( !MMG5_boulen(mesh,mesh->adjt,kk,ii,n) ) {
          ++nf;
          continue;
        }
        memcpy(pxp->n2,n,3*sizeof(double));

        /* compute tangent as intersection of n1 + n2 */
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
        ppt->flag = mesh->base;
        ++nt;
        continue;
      }

      /* compute tgte */
      ppt->flag = mesh->base;
      ++nt;
      if ( !MMG5_boulec(mesh,mesh->adjt,k,i,ppt->n) ) {
        ++nf;
        continue;
      }
      dd = pxp->n1[0]*ppt->n[0] + pxp->n1[1]*ppt->n[1] + pxp->n1[2]*ppt->n[2];
      ppt->n[0] -= dd*pxp->n1[0];
      ppt->n[1] -= dd*pxp->n1[1];
      ppt->n[2] -= dd*pxp->n1[2];
      dd = ppt->n[0]*ppt->n[0] + ppt->n[1]*ppt->n[1] + ppt->n[2]*ppt->n[2];
      if ( dd > MMG5_EPSD2 ) {
        dd = 1.0 / sqrt(dd);
        ppt->n[0] *= dd;
        ppt->n[1] *= dd;
        ppt->n[2] *= dd;
      }
    }
  }
  mesh->nc1 = 0;

  if ( abs(mesh->info.imprim) > 3 && nn+nt > 0 ) {
    if ( nnr )
      fprintf(stdout,"     %d input normals ignored\n",nnr);
    fprintf(stdout,"     %d normals,  %d tangents updated  (%d failed)\n",nn,nt,nf);
  }
  return 1;
}

int PMMG_update_analys(PMMG_pParMesh parmesh) {
  PMMG_pGrp  grp;
  MMG5_pMesh mesh;
  MMG5_Hash  hash;
  int igrp,ip;

  for( igrp = 0; igrp < parmesh->ngrp; igrp++ ) {
    grp = &parmesh->listgrp[igrp];
    mesh = grp->mesh;

    for( ip = 1; ip <= mesh->np; ip++ )
      mesh->point[ip].flag = mesh->base;

    /* rebuild triangles*/
    if ( mesh->tria )
      MMG5_DEL_MEM(mesh,mesh->tria);
    if ( mesh->adjt )
      MMG5_DEL_MEM(mesh,mesh->adjt);
    mesh->nt = 0;

    if ( !MMG5_chkBdryTria(mesh) ) {
      fprintf(stderr,"\n  ## Error: %s: unable to rebuild triangles\n",__func__);
      return -1;
    }

    /* create tetra adjacency */
    if ( !MMG3D_hashTetra(mesh,0) ) {
      fprintf(stderr,"\n  ## Hashing problem (1). Exit program.\n");
      return 0;
    }

    /* create surface adjacency */
    if ( !MMG3D_hashTria(mesh,&hash) ) {
      MMG5_DEL_MEM(mesh,hash.item);
      fprintf(stderr,"\n  ## Hashing problem (2). Exit program.\n");
      return 0;
    }

    if( !PMMG_update_norver(parmesh,mesh) ) return 0;

    if( !PMMG_update_singul(parmesh,mesh) ) return 0;


    /* define geometry for non manifold points */
    if( !PMMG_update_nmgeom(parmesh,mesh) ) return 0;

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

int PMMG_setdhd(PMMG_pParMesh parmesh,MMG5_pMesh mesh,MMG5_HGeom *pHash ) {
  PMMG_pInt_comm int_face_comm,int_edge_comm;
  PMMG_pExt_comm ext_face_comm,ext_edge_comm;
  MMG5_Hash      hash;
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

  /** Hash triangles */
  if ( ! MMG5_hashNew(mesh,&hash,0.51*mesh->nt,1.51*mesh->nt) ) return 0;

  for( k = 1; k <= mesh->nt; k++) {
    ptr = &mesh->tria[k];
    if ( !MMG5_hashFace(mesh,&hash,ptr->v[0],ptr->v[1],ptr->v[2],k) ) {
      MMG5_DEL_MEM(mesh,hash.item);
      return 0;
    }
    ptr->flag = PMMG_UNSET;
  }

  /* Flag triangles with the color of the external process sharing them */
  int_face_comm = parmesh->int_face_comm;
  PMMG_CALLOC(parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,
              "intvalues",return 0);

  for( k = 0; k < parmesh->next_face_comm; k++ ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    nitem         = ext_face_comm->nitem;
    color         = ext_face_comm->color_out;
    for( i = 0; i < nitem; i++ ) {
      idx = ext_face_comm->int_comm_index[i];
      int_face_comm->intvalues[idx] = color;
    }
  }

  PMMG_pGrp grp = &parmesh->listgrp[0];
  MMG5_pTetra pt;
  int ie,ifac,ia,ib,ic;
  for( i = 0; i < grp->nitem_int_face_comm; i++ ) {
    ie   =  grp->face2int_face_comm_index1[i]/12;
    ifac = (grp->face2int_face_comm_index1[i]%12)/3;
    idx  =  grp->face2int_face_comm_index2[i];
    /* Get triangle index from hash table */
    pt = &mesh->tetra[ie];
    ia = pt->v[MMG5_idir[ifac][0]];
    ib = pt->v[MMG5_idir[ifac][1]];
    ic = pt->v[MMG5_idir[ifac][2]];
    k = MMG5_hashGetFace(&hash,ia,ib,ic);
    ptr = &mesh->tria[k];
    /* Flag triangle with the color of the external process */
    ptr->flag = int_face_comm->intvalues[idx];
  }

  PMMG_DEL_MEM(parmesh,int_face_comm->intvalues,int,"intvalues");
  MMG5_DEL_MEM(mesh,hash.item);


  /** Loop on true boundary triangles and store the normal in the edge internal
   *  communicator where the triangle touches a parallel edge. */
  for( k = 1; k <= mesh->nt; k++ ) {
    ptr = &mesh->tria[k];
    if( !MG_EOK(ptr) )  continue;

    /* Skip faces that are just parallel or analyzed by another process */
    tag = ptr->tag[0] & ptr->tag[1] & ptr->tag[2];
    if( (tag & MG_PARBDY) )
      if( !(tag & MG_BDY) )
        if( !(tag & MG_PARBDYBDY && ptr->flag > parmesh->myrank ) )
          continue;

    /* triangle normal */
    MMG5_nortri(mesh,ptr,n1);

    /* Get parallel edge touched by a MG_BDY face and store normal vectors */
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

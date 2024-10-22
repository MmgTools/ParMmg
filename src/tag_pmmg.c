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
 * \file tag_pmmg.c
 * \brief Functions related to boundary and interface tags
 * \author Luca Cirrottola (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "parmmg.h"


/**
 * \param ppt pointer to the point.
 *
 * Tag a node as parallel.
 */
void PMMG_tag_par_node(MMG5_pPoint ppt){

  ppt->tag |= (MG_PARBDY + MG_BDY);
  if( !(ppt->tag & MG_REQ) ) {
    /* do not add the MG_NOSURF tag on a required entity */
    ppt->tag |= (MG_REQ + MG_NOSURF);
  }
}

/**
 * \param pxt pointer to the xtetra.
 * \param j local index of the edge on the tetra.
 *
 * Tag an edge as parallel.
 */
void PMMG_tag_par_edge(MMG5_pxTetra pxt,int j){

  pxt->tag[j] |= (MG_PARBDY + MG_BDY);
  if( !(pxt->tag[j] & MG_REQ) ) {
    /* do not add the MG_NOSURF tag on a required entity */
    pxt->tag[j] |= (MG_REQ + MG_NOSURF);
  }
}

/**
 * \param pt pointer to the tetra.
 * \param hash edges hash table.
 * \param ia local index of the edge on the tetra.
 *
 * Tag an edge as parallel.
 */
int PMMG_tag_par_edge_hash(MMG5_pTetra pt,MMG5_HGeom hash,int ia){
  int      ip0,ip1,getref;
  uint16_t gettag;

  ip0 = pt->v[MMG5_iare[ia][0]];
  ip1 = pt->v[MMG5_iare[ia][1]];
  if( !MMG5_hTag( &hash, ip0, ip1, 0, MG_PARBDY + MG_BDY ) ) return 0;

  if( !MMG5_hGet( &hash, ip0, ip1, &getref, &gettag ) ) return 0;
  if( !(gettag & MG_REQ) ) {
    /* do not add the MG_NOSURF tag on a required entity */
    if( !MMG5_hTag( &hash, ip0, ip1, 0, MG_REQ + MG_NOSURF ) ) return 0;
  }

  return 1;
}

/**
 * \param pxt pointer to the xtetra.
 * \param j local index of the face on the tetra.
 *
 * Tag a face as parallel.
 */
void PMMG_tag_par_face(MMG5_pxTetra pxt,int j){

  pxt->ftag[j] |= (MG_PARBDY + MG_BDY );
  if( !(pxt->ftag[j] & MG_REQ) ) {
    /* do not add the MG_NOSURF tag on a required entity */
    pxt->ftag[j] |= (MG_REQ + MG_NOSURF);
  }
}

/**
 * \param ptt pointer to the tria.
 *
 * Tag the edges on a triangle as parallel.
 */
void PMMG_tag_par_tria(MMG5_pTria ptt){
  int j;

  for( j = 0; j < 3; j++ ) {
    ptt->tag[j] |= MG_PARBDY;
    if( !(ptt->tag[j] & MG_REQ) ) {
      /* do not add the MG_NOSURF tag on a required entity */
      ptt->tag[j] |= (MG_REQ + MG_NOSURF);
    }
  }
}

/**
 * \param pxt pointer to the point.
 *
 * Untag a parallel node.
 * Also delete the MG_PARBDYBDY tag, as it will be put only on nodes which are
 * visible from parallel faces.
 */
void PMMG_untag_par_node(MMG5_pPoint ppt){

  /* Only work on a parallel entity, so that a successive application of this
   * function on the same entity has no effect */
  if ( ppt->tag & MG_PARBDY ) {
    ppt->tag &= ~MG_PARBDY;
    ppt->tag &= ~MG_BDY;
    ppt->tag &= ~MG_REQ;
    /* a truly required entity has had the MG_NOSURF tag erased by the
     * analysis so reapply the MG_REQ tag on it */
    if( ppt->tag & MG_NOSURF ) {
      ppt->tag &= ~MG_NOSURF;
    } else {
      ppt->tag |= MG_REQ;
    }
    ppt->tag &= ~MG_PARBDYBDY;
  }
}

/**
 * \param pxt pointer to the xtetra.
 * \param j local index of the edge on the tetra.
 *
 * Untag a parallel edge.
 * Also delete the MG_PARBDYBDY tag, as it will be put only on edges which are
 * visible from parallel faces.
 */
void PMMG_untag_par_edge(MMG5_pxTetra pxt,int j){

  /* Only work on a parallel entity, so that a successive application of this
   * function on the same entity has no effect */
  if ( pxt->tag[j] & MG_PARBDY ) {
    pxt->tag[j] &= ~MG_PARBDY;
    pxt->tag[j] &= ~MG_BDY;
    pxt->tag[j] &= ~MG_REQ;
    /* a truly required entity has had the MG_NOSURF tag erased by the
     * analysis so reapply the MG_REQ tag on it */
    if( pxt->tag[j] & MG_NOSURF ) {
      pxt->tag[j] &= ~MG_NOSURF;
    } else {
      pxt->tag[j] |= MG_REQ;
    }
    pxt->tag[j] &= ~MG_PARBDYBDY;
  }
}

/**
 * \param pxt pointer to the xtetra.
 * \param j local index of the face on the tetra.
 *
 * Untag a parallel face.
 * Do not delete the MG_PARBDYBDY tag, as it will be used to recognize boundary
 * faces.
 */
void PMMG_untag_par_face(MMG5_pxTetra pxt,int j){

  /* Only work on a parallel entity, so that a successive application of this
   * function on the same entity has no effect */
  if ( pxt->ftag[j] & MG_PARBDY ) {
    pxt->ftag[j] &= ~MG_PARBDY;
    pxt->ftag[j] &= ~MG_BDY;
    pxt->ftag[j] &= ~MG_REQ;
    /* a truly required entity has had the MG_NOSURF tag erased by the
     * analysis so reapply the MG_REQ tag on it */
    if( pxt->ftag[j] & MG_NOSURF ) {
      pxt->ftag[j] &= ~MG_NOSURF;
    } else {
      pxt->ftag[j] |= MG_REQ;
    }
    /* do not delete the MG_PARBDYBDY tag, as it will be used to recognize
     * boundary faces */
  }
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * \return 1 if success, 0 if fail.
 *
 * Count the number of parallel triangle in each tetra and store the info into
 * the mark field of the tetra.
 * Reset the OLDPARBDY tag after mesh adaptation in order to track previous
 * parallel faces during load balancing.
 *
 */
inline int PMMG_resetOldTag(PMMG_pParMesh parmesh) {
  MMG5_pMesh   mesh;
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_pPoint  ppt;
  int          k,i,j,l,ip;

  for ( i=0; i<parmesh->ngrp; ++i ) {
    mesh = parmesh->listgrp[i].mesh;

    if ( !mesh ) continue;

    /* Untag old parallel points */
    for( ip = 1; ip <= mesh->np; ip++ ) {
      ppt = &mesh->point[ip];
      ppt->tag &= ~MG_OLDPARBDY;
    }

    for ( k=1; k<=mesh->ne; ++k ) {
      pt       = &mesh->tetra[k];
      pt->mark = 1;

      if ( (!MG_EOK(pt)) || (!pt->xt) ) continue;
      pxt = &mesh->xtetra[pt->xt];

      for ( j=0; j<4; ++j ) {
        if ( pxt->ftag[j] & MG_PARBDY ) {
          /* Increase counter */
          ++pt->mark;
          /* Mark face as a previously parallel one */
          pxt->ftag[j] |= MG_OLDPARBDY;
          /* Mark its points */
          for( l = 0; l < 3; l++ )
            mesh->point[pt->v[MMG5_idir[j][l]]].tag |= MG_OLDPARBDY;
          /* Check that there is no reference on an old parallel face that is
           * not a true boundary */
          if( !(pxt->ftag[j] & MG_PARBDYBDY) ) assert( !pxt->ref[j] );
        } else if ( pxt->ftag[j] & MG_OLDPARBDY ) {
          /* Untag faces which are not parallel anymore */
          pxt->ftag[j] &= ~MG_OLDPARBDY;
        }
      }
    }
  }

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Update the parallel-related tags on the points and tetra: first, remove
 * obsolete parallel markers; second, re-tag boundary entites; third, tag new
 * parallel interfaces from the internal communicator.
 *
 */
int PMMG_updateTag(PMMG_pParMesh parmesh) {
  PMMG_pGrp       grp;
  MMG5_pMesh      mesh;
  MMG5_pTetra     pt;
  MMG5_pxTetra    pxt;
  MMG5_pPoint     ppt;
  MMG5_HGeom      hash;
  int             *node2int_node_comm0_index1,*face2int_face_comm0_index1;
  int             grpid,iel,ifac,ia,ip0,ip1,k,j,i,getref;
  uint16_t        gettag;
  int8_t          isbdy;

  /* Loop on groups */
  for ( grpid=0; grpid<parmesh->ngrp; grpid++ ) {
    grp                        = &parmesh->listgrp[grpid];
    mesh                       = grp->mesh;
    node2int_node_comm0_index1 = grp->node2int_node_comm_index1;
    face2int_face_comm0_index1 = grp->face2int_face_comm_index1;

    /** Step 1: Loop on xtetras to untag old parallel entities, then build
     * hash table for edges on xtetras. */
    for ( k=1; k<=mesh->ne; k++ ) {
      pt = &mesh->tetra[k];
      if ( !pt->xt ) continue;
      pxt = &mesh->xtetra[pt->xt];
      /* Untag parallel nodes: remove PARBDY, BDY, REQ, PARBDYBDY and NOSURF
       * tag. Point not marked by NOSURF tag are required by the user: re-add
       * the REQ tag. */
      for ( j=0 ; j<4 ; j++ ) {
        ppt = &mesh->point[pt->v[j]];
        PMMG_untag_par_node(ppt);
      }
      /* Untag parallel edges: remove PARBDY, BDY, REQ, PARBDYBDY and NOSURF
       * tag. Point not marked by NOSURF tag are required by the user: re-add
       * the REQ tag. */
      for ( j=0 ; j<6 ; j++ )
        PMMG_untag_par_edge(pxt,j);
      /* Untag parallel faces: remove PARBDY, BDY, REQ and NOSURF tags but nor
       * PARBDYBDY one (used to recognize BDY faces whose BDY tag has been
       * removed). Point not marked by NOSURF tag are required by the user:
       * re-add the REQ tag.  */
      for ( j=0 ; j<4 ; j++ )
        PMMG_untag_par_face(pxt,j);
    }

    /* Create hash table for edges on xtetra, without reference or tags */
    if ( !MMG5_hNew(mesh, &hash, 6*mesh->xt, 8*mesh->xt) ) return 0;
    for ( k=1; k<=mesh->ne; k++ ) {
      pt = &mesh->tetra[k];
      if ( !pt->xt ) continue;
      for ( j=0; j<6; j++ ) {
        ip0 = pt->v[MMG5_iare[j][0]];
        ip1 = pt->v[MMG5_iare[j][1]];
        if( !MMG5_hEdge( mesh, &hash, ip0, ip1, 0, MG_NOTAG ) ) return 0;
      }
    }

    /** Step 2: Re-tag boundary entities starting from xtetra faces.
     *  Just add the appropriate tag to the face, edge (hash) or node. */
    for ( k=1; k<=mesh->ne; k++ ) {
      pt = &mesh->tetra[k];
      if ( !pt->xt ) continue;
      pxt = &mesh->xtetra[pt->xt];
      /* Look for external boundary faces (MG_BDY) or internal boundary faces
       * previously on parallel interfaces (MG_PARBDYBDY), tag their edges and
       * nodes (the BDY tag could have been removed when deleting old parallel
       * interfaces in step 1).*/
      for ( ifac=0 ; ifac<4 ; ifac++ ) {
        if ( pxt->ftag[ifac] & MG_PARBDYBDY ) {
          pxt->ftag[ifac] &= ~MG_PARBDYBDY;
          for( j = 0; j < 3; j++)
            mesh->point[pt->v[MMG5_idir[ifac][j]]].tag &= ~MG_PARBDYBDY;
          pxt->ftag[ifac] |= MG_BDY;
        }
        /* Only a "true" boundary after this line */
        if ( pxt->ftag[ifac] & MG_BDY ) {
          /* Constrain boundary if -nosurf option */
          if( mesh->info.nosurf ) {
            if( !(pxt->ftag[ifac] & MG_REQ) ) {
              /* do not add the MG_NOSURF tag on a required entity */
              pxt->ftag[ifac] |= MG_REQ + MG_NOSURF;
            }
          }
          /* Tag face edges */
          for ( j=0; j<3; j++ ) {
            ia = MMG5_iarf[ifac][j];
            ip0 = pt->v[MMG5_iare[ia][0]];
            ip1 = pt->v[MMG5_iare[ia][1]];
            if( !MMG5_hTag( &hash, ip0, ip1, 0, MG_BDY ) ) return 0;
          }
          /* Tag face nodes */
          for ( j=0 ; j<3 ; j++) {
            ppt = &mesh->point[pt->v[MMG5_idir[ifac][j]]];
            ppt->tag |= MG_BDY;
            /* Constrain boundary if -nosurf option */
            if( mesh->info.nosurf ) {
              if( !(ppt->tag & MG_REQ) ) {
                /* do not add the MG_NOSURF tag on a required entity */
                ppt->tag |= MG_REQ + MG_NOSURF;
              }
            }
          }
        }
      }
    }

    /** Step 3: Tag new parallel interface entities starting from int_face_comm.
     *
     *  This step needs to be done even if the external face communicator is
     *  not allocated (for example, when a centralized mesh is loaded and split
     *  on a single proc: in this case the internal communicator is allocated,
     *  but not the external one).
     */
    for ( i=0; i<grp->nitem_int_face_comm; i++ ) {
      iel  =   face2int_face_comm0_index1[i] / 12;
      ifac = ( face2int_face_comm0_index1[i] % 12 ) / 3;
      pt = &mesh->tetra[iel];
      assert( pt->xt );
      pxt = &mesh->xtetra[pt->xt];
      /* If already boundary, make it recognizable as a "true" boundary */
      if( pxt->ftag[ifac] & MG_BDY ) {
        pxt->ftag[ifac] |= MG_PARBDYBDY;
        for( j = 0; j < 3; j++)
          mesh->point[pt->v[MMG5_idir[ifac][j]]].tag |= MG_PARBDYBDY;
      }
      /* Tag face */
      PMMG_tag_par_face(pxt,ifac);
      /* Tag face edges */
      for ( j=0; j<3; j++ ) {
        ia = MMG5_iarf[ifac][j];
        if( !PMMG_tag_par_edge_hash(pt,hash,ia) ) return 0;
      }
      /* Tag face nodes */
      for ( j=0 ; j<3 ; j++) {
        ppt = &mesh->point[pt->v[MMG5_idir[ifac][j]]];
        PMMG_tag_par_node(ppt);
      }
    }

    /** Step 4: Get edge tag from hash table, add it to the edge tag on the
     *  xtetra,and delete hash table */
    for ( k=1; k<=mesh->ne; k++ ) {
      pt = &mesh->tetra[k];
      if ( !pt->xt ) continue;
      pxt = &mesh->xtetra[pt->xt];
      /* Unreference xtetra that are not on the boundary anymore */
      isbdy = 0;
      for( ifac = 0; ifac < 4; ifac++ ) {
        if( pxt->ftag[ifac] & MG_BDY ) {
          isbdy = 1;
          break;
        }
      }
      if( !isbdy ) {
        pt->xt = 0;
        continue;
      }
      /* update edge tags */
      for ( j=0; j<6; j++ ) {
        ip0 = pt->v[MMG5_iare[j][0]];
        ip1 = pt->v[MMG5_iare[j][1]];
        /* get the tag stored in the hash table to the xtetra edge */
        if( !MMG5_hGet( &hash, ip0, ip1, &getref, &gettag ) ) return 0;
        /* the hash table should only contain boundary/parallel related tags,
         * so remove the MG_NOSURF tag if the edge is truly required */
        if( pxt->tag[j] & MG_REQ )
          gettag &= ~MG_NOSURF;

        /* set edge tag (without NOSURF tag if the edge is required by the
         * user): here we preserve the initial MG_REQ tag of each tetra, thus,
         * potential inconsistencies will not be solved.
         *
         * A xtetra may have an edge that is boundary but doesn't belong to any
         * boundary face:
         *  - if this edge is marked as MG_BDY, the edge tag should be
         *    consistent with edge tag stored from a boundary face and we have
         *    to maintain this consistency;
         *
         * - if this edge is not marked as MG_BDY (tag == 0), we are not able to
         *    know if the edge is ref or required or if it has any other tag so
         *    we are not able to maintain the tag consistency and we have to
         *    preserve the fact that the edge is not MG_BDY.
         *
         */
        if ( (pxt->tag[j] & MG_BDY) ||
             ( (pxt->ftag[MMG5_ifar[j][0]] & MG_BDY) || (pxt->ftag[MMG5_ifar[j][1]] & MG_BDY) ) ) {
               pxt->tag[j] |= gettag;
        }
      }
    }
    PMMG_DEL_MEM( mesh, hash.geom, MMG5_hgeom, "Edge hash table" );

    /** Step 5: Unreference xpoints not on BDY (or PARBDY) */
    for ( i=1; i<=mesh->np; i++ ) {
      ppt = &mesh->point[i];
      if( ppt->tag & MG_BDY ) continue;
      if( ppt->xp ) ppt->xp = 0;
    }

  }

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param mesh pointer to the mesh structure
 *
 * Update the nodes tag with MG_REF if the edge in tetra is also MG_REF
 * \remark When we enter this function, the tag MG_OLDPBDY is correctly assigned to xtetra.ftag only.
 * However, we can have only an edge on the // interface (and no face). Therefore, we need 
 * to loop over all the edges of all the tetra.
 *
 */
void PMMG_updateTagRef_node(PMMG_pParMesh parmesh, MMG5_pMesh mesh) {
  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_int      k,ip0,ip1;
  int           ia;

  for ( k=1; k<=mesh->ne; k++ ) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    if ( !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];
    for ( ia=0 ; ia<6 ; ia++ ) {
      if ( pxt->tag[ia] & MG_REF) {
        ip0 = pt->v[MMG5_iare[ia][0]];
        ip1 = pt->v[MMG5_iare[ia][1]];
        if ( !(mesh->point[ip0].tag & MG_REF) ) {
          mesh->point[ip0].tag |= MG_REF;
        }
        if ( !(mesh->point[ip1].tag & MG_REF) ) {
          mesh->point[ip1].tag |= MG_REF;
        }
      }
    }
  }
}

/**
 * \param parmesh pointer to parmesh structure.
 * \return 0 if fail, 1 if success.
 *
 * Check if faces on a parallel communicator connect elements with different
 * references, and tag them as a "true" boundary (thus PARBDYBDY).
 *
 * \remark: Edge tags are not updated along faces with the PARBDYBDY tag as
 * it is not sufficient to maintain the consistency of PARBDYBDY tags through the mesh.
 * Morover, even if we manage to have consistent tags inside one mesh, we will
 * still have to synchronize the edge tags through the partition interfaces. In consequence,
 * the PARBDYBDY tags may be not consistent throught the entire remeshing process.
 */
int PMMG_parbdySet( PMMG_pParMesh parmesh ) {
  PMMG_pGrp      grp;
  PMMG_pExt_comm ext_face_comm;
  PMMG_pInt_comm int_face_comm;
  MMG5_pMesh     mesh;
  MMG5_pTetra    pt;
  MMG5_pxTetra   pxt;
  MMG5_pPoint    ppt;
  MPI_Status     status;
  int            *face2int_face_comm_index1,*face2int_face_comm_index2;
  int            *seenFace,*intvalues,*itosend,*itorecv;
  int            ngrp,myrank,color,nitem,k,igrp,i,idx,ie,ifac;

  grp    = parmesh->listgrp;
  myrank = parmesh->myrank;
  ngrp   = parmesh->ngrp;

  /* intvalues will be used to store tetra ref */
  int_face_comm = parmesh->int_face_comm;
  PMMG_MALLOC(parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,
              "intvalues",return 0);
  intvalues = parmesh->int_face_comm->intvalues;

  /* seenFace will be used to recognize already visited faces */
  PMMG_CALLOC(parmesh,seenFace,int_face_comm->nitem,int,"seenFace",return 0);

  /** Fill the internal communicator with the first ref found */
  for( igrp = 0; igrp < ngrp; igrp++ ) {
    grp                       = &parmesh->listgrp[igrp];
    mesh                      = grp->mesh;
    face2int_face_comm_index1 = grp->face2int_face_comm_index1;
    face2int_face_comm_index2 = grp->face2int_face_comm_index2;

    /* Remove PARBDYBDY tag from faces that are not parallel (it could have
     * been put by MMG5_bdrySet if all edges on face have been tagged as
     * MG_PARBDYBDY by MMG5_mmgHashTria) */
    for( ie = 1; ie <= mesh->ne; ie++ ) {
      pt = &mesh->tetra[ie];
      if( !MG_EOK(pt) || !pt->xt ) continue;
      pxt = &mesh->xtetra[pt->xt];
      for( ifac = 0; ifac < 4; ifac++ )
        if( pxt->ftag[ifac] & MG_PARBDYBDY && !(pxt->ftag[ifac] & MG_PARBDY) )
          pxt->ftag[ifac] &= ~MG_PARBDYBDY;
    }

    for ( k=0; k<grp->nitem_int_face_comm; ++k ) {
      ie   =  face2int_face_comm_index1[k]/12;
      ifac = (face2int_face_comm_index1[k]%12)/3;
      idx  =  face2int_face_comm_index2[k];
      pt = &mesh->tetra[ie];
      assert( MG_EOK(pt) && pt->xt );
      pxt = &mesh->xtetra[pt->xt];

      /* Tag face as "true" boundary if its second ref is different or if
       * triangle has a non-nul ref in opnbdy mode */
      if( !seenFace[idx] ) {
        intvalues[idx] = pt->ref;
      }
      else if ( (mesh->info.opnbdy && pxt->ref[ifac]>0) || (intvalues[idx] != pt->ref ) ) {
        pxt->ftag[ifac] |= MG_PARBDYBDY;
        for( i = 0; i < 3; i++)
          mesh->point[pt->v[MMG5_idir[ifac][i]]].tag |= MG_PARBDYBDY;
      }

      /* Mark face each time that it's seen */
      seenFace[idx]++;
    }
  }

  /** Send and receive external communicators filled with the tetra ref */
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    nitem         = ext_face_comm->nitem;
    color         = ext_face_comm->color_out;

    PMMG_CALLOC(parmesh,ext_face_comm->itosend,nitem,int,"itosend array",
                return 0);
    itosend = ext_face_comm->itosend;

    PMMG_CALLOC(parmesh,ext_face_comm->itorecv,nitem,int,"itorecv array",
                return 0);
    itorecv = ext_face_comm->itorecv;

    for ( i=0; i<nitem; ++i ) {
      idx            = ext_face_comm->int_comm_index[i];
      itosend[i]     = intvalues[idx];
    }

    // As this function is called only during analysis we can use directly
    // read_comm communicator
    MPI_CHECK(
      MPI_Sendrecv(itosend,nitem,MPI_INT,color,MPI_COMMUNICATORS_REF_TAG,
                   itorecv,nitem,MPI_INT,color,MPI_COMMUNICATORS_REF_TAG,
                   parmesh->info.read_comm,&status),return 0 );

    /* Store the info in intvalues */
    for ( i=0; i<nitem; ++i ) {
      idx            = ext_face_comm->int_comm_index[i];
      intvalues[idx] = itorecv[i];
    }
  }

  /* Check the internal communicator */
  for( igrp = 0; igrp < ngrp; igrp++ ) {
    grp                       = &parmesh->listgrp[igrp];
    mesh                      = grp->mesh;
    face2int_face_comm_index1 = grp->face2int_face_comm_index1;
    face2int_face_comm_index2 = grp->face2int_face_comm_index2;

    for ( k=0; k<grp->nitem_int_face_comm; ++k ) {
      ie   =  face2int_face_comm_index1[k]/12;
      ifac = (face2int_face_comm_index1[k]%12)/3;
      idx  =  face2int_face_comm_index2[k];
      pt = &mesh->tetra[ie];
      assert( MG_EOK(pt) && pt->xt );
      pxt = &mesh->xtetra[pt->xt];

      /* Faces on the external communicator have been visited only once */
      if( seenFace[idx] != 1 ) continue;

      /* Tag face as "true" boundary if its ref is different (or if triangle has
       * a non nul ref in opnbdy mode), delete reference if it is only a
       * parallel boundary. */
      if ( (mesh->info.opnbdy && pxt->ref[ifac]>0) ) {
        pxt->ftag[ifac] |= MG_PARBDYBDY;
        for( i = 0; i < 3; i++)
          mesh->point[pt->v[MMG5_idir[ifac][i]]].tag |= MG_PARBDYBDY;
      } else if ( intvalues[idx] > pt->ref ) {
        /* Tria belongs to my neighbor */
        pxt->ftag[ifac] |= MG_PARBDYBDY;
        for( i = 0; i < 3; i++)
          mesh->point[pt->v[MMG5_idir[ifac][i]]].tag |= MG_PARBDYBDY;
        MG_CLR(pxt->ori,ifac);
      }
      else if ( intvalues[idx] < pt->ref ) {
        /* Tria belongs to me */
        pxt->ftag[ifac] |= MG_PARBDYBDY;
        for( i = 0; i < 3; i++)
          mesh->point[pt->v[MMG5_idir[ifac][i]]].tag |= MG_PARBDYBDY;
        MG_SET(pxt->ori,ifac);
      }
      else {
        pxt->ref[ifac] = PMMG_NUL;
      }
    }
  }

  /* Tag parallel points touched by simple MG_BDY faces as MG_PARBDYBDY
   * (a parallel surface can pinch a regular surface in just one point).
   * The same problem on edges is handled by MMG5_mmgHashTria. */
  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    if( !MG_EOK(pt) || !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];
    for( ifac = 0; ifac < 4; ifac++ ) {
      /* Loop on simple boundary faces */
      if( !(pxt->ftag[ifac] & MG_BDY) || (pxt->ftag[ifac] & MG_PARBDY ) ) continue;
      for( i = 0; i < 3; i++ ) {
        ppt = &mesh->point[pt->v[MMG5_idir[ifac][i]]];
        if( ppt->tag & MG_PARBDY )
          ppt->tag |= MG_PARBDYBDY;
      }
    }
  }


  /* Deallocate and return */
  PMMG_DEL_MEM(parmesh,int_face_comm->intvalues,int,"intvalues");
  PMMG_DEL_MEM(parmesh,seenFace,int,"seenFace");
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    PMMG_DEL_MEM(parmesh,ext_face_comm->itosend,int,"itosend array");
    PMMG_DEL_MEM(parmesh,ext_face_comm->itorecv,int,"itorecv array");
  }

  return 1;
}

/**
 * \param parmesh pointer to parmesh structure.
 * \return 0 if fail, 1 if success.
 *
 * Check if faces on a parallel communicator connect elements with different
 * references, and tag them as a "true" boundary (thus PARBDYBDY).
 */
int PMMG_parbdyTria( PMMG_pParMesh parmesh ) {
  MMG5_Hash      hash;
  PMMG_pGrp      grp = &parmesh->listgrp[0];
  PMMG_pExt_comm ext_face_comm;
  PMMG_pInt_comm int_face_comm;
  MMG5_pMesh     mesh;
  MMG5_pTetra    pt;
  MMG5_pxTetra   pxt;
  MMG5_pTria     ptt;
  MPI_Comm       comm;
  MPI_Status     status;
  int            *face2int_face_comm_index1,*face2int_face_comm_index2;
  int            *seenFace,*intvalues,*itosend,*itorecv;
  int            myrank,color,nitem,k,i,idx,ie,ifac,kt,j,ia,ib,ic;

  comm   = parmesh->comm;
  myrank = parmesh->myrank;
  assert( parmesh->ngrp == 1 );
  mesh = grp->mesh;

  /* intvalues will be used to store tetra ref */
  int_face_comm = parmesh->int_face_comm;
  PMMG_MALLOC(parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,
              "intvalues",return 0);
  intvalues = parmesh->int_face_comm->intvalues;

  /* seenFace will be used to recognize already visited faces */
  PMMG_CALLOC(parmesh,seenFace,int_face_comm->nitem,int,"seenFace",return 0);

  /* Hash triangles */
  if ( ! MMG5_hashNew(mesh,&hash,0.51*mesh->nt,1.51*mesh->nt) ) return 0;

  for (kt=1; kt<=mesh->nt; kt++) {
    ptt = &mesh->tria[kt];
    if ( !MMG5_hashFace(mesh,&hash,ptt->v[0],ptt->v[1],ptt->v[2],kt) ) {
      MMG5_DEL_MEM(mesh,hash.item);
      return 0;
    }
  }

  /* Tag face as "true" boundary if triangle has a non-nul ref in opnbdy mode */
  for( kt = 1; kt <= mesh->nt; kt++ ) {
    ptt = &mesh->tria[kt];
    if( mesh->info.opnbdy && ptt->ref > 0 ) {
      for( j = 0; j < 3; j++)
        ptt->tag[j] |= MG_PARBDYBDY;
      for( j = 0; j < 3; j++)
        mesh->point[ptt->v[j]].tag |= MG_PARBDYBDY;
    }
  }

  /** Fill the internal communicator with the first ref found */
  face2int_face_comm_index1 = grp->face2int_face_comm_index1;
  face2int_face_comm_index2 = grp->face2int_face_comm_index2;

  for ( k=0; k<grp->nitem_int_face_comm; ++k ) {
    ie   =  face2int_face_comm_index1[k]/12;
    ifac = (face2int_face_comm_index1[k]%12)/3;
    idx  =  face2int_face_comm_index2[k];
    pt = &mesh->tetra[ie];
    assert( MG_EOK(pt) );

    if( !seenFace[idx] ) {
      intvalues[idx] = pt->ref;
    }

    /* Mark face each time that it's seen */
    seenFace[idx]++;
  }

  /** Send and receive external communicators filled with the tetra ref */
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    nitem         = ext_face_comm->nitem;
    color         = ext_face_comm->color_out;

    PMMG_CALLOC(parmesh,ext_face_comm->itosend,nitem,int,"itosend array",
                return 0);
    itosend = ext_face_comm->itosend;

    PMMG_CALLOC(parmesh,ext_face_comm->itorecv,nitem,int,"itorecv array",
                return 0);
    itorecv = ext_face_comm->itorecv;

    for ( i=0; i<nitem; ++i ) {
      idx            = ext_face_comm->int_comm_index[i];
      itosend[i]     = intvalues[idx];
    }

    MPI_CHECK(
      MPI_Sendrecv(itosend,nitem,MPI_INT,color,MPI_COMMUNICATORS_REF_TAG,
                   itorecv,nitem,MPI_INT,color,MPI_COMMUNICATORS_REF_TAG,
                   comm,&status),return 0 );

    /* Store the info in intvalues */
    for ( i=0; i<nitem; ++i ) {
      idx            = ext_face_comm->int_comm_index[i];
      intvalues[idx] = itorecv[i];
    }
  }

  /* Check the internal communicator */
  face2int_face_comm_index1 = grp->face2int_face_comm_index1;
  face2int_face_comm_index2 = grp->face2int_face_comm_index2;

  for ( k=0; k<grp->nitem_int_face_comm; ++k ) {
    ie   =  face2int_face_comm_index1[k]/12;
    ifac = (face2int_face_comm_index1[k]%12)/3;
    idx  =  face2int_face_comm_index2[k];
    pt = &mesh->tetra[ie];
    assert( MG_EOK(pt) );

    /* Get triangle index from hash table */
    ia = pt->v[MMG5_idir[ifac][0]];
    ib = pt->v[MMG5_idir[ifac][1]];
    ic = pt->v[MMG5_idir[ifac][2]];
    kt = MMG5_hashGetFace(&hash,ia,ib,ic);
    ptt = &mesh->tria[kt];

    /* Faces on the external communicator have been visited only once */
    if( seenFace[idx] != 1 ) continue;

    /* Tag face as "true" boundary if its ref is different */
    if( intvalues[idx] != pt->ref ) {
      for( j = 0; j < 3; j++)
        ptt->tag[j] |= MG_PARBDYBDY;
      for( j = 0; j < 3; j++)
        mesh->point[ptt->v[j]].tag |= MG_PARBDYBDY;
      /* check orientation: set orientation of triangle outward w.r.t. the
       * highest tetra reference, remove MG_PARBDYBDY from the halo triangle */
#warning Luca: no opnbdy yet
      if( pt->ref > intvalues[idx] ) {
        ptt->v[0] = ia;
        ptt->v[1] = ib;
        ptt->v[2] = ic;
      } else {
        /* The boundary will belong to one partition only, for the other
         * partition, the triangle will be considered as simply parallel.  This
         * means that the PARBDYBDY triangle tags may be not consistent from
         * here.  */
        ptt->tag[0] &= ~MG_PARBDYBDY;
        ptt->tag[1] &= ~MG_PARBDYBDY;
        ptt->tag[2] &= ~MG_PARBDYBDY;
      }
    }
  }


  /* Deallocate and return */
  MMG5_DEL_MEM(mesh,hash.item);
  PMMG_DEL_MEM(parmesh,int_face_comm->intvalues,int,"intvalues");
  PMMG_DEL_MEM(parmesh,seenFace,int,"seenFace");
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    PMMG_DEL_MEM(parmesh,ext_face_comm->itosend,int,"itosend array");
    PMMG_DEL_MEM(parmesh,ext_face_comm->itorecv,int,"itorecv array");
  }

  return 1;
}

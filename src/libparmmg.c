/**
 * \file libparmmg.c
 * \brief Wrapper for the parallel remeshing library.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (InriaSoft)
 * \author Nikos Pattakos (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Mmain library functions (parallel remeshing starting from centralized or
 * distributed data.
 *
 */

#include "parmmg.h"


/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Check the validity of the input mesh data (tetra orientation, solution
 * compatibility with respect to the provided mesh, Mmg options).
 *
 */
int PMMG_check_inputData(PMMG_pParMesh parmesh)
{
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  int        k;

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n  -- PMMG: CHECK INPUT DATA\n");

  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    met  = parmesh->listgrp[k].met;

    /* Check options */
    if ( mesh->info.lag > -1 ) {
      fprintf(stderr,
              "  ## Error: lagrangian mode unavailable (MMG3D_IPARAM_lag):\n");
      return 0;
    } else if ( mesh->info.iso ) {
      fprintf(stderr,"  ## Error: level-set discretisation unavailable"
              " (MMG3D_IPARAM_iso):\n");
      return 0;
    } else if ( mesh->info.optimLES && met->size==6 ) {
      fprintf(stdout,"  ## Error: strong mesh optimization for LES methods"
              " unavailable (MMG3D_IPARAM_optimLES) with an anisotropic metric.\n");
      return 0;
    }
    /* specific meshing */
    if ( met->np ) {
      if ( mesh->info.optim ) {
        printf("\n  ## ERROR: MISMATCH OPTIONS: OPTIM OPTION CAN NOT BE USED"
               " WITH AN INPUT METRIC.\n");
        return 0;
      }

      if ( mesh->info.hsiz>0. ) {
        printf("\n  ## ERROR: MISMATCH OPTIONS: HSIZ OPTION CAN NOT BE USED"
               " WITH AN INPUT METRIC.\n");
        return 0;
      }
    }

    if ( mesh->info.optim &&  mesh->info.hsiz>0. ) {
      printf("\n  ## ERROR: MISMATCH OPTIONS: HSIZ AND OPTIM OPTIONS CAN NOT BE USED"
             " TOGETHER.\n");
      return 0;
    }

    /* load data */
    MMG5_warnOrientation(mesh);

    if ( met->np && (met->np != mesh->np) ) {
      fprintf(stdout,"  ## WARNING: WRONG METRIC NUMBER. IGNORED\n");
      MMG5_DEL_MEM(mesh,met->m);
      met->np = 0;
    } else if ( met->size!=1 && met->size!=6 ) {
      fprintf(stderr,"  ## ERROR: WRONG DATA TYPE.\n");
      return 0;
    }
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \return 0 if failed, 1 if success.
 *
 * Set the triangles references to the tetrahedra faces and edges, and build
 * face communicators.
 *
 * \remark Modeled after the MMG3D_bdrySet function, with additional
 * communicators construction.
 * For the face communicator:
 * - The index of the interface triangles are initially stored in
 *   ext_comm->int_comm_index
 * - Set communicator indices:
 *   - ext communicators are concatenated to index the unique int communicator
 *   - group index2 coincides with the int communicator index for ngrp==1
 *   - index1 temporarily set to the interface triangle index
 * - Invert index1 to get the comm index from the triangle index.
 * - Replace triangle index with tetra-face-node triplet (global node index
 *   is used to sort the face nodes).
 *
 * For the node communicator:
 * - Local node ids are initially stored in ext_comm->itosend
 * - Global node ids are initially stored in ext_comm->itorecv
 * - Store global node ids in point flags (use it to match corresponding
 *   triangles on the two sides of a parallel interface)
 * - The node communicator will be built after this function.
 *
 */
int PMMG_bdrySet_buildComm(PMMG_pParMesh parmesh,MMG5_pMesh mesh) {
  PMMG_pGrp      grp;
  PMMG_pExt_comm ext_face_comm,ext_node_comm;
  PMMG_pInt_comm int_face_comm,int_node_comm;
  MMG5_pTetra   pt,pt1;
  MMG5_pPrism   pp;
  MMG5_pTria    ptt;
  MMG5_pxTetra  pxt;
  MMG5_pxPrism  pxp;
  MMG5_Hash     hash;
  int      ref,*adja,adj,k,kt,ia,ib,ic,j,na,initedg[3];
  int      nitem_int_face_comm,*posInIdx1,iint,iext_comm,iext,iloc,iploc,imax;
  int16_t  tag,inittag[3];
  char     i,i1,i2;

  grp = &parmesh->listgrp[0];


  if ( !mesh->nt )  return 1;

  if ( mesh->xtetra ) {
    if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Error: %s: mesh->xtetra must be freed.\n",__func__);
    }
    return 0;
  }
  if ( mesh->xprism ) {
    if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug ) {
      fprintf(stderr,"\n  ## Error: %s: mesh->xprism must be freed.\n",__func__);
    }
    return 0;
  }

  if ( ! MMG5_hashNew(mesh,&hash,0.51*mesh->nt,1.51*mesh->nt) ) return 0;
  for (k=1; k<=mesh->nt; k++) {
    ptt = &mesh->tria[k];
    if ( !MMG5_hashFace(mesh,&hash,ptt->v[0],ptt->v[1],ptt->v[2],k) ) return 0;
  }
  na = 0;

  mesh->xt     = 0;
  mesh->xtmax  = mesh->ntmax + 2*na;
  assert(mesh->xtmax);

  MMG5_ADD_MEM(mesh,(mesh->xtmax+1)*sizeof(MMG5_xTetra),"boundary tetrahedra",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(mesh->xtetra,mesh->xtmax+1,MMG5_xTetra,return 0);


  /* Step 1: Fill links among communicators, stock local triangle index */
  if( !PMMG_build_faceCommIndex(parmesh) ) return 0;

  /* Put global node enumeration in the point flag:
   * itosend contains the local bdy node index,
   * itorecv contains the global bdy node index.
   * Then, destroy external node communicator.*/
  for( i = 1; i <= mesh->np; i++ )
    mesh->point[i].flag = PMMG_UNSET;
  for( iext_comm = 0; iext_comm < parmesh->next_node_comm; iext_comm++ ) {
    ext_node_comm = &parmesh->ext_node_comm[iext_comm];
    for( iext = 0; iext < ext_node_comm->nitem_to_share; iext++ ) {
      mesh->point[ext_node_comm->itosend[iext]].flag = ext_node_comm->itorecv[iext];
    }
  }
  PMMG_parmesh_ext_comm_free( parmesh,parmesh->ext_node_comm,parmesh->next_node_comm);
  PMMG_DEL_MEM(parmesh, parmesh->ext_node_comm,PMMG_Ext_comm,"ext node comm");
  parmesh->next_node_comm = 0;
 
  /* Step 2: Invert mapping from bdy node list to triangles. */
  PMMG_CALLOC(parmesh,posInIdx1,mesh->nt+1,int,"posInIdx1",return 0);
  for( kt = 1; kt <= mesh->nt; kt++ )
    posInIdx1[kt] = PMMG_UNSET;
  for( iint = 0; iint < parmesh->int_face_comm->nitem; iint++ ) {
    posInIdx1[grp->face2int_face_comm_index1[iint]] = iint;
  }

  
  /* assign references to tetras faces */
  if ( !mesh->info.opnbdy ) {
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;
      adja = &mesh->adja[4*(k-1)+1];
      for (i=0; i<4; i++) {
        adj = adja[i] / 4;
        pt1 = &mesh->tetra[adj];
        if ( !adj || ( pt->ref != pt1->ref) ) {
          ia = pt->v[MMG5_idir[i][0]];
          ib = pt->v[MMG5_idir[i][1]];
          ic = pt->v[MMG5_idir[i][2]];
          kt = MMG5_hashGetFace(&hash,ia,ib,ic);
          assert(kt);
          if ( !pt->xt ) {
            mesh->xt++;
            if ( mesh->xt > mesh->xtmax ) {
              MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                                 "larger xtetra table",
                                 mesh->xt--;
                                 fprintf(stderr,"  Exit program.\n");return 0;);
            }
            pt->xt = mesh->xt;
          }
          ptt = &mesh->tria[kt];
          pxt = &mesh->xtetra[mesh->xt];
          pxt->ref[i]   = ptt->ref;
          pxt->ftag[i] |= MG_BDY;
          pxt->ftag[i] |= (ptt->tag[0] & ptt->tag[1] & ptt->tag[2]);
          /* Step 3: Put face in the internal communicator.
           * Triangle node matching is based on maximum node global index,
           * stored in point->flag */
          if( posInIdx1[kt] != PMMG_UNSET ) {
            pxt->ftag[i] |= MG_PARBDY + MG_REQ + MG_NOSURF;
            iploc = imax = 0;
            for( iloc = 0; iloc<3; iloc++ )
              if( mesh->point[pt->v[MMG5_idir[i][iloc]]].flag > imax ) {
                imax = mesh->point[pt->v[MMG5_idir[i][iloc]]].flag;
                iploc = iloc;
              }
            grp->face2int_face_comm_index1[posInIdx1[kt]] = 12*k+3*i+iploc;
          }
        }
      }
    }
  }
  else {
    /* Internal triangles preservations */
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;

      for (i=0; i<4; i++) {
        ia = pt->v[MMG5_idir[i][0]];
        ib = pt->v[MMG5_idir[i][1]];
        ic = pt->v[MMG5_idir[i][2]];
        kt = MMG5_hashGetFace(&hash,ia,ib,ic);

        if ( !kt ) continue;

        if ( !pt->xt ) {
          mesh->xt++;
          if ( mesh->xt > mesh->xtmax ) {
            MMG5_TAB_RECALLOC(mesh,mesh->xtetra,mesh->xtmax,0.2,MMG5_xTetra,
                               "larger xtetra table",
                               mesh->xt--;
                               fprintf(stderr,"  Exit program.\n");return 0;);
          }
          pt->xt = mesh->xt;
        }
        ptt = &mesh->tria[kt];
        pxt = &mesh->xtetra[mesh->xt];
        pxt->ref[i]   = ptt->ref;
        pxt->ftag[i] |= MG_BDY;
        pxt->ftag[i] |= (ptt->tag[0] & ptt->tag[1] & ptt->tag[2]);
        /* Step 3: Put face in the internal communicator.
         * Triangle node matching is based on maximum node global index,
         * stored in point->flag */
        if( posInIdx1[kt] != PMMG_UNSET ) {
          pxt->ftag[i] |= MG_PARBDY + MG_REQ + MG_NOSURF;
          iploc = imax = 0;
          for( iloc = 0; iloc<3; iloc++ )
            if( mesh->point[pt->v[MMG5_idir[i][iloc]]].flag > imax ) {
              imax = mesh->point[pt->v[MMG5_idir[i][iloc]]].flag;
              iploc = iloc;
            }
          grp->face2int_face_comm_index1[posInIdx1[kt]] = 12*k+3*i+iploc;
        }
      }
    }
  }

  /* Reset point flag */
  for( i = 1; i <= mesh->np; i++ )
    mesh->point[i].flag = PMMG_UNSET;

  if ( !mesh->info.opnbdy ) {
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;
      if ( !pt->xt )  continue;
      pxt = &mesh->xtetra[pt->xt];
      adja = &mesh->adja[4*(k-1)+1];
      for (i=0; i<4; i++) {
        adj = adja[i] / 4;
        pt1 = &mesh->tetra[adj];
        /* Set edge tag */
        if ( pxt->ftag[i] ) {
          if ( adj && (pt->ref == pt1->ref ) ) {
            continue;
          }
          else {
            ia = pt->v[MMG5_idir[i][0]];
            ib = pt->v[MMG5_idir[i][1]];
            ic = pt->v[MMG5_idir[i][2]];
            kt = MMG5_hashGetFace(&hash,ia,ib,ic);
            ptt = &mesh->tria[kt];

            /* Set flag to know if tetra has the same orientation than the
             * triangle: here, we can not suppose that the triangle are oriented
             * with respect to the face orientation because for non manifold
             * cases, setadj may have reoriented the triangles */
            if ( ptt->v[0] == ia && ptt->v[1] == ib && ptt->v[2] == ic ) {
              MG_SET(pxt->ori,i);
              for (j=0; j<3; j++) {
                tag = pxt->ftag[i] | ptt->tag[j];
                if ( tag ) {
                  if ( !MMG5_settag(mesh,k,MMG5_iarf[i][j],tag,ptt->edg[j]) )
                    return 0;
                }
              }
            }
            else
              MG_CLR(pxt->ori,i);
          }
        }
      }
    }
  }
  else {
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )  continue;
      if ( !pt->xt )  continue;
      pxt = &mesh->xtetra[pt->xt];
      adja = &mesh->adja[4*(k-1)+1];
      for (i=0; i<4; i++) {
        adj = adja[i] / 4;
        pt1 = &mesh->tetra[adj];

        ia = pt->v[MMG5_idir[i][0]];
        ib = pt->v[MMG5_idir[i][1]];
        ic = pt->v[MMG5_idir[i][2]];
        kt = MMG5_hashGetFace(&hash,ia,ib,ic);

        if ( !kt ) continue;

        ptt = &mesh->tria[kt];

        /* Set flag to know if tetra has the same orientation than the triangle
         * + forcd the triangle numbering to match the tetra face numbering */
        if ( adj ) {
          for ( j=0; j<3; ++j ) {
            i1 = MMG5_inxt2[j];
            i2 = MMG5_inxt2[i1];
            if (  ptt->v[j]==ia && ptt->v[i1]==ib && ptt->v[i2]==ic )
              break;
          }
          if ( j<3 ) {
            MG_SET(pxt->ori,i);
            if ( j!=0 ) {
              /* Triangle vertices+tag/edg reordering */
              ptt->v[0] = ia;
              ptt->v[1] = ib;
              ptt->v[2] = ic;

              inittag[0]  = ptt->tag[0];
              inittag[1]  = ptt->tag[1];
              inittag[2]  = ptt->tag[2];
              ptt->tag[0] = inittag[j];
              ptt->tag[1] = inittag[i1];
              ptt->tag[2] = inittag[i2];

              initedg[0]  = ptt->edg[0];
              initedg[1]  = ptt->edg[1];
              initedg[2]  = ptt->edg[2];
              ptt->edg[0] = initedg[j];
              ptt->edg[1] = initedg[i1];
              ptt->edg[2] = initedg[i2];
            }
          }
          else {
            MG_CLR(pxt->ori,i);
          }
        }
        else  MG_SET(pxt->ori,i);

        /* Set edge tag */
        if ( pxt->ftag[i] ) {
          if ( adj && ( (pt->ref < pt1->ref) || !MG_GET(pxt->ori,i) ) ) {
            continue;
          }
          else {
            for (j=0; j<3; j++) {
              tag = pxt->ftag[i] | ptt->tag[j];
              if ( tag ) {
                if ( !MMG5_settag(mesh,k,MMG5_iarf[i][j],tag,ptt->edg[j]) )
                  return 0;
              }
            }
          }
        }
      }
    }
  }

  if ( !mesh->nprism ) {
    MMG5_DEL_MEM(mesh,hash.item);
    return 1;
  }

  mesh->xpr     = 0;
  MMG5_ADD_MEM(mesh,(mesh->nprism+1)*sizeof(MMG5_xPrism),"boundary prisms",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(mesh->xprism,mesh->nprism+1,MMG5_xPrism,return 0);

  /* assign references to prism faces */
  for (k=1; k<=mesh->nprism; k++) {
    pp = &mesh->prism[k];
    if ( !MG_EOK(pp) )  continue;
    adja = &mesh->adjapr[5*(k-1)+1];
    for (i=0; i<2; i++) {
      adj = adja[i] / 5;
      if ( adj<0 ) {
        ref = mesh->tetra[abs(adj)].ref;
      } else {
        ref = mesh->prism[adj].ref;
      }
      if ( adj && (pp->ref == ref) ) continue;

      ia = pp->v[MMG5_idir_pr[i][0]];
      ib = pp->v[MMG5_idir_pr[i][1]];
      ic = pp->v[MMG5_idir_pr[i][2]];
      kt = MMG5_hashGetFace(&hash,ia,ib,ic);
      assert(kt);
      if ( !pp->xpr ) {
        mesh->xpr++;
        pp->xpr = mesh->xpr;
      }
      ptt = &mesh->tria[kt];
      pxp = &mesh->xprism[mesh->xpr];
      pxp->ref[i]   = ptt->ref;
      pxp->ftag[i] |= MG_BDY;
      pxp->ftag[i] |= (ptt->tag[0] & ptt->tag[1] & ptt->tag[2]);

      for (j=0; j<3; j++) {
        pxp->tag[MMG5_iarf[i][j]] |= pxp->ftag[i] | ptt->tag[j];
        pxp->edg[MMG5_iarf[i][j]] = ptt->edg[j];
      }
    }
  }
  MMG5_ADD_MEM(mesh,(mesh->xpr-mesh->nprism)*sizeof(MMG5_xPrism),"boundary prisms",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_RECALLOC(mesh->xprism,mesh->nprism+1,mesh->xpr+1,MMG5_xPrism,
                      "boundary prisms",return 0);

  MMG5_DEL_MEM(mesh,hash.item);
  PMMG_DEL_MEM(parmesh,posInIdx1,int,"posInIdx1");
  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param mesh pointer toward the mesh structure
 *
 * \remark Modeled after the MMG3D_analys function, with additional face
 *         communicators construction.
 */
int PMMG_analys_buildComm(PMMG_pParMesh parmesh,MMG5_pMesh mesh) {
  PMMG_pExt_comm ext_node_comm;
  MMG5_Hash      hash;
  int            i,iext_comm,iext;

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

  /* check for ridges */
  if ( mesh->info.dhd > MMG5_ANGLIM && !MMG5_setdhd(mesh) ) {
    fprintf(stderr,"\n  ## Geometry problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /* identify singularities */
  if ( !MMG5_singul(mesh) ) {
    fprintf(stderr,"\n  ## MMG5_Singularity problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  if ( abs(mesh->info.imprim) > 3 || mesh->info.ddebug )
    fprintf(stdout,"  ** DEFINING GEOMETRY\n");

  /* define (and regularize) normals */
  if ( !MMG5_norver(mesh) ) {
    fprintf(stderr,"\n  ## Normal problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    return 0;
  }

  /* set bdry entities to tetra */
  if ( !MMG5_bdrySet(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    MMG5_DEL_MEM(mesh,hash.item);
    MMG5_DEL_MEM(mesh,mesh->xpoint);
    return 0;
  }

  switch( parmesh->info.API_mode ) {
    case PMMG_APIDISTRIB_faces :
      /* Set face communicators indexing, convert tria index into iel face index
       * (it needs a valid cc field in each tria), and tag xtetra face as PARBDY
       * */
      if( !PMMG_build_faceCommIndex( parmesh ) ) return 0;
      PMMG_tria2elmFace_coords( parmesh );
      break;
    case PMMG_APIDISTRIB_nodes :
      if( !PMMG_build_nodeCommIndex( parmesh ) ) return 0;
      break;
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

  /* check subdomains connected by a vertex and mark these vertex as corner and required */
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

  if( parmesh->info.API_mode == PMMG_APIDISTRIB_nodes ) {
    if ( !PMMG_build_faceCommFromNodes(parmesh) ) {
      return PMMG_STRONGFAILURE;
    }
  }

  /* release memory */
  MMG5_DEL_MEM(mesh,mesh->htab.geom);
  MMG5_DEL_MEM(mesh,mesh->adjt);
  MMG5_DEL_MEM(mesh,mesh->tria);
  mesh->nt = 0;

  if ( mesh->nprism ) MMG5_DEL_MEM(mesh,mesh->adjapr);

  return 1;
}

/**
 * \param  parmesh pointer to parmesh structure
 *
 * \return PMMG_SUCCESS if success, PMMG_LOWFAILURE if fail and return an
 * unscaled mesh, PMMG_STRONGFAILURE if fail and return a scaled mesh.
 *
 * Mesh preprocessing: set function pointers, scale mesh, perform mesh
 * analysis and display length and quality histos.
 */
int PMMG_preprocessMesh( PMMG_pParMesh parmesh )
{
  MMG5_pMesh mesh;
  MMG5_pSol  met;

  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  assert ( ( mesh != NULL ) && ( met != NULL ) && "Preprocessing empty args");

  /** Function setters (must be assigned before quality computation) */
  MMG3D_Set_commonFunc();

  /** Mesh scaling and quality histogram */
  if ( !MMG5_scaleMesh(mesh,met) ) {
    return PMMG_LOWFAILURE;
  }

  /** specific meshing */
  if ( mesh->info.optim && !met->np ) {
    if ( !MMG3D_doSol(mesh,met) ) {
      return PMMG_STRONGFAILURE;
    }
    MMG3D_solTruncatureForOptim(mesh,met);
  }

  if ( mesh->info.hsiz > 0. ) {
    if ( !MMG3D_Set_constantSize(mesh,met) ) {
      return PMMG_STRONGFAILURE;
    }
  }

  MMG3D_setfunc(mesh,met);

  if ( !MMG3D_tetraQual( mesh, met, 0 ) ) {
    return PMMG_STRONGFAILURE;
  }

  if ( !PMMG_qualhisto(parmesh,PMMG_INQUA) ) {
    return PMMG_STRONGFAILURE;
  }

  /** Mesh analysis */
  if ( !MMG3D_analys(mesh) ) {
    return PMMG_STRONGFAILURE;
  }

  if ( parmesh->info.imprim > PMMG_VERB_ITWAVES && (!mesh->info.iso) && met->m ) {
    MMG3D_prilen(mesh,met,0);
  }

  /** Mesh unscaling */
  if ( !MMG5_unscaleMesh(mesh,met) ) {
    return PMMG_STRONGFAILURE;
  }

  return PMMG_SUCCESS;
}

/**
 * \param  parmesh pointer to parmesh structure
 *
 * \return PMMG_SUCCESS if success, PMMG_LOWFAILURE if fail and return an
 * unscaled mesh, PMMG_STRONGFAILURE if fail and return a scaled mesh.
 *
 * Mesh preprocessing (distributed input mesh): set function pointers,
 * scale mesh, perform mesh analysis, display length and quality histos,
 * and build communicators.
 */
int PMMG_preprocessMesh_distributed( PMMG_pParMesh parmesh )
{
  MMG5_pMesh mesh;
  MMG5_pSol  met;

  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  assert ( ( mesh != NULL ) && ( met != NULL ) && "Preprocessing empty args");

  /* Set distributed API mode */
  if( parmesh->next_face_comm ) {
    parmesh->info.API_mode = PMMG_APIDISTRIB_faces;
  } else if( parmesh->next_node_comm ) {
    parmesh->info.API_mode = PMMG_APIDISTRIB_nodes;
  } else {
    fprintf(stderr," ## Error: %s: parallel faces or nodes must be set through the API interface\n",__func__);
    return 0;
  }

  /** Function setters (must be assigned before quality computation) */
  MMG3D_Set_commonFunc();

  /** Mesh scaling and quality histogram */
  if ( !MMG5_scaleMesh(mesh,met) ) {
    return PMMG_LOWFAILURE;
  }

  /** specific meshing */
  if ( mesh->info.optim && !met->np ) {
    if ( !MMG3D_doSol(mesh,met) ) {
      return PMMG_STRONGFAILURE;
    }
    MMG3D_solTruncatureForOptim(mesh,met);
  }

  if ( mesh->info.hsiz > 0. ) {
    if ( !MMG3D_Set_constantSize(mesh,met) ) {
      return PMMG_STRONGFAILURE;
    }
  }

  MMG3D_setfunc(mesh,met);

  if ( !MMG3D_tetraQual( mesh, met, 0 ) ) {
    return PMMG_STRONGFAILURE;
  }

  if ( !PMMG_qualhisto(parmesh,PMMG_INQUA) ) {
    return PMMG_STRONGFAILURE;
  }

  /** Mesh analysis and face communicators construction */
  if ( !PMMG_analys_buildComm(parmesh,mesh) ) {
    return PMMG_STRONGFAILURE;
  }

  if ( parmesh->info.imprim > PMMG_VERB_ITWAVES && (!mesh->info.iso) && met->m ) {
    MMG3D_prilen(mesh,met,0);
  }

  /** Mesh unscaling */
  if ( !MMG5_unscaleMesh(mesh,met) ) {
    return PMMG_STRONGFAILURE;
  }

  if( parmesh->info.API_mode == PMMG_APIDISTRIB_faces ) {
    /** Build node communicators (mesh needs to be unscaled) */
    PMMG_parmesh_ext_comm_free( parmesh,parmesh->ext_node_comm,parmesh->next_node_comm);
    PMMG_DEL_MEM(parmesh, parmesh->ext_node_comm,PMMG_Ext_comm,"ext node comm");
    parmesh->next_node_comm = 0;
    PMMG_DEL_MEM(parmesh, parmesh->ext_node_comm,PMMG_Int_comm,"int node comm");
    PMMG_CALLOC(parmesh,parmesh->int_node_comm,1,PMMG_Int_comm,"int node comm",return 0);
    if ( !PMMG_build_nodeCommFromFaces(parmesh) ) {
      return PMMG_STRONGFAILURE;
    }
  }

  return PMMG_SUCCESS;
}

int PMMG_parmmglib_centralized(PMMG_pParMesh parmesh) {
  MMG5_pMesh    mesh;
  MMG5_pSol     met;
  int           ier;
  int           iresult,ierlib;
  long int      tmpmem;
  mytime        ctim[TIMEMAX];
  int8_t        tim;
  char          stim[32];

 if ( parmesh->info.imprim > PMMG_VERB_NO ) {
    fprintf(stdout,"\n  %s\n   MODULE PARMMGLIB_CENTRALIZED: IMB-LJLL : "
            "%s (%s)\n  %s\n",PMMG_STR,PMMG_VER,PMMG_REL,PMMG_STR);
  }

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /** Check input data */
  tim = 1;
  chrono(ON,&(ctim[tim]));

  ier = PMMG_check_inputData( parmesh );
  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !iresult ) return PMMG_LOWFAILURE;

  chrono(OFF,&(ctim[tim]));
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf(stdout,"  -- CHECK INPUT DATA COMPLETED.     %s\n",stim);
  }


  chrono(ON,&(ctim[2]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS AND MESH DISTRIBUTION\n");
  }

  /** Send mesh to other procs */
  tim = 6;
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    chrono(ON,&(ctim[tim]));
    fprintf(stdout,"\n  -- BCAST" );
  }
  ier = PMMG_bcast_mesh( parmesh );
  if ( ier!=1 ) return PMMG_LOWFAILURE;

  if ( parmesh->info.imprim >= PMMG_VERB_STEPS  ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"\n  -- BCAST COMPLETED    %s\n",stim );
  }

  /** Mesh preprocessing: set function pointers, scale mesh, perform mesh
   * analysis and display length and quality histos. */
  tim = 7;
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    chrono(ON,&(ctim[tim]));
    fprintf(stdout,"\n  -- ANALYSIS" );
  }
  ier = PMMG_preprocessMesh( parmesh );
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"\n  -- ANALYSIS COMPLETED    %s\n",stim );
  }

  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;
  if ( (ier==PMMG_STRONGFAILURE) && MMG5_unscaleMesh( mesh, met ) ) {
    ier = PMMG_LOWFAILURE;
  }
  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MAX, parmesh->comm );
  if ( iresult!=PMMG_SUCCESS ) {
    return iresult;
  }

  /** Send mesh partionning to other procs */
  tim = 8;
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    chrono(ON,&(ctim[tim]));
    fprintf(stdout,"\n  -- PARTITIONING" );
  }
  if ( !PMMG_distribute_mesh( parmesh ) ) {
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_LOWFAILURE);
  }
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"\n  -- PARTITIONING COMPLETED    %s\n",stim );
  }

  chrono(OFF,&(ctim[2]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    printim(ctim[2].gdif,stim);
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);
  }

  /** Remeshing */
  tim = 3;
  chrono(ON,&(ctim[tim]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf( stdout,"\n  -- PHASE 2 : %s MESHING\n",
             met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC" );
  }

  ier = PMMG_parmmglib1(parmesh);
  MPI_Allreduce( &ier, &ierlib, 1, MPI_INT, MPI_MAX, parmesh->comm );

  chrono(OFF,&(ctim[tim]));
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
  }
  if ( ierlib == PMMG_STRONGFAILURE ) {
    return ierlib;
  }

  /** Merge all the meshes on the proc 0 */
  tim = 4;
  chrono(ON,&(ctim[tim]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf( stdout,"\n   -- PHASE 3 : MERGE MESHES OVER PROCESSORS\n" );
  }

  iresult = PMMG_merge_parmesh( parmesh );
  if ( !iresult ) {
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_STRONGFAILURE);
  }

  chrono(OFF,&(ctim[tim]));
  if ( parmesh->info.imprim >  PMMG_VERB_VERSION  ) {
    printim(ctim[tim].gdif,stim);
    fprintf( stdout,"   -- PHASE 3 COMPLETED.     %s\n",stim );
  }

  if ( !parmesh->myrank ) {
    /** Boundaries reconstruction */
    tim = 5;
    chrono(ON,&(ctim[tim]));
    if (  parmesh->info.imprim > PMMG_VERB_VERSION ) {
      fprintf( stdout,"\n   -- PHASE 4 : MESH PACKED UP\n" );
    }

    tmpmem = parmesh->memMax - parmesh->memCur;
    parmesh->memMax = parmesh->memCur;
    parmesh->listgrp[0].mesh->memMax += tmpmem;

    mesh = parmesh->listgrp[0].mesh;
    if ( (!MMG3D_hashTetra( mesh, 0 )) || (-1 == MMG3D_bdryBuild( mesh )) ) {
      /** Impossible to rebuild the triangle */
      fprintf(stdout,"\n\n\n  -- IMPOSSIBLE TO BUILD THE BOUNDARY MESH\n\n\n");
      PMMG_CLEAN_AND_RETURN(parmesh,PMMG_LOWFAILURE);
    }

    chrono(OFF,&(ctim[tim]));
    if (  parmesh->info.imprim >  PMMG_VERB_VERSION ) {
      printim(ctim[tim].gdif,stim);
      fprintf( stdout,"   -- PHASE 4 COMPLETED.     %s\n",stim );
    }
  }

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( parmesh->info.imprim >= PMMG_VERB_VERSION ) {
    fprintf(stdout,"\n   PARMMGLIB_CENTRALIZED: ELAPSED TIME  %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE PARMMGLIB_CENTRALIZED: IMB-LJLL \n  %s\n",
            PMMG_STR,PMMG_STR);
  }

  PMMG_CLEAN_AND_RETURN(parmesh,ierlib);
}

int PMMG_parmmglib_distributed(PMMG_pParMesh parmesh) {
  MMG5_pMesh       mesh;
  MMG5_pSol        met;
  int              ier,iresult,ierlib;
  long int         tmpmem;
  mytime           ctim[TIMEMAX];
  int8_t           tim;
  char             stim[32];


  if ( parmesh->info.imprim >= PMMG_VERB_VERSION ) {
    fprintf(stdout,"\n  %s\n   MODULE PARMMGLIB_DISTRIBUTED: IMB-LJLL : "
            "%s (%s)\n  %s\n",PMMG_STR,PMMG_VER,PMMG_REL,PMMG_STR);
  }

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  /** Check input data */
  tim = 1;
  chrono(ON,&(ctim[tim]));

  ier = PMMG_check_inputData( parmesh );
  MPI_CHECK( MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm ),
             return PMMG_LOWFAILURE);
  if ( !iresult ) return PMMG_LOWFAILURE;

  chrono(OFF,&(ctim[tim]));
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf(stdout,"  -- CHECK INPUT DATA COMPLETED.     %s\n",stim);
  }

  tim = 2;
  chrono(ON,&(ctim[tim]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  }

  assert ( parmesh->ngrp < 2 );
  if ( parmesh->ngrp ) {
    /** Mesh preprocessing: set function pointers, scale mesh, perform mesh
     * analysis and display length and quality histos. */
    ier  = PMMG_preprocessMesh_distributed( parmesh );
    mesh = parmesh->listgrp[0].mesh;
    met  = parmesh->listgrp[0].met;
    if ( (ier==PMMG_STRONGFAILURE) && MMG5_unscaleMesh( mesh, met ) ) {
      ier = PMMG_LOWFAILURE;
    }
  }
  else { ier = PMMG_SUCCESS; }

  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MAX, parmesh->comm );
  if ( iresult!=PMMG_SUCCESS ) {
    return iresult;
  }

  chrono(OFF,&(ctim[tim]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"   -- PHASE 1 COMPLETED.     %s\n",stim);
  }

  /** Remeshing */
  tim = 3;
  chrono(ON,&(ctim[tim]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf( stdout,"\n  -- PHASE 2 : %s MESHING\n",
             met->size < 6 ? "ISOTROPIC" : "ANISOTROPIC" );
  }

  ier = PMMG_parmmglib1(parmesh);
  MPI_Allreduce( &ier, &ierlib, 1, MPI_INT, MPI_MAX, parmesh->comm );

  chrono(OFF,&(ctim[tim]));
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
  }
  if ( ierlib == PMMG_STRONGFAILURE ) {
    return ierlib;
  }

  /** Boundaries reconstruction */
  tim = 4;
  chrono(ON,&(ctim[tim]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n   -- PHASE 3 : MESH PACKED UP\n");

  /** All the memory is devoted to the mesh **/
  tmpmem = parmesh->memMax - parmesh->memCur;
  parmesh->memMax = parmesh->memCur;
  parmesh->listgrp[0].mesh->memMax += tmpmem;

  mesh = parmesh->listgrp[0].mesh;
  if ( (!MMG3D_hashTetra( mesh, 0 )) || ( -1 == MMG3D_bdryBuild(parmesh->listgrp[0].mesh) ) ) {
    /** Impossible to rebuild the triangle **/
    fprintf(stdout,"\n\n\n  -- IMPOSSIBLE TO BUILD THE BOUNDARY MESH\n\n\n");
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_LOWFAILURE);
  }

  chrono(OFF,&(ctim[tim]));
  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"\n   -- PHASE 3 COMPLETED.     %s\n",stim);
  }

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( parmesh->info.imprim >= PMMG_VERB_VERSION ) {
    fprintf(stdout,"\n   PARMMGLIB_DISTRIBUTED: ELAPSED TIME  %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE PARMMGLIB_DISTRIBUTED: IMB-LJLL \n  %s\n",
            PMMG_STR,PMMG_STR);
  }

  PMMG_CLEAN_AND_RETURN(parmesh,ierlib);
}

/**
 * \file mergemesh.c
 * \brief Merge the mesh.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 */
#include "parmmg.h"
#include "mpitypes_pmmg.h"

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param grpI pointer toward the group in which we want to merge.
 * \param grpJ pointer toward the group that we merge.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Add the interface points of the group \a grpJ in the group \a grpI (mesh +
 * communicator).
 *
 */
int PMMG_mergeGrpJinI_interfacePoints_addGrpJ( PMMG_pParMesh parmesh,
                                               PMMG_pGrp grpI,PMMG_pGrp grpJ ) {
  MMG5_pMesh     meshI,meshJ;
  MMG5_pSol      metI,metJ;
  MMG5_pPoint    pptI,pptJ;
  MMG5_pxPoint   pxpI,pxpJ;
  int            *intvalues,nitem_int_node_comm;
  int            *node2int_node_comm_index1,*node2int_node_comm_index2;
  int            size,k,poi_id_glo,ip;
  static int     pmmg_warn = 0;

  intvalues = parmesh->int_node_comm->intvalues;

  meshI                     = grpI->mesh;
  metI                      = grpI->met;

  meshJ                     = grpJ->mesh;
  metJ                      = grpJ->met;
  nitem_int_node_comm       = grpJ->nitem_int_node_comm;
  node2int_node_comm_index1 = grpJ->node2int_node_comm_index1;
  node2int_node_comm_index2 = grpJ->node2int_node_comm_index2;

  for ( k=0; k<nitem_int_node_comm; ++k ) {
    poi_id_glo = node2int_node_comm_index2[k];
    assert(   ( 0 <= poi_id_glo )
              && ( poi_id_glo < parmesh->int_node_comm->nitem )
              && "check intvalues indices"  );

    /* Position of the point in grpJ */
    ip =  node2int_node_comm_index1[ k ];
    assert ( ip && ip <= meshJ->np );
    pptJ = &meshJ->point[ node2int_node_comm_index1[ k ] ];

    if ( (!MG_VOK(pptJ)) || pptJ->tmp ) continue;

    /* Point pptJ is not found in the merged mesh. Add it */
    if ( !intvalues[ poi_id_glo ] ) {
      ip = MMG3D_newPt(meshI,pptJ->c,pptJ->tag);
      if ( !ip ) {
        /* reallocation of point table */
        MMG3D_POINT_REALLOC(meshI,metI,ip,meshI->gap,
                            printf("  ## Error: unable to merge group points\n");
                            MMG5_INCREASE_MEM_MESSAGE();
                            metI->np = meshI->np;
                            return 0;,
                            pptJ->c,pptJ->tag);

      }
      assert( (ip <= meshI->npmax) && "run out of points" );

      pptJ->tmp = ip;
      intvalues[ poi_id_glo ] = ip;

      /* Add xpoint if needed */
      if ( pptJ->xp ) {
        pxpJ = &meshJ->xpoint[pptJ->xp];
        pptI = &meshI->point[ip];
        pxpI = &meshI->xpoint[pptI->xp];
        assert( (pptI->xp <= meshI->xpmax) && "increase xpoints" );
        assert( (pptI->xp > 0) && "negative xpoints" );
        memcpy(pxpI,pxpJ,sizeof(MMG5_xPoint));
      }
      if ( metI->m ) {
        size = metI->size;
        if ( !metJ->m ) {
          if ( !pmmg_warn ) {
            pmmg_warn = 1;
            printf("  ## Error: unable to merge metrics:"
                   " group I has a metric while group J don't.\n");
          }
        }
        else {
          assert( ((size==1||size==6) && size==metJ->size) && "met size issues" );
          memcpy(&metI->m[size*ip],&metJ->m[size*node2int_node_comm_index1[ k ]],size*sizeof(double));
        }
      }
    } else {
      /* point already exists in merged mesh. update his tmp field to point to
       * its meshI index */
      pptJ->tmp = abs(intvalues[ poi_id_glo ]);
      intvalues[ poi_id_glo ] *= -1;
    }
  }
  metI->np = meshI->np;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Merge the interface points of all the parmesh groups (and listed inside the
 * internal communicator) into the group 0.
 *
 */
int PMMG_mergeGrps_interfacePoints( PMMG_pParMesh parmesh ) {
  PMMG_pGrp      listgrp,grpI;
  int            *intvalues;
  int            poi_id_int,poi_id_glo,imsh,k;

  listgrp   = parmesh->listgrp;

  grpI      = &listgrp[0];
  intvalues = parmesh->int_node_comm->intvalues;

  /** Use the tmp field of points in meshes to remember the id in the merged mesh
   * of points that have already been added to the merged mesh or 0 if they
   * haven't been merged yet */
  for ( imsh = 0; imsh < parmesh->ngrp; ++imsh )
    for ( k = 1; k <= listgrp[imsh].mesh->np; ++k )
      listgrp[imsh].mesh->point[k].tmp = 0;

  /** Step 1: store the indices of the interface entities of meshI into the
   * internal communicators */
  for ( k = 0; k < grpI->nitem_int_node_comm; ++k ) {
    poi_id_int = grpI->node2int_node_comm_index1[k];
    poi_id_glo = grpI->node2int_node_comm_index2[k];
    assert(   ( 0 <= poi_id_glo )
           && ( poi_id_glo < parmesh->int_node_comm->nitem )
           && "check intvalues indices" );
    intvalues[ poi_id_glo ] = poi_id_int;
  }

  /** Step 2: add points referenced by the rest of the groups meshes' internal
   * communicators to the merging mesh (meshI) and in intvalues */
  for ( imsh=1; imsh<parmesh->ngrp; ++imsh ) {
    if ( !PMMG_mergeGrpJinI_interfacePoints_addGrpJ(parmesh,grpI,&listgrp[imsh]) )
      return 0;
  }

  return 1;
}

/**
 * \param parmesh pointer toward a parmmg parmesh mesh structure.
 * \param metI pointer toward the group in which we wnat to merge.
 * \param grpJ pointer toward the group that we want to merge into grpI
 *
 * \return 0 if fail, 1 otherwise
 *
 * Merge the internal nodes of the group \a grpJ into the group grpI.
 *
 */
int PMMG_mergeGrpJinI_internalPoints( PMMG_pGrp grpI, PMMG_pGrp grpJ ) {
  MMG5_pMesh     meshI,meshJ;
  MMG5_pSol      metI,metJ;
  MMG5_pPoint    pptI,pptJ;
  MMG5_pxPoint   pxpI,pxpJ;
  int            ip,ier,k;
  static int     pmmg_warn = 0;

  meshI = grpI->mesh;
  metI  = grpI->met;

  /** Loop over points and add the ones that are not already in the merged
   * mesh (meshI) */
  meshJ = grpJ->mesh;
  metJ  = grpJ->met;

  for ( k=1; k<=meshJ->np; k++ ) {
    pptJ = &meshJ->point[k];
    if ( !MG_VOK(pptJ) ) continue;
    if ( pptJ->tmp )     continue;

    ip = MMG3D_newPt(meshI,pptJ->c,pptJ->tag);
    if ( !ip ) {
      /* reallocation of point table */
      MMG3D_POINT_REALLOC(meshI,metI,ip,meshI->gap,
                          printf("  ## Error: unable to merge group points\n");
                          MMG5_INCREASE_MEM_MESSAGE();
                          metI->np = meshI->np;
                          return 0;,
                          pptJ->c,pptJ->tag);
    }
    pptJ->tmp = ip;

    /* Add xpoint if needed */
    ier = 1;
    if ( pptJ->xp ) {
      pxpJ = &meshJ->xpoint[pptJ->xp];
      pptI = &meshI->point[ip];
      pxpI = &meshI->xpoint[pptI->xp];
      memcpy(pxpI,pxpJ,sizeof(MMG5_xPoint));
    }
    if ( metI->m ) {
      if ( !metJ->m ) {
        if ( !pmmg_warn ) {
          pmmg_warn = 1;
          printf("  ## Error: unable to merge metrics:"
                 " group I has a metric while group J don't.\n");
        }
      }
      else {
        memcpy(&metI->m[metI->size*ip],&metJ->m[metJ->size*k],metJ->size*sizeof(double));
      }
    }
  }
  metI->np = meshI->np;
  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param grpI pointer toward the group in which we want to merge
 * \param grpJ pointer toward the group that we want to merge
 *
 * \return 0 if fail, 1 otherwise
 *
 * Merge the interfaces tetra of \a grpJ->mesh mesh into the \a grpI->mesh
 * mesh. An interface tetra is a tetra that has at least 1 interface triangle.
 *
 */
int PMMG_mergeGrpJinI_interfaceTetra(PMMG_pParMesh parmesh,PMMG_pGrp grpI,
                                     PMMG_pGrp grpJ) {
  MMG5_pMesh     meshI,meshJ;
  MMG5_pTetra    ptI,ptJ;
  MMG5_pxTetra   pxtI,pxtJ;
  int            nitem_int_face_comm;
  int            *intvalues;
  int            *face2int_face_comm_index1,*face2int_face_comm_index2;
  int            face_id_glo;
  int            k,iel,ie,ier,ifac,iploc,i,newsize;

  intvalues = parmesh->int_face_comm->intvalues;

  meshI     = grpI->mesh;
  ier       = 1;

  /** Add the interfaces tetra to the meshI mesh */
  meshJ                     = grpJ->mesh;
  nitem_int_face_comm       = grpJ->nitem_int_face_comm;
  face2int_face_comm_index1 = grpJ->face2int_face_comm_index1;
  face2int_face_comm_index2 = grpJ->face2int_face_comm_index2;

  ++meshJ->base;
  for ( k=0; k<nitem_int_face_comm; ++k ) {
    face_id_glo = face2int_face_comm_index2[k];
    assert( 0 <= face_id_glo );
    assert(face_id_glo < (parmesh->int_face_comm->nitem));

    /* Index of the interface tetra in the mesh */
    iel   =  face2int_face_comm_index1[k]/12;
    ifac  = (face2int_face_comm_index1[k]%12)/3;
    iploc = (face2int_face_comm_index1[k]%12)%3;

    assert ( iel && iel<=meshJ->ne );
    assert ( 0<=ifac && ifac<4 );

    ptJ       = &meshJ->tetra[iel];

    assert ( MG_EOK(ptJ) );

    /* If the tetra has already been added to meshI */
    if ( ptJ->base == meshJ->base ) {
      /* Store the interface face if it has not been seen from another group */
      assert( ptJ->flag );

      if ( !intvalues[face_id_glo] )
        intvalues[face_id_glo] = 12*ptJ->flag+3*ifac+iploc;
      else intvalues[face_id_glo] *= -1;

      continue;
    }

    /* Else */
    ptJ->base = meshJ->base;

    /* Add the tetra iel to meshI */
    ie = MMG3D_newElt(meshI);
    if ( !ie ) {
      MMG3D_TETRA_REALLOC(meshI,ie,meshI->gap,
                           fprintf(stderr,"  ## Error: unable to merge group elts.\n");
                           MMG5_INCREASE_MEM_MESSAGE();
                           return 0);
    }
    ptI = &meshI->tetra[ie];

    for ( i=0; i<4; ++i ) ptI->v[i] = meshJ->point[ptJ->v[i]].tmp;
    ptI->ref  = ptJ->ref;
    ptI->qual = ptJ->qual;
    ptI->mark = ptJ->mark;

    ptJ->flag = ie;

    /* Store the interface face if it has not been seen from another group */
    if ( !intvalues[face_id_glo] ) intvalues[face_id_glo] = 12*ie+3*ifac+iploc;
    else intvalues[face_id_glo] *= -1;

    /** Add xtetra if needed */
    if ( ptJ->xt ) {
      pxtJ = &meshJ->xtetra[ptJ->xt];
      meshI->xt++;
      if ( meshI->xt > meshI->xtmax ) {
        newsize = MG_MAX((1+meshI->gap) * meshI->xtmax,meshI->xtmax+1);
        PMMG_RECALLOC(meshI, meshI->xtetra, newsize+1,
                      meshI->xtmax+1, MMG5_xTetra,
                      "larger xtetra table", meshI->xt--; ier = 0);
        if ( ier )
          meshI->xtmax = newsize;
      }
      if ( ier ) {
        ptI->xt = meshI->xt;
        pxtI = &meshI->xtetra[ptI->xt];
        memcpy(pxtI,pxtJ,sizeof(MMG5_xTetra));
      }
      else break;
    }
  }

  return ier;
}

/**
 * \param parmesh pointer toward a PMMG parmesh structure.
 * \param grpI pointer toward the group in which we want to merge.
 * \param grpJ pointer toward the group that we want to merge into the group
 * \a grpI
 *
 * \return 0 if fail, 1 otherwise
 *
 * Merge the tetra of the \a grpJ->mesh mesh into the \a grpI->mesh mesh.
 *
 */
int PMMG_mergeGrpJinI_internalTetra( PMMG_pGrp grpI, PMMG_pGrp grpJ ) {
  MMG5_pMesh     meshI,meshJ;
  MMG5_pTetra    ptI,ptJ;
  MMG5_pxTetra   pxtI,pxtJ;
  int            k,ie,ier,i,newsize;

  meshI = grpI->mesh;

  /** Add current meshs' tetras to the merged mesh (meshI) */
  meshJ = grpJ->mesh;
  ier   = 1;

  for ( k=1; k<=meshJ->ne; k++ ) {
    ptJ  = &meshJ->tetra[k];

    if ( (!MG_EOK(ptJ)) || ptJ->base==meshJ->base ) continue;

    ptJ->base = meshJ->base;

    ie  = MMG3D_newElt(meshI);
    if ( !ie ) {
      MMG3D_TETRA_REALLOC(meshI,ie,meshI->gap,
                           fprintf(stderr,"  ## Error: unable to merge group elts.\n");
                           MMG5_INCREASE_MEM_MESSAGE();
                           return 0);
    }
    ptI = &meshI->tetra[ie];

    for ( i=0; i<4; ++i ) ptI->v[i] = meshJ->point[ptJ->v[i]].tmp;
    ptI->ref  = ptJ->ref;
    ptI->qual = ptJ->qual;
    ptI->mark = ptJ->mark;

    /** Add xtetra if needed */
    if ( ptJ->xt ) {
      pxtJ = &meshJ->xtetra[ptJ->xt];
      meshI->xt++;

      if ( meshI->xt > meshI->xtmax ) {
        newsize =  MG_MAX((1.+meshI->gap) * meshI->xtmax,meshI->xtmax+1);
        PMMG_RECALLOC(meshI, meshI->xtetra,newsize+1,
                      meshI->xtmax+1, MMG5_xTetra,
                      "larger xtetra table", meshI->xt--; ier = 0;);
        if ( ier )
          meshI->xtmax = newsize;
      }
      if ( ier ) {
        ptI->xt = meshI->xt;
        pxtI = &meshI->xtetra[ptI->xt];
        memcpy(pxtI,pxtJ,sizeof(MMG5_xTetra));
      }
      else break;
    }
  }

  return ier;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param grpI pointer toward the group in which we want to merge.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Update the node communicators when merging the groups of listgroups into 1
 * group.
 *
 */
static inline
int PMMG_mergeGrps_nodeCommunicators( PMMG_pParMesh parmesh,PMMG_pGrp grpI ) {
  MMG5_pMesh     meshI;
  MMG5_pPoint    ppt;
  PMMG_pExt_comm ext_node_comm;
  PMMG_pInt_comm int_node_comm;
  int            *intvalues;
  int           *node2int_node_comm0_index1;
  int           *node2int_node_comm0_index2;
  int            poi_id_int,poi_id_glo,idx,k,i,ip;
  int            new_nitem_int_node_comm;

  int_node_comm              = parmesh->int_node_comm;
  intvalues                  = int_node_comm->intvalues;
  meshI                      = grpI->mesh;
  node2int_node_comm0_index1 = grpI->node2int_node_comm_index1;
  node2int_node_comm0_index2 = grpI->node2int_node_comm_index2;

  /** Reset the tmp field of the point: it will be used to store the position of
   * a point in the internal communicator */
  for ( k=1; k<=meshI->np; k++ )
    meshI->point[k].tmp = PMMG_UNSET;

  /** Travel through the external communicators and udpate all the communicators */
  poi_id_int = 0;
  for ( k=0; k<parmesh->next_node_comm; ++k ) {

    /* currently working external communicator */
    ext_node_comm = &parmesh->ext_node_comm[k];

    poi_id_glo = 0;
    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx = ext_node_comm->int_comm_index[i];
      assert( (0<=idx ) && (idx<parmesh->int_node_comm->nitem) &&
              "check intvalues indices" );
      ip  = abs(intvalues[idx]);
      assert ( ip && ip<=meshI->np );

      /* The point belong to the merged mesh */
      ppt = &meshI->point[ip];

      /* New point in the internal communicator */
      if ( ppt->tmp<0 ) {
        if ( poi_id_int == grpI->nitem_int_node_comm ) {
          new_nitem_int_node_comm = (int)((1.+PMMG_GAP)*grpI->nitem_int_node_comm) + 1;
          PMMG_REALLOC(parmesh,grpI->node2int_node_comm_index1,
                       new_nitem_int_node_comm,
                       grpI->nitem_int_node_comm,int,
                       "(mergeGrps) node2int_node_comm_index1",
                       grpI->nitem_int_node_comm = new_nitem_int_node_comm;
                       return 0);
          PMMG_REALLOC(parmesh,grpI->node2int_node_comm_index2,
                       new_nitem_int_node_comm,
                       grpI->nitem_int_node_comm,int,
                       "(mergeGrps) node2int_node_comm_index2",
                       grpI->nitem_int_node_comm = new_nitem_int_node_comm;
                       return 0);
          grpI->nitem_int_node_comm = new_nitem_int_node_comm;
          node2int_node_comm0_index1 = grpI->node2int_node_comm_index1;
          node2int_node_comm0_index2 = grpI->node2int_node_comm_index2;
        }
        node2int_node_comm0_index1[ poi_id_int ]     = ip;
        node2int_node_comm0_index2[ poi_id_int ]     = poi_id_int;

        ext_node_comm->int_comm_index[ poi_id_glo++ ]= poi_id_int;
        ppt->tmp                                     = poi_id_int + k*meshI->np;
        poi_id_int++;
      } else {
        /* The point has already a position in the internal comm */
        if ( ppt->tmp < k * meshI->np ) {
          /* The point has been stored in the internal comm by another external
           * comm: update its position in our external comm */
          ext_node_comm->int_comm_index[ poi_id_glo++] = ppt->tmp%meshI->np;
        }
      }
    }
    assert( (poi_id_glo != 0)  && "empty communicator?????" );
    assert(poi_id_glo==ext_node_comm->nitem);
  }

  PMMG_REALLOC(parmesh,grpI->node2int_node_comm_index1,poi_id_int,
               grpI->nitem_int_node_comm,int,
               "(mergeGrps) node2int_node_comm_index1",return 0);
  PMMG_REALLOC(parmesh,grpI->node2int_node_comm_index2,poi_id_int,
               grpI->nitem_int_node_comm,int,
               "(mergeGrps) node2int_node_comm_index2",return 0);
  grpI->nitem_int_node_comm = poi_id_int;
  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,int,
               "free int_node_comm intvalues");
  int_node_comm->nitem       = poi_id_int;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Update the face communicators when merging the groups of listgroups into 1
 * group.
 *
 */
static inline
int PMMG_mergeGrps_faceCommunicators(PMMG_pParMesh parmesh) {
  PMMG_pGrp      grp;
  PMMG_pExt_comm ext_face_comm;
  PMMG_pInt_comm int_face_comm;
  int            *intvalues;
  int           *face2int_face_comm0_index1;
  int           *face2int_face_comm0_index2;
  int            face_id_int,idx,k,i,iel;
  int            new_nitem_int_face_comm;

  grp  = parmesh->listgrp;

  int_face_comm              = parmesh->int_face_comm;
  intvalues                  = int_face_comm->intvalues;
  face2int_face_comm0_index1 = grp[0].face2int_face_comm_index1;
  face2int_face_comm0_index2 = grp[0].face2int_face_comm_index2;

  /** Travel through the external communicators and udpate all the communicators */
  face_id_int = 0;
  for ( k=0; k<parmesh->next_face_comm; ++k ) {

    /* currently working external communicator */
    ext_face_comm = &parmesh->ext_face_comm[k];

    for ( i=0; i<ext_face_comm->nitem; ++i ) {
      idx = ext_face_comm->int_comm_index[i];
      iel = abs(intvalues[idx]);
      assert(iel);

      /* Add this face to the face communicators */
      if ( face_id_int == grp[0].nitem_int_face_comm ) {
        new_nitem_int_face_comm = (int)((1.+PMMG_GAP)*grp[0].nitem_int_face_comm) + 1;

        PMMG_REALLOC(parmesh,grp[0].face2int_face_comm_index1,
                     new_nitem_int_face_comm,
                     grp[0].nitem_int_face_comm,int,
                     "(mergeGrps) face2int_face_comm_index1",
                     grp[0].nitem_int_face_comm = new_nitem_int_face_comm;
                     return 0);
        PMMG_REALLOC(parmesh,grp[0].face2int_face_comm_index2,
                     new_nitem_int_face_comm,
                     grp[0].nitem_int_face_comm,int,
                     "(mergeGrps) face2int_face_comm_index2",
                     grp[0].nitem_int_face_comm = new_nitem_int_face_comm;
                     return 0);
        grp[0].nitem_int_face_comm = new_nitem_int_face_comm;
        face2int_face_comm0_index1 = grp[0].face2int_face_comm_index1;
        face2int_face_comm0_index2 = grp[0].face2int_face_comm_index2;
      }
      face2int_face_comm0_index1[ face_id_int ] = iel;
      face2int_face_comm0_index2[ face_id_int ] = face_id_int;
      ext_face_comm->int_comm_index[ i ]        = face_id_int;
      face_id_int++;
    }
  }

  PMMG_REALLOC(parmesh,grp[0].face2int_face_comm_index1,face_id_int,
               grp[0].nitem_int_face_comm,int,
               "(mergeGrps) face2int_face_comm_index1",return 0);
  PMMG_REALLOC(parmesh,grp[0].face2int_face_comm_index2,face_id_int,
               grp[0].nitem_int_face_comm,int,
               "(mergeGrps) face2int_face_comm_index2",return 0);
  grp[0].nitem_int_face_comm = face_id_int;
  PMMG_DEL_MEM(parmesh,int_face_comm->intvalues,int,
               "free int_face_comm intvalues");
  int_face_comm->nitem       = face_id_int;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Update the communicators when merging the groups of listgroups into 1 group.
 *
 */
static inline
int PMMG_mergeGrps_communicators(PMMG_pParMesh parmesh) {

  if ( !PMMG_mergeGrps_nodeCommunicators(parmesh,&parmesh->listgrp[0]) )
    return 0;

  if ( !PMMG_mergeGrps_faceCommunicators(parmesh) ) return 0;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Update the tag on the points and tetra
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
  size_t          available,oldMemMax;

  /* Compute available memory (previously given to the communicators) */
  PMMG_TRANSFER_AVMEM_TO_PARMESH(parmesh,available,oldMemMax);

  /* Loop on groups */
  for ( grpid=0; grpid<parmesh->ngrp; grpid++ ) {
    grp                        = &parmesh->listgrp[grpid];
    mesh                       = grp->mesh;
    node2int_node_comm0_index1 = grp->node2int_node_comm_index1;
    face2int_face_comm0_index1 = grp->face2int_face_comm_index1;

    PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh,available,oldMemMax);

    /** Step 1: Loop on xtetras to untag old parallel entities, then build
     * hash table for edges on xtetras. */
    for ( k=1; k<=mesh->ne; k++ ) {
      pt = &mesh->tetra[k];
      if ( !pt->xt ) continue;
      pxt = &mesh->xtetra[pt->xt];
      /* Untag parallel nodes */
      for ( j=0 ; j<4 ; j++ ) {
        ppt = &mesh->point[pt->v[j]];
        if ( ppt->tag & MG_PARBDY ) {
          ppt->tag &= ~MG_PARBDY;
          if ( ppt->tag & MG_BDY )    ppt->tag &= ~MG_BDY;
          if ( ppt->tag & MG_REQ )    ppt->tag &= ~MG_REQ;
          if ( ppt->tag & MG_NOSURF ) ppt->tag &= ~MG_NOSURF;
        }
      }
      /* Untag parallel edges */
      for ( j=0 ; j<6 ; j++ )
        if ( pxt->tag[j] & MG_PARBDY ) {
          pxt->tag[j] &= ~MG_PARBDY;
          if ( pxt->tag[j] & MG_BDY)    pxt->tag[j] &= ~MG_BDY;
          if ( pxt->tag[j] & MG_REQ)    pxt->tag[j] &= ~MG_REQ;
          if ( pxt->tag[j] & MG_NOSURF) pxt->tag[j] &= ~MG_NOSURF;
        }
      /* Untag parallel faces */
      for ( j=0 ; j<4 ; j++ )
        if ( pxt->ftag[j] & MG_PARBDY ) {
          pxt->ftag[j] &= ~MG_PARBDY;
          if ( pxt->ftag[j] & MG_BDY)    pxt->ftag[j] &= ~MG_BDY;
          if ( pxt->ftag[j] & MG_REQ)    pxt->ftag[j] &= ~MG_REQ;
          if ( pxt->ftag[j] & MG_NOSURF) pxt->ftag[j] &= ~MG_NOSURF;
        }
    }

    /* Create hash table for edges */
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

    /** Step 2: Re-tag boundary entities starting from xtetra faces. */
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
          pxt->ftag[ifac] |= MG_BDY;
        }
        /* Only a "true" boundary after this line */
        if ( pxt->ftag[ifac] & MG_BDY ) {
          /* Constrain boundary if -nosurf option */
          if( mesh->info.nosurf ) pxt->ftag[ifac] |= MG_REQ + MG_NOSURF;
          /* Tag face edges */
          for ( j=0; j<3; j++ ) {
            ia = MMG5_iarf[ifac][j];
            ip0 = pt->v[MMG5_iare[ia][0]];
            ip1 = pt->v[MMG5_iare[ia][1]];
            if( !MMG5_hTag( &hash, ip0, ip1, 0, MG_BDY ) ) return 0;
            /* Constrain boundary if -nosurf option */
            if( mesh->info.nosurf )
              if( !MMG5_hTag( &hash, ip0, ip1, 0, MG_REQ + MG_NOSURF ) ) return 0;
          }
          /* Tag face nodes */
          for ( j=0 ; j<3 ; j++) {
            ppt = &mesh->point[pt->v[MMG5_idir[ifac][j]]];
            ppt->tag |= MG_BDY;
            /* Constrain boundary if -nosurf option */
            if( mesh->info.nosurf ) ppt->tag |= MG_REQ + MG_NOSURF;
          }
        }
      }
    }

    /** Step 3: if communicators are allocated: tag new parallel interface
     * entities starting from int_face_comm. */
    if ( parmesh->ext_face_comm ) {
      for ( i=0; i<grp->nitem_int_face_comm; i++ ) {
        iel  =   face2int_face_comm0_index1[i] / 12;
        ifac = ( face2int_face_comm0_index1[i] % 12 ) / 3;
        pt = &mesh->tetra[iel];
        assert( pt->xt );
        pxt = &mesh->xtetra[pt->xt];
        /* If already boundary, make it recognizable as a "true" boundary */
        if( pxt->ftag[ifac] & MG_BDY ) pxt->ftag[ifac] |= MG_PARBDYBDY;
        /* Tag face */
        pxt->ftag[ifac] |= (MG_PARBDY + MG_BDY + MG_REQ + MG_NOSURF);
        /* Tag face edges */
        for ( j=0; j<3; j++ ) {
          ia = MMG5_iarf[ifac][j];
          ip0 = pt->v[MMG5_iare[ia][0]];
          ip1 = pt->v[MMG5_iare[ia][1]];
          if( !MMG5_hTag( &hash, ip0, ip1, 0,
                          MG_PARBDY + MG_BDY + MG_REQ + MG_NOSURF ) ) return 0;
        }
        /* Tag face nodes */
        for ( j=0 ; j<3 ; j++) {
          ppt = &mesh->point[pt->v[MMG5_idir[ifac][j]]];
          ppt->tag |= (MG_PARBDY + MG_BDY + MG_REQ + MG_NOSURF);
        }
      }
    }

    /** Step 4: Get edge tag and delete hash table */
    for ( k=1; k<=mesh->ne; k++ ) {
      pt = &mesh->tetra[k];
      if ( !pt->xt ) continue;
      pxt = &mesh->xtetra[pt->xt];
      for ( j=0; j<6; j++ ) {
        ip0 = pt->v[MMG5_iare[j][0]];
        ip1 = pt->v[MMG5_iare[j][1]];
        /* Put the tag stored in the hash table on the xtetra edge */
        if( !MMG5_hGet( &hash, ip0, ip1, &getref, &pxt->tag[j] ) ) return 0;
      }
    }
    PMMG_DEL_MEM( mesh, hash.geom, MMG5_hgeom, "Edge hash table" );

    /** Step 5: Unreference xpoints not on BDY (or PARBDY) */
    for ( i=1; i<=mesh->np; i++ ) {
      ppt = &mesh->point[i];
      if( ppt->tag & MG_BDY ) continue;
      if( ppt->xp ) ppt->xp = 0;
    }

    PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,mesh,available,oldMemMax);
  }

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if fail, 1 if success
 *
 * Merge all meshes (mesh elements + internal communicator) of a group into the
 * first mesh of the group. This function free the adjacency array.
 *
 * \remark the tetra must be packed.
 *
 */
int PMMG_merge_grps( PMMG_pParMesh parmesh )
{
  PMMG_pGrp      listgrp,grp;
  MMG5_pMesh     mesh0;
  PMMG_pInt_comm int_node_comm,int_face_comm;
  size_t         available,oldMemMax;
  int            *face2int_face_comm_index1,*face2int_face_comm_index2;
  int            imsh,k,iel;

  if ( !parmesh->ngrp ) return 1;

  listgrp  = parmesh->listgrp;

  /** Free the adjacency array: a possible improvement is to update it */
  mesh0 = listgrp[0].mesh;

  if ( !mesh0 ) return 1;

  if ( mesh0->adja )
    PMMG_DEL_MEM(mesh0, mesh0->adja,int, "adjacency table" );

  if ( parmesh->ngrp == 1 ) return 1;

  /* Give the memory to the parmesh */
  PMMG_TRANSFER_AVMEM_TO_PARMESH(parmesh,available,oldMemMax);

  /** Use the internal communicators to store the interface entities indices */
  int_node_comm = parmesh->int_node_comm;
  PMMG_CALLOC(parmesh,int_node_comm->intvalues,int_node_comm->nitem,int,
              "node communicator",return 0);

  int_face_comm = parmesh->int_face_comm;
  PMMG_CALLOC(parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,
              "face communicator",goto fail_ncomm);

  //DEBUGGING:
  //saveGrpsToMeshes(listgrp,parmesh->ngrp,parmesh->myrank,"BeforeMergeGrp");

  /* Give all the memory to mesh0 */
  PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh0,available,oldMemMax);

  /** Step 0: Store the indices of the interface faces of mesh0 into the
   * internal face communicator */
  face2int_face_comm_index1 = listgrp[0].face2int_face_comm_index1;
  face2int_face_comm_index2 = listgrp[0].face2int_face_comm_index2;
  for ( k=0; k<listgrp[0].nitem_int_face_comm; ++k ) {
    iel = face2int_face_comm_index1[k];
    assert(   ( 0 <= face2int_face_comm_index2[k] )
              && ( face2int_face_comm_index2[k] < parmesh->int_face_comm->nitem )
              && "check intvalues indices" );
    parmesh->int_face_comm->intvalues[face2int_face_comm_index2[k]] = iel;
  }

  /** Step 1: Merge interface points from all meshes into mesh0->points */
  if ( !PMMG_mergeGrps_interfacePoints(parmesh) ) goto fail_comms;

  for ( imsh=1; imsh<parmesh->ngrp; ++imsh ) {
    grp = &listgrp[imsh];

    /** Step 2: Merge internal points of the mesh mesh into the mesh0 mesh */
    if ( !PMMG_mergeGrpJinI_internalPoints(&listgrp[0],grp) )
      goto fail_comms;

    /* Step 3: Add the interfaces tetra of the imsh mesh to the mesh0 mesh */
    if ( !PMMG_mergeGrpJinI_interfaceTetra(parmesh,&listgrp[0],grp) )
      goto fail_comms;

    /** Step 4: Merge internal tetras of the imsh mesh into the mesh0 mesh */
    if ( !PMMG_mergeGrpJinI_internalTetra(&listgrp[0],grp) )
      goto fail_comms;

    /* Free merged mesh and increase mesh0->memMax*/
    mesh0->memMax += grp->mesh->memCur;
    PMMG_grp_free(parmesh,grp);
  }
  assert ( mesh0->memMax+parmesh->memMax<=parmesh->memGloMax );

  /** Step 5: Update the communicators */
  /* Give all the memory to the communicators */
  PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,mesh0,available,oldMemMax);

  if ( !PMMG_mergeGrps_communicators(parmesh) ) goto fail_comms;

  PMMG_REALLOC(parmesh,parmesh->listgrp,1,parmesh->ngrp,PMMG_Grp,"listgrp",return 0;);
  parmesh->ngrp = 1;

  /** Step 6: Update tag on points, tetra */
  if ( !PMMG_updateTag(parmesh) ) goto fail_comms;

  return 1;

fail_comms:
  PMMG_DEL_MEM(parmesh,int_face_comm->intvalues,int,"face communicator");

fail_ncomm:
  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,int,"node communicator");

  return 0;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param rcv_point Buffer to gather points
 * \param np_tot number of points in \a rcv_point
 * \param rcv_xpoint Buffer to gather xPoints
 * \param xp_tot number of xpoints in \a rcv_xpoint
 * \param rcv_tetra Buffer to gather tetra
 * \param ne_tot number of tetra in \a rcv_tetra
 * \param rcv_xtetra Buffer to gather xtetra
 * \param xt_tot number of xtetra in \a rcv_xtetra
 * \param rcv_met buffer to gather metric
 * \param nmet_tot number of metrics in \a rcv_met
 * \param rcv_isMet met->size if we are waiting a metric
 * \param rcv_intvalues Buffer to gather the intvalue array of the internal comm
 * \param nitem_int_node_comm_tot number of items in \a  rcv_intvalues
 * \param rcv_nitem_ext_tab Buffer to gather the number of item in the ext comm
 * \param ext_comm_displs_tot number of items in \a rcv_nitem_ext_tab
 * \param rcv_color_in_tab Buffer to gather the color_in field of the ext comm
 * \param rcv_color_out_tab Buffer to gather the color_out field of the ext comm
 * \param rcv_node2int_node_comm_index1 Buffer to gather node2int_node_comm_index1
 * \param rcv_node2int_node_comm_index2 Buffer to gather node2int_node_comm_index2
 * \param point_displs Position of the 1st point of each mesh in rcv_point
 * \param xpoint_displs Position of the 1st xpoint of each mesh in rcv_xpoint
 * \param tetra_displs Position of the 1st tetra of each mesh in rcv_tetra
 * \param xtetra_displs Position of the 1st xtetra of each mesh in rcv_xtetra
 * \param met_displs Position of the 1st metric of each mesh in rcv_met
 * \param intval_displs Position of the 1st data of each internal comm in rcv_intvalues
 * \param ext_comm_displs Position of the 1st data of each external comm in arrays
 * related to external comm
 * \param int_comm_index_displs Position of the 1st data of each internal comm in
 * the rcv_node2int_node_comm_index arrays
 * \param rcv_np Buffer to gather the number of points
 * \param rcv_xp Buffer to gather the number of xPoints
 * \param rcv_ne Buffer to gather the number of tetra
 * \param rcv_xt Buffer to gather the number of xtetra
 * \param rcv_nmet Buffer to gather the number of metrics
 * \param rcv_int_comm_index Buffer to gather the internal comm sizes
 * \param nitem_icidx_tot number of items in \a rcv_int_comm_index
 * \param rcv_next_node_comm Buffer to gather the numbers of external comm
 * \param rcv_nitem_int_node_comm Buffer to gather the node2int_node_comm arrays
 * sizes
 *
 * \return 0 if fail, 1 otherwise (on all procs)
 *
 * Gather the parmeshes on the proc 0.
 *
 * \warning We must have 1 group per parmesh
 *
 */
static inline
int PMMG_gather_parmesh( PMMG_pParMesh parmesh,
                         MMG5_pPoint  *rcv_point, int *np_tot,
                         MMG5_pxPoint *rcv_xpoint,int *xp_tot,
                         MMG5_pTetra  *rcv_tetra, int *ne_tot,
                         MMG5_pxTetra *rcv_xtetra,int *xt_tot,
                         double      **rcv_met,   int *nmet_tot,int *rcv_isMet,
                         int **rcv_intvalues,int *nitem_int_node_comm_tot,
                         int **rcv_nitem_ext_tab,int *ext_comm_displs_tot,
                         int **rcv_color_in_tab,int **rcv_color_out_tab,
                         int **rcv_node2int_node_comm_index1,
                         int **rcv_node2int_node_comm_index2,int **point_displs,
                         int **xpoint_displs,int **tetra_displs,
                         int **xtetra_displs,int **met_displs,
                         int **intval_displs,int **ext_comm_displs,
                         int **int_comm_index_displs,int **rcv_np,int **rcv_xp,
                         int **rcv_ne,int **rcv_xt,int **rcv_nmet,
                         int **rcv_int_comm_index,int *nitem_icidx_tot,
                         int** rcv_next_node_comm,
                         int **rcv_nitem_int_node_comm ) {
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_pSol      met;
  PMMG_pInt_comm int_node_comm;
  PMMG_pExt_comm ext_node_comm;
  MPI_Comm       comm;
  MPI_Datatype   mpi_point,mpi_xpoint,mpi_tetra,mpi_xtetra;
  int            *color_in_tab,*color_out_tab,*nitem_ext_tab;
  int            *int_comm_index,*nitems_ext_idx,nitem_ext_tot;
  int            nprocs,root,k,i,idx;
  int            isMet,np,xp,ne,xt,size2send;
  int            ier,ieresult;
  void           *ptr;

  nprocs        = parmesh->nprocs;
  root          = parmesh->info.root;
  grp           = &parmesh->listgrp[0];
  mesh          = grp ? grp->mesh : NULL;
  met           = grp ? grp->met  : NULL;
  isMet         = (met && (met->m) )? met->size : 0;
  comm          = parmesh->comm;
  int_node_comm = parmesh->int_node_comm;
  ier           = 1;

  (*rcv_np)        = NULL;
  (*point_displs)  = NULL;
  (*rcv_point)     = NULL;

  (*rcv_xp)        = NULL;
  (*xpoint_displs) = NULL;
  (*rcv_xpoint)    = NULL;

  (*rcv_ne)        = NULL;
  (*tetra_displs)  = NULL;
  (*rcv_tetra)     = NULL;

  (*rcv_xt)        = NULL;
  (*xtetra_displs) = NULL;
  (*rcv_xtetra)    = NULL;

  (*rcv_nmet)      = NULL;
  (*met_displs)    = NULL;
  (*rcv_met)       = NULL;

  (*rcv_intvalues)                 = NULL;
  (*rcv_node2int_node_comm_index1) = NULL;
  (*rcv_node2int_node_comm_index2) = NULL;
  (*rcv_int_comm_index)            = NULL;
  (*rcv_next_node_comm)            = NULL;
  (*rcv_nitem_ext_tab)             = NULL;

  (*rcv_color_in_tab)  = NULL;
  (*rcv_color_out_tab) = NULL;
  (*ext_comm_displs )  = NULL;

  color_in_tab   = NULL;
  color_out_tab  = NULL;
  nitem_ext_tab  = NULL;
  nitems_ext_idx = NULL;
  int_comm_index = NULL;

  /** Memory alloc */
  /* 1: Mesh data */
  if ( !parmesh->myrank ) {
    PMMG_CALLOC( parmesh, (*rcv_np)        ,nprocs,int,"rcv_np"       ,ier=0);
    PMMG_CALLOC( parmesh, (*point_displs)  ,nprocs,int,"point_displs" ,ier=0);

    PMMG_CALLOC( parmesh, (*rcv_xp)        ,nprocs,int,"rcv_xp"       ,ier=3);
    PMMG_CALLOC( parmesh, (*xpoint_displs) ,nprocs,int,"xpoint_displs",ier=3);

    PMMG_CALLOC( parmesh, (*rcv_ne)        ,nprocs,int,"rcv_ne"       ,ier=3);
    PMMG_CALLOC( parmesh, (*tetra_displs)  ,nprocs,int,"tetra_displs" ,ier=3);

    PMMG_CALLOC( parmesh, (*rcv_xt)        ,nprocs,int,"rcv_xt"       ,ier=3);
    PMMG_CALLOC( parmesh, (*xtetra_displs) ,nprocs,int,"xtetra_displs",ier=3);
  }

  /* 2: Communicators */
  PMMG_CALLOC( parmesh, (*intval_displs)          ,nprocs,int,"intval_displs"        ,ier=3);
  PMMG_CALLOC( parmesh, (*rcv_nitem_int_node_comm),nprocs,int,"rcv_nitem_inc0"       ,ier=3);
  PMMG_CALLOC( parmesh, (*int_comm_index_displs)  ,nprocs,int,"int_comm_index_displs",ier=3);
  PMMG_CALLOC( parmesh, (*rcv_next_node_comm)     ,nprocs,int,"rcv_next_node_comm"   ,ier=3);
  PMMG_CALLOC( parmesh, (*ext_comm_displs)        ,nprocs,int,"ext_comm_displs"      ,ier=3);

  /** Gather parmesh size infos on proc 0 */
  np = mesh ? mesh->np : 0;
  xp = mesh ? mesh->xp : 0;
  ne = mesh ? mesh->ne : 0;
  xt = mesh ? mesh->xt : 0;

  MPI_CHECK( MPI_Gather(&np,1,MPI_INT,(*rcv_np),1,MPI_INT,root,comm),ier = 3);
  MPI_CHECK( MPI_Gather(&ne,1,MPI_INT,(*rcv_ne),1,MPI_INT,root,comm),ier = 3);
  MPI_CHECK( MPI_Gather(&xp,1,MPI_INT,(*rcv_xp),1,MPI_INT,root,comm),ier = 3);
  MPI_CHECK( MPI_Gather(&xt,1,MPI_INT,(*rcv_xt),1,MPI_INT,root,comm),ier = 3);

  size2send = grp ? grp->nitem_int_node_comm : 0;
  MPI_CHECK( MPI_Gather(&size2send,1,MPI_INT,
                        (*rcv_nitem_int_node_comm),1,MPI_INT,root,comm),ier = 3);
  MPI_CHECK( MPI_Gather(&parmesh->next_node_comm,1,MPI_INT,
                        (*rcv_next_node_comm),1,MPI_INT,root,comm),ier = 3);
  MPI_CHECK( MPI_Allreduce(&isMet,rcv_isMet,1,MPI_INT,MPI_MAX,comm),ier = 3);

  if ( (!parmesh->myrank) && rcv_isMet ) {
    PMMG_CALLOC( parmesh, (*rcv_nmet)    ,nprocs,int,"rcv_nmet"     ,ier=3);
    PMMG_CALLOC( parmesh, (*met_displs)  ,nprocs,int,"met_displs"   ,ier=3);
  }

  MPI_CHECK( MPI_Allreduce(&ier,&ieresult,1,MPI_INT,MPI_MAX,comm),ieresult=3);
  if ( ieresult>1 ) goto end;

  /** Alloc arrays */
  /* Creation of MPI types for the mesh comm */
  PMMG_create_MPI_Point ( &mpi_point  );
  PMMG_create_MPI_xPoint( &mpi_xpoint );
  PMMG_create_MPI_Tetra ( &mpi_tetra  );
  PMMG_create_MPI_xTetra( &mpi_xtetra );

  /* Points */
  if ( mesh )
    for ( k=1; k<=mesh->np; ++k ) mesh->point[k].tmp = 0;

  if ( !parmesh->myrank ) {
    (*point_displs)[0] = 0;
    for ( k=1; k<nprocs; ++k ) {
      (*point_displs)[k] += (*point_displs)[k-1] + (*rcv_np)[k-1];
    }
    *np_tot       = (*point_displs)[nprocs-1]+(*rcv_np)[nprocs-1];
    PMMG_CALLOC( parmesh,(*rcv_point),*np_tot+1,MMG5_Point,"rcv_point",ier=2);
  }
  /* xPoints */
  if ( !parmesh->myrank ) {
    (*xpoint_displs)[0] = 0;
    for ( k=1; k<nprocs; ++k ) {
      (*xpoint_displs)[k] += (*xpoint_displs)[k-1] + (*rcv_xp)[k-1];
    }
    *xp_tot        = (*xpoint_displs)[nprocs-1]+(*rcv_xp)[nprocs-1];
    PMMG_CALLOC( parmesh,(*rcv_xpoint),*xp_tot+1,MMG5_xPoint,"rcv_xpoint",ier=2);
  }
  /* Tetra */
  if ( !parmesh->myrank ) {
    (*tetra_displs)[0] = 0;
    for ( k=1; k<nprocs; ++k ) {
      (*tetra_displs)[k] += (*tetra_displs)[k-1] + (*rcv_ne)[k-1];
    }
    *ne_tot       = (*tetra_displs)[nprocs-1]+(*rcv_ne)[nprocs-1];
    PMMG_CALLOC( parmesh,(*rcv_tetra),*ne_tot+1,MMG5_Tetra,"rcv_tetra",ier=2);
  }
  /* xTetra */
  if ( !parmesh->myrank ) {
    (*xtetra_displs)[0] = 0;
    for ( k=1; k<nprocs; ++k ) {
      (*xtetra_displs)[k] += (*xtetra_displs)[k-1] + (*rcv_xt)[k-1];
    }
    *xt_tot        = (*xtetra_displs)[nprocs-1]+(*rcv_xt)[nprocs-1];
    PMMG_CALLOC( parmesh,(*rcv_xtetra),*xt_tot+1,MMG5_xTetra,"rcv_xtetra",ier=2);
  }
  /* Solutions */
  if ( *rcv_isMet && !parmesh->myrank ) {
    for ( k=0; k<nprocs; ++k ) {
      /* isMet = met->size here */
      (*rcv_nmet)[k] = *rcv_isMet*(*rcv_np)[k];
    }

    (*met_displs)[0] = 0;
    for ( k=1; k<nprocs; ++k ) {
      (*met_displs)[k] += (*met_displs)[k-1] + (*rcv_nmet)[k-1];
    }
    *nmet_tot   = (*met_displs)[nprocs-1]+(*rcv_nmet)[nprocs-1]+*rcv_isMet;
    PMMG_CALLOC( parmesh,(*rcv_met),*nmet_tot,double,"rcv_met",ier=2);
  }

  /* Internal communicator */
  (*intval_displs)[0] = 0;
  for ( k=1; k<nprocs; ++k ) {
    (*intval_displs)[k] += (*intval_displs)[k-1] + (*rcv_nitem_int_node_comm)[k-1];
  }
  *nitem_int_node_comm_tot = (*intval_displs)[nprocs-1]
    + (*rcv_nitem_int_node_comm)[nprocs-1];

  PMMG_MALLOC(parmesh,(*rcv_intvalues)                ,*nitem_int_node_comm_tot,
                int,"rcv_intvalues",ier=2 );
  PMMG_MALLOC(parmesh,(*rcv_node2int_node_comm_index1),*nitem_int_node_comm_tot,
                int,"rcv_n2inc_idx1",ier=2);
  PMMG_MALLOC(parmesh,(*rcv_node2int_node_comm_index2),*nitem_int_node_comm_tot,
                int,"rcv_n2inc_idx2",ier=2);

  /* External communicator */
  /* For each proc, we concatenate the external communicators fields into arrays */
  PMMG_MALLOC(parmesh,color_in_tab ,parmesh->next_node_comm,int,"color_in" ,ier=2);
  PMMG_MALLOC(parmesh,color_out_tab,parmesh->next_node_comm,int,"color_out",ier=2);
  PMMG_MALLOC(parmesh,nitem_ext_tab,parmesh->next_node_comm,int,"nitem_ext",ier=2);

  nitem_ext_tot = 0;
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm    = &parmesh->ext_node_comm[k];
    color_in_tab[k ] = ext_node_comm->color_in;
    color_out_tab[k] = ext_node_comm->color_out;
    nitem_ext_tab[k] = ext_node_comm->nitem;
    nitem_ext_tot += ext_node_comm->nitem;
  }

  PMMG_MALLOC(parmesh,int_comm_index,nitem_ext_tot,int,"int_comm_idx",ier=2);
  idx = 0;
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    memcpy(&int_comm_index[idx],ext_node_comm->int_comm_index,
           ext_node_comm->nitem*sizeof(int));
    idx += ext_node_comm->nitem;
  }

  (*ext_comm_displs)[0] = 0;
  for ( k=1; k<nprocs; ++k ) {
    (*ext_comm_displs)[k] += (*ext_comm_displs)[k-1] + (*rcv_next_node_comm)[k-1];
  }
  *ext_comm_displs_tot = (*ext_comm_displs)[nprocs-1] + (*rcv_next_node_comm)[nprocs-1];

  PMMG_MALLOC(parmesh,(*rcv_color_in_tab) ,*ext_comm_displs_tot,int,"rcv_color_in" ,ier=2);
  PMMG_MALLOC(parmesh,(*rcv_color_out_tab),*ext_comm_displs_tot,int,"rcv_color_out",ier=2);
  PMMG_MALLOC(parmesh,(*rcv_nitem_ext_tab),*ext_comm_displs_tot,int,"rcv_nitem_ext",ier=2);
  PMMG_MALLOC(parmesh,nitems_ext_idx,nprocs,int,"nitems_ext_idx",ier=2);

  MPI_CHECK( MPI_Gatherv(nitem_ext_tab,parmesh->next_node_comm,MPI_INT,
                         (*rcv_nitem_ext_tab),(*rcv_next_node_comm),
                         (*ext_comm_displs),MPI_INT,root,comm),ier=2);

  (*int_comm_index_displs)[0] = 0;
  nitems_ext_idx[0] = 0;
  idx = 0;
  for ( i=0; i<(*rcv_next_node_comm)[0]; ++i ) {
    nitems_ext_idx[0] += (*rcv_nitem_ext_tab)[idx++];
  }
  for ( k=1; k<nprocs; ++k ) {
    nitems_ext_idx[k] = 0;
    for ( i=0; i<(*rcv_next_node_comm)[k]; ++i ) {
      nitems_ext_idx[k] += (*rcv_nitem_ext_tab)[idx++];
    }
    (*int_comm_index_displs)[k] = (*int_comm_index_displs)[k-1] + nitems_ext_idx[k-1];
  }
  *nitem_icidx_tot = (*int_comm_index_displs)[nprocs-1]+nitems_ext_idx[nprocs-1];

  PMMG_MALLOC(parmesh,(*rcv_int_comm_index),*nitem_icidx_tot,int,"rcv_int_comm_idx",ier=2);

  MPI_CHECK( MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MAX, parmesh->comm ),
             ieresult=2 );

  if ( ieresult > 1 ) goto end;

  /** Gathers data on root */


  /* Points */
  if ( !mesh ) {
    ptr       = NULL;
    size2send = 0;
  }
  else {
    ptr = &mesh->point[1];
    size2send = mesh->np;
  }

  MPI_CHECK( MPI_Gatherv(ptr,size2send,mpi_point,&(*rcv_point)[1],
                         (*rcv_np),(*point_displs),mpi_point,root,comm),ier=2);

  /* xPoints */
  if ( !mesh ) {
    ptr       = NULL;
    size2send = 0;
  }
  else {
    ptr       = &mesh->xpoint[1];
    size2send = mesh->xp;
  }

  MPI_CHECK( MPI_Gatherv(ptr,size2send,mpi_xpoint,&(*rcv_xpoint)[1],
                         (*rcv_xp),(*xpoint_displs),mpi_xpoint,root,comm),ier=2);
  /* Tetra */
  if ( !mesh ) {
    ptr       = NULL;
    size2send = 0;
  }
  else {
    ptr       = &mesh->tetra[1];
    size2send = mesh->ne;
  }

  MPI_CHECK( MPI_Gatherv(ptr,size2send,mpi_tetra,&(*rcv_tetra)[1],
                         (*rcv_ne),(*tetra_displs),mpi_tetra,root,comm),ier=2);
  /* xTetra */
  if ( !mesh ) {
    ptr       = NULL;
    size2send = 0;
  }
  else {
    ptr       = &mesh->xtetra[1];
    size2send = mesh->xt;
  }

  MPI_CHECK( MPI_Gatherv(ptr,size2send,mpi_xtetra,&(*rcv_xtetra)[1],
                         (*rcv_xt),(*xtetra_displs),mpi_xtetra,root,comm),ier=2);
  /* Solutions */
  if ( *rcv_isMet ) {
    if ( !met ) {
      ptr       = NULL;
      size2send = 0;
    }
    else {
      ptr = &met->m[met->size];
      size2send = mesh->np*met->size;
    }

    MPI_CHECK( MPI_Gatherv(ptr,size2send,MPI_DOUBLE,
                           &(*rcv_met)[*rcv_isMet],(*rcv_nmet),(*met_displs),
                           MPI_DOUBLE,root,comm),ier=2);
  }

  /* Internal communicator */
  MPI_CHECK( MPI_Gatherv(int_node_comm->intvalues,int_node_comm->nitem,MPI_INT,
                         (*rcv_intvalues),(*rcv_nitem_int_node_comm),
                         (*intval_displs),MPI_INT,root,comm),ier=2);

  if ( !grp ) {
    ptr       = NULL;
    size2send = 0;
  }
  else {
    ptr       = grp->node2int_node_comm_index1;
    size2send = int_node_comm->nitem;
  }
  MPI_CHECK( MPI_Gatherv(ptr,size2send,MPI_INT,
                         (*rcv_node2int_node_comm_index1),(*rcv_nitem_int_node_comm),
                         (*intval_displs),MPI_INT,root,comm),ier=2);

  if ( !grp ) {
    ptr       = NULL;
    size2send = 0;
  }
  else {
    ptr       = grp->node2int_node_comm_index2;
    size2send = int_node_comm->nitem;
  }
  MPI_CHECK( MPI_Gatherv(ptr,size2send,MPI_INT,
                         (*rcv_node2int_node_comm_index2),(*rcv_nitem_int_node_comm),
                         (*intval_displs),MPI_INT,root,comm),ier=2);


  /* External communicator */
  /* For each proc, we concatenate the external communicators fields into arrays */
  MPI_CHECK( MPI_Gatherv(color_in_tab,parmesh->next_node_comm,MPI_INT,
                         (*rcv_color_in_tab),(*rcv_next_node_comm),
                         (*ext_comm_displs),MPI_INT,root,comm),ier=2);
  MPI_CHECK( MPI_Gatherv(color_out_tab,parmesh->next_node_comm,MPI_INT,
                         (*rcv_color_out_tab),(*rcv_next_node_comm),
                         (*ext_comm_displs),MPI_INT,root,comm),ier=2);
  MPI_CHECK( MPI_Gatherv(int_comm_index,nitem_ext_tot,MPI_INT,
                         (*rcv_int_comm_index),nitems_ext_idx,(*int_comm_index_displs),MPI_INT,
                         root,comm),ier=2);
  MPI_CHECK( MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MAX, parmesh->comm ),
             ieresult=2 );

end:
  /** Free the memory */
  if ( ieresult < 4 ) {
    /* Free temporary arrays */
    PMMG_DEL_MEM(parmesh,color_in_tab  ,int,"color_in");
    PMMG_DEL_MEM(parmesh,color_out_tab ,int,"color_out");
    PMMG_DEL_MEM(parmesh,nitem_ext_tab ,int,"nitem_ext");
    PMMG_DEL_MEM(parmesh,nitems_ext_idx,int,"nitems_ext_idx");
    PMMG_DEL_MEM(parmesh,int_comm_index,int,"int_comm_idx");

    /* Free useless groups */
    /* 1: mesh */
    if ( mesh ) {
      PMMG_DEL_MEM(mesh,mesh->point ,MMG5_Point ,"point");
      PMMG_DEL_MEM(mesh,mesh->tetra ,MMG5_Tetra ,"tetra");
      PMMG_DEL_MEM(mesh,mesh->xpoint,MMG5_xPoint,"xpoint");
      PMMG_DEL_MEM(mesh,mesh->xtetra,MMG5_xTetra,"xtetra");
      if ( mesh->adja ) {
        PMMG_DEL_MEM( mesh,mesh->adja,int,"adja");
      }
      if ( isMet )
        PMMG_DEL_MEM(mesh,met->m,double,"met");

      mesh->np = mesh->npmax = 0;
      mesh->ne = mesh->nemax = 0;
      mesh->xp = mesh->xpmax = 0;
      mesh->xt = mesh->xtmax = 0;
      met->np = met->npmax = 0;
    }

    /* 2: communicators */
    if ( grp ) {
      PMMG_DEL_MEM(parmesh,grp->edge2int_edge_comm_index1,int,"e2iec_idx1");
      PMMG_DEL_MEM(parmesh,grp->edge2int_edge_comm_index2,int,"e2iec_idx2");
    }

    PMMG_DEL_MEM(parmesh,parmesh->int_node_comm->intvalues,int,"intval");

    for ( i=0; i<parmesh->next_node_comm; ++i ) {
      PMMG_DEL_MEM(parmesh,parmesh->ext_node_comm[i].int_comm_index,int,"int_comm_idx");
    }
    PMMG_DEL_MEM(parmesh,parmesh->ext_node_comm,PMMG_Ext_comm,"ext_node_comm");

    if ( ieresult < 3 ) {
      /* Free MPI types */
      MPI_Type_free(&mpi_point);
      MPI_Type_free(&mpi_xpoint);
      MPI_Type_free(&mpi_tetra);
      MPI_Type_free(&mpi_xtetra);
    }
  }
  return ieresult==1;
}


/**
 * \param parmesh pointer toward the parmesh structure.
 * \param rcv_point Buffer that gathers points
 * \param rcv_xpoint Buffer that gathers xPoints
 * \param rcv_tetra Buffer that gathers tetra
 * \param rcv_xtetra Buffer that gathers xtetra
 * \param rcv_met buffer that gathers metric
 * \param rcv_intvalues Buffer that gathers the intvalue array of the internal comm
 * \param rcv_nitem_ext_tab Buffer that gathers the number of item in the ext comm
 * \param rcv_color_in_tab Buffer that gathers the color_in field of the ext comm
 * \param rcv_color_out_tab Buffer that gathers the color_out field of the ext comm
 * \param point_displs Position of the 1st point of each mesh in rcv_point
 * \param xpoint_displs Position of the 1st xpoint of each mesh in rcv_xpoint
 * \param tetra_displs Position of the 1st tetra of each mesh in rcv_tetra
 * \param xtetra_displs Position of the 1st xtetra of each mesh in rcv_xtetra
 * \param met_displs Position of the 1st metric of each mesh in rcv_met
 * \param intval_displs Position of the 1st data of each internal comm in rcv_intvalues
 * \param ext_comm_displs Position of the 1st data of each external comm in arrays
 * related to external comm
 * \param int_comm_index_displs Position of the 1st data of each internal comm in
 * the rcv_node2int_node_comm_index arrays
 * \param rcv_np Buffer that gathers the number of points
 * \param rcv_ne Buffer that gathers the number of tetra
 * \param rcv_xt Buffer that gathers the number of xtetra
 * \param rcv_int_comm_index Buffer that gathers the internal comm sizes
 * \param rcv_next_node_comm Buffer that gathers the numbers of external comm
 * sizes
 *
 * \return 0 if fail, 1 otherwise
 *
 * Merge the parmeshes data contained in the rcv_* arrays into 1 parmesh with 1
 * group.
 *
 */
static inline
int PMMG_mergeParmesh_rcvParMeshes(PMMG_pParMesh parmesh,MMG5_pPoint rcv_point,
                                   MMG5_pxPoint rcv_xpoint,MMG5_pTetra rcv_tetra,
                                   MMG5_pxTetra rcv_xtetra,double *rcv_met,
                                   int *rcv_intvalues,int *rcv_nitem_ext_tab,
                                   int *rcv_color_in_tab,int *rcv_color_out_tab,
                                   int *point_displs,
                                   int *xpoint_displs,int *tetra_displs,
                                   int *xtetra_displs,int *met_displs,
                                   int *intval_displs,int *ext_comm_displs,
                                   int *int_comm_index_displs,int *rcv_np,
                                   int *rcv_ne,int *rcv_xt,
                                   int *rcv_int_comm_index,int* rcv_next_node_comm) {
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_pPoint    point_1,point_2,ppt;
  MMG5_pxPoint   xpoint,pxp;
  MMG5_pTetra    tetra,pt;
  MMG5_pxTetra   xtetra,pxt;
  MMG5_pSol      met;
  double         *met_1;
  size_t         memAv;
  int            *int_comm_index,*int_comm_index_2;
  int            *intvalues_1,*intvalues_2,nitems_1,nitems_2;
  int            nprocs,k,i,j,idx,idx_2,cursor,color_in,color_out;
  int            np,ne,ne_tot,xt_tot,nnpar;

  nprocs = parmesh->nprocs;

  if ( parmesh->myrank != parmesh->info.root ) return 1;

  /** Give all the memory to the mesh */
  parmesh->memMax  = parmesh->memCur;

  if ( !parmesh->listgrp ) {

    if ( parmesh->memMax + sizeof(PMMG_Grp) <= parmesh->memGloMax )
      parmesh->memMax += sizeof(PMMG_Grp);

    PMMG_CALLOC(parmesh,parmesh->listgrp,1,PMMG_Grp,"listgrp", return 0);

    MMG3D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &parmesh->listgrp[0].mesh,
                     MMG5_ARG_ppMet, &parmesh->listgrp[0].met, MMG5_ARG_end );

  }

  grp    = &parmesh->listgrp[0];
  mesh   = grp->mesh;
  met    = grp->met ;

  /** Give all the memory to the mesh */
  memAv            = parmesh->memGloMax - parmesh->memMax - mesh->memMax;
  mesh->memMax    += memAv;

  np = 0;

  /** Mesh renumbering to have the same indices at points shared by multiple
   * processors; The new indices are stored in the tmp field of the MMG5_Point
   * structure. */
  for ( k=0; k<nprocs; ++k ) {
    /* Points and internal communicators for proc \a k */
    point_1     = &rcv_point[point_displs[k]];
    cursor      = intval_displs[k];
    intvalues_1 = &rcv_intvalues[cursor];

    /* Travel through the external communicators that lists the points at the
     * interface of the procs color_in and color_out: if color_in<color_out,
     * create new indices for the points over the color_in proc. If
     * color_out<color_in, get the indices from the matching points over the
     * proc color_out. */
    idx = int_comm_index_displs[k];

    for ( i=0; i<rcv_next_node_comm[k]; ++i ) {
      cursor   = ext_comm_displs[k]+i;
      nitems_1 = rcv_nitem_ext_tab[cursor];
      int_comm_index = &rcv_int_comm_index[idx];
      idx     += nitems_1;

      /* External communicator k->color_out */
      color_in    = rcv_color_in_tab[cursor];
      color_out   = rcv_color_out_tab[cursor];

      assert( color_in==k );
      assert( color_in!=color_out );

      if ( color_in < color_out ) {
        /* New point */
        for ( j=0; j<nitems_1; ++j ) {
          if ( !point_1[intvalues_1[int_comm_index[j]]].tmp ) {
            point_1[intvalues_1[int_comm_index[j]]].tmp = ++np;
          }
        }
      }
      else {
        /* Get the point index from the neighbouring proc */

        /* Find the matching external communicator over proc color_out (->k) */
        idx_2 =  int_comm_index_displs[color_out];
        for ( j=0; j<rcv_next_node_comm[color_out]; ++j ) {
          cursor   = ext_comm_displs[color_out]+j;
          nitems_2 = rcv_nitem_ext_tab[cursor];
          int_comm_index_2 = &rcv_int_comm_index[idx_2];
          idx_2   += nitems_2;

          if ( rcv_color_in_tab[cursor]  == color_out &&
               rcv_color_out_tab[cursor] == color_in  )  break;
        }
        assert(j<rcv_next_node_comm[color_out]);
        assert(nitems_1 == nitems_2);


        /* Points and internal communicators for proc \a color_out */
        point_2     = &rcv_point[point_displs[color_out]];

        cursor      = intval_displs[color_out];
        intvalues_2 = &rcv_intvalues[cursor];

        /* Update point indices (stored in the tmp field) */
        for ( j=0; j<nitems_1; ++j ) {
          point_1[intvalues_1[int_comm_index[j]]].tmp =
            point_2[intvalues_2[int_comm_index_2[j]]].tmp;
        }
      }
    }

    /* Create new indices for the points that haven't been seen. */
    for ( i=1; i<=rcv_np[k]; ++i ) {
      if ( !point_1[i].tmp ) {
        point_1[i].tmp = ++np;
      }
    }
  }

  /** Tetra + xTetra */
  ne_tot       = tetra_displs[nprocs-1]+rcv_ne[nprocs-1];
  xt_tot       = xtetra_displs[nprocs-1]+rcv_xt[nprocs-1];

  mesh->nemax  = mesh->ne = ne_tot;
  mesh->nenil  = 0;
  mesh->xtmax  = mesh->xt = xt_tot;

  MMG5_ADD_MEM(mesh,(mesh->nemax+1)*sizeof(MMG5_Tetra),"tetra",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_ADD_MEM(mesh,(mesh->xtmax+1)*sizeof(MMG5_xTetra),"xtetra",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(mesh->xtetra,mesh->xtmax+1,MMG5_xTetra,return 0);
  MMG5_SAFE_CALLOC(mesh->tetra,mesh->nemax+1,MMG5_Tetra,return 0);

  ne = idx = 0;
  for ( k=0; k<nprocs; ++k ) {
    xtetra      = &rcv_xtetra[xtetra_displs[k]];
    tetra       = &rcv_tetra[tetra_displs[k]];
    point_1     = &rcv_point[point_displs[k]];

    for ( i=1; i<=rcv_ne[k]; ++i ) {
      pt = &mesh->tetra[tetra_displs[k]+i];

      for ( j=0; j<4; ++j ) {
        tetra[i].v[j] = point_1[tetra[i].v[j]].tmp;
      }

      memcpy(pt,&tetra[i],sizeof(MMG5_Tetra));

      if ( tetra[i].xt ) {
        nnpar = 0;
        pxt = &xtetra[tetra[i].xt];
        memcpy(&mesh->xtetra[++ne],pxt,sizeof(MMG5_xTetra));
        pt->xt = ne;
      }
    }
  }
  mesh->xt=ne;

  /** Points and solutions */
  mesh->np = met->np = np;
  mesh->npmax = met->npmax = mesh->np;
  mesh->npnil = 0;
  MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(MMG5_Point),"merge point",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(mesh->point,mesh->npmax+1,MMG5_Point,return 0);

  if ( rcv_met ) {
    MMG5_ADD_MEM(mesh,(met->npmax+1)*met->size*sizeof(double),"merge met",
                  fprintf(stderr,"  Exit program.\n");
                  return 0);
    MMG5_SAFE_CALLOC(met->m,(met->npmax+1)*met->size,double,return 0);
  }

  for ( i=1; i<=mesh->np; ++i ) mesh->point[i].tag = MG_NUL;

  np = 0;
  for ( k=0; k<nprocs; ++k ) {
    point_1     = &rcv_point[point_displs[k]];
    if ( rcv_met )
      met_1       = &rcv_met[met_displs[k]];

    for ( i=1; i<=rcv_np[k]; ++i ) {
      idx = point_1[i].tmp;
      assert(idx);
      ppt = &mesh->point[idx];

      if ( MG_VOK(ppt) ) {
        point_1[i].tmp = 0;
        continue;
      }

      memcpy(ppt,&point_1[i],sizeof(MMG5_Point));
      ppt->tmp = 0;

      if ( rcv_met )
        memcpy( &met->m[idx*met->size],
                &met_1[i*met->size],
                met->size*sizeof(double) );

      if ( point_1[i].xp ) ++np;
    }
  }

  /** xPoints */
  mesh->xpmax = mesh->xp = np;
  MMG5_ADD_MEM(mesh,(mesh->xpmax+1)*sizeof(MMG5_xPoint),"merge xPoint",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(mesh->xpoint,mesh->xpmax+1,MMG5_xPoint,return 0);
  np = 0;
  for ( k=0; k<nprocs; ++k ) {
    point_1     = &rcv_point[point_displs[k]];
    xpoint      = &rcv_xpoint[xpoint_displs[k]];

    for ( i=1; i<=rcv_np[k]; ++i ) {
      idx = point_1[i].tmp;

      if ( !idx ) continue;
      ppt = &mesh->point[idx];

      if ( !point_1[i].xp ) continue;

      pxp = &mesh->xpoint[++np];
      memcpy(pxp, &xpoint[point_1[i].xp],sizeof(MMG5_xPoint));
      ppt->xp = np;
    }
  }

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if fail, 1 if successc(on all procs)
 *
 *  merge all meshes to a single mesh in P0's parmesh
 */
int PMMG_merge_parmesh( PMMG_pParMesh parmesh ) {
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_pPoint    rcv_point,ppt;
  MMG5_pxPoint   rcv_xpoint;
  MMG5_pTetra    rcv_tetra;
  MMG5_pxTetra   rcv_xtetra;
  PMMG_pInt_comm int_node_comm;
  double         *rcv_met;
  size_t         available;
  int            *rcv_np,np_tot,*rcv_ne,ne_tot,*rcv_xp,xp_tot,*rcv_xt,xt_tot;
  int            *rcv_nmet,nmet_tot;
  int            *point_displs,*xpoint_displs,*tetra_displs,*xtetra_displs;
  int            *met_displs,*intval_displs,*ext_comm_displs;
  int            *int_comm_index_displs;
  int            *rcv_intvalues,nitem_inc_tot;
  int            *rcv_int_comm_index,nitem_icidx_tot,ext_comm_displs_tot;
  int            *rcv_node2int_node_comm_index1,*rcv_node2int_node_comm_index2;
  int            *rcv_nitem_int_node_comm,*rcv_next_node_comm,rcv_isMet;
  int            *rcv_color_in_tab,*rcv_color_out_tab,*rcv_nitem_ext_tab;
  int            k,idx,ier,ieresult;

  ier = 1;
  grp = &parmesh->listgrp[0];

  assert ( parmesh->ngrp <= 1 );

  //DEBUGGING:
  // saveGrpsToMeshes( PMMG_pGrp listgrp, 0, parmesh->myrank, "mesh" );

  /** Step 1: Allocate internal communicator buffer and fill it: the
   *  intvalues array contains the indices of the matching nodes on the proc. */

  /* Give all the memory to the communicators */
  available = parmesh->memGloMax;
  if ( grp ) {
    grp->mesh->memMax = grp->mesh->memCur;
    available -= grp->mesh->memMax;
  }

  assert ( available >= 0 );
  parmesh->memMax += available;

  /* Internal comm allocation */
  int_node_comm = parmesh->int_node_comm;

  PMMG_CALLOC(parmesh,int_node_comm->intvalues,int_node_comm->nitem,int,"intvalues",ier=0);
  MPI_CHECK( MPI_Allreduce(&ier,&ieresult,1,MPI_INT,MPI_MIN,parmesh->comm),ieresult=0);

  if ( !ieresult ) return 0;

#ifndef NDEBUG
  if ( grp )
    assert(int_node_comm->nitem == grp->nitem_int_node_comm);
#endif

  if ( grp ) {
    for ( k=0; k<grp->nitem_int_node_comm; ++k ) {
      idx = grp->node2int_node_comm_index2[k];
      int_node_comm->intvalues[idx] = grp->node2int_node_comm_index1[k];
    }
  }

  /** Step 2: Procs send their parmeshes to Proc 0 and Proc 0 recieve the data */
  ier = PMMG_gather_parmesh(parmesh,&rcv_point,&np_tot,&rcv_xpoint,&xp_tot,
                            &rcv_tetra,&ne_tot,&rcv_xtetra,&xt_tot,
                            &rcv_met,&nmet_tot,&rcv_isMet,&rcv_intvalues,
                            &nitem_inc_tot,
                            &rcv_nitem_ext_tab,&ext_comm_displs_tot,
                            &rcv_color_in_tab,&rcv_color_out_tab,
                            &rcv_node2int_node_comm_index1,
                            &rcv_node2int_node_comm_index2,&point_displs,
                            &xpoint_displs,&tetra_displs,&xtetra_displs,
                            &met_displs,&intval_displs,&ext_comm_displs,
                            &int_comm_index_displs,&rcv_np,&rcv_xp,&rcv_ne,
                            &rcv_xt,&rcv_nmet,&rcv_int_comm_index,&nitem_icidx_tot,
                            &rcv_next_node_comm,&rcv_nitem_int_node_comm);


  if ( ier ) {
    /** Step 3: Proc 0 merges the meshes: We travel through the external
     * communicators to recover the numbering of the points shared with a lower
     * proc. The other points are concatenated with the proc 0. */

    /* Give all the memory to mesh0 */
    parmesh->memMax = parmesh->memCur;
    available = parmesh->memGloMax - parmesh->memMax;
    if ( grp ) {
      available -= grp->mesh->memMax;
      assert ( available >= 0 );

      grp->mesh->memMax +=  available;
    }

    ier = PMMG_mergeParmesh_rcvParMeshes(parmesh,rcv_point,rcv_xpoint,rcv_tetra,
                                         rcv_xtetra,rcv_met,rcv_intvalues,rcv_nitem_ext_tab,
                                         rcv_color_in_tab,rcv_color_out_tab,point_displs,
                                         xpoint_displs,tetra_displs,xtetra_displs,
                                         met_displs,intval_displs,ext_comm_displs,
                                         int_comm_index_displs,rcv_np,rcv_ne,rcv_xt,
                                         rcv_int_comm_index,rcv_next_node_comm);
  }

  MPI_CHECK( MPI_Allreduce(&ier,&ieresult,1,MPI_INT,MPI_MIN,parmesh->comm),ieresult=0);

  /** Step 4: Free memory */
  /* 1: Mesh data */
  if ( parmesh->myrank == parmesh->info.root ) {
    PMMG_DEL_MEM(parmesh,rcv_np,int,"rcv_np");
    PMMG_DEL_MEM(parmesh,rcv_point,MMG5_Point,"rcv_pt");
    PMMG_DEL_MEM(parmesh,point_displs,int,"pt_displs");
    PMMG_DEL_MEM(parmesh,rcv_xp,int,"rcv_xp");
    PMMG_DEL_MEM(parmesh,rcv_xpoint,MMG5_xPoint,"rcv_xpt");
    PMMG_DEL_MEM(parmesh,xpoint_displs,int,"xp_displs");

    PMMG_DEL_MEM(parmesh,rcv_ne,int,"rcv_ne");
    PMMG_DEL_MEM(parmesh,rcv_tetra,MMG5_Tetra,"rcv_tetra");
    PMMG_DEL_MEM(parmesh,tetra_displs,int,"ne_displs");

    PMMG_DEL_MEM(parmesh,rcv_xt,int,"rcv_xt");
    PMMG_DEL_MEM(parmesh,rcv_xtetra,MMG5_xTetra,"rcv_xtetra");
    PMMG_DEL_MEM(parmesh,xtetra_displs,int,"xt_displs");

    if ( rcv_isMet ) {
      PMMG_DEL_MEM(parmesh,rcv_nmet,int,"rcv_nmet");
      PMMG_DEL_MEM(parmesh,rcv_met,double,"rcv_met");
      PMMG_DEL_MEM(parmesh,met_displs,int,"met_displs");
    }

    parmesh->ngrp = 1;
  }
  else parmesh->ngrp = 0;


  /* 2: communicators data */
  PMMG_DEL_MEM(parmesh,rcv_int_comm_index,int,"rcv_ic_idx");
  PMMG_DEL_MEM(parmesh,rcv_intvalues,int,"rcv_intvalues");
  PMMG_DEL_MEM(parmesh,intval_displs,int,"intval_displs");

  PMMG_DEL_MEM(parmesh,rcv_next_node_comm,int,"rcv_next_nc");
  PMMG_DEL_MEM(parmesh,rcv_nitem_ext_tab,int,"rcv_nitem_ext");
  PMMG_DEL_MEM(parmesh,rcv_color_in_tab,int,"rcv_color_in");
  PMMG_DEL_MEM(parmesh,rcv_color_out_tab,int,"rcv_colour_out");
  PMMG_DEL_MEM(parmesh,ext_comm_displs,int,"ext_comm_displs");

  PMMG_DEL_MEM(parmesh,rcv_nitem_int_node_comm,int,"rcv_nitem_inc");
  PMMG_DEL_MEM(parmesh,rcv_node2int_node_comm_index1,int,"rcv_inc_idx1");
  PMMG_DEL_MEM(parmesh,rcv_node2int_node_comm_index2,int,"rcv_inc_idx2");
  PMMG_DEL_MEM(parmesh,int_comm_index_displs,int,"icidx_displs");

  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,int,"intval");

  PMMG_parmesh_Free_Comm(parmesh);

  /** Step 5: Update tag on points, tetra */
  if ( ieresult > 0 ) {
    ieresult = PMMG_updateTag(parmesh);

    /** Step 6: In nosurf mode, the updateTag function has added nosurf + required
     * tag to the boundary mesh : remove its */
    mesh = parmesh->listgrp[0].mesh;

    if ( mesh->info.nosurf ) {

      MMG3D_unset_reqBoundaries(mesh);

      for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) )  continue;

        if ( ppt->tag & MG_NOSURF ) {
          ppt->tag &= ~MG_NOSURF;
          ppt->tag &= ~MG_REQ;
        }
      }
    }
  }

  return ieresult;
}

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
 * \param mesh0 pointer toward a MMG5 mesh structure.
 * \param met0 pointer toward a MMG5 solution structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Merge the points shared between groups (and listed inside the internal
 * communicator) into the mesh0 mesh.
 *
 */
static inline
int PMMG_mergeGrps_interfacePoints(PMMG_pParMesh parmesh,MMG5_pMesh mesh0,
                                   MMG5_pSol met0) {
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_pSol      met;
  MMG5_pPoint    ppt0,ppt;
  MMG5_pxPoint   pxp0,pxp;
  int            nitem_int_node_comm;
  int            *intvalues;
  int            *node2int_node_comm_index1,*node2int_node_comm_index2;
  int            size,poi_id_int,poi_id_glo,imsh,k,ip;

  grp       = parmesh->listgrp;
  intvalues = parmesh->int_node_comm->intvalues;

  /** Use the tmp field of points in meshes to remember the id in the merged mesh
   * of points that have already been added to the merged mesh or 0 if they
   * haven't been merged yet */
  for ( imsh = 0; imsh < parmesh->ngrp; ++imsh )
    for ( k = 1; k <= grp[imsh].mesh->np; ++k )
      grp[imsh].mesh->point[k].tmp = 0;

  /** Step 1: store the indices of the interface entities of mesh0 into the
   * internal communicators */
  for ( k = 0; k < grp[0].nitem_int_node_comm; ++k ) {
    poi_id_int = grp[0].node2int_node_comm_index1[k];
    poi_id_glo = grp[0].node2int_node_comm_index2[k];
    assert(   ( 0 <= poi_id_glo )
           && ( poi_id_glo < parmesh->int_node_comm->nitem )
           && "check intvalues indices" );
    intvalues[ poi_id_glo ] = poi_id_int;
  }

  /** Step 2: add points referenced by the rest of the groups meshes' internal
   * communicators to the merging mesh (mesh0) and in intvalues */
  for ( imsh=1; imsh<parmesh->ngrp; ++imsh ) {
    mesh                      = grp[imsh].mesh;
    met                       = grp[imsh].met;
    nitem_int_node_comm       = grp[imsh].nitem_int_node_comm;
    node2int_node_comm_index1 = grp[imsh].node2int_node_comm_index1;
    node2int_node_comm_index2 = grp[imsh].node2int_node_comm_index2;

    for ( k=0; k<nitem_int_node_comm; ++k ) {
      poi_id_glo = node2int_node_comm_index2[k];
      assert(   ( 0 <= poi_id_glo )
             && ( poi_id_glo < parmesh->int_node_comm->nitem )
             && "check intvalues indices"  );

      // location in currently working mesh where point data actually are
      assert ( node2int_node_comm_index1[ k ] <= mesh->np );
      ppt = &mesh->point[ node2int_node_comm_index1[ k ] ];

      if ( !MG_VOK(ppt) ) continue;

      /* point ppt is not found in the merged mesh. add it */
      if ( !intvalues[ poi_id_glo ] ) {
        ip = _MMG3D_newPt(mesh0,ppt->c,ppt->tag);
        if ( !ip ) {
          /* reallocation of point table */
          _MMG5_POINT_REALLOC(mesh0,met0,ip,mesh0->gap,
                              printf("  ## Error: unable to merge group points\n");
                              _MMG5_INCREASE_MEM_MESSAGE();
                              return 0;,
                              ppt->c,ppt->tag,0);

        }
        assert( (ip < grp[ 0 ].mesh->npmax) && "run out of points" );

        //NIKOS TODO: Remember location in merged mesh in mesh, mesh0 and
        //communicator tmp fields
        ppt->tmp = ip;
        intvalues[ poi_id_glo ] = ip;

        /* Add xpoint if needed */
        if ( ppt->xp ) {
          pxp  = &mesh->xpoint[ppt->xp];
          ppt0 = &mesh0->point[ip];
          pxp0 = &mesh0->xpoint[ppt0->xp];
          assert( (ppt0->xp <= mesh0->xpmax) && "increase xpoints" );
          assert( (ppt0->xp > 0) && "negative xpoints?!?!?!?!?" );
          memcpy(pxp0,pxp,sizeof(MMG5_xPoint));
        }
        if ( met0->m ) {
          size = met0->size;
          assert(met->m);
          assert( ((size==1||size==6) && size==met->size) && "met size issues" );
          memcpy(&met0->m[size*ip],&met->m[size*k],size*sizeof(double));
        }
      } else {
        /* point already exists in merged mesh. update his tmp field to point to
         * its mesh0 index */
        ppt->tmp = intvalues[ poi_id_glo ];
      }
    }
  }
  return 1;
}

/**
 * \param parmesh pointer toward a parmmg parmesh mesh structure.
 * \param mesh0 pointer toward a MMG5 mesh structure.
 * \param met0 pointer toward a MMG5 solution structure.
 * \param imsh index in listgrp of the group that we want to merge into mesh0
 *
 * \return 0 if fail, 1 otherwise
 *
 * Merge the non internal nodes of the lisgrp[imsh].mesh mesh into the mesh0
 * mesh.
 *
 */
static inline
int PMMG_mergeGrps_internalPoints(PMMG_pParMesh parmesh,MMG5_pMesh mesh0,
                                  MMG5_pSol met0,int imsh) {
  MMG5_pMesh     mesh;
  MMG5_pSol      met;
  MMG5_pPoint    ppt0,ppt;
  MMG5_pxPoint   pxp0,pxp;
  int            np,ip,k;

  np = mesh0->np;

  /** Loop over points and add the ones that are not already in the merged
   * mesh (mesh0) */
  mesh = parmesh->listgrp[imsh].mesh;
  met   = parmesh->listgrp[imsh].met;

  for ( k=1; k<=mesh->np; k++ ) {
    ppt = &mesh->point[k];
    if ( !MG_VOK(ppt) ) continue;
    if ( ppt->tmp )     continue;

    ip = _MMG3D_newPt(mesh0,ppt->c,ppt->tag);
    if ( !ip ) {
      /* reallocation of point table */
      _MMG5_POINT_REALLOC(mesh0,met0,ip,mesh0->gap,
                          printf("  ## Error: unable to merge group points\n");
                          _MMG5_INCREASE_MEM_MESSAGE();
                          return 0;,
                          ppt->c,ppt->tag,0);
    }
    ppt->tmp = ip;

    /* Add xpoint if needed */
    if ( ppt->xp ) {
#warning add only the "true" xpoints (not those linked to the interfaces)
      pxp  = &mesh->xpoint[ppt->xp];
      ppt0 = &mesh0->point[ip];
      pxp0 = &mesh0->xpoint[ppt0->xp];
      memcpy(pxp0,pxp,sizeof(MMG5_xPoint));
    }
    if ( met0->m ) {
      assert(met->m);
      memcpy(&met0->m[met0->size*ip],&met->m[met->size*k],met->size*sizeof(double));
    }
  }
  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param mesh0 pointer toward a MMG5 mesh structure.
 * \param imsh index in listgrp of the group that we want to merge
 *
 * \return 0 if fail, 1 otherwise
 *
 * Merge the interfaces tetra of the lisgrp[imsh].mesh mesh into the mesh0
 * mesh. An interface tetra is a tetra that has at least 1 interface triangle.
 *
 */
static inline
int PMMG_mergeGrps_interfaceTetra(PMMG_pParMesh parmesh,MMG5_pMesh mesh0,int imsh) {
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_pTetra    pt0,pt;
  MMG5_pxTetra   pxt0,pxt;
  int            nitem_int_face_comm;
  int            *intvalues;
  int            *face2int_face_comm_index1,*face2int_face_comm_index2;
  int            face_id_glo;
  int            k,iel,ie,ifac,i;

  grp       = parmesh->listgrp;
  intvalues = parmesh->int_face_comm->intvalues;

  /** Add the interfaces tetra to the mesh0 mesh */
  mesh                      = grp[imsh].mesh;
  nitem_int_face_comm       = grp[imsh].nitem_int_face_comm;
  face2int_face_comm_index1 = grp[imsh].face2int_face_comm_index1;
  face2int_face_comm_index2 = grp[imsh].face2int_face_comm_index2;

  ++mesh->base;
  for ( k=0; k<nitem_int_face_comm; ++k ) {
    face_id_glo = face2int_face_comm_index2[k];
    assert( 0 <= face_id_glo );
    assert(face_id_glo < (parmesh->int_face_comm->nitem));

    /* Index of the interface tetra in the mesh */
    iel  = face2int_face_comm_index1[k]/4;
    ifac = face2int_face_comm_index1[k]%4;
    assert ( iel && iel<=mesh->ne );
    pt       = &mesh->tetra[iel];

    assert ( MG_EOK(pt) );

    /* If the tetra has already been added to mesh0 */
    if ( pt->base == mesh->base ) {
      /* Store the interface face if it has not been seen from another group */
      assert( pt->flag );

      if ( !intvalues[face_id_glo] )
        intvalues[face_id_glo] = 4*pt->flag+ifac;
      continue;
    }

    /* Else */
    pt->base = mesh->base;

    /* Add the tetra iel to mesh0 */
    ie = _MMG3D_newElt(mesh0);
    if ( !ie ) {
      _MMG5_TETRA_REALLOC(mesh0,ie,mesh0->gap,
                          fprintf(stderr,"  ## Error: unable to merge group elts.\n");
                          _MMG5_INCREASE_MEM_MESSAGE();
                          return 0;,
                          0);
    }
    pt0 = &mesh0->tetra[ie];

    for ( i=0; i<4; ++i ) pt0->v[i] = mesh->point[pt->v[i]].tmp;
    pt0->ref  = pt->ref;
    pt->flag  = ie;

    /* Store the interface face if it has not been seen from another group */
    if ( !intvalues[face_id_glo] )
      intvalues[face_id_glo] = 4*ie+ifac;

    /** Add xtetra if needed */
    if ( pt->xt ) {
#warning add only the "true" xtetras (not those linked to the interfaces)
      pxt = &mesh->xtetra[pt->xt];
      mesh0->xt++;
      if ( mesh0->xt > mesh0->xtmax ) {
        PMMG_RECALLOC(mesh0, mesh0->xtetra, 1.2 * mesh0->xtmax + 1,
                      mesh0->xtmax + 1, MMG5_xTetra,
                      "larger xtetra table", mesh0->xt--; goto fail_ncomm);
        mesh0->xtmax = 1.2 * mesh0->xtmax;
      }
      pt0->xt = mesh0->xt;
      pxt0 = &mesh0->xtetra[pt0->xt];
      memcpy(pxt0,pxt,sizeof(MMG5_xTetra));
    }
  }

  return 1;
}

/**
 * \param parmesh pointer toward a PMMG parmesh structure.
 * \param mesh0 pointer toward a MMG5 mesh structure.
 * \param imsh index in listgrp of the group that we want to merge into mesh0
 *
 * \return 0 if fail, 1 otherwise
 *
 * Merge the tetra of the lisgrp[imsh].mesh mesh into the mesh0 mesh.
 *
 */
static inline
int PMMG_mergeGrps_internalTetra(PMMG_pParMesh parmesh,MMG5_pMesh mesh0,int imsh) {
  MMG5_pMesh     mesh;
  MMG5_pTetra    pt0,pt;
  MMG5_pxTetra   pxt0,pxt;
  int            k,ie,i;

  /** Add current meshs' tetras to the merged mesh (mesh0) */
  mesh = parmesh->listgrp[imsh].mesh;
  for ( k=1; k<=mesh->ne; k++ ) {
    pt  = &mesh->tetra[k];

    if ( (!MG_EOK(pt)) || pt->base==mesh->base ) continue;

    pt->base = mesh->base;

    ie  = _MMG3D_newElt(mesh0);
    if ( !ie ) {
      _MMG5_TETRA_REALLOC(mesh0,ie,mesh0->gap,
                          fprintf(stderr,"  ## Error: unable to merge group elts.\n");
                          _MMG5_INCREASE_MEM_MESSAGE();
                          return 0;,
                          0);
    }
    pt0 = &mesh0->tetra[ie];

    for ( i=0; i<4; ++i ) pt0->v[i] = mesh->point[pt->v[i]].tmp;
    pt0->ref = pt->ref;

#warning need to update the adjacents?

    /** Add xtetra if needed */
    if ( pt->xt ) {
      pxt = &mesh->xtetra[pt->xt];
      mesh0->xt++;
      if ( mesh0->xt > mesh0->xtmax ) {
        PMMG_RECALLOC(mesh0, mesh0->xtetra, 1.2 * mesh0->xtmax + 1,
                      mesh0->xtmax + 1, MMG5_xTetra,
                      "larger xtetra table", mesh0->xt--; goto fail_ncomm);
        mesh0->xtmax = 1.2 * mesh0->xtmax;
      }
      pt0->xt = mesh0->xt;
      pxt0 = &mesh0->xtetra[pt0->xt];
      memcpy(pxt0,pxt,sizeof(MMG5_xTetra));
    }
  }

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Update the node communicators when merging the groups of listgroups into 1
 * group.
 *
 */
static inline
int PMMG_mergeGrps_nodeCommunicators(PMMG_pParMesh parmesh) {
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh0;
  MMG5_pSol      met0;
  MMG5_pPoint    ppt;
  PMMG_pext_comm ext_node_comm;
  PMMG_pint_comm int_node_comm;
  int            nitem_int_node_comm0,*intvalues;
  int           *node2int_node_comm0_index1;
  int           *node2int_node_comm0_index2;
  int            poi_id_int,poi_id_glo,idx,k,i,ip;
  int            new_nitem_int_node_comm;

  grp  = parmesh->listgrp;

  int_node_comm              = parmesh->int_node_comm;
  intvalues                  = int_node_comm->intvalues;
  mesh0                      = grp[0].mesh;
  met0                       = grp[0].met;
  nitem_int_node_comm0       = grp[0].nitem_int_node_comm;
  node2int_node_comm0_index1 = grp[0].node2int_node_comm_index1;
  node2int_node_comm0_index2 = grp[0].node2int_node_comm_index2;

  /** Reset the tmp field of the point: it will be used to store the position of
   * a point in the internal communicator */
  for ( k=1; k<=mesh0->np; k++ )
    mesh0->point[k].tmp = 0;

  /** Travel through the external communicators and udpate all the communicators */
  poi_id_int = 0;
  for ( k=0; k<parmesh->next_node_comm; ++k ) {

    /* currently working external communicator */
    ext_node_comm = &parmesh->ext_node_comm[k];

    poi_id_glo = 0;
    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx = ext_node_comm->int_comm_index[i];
      assert( ( 0 <= idx ) && ( idx < parmesh->int_node_comm->nitem ) && "check intvalues indices" );
      ip  = intvalues[idx];
      ppt = &mesh0->point[ip];

      /* New point in the internal communicator */
      if ( !ppt->tmp ) {
        if ( poi_id_int == grp[0].nitem_int_node_comm ) {
          new_nitem_int_node_comm = (int)(1.2*grp[0].nitem_int_node_comm);
          PMMG_REALLOC(parmesh,grp[0].node2int_node_comm_index1,
                       new_nitem_int_node_comm,
                       grp[0].nitem_int_node_comm,int,
                       "(mergeGrps) node2int_node_comm_index1",
                       grp[0].nitem_int_node_comm = new_nitem_int_node_comm;
                       return 0);
          PMMG_REALLOC(parmesh,grp[0].node2int_node_comm_index2,
                       new_nitem_int_node_comm,
                       grp[0].nitem_int_node_comm,int,
                       "(mergeGrps) node2int_node_comm_index2",
                       grp[0].nitem_int_node_comm = new_nitem_int_node_comm;
                       return 0);
          grp[0].nitem_int_node_comm = new_nitem_int_node_comm;
          node2int_node_comm0_index1 = grp[0].node2int_node_comm_index1;
          node2int_node_comm0_index2 = grp[0].node2int_node_comm_index2;
        }
        node2int_node_comm0_index1[ poi_id_int ]     = ip;
        node2int_node_comm0_index2[ poi_id_int ]     = poi_id_int;

        ext_node_comm->int_comm_index[ poi_id_glo++ ]= poi_id_int;
        ppt->tmp                                     = poi_id_int + k*mesh0->np;
        poi_id_int++;
      } else {
        /* The point has already a position in the internal comm */
        if ( ppt->tmp < k * mesh0->np ) {
          /* The point has been stored in the internal comm by another external
           * comm: update its position in our external comm */
          ext_node_comm->int_comm_index[ poi_id_glo++] = ppt->tmp%mesh0->np;
        }
      }
    }
    assert( (poi_id_glo != 0)  && "empty communicator?????" );
    assert(poi_id_glo==ext_node_comm->nitem);

    // To remove if the assertion is always true.
    // PMMG_REALLOC(parmesh,ext_node_comm->int_comm_index, poi_id_glo,
    //              ext_node_comm->nitem,int,"(mergeGrps) ext_node_comm",
    //              ext_node_comm->nitem = poi_id_glo;return 0);
    // ext_node_comm->nitem = poi_id_glo;
  }

  PMMG_REALLOC(parmesh,grp[0].node2int_node_comm_index1,poi_id_int,
               grp[0].nitem_int_node_comm,int,
               "(mergeGrps) node2int_node_comm_index1",return 0);
  PMMG_REALLOC(parmesh,grp[0].node2int_node_comm_index2,poi_id_int,
               grp[0].nitem_int_node_comm,int,
               "(mergeGrps) node2int_node_comm_index2",return 0);
  grp[0].nitem_int_node_comm = poi_id_int;
  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,int_node_comm->nitem,int,
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
  MMG5_pMesh     mesh0;
  MMG5_pSol      met0;
  PMMG_pext_comm ext_face_comm;
  PMMG_pint_comm int_face_comm;
  int            nitem_int_face_comm0,*intvalues;
  int           *face2int_face_comm0_index1;
  int           *face2int_face_comm0_index2;
  int            face_id_int,idx,k,i,iel;
  int            new_nitem_int_face_comm;

  grp  = parmesh->listgrp;

  int_face_comm              = parmesh->int_face_comm;
  intvalues                  = int_face_comm->intvalues;
  mesh0                      = grp[0].mesh;
  met0                       = grp[0].met;
  nitem_int_face_comm0       = grp[0].nitem_int_face_comm;
  face2int_face_comm0_index1 = grp[0].face2int_face_comm_index1;
  face2int_face_comm0_index2 = grp[0].face2int_face_comm_index2;

  /** Travel through the external communicators and udpate all the communicators */
  face_id_int = 0;
  for ( k=0; k<parmesh->next_face_comm; ++k ) {

    /* currently working external communicator */
    ext_face_comm = &parmesh->ext_face_comm[k];

    for ( i=0; i<ext_face_comm->nitem; ++i ) {
      idx = ext_face_comm->int_comm_index[i];
      iel = intvalues[idx];
      assert(iel);

      /* Add this face to the face communicators */
      if ( face_id_int == grp[0].nitem_int_face_comm ) {
        new_nitem_int_face_comm = (int)(1.2*grp[0].nitem_int_face_comm);
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
  PMMG_DEL_MEM(parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,
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

  if ( !PMMG_mergeGrps_nodeCommunicators(parmesh) ) return 0;

  if ( !PMMG_mergeGrps_faceCommunicators(parmesh) ) return 0;

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
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh0,mesh;
  MMG5_pSol      met0,met;
  PMMG_pint_comm int_node_comm,int_face_comm;
  int            *face2int_face_comm_index1,*face2int_face_comm_index2;
  int            imsh,k,iel;

  grp  = parmesh->listgrp;

  /** Free the adjacency array: a possible improvement is to update it */
  mesh0 = grp[0].mesh;
  met0  = grp[0].met;

  if ( mesh0->adja )
    PMMG_DEL_MEM(mesh0, mesh0->adja, 4*mesh0->nemax+5, int, "adjacency table" );

  if ( parmesh->ngrp == 1 ) return 1;

  /** Use the internal communicators to store the interface entities indices */
  int_node_comm              = parmesh->int_node_comm;
  PMMG_CALLOC(parmesh,int_node_comm->intvalues,int_node_comm->nitem,int,
              "node communicator",return 0);

  int_face_comm              = parmesh->int_face_comm;
  PMMG_CALLOC(parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,
              "face communicator",goto fail_ncomm);

  //DEBUGGING:
  //saveGrpsToMeshes(parmesh->listgrp,parmesh->ngrp,parmesh->myrank,"BeforeMergeGrp");

   /** Step 1: Merge interface points from all meshes into mesh0->points */
  if ( !PMMG_mergeGrps_interfacePoints(parmesh,mesh0,met0) ) goto fail_comms;

  for ( imsh=1; imsh<parmesh->ngrp; ++imsh ) {
    /** Step 2: Merge internal points of the mesh mesh into the mesh0 mesh */
    if ( !PMMG_mergeGrps_internalPoints(parmesh,mesh0,met0,imsh) )
      goto fail_comms;

    /** Step 3: Merge interfaces tetras of the imsh mesh into the mesh0 mesh */
    /* 1) Store the indices of the interface faces of mesh0 into the internal
     * face communicator */
    face2int_face_comm_index1 = parmesh->listgrp[0].face2int_face_comm_index1;
    face2int_face_comm_index2 = parmesh->listgrp[0].face2int_face_comm_index2;
    for ( k=0; k<grp->nitem_int_face_comm; ++k ) {
      iel = face2int_face_comm_index1[k];
      assert(   ( 0 <= face2int_face_comm_index2[k] )
             && ( face2int_face_comm_index2[k] < parmesh->int_node_comm->nitem )
             && "check intvalues indices" );
      parmesh->int_face_comm->intvalues[face2int_face_comm_index2[k]] = iel;
    }

    /* 2) Add the interfaces tetra of the imsh mesh to the mesh0 mesh */
    if ( !PMMG_mergeGrps_interfaceTetra(parmesh,mesh0,imsh) ) goto fail_comms;

    /** Step 4: Merge internal tetras of the mesh mesh into the mesh0 mesh */
    if ( !PMMG_mergeGrps_internalTetra(parmesh,mesh0,imsh) ) goto fail_comms;
  }
  /* add the merged groups' memory allowance to the merged groupd and free them */
  for ( imsh=1; imsh<parmesh->ngrp; ++imsh ) {
    mesh = parmesh->listgrp[imsh].mesh;
    met  = parmesh->listgrp[imsh].met;
    mesh0->memMax += mesh->memMax;
    MMG3D_Free_all(MMG5_ARG_start,
                   MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,
                   MMG5_ARG_end);
  }
  assert ( mesh0->memMax+parmesh->memMax<=parmesh->memGloMax );

  /** Step 5: Update the communicators */
  if ( !PMMG_mergeGrps_communicators(parmesh) ) goto fail_comms;

  _MMG5_SAFE_REALLOC(grp,1,PMMG_Grp,"(mergeGrps) listgrp",0);
  parmesh->listgrp = grp;

  parmesh->ngrp = 1;

  return 1;

fail_comms:
  PMMG_DEL_MEM(parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,
               "face communicator");

fail_ncomm:
  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,int_node_comm->nitem,int,
               "node communicator");

  return 0;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param rcv_point Buffer to gather points
 * \param rcv_xpoint Buffer to gather xPoints
 * \param rcv_tetra Buffer to gather tetra
 * \param rcv_xtetra Buffer to gather xtetra
 * \param rcv_met buffer to gather metric
 * \param rcv_intvalues Buffer to gather the intvalue array of the internal comm
 * \param rcv_nitem_ext_tab Buffer to gather the number of item in the ext comm
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
 * \param rcv_next_node_comm Buffer to gather the numbers of external comm
 * \param rcv_nitem_int_node_comm Buffer to gather the node2int_node_comm arrays
 * sizes
 *
 * \return 0 if fail, 1 otherwise
 *
 * Gather the parmeshes on the proc 0.
 *
 * \warning We must have 1 group per parmesh
 *
 */
static inline
int PMMG_gather_parmesh( PMMG_pParMesh parmesh,MMG5_pPoint *rcv_point,
                         MMG5_pxPoint *rcv_xpoint,MMG5_pTetra *rcv_tetra,
                         MMG5_pxTetra *rcv_xtetra,double **rcv_met,
                         int **rcv_intvalues,int **rcv_nitem_ext_tab,
                         int **rcv_color_in_tab,int **rcv_color_out_tab,
                         int **rcv_node2int_node_comm_index1,
                         int **rcv_node2int_node_comm_index2,int **point_displs,
                         int **xpoint_displs,int **tetra_displs,
                         int **xtetra_displs,int **met_displs,
                         int **intval_displs,int **ext_comm_displs,
                         int **int_comm_index_displs,int **rcv_np,int **rcv_xp,
                         int **rcv_ne,int **rcv_xt,int **rcv_nmet,
                         int **rcv_int_comm_index,int** rcv_next_node_comm,
                         int **rcv_nitem_int_node_comm ) {
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_pSol      met;
  PMMG_pint_comm int_node_comm;
  PMMG_pext_comm ext_node_comm;
  MPI_Comm       comm;
  MPI_Datatype   mpi_point,mpi_xpoint,mpi_tetra,mpi_xtetra;
  int            nmet_tot,ne_tot,nitem_int_node_comm_tot;
  int            *color_in_tab,*color_out_tab,*nitem_ext_tab;
  int            *int_comm_index,*nitems_ext_idx;
  int            ext_comm_displs_tot,nitem_ext_tot;
  int            nprocs,k,i,idx;
  int            xt_tot,np_tot,xp_tot,isMet;
  //  char           filename[32];

  nprocs        = parmesh->nprocs;
  grp           = &parmesh->listgrp[0];
  mesh          = grp->mesh;
  met           = grp->met;
  isMet         = (met->m)? 1 : 0;
  comm          = parmesh->comm;
  int_node_comm = parmesh->int_node_comm;

  /** Memory alloc */
  /* 1: Mesh data */
  _MMG5_SAFE_CALLOC( (*rcv_np)        ,nprocs,int,0);
  _MMG5_SAFE_CALLOC( (*point_displs)  ,nprocs,int,0);

  _MMG5_SAFE_CALLOC( (*rcv_xp)        ,nprocs,int,0);
  _MMG5_SAFE_CALLOC( (*xpoint_displs) ,nprocs,int,0);

  _MMG5_SAFE_CALLOC( (*rcv_ne)        ,nprocs,int,0);
  _MMG5_SAFE_CALLOC( (*tetra_displs)  ,nprocs,int,0);

  _MMG5_SAFE_CALLOC( (*rcv_xt)        ,nprocs,int,0);
  _MMG5_SAFE_CALLOC( (*xtetra_displs) ,nprocs,int,0);

  if(isMet) {
    _MMG5_SAFE_CALLOC( (*rcv_nmet)    ,nprocs,int,0);
    _MMG5_SAFE_CALLOC( (*met_displs)  ,nprocs,int,0);
  }

  /* 2: Communicators */
  _MMG5_SAFE_CALLOC( (*intval_displs)          ,nprocs,int,0);
  _MMG5_SAFE_CALLOC( (*rcv_nitem_int_node_comm),nprocs,int,0);
  _MMG5_SAFE_CALLOC( (*int_comm_index_displs)  ,nprocs,int,0);
  _MMG5_SAFE_CALLOC( (*rcv_next_node_comm)     ,nprocs,int,0);
  _MMG5_SAFE_CALLOC( (*ext_comm_displs)        ,nprocs,int,0);

  /* if ( parmesh->ddebug ) { */
  /*   sprintf(filename,"Before_Gather_proc%d.mesh",parmesh->myrank); */
  /*   _MMG3D_bdryBuild(parmesh->listgrp[0].mesh); */
  /*   MMG3D_saveMesh(mesh,filename); */
  /*   MMG3D_saveSol(mesh,met,filename); */
  /* } */

  /** Gather parmesh size infos on proc 0 */
  MPI_Gather(&mesh->np,1,MPI_INT,(*rcv_np),1,MPI_INT,0,comm);
  MPI_Gather(&mesh->ne,1,MPI_INT,(*rcv_ne),1,MPI_INT,0,comm);
  MPI_Gather(&mesh->xp,1,MPI_INT,(*rcv_xp),1,MPI_INT,0,comm);
  MPI_Gather(&mesh->xt,1,MPI_INT,(*rcv_xt),1,MPI_INT,0,comm);
  MPI_Gather(&grp->nitem_int_node_comm,1,MPI_INT,
              (*rcv_nitem_int_node_comm),1,MPI_INT,0,comm);
  MPI_Gather(&parmesh->next_node_comm,1,MPI_INT,
              (*rcv_next_node_comm),1,MPI_INT,0,comm);

  /** Gather meshes on proc 0 */
  /* Creation of MPI types for the mesh comm */
  PMMG_create_MPI_Point (mesh->point,  &mpi_point );
  PMMG_create_MPI_xPoint(mesh->xpoint, &mpi_xpoint);
  PMMG_create_MPI_Tetra (mesh->tetra,  &mpi_tetra );
  PMMG_create_MPI_xTetra(mesh->xtetra, &mpi_xtetra);

  /* Points */
  for ( k=1; k<=mesh->np; ++k ) mesh->point[k].tmp = 0;

  (*point_displs)[0] = 0;
  for ( k=1; k<nprocs; ++k ) {
    (*point_displs)[k] += (*point_displs)[k-1] + (*rcv_np)[k-1];
  }
  np_tot       = (*point_displs)[nprocs-1]+(*rcv_np)[nprocs-1];
  (*rcv_point) = (MMG5_pPoint)calloc(np_tot+1,sizeof(MMG5_Point));

  MPI_Gatherv(&mesh->point[1],mesh->np,mpi_point,&(*rcv_point)[1],
               (*rcv_np),(*point_displs),mpi_point,0,comm);

  /* xPoints */
  (*xpoint_displs)[0] = 0;
  for ( k=1; k<nprocs; ++k ) {
    (*xpoint_displs)[k] += (*xpoint_displs)[k-1] + (*rcv_xp)[k-1];
  }
  xp_tot        = (*xpoint_displs)[nprocs-1]+(*rcv_xp)[nprocs-1];
  (*rcv_xpoint) = (MMG5_pxPoint)calloc(xp_tot+1,sizeof(MMG5_xPoint));

  MPI_Gatherv(&mesh->xpoint[1],mesh->xp,mpi_xpoint,&(*rcv_xpoint)[1],
              (*rcv_xp),(*xpoint_displs),mpi_xpoint,0,comm);

  /* Tetra */
  (*tetra_displs)[0] = 0;
  for ( k=1; k<nprocs; ++k ) {
    (*tetra_displs)[k] += (*tetra_displs)[k-1] + (*rcv_ne)[k-1];
  }
  ne_tot       = (*tetra_displs)[nprocs-1]+(*rcv_ne)[nprocs-1];
  (*rcv_tetra) = (MMG5_pTetra)calloc(ne_tot+1,sizeof(MMG5_Tetra));
  MPI_Gatherv(&mesh->tetra[1],mesh->ne,mpi_tetra,&(*rcv_tetra)[1],
              (*rcv_ne),(*tetra_displs),mpi_tetra,0,comm);

  /* xTetra */
  (*xtetra_displs)[0] = 0;
  for ( k=1; k<nprocs; ++k ) {
    (*xtetra_displs)[k] += (*xtetra_displs)[k-1] + (*rcv_xt)[k-1];
  }
  xt_tot        = (*xtetra_displs)[nprocs-1]+(*rcv_xt)[nprocs-1];
  (*rcv_xtetra) = (MMG5_pxTetra)calloc(xt_tot+1,sizeof(MMG5_xTetra));

  MPI_Gatherv(&mesh->xtetra[1],mesh->xt,mpi_xtetra,&(*rcv_xtetra)[1],
              (*rcv_xt),(*xtetra_displs),mpi_xtetra,0,comm);

  /* Solutions */
  *rcv_met = NULL;
  if(isMet) {
    for ( k=0; k<nprocs; ++k ) {
      (*rcv_nmet)[k] = met->size*(*rcv_np)[k];
    }

    (*met_displs)[0] = 0;
    for ( k=1; k<nprocs; ++k ) {
      (*met_displs)[k] += (*met_displs)[k-1] + (*rcv_nmet)[k-1];
    }
    nmet_tot   = (*met_displs)[nprocs-1]+(*rcv_nmet)[nprocs-1];
    (*rcv_met) = (double*)calloc(nmet_tot+met->size,sizeof(double));

    MPI_Gatherv(&met->m[met->size],mesh->np*met->size,MPI_DOUBLE,
                &(*rcv_met)[met->size],(*rcv_nmet),(*met_displs),MPI_DOUBLE,0,
                comm);
  }

  /* Internal communicator */
  (*intval_displs)[0] = 0;
  for ( k=1; k<nprocs; ++k ) {
    (*intval_displs)[k] += (*intval_displs)[k-1] + (*rcv_nitem_int_node_comm)[k-1];
  }
  nitem_int_node_comm_tot = (*intval_displs)[nprocs-1]
    + (*rcv_nitem_int_node_comm)[nprocs-1];

  (*rcv_intvalues)                 = malloc(nitem_int_node_comm_tot*sizeof(int));
  (*rcv_node2int_node_comm_index1) = malloc(nitem_int_node_comm_tot*sizeof(int));
  (*rcv_node2int_node_comm_index2) = malloc(nitem_int_node_comm_tot*sizeof(int));

  MPI_Gatherv(int_node_comm->intvalues,int_node_comm->nitem,MPI_INT,
               (*rcv_intvalues),(*rcv_nitem_int_node_comm),(*intval_displs),MPI_INT,
              0,comm);
  MPI_Gatherv(grp->node2int_node_comm_index1,int_node_comm->nitem,MPI_INT,
              (*rcv_node2int_node_comm_index1),(*rcv_nitem_int_node_comm),
              (*intval_displs),MPI_INT,0,comm);
  MPI_Gatherv(grp->node2int_node_comm_index2,int_node_comm->nitem,MPI_INT,
              (*rcv_node2int_node_comm_index2),(*rcv_nitem_int_node_comm),
              (*intval_displs),MPI_INT,0,comm);


  /* External communicator */
  /* For each proc, we concatenate the external communicators fields into arrays */
  color_in_tab  =(int*)malloc(parmesh->next_node_comm*sizeof(int));
  color_out_tab =(int*)malloc(parmesh->next_node_comm*sizeof(int));
  nitem_ext_tab =(int*)malloc(parmesh->next_node_comm*sizeof(int));

  nitem_ext_tot = 0;
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm    = &parmesh->ext_node_comm[k];
    color_in_tab[k ] = ext_node_comm->color_in;
    color_out_tab[k] = ext_node_comm->color_out;
    nitem_ext_tab[k] = ext_node_comm->nitem;
    nitem_ext_tot += ext_node_comm->nitem;
  }

  int_comm_index = (int*)malloc(nitem_ext_tot*sizeof(int));
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
  ext_comm_displs_tot = (*ext_comm_displs)[nprocs-1] + (*rcv_next_node_comm)[nprocs-1];

  (*rcv_color_in_tab)  = (int*)malloc(ext_comm_displs_tot*sizeof(int));
  (*rcv_color_out_tab) = (int*)malloc(ext_comm_displs_tot*sizeof(int));
  (*rcv_nitem_ext_tab) = (int*)malloc(ext_comm_displs_tot*sizeof(int));

  MPI_Gatherv(color_in_tab,parmesh->next_node_comm,MPI_INT,
              (*rcv_color_in_tab),(*rcv_next_node_comm),(*ext_comm_displs),MPI_INT,
              0,comm);
  MPI_Gatherv(color_out_tab,parmesh->next_node_comm,MPI_INT,
              (*rcv_color_out_tab),(*rcv_next_node_comm),(*ext_comm_displs),MPI_INT,
              0,comm);
  MPI_Gatherv(nitem_ext_tab,parmesh->next_node_comm,MPI_INT,
              (*rcv_nitem_ext_tab),(*rcv_next_node_comm),(*ext_comm_displs),MPI_INT,
              0,comm);

#warning NIKOS: the  return 0 in this macro results in leaking memory(if ever executed): int_comm_index
  _MMG5_SAFE_CALLOC(nitems_ext_idx,nprocs,int,0);

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

  (*rcv_int_comm_index) = (int*)malloc(((*int_comm_index_displs)[nprocs-1]+
                                        nitems_ext_idx[nprocs-1])*sizeof(int));
  MPI_Gatherv(int_comm_index,nitem_ext_tot,MPI_INT,
              (*rcv_int_comm_index),nitems_ext_idx,(*int_comm_index_displs),MPI_INT,
              0,comm);

  /** Free the memory */
  /* Free MPI types */
  MPI_Type_free(&mpi_point);
  MPI_Type_free(&mpi_xpoint);
  MPI_Type_free(&mpi_tetra);
  MPI_Type_free(&mpi_xtetra);

  /* Free temporary arrays */
  _MMG5_SAFE_FREE(color_in_tab);
  _MMG5_SAFE_FREE(color_out_tab);
  _MMG5_SAFE_FREE(nitem_ext_tab);
  _MMG5_SAFE_FREE(nitems_ext_idx);
  _MMG5_SAFE_FREE(int_comm_index);

  /* Free useless groups */
  /* 1: mesh */
  _MMG5_DEL_MEM(mesh,mesh->point,(mesh->npmax+1)*sizeof(MMG5_Point));
  _MMG5_DEL_MEM(mesh,mesh->tetra,(mesh->nemax+1)*sizeof(MMG5_Tetra));
  _MMG5_DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(MMG5_xPoint));
  _MMG5_DEL_MEM(mesh,mesh->xtetra,(mesh->xtmax+1)*sizeof(MMG5_Tetra));
  if(isMet)
    _MMG5_DEL_MEM(mesh,met->m,(met->npmax+1)*met->size*sizeof(double));

  /* 2: communicators */
  _MMG5_SAFE_FREE(grp->edge2int_edge_comm_index1);
  _MMG5_SAFE_FREE(grp->edge2int_edge_comm_index2);
  _MMG5_SAFE_FREE(parmesh->int_node_comm->intvalues);

  for ( i=0; i<parmesh->next_node_comm; ++i ) {
    _MMG5_SAFE_FREE(parmesh->ext_node_comm->int_comm_index);
  }
  _MMG5_SAFE_FREE(parmesh->ext_node_comm);

  return 1;
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
 * \param rcv_node2int_node_comm_index1 Buffer that gathers node2int_node_comm_index1
 * \param rcv_node2int_node_comm_index2 Buffer that gathers node2int_node_comm_index2
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
 * \param rcv_xp Buffer that gathers the number of xPoints
 * \param rcv_ne Buffer that gathers the number of tetra
 * \param rcv_xt Buffer that gathers the number of xtetra
 * \param rcv_nmet Buffer that gathers the number of metrics
 * \param rcv_int_comm_index Buffer that gathers the internal comm sizes
 * \param rcv_next_node_comm Buffer that gathers the numbers of external comm
 * \param rcv_nitem_int_node_comm Buffer that gathers the node2int_node_comm arrays
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
                         int *rcv_node2int_node_comm_index1,
                         int *rcv_node2int_node_comm_index2,int *point_displs,
                         int *xpoint_displs,int *tetra_displs,
                         int *xtetra_displs,int *met_displs,
                         int *intval_displs,int *ext_comm_displs,
                         int *int_comm_index_displs,int *rcv_np,int *rcv_xp,
                         int *rcv_ne,int *rcv_xt,int *rcv_nmet,
                         int *rcv_int_comm_index,int* rcv_next_node_comm,
                         int *rcv_nitem_int_node_comm  ) {
  MMG5_pMesh     mesh;
  MMG5_pPoint    point_1,point_2,ppt;
  MMG5_pxPoint   xpoint,pxp;
  MMG5_pTetra    tetra,pt;
  MMG5_pxTetra   xtetra,pxt;
  MMG5_pSol      met;
  double         *met_1;
  int            *node2int_node_comm_index1  ,*node2int_node_comm_index2;
  int            *node2int_node_comm_index1_2,*node2int_node_comm_index2_2;
  int            *int_comm_index,*int_comm_index_2;
  int            *intvalues_1,*intvalues_2,nitems_1,nitems_2;
  int            nprocs,k,i,j,idx,idx_2,cursor,color_in,color_out;
  int            np,ne,ne_tot,xt_tot,nnpar,l;

  nprocs = parmesh->nprocs;
  mesh   = parmesh->listgrp[0].mesh;
  met    = parmesh->listgrp[0].met;

  if ( !parmesh->myrank ) {
    np = 0;

    /** Mesh renumbering to have the same indices at points shared by multiple
     * processors; The new indices are stored in the tmp field of the MMG5_Point
     * structure. */
    for ( k=0; k<nprocs; ++k ) {
      /* Points and internal communicators for proc \a k */
      point_1     = &rcv_point[point_displs[k]];
      cursor      = intval_displs[k];
      intvalues_1 = &rcv_intvalues[cursor];
      node2int_node_comm_index1 = &rcv_node2int_node_comm_index1[cursor];
      node2int_node_comm_index2 = &rcv_node2int_node_comm_index2[cursor];

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
          node2int_node_comm_index1_2 = &rcv_node2int_node_comm_index1[cursor];
          node2int_node_comm_index2_2 = &rcv_node2int_node_comm_index2[cursor];

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
    _MMG5_SAFE_CALLOC(mesh->xtetra,xt_tot+1,MMG5_xTetra,0);
    _MMG5_SAFE_CALLOC(mesh->tetra,mesh->nemax+1,MMG5_Tetra,0);

    mesh->xt = xt_tot;
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
          for ( l=0; l<4; ++l ) {
            if ( pxt->ftag[l] && !(pxt->ftag[l] & MG_PARBDY) ) ++nnpar;
            if ( (pxt->ftag[l] & MG_PARBDY) && (pxt->ftag[l] & MG_REQ) )
              pxt->ftag[l] &= ~MG_REQ;
          }

          if ( !nnpar ) {
            pt->xt = 0;
          }
          else {
            memcpy(&mesh->xtetra[++ne],pxt,sizeof(MMG5_xTetra));
            pt->xt = ne;
          }
        }
      }
    }
    mesh->xt=ne;

    /** Points and solutions */
    mesh->np = met->np = np;
    mesh->npmax = met->npmax = mesh->np;
    mesh->npnil = 0;
    _MMG5_SAFE_CALLOC(mesh->point,mesh->npmax+1,MMG5_Point,0);
    if ( rcv_met )
      _MMG5_SAFE_CALLOC(met->m,(met->npmax+1)*met->size,double,0);

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
    _MMG5_SAFE_CALLOC(mesh->xpoint,np+1,MMG5_xPoint,0);
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
  }

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if fail, 1 if success
 *
 *  TODO what?
 */
int PMMG_merge_parmesh( PMMG_pParMesh parmesh ) {
  PMMG_pGrp      grp;
  MMG5_pPoint    rcv_point;
  MMG5_pxPoint   rcv_xpoint;
  MMG5_pTetra    rcv_tetra;
  MMG5_pxTetra   rcv_xtetra;
  PMMG_pint_comm int_node_comm;
  double         *rcv_met;
  int            *rcv_np,*rcv_ne,*rcv_xp,*rcv_xt,*rcv_nmet;
  int            *point_displs,*xpoint_displs,*tetra_displs,*xtetra_displs;
  int            *met_displs,*intval_displs,*ext_comm_displs;
  int            *int_comm_index_displs;
  int            *rcv_intvalues,*rcv_int_comm_index;
  int            *rcv_node2int_node_comm_index1,*rcv_node2int_node_comm_index2;
  int            *rcv_nitem_int_node_comm,*rcv_next_node_comm;
  int            *rcv_color_in_tab,*rcv_color_out_tab,*rcv_nitem_ext_tab;
  int            k,idx,ier;

  ier = 1;
  grp = &parmesh->listgrp[0];

  assert ( parmesh->ngrp==1 );

  //DEBUGGING:
  // saveGrpsToMeshes( PMMG_pGrp listgrp, 0, parmesh->myrank, "mesh" );

  /** Step 1: Allocate internal communicator buffer and fill it: the
   *  intvalues array contains the indices of the matching nodes on the proc. */
  int_node_comm = parmesh->int_node_comm;

  if ( !int_node_comm->nitem ) return 1;

  _MMG5_SAFE_CALLOC(int_node_comm->intvalues,int_node_comm->nitem,int,0);
  assert(int_node_comm->nitem == grp->nitem_int_node_comm);

  for ( k=0; k<grp->nitem_int_node_comm; ++k ) {
    idx = grp->node2int_node_comm_index2[k];
    int_node_comm->intvalues[idx] = grp->node2int_node_comm_index1[k];
  }

  /** Step 2: Procs send their parmeshes to Proc 0 and Proc 0 recieve the data */
  if ( !PMMG_gather_parmesh(parmesh,&rcv_point,&rcv_xpoint,&rcv_tetra,
                            &rcv_xtetra,&rcv_met,&rcv_intvalues,&rcv_nitem_ext_tab,
                            &rcv_color_in_tab,&rcv_color_out_tab,
                            &rcv_node2int_node_comm_index1,
                            &rcv_node2int_node_comm_index2,&point_displs,
                            &xpoint_displs,&tetra_displs,&xtetra_displs,
                            &met_displs,&intval_displs,&ext_comm_displs,
                            &int_comm_index_displs,&rcv_np,&rcv_xp,&rcv_ne,
                            &rcv_xt,&rcv_nmet,&rcv_int_comm_index,
                            &rcv_next_node_comm,&rcv_nitem_int_node_comm) ) {
    ier = 0;
    goto fail;
  }

  /** Step 3: Proc 0 merges the meshes: We travel through the external
   * communicators to recover the numbering of the points shared with a lower
   * proc. The other points are concatenated with the proc 0. */
  if ( !PMMG_mergeParmesh_rcvParMeshes(parmesh,rcv_point,rcv_xpoint,rcv_tetra,
                            rcv_xtetra,rcv_met,rcv_intvalues,rcv_nitem_ext_tab,
                            rcv_color_in_tab,rcv_color_out_tab,
                            rcv_node2int_node_comm_index1,
                            rcv_node2int_node_comm_index2,point_displs,
                            xpoint_displs,tetra_displs,xtetra_displs,
                            met_displs,intval_displs,ext_comm_displs,
                            int_comm_index_displs,rcv_np,rcv_xp,rcv_ne,
                            rcv_xt,rcv_nmet,rcv_int_comm_index,
                            rcv_next_node_comm,rcv_nitem_int_node_comm) ) {
    ier = 0;
  }

fail:
  /* Free memory */
  /* 1: Mesh data */
  _MMG5_SAFE_FREE(rcv_np);
  _MMG5_SAFE_FREE(rcv_point);
  _MMG5_SAFE_FREE(point_displs);

  _MMG5_SAFE_FREE(rcv_xp);
  _MMG5_SAFE_FREE(rcv_xpoint);
  _MMG5_SAFE_FREE(xpoint_displs);

  _MMG5_SAFE_FREE(rcv_ne);
  _MMG5_SAFE_FREE(rcv_tetra);
  _MMG5_SAFE_FREE(tetra_displs);

  _MMG5_SAFE_FREE(rcv_xt);
  _MMG5_SAFE_FREE(rcv_xtetra);
  _MMG5_SAFE_FREE(xtetra_displs);

  _MMG5_SAFE_FREE(rcv_nmet);
  _MMG5_SAFE_FREE(rcv_met);
  _MMG5_SAFE_FREE(met_displs);

  /* 2: communicators data */
  _MMG5_SAFE_FREE(rcv_int_comm_index);
  _MMG5_SAFE_FREE(rcv_intvalues);
  _MMG5_SAFE_FREE(intval_displs);

  _MMG5_SAFE_FREE(rcv_next_node_comm);
  _MMG5_SAFE_FREE(rcv_nitem_ext_tab);
  _MMG5_SAFE_FREE(rcv_color_in_tab);
  _MMG5_SAFE_FREE(rcv_color_out_tab);
  _MMG5_SAFE_FREE(ext_comm_displs);

  _MMG5_SAFE_FREE(rcv_nitem_int_node_comm);
  _MMG5_SAFE_FREE(rcv_node2int_node_comm_index1);
  _MMG5_SAFE_FREE(rcv_node2int_node_comm_index2);
  _MMG5_SAFE_FREE(int_comm_index_displs);

  _MMG5_SAFE_FREE(int_node_comm->intvalues);

  return(ier);
}

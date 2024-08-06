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
 * \file mergemesh.c
 * \brief Merge the mesh.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 */
#include "parmmg.h"
#include "mpipack_pmmg.h"
#include "mpiunpack_pmmg.h"
#include "moveinterfaces_pmmg.h"

/**
 * \param solI pointer toward the destination solution structure
 * \param solJ pointer toward the source solution structure
 * \param ip index of destination point
 * \param k pos in the group communicator (if provided) or idx (if not) of the src point
 * \param ni_n_c_idx1 group communicator (facultative)
 * \param warn 0 if no errors for now
 * \param whoisssol name of the solution structure on which qe are working (met, ls..)
 *
 * Copy the value of the solution J at point k or n2i_n_c_idx1[k] into the point \a
 * ip of the solution I.
 *
 */
static inline
void PMMG_grpJinI_copySol(MMG5_pSol solI,MMG5_pSol solJ,int ip,int k,
                          int* n2i_n_c_idx1,int8_t *warn,char *whoissol) {
  int size,ip2;

  if ( solI && solI->m ) {
    size = solJ->size;
    if ( (!solJ) || (!solJ->m) ) {
      if ( !(*warn) ) {
        (*warn) = 1;
        printf("  ## Error: unable to merge %s:"
               " group I has a metric while group J don't.\n",whoissol);
      }
    }
    else {
      assert( (size==solJ->size) && "size issues" );
      if ( n2i_n_c_idx1 ) {
        ip2 = n2i_n_c_idx1[k];
      }
      else {
        ip2 = k;
      }
      memcpy(&solI->m[size*ip],&solJ->m[size*ip2],size*sizeof(double));
    }
  }
}

/**
 * \param mesh pointer toward a mesh structure
 * \param met pointer toward a solution structure storing a metric
 * \param ls pointer toward a solution structure storing a level-set
 * \param disp pointer toward a solution structure storing a displacement
 * \param field pointer toward an array of solutions fields
 *
 * \return 1 if the np and npmax field of all the structures are matching, 0 otherwise
 *
 * Check the consistency between the np and npmax of the different structures.
 *
 */
static inline
int PMMG_check_solsNp(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol ls,
                      MMG5_pSol disp,MMG5_pSol field) {
  int i;

  if ( met && met->m ) {
    if ( met->np    != mesh->np    ) {
      return 0;
    }
    if ( met->npmax != mesh->npmax ) {
      return 0;
    }
  }

  if ( ls && ls->m ) {
    if ( ls->np    != mesh->np    ) {
      return 0;
    }
    if ( ls->npmax != mesh->npmax ) {
      return 0;
    }
  }

  if ( disp && disp->m ) {
    if ( disp->np    != mesh->np    ) {
      return 0;
    }
    if ( disp->npmax != mesh->npmax ) {
      return 0;
    }
  }

  if ( field ) {
    for ( i=0; i < mesh->nsols; ++i ) {
      assert ( field[i].np    == mesh->np    );
      assert ( field[i].npmax == mesh->npmax );
    }
  }

  return 1;
}

/**
 * \param mesh pointer toward a mesh structure
 * \param met pointer toward a solution structure storing a metric
 * \param ls pointer toward a solution structure storing a level-set
 * \param disp pointer toward a solution structure storing a displacement
 * \param field pointer toward an array of solutions fields
 *
 * \return 1 if the npmax field of all the structures are matching, 0 otherwise
 *
 * Check the consistency between the npmax values of the different structures.
 *
 */
static inline
int PMMG_check_solsNpmax(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol ls,
                      MMG5_pSol disp,MMG5_pSol field) {
  int i;

  if ( met && met->m ) {
    if ( met->npmax != mesh->npmax ) {
      return 0;
    }
  }

  if ( ls && ls->m ) {
    if ( ls->npmax != mesh->npmax ) {
      return 0;
    }
  }

  if ( disp && disp->m ) {
    if ( disp->npmax != mesh->npmax ) {
      return 0;
    }
  }

  if ( field ) {
    for ( i=0; i < mesh->nsols; ++i ) {
      assert ( field[i].npmax == mesh->npmax );
    }
  }

  return 1;
}

/**
 * \param mesh pointer toward a mesh structure
 * \param met pointer toward a solution structure storing a metric
 * \param ls pointer toward a solution structure storing a level-set
 * \param disp pointer toward a solution structure storing a displacement
 * \param field pointer toward an array of solutions fields
 * \param c coordinates of the new point to create (\a ip)
 * \param tag tag of the new point to create (\a ip)
 *
 * \return 0 if fail, \a ip, the index of the new point if success.
 *
 * Realloc all the array linked to the point array (included) in order to be
 * able to insert a new point. Perform the point creation if success.
 *
 */
static inline
int PMMG_realloc_pointAndSols(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol ls,
                              MMG5_pSol disp,MMG5_pSol field,double *c,int16_t tag,int src) {
  MMG5_pSol psl;
  int ip       = 0;
  int oldnpmax = mesh->npmax;
  int oldnp    = mesh->np;
  int fail     = 0;
  int i;

  assert ( PMMG_check_solsNpmax(mesh,met,ls,disp,field) );

  MMG3D_POINT_REALLOC(mesh,met,ip,mesh->gap,
                     printf("  ## Error: unable to merge group points\n");
                     MMG5_INCREASE_MEM_MESSAGE();
                     met->np = mesh->np; met->npmax = mesh->npmax;
                     ip = 0,
                     c,tag,src);

  if ( !ip ) {
    assert ( PMMG_check_solsNpmax(mesh,met,ls,disp,field) );
    return 0;
  }

  /* Reallocation of solution structures */
  /* ls field */
  if ( ls && ls->m ) {
    PMMG_REALLOC(mesh,ls->m,ls->size*(mesh->npmax+1),ls->size*(ls->npmax+1),
                 double,"larger level-set",fail = 1);
    ls->npmax = mesh->npmax;
    ls->np    = mesh->np;
  }

  /* displacement field */
  if ( disp && disp->m && (!fail) ) {
    PMMG_REALLOC(mesh,disp->m,disp->size*(mesh->npmax+1),disp->size*(disp->npmax+1),
                 double,"larger displacement",fail = 2);
    disp->npmax = mesh->npmax;
    disp->np    = mesh->np;
  }

  if ( field && (!fail) ) {
    for ( i=0; i < mesh->nsols; ++i ) {
      psl = &field[i];
      if ( psl->m ) {
        PMMG_REALLOC(mesh,psl->m,psl->size*(mesh->npmax+1),psl->size*(psl->npmax+1),
                     double,"larger field",fail = 3+i);
        psl->npmax = mesh->npmax;
        psl->np    = mesh->np;
      }
    }
  }

  if ( !fail )  return ip;

  for ( i=0; i < fail-3; ++i ) {
    psl = &field[i];
    PMMG_REALLOC(mesh,psl->m,psl->size*(oldnpmax+1),psl->size*(psl->npmax+1),
                 double,"smaller field",);
    psl->npmax = oldnpmax;
    psl->np    = oldnp;
  }

  if ( disp && disp->m ) {
    PMMG_REALLOC(mesh,disp->m,disp->size*(oldnpmax+1),disp->size*(disp->npmax+1),
                 double,"smaller displacement",);
    disp->npmax = oldnpmax;
    disp->np    = oldnp;
  }

  if ( ls && ls->m ) {
    PMMG_REALLOC(mesh,ls->m,ls->size*(oldnpmax+1),ls->size*(ls->npmax+1),
                 double,"smaller level-set",);
    ls->npmax = oldnpmax;
    ls->np    = oldnp;
  }

  if ( met && met->m ) {
    PMMG_REALLOC(mesh,met->m,met->size*(oldnpmax+1),met->size*(met->npmax+1),
                 double,"smaller metric",);
    met->npmax = oldnpmax;
    met->np    = oldnp;
  }

  MMG3D_delPt(mesh,ip);

  PMMG_REALLOC(mesh,mesh->point,oldnpmax+1,mesh->npmax+1,
               MMG5_Point,"smaller point array",);
  mesh->npmax = oldnpmax;
  mesh->np    = oldnp;
  mesh->npnil = 0;

  assert ( PMMG_check_solsNp(mesh,met,ls,disp,field) );

  printf("  ## Error: unable to merge group points\n");
  MMG5_INCREASE_MEM_MESSAGE();

  return 0;
}


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
  MMG5_pSol      metI,metJ,lsI,lsJ,dispI,dispJ;
  MMG5_pSol      fieldI,fieldJ,pslI,pslJ;
  MMG5_pPoint    pptI,pptJ;
  MMG5_pxPoint   pxpI,pxpJ;
  int            *intvalues,nitem_int_node_comm;
  int            *node2int_node_comm_index1,*node2int_node_comm_index2;
  int            k,poi_id_glo,ip,is,src;
  static int8_t  warnMet=0,warnLs=0,warnDisp=0,warnField=0;

  intvalues = parmesh->int_node_comm->intvalues;

  meshI                     = grpI->mesh;
  metI                      = grpI->met;
  lsI                       = grpI->ls;
  dispI                     = grpI->disp;
  fieldI                    = grpI->field;
  meshJ                     = grpJ->mesh;
  metJ                      = grpJ->met;
  lsJ                       = grpJ->ls;
  dispJ                     = grpJ->disp;
  fieldJ                    = grpJ->field;
  nitem_int_node_comm       = grpJ->nitem_int_node_comm;
  node2int_node_comm_index1 = grpJ->node2int_node_comm_index1;
  node2int_node_comm_index2 = grpJ->node2int_node_comm_index2;


  assert ( PMMG_check_solsNpmax(meshI,metI,lsI,dispI,fieldI) );
  assert ( PMMG_check_solsNpmax(meshJ,metJ,lsJ,dispJ,fieldJ) );

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
#ifdef USE_POINTMAP
    src = pptJ->src;
#else
    src = 1;
#endif
    /* Point pptJ is not found in the merged mesh. Add it */
    if ( !intvalues[ poi_id_glo ] ) {
      ip = MMG3D_newPt(meshI,pptJ->c,pptJ->tag,src);
      if ( !ip ) {
        /* reallocation of point table and associated solutions structures*/
        ip = PMMG_realloc_pointAndSols(meshI,metI,lsI,dispI,fieldI,pptJ->c,pptJ->tag,src);
        if ( !ip ) {
          return 0;
        }
        metI   = grpI->met;
        lsI    = grpI->ls;
        dispI  = grpI->disp;
        fieldI = grpI->field;
      }
      assert( (ip <= meshI->npmax) && "run out of points" );

      pptJ->tmp = ip;
      intvalues[ poi_id_glo ] = ip;
      pptI = &meshI->point[ip];
      pptI->n[0] = pptJ->n[0];
      pptI->n[1] = pptJ->n[1];
      pptI->n[2] = pptJ->n[2];
      pptI->ref = pptJ->ref;

      /* Add xpoint if needed */
      if ( pptJ->xp ) {
        pxpJ = &meshJ->xpoint[pptJ->xp];
        pptI = &meshI->point[ip];
        pxpI = &meshI->xpoint[pptI->xp];
        assert( (pptI->xp <= meshI->xpmax) && "increase xpoints" );
        assert( (pptI->xp > 0) && "negative xpoints" );
        memcpy(pxpI,pxpJ,sizeof(MMG5_xPoint));
      }

      PMMG_grpJinI_copySol(metI ,metJ ,ip,k,node2int_node_comm_index1,&warnMet ,"metric");
      PMMG_grpJinI_copySol(lsI  ,lsJ  ,ip,k,node2int_node_comm_index1,&warnLs  ,"ls");
      PMMG_grpJinI_copySol(dispI,dispJ,ip,k,node2int_node_comm_index1,&warnDisp,"displacement");

      if ( fieldI ) {
        for ( is=0; is<meshI->nsols; ++is ) {
          pslI = &fieldI[is];
          pslJ = &fieldJ[is];
          PMMG_grpJinI_copySol(pslI,pslJ,ip,k,node2int_node_comm_index1,&warnField,"displacement");
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
  if ( lsI && lsI->m) {
    lsI->np = lsI->np;
  }

  if ( dispI && dispI->m) {
    dispI->np = meshI->np;
  }

  if ( fieldI ) {
    for ( is=0; is<meshI->nsols; ++is ) {
      pslI = &fieldI[is];
      pslI->np = meshI->np;
    }
  }

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
  MMG5_pSol      lsI,lsJ;
  MMG5_pSol      dispI,dispJ;
  MMG5_pSol      fieldI,fieldJ,pslI,pslJ;
  MMG5_pPoint    pptI,pptJ;
  MMG5_pxPoint   pxpI,pxpJ;
  int            ip,ier,k,is,src;
  static int8_t  warnMet=0,warnLs=0,warnDisp=0,warnField=0;

  meshI  = grpI->mesh;
  metI   = grpI->met;
  lsI    = grpI->ls;
  dispI  = grpI->disp;
  fieldI = grpI->field;

  /** Loop over points and add the ones that are not already in the merged
   * mesh (meshI) */
  meshJ  = grpJ->mesh;
  metJ   = grpJ->met;
  lsJ    = grpJ->ls;
  dispJ  = grpJ->disp;
  fieldJ = grpJ->field;

  for ( k=1; k<=meshJ->np; k++ ) {
    pptJ = &meshJ->point[k];
    if ( !MG_VOK(pptJ) ) continue;
    if ( pptJ->tmp )     continue;
#ifdef USE_POINTMAP
    src = pptJ->src;
#else
    src = 1;
#endif
    ip = MMG3D_newPt(meshI,pptJ->c,pptJ->tag,src);
    if ( !ip ) {
      /* reallocation of point table and associated solutions structures*/
      ip = PMMG_realloc_pointAndSols(meshI,metI,lsI,dispI,fieldI,pptJ->c,pptJ->tag,src);
      if ( !ip ) {
        return 0;
      }
      metI   = grpI->met;
      lsI    = grpI->ls;
      dispI  = grpI->disp;
      fieldI = grpI->field;
    }
    pptJ->tmp = ip;
    pptI = &meshI->point[ip];
    pptI->n[0] = pptJ->n[0];
    pptI->n[1] = pptJ->n[1];
    pptI->n[2] = pptJ->n[2];
    pptI->ref = pptJ->ref;

    /* Add xpoint if needed */
    ier = 1;
    if ( pptJ->xp ) {
      pxpJ = &meshJ->xpoint[pptJ->xp];
      pptI = &meshI->point[ip];
      pxpI = &meshI->xpoint[pptI->xp];
      memcpy(pxpI,pxpJ,sizeof(MMG5_xPoint));
    }
    PMMG_grpJinI_copySol(metI ,metJ ,ip,k,NULL,&warnMet ,"metric");
    PMMG_grpJinI_copySol(lsI  ,lsJ  ,ip,k,NULL,&warnLs  ,"ls");
    PMMG_grpJinI_copySol(dispI,dispJ,ip,k,NULL,&warnDisp,"displacement");

    if ( fieldI ) {
      for ( is=0; is<meshI->nsols; ++is ) {
        pslI = &fieldI[is];
        pslJ = &fieldJ[is];
        PMMG_grpJinI_copySol(pslI,pslJ,ip,k,NULL,&warnField,"displacement");
      }
    }
  }

  metI->np = meshI->np;
  if ( lsI && lsI->m) {
    lsI->np = lsI->np;
  }

  if ( dispI && dispI->m) {
    dispI->np = meshI->np;
  }

  if ( fieldI ) {
    for ( is=0; is<meshI->nsols; ++is ) {
      pslI = &fieldI[is];
      pslI->np = meshI->np;
    }
  }

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
 * \return 0 if fail, 1 if success
 *
 * Merge all meshes (mesh elements + internal communicator) of a group into the
 * first mesh of the group. This function free the adjacency array.
 *
 * \remark the tetra must be packed.
 *
 */
int PMMG_merge_grps( PMMG_pParMesh parmesh,int target )
{
  PMMG_pGrp      listgrp,grp;
  MMG5_pMesh     mesh0;
  PMMG_pInt_comm int_node_comm,int_face_comm;
  int            *face2int_face_comm_index1,*face2int_face_comm_index2;
  int            imsh,k,iel;

  if ( !parmesh->ngrp ) return 1;

  listgrp  = parmesh->listgrp;

  /** Free the adjacency array: a possible improvement is to update it */
  mesh0 = listgrp[0].mesh;

  if ( !mesh0 ) return 1;

  /* Use mark field to store previous grp index */
  if( target == PMMG_GRPSPL_DISTR_TARGET ) PMMG_set_color_tetra( parmesh,0 );

  if ( mesh0->adja )
    PMMG_DEL_MEM(mesh0, mesh0->adja,int, "adjacency table" );

  if ( parmesh->ngrp == 1 ) return 1;

  /** Use the internal communicators to store the interface entities indices */
  int_node_comm = parmesh->int_node_comm;
  PMMG_CALLOC(parmesh,int_node_comm->intvalues,int_node_comm->nitem,int,
              "node communicator",return 0);

  int_face_comm = parmesh->int_face_comm;
  PMMG_CALLOC(parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,
              "face communicator",goto fail_ncomm);

  //DEBUGGING:
  //saveGrpsToMeshes(listgrp,parmesh->ngrp,parmesh->myrank,"BeforeMergeGrp");

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

    /* Use mark field to store previous grp index */
    if( target == PMMG_GRPSPL_DISTR_TARGET )
      PMMG_set_color_tetra( parmesh,imsh );

    /** Step 2: Merge internal points of the mesh mesh into the mesh0 mesh */
    if ( !PMMG_mergeGrpJinI_internalPoints(&listgrp[0],grp) )
      goto fail_comms;

    /* Step 3: Add the interfaces tetra of the imsh mesh to the mesh0 mesh */
    if ( !PMMG_mergeGrpJinI_interfaceTetra(parmesh,&listgrp[0],grp) )
      goto fail_comms;

    /** Step 4: Merge internal tetras of the imsh mesh into the mesh0 mesh */
    if ( !PMMG_mergeGrpJinI_internalTetra(&listgrp[0],grp) )
      goto fail_comms;

    mesh0->npi = mesh0->np;
    mesh0->nei = mesh0->ne;

    /* Free merged mesh */
    PMMG_grp_free(parmesh,grp);
  }


  /** Step 5: Update the communicators */

  if ( !PMMG_mergeGrps_communicators(parmesh) ) goto fail_comms;

  PMMG_REALLOC(parmesh,parmesh->listgrp,1,parmesh->ngrp,PMMG_Grp,"listgrp",return 0;);
  parmesh->ngrp = 1;

  /** Step 6: Update tag on points, tetra */
  if ( !PMMG_updateTag(parmesh) ) goto fail_comms;

  assert ( PMMG_check_solsNp ( parmesh->listgrp[0].mesh,
                               parmesh->listgrp[0].met,
                               parmesh->listgrp[0].ls,
                               parmesh->listgrp[0].disp,
                               parmesh->listgrp[0].field) );

  return 1;

fail_comms:
  PMMG_DEL_MEM(parmesh,int_face_comm->intvalues,int,"face communicator");

fail_ncomm:
  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,int,"node communicator");
  parmesh->listgrp[0].mesh->npi = parmesh->listgrp[0].mesh->np;
  parmesh->listgrp[0].mesh->nei = parmesh->listgrp[0].mesh->ne;

  return 0;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param rcv_grps array of groups to allocate and fill.
 * \param rcv_int_node_comm array of internal communicators to allocate and fill.
 * \param rcv_next_node_comm array of number external comm to alloc and fill.
 * \param rcv_ext_node_comm array of external comms to allocate and fill.
 *
 * \return 0 if fail, 1 otherwise (no allreduce over procs)
 *
 * Gather the parmeshes on the proc 0.
 *
 * The data received from processor \a k
 * are strored in the \a k index of the arrays \a rcv_grps, \a
 * rcv_int_node_comm, \a rcv_next_node_comm and \a rcv_ext_node_comm.
 *
 *
 * \warning We must have at most 1 group per parmesh
 *
 */
static inline
int PMMG_gather_parmesh( PMMG_pParMesh parmesh,
                         PMMG_pGrp *rcv_grps,
                         PMMG_pInt_comm *rcv_int_node_comm,
                         int **rcv_next_node_comm,
                         PMMG_pExt_comm **rcv_ext_node_comm ) {

  size_t     pack_size_tot,next_disp,*displs,buf_idx;
  int        *rcv_pack_size,ier,ier_glob,k,ier_pack;
  int        nprocs,root,pack_size;
  char       *rcv_buffer,*ptr_to_free,*buffer;

  nprocs        = parmesh->nprocs;
  root          = parmesh->info.root;

  ier           = ier_glob = 1;

  (*rcv_grps)           = NULL;
  (*rcv_int_node_comm)  = NULL;
  (*rcv_next_node_comm) = NULL;
  (*rcv_ext_node_comm)  = NULL;
  rcv_pack_size         = NULL;
  displs                = NULL;
  rcv_buffer            = NULL;
  buffer                = NULL;

  /** 1: Memory alloc */
  if ( parmesh->myrank == root ) {
    PMMG_MALLOC( parmesh, rcv_pack_size        ,nprocs,int,"rcv_pack_size",ier=0);
    PMMG_MALLOC( parmesh, displs               ,nprocs,size_t,"displs for gatherv",ier=0);
    PMMG_CALLOC( parmesh, (*rcv_grps)          ,nprocs,PMMG_Grp,"rcv_grps",ier=0);
    PMMG_MALLOC( parmesh, (*rcv_int_node_comm) ,nprocs,PMMG_Int_comm,"rcv_int_comm" ,ier=0);
    PMMG_MALLOC( parmesh, (*rcv_next_node_comm),nprocs,int,"rcv_next_comm" ,ier=0);
    PMMG_MALLOC( parmesh, (*rcv_ext_node_comm) ,nprocs,PMMG_pExt_comm,"rcv_ext_comm" ,ier=0);
  }

  /** 2: Gather pack size of parmeshes on proc 0 */
  pack_size = PMMG_mpisizeof_parmesh ( parmesh );

#ifndef NDEBUG
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);
  if ( !ier_glob ) {
    return ier_glob;
  }
#endif

  MPI_CHECK(MPI_Gather(&pack_size,1,MPI_INT,rcv_pack_size,1,MPI_INT,root,parmesh->comm),ier = 0);

  /** 3: Gather compressed parmeshes */
  /* Compute data for gatherv: displacement array and receiver buffer size */
  if ( parmesh->myrank == root ) {
    displs[0] = 0;
    for ( k=1; k<nprocs; ++k ) {
      next_disp = displs[k-1] + rcv_pack_size[k-1];
      displs[k] = next_disp;
    }

    /* On root, we will gather all the meshes in rcv_buffer so we have to
     * compute the total pack size */
    pack_size_tot        = (size_t)(displs[nprocs-1])+(size_t)(rcv_pack_size[nprocs-1]);
    assert ( pack_size_tot < SIZE_MAX && "SIZE_MAX overflow" );

    /* root will write directly in the suitable position of rcv_buffer */
    buf_idx = displs[root];
  }
  else {
    /* on ranks other than root we just need to store the local mesh so buffer
     * will be of size pack_size */
    pack_size_tot = pack_size;
    /* we will write the mesh at the starting position */
    buf_idx = 0;
  }

  PMMG_MALLOC( parmesh,rcv_buffer,pack_size_tot,char,"rcv_buffer",ier=0);

  /* Parmesh compression */
  buffer = &rcv_buffer[buf_idx];

  /* Save input allocated address to avoid arrors at unalloc */
  ptr_to_free = rcv_buffer;

#ifndef NDEBUG
  /* Remark: in release mode, a non allocated buffer used in gatherv creates a
     segfault over all procs */
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);
  if ( !ier_glob ) {
    return ier_glob;
  }
#endif

  /* /!\ mpipack_parmesh and mpiunpack_parmesh are modifying the buffer pointer
   * making it not valid for realloc / unalloc */

  /* Save adress of buffer because it will be set the the end of the char array
  by the \a PMMG_mpipack_parmesh function */
  char *buffer_to_send = buffer;
  ier_pack = PMMG_mpipack_parmesh ( parmesh ,&buffer );

  /* Do not use \a buffer pointer after this call: it points toward the end of
   * the packed array which is useless */
  buffer = NULL;

  assert ( ier_pack );

  /* Gather the packed parmeshes */
  ier = MG_MIN ( ier, ier_pack );

  /* Here the gatherv call has been replaced by a send/recv to avoid errors when
   * displacements overflow the INT_MAX value */
  if (parmesh->myrank == root) {
    int i;
    for ( i = 0; i < nprocs; ++i ) {
      if ( i != root ) {
        MPI_CHECK(
          MPI_Recv(rcv_buffer + displs[i], rcv_pack_size[i], MPI_CHAR, i,
                   MPI_MERGEMESH_TAG, parmesh->comm, MPI_STATUS_IGNORE),
          ier = 0);
      }
    }
  } else {
    MPI_CHECK(
      MPI_Send(buffer_to_send, pack_size, MPI_CHAR, root, MPI_MERGEMESH_TAG,parmesh->comm),
      ier = 0);
  }

  /** 4: Unpack parmeshes */
#ifndef NDEBUG
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);
  if ( !ier_glob ) {
    return ier_glob;
  }
#endif

  if ( parmesh->myrank == root ) {
    for ( k=0; k<nprocs; ++k ) {
      ier_pack = PMMG_mpiunpack_parmesh ( parmesh,(*rcv_grps),k,(*rcv_int_node_comm)+k,
                                          (*rcv_next_node_comm)+k,(*rcv_ext_node_comm)+k,
                                          &rcv_buffer );
      ier = MG_MIN(ier_pack,ier);
    }
  }

  /** Free the memory */
  /* Free temporary arrays */
  PMMG_DEL_MEM(parmesh,rcv_pack_size,int,"rcv_pack_size");
  PMMG_DEL_MEM(parmesh,displs,int,"displs");
  /* the address of rcv_buffer is modified by packing/unpacking so it is needed
   * to send the initially allocated address stored in to the unalloc macro */
  PMMG_DEL_MEM(parmesh,ptr_to_free,char,"rcv_buffer");

  return ier;
}


/**
 * \param parmesh pointer toward the parmesh structure.
 * \param rcv_grps array containing groups data.
 * \param rcv_int_node_comm array containing internal communicators.
 * \param rcv_next_node_comm array of number external comm.
 * \param rcv_ext_node_comm array containing external comms.
 *
 * \return 0 if fail, 1 otherwise (no allreduce over procs)
 *
 * Merge the parmeshes data contained in the rcv_* arrays into 1 parmesh with 1
 * group.
 *
 */
static inline
int PMMG_mergeParmesh_rcvParMeshes ( PMMG_pParMesh parmesh,PMMG_pGrp rcv_grps,
                                     PMMG_pInt_comm rcv_int_node_comm,
                                     int *rcv_next_node_comm,
                                     PMMG_pExt_comm *rcv_ext_node_comm ) {
  PMMG_pGrp      grp,grp_1,grp_2;
  PMMG_pInt_comm int_node_comm_1,int_node_comm_2;
  PMMG_pExt_comm ext_node_comm_1,ext_node_comm_2;
  MMG5_pMesh     mesh;
  MMG5_pPoint    point_1,point_2,ppt;
  MMG5_pxPoint   xpoint_1,pxp;
  MMG5_pTetra    tetra_1,pt;
  MMG5_pxTetra   xtetra_1,pxt;
  MMG5_pSol      met,ls,disp,psl,met_1,ls_1,disp_1,psl_1;
  size_t         memAv;
  int            is,*int_comm_index_1,*int_comm_index_2;
  int            next_node_comm_1,next_node_comm_2;
  int            *intvalues_1,*intvalues_2,nitem_1,nitem_2;
  int            nprocs,k,i,j,idx,color_in,color_out;
  int            np,ne,xt,ne_tot,xt_tot,ismet,isls,isdisp;
  int            type[MMG5_NSOLS_MAX];

  nprocs = parmesh->nprocs;

  if ( parmesh->myrank != parmesh->info.root ) return 1;

  assert ( !parmesh->listgrp );

  PMMG_CALLOC(parmesh,parmesh->listgrp,1,PMMG_Grp,"listgrp", return 0);

  MMG3D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &parmesh->listgrp[0].mesh,
                   MMG5_ARG_ppMet, &parmesh->listgrp[0].met, MMG5_ARG_end );


  grp    = &parmesh->listgrp[0];
  mesh   = grp->mesh;
  met    = grp->met ;

  /* Set maximum memory */
  mesh->memMax = parmesh->memGloMax;

  np = 0;

  /** Mesh renumbering to have the same indices at points shared by multiple
   * processors; The new indices are stored in the tmp field of the MMG5_Point
   * structure. */
  for ( k=0; k<nprocs; ++k ) {
    grp_1            = &rcv_grps[k];
    int_node_comm_1  = &rcv_int_node_comm[k];
    ext_node_comm_1  = rcv_ext_node_comm[k];
    next_node_comm_1 = rcv_next_node_comm[k];

    /* Points and internal communicators for proc \a k */
    point_1     = grp_1->mesh->point;
    intvalues_1 = int_node_comm_1->intvalues;

    /* Travel through the external communicators that lists the points at the
     * interface of the procs color_in and color_out: if color_in<color_out,
     * create new indices for the points over the color_in proc. If
     * color_out<color_in, get the indices from the matching points over the
     * proc color_out. */
    for ( i=0; i<next_node_comm_1; ++i ) {
      nitem_1 = ext_node_comm_1[i].nitem;

      int_comm_index_1 = ext_node_comm_1[i].int_comm_index;

      /* External communicator k->color_out */
      color_in    = ext_node_comm_1[i].color_in;
      color_out   = ext_node_comm_1[i].color_out;

      assert( color_in==k );
      assert( color_in!=color_out );

      if ( color_in < color_out ) {
        /* New point */
        for ( j=0; j<nitem_1; ++j ) {
          if ( !point_1[intvalues_1[int_comm_index_1[j]]].tmp ) {
            point_1[intvalues_1[int_comm_index_1[j]]].tmp = ++np;
          }
        }
      }
      else {
        /* Get the point index from the neighbouring proc */

        /* Find the matching external communicator over proc color_out (->k) */
        grp_2            = &rcv_grps[color_out];
        int_node_comm_2  = &rcv_int_node_comm[color_out];
        ext_node_comm_2  = rcv_ext_node_comm[color_out];
        next_node_comm_2 = rcv_next_node_comm[color_out];

        for ( j=0; j<next_node_comm_2; ++j ) {
          nitem_2 = ext_node_comm_2[j].nitem;
          int_comm_index_2 = ext_node_comm_2[j].int_comm_index;

          if ( ext_node_comm_2[j].color_in  == color_out &&
               ext_node_comm_2[j].color_out == color_in  )  break;
        }
        assert ( j < next_node_comm_2 );
        assert ( nitem_1 == nitem_2 );


        /* Points and internal communicators for proc \a color_out */
        point_2     = grp_2->mesh->point;
        intvalues_2 = int_node_comm_2->intvalues;

        /* Update point indices (stored in the tmp field) */
        for ( j=0; j<nitem_1; ++j ) {
          point_1[intvalues_1[int_comm_index_1[j]]].tmp =
            point_2[intvalues_2[int_comm_index_2[j]]].tmp;
        }
      }
    }

    /* Create new indices for the points that haven't been seen. */
    for ( i=1; i <= grp_1->mesh->np; ++i ) {
      if ( !point_1[i].tmp ) {
        point_1[i].tmp = ++np;
      }
    }
  }

  /** Tetra + xTetra */
  ne_tot = xt_tot = 0;
  for ( k=0; k<nprocs; ++k ) {
    ne_tot += rcv_grps[k].mesh->ne;
    xt_tot += rcv_grps[k].mesh->xt;
  }

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

  ne = xt = idx = 0;
  for ( k=0; k<nprocs; ++k ) {
    grp_1       = &rcv_grps[k];
    xtetra_1    = grp_1->mesh->xtetra;
    tetra_1     = grp_1->mesh->tetra;
    point_1     = grp_1->mesh->point;

    for ( i=1; i<=grp_1->mesh->ne; ++i ) {
      pt = &mesh->tetra[++ne];

      for ( j=0; j<4; ++j ) {
        tetra_1[i].v[j] = point_1[tetra_1[i].v[j]].tmp;
      }

      memcpy(pt,&tetra_1[i],sizeof(MMG5_Tetra));

      if ( tetra_1[i].xt ) {
        pxt = &xtetra_1[tetra_1[i].xt];
        memcpy(&mesh->xtetra[++xt],pxt,sizeof(MMG5_xTetra));
        pt->xt = xt;
      }
    }
  }
  mesh->xt=xt;

  /** Points and solutions */
  mesh->np = np;
  mesh->npmax = mesh->np;
  mesh->npnil = 0;
  MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(MMG5_Point),"merge point",
                fprintf(stderr,"  Exit program.\n");
                return 0);
  MMG5_SAFE_CALLOC(mesh->point,mesh->npmax+1,MMG5_Point,return 0);

  /* Check the presence of solution structures */
  ismet  = 0;
  isls   = 0;
  isdisp = 0;
  for ( k=0; k<nprocs; ++k ) {
    if ( rcv_grps[k].mesh->np ) {
      ismet  = (rcv_grps[k].met  && rcv_grps[k].met->m  )? rcv_grps[k].met->size:0;
      isls   = (rcv_grps[k].ls   && rcv_grps[k].ls->m   )? rcv_grps[k].ls->size:0;
      isdisp = (rcv_grps[k].disp && rcv_grps[k].disp->m )? rcv_grps[k].disp->size:0;
      mesh->nsols = rcv_grps[k].mesh->nsols;

      /* Allocation of solution structures */
      if ( ismet ) {
        MMG3D_Set_solSize(mesh,met,MMG5_Vertex,mesh->np,rcv_grps[k].met->type);
      }
      if ( isls ) {
        PMMG_CALLOC(mesh,grp->ls,1,MMG5_Sol,"ls",return 0);
        MMG3D_Set_solSize(mesh,grp->ls,MMG5_Vertex,mesh->np,rcv_grps[k].ls->type);
      }
      if ( isdisp ) {
        PMMG_CALLOC(mesh,grp->disp,1,MMG5_Sol,"disp",return 0);
        MMG3D_Set_solSize(mesh,grp->disp,MMG5_Vertex,mesh->np,rcv_grps[k].disp->type);
      }

      if ( mesh->nsols ) {
        assert ( !grp->field );
        for ( i=0; i<mesh->nsols; ++i ) {
          type[i] = rcv_grps[k].field[i].type;
        }
        MMG3D_Set_solsAtVerticesSize( mesh,&grp->field,mesh->nsols,mesh->np,type);
      }
      break;
    }
  }
  ls    = grp->ls;
  disp  = grp->disp;

  for ( i=1; i<=mesh->np; ++i ) mesh->point[i].tag = MG_NUL;

  np = 0;
  for ( k=0; k<nprocs; ++k ) {
    grp_1   = &rcv_grps[k];
    point_1 = grp_1->mesh->point;
    met_1   = grp_1->met;
    ls_1    = grp_1->ls;
    disp_1  = grp_1->disp;

    for ( i=1; i<=grp_1->mesh->np; ++i ) {
      idx = point_1[i].tmp;
      assert(idx);
      ppt = &mesh->point[idx];

      if ( MG_VOK(ppt) ) {
        point_1[i].tmp = 0;
        continue;
      }

      memcpy(ppt,&point_1[i],sizeof(MMG5_Point));
      ppt->tmp = 0;

      /* Copy solution structures */
      if ( ismet ) {
        assert ( met && met->m );
        assert ( met_1->size == met->size );
        memcpy ( &met->m[idx*met->size],&met_1->m[i*met_1->size],
                 met->size*sizeof(double) );
      }
      if ( isls ) {
        assert ( ls && ls->m );
        assert ( ls_1->size == ls->size );
        memcpy ( &ls->m[idx*ls->size],&ls_1->m[i*ls_1->size],
                 ls->size*sizeof(double) );
      }
      if ( isdisp ) {
        assert ( disp && disp->m );
        assert ( disp_1->size == disp->size );
        memcpy ( &disp->m[idx*disp->size],&disp_1->m[i*disp_1->size],
                 disp->size*sizeof(double) );
      }
      if ( mesh->nsols ) {
        for ( is=0; is<mesh->nsols; ++is ) {
          psl   = &grp->field[is];
          psl_1 = &grp_1->field[is];
          assert ( psl && psl->m );
          assert ( psl_1->size == psl->size );
          memcpy ( &psl->m[idx*psl->size],&psl_1->m[i*psl_1->size],
                   psl->size*sizeof(double) );
        }
      }

      /* Count xpoints */
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
    grp_1       = &rcv_grps[k];
    point_1     = grp_1->mesh->point;
    xpoint_1    = grp_1->mesh->xpoint;

    for ( i=1; i<=grp_1->mesh->np; ++i ) {
      idx = point_1[i].tmp;

      if ( !idx ) continue;
      ppt = &mesh->point[idx];

      if ( !point_1[i].xp ) continue;

      pxp = &mesh->xpoint[++np];
      memcpy(pxp, &xpoint_1[point_1[i].xp],sizeof(MMG5_xPoint));
      ppt->xp = np;
    }
  }

  /** Recover mesh infos */
  grp_1 = &rcv_grps[0];
  assert ( grp_1->mesh );
  if ( !PMMG_copy_mmgInfo ( &grp_1->mesh->info,&mesh->info ) ) return 0;

  /** Recover mesh name */

  if ( parmesh->meshin ) {
    MMG3D_Set_inputMeshName (mesh, parmesh->meshin);
  }
  if ( parmesh->meshout ) {
    MMG3D_Set_outputMeshName(mesh, parmesh->meshout);
  }

  if ( met ) {
    if ( parmesh->metin ) {
      MMG3D_Set_inputSolName (mesh,met, parmesh->metin);
    }
    if ( parmesh->metout ) {
      MMG3D_Set_outputSolName(mesh,met, parmesh->metout);
    }
  }
  if ( ls ) {
    if ( parmesh->lsin ) {
      MMG3D_Set_inputSolName (mesh,ls, parmesh->lsin);
    }
  }
  if ( disp ) {
    if ( parmesh->dispin ) {
      MMG3D_Set_inputSolName (mesh,disp, parmesh->dispin);
    }
  }
  if ( mesh->nsols ) {
    for ( is=0; is < mesh->nsols; ++is ) {
      psl = &parmesh->listgrp[0].field[is];
      if ( parmesh->fieldin ) {
        MMG3D_Set_inputSolName (mesh, psl,parmesh->fieldin);
      }
      if ( parmesh->fieldout ) {
        MMG3D_Set_outputSolName(mesh, psl,parmesh->fieldout);
      }
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
  PMMG_pGrp      grp,rcv_grps;
  PMMG_pInt_comm int_node_comm,rcv_int_node_comm;
  PMMG_pExt_comm *rcv_ext_node_comm,enc;
  MMG5_pMesh     mesh;
  MMG5_pPoint    ppt;
  int            *rcv_next_node_comm;
  int            k,idx,ier,ieresult;

  ier = 1;
  grp = &parmesh->listgrp[0];

  assert ( parmesh->ngrp <= 1 );

  //DEBUGGING:
  // saveGrpsToMeshes( PMMG_pGrp listgrp, 0, parmesh->myrank, "mesh" );

  /** Step 1: Allocate internal communicator buffer and fill it: the
   *  intvalues array contains the indices of the matching nodes on the proc. */

#warning MEMORY: small inconsistency

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
  ier = PMMG_gather_parmesh ( parmesh,&rcv_grps,&rcv_int_node_comm,
                              &rcv_next_node_comm,&rcv_ext_node_comm );
  MPI_CHECK( MPI_Allreduce(&ier,&ieresult,1,MPI_INT,MPI_MIN,parmesh->comm),ieresult=0);

  if ( ieresult ) {
    /* Free useless data of the parmesh (rcv_* arrays contains all the needed data) */
    /* 1: groups */
    PMMG_listgrp_free( parmesh, &parmesh->listgrp, parmesh->ngrp );

    /* 2: communicators */
    PMMG_parmesh_Free_Comm(parmesh);

    /** Step 3: Proc 0 merges the meshes: We travel through the external
     * communicators to recover the numbering of the points shared with a lower
     * proc. The other points are concatenated with the proc 0. */
    ieresult = PMMG_mergeParmesh_rcvParMeshes(parmesh,rcv_grps,rcv_int_node_comm,
                                              rcv_next_node_comm,rcv_ext_node_comm);
    /* Only root may return a non-zero value */
    MPI_CHECK( MPI_Bcast(&ieresult,1,MPI_INT,parmesh->info.root,parmesh->comm),
               ieresult=0 );
  }

  /* Free useless receivers */
  if ( rcv_grps ) {
    PMMG_listgrp_free( parmesh, &rcv_grps, parmesh->nprocs );

    for ( k=0; k<parmesh->nprocs; ++k ) {
      PMMG_DEL_MEM(parmesh,rcv_int_node_comm[k].intvalues,int,"intvalues");
      enc = rcv_ext_node_comm[k];
      for ( idx=0; idx<rcv_next_node_comm[k]; ++idx ) {
        PMMG_DEL_MEM(parmesh,enc[idx].int_comm_index,int,"int_comm_index");
      }
      PMMG_DEL_MEM(parmesh,enc,PMMG_Ext_comm,"ext_node_comm");
    }

    PMMG_DEL_MEM(parmesh,rcv_int_node_comm,PMMG_Int_comm,"rcv_int_node_comm");
    PMMG_DEL_MEM(parmesh,rcv_next_node_comm,int,"rcv_next_node_comm");
    PMMG_DEL_MEM(parmesh,rcv_ext_node_comm,PMMG_Ext_comm,"rcv_ext_node_comm");
  }

  if ( !ieresult ) {
    fprintf ( stderr, " ## Warning: unable to merge meshes on one proc.\n"
              "            Try to save parallel mesh.");
    return ieresult;
  }

  if ( parmesh->myrank != parmesh->info.root ) {
    /** Empty normally */
    parmesh->memGloMax = parmesh->memCur;
    parmesh->ngrp = 0;
  }
  else {
    /** One group only */
    parmesh->ngrp = 1;

    /** Step 5: Update tag on points, tetra */
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

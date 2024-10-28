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
 * \file ls_pmmg.c
 * \brief Create implicit surface in distribuited mesh.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (InriaSoft)
 * \author Laetitia Mottet (UBordeaux)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Functions to perform the level-set discretization in parallel
 *
 */

#include "parmmg.h"
#include "mmgexterns_private.h"

/**
 * \param parmesh pointer toward a parmesh structure
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set values.
 * \param met pointer toward a metric (non-mandatory).
 *
 * \return 1 if success, 0 otherwise.
 *
 * \todo Multimaterial isSplit() on // interface
 *
 * Proceed to discretization of the implicit function carried by sol into mesh,
 * once values of sol have been snapped/checked
 *
 */
int PMMG_cuttet_ls(PMMG_pParMesh parmesh){
  MMG5_pTetra  pt,pt0;
  MMG5_pxTetra pxt,pxt0;
  MMG5_pPoint  p0,p1;
  MMG5_Hash    hash;

  PMMG_pExt_comm ext_node_comm,ext_edge_comm,ext_face_comm;
  PMMG_pInt_comm int_node_comm, int_edge_comm, int_face_comm;
  PMMG_pGrp      grp;
  MMG5_pSol      field;
  MMG5_pSol      psl;
  MMG5_pMat      mat;
  MMG5_pMesh     mesh;
  MMG5_pSol      met,sol;

  MMG5_int ne_init,ne_tmp;
  MMG5_int i,j,k,k0;
  MMG5_int ie,idx_tmp;
  MMG5_int pos,pos_edge,pos_node,pos_face;
  MMG5_int ip0,ip1,np,nb,ns;
  MMG5_int src;
  MMG5_int refext,refint,ref;
  MMG5_int vGlobNum[4],vx[6];
  MMG5_int tetra_sorted[3], node_sorted[3];
  MMG5_int *ne_tmp_tab,*vGlobNum_tab;

  static int8_t  mmgWarn = 0;
  int8_t ia;
  int8_t npneg,nface_added;
  int8_t already_split;

  const uint8_t *taued=NULL;
  uint8_t        tau[4];
  uint8_t        imin0,imin2;

  double c[3],v0,v1,s;

  int i_commf;
  int ifac,iploc,val_face;
  int flag;
  int ier;

  int nitem_grp_node_firstalloc,nitem_grp_face_firstalloc;
  int nitem_ext_face_init;
  int next_node_comm,next_face_comm,next_edge_comm;

  int nitem_int_node, nitem_grp_node, nitem_ext_node;
  int nitem_int_edge, nitem_grp_edge, nitem_ext_edge;
  int nitem_int_face, nitem_grp_face, nitem_ext_face;
  int nitem_grp_face_tmp;

  int color_in_node,color_out_node;
  int color_in_edge,color_out_edge;

  if ( parmesh->myrank == parmesh->info.root )
    fprintf(stdout,"\n      ## PMMG_cuttet_ls: Multimaterial not fully supported yet.\n");

  /* Ensure only one group on each proc */
  assert ( (parmesh->ngrp == 1 || parmesh->ngrp == 0) &&
           "Implemented for 1 group per rank" );

  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;
  sol  = parmesh->listgrp[0].ls;

  /* For now, does not support `nosplit` option in multimat  */
  // To be removed when supported
  for (i=0; i<mesh->info.nmat; i++) {
    mat    = &mesh->info.mat[i];
    ref    = mat->ref;
    refint = mat->rin;
    refext = mat->rex;
    if ( (ref == refint) && (ref == refext) ) {
      if ( parmesh->myrank == parmesh->info.root )
        fprintf(stderr,"\n      -- ERROR: The option `nosplit` in multimat is not supported yet.\n");
      return 0;
    }
  }

  /* Initialization */
  grp = &parmesh->listgrp[0];
  field = grp->field;
  next_node_comm = parmesh->next_node_comm;  // Number of communicator for nodes
  next_edge_comm = parmesh->next_edge_comm;  // Number of communicator for edges
  next_face_comm = parmesh->next_face_comm;  // Number of communicator for faces
  nitem_grp_node = grp->nitem_int_node_comm; // Number of initial total nodes in internal node communicator
  nitem_grp_edge = grp->nitem_int_edge_comm;
  nitem_grp_face = grp->nitem_int_face_comm; // Number of initial total faces in internal node communicator
  int_node_comm = parmesh->int_node_comm; // Internal node communicator
  int_edge_comm = parmesh->int_edge_comm; // Internal edge communicator
  int_face_comm = parmesh->int_face_comm; // Internal face communicator

  nitem_int_node = int_node_comm->nitem;
  nitem_int_edge = int_edge_comm->nitem;
  nitem_int_face = int_face_comm->nitem;

  ne_init = mesh->ne; // Initial number of tetra - before ls - needed in step 6.3

  /** STEP 1 - Reset flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  for (k=1; k<=mesh->ne; k++)
    mesh->tetra[k].flag = 0;

  /** STEP 2 - Approximate the number nb of intersection points on edges */
  nb = 0;

  /* Loop over tetra */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];

    if ( !MG_EOK(pt) ) {
      continue;
    }

    /* Loop over edges */
    for (ia=0; ia<6; ia++) {

      /* Get the points defining the edges */
      ip0 = pt->v[MMG5_iare[ia][0]];
      ip1 = pt->v[MMG5_iare[ia][1]];
      p0  = &mesh->point[ip0];
      p1  = &mesh->point[ip1];

      /* If both points have flag, then pass as these points have been treated */
      if ( p0->flag && p1->flag )  continue;

      /* Otherwise take the values at these points */
      v0  = sol->m[ip0];
      v1  = sol->m[ip1];

      /* If the points are not already exactly on the level-set
         and does not have the same sign, then this edge needs to be split */
      if ( fabs(v0) > MMG5_EPSD2 && fabs(v1) > MMG5_EPSD2 && v0*v1 < 0.0 ) {
        /* If the points have not been treated yet, assign a new flag and increase nb value */
        if ( !p0->flag ) {
          p0->flag = ++nb;
        }
        if ( !p1->flag ) {
          p1->flag = ++nb;
        }
      }
    }
  }
  if ( ! nb )  return 1;

  /* TODO:: test if the number of point proc by proc is correct */
  // Cannot be done here as it is an approximation. Otherwise, need to robustify step 2 above.
#ifndef NDEBUG
  /* TODO */
#endif

  /** STEP 3 - Memory allocation */
  /* STEP 3.1 - Initialize hash table for edges */
  if ( !MMG5_hashNew(mesh,&hash,nb,7*nb) ) return 0;

  /* STEP 3.2 - Realloc internal node communicators */
  PMMG_REALLOC(parmesh, grp->node2int_node_comm_index1,
               nitem_grp_node+nitem_grp_edge,
               nitem_grp_node,
               int,"Allocation of node2int_node_comm_index1", return 0);
  PMMG_REALLOC(parmesh, grp->node2int_node_comm_index2,
               nitem_grp_node+nitem_grp_edge,
               nitem_grp_node,
               int,"Allocation of node2int_node_comm_index2", return 0);
  nitem_grp_node_firstalloc = nitem_grp_node+nitem_grp_edge;

  /* STEP 3.3 - Realloc internal face communicators */
  PMMG_REALLOC(parmesh, grp->face2int_face_comm_index1,
               3*nitem_grp_face,
               nitem_grp_face,
               int,"Allocation of face2int_face_comm_index1", return 0);
  PMMG_REALLOC(parmesh, grp->face2int_face_comm_index2,
               3*nitem_grp_face,
               nitem_grp_face,
               int,"Allocation of face2int_face_comm_index2", return 0);
  nitem_grp_face_firstalloc = 3*nitem_grp_face;

  /* STEP 3.4 - Realloc external node communicator */
  for ( k=0; k<parmesh->next_node_comm; k++ ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    nitem_ext_node = ext_node_comm->nitem;
    PMMG_REALLOC(parmesh,ext_node_comm->int_comm_index,
                 nitem_ext_node+2*nb,
                 nitem_ext_node,
                 int,"Allocation of external node communicator",return 0);
    ext_node_comm->nitem_to_share = nitem_ext_node+2*nb;
  }

  /* STEP 3.5 - Realloc external face communicator */
  for ( k=0; k<parmesh->next_face_comm; k++ ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    nitem_ext_face = ext_face_comm->nitem;
    PMMG_REALLOC(parmesh,ext_face_comm->int_comm_index,
                 nitem_ext_face+2*nitem_ext_face,
                 nitem_ext_face,
                 int,"Allocation of external face communicator",return 0);
    ext_face_comm->nitem_to_share = nitem_ext_face+2*nitem_ext_face;
  }

  /* STEP 3.6 - Allocate all the other variables needed to be allocated ! */
  PMMG_CALLOC( parmesh,vGlobNum_tab,4*(nitem_grp_face),MMG5_int,"vGlobNum_tab",return 0 );
  PMMG_CALLOC( parmesh,ne_tmp_tab,nitem_grp_face+1,MMG5_int,"ne_tmp_tab",return 0 );

  /** STEP 4 - Identify required edges. Put hash.item[key].k = -1 */
  /* Loop over tetra */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    /* Check whether the tetra with reference ref should be split */
    // For now, the option nosplit is not supported
    // So we assume all the elements should be split
    // if ( !MMG5_isSplit(mesh,pt->ref,&refint,&refext) ) continue;

    /** Step 4.1 - Identification of edges belonging to a required tet */
    /* If the tetra is required MG_REQ  */
    if ( pt->tag & MG_REQ ) {
      /* Loop over the edges */
      for (ia=0; ia<6; ia++) {
        ip0 = pt->v[MMG5_iare[ia][0]];
        ip1 = pt->v[MMG5_iare[ia][1]];
        np  = -1;
        /* Add an edge to the edge table with hash.item[key].k = -1 */
        if ( !MMG5_hashEdge(mesh,&hash,ip0,ip1,np) )  return -1;
      }
      continue;
    }

    /** Step 4.2 - Identification of edges belonging to a (par)boundary or being explicitely required */
    /* If the xtetra associated to this tetra exists */
    if ( !pt->xt ) continue;

    /* Point towards the xtetra corresponding to the tetra...  */
    pxt = &mesh->xtetra[pt->xt];
    /* ... then loop over the faces */
    for (ia=0; ia<4; ia++) {

      /* (a) If the face is not a boundary MG_BDY, then continue */
      if ( !(pxt->ftag[ia] & MG_BDY) ) continue;

      /* (a) otherwise loop over the edges */
      for (j=0; j<3; j++) {

        /* (b) If the edges is not required, then continue */
        if ( !(pxt->tag[ MMG5_iarf[ia][j] ] & MG_REQ) ) continue;

        /* (b) otherwise get the extremity of the edges ... */
        ip0 = pt->v[MMG5_idir[ia][MMG5_inxt2[j]]];
        ip1 = pt->v[MMG5_idir[ia][MMG5_iprv2[j]]];
        np  = -1;

        /* (c) ... and add an edge to the edge table with hash.item[key].k = -1 */
        if ( !MMG5_hashEdge(mesh,&hash,ip0,ip1,np) )  return -1;
      }
    }
  }

  /** STEP 5 - Create points at iso-value. Fill or update edge hash table */
  /** STEP 5.1 - Create new points located on parallel interfaces */
  /* Internal edge communicator - intvalues stores:
    - point position as in node2int_node_comm_index2 if the edge is split
    - otherwise, -1 */
  PMMG_CALLOC(parmesh,int_edge_comm->intvalues,nitem_int_edge,int,"int_edge_comm intvalues",return 0);

  /* Loop on the internal edge communicator */
  for (i=0; i < nitem_grp_edge; i++) {

    ie  = grp->edge2int_edge_comm_index1[i]; // id of edge
    pos = grp->edge2int_edge_comm_index2[i]; // position in int_edge_comm->intvalues
    int_edge_comm->intvalues[pos] = -1;      // initialization at -1

    /* Find extremities of this edge */
    ip0 = mesh->edge[ie].a;
    ip1 = mesh->edge[ie].b;

    // TODO:: Multimaterial - Check whether an entity with reference ref should be split
    // Pb1: This function does not work here because we come from an edge and not a tetra
    // Need to find a way to know if we need to split or not this edge...
    // Pb2: even if I find a way to know if this edge need to be split on this partition
    // still no way to know if it is split from another partition if we are in
    // the case where on partition 1/material 1 and partition 2/material 2 and
    // we split mat1 but we do not split mat2. In that case, on partition 2, we will never know
    // that edges on // interface need to be split.
    // if ( !MMG5_isSplit(mesh,pt->ref,&refint,&refext) ) continue;

    /* STEP 5.1.1 - Create a new point if this edge needs to be split */
    /* Check the ls value at the edge nodes */
    p0 = &mesh->point[ip0];
    p1 = &mesh->point[ip1];
    v0 = sol->m[ip0];
    v1 = sol->m[ip1];

    /* Check if the edge should be split */
    /* If one of the points is exactly on the level set, the point exists already, pass */
    if ( fabs(v0) < MMG5_EPSD2 || fabs(v1) < MMG5_EPSD2 ) continue;
    /* If the points have the same sign, no need to split the edge, pass */
    else if ( MG_SMSGN(v0,v1) ) continue;
    /* If one or the other point has never been treated, pass */
    else if ( !p0->flag || !p1->flag ) continue;

    /* Define the weighting factor */
    s = v0 / (v0-v1);
    s = MG_MAX(MG_MIN(s,1.0-MMG5_EPS),MMG5_EPS);

    /* Find the coordinates of the new points using the weighting factor */
    c[0] = p0->c[0] + s*(p1->c[0]-p0->c[0]);
    c[1] = p0->c[1] + s*(p1->c[1]-p0->c[1]);
    c[2] = p0->c[2] + s*(p1->c[2]-p0->c[2]);

    /* Create a new point with coordinates c, tags MG_PARBDY+MG_NOSURF+MG_REQ
        and source src. Return the new number of points np in this partition    */
#ifdef USE_POINTMAP
    src = p0->src;
#else
    src = 1;
#endif
    np = MMG3D_newPt(mesh,c,MG_PARBDY+MG_NOSURF+MG_REQ,src);

    /* Update internal communicators of node and edge. For int edge comm
       intvalues stores the point position as in node2int_node_comm_index2 */
    grp->node2int_node_comm_index1[nitem_grp_node]=np; // Add this new point in int node comm index1
    grp->node2int_node_comm_index2[nitem_grp_node]=nitem_int_node; // Position in int node comm
    int_edge_comm->intvalues[pos] = nitem_int_node; // In int edge comm, assign position of node in int comm
    nitem_int_node += 1;
    nitem_grp_node += 1;
    grp->nitem_int_node_comm = nitem_grp_node;
    int_node_comm->nitem     = nitem_int_node;

    /* STEP 5.1.2 - Update hash table, met, sol and field for the new point */
    /* Memory allocation for sol and met */
    if ( !np ) {
      MMG5_int oldnpmax = mesh->npmax;
      MMG3D_POINT_REALLOC(mesh,sol,np,MMG5_GAP,
                            fprintf(stderr,"\n  ## Error: %s: unable to"
                                    " allocate a new point\n",__func__);
                            MMG5_INCREASE_MEM_MESSAGE();
                            return 0
                            ,c,0,src);
      if( met ) {
        if( met->m ) {
          MMG5_ADD_MEM(mesh,(met->size*(mesh->npmax-met->npmax))*sizeof(double),
                        "larger solution",
                        MMG5_SAFE_RECALLOC(mesh->point,mesh->npmax+1,oldnpmax+1,MMG5_Point,,);
                        mesh->memCur -= (mesh->npmax - oldnpmax)*sizeof(MMG5_Point);
                        mesh->npmax = oldnpmax;
                        mesh->np = mesh->npmax-1;
                        mesh->npnil = 0;
                        return 0);
          MMG5_SAFE_REALLOC(met->m,met->size*(met->npmax+1),
                            met->size*(mesh->npmax+1),
                            double,"larger solution",
                            MMG5_SAFE_RECALLOC(mesh->point,mesh->npmax+1,oldnpmax+1,MMG5_Point,,);
                            mesh->memCur -= (mesh->npmax - oldnpmax)*sizeof(MMG5_Point);
                            mesh->npmax = oldnpmax;
                            mesh->np = mesh->npmax-1;
                            mesh->npnil = 0;
                            return 0);
        }
        met->npmax = mesh->npmax;
      }
    }

    /* For this new point, add the value of the solution, i.e. the isovalue 0 */
    sol->m[np] = 0;

    /* If user provide a metric, interpolate it at the new point */
    if ( met && met->m ) {
      if ( met->size > 1 ) {
        ier = MMG3D_intmet33_ani_edge(met,ip0,ip1,np,s);
      }
      else {
        ier = MMG5_intmet_iso_edge(met,ip0,ip1,np,s);
      }
      if ( ier <= 0 ) {
        /* Unable to compute the metric */
        fprintf(stderr,"\n  ## Error: %s: unable to"
                " interpolate the metric during the level-set"
                " discretization\n",__func__);
        return 0;
      }
    }

    /* If user provide fields, interpolate them at the new point */
    if ( mesh->nsols ) {
      for ( j=0; j<mesh->nsols; ++j ) {
        psl    = field + j;
        if ( field->size > 1 ) {
          ier = MMG3D_intmet33_ani_edge(psl,ip0,ip1,np,s);
        }
        else {
          ier = MMG5_intmet_iso_edge(psl,ip0,ip1,np,s);
        }
        if ( ier <= 0 ) {
          /* Unable to compute fields */
          fprintf(stderr,"\n  ## Error: %s: unable to"
                  " interpolate fields during the level-set"
                  " discretization\n",__func__);
          return 0;
        }
      }
    }

    /* Update hash table */
    MMG5_hashUpdate(&hash,ip0,ip1,np);

  }

  /** STEP 5.2 -  Update external node communicator  */
  /* Loop over the external edge comm */
  for (i=0; i < next_edge_comm; i++) {

    /* Get external edge communicator information */
    ext_edge_comm  = &parmesh->ext_edge_comm[i]; // External edge communicator
    color_in_edge  = ext_edge_comm->color_in;    // Color of the hosting proc - this proc
    color_out_edge = ext_edge_comm->color_out;   // Color of the remote  proc - the proc to exchange with
    nitem_ext_edge = ext_edge_comm->nitem;       // Nbr of edges in common between these 2 procs

    /* Loop over the edges in the external edge communicator */
    for (j=0; j < nitem_ext_edge; j++) {

      /* Get the position of the edge and node in internal communicators */
      pos_edge = ext_edge_comm->int_comm_index[j];
      pos_node = int_edge_comm->intvalues[pos_edge];

      /* If pos_node < 0, this edge j is not split, so ignore it */
      if (pos_node < 0) continue;

      /* Update the external node communicator */
      /* Loop over the external node comm to find the appropriate one to be updated */
      for (k=0; k < next_node_comm; k++) {
        ext_node_comm  = &parmesh->ext_node_comm[k]; // External node communicator
        color_in_node  = ext_node_comm->color_in;    // Color of the hosting proc - this proc
        color_out_node = ext_node_comm->color_out;   // Color of the remote  proc - the proc to exchange with
        assert(color_in_node == color_in_edge);      // Ensure that the hosting proc is the same

        /* While color_out_node and color_out_edge are different, continue */
        if (color_out_node != color_out_edge) continue;

        /* If color_out of edge and node comm are the same - Update external node communicator */
        nitem_ext_node = ext_node_comm->nitem;                    // Initial nbr of nodes in common between these 2 procs
        ext_node_comm->int_comm_index[nitem_ext_node] = pos_node; // Add the node to the external node comm
        ext_node_comm->nitem = nitem_ext_node + 1;                // Updated nbr of nodes in common between these 2 procs
        break;
      }
    }
  }

  /** STEP 5.3 - Create all the other new points located elsewhere and update hash table */
  /* Loop over tetra k */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    /* Loop over the edges ia */
    for (ia=0; ia<6; ia++) {
      /* Find extremities of this edge and np (value stored in hash.item[key].k) */
      ip0 = pt->v[MMG5_iare[ia][0]];
      ip1 = pt->v[MMG5_iare[ia][1]];
      np  = MMG5_hashGet(&hash,ip0,ip1);

      /* If np>0 (i.e. hash.item[key].k != [0;-1]), this edge has already been split, pass to the next edge */
      if ( np>0 ) continue;

      /* Check whether an entity with reference ref should be split */
      // For now, the option nosplit is not supported
      // So we assume all the elements should be split
      // if ( !MMG5_isSplit(mesh,pt->ref,&refint,&refext) ) continue;

      /* STEP 5.3.1 - Create a new point if this edge needs to be split */
      /* Check the ls value at the edge nodes */
      p0 = &mesh->point[ip0];
      p1 = &mesh->point[ip1];
      v0 = sol->m[ip0];
      v1 = sol->m[ip1];

      /* Check if the edge should be split */
      /* If one of the points is exactly on the level set, the point exists already, pass */
      if ( fabs(v0) < MMG5_EPSD2 || fabs(v1) < MMG5_EPSD2 ) continue;
      /* If the points have the same sign, no need to split the edge, pass */
      else if ( MG_SMSGN(v0,v1) ) continue;
      /* If one or the other point has never been treated, pass */
      else if ( !p0->flag || !p1->flag ) continue;

      /* If np is = -1; then npneg is = 1  */
      npneg = (np<0);

      /* Define the weighting factor */
      s = v0 / (v0-v1);
      s = MG_MAX(MG_MIN(s,1.0-MMG5_EPS),MMG5_EPS);

      /* Find the coordinates of the new points using the weighting factor */
      c[0] = p0->c[0] + s*(p1->c[0]-p0->c[0]);
      c[1] = p0->c[1] + s*(p1->c[1]-p0->c[1]);
      c[2] = p0->c[2] + s*(p1->c[2]-p0->c[2]);

      /* Create a new point with coordinates c, tag 0 and source src.
         Return the new number of points np in this partition  */
#ifdef USE_POINTMAP
      src = p0->src;
#else
      src = 1;
#endif
      np = MMG3D_newPt(mesh,c,0,src);

      /* STEP 5.3.2 - Update of met, sol and field for the new point */
      /* Memory allocation for sol and met */
      if ( !np ) {
        MMG5_int oldnpmax = mesh->npmax;
        MMG3D_POINT_REALLOC(mesh,sol,np,MMG5_GAP,
                             fprintf(stderr,"\n  ## Error: %s: unable to"
                                     " allocate a new point\n",__func__);
                             MMG5_INCREASE_MEM_MESSAGE();
                             return 0
                             ,c,0,src);
        if( met ) {
          if( met->m ) {
            MMG5_ADD_MEM(mesh,(met->size*(mesh->npmax-met->npmax))*sizeof(double),
                         "larger solution",
                         MMG5_SAFE_RECALLOC(mesh->point,mesh->npmax+1,oldnpmax+1,MMG5_Point,,);
                         mesh->memCur -= (mesh->npmax - oldnpmax)*sizeof(MMG5_Point);
                         mesh->npmax = oldnpmax;
                         mesh->np = mesh->npmax-1;
                         mesh->npnil = 0;
                         return 0);
            MMG5_SAFE_REALLOC(met->m,met->size*(met->npmax+1),
                              met->size*(mesh->npmax+1),
                              double,"larger solution",
                              MMG5_SAFE_RECALLOC(mesh->point,mesh->npmax+1,oldnpmax+1,MMG5_Point,,);
                              mesh->memCur -= (mesh->npmax - oldnpmax)*sizeof(MMG5_Point);
                              mesh->npmax = oldnpmax;
                              mesh->np = mesh->npmax-1;
                              mesh->npnil = 0;
                              return 0);
          }
          met->npmax = mesh->npmax;
        }
      }

      /* For this new point, add the value of the solution, i.e. the isovalue 0 */
      sol->m[np] = 0;

      /* If user provide a metric, interpolate it at the new point */
      if ( met && met->m ) {
        if ( met->size > 1 ) {
          ier = MMG3D_intmet33_ani(mesh,met,k,ia,np,s);
        }
        else {
          ier = MMG5_intmet_iso(mesh,met,k,ia,np,s);
        }
        if ( ier <= 0 ) {
          /* Unable to compute the metric */
          fprintf(stderr,"\n  ## Error: %s: unable to"
                  " interpolate the metric during the level-set"
                  " discretization\n",__func__);
          return 0;
        }
      }

      /* If user provide fields, interpolate them at the new point */
      if ( mesh->nsols ) {
        for ( j=0; j<mesh->nsols; ++j ) {
          psl    = field + j;
          if ( field->size > 1 ) {
            ier = MMG3D_intmet33_ani(mesh,psl,k,ia,np,s);
          }
          else {
            ier = MMG5_intmet_iso(mesh,psl,k,ia,np,s);
          }
          if ( ier <= 0 ) {
            /* Unable to compute fields */
            fprintf(stderr,"\n  ## Error: %s: unable to"
                    " interpolate fields during the level-set"
                    " discretization\n",__func__);
            return 0;
          }
        }
      }

      /* STEP 5.2.3 - Update edge hash table */
      /* If this edge is required, then inform the user it is split anyway
         and update the hash: hash.item[key].k = - 1 becomes = np */
      if ( npneg ) {
        /* We split a required edge */
        if ( !mmgWarn ) {
          mmgWarn = 1;
          if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
            fprintf(stderr,"  ## Warning: %s: the level-set intersect at least"
                    " one required entity. Required entity ignored.\n\n",__func__);
          }
        }
        MMG5_hashUpdate(&hash,ip0,ip1,np);
      }
      /* Otherwise add the edge to be split into hash table */
      else {
        MMG5_hashEdge(mesh,&hash,ip0,ip1,np);
      }
    }
  }

#ifndef NDEBUG
  /** Check the internal node communicator */
  assert( PMMG_check_intNodeComm( parmesh ) );

  /** Check the external node communicator */
  assert( PMMG_check_extNodeComm( parmesh,parmesh->info.read_comm ) );
#endif

  /** STEP 6 - Split according to tets flags */
  /** STEP 6.1 - Compute global node vertices */
  if ( !PMMG_Compute_verticesGloNum( parmesh,parmesh->comm ) ) {
    fprintf(stderr,"\n\n\n  -- WARNING: IMPOSSIBLE TO COMPUTE NODE GLOBAL NUMBERING\n\n\n");
    return 0;
  }

  /** STEP 6.2 - Do the splitting for tetra on parallel interface */
  ns  = 0;     // Number of total split on this proc
  ier = 1;     // Error
  idx_tmp = 0; // Index of an already split tetra to recover info stored in ne_tmp_tab and vGlobNum_tab
  nitem_grp_face_tmp = nitem_grp_face;

  /* Allocate internal face comm */
  PMMG_CALLOC(parmesh,int_face_comm->intvalues,3*nitem_int_face,int,"int_face_comm intvalues",return 0);

  /* Loop over the internal faces communicator */
  for (i=0; i < nitem_grp_face; i++) {

    /* Get the position of the face in internal communicator and the value of the face */
    val_face = grp->face2int_face_comm_index1[i];
    pos      = grp->face2int_face_comm_index2[i];

    /* Initialize interval face comm at -1 */
    int_face_comm->intvalues[3*pos]   = -1;
    int_face_comm->intvalues[3*pos+1] = -1;
    int_face_comm->intvalues[3*pos+2] = -1;

    /* Find the local tetra, the face and node associated */
    k     = val_face/12;     // Index of the tetra on this proc
    ifac  = (val_face%12)/3; // Index of the face
    iploc = (val_face%12)%3; // Index of the node

    /* Get the tetra k and xtetra associated */
    pt  = &mesh->tetra[k];
    assert ( MG_EOK(pt) && "Invalid tetra stored in the communicator" );

    pxt = &mesh->xtetra[pt->xt];

    /* STEP 6.2.1 - Find global numbering and split pattern of the tetra */
    /* If the tetra has already a flag (flag !=0) - then it has already been processed */
    already_split = 0;
    if (pt->flag) {
      /* If flag>0, we need to update the face communicators; */
      if (pt->flag > 0) {
        already_split = 1;             // Identify the tetra as already split
        idx_tmp = pt->mark;            // Index of ne_tmp stored in ne_tmp_tab
        ne_tmp  = ne_tmp_tab[idx_tmp]; // Old total number of tetras just after the split of this old tetra
        /* Get global num of the old tetra vertices before the split */
        for (j=0; j<4; j++) {
          vGlobNum[j] = vGlobNum_tab[idx_tmp*4+j];
        }
      }
      /* Otherwise, if flag equals to -1 (flag<0), the tetra does not need to be split. Pass to the next tetra. */
      else if (pt->flag < 0) continue;
    }
    /* Otherwise, if flag=0, the tetra has never been processed. */
    else {
      /* Get the split pattern: loop over the edges, get hash.item[key].k */
      memset(vx,0,6*sizeof(MMG5_int));
      for (ia=0; ia<6; ia++) {
        vx[ia] = MMG5_hashGet(&hash,pt->v[MMG5_iare[ia][0]],pt->v[MMG5_iare[ia][1]]);
        if ( vx[ia] > 0 )
          MG_SET(pt->flag,ia);
      }

      /* Get and store global num of the tetra vertices */
      for (j=0; j<4; j++) {
        vGlobNum[j] = mesh->point[pt->v[j]].tmp;
        vGlobNum_tab[ns*4+j] = vGlobNum[j];
      }
    }

    /* Get the split pattern stored in flag */
    flag = pt->flag;

    /* Initialize tetra_sorted and node_sorted at -1 */
    for ( j=0; j<3; ++j ) {
      tetra_sorted[j] = -1;
      node_sorted[j] = -1;
    }

    /* STEP 6.2.2 - If not already done, split the tetra according to the flag */
    switch (flag) {
    case 1: case 2: case 4: case 8: case 16: case 32: // 1 edge split
      if (!already_split) {
        ier = MMG5_split1(mesh,met,k,vx,1);
        pt->flag = flag;           // Re-flag tetra k as the flag has been reset in the split
        pt->mark = ns;             // Split number to recover info later if needed
        ne_tmp_tab[ns] = mesh->ne; // Total number of tetras after this split
        ne_tmp         = mesh->ne; // Total number of tetras after this split
        ns++;                      // Incremente the total number of split
      }

      /* Find the tetras and nodes defining the face ifac of tetra k */
      MMG3D_split1_cfg(flag,tau,&taued); // Compute tau
      PMMG_split1_sort(mesh,k,ifac,tau,ne_tmp,tetra_sorted,node_sorted); // Find tetra_sorted and node_sorted
      break;

    case 48: case 24: case 40: case 6: case 34: case 36: // 2 edges (same face) split
    case 20: case 5:  case 17: case 9: case 3:  case 10:
      if (!already_split) {
        ier = MMG5_split2sf_globNum(mesh,met,k,vx,vGlobNum,1);
        pt->flag = flag;           // Re-flag tetra k as the flag has been reset in the split
        pt->mark = ns;             // Split number to recover info later if needed
        ne_tmp_tab[ns] = mesh->ne; // Total number of tetras after this split
        ne_tmp         = mesh->ne; // Total number of tetras after this split
        ns++;                      // Incremente the total number of split
      }

      /* Find the tetras and nodes defining the face ifac of tetra k */
      MMG3D_split2sf_cfg(flag,vGlobNum,tau,&taued,&imin0); // Compute tau and imin0
      PMMG_split2sf_sort(mesh,k,ifac,tau,imin0,ne_tmp,tetra_sorted,node_sorted); // Find tetra_sorted and node_sorted
      break;

    case 7: case 25: case 42: case 52: // 3 edges on conic configuration split
      if (!already_split) {
        ier = MMG5_split3cone_globNum(mesh,met,k,vx,vGlobNum,1);
        pt->flag = flag;           // Re-flag tetra k as the flag has been reset in the split
        pt->mark = ns;             // Split number to recover info later if needed
        ne_tmp_tab[ns] = mesh->ne; // Total number of tetras after this split
        ne_tmp         = mesh->ne; // Total number of tetras after this split
        ns++;                      // Incremente the total number of split
      }

      /* Find the tetras and nodes defining the face ifac of tetra k */
      MMG3D_split3cone_cfg(flag,vGlobNum,tau,&taued,&imin0,&imin2); // Compute tau, imin0 and imin2
      PMMG_split3cone_sort(mesh,k,ifac,tau,imin0,imin2,ne_tmp,tetra_sorted,node_sorted); // Find tetra_sorted and node_sorted
      break;

    case 30: case 45: case 51: // 4 edges on opposite configuration split
      if (!already_split) {
        ier = MMG5_split4op_globNum(mesh,met,k,vx,vGlobNum,1);
        pt->flag = flag;           // Re-flag tetra k as the flag has been reset in the split
        pt->mark = ns;             // Split number to recover info later if needed
        ne_tmp_tab[ns] = mesh->ne; // Total number of tetras after this split
        ne_tmp         = mesh->ne; // Total number of tetras after this split
        ns++;                      // Incremente the total number of split
      }

      /* Find the tetras and nodes defining the face ifac of tetra k */
      MMG3D_split4op_cfg(flag,vGlobNum,tau,&taued,&imin0,&imin2); // Compute tau, imin0 and imin2
      PMMG_split4op_sort(mesh,k,ifac,tau,imin0,imin2,ne_tmp,tetra_sorted,node_sorted); // Find tetra_sorted and node_sorted
      break;

    default: // This tetra does not need to be split and is processed for the first time
      assert(pt->flag == 0);
      pt->flag = -1; // Put this flag to -1 to specify that this tetra has been processed
      PMMG_nosplit_sort(mesh,k,ifac,tetra_sorted,node_sorted);
      break;
    }

    if ( !ier ) return 0;

    if (pt->flag != -1) {
      /* STEP 6.2.3 - Update tag of edges in xtetra with MG_PARBDY */
      for (j=0; j<3; j++) {
        k0 = tetra_sorted[j];
        if (k0 != -1) {
          pt0  = &mesh->tetra[k0];
          pxt0 = &mesh->xtetra[pt0->xt];
          for (k=0; k<3; k++) {
            ia = MMG5_iarf[ifac][k];
            pxt0->tag[ia] |= MG_PARBDY;
          }
        }
      }

      /* STEP 6.2.4 - Update internal face communicators */
      /* (a) Update the first face located at i - Modify only index1 - index2 stays the same */
      grp->face2int_face_comm_index1[i] = 12*tetra_sorted[0]+3*ifac+node_sorted[0];
      int_face_comm->intvalues[3*pos]   = pos;

      /* (b) Update the communicators for the potential 2 other faces */
      nface_added = 0;
      for (j=0; j<2; j++) {
        if ( tetra_sorted[j+1] != -1) {
          grp->face2int_face_comm_index1[nitem_grp_face_tmp+j] = 12*tetra_sorted[j+1]+3*ifac+node_sorted[j+1];
          grp->face2int_face_comm_index2[nitem_grp_face_tmp+j] = nitem_grp_face_tmp+j;
          int_face_comm->intvalues[3*pos+j+1] = nitem_grp_face_tmp+j;
          nface_added += 1;
        }
      }

      /* (c) Update the total number of faces */
      nitem_grp_face_tmp += nface_added;
      int_face_comm->nitem = nitem_grp_face_tmp;
    }
    else if (pt->flag == -1) {
      /* STEP 6.2.5 - Update internal face communicators for tetra not split */
      /* As the tetra is not split, the tetra index has not changed.
          Moreover, the node of face ifac is chosen to be the node with highest
          coordinates, so it should not have been changed either. However, to be
          sure we are not missing a case, we still update face2int_face_comm_index1 */
      grp->face2int_face_comm_index1[i] = 12*tetra_sorted[0]+3*ifac+node_sorted[0];
    }
  }

  /** STEP 6.3 - Update internal and external face communicator */
  /* Update number of element in internal face comm*/
  nitem_grp_face = nitem_grp_face_tmp;
  grp->nitem_int_face_comm = nitem_grp_face;
  int_face_comm->nitem = nitem_grp_face;

  /* Update external face comm */
  for (k=0; k < next_face_comm; k++) {

    /* Get current external face communicator */
    ext_face_comm  = &parmesh->ext_face_comm[k]; // External face communicator
    nitem_ext_face = ext_face_comm->nitem;       // Nbr of faces in common between these 2 procs
    nitem_ext_face_init = ext_face_comm->nitem;  // Initial nbr of faces in common between these 2 procs

    /* Loop over the faces in the external face communicator */
    for (i=0; i < nitem_ext_face_init; i++) {
      pos = ext_face_comm->int_comm_index[i];

      /* Loop over the potential 3 faces created after the split */
      for (j=0; j < 3; j++) {
        pos_face = int_face_comm->intvalues[3*pos+j];
        if (pos_face < 0) continue; // If pos_face=-1, there is not extra face to add
        /* The first face is located at i in the ext face comm */
        if (j==0) {
          ext_face_comm->int_comm_index[i] = pos_face;
        }
        /* The next faces are added at the end of the ext face comm */
        else {
          ext_face_comm->int_comm_index[nitem_ext_face] = pos_face;
          nitem_ext_face += 1;
        }
      }
    }
    ext_face_comm->nitem = nitem_ext_face; // Update nbr of face in ext face comm
  }

#ifndef NDEBUG
  /** Check the internal face communicator */
  assert( PMMG_check_intFaceComm( parmesh ) );

  /** Check the external face communicator */
  assert( PMMG_check_extFaceComm( parmesh,parmesh->info.read_comm ) );
#endif

  /** STEP 6.4 - Do the splitting for tetra located elsewhere */
  /* Loop over tetra */
  for (k=1; k<=ne_init; k++) {

    /* Get the tetra k and xtetra associated */
    pt  = &mesh->tetra[k];
    pxt = &mesh->xtetra[pt->xt];
    if ( !MG_EOK(pt) )  continue;

    /* If the tetra has already a flag (flag !=0) - it has already been processed. Pass to the next tetra. */
    if (pt->flag) continue;

    /* STEP 6.4.1 - Find global numbering and split pattern of the tetra */
    /* Get the split pattern: loop over the edges, get hash.item[key].k */
    memset(vx,0,6*sizeof(MMG5_int));
    for (ia=0; ia<6; ia++) {
      vx[ia] = MMG5_hashGet(&hash,pt->v[MMG5_iare[ia][0]],pt->v[MMG5_iare[ia][1]]);
      if ( vx[ia] > 0 )
        MG_SET(pt->flag,ia);
    }

    /* Get global num of the tetra vertices */
    for (j=0; j<4; j++) {
      vGlobNum[j] = mesh->point[pt->v[j]].tmp;
    }

    /* Get the split pattern stored in flag */
    flag = pt->flag;

    /* STEP 6.4.2 - If not already done, split the tetra according to the flag */
    switch (flag) {
    case 1: case 2: case 4: case 8: case 16: case 32: // 1 edge split
      ier = MMG5_split1(mesh,met,k,vx,1);
      pt->flag = flag; // Re-flag tetra k as the flag has been reset in the split
      ns++;            // Incremente the total number of split
      break;

    case 48: case 24: case 40: case 6: case 34: case 36: // 2 edges (same face) split
    case 20: case 5:  case 17: case 9: case 3:  case 10:
      ier = MMG5_split2sf_globNum(mesh,met,k,vx,vGlobNum,1);
      pt->flag = flag; // Re-flag tetra k as the flag has been reset in the split
      ns++;            // Incremente the total number of split
      break;

    case 7: case 25: case 42: case 52: // 3 edges on conic configuration split
      ier = MMG5_split3cone_globNum(mesh,met,k,vx,vGlobNum,1);
      pt->flag = flag; // Re-flag tetra k as the flag has been reset in the split
      ns++;            // Incremente the total number of split
      break;

    case 30: case 45: case 51: // 4 edges on opposite configuration split
      ier = MMG5_split4op_globNum(mesh,met,k,vx,vGlobNum,1);
      pt->flag = flag; // Re-flag tetra k as the flag has been reset in the split
      ns++;            // Incremente the total number of split
      break;

    default:
      assert(pt->flag == 0);
      pt->flag = -1; // Put this flag to -1 to specify that this tetra has been processed
      break;
    }

    if ( !ier ) return 0;
  }

  /** STEP 7 - Deallocation/Allocation of memory and reset of some fields */
  /* Delete the tables storing imin0, imin2 and ne_tmp_tab */
  PMMG_DEL_MEM(parmesh,vGlobNum_tab,MMG5_int,"vGlobNum_tab");
  PMMG_DEL_MEM(parmesh,ne_tmp_tab,MMG5_int,"ne_tmp_tab");

  /* Delete the edges hash table */
  PMMG_DEL_MEM(parmesh,int_edge_comm->intvalues,int,"edge intvalues");
  PMMG_DEL_MEM(parmesh,int_face_comm->intvalues,int,"face intvalues");

  /* Delete internal communicators */
  MMG5_DEL_MEM(mesh,hash.item);

  /* Realloc internal node communicators to exact final size */
  PMMG_REALLOC(parmesh, grp->node2int_node_comm_index1,
               nitem_grp_node,
               nitem_grp_node_firstalloc,
               int,"Allocation of node2int_node_comm_index1", return 0);
  PMMG_REALLOC(parmesh, grp->node2int_node_comm_index2,
               nitem_grp_node,
               nitem_grp_node_firstalloc,
               int,"Allocation of node2int_node_comm_index2", return 0);

  /* Realloc internal face communicators to exact final size */
  PMMG_REALLOC(parmesh, grp->face2int_face_comm_index1,
               nitem_grp_face,
               nitem_grp_face_firstalloc,
               int,"Allocation of face2int_face_comm_index1", return 0);
  PMMG_REALLOC(parmesh, grp->face2int_face_comm_index2,
               nitem_grp_face,
               nitem_grp_face_firstalloc,
               int,"Allocation of face2int_face_comm_index2", return 0);

  /* Realloc external node communicator to exact final size */
  for ( k=0; k<parmesh->next_node_comm; k++ ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    PMMG_REALLOC(parmesh,ext_node_comm->int_comm_index,
                 ext_node_comm->nitem,
                 ext_node_comm->nitem_to_share,
                 int,"Re-allocation of external node communicator",return 0);
    ext_node_comm->nitem_to_share = 0;
  }

  /* Realloc external face communicator to exact final size */
  for ( k=0; k<parmesh->next_face_comm; k++ ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    PMMG_REALLOC(parmesh,ext_face_comm->int_comm_index,
                 ext_face_comm->nitem,
                 ext_face_comm->nitem_to_share,
                 int,"Re-allocation of external face communicator",return 0);
    ext_face_comm->nitem_to_share = 0;
  }

  /* Reset mark and flag in mesh->tetra */
  for (k=1; k<=ne_init; k++) {
    mesh->tetra[k].mark = 0;
    mesh->tetra[k].flag = 0;
  }

  return ns;
}

/**
 * \param mesh   pointer toward the mesh structure
 * \param k      index of the tetra that we do not split
 * \param ifac   local index of the face located on a parallel boundary
 * \param tetra_sorted indices of tetra
 *                     sorted by increasing order of their global node index
 * \param node_sorted  for each tetras in tetra_sorted: local index of the node
 *                     on ifac having the minimum global node index
 *
 * Find the node on ifac with the minimum global node index
 * for the tetra that we do not split and store the local index in node_sorted.
 *
 */
void PMMG_nosplit_sort(MMG5_pMesh mesh,MMG5_int k,int ifac,MMG5_int *tetra_sorted,MMG5_int *node_sorted) {
  MMG5_int v_t0[3], v_t1[3], v_t2[3];

  /* STEP 1 - Find the indices of the tetra */
  /* The tetra k is not split
     Note that 2 faces of the initial tetra are divided into 2
           and 2 faces are not divided.                        */

  /* STEP 1.1 - Index of the first tetra */
  /* Tetra #0 created by MMG5_split1 */
  tetra_sorted[0] = k;

  /* Store index of the node with highest coordinates in node_sorted[0]
     and sort the vertices by increasing order in v_t0 */
  node_sorted[0]=PMMG_sort_vertices(mesh,tetra_sorted[0],v_t0,ifac);

  /* STEP 1.2 - Index of the second tetra */
  tetra_sorted[1] = -1;
  v_t1[0] = v_t1[1] = v_t1[2] = -1;

  /* STEP 1.3 - Index of the third tetra */
  /* There is no third tetra created by MMG5_split1 - we assign -1 */
  tetra_sorted[2] = -1;
  v_t2[0] = v_t2[1] = v_t2[2] = -1;

  /* STEP 2 - Sort these tetras by their global indices */
  PMMG_sort_tetra(tetra_sorted,node_sorted,v_t0,v_t1,v_t2);

  return;
}

/**
 * \param mesh   pointer toward the mesh structure
 * \param k      index of the tetra that we have split with MMG5_split1
 * \param ifac   local index of the face located on a parallel boundary
 *               that we have split with MMG5_split1
 * \param tau    vertices permutation
 * \param ne_tmp number of tetra created after MMG5_split1
 * \param tetra_sorted indices of tetras (defining ifac after the split MMG5_split1)
 *                     sorted by increasing order of their global node index
 * \param node_sorted  for each tetras in tetra_sorted: local index of the node
 *                     on ifac having the minimum global node index
 *
 * Sort the tetras created by MMG5_split1 defining the face ifac.
 * Find the node on ifac with the minimum global node index
 * for each tetra in tetra_sorted and store the local index in node_sorted.
 *
 */
void PMMG_split1_sort(MMG5_pMesh mesh,MMG5_int k,int ifac,uint8_t tau[4],
                      MMG5_int ne_tmp,MMG5_int *tetra_sorted,MMG5_int *node_sorted) {
  MMG5_int v_t0[3], v_t1[3], v_t2[3];

  /* STEP 1 - Find the indices of the new tetras defining the face ifac */
  /* The 2 tetras created by MMG5_split1 are at
           mesh.tetra[k] (tetra #0) and [ne_tmp] (tetra #1)
     Note that 2 faces of the initial tetra are divided into 2
           and 2 faces are not divided.                        */

  /* STEP 1.1 - Index of the first tetra */
  /* Tetra #0 created by MMG5_split1 */
  tetra_sorted[0] = k;
  /* Except for the following: treta #1 created by MMG5_split1 */
  if ( ifac==tau[0] ) tetra_sorted[0] = ne_tmp;

  /* Store index of the node with highest coordinates in node_sorted[0]
     and sort the vertices by increasing order in v_t0 */
  node_sorted[0]=PMMG_sort_vertices(mesh,tetra_sorted[0],v_t0,ifac);

  /* STEP 1.2 - Index of the second tetra */
  /* Tetra #1 created by MMG5_split1 */
  tetra_sorted[1] = ne_tmp;
  /* Except for the following: no more tetra define ifac - we assign -1 */
  if ( ifac == tau[0] ) tetra_sorted[1] = -1;
  if ( ifac == tau[1] ) tetra_sorted[1] = -1;

  if ( tetra_sorted[1] != -1 ) {
    /* Store index of the node with highest coordinates in node_sorted[1]
      and sort the vertices by increasing order in v_t1 */
    node_sorted[1]=PMMG_sort_vertices(mesh,tetra_sorted[1],v_t1,ifac);
  }
  else {
    v_t1[0] = v_t1[1] = v_t1[2] = -1;
  }

  /* STEP 1.3 - Index of the third tetra */
  /* There is no third tetra created by MMG5_split1 - we assign -1 */
  tetra_sorted[2] = -1;
  v_t2[0] = v_t2[1] = v_t2[2] = -1;

  /* STEP 2 - Sort these tetras by their global indices */
  PMMG_sort_tetra(tetra_sorted,node_sorted,v_t0,v_t1,v_t2);

  return;
}

/**
 * \param mesh   pointer toward the mesh structure
 * \param k      index of the tetra that we have split with MMG5_split2sf_globNum
 * \param ifac   local index of the face located on a parallel boundary
 *               that we have split with MMG5_split2sf_globNum
 * \param tau    vertices permutation
 * \param imin   minimal index of vertices \a tau[1] and \a tau[2]
 * \param ne_tmp number of tetra created after MMG5_split2sf_globNum
 * \param tetra_sorted indices of tetras (defining ifac after the split MMG5_split2sf_globNum)
 *                     sorted by increasing order of their global node index
 * \param node_sorted  for each tetras in tetra_sorted: local index of the node on ifac having
 *                     the minimum global node index
 *
 * Sort the tetras created by MMG5_split2sf_globNum defining the face ifac.
 * Find the node on ifac with the minimum global node index
 * for each tetra in tetra_sorted and store the local index in node_sorted.
 *
 */
void PMMG_split2sf_sort(MMG5_pMesh mesh,MMG5_int k,int ifac,uint8_t tau[4],int imin,
                        MMG5_int ne_tmp,MMG5_int *tetra_sorted,MMG5_int *node_sorted) {
  MMG5_int v_t0[3], v_t1[3], v_t2[3];

  /* STEP 1 - Find the indices of the new tetras defining the face ifac */
  /* The 3 tetras created by MMG3D_split2sf_globNum are at
        mesh.tetra[k] (tetra #0), [ne_tmp-1] (tetra #1) and [ne_tmp] (tetra #2)
     Note that 1 face of the initial tetra is divided into 3;
               2 faces are divided into 2 and 1 face is not divided.            */

  /* STEP 1.1 - Index of the first tetra */
  /* Tetra #0 created by MMG5_split2sf_globNum */
  tetra_sorted[0] = k;
  /* Except for the following: treta #1 or #2 created by MMG5_split2sf_globNum */
  if ( (imin==tau[1]) && (ifac==tau[3]) ) tetra_sorted[0] = ne_tmp;
  if ( (imin==tau[2]) && (ifac==tau[3]) ) tetra_sorted[0] = ne_tmp-1;

  /* Store index of the node with highest coordinates in node_sorted[0]
     and sort the vertices by increasing order in v_t0 */
  node_sorted[0]=PMMG_sort_vertices(mesh,tetra_sorted[0],v_t0,ifac);

  /* STEP 1.2 - Index of the second tetra */
  /* Tetra #1 created by MMG5_split2sf_globNum */
  tetra_sorted[1] = ne_tmp-1;
  /* Except for the following: treta #2 created by MMG5_split2sf_globNum or no more tetra */
  if ( ifac == tau[1] ) tetra_sorted[1] = ne_tmp;
  if ( ifac == tau[3] ) tetra_sorted[1] = -1;

  if ( tetra_sorted[1] != -1 ) {
    /* Store index of the node with highest coordinates in node_sorted[1]
      and sort the vertices by increasing order in v_t1 */
    node_sorted[1]=PMMG_sort_vertices(mesh,tetra_sorted[1],v_t1,ifac);
  }
  else {
    v_t1[0] = v_t1[1] = v_t1[2] = -1;
  }

  /* STEP 1.3 - Index of the third tetra */
  /* Tetra #2 created by MMG5_split2sf_globNum */
  tetra_sorted[2] = ne_tmp;
  /* Except for the following: no more tetra define ifac */
  if ( ifac != tau[0] ) tetra_sorted[2] = -1;

  if ( tetra_sorted[2] != -1 ) {
    /* Store index of the node with highest coordinates in node_sorted[2]
      and sort the vertices by increasing order in v_t2 */
    node_sorted[2]=PMMG_sort_vertices(mesh,tetra_sorted[2],v_t2,ifac);
  }
  else{
    v_t2[0] = v_t2[1] = v_t2[2] = -1;
  }

  /* STEP 2 - Sort these tetras by their global indices */
  PMMG_sort_tetra(tetra_sorted,node_sorted,v_t0,v_t1,v_t2);

  return;
}

/**
 * \param mesh   pointer toward the mesh structure
 * \param k      index of the tetra that we have split with MMG5_split3cone_globNum
 * \param ifac   local index of the face located on a parallel boundary
 *               that we have split with MMG5_split3cone_globNum
 * \param tau    vertices permutation
 * \param ia     first  condition to choose the appropriate split in MMG5_split3cone_globNum
 * \param ib     second condition to choose the appropriate split in MMG5_split3cone_globNum
 * \param ne_tmp number of tetra created after MMG5_split3cone_globNum
 * \param tetra_sorted indices of tetras (defining ifac after the split MMG5_split3cone_globNum)
 *                     sorted by increasing order of their global node index
 * \param node_sorted  for each tetras in tetra_sorted: local index of the node on ifac having
 *                     the minimum global node index
 *
 * Sort the tetras created by MMG5_split3cone_globNum defining the face ifac.
 * Find the node on ifac with the minimum global node index
 * for each tetra in tetra_sorted and store the local index in node_sorted.
 *
 */
void PMMG_split3cone_sort(MMG5_pMesh mesh, MMG5_int k,int ifac,uint8_t tau[4],int ia,int ib,
                          MMG5_int ne_tmp,MMG5_int *tetra_sorted,MMG5_int *node_sorted) {
  MMG5_int v_t0[3], v_t1[3], v_t2[3];

  /* STEP 1 - Find the indices of the new tetras defining the face ifac */
  /* The 4 tetras created by MMG3D_split3cone_globNum are at
           mesh.tetra[k] (tetra #0),  [ne_tmp-2] (tetra #1),
              [ne_tmp-1] (tetra #2) and [ne_tmp] (tetra #3)
     Note that 3 faces of the initial tetra are divided into 3
           and 1 face is not divided.                          */

  /* STEP 1.1 - Index of the first tetra */
  /* Tetra #0 created by MMG5_split3cone_globNum */
  tetra_sorted[0] = k;
  /* Except for the following: treta #3 created by MMG5_split3cone_globNum */
  if ( ifac == tau[0] ) tetra_sorted[0] = ne_tmp;

  /* Store index of the node with highest coordinates in node_sorted[0]
     and sort the vertices by increasing order in v_t0 */
  node_sorted[0]=PMMG_sort_vertices(mesh,tetra_sorted[0],v_t0,ifac);

  /* STEP 1.2 - Index of the second tetra */
  /* Tetra #1 created by MMG5_split3cone_globNum */
  tetra_sorted[1] = ne_tmp-2;
  /* Except for the following: treta #2 created by MMG5_split3cone_globNum or no more tetra */
  if ( (ia==tau[1]) && (ifac == tau[1]) ) tetra_sorted[1] = ne_tmp-1;
  if ( (ia==tau[2]) && (ifac == tau[2]) ) tetra_sorted[1] = ne_tmp-1;
  if ( (ia==tau[3]) && (ifac == tau[3]) ) tetra_sorted[1] = ne_tmp-1;
  if ( ifac == tau[0] ) tetra_sorted[1] = -1;

  if ( tetra_sorted[1] != -1 ) {
    /* Store index of the node with highest coordinates in node_sorted[1]
      and sort the vertices by increasing order in v_t1 */
    node_sorted[1]=PMMG_sort_vertices(mesh,tetra_sorted[1],v_t1,ifac);
  }
  else {
    v_t1[0] = v_t1[1] = v_t1[2] = -1;
  }

  /* STEP 1.3 - Index of the third tetra */
  /* Tetra #3 created by MMG5_split3cone_globNum */
  tetra_sorted[2] = ne_tmp;
  /* Except for the following: treta #2 by MMG5_split3cone_globNum or no more tetra */
  if ( (ia==tau[1]) && (ib==tau[2]) && (ifac == tau[3]) ) tetra_sorted[2] = ne_tmp-1;
  if ( (ia==tau[1]) && (ib==tau[3]) && (ifac == tau[2]) ) tetra_sorted[2] = ne_tmp-1;
  if ( (ia==tau[2]) && (ib==tau[1]) && (ifac == tau[3]) ) tetra_sorted[2] = ne_tmp-1;
  if ( (ia==tau[2]) && (ib==tau[3]) && (ifac == tau[1]) ) tetra_sorted[2] = ne_tmp-1;
  if ( (ia==tau[3]) && (ib==tau[1]) && (ifac == tau[2]) ) tetra_sorted[2] = ne_tmp-1;
  if ( (ia==tau[3]) && (ib==tau[2]) && (ifac == tau[1]) ) tetra_sorted[2] = ne_tmp-1;
  if ( ifac == tau[0] ) tetra_sorted[2] = -1;

  if ( tetra_sorted[2] != -1 ) {
    /* Store index of the node with highest coordinates in node_sorted[2]
      and sort the vertices by increasing order in v_t2 */
    node_sorted[2]=PMMG_sort_vertices(mesh,tetra_sorted[2],v_t2,ifac);
  }
  else{
    v_t2[0] = v_t2[1] = v_t2[2] = -1;
  }

  /* STEP 2 - Sort these tetras by their global indices */
  PMMG_sort_tetra(tetra_sorted,node_sorted,v_t0,v_t1,v_t2);

  return;
}


/**
 * \param mesh   pointer toward the mesh structure
 * \param k      index of the tetra that we have split with MMG5_split4op_globNum
 * \param ifac   local index of the face located on a parallel boundary
 *               that we have split with MMG5_split4op_globNum
 * \param tau    vertices permutation
 * \param imin01 minimal index of vertices \a tau[0] and \a tau[1]
 * \param imin23 minimal index of vertices \a tau[2] and \a tau[3]
 * \param ne_tmp number of tetra created after MMG5_split4op_globNum
 * \param tetra_sorted indices of tetras (defining ifac after the split MMG5_split4op_globNum)
 *                     sorted by increasing order of their global node index
 * \param node_sorted  for each tetras in tetra_sorted: local index of the node on ifac having
 *                     the minimum global node index
 *
 * Sort the tetras created by MMG5_split4op_globNum defining the face ifac.
 * Find the node on ifac with the minimum global node index
 * for each tetra in tetra_sorted and store the local index in node_sorted.
 *
 */
void PMMG_split4op_sort(MMG5_pMesh mesh,MMG5_int k,int ifac,uint8_t tau[4],int imin01,int imin23,
                        MMG5_int ne_tmp,MMG5_int *tetra_sorted,MMG5_int *node_sorted) {
  MMG5_int v_t0[3], v_t1[3], v_t2[3];

  /* STEP 1 - Find the indices of the new tetras defining the face ifac */
  /* The 6 tetras created by MMG5_split4op_globNum are at
        mesh.tetra[k] (tetra #0), [ne_tmp-4] (tetra #1),  [ne_tmp-3] (tetra #2),
           [ne_tmp-2] (tetra #3), [ne_tmp-1] (tetra #4) and [ne_tmp] (tetra #5)
     Note that all the 4 faces of the initial tetra are divided into 3.          */

  /* STEP 1.1 - Index of the first tetra */
  /* Tetra #0 created by MMG5_split4op_globNum */
  tetra_sorted[0] = k;
  /* Except for the following: treta #2 created by MMG5_split4op_globNum */
  if ( (imin01==tau[0]) && (imin23==tau[2]) && (ifac == tau[1]) ) tetra_sorted[0] = ne_tmp-3;
  if ( (imin01==tau[1]) && (imin23==tau[2]) && (ifac == tau[0]) ) tetra_sorted[0] = ne_tmp-3;
  if ( (imin01==tau[0]) && (imin23==tau[3]) && (ifac == tau[1]) ) tetra_sorted[0] = ne_tmp-3;
  if ( (imin01==tau[1]) && (imin23==tau[3]) && (ifac == tau[0]) ) tetra_sorted[0] = ne_tmp-3;

  /* Store index of the node with highest coordinates in node_sorted[0]
     and sort the vertices by increasing order in v_t0 */
  node_sorted[0]=PMMG_sort_vertices(mesh,tetra_sorted[0],v_t0,ifac);

  /* STEP 1.2 - Index of the second tetra */
  /* Tetra #3 created by MMG5_split4op_globNum */
  tetra_sorted[1] = ne_tmp-2;
  /* Except for the following: treta #1 or #2 created by MMG5_split4op_globNum */
  if (ifac == tau[2]) {
    tetra_sorted[1] = ne_tmp-4;
    if ( imin01 == tau[1]) tetra_sorted[1] = ne_tmp-3;
  }
  else if (ifac==tau[3]) {
    tetra_sorted[1] = ne_tmp-3;
    if ( imin01 == tau[1]) tetra_sorted[1] = ne_tmp-4;
  }

  /* Store index of the node with highest coordinates in node_sorted[1]
     and sort the vertices by increasing order in v_t1 */
  node_sorted[1]=PMMG_sort_vertices(mesh,tetra_sorted[1],v_t1,ifac);

  /* STEP 1.3 - Index of the third tetra */
  /* Tetra #5 created by MMG5_split4op_globNum */
  tetra_sorted[2] = ne_tmp;
  /* Except for the following: treta #3 or #4 created by MMG5_split4op_globNum */
  if ( (imin23==tau[2]) && (ifac == tau[0]) ) tetra_sorted[2] = ne_tmp-1;
  if ( (imin23==tau[3]) && (ifac == tau[1]) ) tetra_sorted[2] = ne_tmp-1;
  if ( (imin23==tau[2]) && (ifac == tau[2]) ) tetra_sorted[2] = ne_tmp-2;
  if ( (imin23==tau[3]) && (ifac == tau[3]) ) tetra_sorted[2] = ne_tmp-2;

  /* Store index of the node with highest coordinates in node_sorted[2]
     and sort the vertices by increasing order in v_t2 */
  node_sorted[2]=PMMG_sort_vertices(mesh,tetra_sorted[2],v_t2,ifac);

  /* STEP 2 - Sort these tetras by their global indices */
  PMMG_sort_tetra(tetra_sorted,node_sorted,v_t0,v_t1,v_t2);

  return;
}

/**
 * \param mesh     pointer toward the mesh structure
 * \param k        index of the tetra
 * \param v_t      table of a triangle vertices
 * \param ifac     local index of the face located on a parallel boundary
 *
 * \return iploc index (0,1 or 2) the node with highest coordinates defining triangle ifac
 *
 * Find the node index (0,1 or 2) with highest coordinates defining triangle ifac
 * and sort the vertices by increasing order.
 *
 */
int PMMG_sort_vertices(MMG5_pMesh mesh,MMG5_int k,MMG5_int *v_t,int ifac) {
  MMG5_pTetra pt;
  int         iploc;

  /* Pointer to the tetra structure */
  pt = &mesh->tetra[k];

  /* Local node indices of the 3 vertices on ifac */
  v_t[0] = pt->v[MMG5_idir[ifac][0]];
  v_t[1] = pt->v[MMG5_idir[ifac][1]];
  v_t[2] = pt->v[MMG5_idir[ifac][2]];

  /* Index [0,1,2] of the node with highest coordinates defining triangle ifac */
  iploc = PMMG_tria_highestcoord(mesh,v_t);

  /* Global node indices of the 3 vertices on ifac */
  v_t[0] = mesh->point[pt->v[MMG5_idir[ifac][0]]].tmp;
  v_t[1] = mesh->point[pt->v[MMG5_idir[ifac][1]]].tmp;
  v_t[2] = mesh->point[pt->v[MMG5_idir[ifac][2]]].tmp;

  /* Sort the vertices by increasing order of the global node indices */
  PMMG_swap_vertices(v_t);

  return iploc;

}

/**
 * \param tetra Indices of the tetras to be sorted
 * \param node  Indices of the nodes associated with the tetras
 * \param v_t0  First  tetra: indices of triangle vertices
 * \param v_t1  Second tetra: indices of triangle vertices
 * \param v_t2  Third  tetra: indices of triangle vertices
 *
 * Sort the tetras in increasing order based on vertices indices
 * Sort accordingly to the tetra sorting, the array storing the indices of nodes
 *
 */
void PMMG_sort_tetra(MMG5_int *tetra,MMG5_int *node,MMG5_int *v_t0,MMG5_int *v_t1,MMG5_int *v_t2) {
  /* Sorting using conditional statements */
  if ( v_t1[0] != -1 ) {
    if (PMMG_compare_3ints_array(v_t0, v_t1) > 0) {
      PMMG_swap_ints(&tetra[0], &tetra[1]);
      PMMG_swap_ints(&node[0], &node[1]);
      PMMG_swap_3int_arrays(v_t0, v_t1);
    }
    if ( v_t2[0] != -1 ) {
      if (PMMG_compare_3ints_array(v_t1, v_t2) > 0) {
        PMMG_swap_ints(&tetra[1], &tetra[2]);
        PMMG_swap_ints(&node[1], &node[2]);
        PMMG_swap_3int_arrays(v_t1, v_t2);
      }
      if (PMMG_compare_3ints_array(v_t0, v_t1) > 0) {
        PMMG_swap_ints(&tetra[0], &tetra[1]);
        PMMG_swap_ints(&node[0], &node[1]);
        PMMG_swap_3int_arrays(v_t0, v_t1);
      }
    }
  }
}

/**
 * \param v_t table of a triangle vertices
 *
 * Sort the vertices of the triangle in increasing order
 *
 */
void PMMG_swap_vertices(MMG5_int *v_t) {
  /* Sorting using conditional statements */
  if (v_t[0] > v_t[1]) {
    PMMG_swap_ints(&v_t[0], &v_t[1]);
  }
  if (v_t[1] > v_t[2]) {
    PMMG_swap_ints(&v_t[1], &v_t[2]);
  }
  if (v_t[0] > v_t[1]) {
    PMMG_swap_ints(&v_t[0], &v_t[1]);
  }
}

/**
 * \param a table of the 3 vertices of first  triangle
 * \param b table of the 3 vertices of second triangle
 *
 * \return -1 if a < b
 * \return  0 if a = b
 * \return +1 if a > b
 *
 * Compare vertices of 2 triangles to sort the triangle by increasing order
 * of their vertices indices
 *
 */
int PMMG_compare_3ints_array(int *a, int *b) {
  MMG5_int result;
  if (a[0] > b[0]) return 1;
  else if (a[0] < b[0]) return -1;

  if (a[1] > b[1]) return 1;
  else if (a[1] < b[1]) return -1;

  if (a[2] > b[2]) return 1;
  else if (a[2] < b[2]) return -1;

  return 0;
}

/**
 * \param a first  integer to swap
 * \param b second integer to swap
 *
 * Swap the integer a and b
 *
 */
void PMMG_swap_ints(int *a, int *b) {
  MMG5_int temp = *a;
  *a = *b;
  *b = temp;
}

/**
 * \param a first  array of 3 integers to swap
 * \param b second array of 3 integers to swap
 *
 * Swap the array of 3 integers a and b
 *
 */
void PMMG_swap_3int_arrays(int *a, int *b) {
  for ( int i = 0; i < 3; i++ ) {
      MMG5_int temp = a[i];
      a[i] = b[i];
      b[i] = temp;
  }
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param mesh pointer toward the mesh
 * \param sol pointer toward the level-set
 *
 * \return 1 if success, 0 otherwise
 *
 * \todo Fill the funtion
 *
 * Removal of small parasitic components (bubbles of material, etc) with volume
 * less than mesh->info.rmc (default VOLFRAC) * volume of the mesh.
 *
 */
int PMMG_rmc(PMMG_pParMesh parmesh,MMG5_pMesh mesh,MMG5_pSol sol){
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n      ## TODO:: PMMG_rmc.\n");
  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set function.
 * \param comm MPI communicator for ParMmg
 *
 * \return 1 if success, 0 if fail.
 *
 * \todo Fill the funtion
 *
 * Snap values of the level set function very close to 0 to exactly 0,
 * and prevent nonmanifold patterns from being generated.
 *
 * \todo all MPI_Abort have to be removed and replaced by a clean error handling
 * without deadlocks.
 *
 */
int PMMG_snpval_ls(PMMG_pParMesh parmesh,MPI_Comm comm) {

  MMG5_pTetra   pt;
  MMG5_pPoint   p0;
  PMMG_pInt_comm int_comm;
  PMMG_pExt_comm ext_comm;
  PMMG_pGrp      grp;
  MPI_Status     status;
  MMG5_pMesh     mesh;
  MMG5_pSol      sol;

  int nitem_ext,next_comm,nitem_ToShare_ToSend,nitem_ToShare_ToRecv;
  int color_in,color_out;
  int idx_ext,idx_int;
  int icomm,i,ip,idx;
  double *tmp;
  MMG5_int k,nc,ns,ncg;
  double         *rtosend,*rtorecv,*doublevalues;
  int            *itosend,*itorecv,*intvalues;
  int ier = 1;       // Initialize error

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n      ## TODO:: PMMG_snpval_ls.\n");

  /* Ensure only one group on each proc */
  assert(parmesh->ngrp == 1 && "more than one group per rank not implemented");

  /* Get external node communicator information */
  grp  = &parmesh->listgrp[0];
  mesh =  parmesh->listgrp[0].mesh;
  sol  = parmesh->listgrp[0].ls;
  next_comm = parmesh->next_node_comm; // Nbr of external node communicators
  int_comm  = parmesh->int_node_comm;  // Internal node communicator

  /** STEP 1 - Create tetra adjacency */
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Error: %s: hashing problem (1). Exit program.\n",
      __func__);
    ier = 0;
  }

  /* Reset point flags and s */
  for (k=1; k<=mesh->np; k++) {
    mesh->point[k].flag =  0;
    mesh->point[k].s    = -1;
  }

  /* Allocation memory */
  PMMG_CALLOC(parmesh,int_comm->intvalues,int_comm->nitem,int,"intvalues",ier = 0);
  MMG5_ADD_MEM(mesh,(mesh->npmax+1)*sizeof(double),"temporary table",
                fprintf(stderr,"  Exit program.\n"); ier = 0);

  if ( ier ) {
    MMG5_SAFE_CALLOC(tmp,mesh->npmax+1,double,ier = 0);
  }

  if ( !ier ) {
    /* Comms of step 6 will fail */
    MPI_Abort(parmesh->comm, PMMG_TMPFAILURE);
  }

  /** STEP 2 - Identify proc owner of interface points */
  /* Store point index in internal communicator intvalues */
  for( i = 0; i < grp->nitem_int_node_comm; i++ ){
    ip   = grp->node2int_node_comm_index1[i];
    idx  = grp->node2int_node_comm_index2[i];
    int_comm->intvalues[idx] = ip;
  }

  /* The highest rank to which this point belong store in mesh->point.s */
  for (icomm=0; icomm<next_comm; icomm++) {
    ext_comm  = &parmesh->ext_node_comm[icomm]; // External node communicator
    color_in  = ext_comm->color_in;             // Color of this partition Pcolor_in
    color_out = ext_comm->color_out;            // Color of the remote partition Pcolor_out
    nitem_ext = ext_comm->nitem;                // Nbr of nodes in common between Pcolor_in and Pcolor_out

    /* Loop over the nodes in the external node communicator **/
    for (i=0; i < nitem_ext; i++) {
      /* Get the indices of the nodes in internal communicators */
      idx = ext_comm->int_comm_index[i];
      ip  = int_comm->intvalues[idx];
      p0  = &mesh->point[ip];

      if (color_out > p0->s)
        p0->s = color_out;
      if (color_in > p0->s)
        p0->s = color_in;
    }
  }

  /** STEP 3 - Include tetras with very poor quality that are connected to the negative part */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !pt->v[0] ) continue;
    if ( pt->tag & MG_OVERLAP) continue; // Ignore overlap tetra
    if ( pt->qual < MMG5_EPS ) {
      fprintf(stdout, " PROC %d - Bad qual=%f < 1e-6 at tetra=%d \n", parmesh->myrank,pt->qual,k);
      for (i=0; i<4; i++) {
        ip = pt->v[i];
        if ( sol->m[ip] < 1000.0*MMG5_EPS ) break;
      }
      if ( i < 4 ) {
        for (i=0; i<4; i++) {
          ip = pt->v[i];
          sol->m[ip] = -1000.0*MMG5_EPS;
        }
      }
    }
  }

  /** STEP 4 - Snap values of sol that are close to 0 to 0 exactly */
  ns = 0;
  for (k=1; k<=mesh->np; k++) {
    p0 = &mesh->point[k];
    if ( !MG_VOK(p0) ) continue;
    if (p0->tag & MG_OVERLAP) continue; // Ignore overlap points
    /* Snap points in the interior of the partition and
       interface points with proc owner being this proc color */
    if ( (p0->s == -1) || (p0->s == parmesh->myrank))  {
      if ( fabs(sol->m[k]) < MMG5_EPS ) {
        if ( mesh->info.ddebug )
          fprintf(stderr,"  ## Warning: %s: snapping value at vertex %" MMG5_PRId "; "
                  "previous value: %E.\n",__func__,k,fabs(sol->m[k]));

        tmp[k] = ( fabs(sol->m[k]) < MMG5_EPSD ) ?
          (-100.0*MMG5_EPS) : sol->m[k];

        if ( parmesh->ddebug ) {
          fprintf(stderr, "  ## Warning: %s: rank %d - snapping value at "
                  "vertex %d, s=%d, tmp=%f, sol=%f \n",__func__,
                  parmesh->myrank,k,p0->s,tmp[k],sol->m[k]);
        }

        p0->flag = 1;
        sol->m[k] = 0;
        ns++;
      }
    }
  }

  /** STEP 5 - Check snapping did not lead to a nonmanifold situation */
  ncg = 0;
  do {
    nc = 0;
    /* Check snapping did not lead to a nonmanifold situation */
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) ) continue;
      if (pt->tag & MG_OVERLAP) continue; // Ignore overlap tetra
      for (i=0; i<4; i++) {
        ip = pt->v[i];
        p0 = &mesh->point[ip];
        // if (p0->tag & MG_PARBDY) continue;
        if ( p0->flag == 1 ) {
          if ( parmesh->ddebug ) {
            fprintf(stdout, "  ## Info: %s: rank %d - call MMG3D_ismaniball:\n"
                            "                              Tetra=%d, Point=%d, maniball=%d \n",
                    __func__,parmesh->myrank,k,ip,MMG3D_ismaniball(mesh,sol,k,i));
          }
          if ( !MMG3D_ismaniball(mesh,sol,k,i) ) {
            if ( tmp[ip] < 0.0 )
              sol->m[ip] = -100.0*MMG5_EPS;
            else
              sol->m[ip] = +100.0*MMG5_EPS;

            p0->flag = -1;
            nc++;
          }
        }
      }
    }
    ncg += nc;
  }
  while ( nc );

  /** TODO :: STEP 6 - Transfer data of snap point to the other proc */
  // nitem_ToShare_ToSend = 0; // Nbr of nodes in common between Pcolor_in and Pcolor_out
  // PMMG_CALLOC(parmesh,intvalues,    mesh->np, int,    "intvalues",    ier = 0);
  // PMMG_CALLOC(parmesh,doublevalues, mesh->np, double, "doublevalues", ier = 0);

  // if ( !ier ) {
  //  /* Avoid deadlock in comms */
  //  MPI_Abort(parmesh->comm, PMMG_TMPFAILURE);
  // }

  // for (k=1; k<=mesh->np; k++) {
  //   p0 = &mesh->point[k];
  //   if ( !MG_VOK(p0) ) continue;
  //   if ( !(p0->tag & MG_PARBDY)  ) continue; // If the point is not MG_PARBDY,  ignore it
  //   if ( !(p0->tag & MG_OVERLAP) ) continue; // If the point is not MG_OVERLAP, ignore it

  //   /* If this point is MG_PARBDY or MG_OVERLAP and has been modified on this partition (p0->flag==-1) */
  //   if ( (p0->flag == -1) )  {
  //     intvalues[nitem_ToShare_ToSend] = ip;
  //     doublevalues[nitem_ToShare_ToSend] = sol->m[ip];
  //     nitem_ToShare_ToSend +=1;

  //   }
  // }

  // for (icomm=0; icomm<next_comm; icomm++) {
  //   ext_comm  = &parmesh->ext_node_comm[icomm]; // External node communicator
  //   color_in  = ext_comm->color_in;             // Color of this partition Pcolor_in
  //   color_out = ext_comm->color_out;            // Color of the remote partition Pcolor_out
  //   nitem_ext = ext_comm->nitem;                // Nbr of nodes in common between Pcolor_in and Pcolor_out
  //   nitem_ToShare_ToSend = 0; // Nbr of nodes in common between Pcolor_in and Pcolor_out

  //   itosend   = ext_comm->itosend;
  //   rtosend   = ext_comm->rtosend;
  //   itorecv   = ext_comm->itorecv;
  //   rtorecv   = ext_comm->rtorecv;

  //   PMMG_CALLOC(parmesh,itosend, nitem_ext, int,    "itosend", ier = 0);
  //   PMMG_CALLOC(parmesh,rtosend, nitem_ext, double, "rtosend", ier = 0);

  //   if ( !ier ) {
  //     /* Avoid deadlock in comms */
  //     MPI_Abort(parmesh->comm, PMMG_TMPFAILURE);
  //   }

  //   /* Loop over the nodes in the external node communicator **/
  //   for (i=0; i < nitem_ext; i++) {
  //     /* Get the indices of the nodes in internal communicators */
  //     idx_ext = ext_comm->int_comm_index[i];
  //     idx_int = grp->node2int_node_comm_index2[idx_ext];
  //     ip      = grp->node2int_node_comm_index1[idx_int];
  //     p0      = &mesh->point[ip];

  //     /* If this point has been treated on this partition, send the data to the other */
  //     if ( (p0->s == color_in) & (p0->flag == -1))  {
  //       itosend[nitem_ToShare_ToSend] = ip;
  //       rtosend[nitem_ToShare_ToSend] = sol->m[ip];
  //       nitem_ToShare_ToSend +=1;

  //     }
  //   }

  //   /* Communication */
  //   // PMMG_CALLOC(parmesh,ext_comm->nitem_to_share,1,int,"nitem_to_share",ier = 0);
  //   MPI_CHECK(
  //     MPI_Sendrecv(&nitem_ToShare_ToSend,1,MPI_INT,color_out,MPI_LS_TAG+2,
  //                  &nitem_ToShare_ToRecv,1,MPI_INT,color_out,MPI_LS_TAG+2,
  //                  comm,&status),MPI_Abort(parmesh->comm, PMMG_TMPFAILURE) );

  //   ext_comm->nitem_to_share = nitem_ToShare_ToRecv;

  //   PMMG_CALLOC(parmesh,itorecv, nitem_ToShare_ToRecv, int,    "itorecv", ier = 0);
  //   PMMG_CALLOC(parmesh,rtorecv, nitem_ToShare_ToRecv, double, "rtorecv", ier = 0);

  //   if ( !ier ) {
  //     /* Avoid deadlock in comms */
  //     MPI_Abort(parmesh->comm, PMMG_TMPFAILURE);
  //   }

  //   MPI_CHECK(
  //     MPI_Sendrecv(itosend,nitem_ToShare_ToSend,MPI_INT,color_out,MPI_LS_TAG,
  //                  itorecv,nitem_ToShare_ToRecv,MPI_INT,color_out,MPI_LS_TAG,
  //                  comm,&status), MPI_Abort(parmesh->comm, PMMG_TMPFAILURE) );
  //   MPI_CHECK(
  //     MPI_Sendrecv(rtosend,nitem_ToShare_ToSend,MPI_DOUBLE,color_out,MPI_LS_TAG+1,
  //                  rtorecv,nitem_ToShare_ToRecv,MPI_DOUBLE,color_out,MPI_LS_TAG+1,
  //                  comm,&status), MPI_Abort(parmesh->comm, PMMG_TMPFAILURE) );

  // }

  // if ( (abs(mesh->info.imprim) > 5 || mesh->info.ddebug) && ns+ncg > 0 )
    fprintf(stdout,"  PROC %d -   %8" MMG5_PRId " points snapped, %" MMG5_PRId " corrected\n",parmesh->myrank,ns,ncg);

  /* Reset point flags */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /* memory free */
  MMG5_DEL_MEM(mesh,mesh->adja);
  MMG5_DEL_MEM(mesh,tmp);

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set.
 * \param met pointer toward  a metric (optionnal).
 *
 * \return 0 if fail, 1 otherwise.
 *
 * Create implicit surface in mesh.
 *
 * \todo all MPI_Abort have to be removed and replaced by a clean error handling
 * without deadlocks.
 *
 */
int PMMG_ls(PMMG_pParMesh parmesh) {
  char str[16]="";
  MMG5_HGeom hpar;
  MMG5_pMesh mesh;
  MMG5_pSol  met,sol;
  MMG5_int k;
  int      ier = 1;

  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;
  sol  = parmesh->listgrp[0].ls;

  /* Set function pointers */
  /** \todo TODO :: Surface ls and alias functions */
  if ( mesh->info.isosurf ) {
    fprintf(stderr," ## Error: Splitting boundaries on isovalue not yet"
            " implemented. Exit program.\n");
    ier = 0;
  }

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"  ** ISOSURFACE EXTRACTION %s\n",str);

  if ( mesh->nprism || mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: Isosurface extraction not available with"
            " hybrid meshes. Exit program.\n");
    ier = 0;
  }

  /* Modify the value of the level-set to work with the 0 level-set  */
  for (k=1; k<= sol->np; k++)
    sol->m[k] -= mesh->info.ls;

  /* Create overlap */
  if ( !ier ) {
    MPI_Abort(parmesh->comm,PMMG_TMPFAILURE);
  }

  if ( !PMMG_create_overlap(parmesh,parmesh->info.read_comm) ) {
    /* To avoid deadlocks in snpval_ls */
    MPI_Abort(parmesh->comm,PMMG_TMPFAILURE);
  }

  /** \todo TODO :: Snap values of level set function if needed */
  if ( !PMMG_snpval_ls(parmesh,parmesh->info.read_comm) ) {
    fprintf(stderr,"\n  ## Problem with implicit function. Exit program.\n");
    /* To avoid deadlocks in parbdyTria */
    ier = 0;
  }

  /* Delete overlap */
  if ( !PMMG_delete_overlap(parmesh,parmesh->info.read_comm) ) {
    fprintf(stderr,"\n  ## Impossible to delete overlap. Exit program.\n");
    ier = 0;
  }

  /* Create table of adjacency for tetra */
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    ier = 0;
  }

  /* Reset the mesh->info.isoref field everywhere */
  if ( ier && !MMG3D_resetRef_ls(mesh) ) {
    fprintf(stderr,"\n  ## Problem in resetting references. Exit program.\n");
    ier = 0;
  }

  /* Tag parallel triangles on material interfaces as boundary. */
  if ( !ier ) {
    /* Avoid deadlock in comms in parbdyTria */
    MPI_Abort(parmesh->comm,PMMG_TMPFAILURE);
  }

  if( !PMMG_parbdyTria( parmesh ) ) {
    fprintf(stderr,"\n  ## Unable to recognize parallel triangles on material interfaces."
            " Exit program.\n");
    ier =  0;
  }

  /* Check the compatibility of triangle orientation with tetra faces */
  if ( !MMG5_bdryPerm(mesh) ) {
    fprintf(stderr,"\n  ## Boundary orientation problem. Exit program.\n");
    ier = 0;
  }

  /* Identify surface mesh. Clean triangle array: remove useless or double
     triangles and add the missing ones.  Remark: spurious boundary triangles
     across parallel interface cannot be removed by the serial function but will
     not be stored inthe xtetra by the MMG5_bdrySet function during analysis.
     This may create inconsistencies between edge and point tags.
  */
  if ( ier && !MMG5_chkBdryTria(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    ier = 0;
  }

  /* Build hash table for initial edges: gather tag infos from edges and
   * triangles and store these infos in tria. Skip non PARBDYBDY // edges. */
  if ( ier && !MMG5_hGeom(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem (0). Exit program.\n");
    ier = 0;
  }

  /* Set the triangles references to the tetrahedra faces and edges */
  if ( ier && !MMG5_bdrySet(mesh) ) {
    fprintf(stderr,"\n  ## Problem in setting boundary. Exit program.\n");
    ier = 0;
  }

  /** \todo TODO :: Removal of small parasitic components */
  if ( mesh->info.rmc > 0 && ier ) {
    PMMG_rmc(parmesh,mesh,sol);
    fprintf(stdout,"\n  ## Warning: rmc option not implemented yet for ParMmg.\n");
    ier = 0;
  }

#ifdef USE_POINTMAP
  /* OK - Initialize source point with input index */
  MMG5_int ip;
  for( ip = 1; ip <= mesh->np; ip++ ) {
      if ( (!MG_VOK(&mesh->point[ip])) ) continue;
      mesh->point[ip].src = ip;
  }
#endif

  /* Compute vertices global numerotation
     This step is needed to compute the edge communicator */

  if ( !ier ) {
    /* Avoid deadlock in comms in compute_verticesGloNum */
    MPI_Abort(parmesh->comm,PMMG_TMPFAILURE);
  }

  if ( !PMMG_Compute_verticesGloNum( parmesh,parmesh->info.read_comm ) ) {
    fprintf(stderr,"\n  ## Warning: impossible to compute node global numbering.\n");
    ier = 0;
  }

  /* Hash parallel edges
     This step is needed to compute the edge communicator */
  if( ier && (PMMG_hashPar_fromFaceComm( parmesh,&hpar ) != PMMG_SUCCESS) ) {
    fprintf(stderr,"\n  ## Warning: impossible to compute the hash parallel edge.\n");
    ier = 0;
  }

  if ( !ier ) {
    /* Avoid deadlock in comms in build_edgeComm */
    MPI_Abort(parmesh->comm,PMMG_TMPFAILURE);
  }

  /* Build edge communicator */
  if( !PMMG_build_edgeComm( parmesh,mesh,&hpar,parmesh->info.read_comm ) ) {
    fprintf(stderr,"\n  ## Warning: Impossible to build edge communicator.\n");
    MPI_Abort(parmesh->comm,PMMG_TMPFAILURE);
  }

  assert ( PMMG_check_extEdgeComm ( parmesh,parmesh->info.read_comm ) );

  /** Discretization of the implicit function - Cut tetra */
  if ( !PMMG_cuttet_ls(parmesh) ) {
    fprintf(stderr,"\n  ## Problem in discretizing implicit function. Exit program.\n");
    ier = 0;
  }

  /* Delete outdated arrays */
  MMG5_DEL_MEM(mesh,mesh->adja);
  MMG5_DEL_MEM(mesh,mesh->adjt);
  MMG5_DEL_MEM(mesh,mesh->tria);

  /* The function MMG5_hGeom is in charge to set mesh->na=0 if not already.
     Here mesh->na is modified by the creation of the edge comm PMMG_build_edgeComm and is
     set equal to the number of // edges.
     Unfortunately here PMMG_build_edgeComm cannot be called before MMG5_hGeom.
     Hence, mesh->na is then set equal to 0 here because:
        1. later in the analysis, mesh->na needs to be equal to 0 and;
        2. at this stage, the edge comm does not exist anymore, so mesh->na should be 0 */
  mesh->na = 0;
  mesh->nt = 0;

  /* Update mesh->npi and mesh->nei to be equal to mesh->np and mesh->ne, respectively */
  mesh->npi = mesh->np;
  mesh->nei = mesh->ne;

  /* Set ref to tetra according to the sign of the level-set */
  if ( ier && !MMG3D_setref_ls(mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem in setting references. Exit program.\n");
    ier = 0;
  }

  /* Clean old bdy analysis */
  for ( k=1; k<=mesh->np; ++k ) {
    if ( mesh->point[k].tag & MG_BDY ) {
      mesh->point[k].tag &= ~MG_BDY;
    }
    if ( mesh->point[k].tag & MG_PARBDYBDY ) {
      mesh->point[k].tag &= ~MG_PARBDYBDY;
    }
  }

  /* Clean memory */
  MMG5_DEL_MEM(mesh,sol->m);

  if ( !ier ) {
    /* Avoid deadlock in comms in build_edgeComm */
    MPI_Abort(parmesh->comm,PMMG_TMPFAILURE);
  }

  /* Check communicators */
  assert ( PMMG_check_extFaceComm ( parmesh,parmesh->info.read_comm ) );
  assert ( PMMG_check_intFaceComm ( parmesh ) );
  assert ( PMMG_check_extNodeComm ( parmesh,parmesh->info.read_comm ) );
  assert ( PMMG_check_intNodeComm ( parmesh ) );

  /* Dealloc edge comm  as it is not up-to-date */
  MMG5_DEL_MEM(mesh,hpar.geom);
  PMMG_edge_comm_free( parmesh );

  return 1;
}

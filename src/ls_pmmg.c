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
 * \file isovalue_pmmg.c
 * \brief Create implicit surface in distribuited mesh.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (InriaSoft)
 * \author Laetitia Mottet (UBordeaux)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Main library functions (parallel remeshing starting from centralized or
 * distributed data.
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
 * \todo Fill the funtion
 *
 * Proceed to discretization of the implicit function carried by sol into mesh,
 * once values of sol have been snapped/checked
 *
 */
int PMMG_cuttet_ls(PMMG_pParMesh parmesh, MMG5_pMesh mesh, MMG5_pSol sol, MMG5_pSol met){
  MMG5_pTetra  pt,pt0;
  MMG5_pxTetra pxt,pxt0;
  MMG5_pPoint  p0,p1;
  MMG5_Hash    hash;

  PMMG_pExt_comm ext_node_comm,ext_edge_comm,ext_face_comm;
  PMMG_pGrp      grp;

  MMG5_int ne_init,ne_tmp;
  MMG5_int k,k0;
  MMG5_int ip0,ip1,np,nb,ns,src,refext,refint;
  MMG5_int vGlobNum[4],vx[6];
  MMG5_int tetra_sorted[3], node_sorted[3];
  MMG5_int *ne_tmp_tab,*vGlobNum_tab;

  static int8_t  mmgWarn = 0;
  int8_t         ia;
  int8_t         i,j;
  int8_t         npneg,nface_added;

  const uint8_t *taued=NULL;
  uint8_t        tau[4];
  uint8_t        imin0,imin2;

  double c[3],v0,v1,s;

  int already_split;
  int idx_tmp;
  int i_commn,i_comme,i_commf;
  int iedge,iface;
  int ifac,iploc;
  int flag;
  int ier;

  int nitem_int_node,nitem_int_face;
  int nitem_ext_node,nitem_ext_edge,nitem_ext_face;
  int nitem_ext_face_init;
  int next_node_comm,next_face_comm,next_edge_comm;

  int color_in_node,color_out_node;
  int color_in_edge,color_out_edge;

  int idx_edge_ext,idx_edge_int,idx_edge_mesh;
  int idx_face_ext,idx_face_int,val_face;

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n      ## PMMG_cuttet_ls:: Under development.\n");

  /* Ensure only one group on each proc */
  assert(parmesh->ngrp == 1);

  /* Initialization */
  grp = &parmesh->listgrp[0];
  next_node_comm = parmesh->next_node_comm;  // Number of communicator for nodes
  next_edge_comm = parmesh->next_edge_comm;  // Number of communicator for edges
  next_face_comm = parmesh->next_face_comm;  // Number of communicator for faces
  nitem_int_node = grp->nitem_int_node_comm; // Number of initial total nodes in internal node communicator
  nitem_int_face = grp->nitem_int_face_comm; // Number of initial total faces in internal node communicator

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
               nitem_int_node+2*nb,
               nitem_int_node,
               int,"Allocation of node2int_node_comm_index1", return 0);
  PMMG_REALLOC(parmesh, grp->node2int_node_comm_index2,
               nitem_int_node+2*nb,
               nitem_int_node,
               int,"Allocation of node2int_node_comm_index2", return 0);

  /* STEP 3.3 - Realloc internal face communicators */
  PMMG_REALLOC(parmesh, grp->face2int_face_comm_index1,
               3*nitem_int_face,
               nitem_int_face,
               int,"Allocation of face2int_face_comm_index1", return 0);
  PMMG_REALLOC(parmesh, grp->face2int_face_comm_index2,
               3*nitem_int_face,
               nitem_int_face,
               int,"Allocation of face2int_face_comm_index2", return 0);

  /* STEP 3.4 - Realloc external node communicator */
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    PMMG_REALLOC(parmesh,ext_node_comm->int_comm_index,
                 ext_node_comm->nitem+2*nb,
                 ext_node_comm->nitem,
                 int,"Allocation of external node communicator",return 0);
  }

  /* STEP 3.5 - Realloc external face communicator */
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    PMMG_REALLOC(parmesh,ext_face_comm->int_comm_index,
                 ext_face_comm->nitem+2*ext_face_comm->nitem,
                 ext_face_comm->nitem,
                 int,"Allocation of external face communicator",return 0);
  }

  /* STEP 3.6 - Allocate all the other variables needed to be allocated ! */
  PMMG_CALLOC( parmesh,vGlobNum_tab,4*(nitem_int_face),MMG5_int,"vGlobNum_tab",return 0 );
  PMMG_CALLOC( parmesh,ne_tmp_tab,nitem_int_face+1,MMG5_int,"ne_tmp_tab",return 0 );

  /** STEP 4 - Identify required edges. Put hash.item[key].k = -1 */
  /* Loop over tetra */
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

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
  /* Loop over the number of edge communicator */
  for (i_comme=0; i_comme < next_edge_comm; i_comme++) {

    /* Get external edge communicator information */
    ext_edge_comm  = &parmesh->ext_edge_comm[i_comme]; // External edge communicator
    color_in_edge  = ext_edge_comm->color_in;          // Color of the hosting proc - this proc
    color_out_edge = ext_edge_comm->color_out;         // Color of the remote  proc - the proc to exchange with
    nitem_ext_edge = ext_edge_comm->nitem;             // Number of edges in common between these 2 procs

    /* Loop over the edges in the external edge communicator */
    for (iedge=0; iedge < nitem_ext_edge; iedge++) {

      /* Get the indices of the edge in internal communicators and mesh->edge */
      idx_edge_ext  = ext_edge_comm->int_comm_index[iedge];
      idx_edge_int  = grp->edge2int_edge_comm_index2[idx_edge_ext];
      idx_edge_mesh = grp->edge2int_edge_comm_index1[idx_edge_int];

      /* Find extremities of this edge */
      ip0 = mesh->edge[idx_edge_mesh].a;
      ip1 = mesh->edge[idx_edge_mesh].b;

      /* Get np (i.e. hash.item[key].k) and ns (i.e. hash.item[key].s) associated with this edge.
         np: index of the node along this edge in mesh->point. If np>0, this node already exists.
         ns: index of the node along this edge in internal node comm index2. */
      PMMG_hashGet_all(&hash,ip0,ip1,&np,&ns);

      /* STEP 5.1.1 - Update external node communicators if edge already split */
      /*      If np>0 (i.e. hash.item[key].k != [0;-1]), this edge has already been split;
              and new node on this edge is already in the internal node comm.
              Update external node comm if needed, then pass to the next edge. */
      if ( np>0 ) {
        /* If mesh->point[np].s != (i_comme+1), add this node to the external node communicator.
           Otherwise, this node has already been added to the external node comm - do nothing.  */
        if ( mesh->point[np].s != (i_comme+1) ) {
          /* (a) Find the appropriate external node comm */
          for (i_commn=0; i_commn < next_node_comm; i_commn++) {
            ext_node_comm  = &parmesh->ext_node_comm[i_commn]; // External node communicator
            color_in_node  = ext_node_comm->color_in;          // Color of the hosting proc - this proc
            color_out_node = ext_node_comm->color_out;         // Color of the remote  proc - the proc to exchange with
            assert(color_in_node == color_in_edge);            // Ensure that the hosting proc is the same
            /* (b) If color_out of edge and node comm are the same - Update external node communicator */
            if (color_out_node == color_out_edge) {
              nitem_ext_node = ext_node_comm->nitem;              // Initial nbr of nodes in common between these 2 procs
              ext_node_comm->int_comm_index[nitem_ext_node] = ns; // Add the node to the external node comm
              ext_node_comm->nitem = nitem_ext_node + 1;          // Updated nbr of nodes in common between these 2 procs
              break;
            }
          }
          /* (c) Update mesh->point[np].s to specify that this node already exists in the external node comm */
          mesh->point[np].s = i_comme+1;
        }
        continue;
      }

      // TODO:: Multimaterial - Check whether an entity with reference ref should be split

      /* STEP 5.1.2 - Create a new point if this edge needs to be split */
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

      /* STEP 5.1.3 - Update of met, sol and field for the new point */
      // TODO:: Add field interpolation for the new point
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

      /* STEP 5.1.4 - Update node communicators and edge hash table */
      /* (a) Update the internal node communicators */
      grp->node2int_node_comm_index1[nitem_int_node] = np;
      grp->node2int_node_comm_index2[nitem_int_node] = nitem_int_node;

      /* (b) Update the external node communicator */
      /* Find the appropriate external node comm */
      for (i_commn=0; i_commn < next_node_comm; i_commn++) {
        ext_node_comm  = &parmesh->ext_node_comm[i_commn]; // External node communicator
        color_in_node  = ext_node_comm->color_in;          // Color of the hosting proc - this proc
        color_out_node = ext_node_comm->color_out;         // Color of the remote  proc - the proc to exchange with
        assert(color_in_node == color_in_edge);            // Ensure that the hosting proc is the same

        /* If color_out of edge and node comm are the same - Update external node communicator */
        if (color_out_node == color_out_edge) {
          nitem_ext_node = ext_node_comm->nitem;                          // Initial nbr of nodes in common between these 2 procs
          ext_node_comm->int_comm_index[nitem_ext_node] = nitem_int_node; // Add the node to the external node comm
          ext_node_comm->nitem = nitem_ext_node + 1;                      // Updated nbr of nodes in common between these 2 procs
          break;
        }
      }

      /* (c) Add the index of this edge comm into mesh->point[np].s to specify */
      /*     this node is already in this external node communicator */
      mesh->point[np].s = i_comme+1;

      /* (d) Update hash hash.item[key].k = -1 becomes = np and         */
      /*                 hash.item[key].s = -1 becomes = nitem_int_node */
      PMMG_hashUpdate_all(&hash,ip0,ip1,np,nitem_int_node);

      /* (e) Update the total number of nodes in internal node communicator */
      nitem_int_node += 1;
      parmesh->int_node_comm->nitem = nitem_int_node;
      grp->nitem_int_node_comm      = nitem_int_node;
    }
  }

  /** STEP 5.2 - Create all the other new points located elsewhere */
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

      /* If np>0 (i.e. hash.item[key].k != [0;-1]), this edge has already been split, pass to the next edge. */
      if ( np>0 ) continue;

      /* Check whether an entity with reference ref should be split */
      if ( !MMG5_isSplit(mesh,pt->ref,&refint,&refext) ) continue;

      /* STEP 5.2.1 - Create a new point if this edge needs to be split */
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

      /* STEP 5.2.2 - Update of met, sol and field for the new point */
      // TODO:: Add field interpolation for the new point
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

  /** STEP 6 - Split according to tets flags */
  /** STEP 6.1 - Compute global node vertices */
  if ( !PMMG_Compute_verticesGloNum( parmesh,parmesh->comm ) ) {
    if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
      fprintf(stdout,"\n\n\n  -- WARNING: IMPOSSIBLE TO COMPUTE NODE GLOBAL NUMBERING\n\n\n");
      PMMG_RETURN_AND_FREE( parmesh, PMMG_LOWFAILURE );
    }
  }

  /** STEP 6.2 - Do the splitting for tetra on parallel interface */
  ns  = 0;     // Number of total split on this proc
  ier = 1;     // Error
  idx_tmp = 0; // Index of an already split tetra to recover info stored in ne_tmp_tab and vGlobNum_tab

  /* Loop over the number of faces communicator */
  for (i_commf=0; i_commf < next_face_comm; i_commf++) {

    /* Get current external face communicator */
    ext_face_comm  = &parmesh->ext_face_comm[i_commf]; // External face communicator
    nitem_ext_face = ext_face_comm->nitem;             // Number of faces in common between these 2 procs
    nitem_ext_face_init = ext_face_comm->nitem;        // Initial number of faces in common between these 2 procs

    /* Loop over the faces in the external face communicator */
    for (iface=0; iface < nitem_ext_face_init; iface++) {

      /* Get the index of the face in internal communicator and the value of the face */
      idx_face_ext = ext_face_comm->int_comm_index[iface];
      idx_face_int = grp->face2int_face_comm_index2[idx_face_ext];
      val_face     = grp->face2int_face_comm_index1[idx_face_int];

      /* Find the local tetra, the face and node associated */
      k     = val_face/12;     // Index of the tetra on this proc
      ifac  = (val_face%12)/3; // Index of the face
      iploc = (val_face%12)%3; // Index of the node

      /* Get the tetra k and xtetra associated */
      pt  = &mesh->tetra[k];
      pxt = &mesh->xtetra[pt->xt];
      if ( !MG_EOK(pt) )  continue;

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
      memset(tetra_sorted,-1,3*sizeof(MMG5_int));
      memset(node_sorted, -1,3*sizeof(MMG5_int));

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
        imin0=MMG3D_split2sf_cfg(flag,vGlobNum,tau,&taued); // Compute tau and imin0
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
        break;
      }

      if ( !ier ) return 0;

      /* STEP 6.2.3 - Update tag of edges in xtetra with MG_PARBDY */
      for (j=0; j<3; j++) {
        k0 = tetra_sorted[j];
        if (k0 != -1) {
          pt0  = &mesh->tetra[k0];
          pxt0 = &mesh->xtetra[pt0->xt];
          for (i=0; i<3; i++) {
            ia = MMG5_iarf[ifac][i];
            if ( !(pxt0->tag[ia] && MG_PARBDY) ) {
              pxt0->tag[ia] |= MG_PARBDY;
            }
          }
        }
      }

      /* STEP 6.2.4 - Update face communicators */
      nitem_ext_face = ext_face_comm->nitem; // Number of faces in common between these 2 procs

      /* (a) Update the first face located at idx_face_int - Modify only index1 - index2 stays the same */
      grp->face2int_face_comm_index1[idx_face_int] = 12*tetra_sorted[0]+3*ifac+node_sorted[0];

      /* (b) Update the communicators for the potential 2 other faces */
      nface_added = 0;
      for (j=0; j<2; j++) {
        if ( tetra_sorted[j+1] != -1) {
          grp->face2int_face_comm_index1[nitem_int_face+j] = 12*tetra_sorted[j+1]+3*ifac+node_sorted[j+1];
          grp->face2int_face_comm_index2[nitem_int_face+j] = nitem_int_face+j;
          ext_face_comm->int_comm_index[nitem_ext_face+j]  = nitem_int_face+j;
          nface_added += 1;
        }
      }

      /* (c) Update the total number of faces */
      nitem_int_face += nface_added;
      nitem_ext_face += nface_added;
      parmesh->int_face_comm->nitem = nitem_int_face;
      grp->nitem_int_face_comm      = nitem_int_face;
      ext_face_comm->nitem          = nitem_ext_face;
    }
  }

  /** STEP 6.3 - Do the splitting for tetra located elsewhere */
  /* Loop over tetra */
  for (k=1; k<=ne_init; k++) {

    /* Get the tetra k and xtetra associated */
    pt  = &mesh->tetra[k];
    pxt = &mesh->xtetra[pt->xt];
    if ( !MG_EOK(pt) )  continue;

    /* If the tetra has already a flag (flag !=0) - it has already been processed. Pass to the next tetra. */
    if (pt->flag) continue;

    /* STEP 6.3.1 - Find global numbering and split pattern of the tetra */
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

    /* STEP 6.3.2 - If not already done, split the tetra according to the flag */
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

  /** STEP 7 - Deallocation of memory and reset of some fields */
  /* Dealloc edge comm  as it is not up-to-date */
  PMMG_edge_comm_free( parmesh );

  /* Delete the tables storing imin0, imin2 and ne_tmp_tab */
  PMMG_DEL_MEM(parmesh,vGlobNum_tab,MMG5_int,"vGlobNum_tab");
  PMMG_DEL_MEM(parmesh,ne_tmp_tab,MMG5_int,"ne_tmp_tab");

  /* Delete the edges hash table */
  MMG5_DEL_MEM(mesh,hash.item);

  /* Reset mark and flag in mesh->tetra */
  for (k=1; k<=ne_init; k++) {
    mesh->tetra[k].mark = 0;
    mesh->tetra[k].flag = 0;
  }

  /* Reset s in mesh->point */
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].s = 0;

  return ns;
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

  MMG5_pTetra pt;
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
  if ( (ifac==tau[0]) ) tetra_sorted[0] = ne_tmp;

  /* Global node indices of the 3 vertices on ifac */
  pt     = &mesh->tetra[tetra_sorted[0]];
  v_t0[0] = mesh->point[pt->v[MMG5_idir[ifac][0]]].tmp;
  v_t0[1] = mesh->point[pt->v[MMG5_idir[ifac][1]]].tmp;
  v_t0[2] = mesh->point[pt->v[MMG5_idir[ifac][2]]].tmp;

  /* Index of the minimum global node index defining ifac */
  node_sorted[0]=PMMG_find_node(v_t0);

  /* Sort the vertices by increasing order */
  PMMG_sort_vertices(v_t0);

  /* STEP 1.2 - Index of the second tetra */
  /* Tetra #1 created by MMG5_split1 */
  tetra_sorted[1] = ne_tmp;
  /* Except for the following: no more tetra define ifac - we assign -1 */
  if ( ifac == tau[0] ) tetra_sorted[1] = -1;
  if ( ifac == tau[1] ) tetra_sorted[1] = -1;

  if ( tetra_sorted[1] != -1 ) {
    /* Global node indices of the 3 vertices on ifac */
    pt     = &mesh->tetra[tetra_sorted[1]];
    v_t1[0] = mesh->point[pt->v[MMG5_idir[ifac][0]]].tmp;
    v_t1[1] = mesh->point[pt->v[MMG5_idir[ifac][1]]].tmp;
    v_t1[2] = mesh->point[pt->v[MMG5_idir[ifac][2]]].tmp;

    /* Index of the minimum global node index defining ifac */
    node_sorted[1]=PMMG_find_node(v_t1);

    /* Sort the vertices by increasing order */
    PMMG_sort_vertices(v_t1);
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

  MMG5_pTetra pt;
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

  /* Global node indices of the 3 vertices on ifac */
  pt     = &mesh->tetra[tetra_sorted[0]];
  v_t0[0] = mesh->point[pt->v[MMG5_idir[ifac][0]]].tmp;
  v_t0[1] = mesh->point[pt->v[MMG5_idir[ifac][1]]].tmp;
  v_t0[2] = mesh->point[pt->v[MMG5_idir[ifac][2]]].tmp;

  /* Index of the minimum global node index defining ifac */
  node_sorted[0]=PMMG_find_node(v_t0);

  /* Sort the vertices by increasing order */
  PMMG_sort_vertices(v_t0);

  /* STEP 1.2 - Index of the second tetra */
  /* Tetra #1 created by MMG5_split2sf_globNum */
  tetra_sorted[1] = ne_tmp-1;
  /* Except for the following: treta #2 created by MMG5_split2sf_globNum or no more tetra */
  if ( ifac == tau[1] ) tetra_sorted[1] = ne_tmp;
  if ( ifac == tau[3] ) tetra_sorted[1] = -1;

  if ( tetra_sorted[1] != -1 ) {
    /* Global node indices of the 3 vertices on ifac */
    pt     = &mesh->tetra[tetra_sorted[1]];
    v_t1[0] = mesh->point[pt->v[MMG5_idir[ifac][0]]].tmp;
    v_t1[1] = mesh->point[pt->v[MMG5_idir[ifac][1]]].tmp;
    v_t1[2] = mesh->point[pt->v[MMG5_idir[ifac][2]]].tmp;

    /* Index of the minimum global node index defining ifac */
    node_sorted[1]=PMMG_find_node(v_t1);

    /* Sort the vertices by increasing order */
    PMMG_sort_vertices(v_t1);
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
    /* Global node indices of the 3 vertices on ifac */
    pt     = &mesh->tetra[tetra_sorted[2]];
    v_t2[0] = mesh->point[pt->v[MMG5_idir[ifac][0]]].tmp;
    v_t2[1] = mesh->point[pt->v[MMG5_idir[ifac][1]]].tmp;
    v_t2[2] = mesh->point[pt->v[MMG5_idir[ifac][2]]].tmp;

    /* Index of the minimum global node index defining ifac */
    node_sorted[2]=PMMG_find_node(v_t2);

    /* Sort the vertices by increasing order */
    PMMG_sort_vertices(v_t2);
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

  MMG5_pTetra pt;
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

  /* Global node indices of the 3 vertices on ifac */
  pt     = &mesh->tetra[tetra_sorted[0]];
  v_t0[0] = mesh->point[pt->v[MMG5_idir[ifac][0]]].tmp;
  v_t0[1] = mesh->point[pt->v[MMG5_idir[ifac][1]]].tmp;
  v_t0[2] = mesh->point[pt->v[MMG5_idir[ifac][2]]].tmp;

  /* Index of the minimum global node index defining ifac */
  node_sorted[0]=PMMG_find_node(v_t0);

  /* Sort the vertices by increasing order */
  PMMG_sort_vertices(v_t0);

  /* STEP 1.2 - Index of the second tetra */
  /* Tetra #1 created by MMG5_split3cone_globNum */
  tetra_sorted[1] = ne_tmp-2;
  /* Except for the following: treta #2 created by MMG5_split3cone_globNum or no more tetra */
  if ( (ia==tau[1]) && (ifac == tau[1]) ) tetra_sorted[1] = ne_tmp-1;
  if ( (ia==tau[2]) && (ifac == tau[2]) ) tetra_sorted[1] = ne_tmp-1;
  if ( (ia==tau[3]) && (ifac == tau[3]) ) tetra_sorted[1] = ne_tmp-1;
  if ( ifac == tau[0] ) tetra_sorted[1] = -1;

  if ( tetra_sorted[1] != -1 ) {
    /* Global node indices of the 3 vertices on ifac */
    pt      = &mesh->tetra[tetra_sorted[1]];
    v_t1[0] = mesh->point[pt->v[MMG5_idir[ifac][0]]].tmp;
    v_t1[1] = mesh->point[pt->v[MMG5_idir[ifac][1]]].tmp;
    v_t1[2] = mesh->point[pt->v[MMG5_idir[ifac][2]]].tmp;

    /* Index of the minimum global node index defining ifac */
    node_sorted[1]=PMMG_find_node(v_t1);

    /* Sort the vertices by increasing order */
    PMMG_sort_vertices(v_t1);
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
    /* Global node indices of the 3 vertices on ifac */
    pt     = &mesh->tetra[tetra_sorted[2]];
    v_t2[0] = mesh->point[pt->v[MMG5_idir[ifac][0]]].tmp;
    v_t2[1] = mesh->point[pt->v[MMG5_idir[ifac][1]]].tmp;
    v_t2[2] = mesh->point[pt->v[MMG5_idir[ifac][2]]].tmp;

    /* Index of the minimum global node index defining ifac */
    node_sorted[2]=PMMG_find_node(v_t2);

    /* Sort the vertices by increasing order */
    PMMG_sort_vertices(v_t2);
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

  MMG5_pTetra pt;
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

  /* Global node indices of the 3 vertices on ifac */
  pt     = &mesh->tetra[tetra_sorted[0]];
  v_t0[0] = mesh->point[pt->v[MMG5_idir[ifac][0]]].tmp;
  v_t0[1] = mesh->point[pt->v[MMG5_idir[ifac][1]]].tmp;
  v_t0[2] = mesh->point[pt->v[MMG5_idir[ifac][2]]].tmp;

  /* Index of the minimum index of the nodes defining ifac */
  node_sorted[0]=PMMG_find_node(v_t0);

  /* Sort the vertices by increasing order */
  PMMG_sort_vertices(v_t0);

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

  /* Global indices of points defining ifac */
  pt     = &mesh->tetra[tetra_sorted[1]];
  v_t1[0] = mesh->point[pt->v[MMG5_idir[ifac][0]]].tmp;
  v_t1[1] = mesh->point[pt->v[MMG5_idir[ifac][1]]].tmp;
  v_t1[2] = mesh->point[pt->v[MMG5_idir[ifac][2]]].tmp;

  /* Index of the minimum index of the nodes defining ifac */
  node_sorted[1]=PMMG_find_node(v_t1);

  /* Sort the vertices by increasing order */
  PMMG_sort_vertices(v_t1);

  /* STEP 1.3 - Index of the third tetra */
  /* Tetra #5 created by MMG5_split4op_globNum */
  tetra_sorted[2] = ne_tmp;
  /* Except for the following: treta #3 or #4 created by MMG5_split4op_globNum */
  if ( (imin23==tau[2]) && (ifac == tau[0]) ) tetra_sorted[2] = ne_tmp-1;
  if ( (imin23==tau[3]) && (ifac == tau[1]) ) tetra_sorted[2] = ne_tmp-1;
  if ( (imin23==tau[2]) && (ifac == tau[2]) ) tetra_sorted[2] = ne_tmp-2;
  if ( (imin23==tau[3]) && (ifac == tau[3]) ) tetra_sorted[2] = ne_tmp-2;

  /* Global node indices of the 3 vertices on ifac */
  pt     = &mesh->tetra[tetra_sorted[2]];
  v_t2[0] = mesh->point[pt->v[MMG5_idir[ifac][0]]].tmp;
  v_t2[1] = mesh->point[pt->v[MMG5_idir[ifac][1]]].tmp;
  v_t2[2] = mesh->point[pt->v[MMG5_idir[ifac][2]]].tmp;

  /* Index of the minimum index of the nodes defining ifac */
  node_sorted[2]=PMMG_find_node(v_t2);

  /* Sort the vertices by increasing order */
  PMMG_sort_vertices(v_t2);

  /* STEP 2 - Sort these tetras by their global indices */
  PMMG_sort_tetra(tetra_sorted,node_sorted,v_t0,v_t1,v_t2);

  return;
}

/**
 * \param v_t Indices of a triangle vertices
 *
 * Find the local index of the minimum of the 3 indices defining the triangle
 *
 * \return 0, 1 or 2 the local index of the minimum
 *
 */
int PMMG_find_node(MMG5_int *v_t) {
  MMG5_int min = MG_MIN(v_t[0],MG_MIN(v_t[1],v_t[2]));
  if (min==v_t[0]) return 0;
  else if (min==v_t[1]) return 1;
  else return 2;
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
void PMMG_sort_tetra(MMG5_int *tetra, MMG5_int *node, MMG5_int *v_t0, MMG5_int *v_t1, MMG5_int *v_t2) {
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
void PMMG_sort_vertices(MMG5_int *v_t) {
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
 *
 * \return 1 if success, 0 if fail.
 *
 * \todo Fill the funtion
 *
 * Snap values of the level set function very close to 0 to exactly 0,
 * and prevent nonmanifold patterns from being generated.
 *
 */
int PMMG_snpval_ls(PMMG_pParMesh parmesh,MMG5_pMesh mesh,MMG5_pSol sol) {
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n      ## TODO:: PMMG_snpval_ls.\n");
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
 */
int PMMG_ls(PMMG_pParMesh parmesh, MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pSol met) {
  char str[16]="";
  MMG5_HGeom hpar;

  /* Set function pointers */
  /** \todo TODO :: Surface ls and alias functions */
  if ( mesh->info.isosurf ) {
    fprintf(stderr," ## Error: Splitting boundaries on isovalue not yet"
            " implemented. Exit program.\n");
    return 0;
  }

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"  ** ISOSURFACE EXTRACTION %s\n",str);

  if ( mesh->nprism || mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: Isosurface extraction not available with"
            " hybrid meshes. Exit program.\n");
    return 0;
  }

  /* Modify the value of the level-set to work with the 0 level-set  */
  MMG5_int k;
  for (k=1; k<= sol->np; k++)
    sol->m[k] -= mesh->info.ls;

  /** \todo TODO :: Snap values of level set function if needed */
  if ( !PMMG_snpval_ls(parmesh,mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem with implicit function. Exit program.\n");
    return 0;
  }

  /* OK - Create table of adjacency for tetra */
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  /* OK - Check the compatibility of triangle orientation with tetra faces */
  if ( !MMG5_bdryPerm(mesh) ) {
    fprintf(stderr,"\n  ## Boundary orientation problem. Exit program.\n");
    return 0;
  }

  /* TO BE CHECKED :: Check behaviour with PMMG_APIDISTRIB_nodes
    Identify surface mesh
    Clean triangle array - remove useless or double triangles
    and add the missing ones */
  if ( !MMG5_chkBdryTria(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    return 0;
  }

  /* OK - Build hash table for initial edges */
  if ( !MMG5_hGeom(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem (0). Exit program.\n");
    return 0;
  }

  /* OK - Set the triangles references to the tetrahedra faces and edges */
  if ( !MMG5_bdrySet(mesh) ) {
    fprintf(stderr,"\n  ## Problem in setting boundary. Exit program.\n");
    return 0;
  }

  /* OK - Reset the mesh->info.isoref field everywhere */
  if ( !MMG3D_resetRef_ls(mesh) ) {
    fprintf(stderr,"\n  ## Problem in resetting references. Exit program.\n");
    return 0;
  }

  /** \todo TODO :: Removal of small parasitic components */
  if ( mesh->info.rmc > 0 ) {
    PMMG_rmc(parmesh,mesh,sol);
    fprintf(stdout,"\n  ## Warning: rmc option not implemented yet for ParMmg\n");
    return 0;
  }

#ifdef USE_POINTMAP
  /* OK - Initialize source point with input index */
  MMG5_int ip;
  for( ip = 1; ip <= mesh->np; ip++ )
    mesh->point[ip].src = ip;
#endif

  /* OK - Compute vertices and triangles global numerotation */
  if ( !PMMG_Compute_verticesGloNum( parmesh,parmesh->comm ) ) {
    if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
      fprintf(stdout,"\n\n\n  -- WARNING: IMPOSSIBLE TO COMPUTE NODE GLOBAL NUMBERING\n\n\n");
      PMMG_RETURN_AND_FREE( parmesh, PMMG_LOWFAILURE );
    }
  }

  /* Hash parallel edges */
  if( PMMG_hashPar_pmmg( parmesh,&hpar ) != PMMG_SUCCESS ) {
    if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
      fprintf(stdout,"\n\n\n  -- WARNING: Impossible to compute the hash parallel edge \n\n\n");
      PMMG_RETURN_AND_FREE( parmesh, PMMG_LOWFAILURE );
    }
  }

  /* Build edge communicator */
  if( !PMMG_build_edgeComm( parmesh,mesh,&hpar,parmesh->comm ) ) {
    if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
      fprintf(stdout,"\n\n\n  -- WARNING: Impossible to build edge communicator \n\n\n");
      PMMG_RETURN_AND_FREE( parmesh, PMMG_LOWFAILURE );
    }
  }

  /** \todo TODO :: Discretization of the implicit funtion - Cut tetra */
  if ( !PMMG_cuttet_ls(parmesh,mesh,sol,met) ) {
    fprintf(stderr,"\n  ## Problem in discretizing implicit function. Exit program.\n");
    return 0;
  }

  /* Not sure which function to be used to deallocate memory */
  MMG5_DEL_MEM(mesh,mesh->adja);
  MMG5_DEL_MEM(mesh,mesh->adjt);
  MMG5_DEL_MEM(mesh,mesh->tria);

  mesh->nt = 0;

  /* OK - Set ref to tetra according to the sign of the level-set */
  if ( !MMG3D_setref_ls(mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem in setting references. Exit program.\n");
    return 0;
  }

  /* Clean old bdy analysis */
  for ( MMG5_int k=1; k<=mesh->np; ++k ) {
    if ( mesh->point[k].tag & MG_BDY ) {
      mesh->point[k].tag &= ~MG_BDY;
    }
  }

  /* Clean memory */
  MMG5_DEL_MEM(mesh,sol->m);

  return 1;
}

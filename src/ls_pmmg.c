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
int PMMG_split_ls(PMMG_pParMesh parmesh, MMG5_pMesh mesh){

}



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
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n      ## PMMG_cuttet_ls:: Under development.\n");

  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   p0,p1;
  MMG5_Hash     hash;
  double        c[3],v0,v1,s;
  int           ier;
  int           bdy_tetra;
  MMG5_int      vx[6],vGlobNum[4],p_tmp;
  MMG5_int      k,ip0,ip1,np,nb,ns,src,refext,refint;
  int8_t        ia,j,npneg;
  static int8_t mmgWarn = 0;

  // Local variables added by LM
  int i_commn,i_comme,i_commf;
  int inode,iedge,iface;

  int nitem_int_node;
  int nitem_ext_node,nitem_ext_edge,nitem_ext_face;
  int next_node_comm,next_face_comm,next_edge_comm;

  int color_in_node,color_out_node;
  int color_in_edge,color_out_edge;
  int color_in_face,color_out_face;

  int pos_edge_ext,pos_edge_int,val_edge;
  int pos_face_ext,pos_face_int,val_face;

  int ip0_g, ip1_g;
  int ifac, iploc;
  int sign;

  PMMG_pExt_comm ext_face_comm, ext_node_comm,ext_edge_comm;
  PMMG_pGrp      grp;

  // Ensure only one group on each proc
  assert(parmesh->ngrp == 1);

  // Initialization
  // base = ++mesh->base;
  grp = &parmesh->listgrp[0];
  next_node_comm = parmesh->next_node_comm;  // Number of communicator for nodes
  next_face_comm = parmesh->next_face_comm;  // Number of communicator for edges
  next_edge_comm = parmesh->next_edge_comm;  // Number of communicator for faces
  nitem_int_node = grp->nitem_int_node_comm; // Number of initial total nodes in internal node communicator

  /*************************/
  /* STEP 1 :: Reset flags */
  /*************************/
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n          --> STEP 1 :: Reset flags ...");

  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  for (k=1; k<=mesh->ne; k++)
    mesh->tetra[k].flag = 0;

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n              ... STEP 1 :: Done. \n");

  /***********************************************************************/
  /* STEP 2 :: Approximate the number nb of intersection points on edges */
  /***********************************************************************/
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n          --> STEP 2 :: Approximate nbr of intersection points ...");

  nb = 0;

  // Loop over tetra
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];

    // Loop over edges
    for (ia=0; ia<6; ia++) {

      // Get the points defining the edges
      ip0 = pt->v[MMG5_iare[ia][0]];
      ip1 = pt->v[MMG5_iare[ia][1]];
      p0  = &mesh->point[ip0];
      p1  = &mesh->point[ip1];

      // If both points have flag, then pass as these points have been treated
      if ( p0->flag && p1->flag )  continue;

      // Otherwise take the values at these points
      v0  = sol->m[ip0];
      v1  = sol->m[ip1];

      // If the points are not already exactly on the level-set
      // and does not have the same sign, then this edge needs to be split
      if ( fabs(v0) > MMG5_EPSD2 && fabs(v1) > MMG5_EPSD2 && v0*v1 < 0.0 ) {
        // If the points have not been treated yet, assign a new flag and increase nb value
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

  if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
    fprintf(stdout,"\n              ... PROC %d :: Approximate nbr of intersection points nb ~ %d. \n",parmesh->myrank,nb);
    fprintf(stdout,"\n              ... STEP 2 :: Done. \n");
  }

  /* TODO :: test if the number of point proc by proc is correct */
  // Cannot be done here as it is an approximation. Otherwise, need to robustify step 2 above.
#ifndef NDEBUG
  /* TODO */
#endif

  /*******************************/
  /* STEP 3 :: Memory allocation */
  /*******************************/
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n          --> STEP 3 :: Memory allocation ...");

  // STEP 3.1 :: Initialize hash table for edges //
  if ( !MMG5_hashNew(mesh,&hash,nb,7*nb) ) return 0;

  // STEP 3.2 :: Realloc internal node communicator //
  PMMG_REALLOC(parmesh, grp->node2int_node_comm_index1,
               nitem_int_node+2*nb,
               nitem_int_node, int,"node2int_node_comm_index1",
               return 0);
  PMMG_REALLOC(parmesh, grp->node2int_node_comm_index2,
               nitem_int_node+2*nb,
               nitem_int_node, int,"node2int_node_comm_index2",
               return 0);

  // STEP 3.3 :: Realloc internal face communicator //
  // TODO

  // STEP 3.4 :: Realloc external node communicator //
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    PMMG_REALLOC(parmesh,ext_node_comm->int_comm_index,ext_node_comm->nitem+2*nb,
                  ext_node_comm->nitem,int,"external communicator",return 0);
  }

  // STEP 3.5 :: Realloc external face communicator //
  // TODO

  // STEP 3.6 :: Allocate all the other variables that need to be allocated !

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n              ... STEP 3 :: Done. \n");

  /****************************************************************/
  /* STEP 4 :: Identify required edges. Put hash.item[key].k = -1 */
  /****************************************************************/
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n          --> STEP 4 :: Identify required edges and put hash.item[key].k = -1 ...");

  // Loop over tetra
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    //-----------------------------------------------------------------//
    // Step 4.1 :: Identification of edges belonging to a required tet //
    //-----------------------------------------------------------------//
    // If the tetra is required MG_REQ
    if ( pt->tag & MG_REQ ) {
      // Loop over the edges
      for (ia=0; ia<6; ia++) {
        ip0 = pt->v[MMG5_iare[ia][0]];
        ip1 = pt->v[MMG5_iare[ia][1]];
        np  = -1;
        // Add an edge to the edge table with hash.item[key].k = -1
        if ( !MMG5_hashEdge(mesh,&hash,ip0,ip1,np) )  return -1;
      }
      continue;
    }

    //------------------------------------------------------------------------------------------------//
    // Step 4.2 :: Identification of edges belonging to a (par)boundary or being explicitely required //
    //------------------------------------------------------------------------------------------------//
    // If the xtetra associated to this tetra exists
    if ( !pt->xt ) continue;

    // Point towards the xtetra corresponding to the tetra...
    pxt = &mesh->xtetra[pt->xt];
    // ... then loop over the faces
    for (ia=0; ia<4; ia++) {

      // (a) If the face is not a boundary MG_BDY, then continue
      if ( !(pxt->ftag[ia] & MG_BDY) ) continue;

      // (a) otherwise loop over the edges
      for (j=0; j<3; j++) {

        // (b) If the edges is not required, then continue
        if ( !(pxt->tag[ MMG5_iarf[ia][j] ] & MG_REQ) ) continue;

        // (b) otherwise get the extremity of the edges ...
        ip0 = pt->v[MMG5_idir[ia][MMG5_inxt2[j]]];
        ip1 = pt->v[MMG5_idir[ia][MMG5_iprv2[j]]];
        np  = -1;

        // (c) ... and add an edge to the edge table with hash.item[key].k = -1
        if ( !MMG5_hashEdge(mesh,&hash,ip0,ip1,np) )  return -1;
      }
    }
  }

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n              ... STEP 4 :: Done. \n");

  /************************************************************************/
  /* STEP 5 :: Create points at iso-value. Fill or update edge hash table */
  /************************************************************************/
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n          --> STEP 5 :: Create points at iso-value. Fill or update edge hash table...");

  //---------------------------------------------------------------//
  // STEP 5.1 :: Create new points located on parallel interfaces  //
  //---------------------------------------------------------------//
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n                 STEP 5.1 :: Create points located on parallel interfaces...");

  // TO BE ADDED - If Multimat; allocation; metric; MG_EOK
  // NOTE :: Here we will split only the edges of the faces on parallel interfaces

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n                  Number of node comm :: %d ; edge comm ::: %d & face comm :: %d",
                   next_node_comm,next_edge_comm,next_face_comm);

  // Loop over the number of edge communicator
  for (i_comme=0; i_comme < next_edge_comm; i_comme++) {

    // Get current external edge communicator
    ext_edge_comm  = &parmesh->ext_edge_comm[i_comme]; // Current external edge communicator
    color_in_edge  = ext_edge_comm->color_in;          // Color of the input proc - this proc
    color_out_edge = ext_edge_comm->color_out;         // Color of the ouput proc - the proc to exchange with
    nitem_ext_edge = ext_edge_comm->nitem;             // Number of edges in common between these 2 procs

    // Loop over the edges in the external communicator
    for (iedge=0; iedge < nitem_ext_edge; iedge++) {

      // Get the index of the edge in internal comm
      pos_edge_ext = ext_edge_comm->int_comm_index[iedge];
      pos_edge_int = grp->edge2int_edge_comm_index2[pos_edge_ext];
      val_edge     = grp->edge2int_edge_comm_index1[pos_edge_int];

      // Find extremities of this edge
      ip0 = mesh->edge[val_edge].a;
      ip1 = mesh->edge[val_edge].b;

      MMG5_hashGet_all(&hash,ip0,ip1,&np,&ns);

      ip0_g = mesh->point[ip0].tmp;
      ip1_g = mesh->point[ip1].tmp;

      // If hash.item[key].k is not 0 or -1 then pass, it means this edge has already been split
      if ( np>0 ) {
        // 1. Update external communicator
        // 1.1. Find the appropriate external node comm
        if ( mesh->point[np].s != (i_comme+1) ) {
          if ( parmesh->info.imprim > PMMG_VERB_VERSION )
          fprintf(stdout,"\n                    MyRank %d :: Proc %d -> %d - Edge tetra %d, "
                        "ip0 %d, ip1 %d, ip0_g %d, ip1_g %d :: Already split at pos %d in int comm - Add into this comm \n",
                        parmesh->myrank, color_in_edge, color_out_edge, iedge,
                        ip0, ip1, ip0_g, ip1_g,ns);
          for (i_commn=0; i_commn < next_node_comm; i_commn++) {
            ext_node_comm  = &parmesh->ext_node_comm[i_commn]; // Current external node communicator
            color_out_node = ext_node_comm->color_out;         // Proc to exchange with
            if (color_out_node == color_out_edge) {
              // 1.2. Update the number of nodes
              nitem_ext_node = ext_node_comm->nitem; // Number of nodes in common between these 2 procs
              ext_node_comm->int_comm_index[nitem_ext_node] = ns;
              ext_node_comm->nitem = nitem_ext_node + 1;
              break;
            }
          }
          mesh->point[np].s = i_comme+1;
        }
        else {
          if ( parmesh->info.imprim > PMMG_VERB_VERSION )
          fprintf(stdout,"\n                    MyRank %d :: Proc %d -> %d - Edge tetra %d, "
                        "ip0 %d, ip1 %d, ip0_g %d, ip1_g %d :: Already split at pos %d in int comm - Already in this comm \n",
                        parmesh->myrank, color_in_edge, color_out_edge, iedge,
                        ip0, ip1, ip0_g, ip1_g,ns);
        }
        continue;
      }

      // Check the ls value at the edge nodes
      p0 = &mesh->point[ip0];
      p1 = &mesh->point[ip1];
      v0 = sol->m[ip0];
      v1 = sol->m[ip1];

      // Check if the edge should be split
      // If one of the points is exactly on the level set, the point exists already, pass
      if ( fabs(v0) < MMG5_EPSD2 || fabs(v1) < MMG5_EPSD2 )
        continue;
      // If the points have the same sign, no need to split the edge, pass
      else if ( MG_SMSGN(v0,v1) )
        continue;
      // If one or the other point has never been treated, pass
      else if ( !p0->flag || !p1->flag )
        continue;

      // Define the weighting factor
      s = v0 / (v0-v1);
      s = MG_MAX(MG_MIN(s,1.0-MMG5_EPS),MMG5_EPS);

      // Find the coordinates of the new points using the weighting factor
      c[0] = p0->c[0] + s*(p1->c[0]-p0->c[0]);
      c[1] = p0->c[1] + s*(p1->c[1]-p0->c[1]);
      c[2] = p0->c[2] + s*(p1->c[2]-p0->c[2]);

      // Create a new point with coordinates c, tag 0 and source src
      // Return the new number of points in the partition np
#ifdef USE_POINTMAP
      src = p0->src;
#else
      src = 1;
#endif
      np = MMG3D_newPt(mesh,c,0,src);

      if ( parmesh->info.imprim > PMMG_VERB_VERSION )
        fprintf(stdout,"\n                    MyRank %d :: Proc %d -> %d - "
                      "Edge tetra %d, ip0 %d, ip1 %d, ip0_g %d, ip1_g %d :: Split edge np=%d \n",
                      parmesh->myrank, color_in_edge, color_out_edge, iedge,
                      ip0, ip1, ip0_g, ip1_g,np);

      // Update node communicator
      // 1. Update the internal node comm
      grp->node2int_node_comm_index1[nitem_int_node] = np;
      grp->node2int_node_comm_index2[nitem_int_node] = nitem_int_node;

      // 2. Update external communicator
      // 2.1. Find the appropriate external node comm
      for (i_commn=0; i_commn < next_node_comm; i_commn++) {
        ext_node_comm  = &parmesh->ext_node_comm[i_commn]; // Current external node communicator
        color_out_node = ext_node_comm->color_out;         // Proc to exchange with
        if (color_out_node == color_out_edge) {
          // 1.2. Update the number of nodes
          nitem_ext_node = ext_node_comm->nitem; // Number of nodes in common between these 2 procs
          ext_node_comm->int_comm_index[nitem_ext_node] = nitem_int_node;
          ext_node_comm->nitem = nitem_ext_node + 1;
          break;
        }
      }

      mesh->point[np].s = i_comme+1;

      // For this new point, add the value of the solution, i.e. the isovalue 0
      sol->m[np] = 0;

      // Update the hash hash.item[key].k = -1 becomes = np
      MMG5_hashUpdate_all(&hash,ip0,ip1,np,nitem_int_node);

      // 3. Update the total number of nodes in internal node comm
      nitem_int_node += 1;
      parmesh->int_node_comm->nitem = nitem_int_node;
      grp->nitem_int_node_comm      = nitem_int_node;
    }
  }

  // Reset s in hash table
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].s = 0;

  //---------------------------------------------------------------//
  // STEP 5.3 :: Create all the other new points located elsewhere //
  //---------------------------------------------------------------//
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n                 STEP 5.2 :: Create points located elsewhere...\n");

  // TO BE ADDED - If Multimat; allocation; metric; warning; MG_EOK

  // Loop over tetra k
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];

    // Loop over the edges ia
    for (ia=0; ia<6; ia++) {
      // Get the points of the edge and np (value stored in hash.item[key].k)
      ip0 = pt->v[MMG5_iare[ia][0]];
      ip1 = pt->v[MMG5_iare[ia][1]];
      np  = MMG5_hashGet(&hash,ip0,ip1);

      // If hash.item[key].k is not 0 or -1 then pass, it means this edge has already been split
      if ( np>0 ) {
        if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
          fprintf(stdout,"                        Edge already split ia %d - ",ia);
          fprintf(stdout," Points %d & %d - ",ip0,ip1);
          fprintf(stdout," np %d \n",np);
        }
        continue;
      }

      // Check the ls value at the edge nodes
      p0 = &mesh->point[ip0];
      p1 = &mesh->point[ip1];
      v0 = sol->m[ip0];
      v1 = sol->m[ip1];

      // Check if the edge should be split
      // If one of the points is exactly on the level set, the point exists already, pass
      if ( fabs(v0) < MMG5_EPSD2 || fabs(v1) < MMG5_EPSD2 )
        continue;
      // If the points have the same sign, no need to split the edge, pass
      else if ( MG_SMSGN(v0,v1) )
        continue;
      // If one or the other point has never been treated, pass
      else if ( !p0->flag || !p1->flag )
        continue;

      // What is this ?
      npneg = (np<0);

      // Define the weighting factor
      s = v0 / (v0-v1);
      s = MG_MAX(MG_MIN(s,1.0-MMG5_EPS),MMG5_EPS);

      // Find the coordinates of the new points using the weighting factor
      c[0] = p0->c[0] + s*(p1->c[0]-p0->c[0]);
      c[1] = p0->c[1] + s*(p1->c[1]-p0->c[1]);
      c[2] = p0->c[2] + s*(p1->c[2]-p0->c[2]);

      // Create a new point with coordinates c, tag 0 and source src
      // Return the new number of points in the partition np
#ifdef USE_POINTMAP
      src = p0->src;
#else
      src = 1;
#endif
      np = MMG3D_newPt(mesh,c,0,src);

      if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
        fprintf(stdout,"                        Split edge ia %d - ",ia);
        fprintf(stdout," Points %d & %d - ",ip0,ip1);
        fprintf(stdout," Weighting factor %f - ",s);
        fprintf(stdout," New coord %f %f %f - ",c[0],c[1],c[2]);
        fprintf(stdout," np %d \n",np);
      }

      // For this new point, add the value of the solution, i.e. the isovalue 0
      sol->m[np] = 0;

      // If this edge is required, then inform the user we split it anyway
      // and update the hash hash.item[key].k = - 1 becomes = np
      // Otherwise add the edge to be split into hash table
      if ( npneg ) {
        /* We split a required edge */
        MMG5_hashUpdate(&hash,ip0,ip1,np);
      }
      else {
        MMG5_hashEdge(mesh,&hash,ip0,ip1,np);
      }
    }
  }

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n              ... STEP 5 :: Done. \n");

  /*******************************************/
  /* STEP 6 :: Split according to tets flags */
  /*******************************************/
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n          --> STEP 6 :: Split according to tets flags...\n");

  //------------------------------------------//
  // STEP 6.1 :: Compute global node vertices //
  //------------------------------------------//
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n                 STEP 6.1 :: Compute global node vertives...\n");

  if ( !PMMG_Compute_verticesGloNum( parmesh,parmesh->comm ) ) {
    if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
      fprintf(stdout,"\n\n\n  -- WARNING: IMPOSSIBLE TO COMPUTE NODE GLOBAL NUMBERING\n\n\n");
      PMMG_RETURN_AND_FREE( parmesh, PMMG_LOWFAILURE );
    }
  }

  //--------------------------------------------------------------//
  // STEP 6.2 :: Do the splitting for tetra on parallel interface //
  //--------------------------------------------------------------//
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n                 STEP 6.2 :: Split tetra on parallel interface...\n");

  ns  = 0; // Number of total split on this proc
  ier = 1;

  // Loop over the number of faces communicator
  for (i_commf=0; i_commf < next_face_comm; i_commf++) {

    // Get current external face communicator
    ext_face_comm  = &parmesh->ext_face_comm[i_commf]; // Current external face communicator
    color_in_face  = ext_face_comm->color_in;          // Color of the input proc - this proc
    color_out_face = ext_face_comm->color_out;         // Color of the ouput proc - the proc to exchange with
    nitem_ext_face = ext_face_comm->nitem;             // Number of faces in common between these 2 procs

    // Process edges in the same order on both MPI process
    // On proc with higher ID, value 1; otherwise value -1
    // sign = (color_in_face > color_out_face) ? -1 : 1 ;

    // Loop over the faces in the external communicator
    for (iface=0; iface < nitem_ext_face; iface++) {

      // Get the value of the face in internal comm
      pos_face_ext = ext_face_comm->int_comm_index[iface];
      pos_face_int = grp->face2int_face_comm_index2[pos_face_ext];
      val_face     = grp->face2int_face_comm_index1[pos_face_int];

      // Find the local tetra, the face and node associated
      k     = val_face/12;     // Index of the tetra on this proc
      ifac  = (val_face%12)/3; // Index of the face
      iploc = (val_face%12)%3; // Index of the node

      // Get the tetra k and xtetra associated
      pt  = &mesh->tetra[k];
      pxt = &mesh->xtetra[pt->xt];
      if ( !MG_EOK(pt) )  continue;

      // Flag the tetra to 0
      pt->flag = 0;
      memset(vx,0,6*sizeof(MMG5_int));

      // Loop over the edges, get hash.item[key].k
      // and set flags to the tetra that need to be split
      bdy_tetra = 0;
      for (ia=0; ia<6; ia++) {
        vx[ia] = MMG5_hashGet(&hash,pt->v[MMG5_iare[ia][0]],pt->v[MMG5_iare[ia][1]]);
        if ( vx[ia] > 0 )
          MG_SET(pt->flag,ia);
        // Get global num of the tetra
        if ( pxt->tag[ia] & MG_PARBDY )
          bdy_tetra = 1;
      }

      // Get global num of the tetra points
      for (j=0; j<4; j++) {
        p_tmp = pt->v[j];
        vGlobNum[j] = mesh->point[p_tmp].tmp;
      }

      // Do the actual split
      switch (pt->flag) {
      case 1: case 2: case 4: case 8: case 16: case 32: /* 1 edge split */
        if ( parmesh->info.imprim > PMMG_VERB_VERSION )
          fprintf(stdout,"\n                  ---------------------------------------");

        if ( parmesh->info.imprim > PMMG_VERB_VERSION )
          fprintf(stdout,"\n                  MyRank %d, Tetra k %d, MMG5_split1, Flag %d, vGlobNum %d-%d-%d-%d, vx %d-%d-%d-%d-%d-%d\n",
                  parmesh->myrank,k,pt->flag,vGlobNum[0],vGlobNum[1],vGlobNum[2],vGlobNum[3],vx[0],vx[1],vx[2],vx[3],vx[4],vx[5]);
          //-- TODO
        ier = MMG5_split1_GlobNum(mesh,met,k,vx,vGlobNum,1,parmesh->myrank);
        ns++;
        break;

      case 48: case 24: case 40: case 6: case 34: case 36:
      case 20: case 5: case 17: case 9: case 3: case 10: /* 2 edges (same face) split */
        if ( parmesh->info.imprim > PMMG_VERB_VERSION )
          fprintf(stdout,"\n                  ---------------------------------------");
        if ( parmesh->info.imprim > PMMG_VERB_VERSION )
          fprintf(stdout,"\n\n                  MyRank %d, Tetra k %d, MMG5_split2sf, Flag %d, vGlobNum %d-%d-%d-%d, vx %d-%d-%d-%d-%d-%d\n",
                  parmesh->myrank,k,pt->flag,vGlobNum[0],vGlobNum[1],vGlobNum[2],vGlobNum[3],vx[0],vx[1],vx[2],vx[3],vx[4],vx[5]);
        //-- TODO
        ier = MMG5_split2sf_GlobNum(mesh,met,k,vx,vGlobNum,1,parmesh->myrank);
        ns++;
        break;

      case 7: case 25: case 42: case 52: /* 3 edges on conic configuration split */
        if ( parmesh->info.imprim > PMMG_VERB_VERSION )
          fprintf(stdout,"\n                  ---------------------------------------");
        if ( parmesh->info.imprim > PMMG_VERB_VERSION )
          fprintf(stdout,"\n                  MyRank %d, Tetra k %d, MMG5_split3cone_GlobNum, Flag %d, vGlobNum %d-%d-%d-%d, vx %d-%d-%d-%d-%d-%d\n",
                  parmesh->myrank,k,pt->flag,vGlobNum[0],vGlobNum[1],vGlobNum[2],vGlobNum[3],vx[0],vx[1],vx[2],vx[3],vx[4],vx[5]);
        //-- TODO
        ier = MMG5_split3cone_GlobNum(mesh,met,k,vx,vGlobNum,1,parmesh->myrank);
        // ier = MMG5_split3cone(mesh,met,k,vx,1);
        ns++;
        break;

      case 30: case 45: case 51:
        if ( parmesh->info.imprim > PMMG_VERB_VERSION )
          fprintf(stdout,"\n                  ---------------------------------------");
        if ( parmesh->info.imprim > PMMG_VERB_VERSION )
          fprintf(stdout,"\n                  MyRank %d, Tetra k %d, MMG5_split4op_GlobNum, Flag %d, vGlobNum %d-%d-%d-%d, vx (%d)-(%d)-(%d)-(%d)-(%d)-(%d)\n",
                  parmesh->myrank,k,pt->flag,vGlobNum[0],vGlobNum[1],vGlobNum[2],vGlobNum[3],vx[0],vx[1],vx[2],vx[3],vx[4],vx[5]);
          //-- TODO
          ier = MMG5_split4op_GlobNum(mesh,met,k,vx,vGlobNum,1,parmesh->myrank);
        ns++;
        break;

      default :
        if ( parmesh->info.imprim > PMMG_VERB_VERSION )
          fprintf(stdout,"\n                  ---------------------------------------");
        assert(pt->flag == 0);
        if ( parmesh->info.imprim > PMMG_VERB_VERSION )
          fprintf(stdout,"\n                  MyRank %d, Tetra k %d, Flag %d, NO SPLIT",parmesh->myrank,k,pt->flag);
        break;
      }
      if ( !ier ) return 0;
    }

    // TODO - Update face communicator


  }


  //----------------------------------------------------------//
  // STEP 6.3 :: Do the splitting for tetra located elsewhere //
  //----------------------------------------------------------//
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n                 STEP 6.3 :: Split tetra located elsewhere...\n");

  // ACHTUNG / TODO :: ne - total nbr of tetra in mesh has changed like we have already
  // split some tetra, this has created nbew tetra.

  // Loop over tetra
  for (k=1; k<=mesh->ne; k++) {

    // Get the tetra k and xtetra associated
    pt  = &mesh->tetra[k];
    pxt = &mesh->xtetra[pt->xt];
    if ( !MG_EOK(pt) )  continue;

    // Flag the tetra to 0
    pt->flag = 0;
    memset(vx,0,6*sizeof(MMG5_int));

    // If this tetra has flag = 0 - then it is either
    // (i)  a tetra that has not been treated - a tetra in the volum
    // (ii) a tetra that has been treated but that does not need to be split - a tetra on parallel interface
    // For now if we are in case (ii), we estimate it is not costly to recompute vx - which is a priori wrong
    // If flag !=0, the tetra has already been split. We pass and do nothing.
    if (pt->flag) {
      if ( parmesh->info.imprim > PMMG_VERB_VERSION )
        fprintf(stdout,"                        Tetra already split k %d - ",k);
      continue;
    }
    // else {
    //   if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    //     fprintf(stdout,"                        Tetra to be split k %d - ",k);
    // }

    // Loop over the edges, get hash.item[key].k
    // and set flags to the tetra that need to be split
    bdy_tetra = 0;
    for (ia=0; ia<6; ia++) {
      vx[ia] = MMG5_hashGet(&hash,pt->v[MMG5_iare[ia][0]],pt->v[MMG5_iare[ia][1]]);
      if ( vx[ia] > 0 )
        MG_SET(pt->flag,ia);
      // Get global num of the tetra
      if ( pxt->tag[ia] & MG_PARBDY )
        bdy_tetra = 1;
    }

    // Get global num of the tetra points
    for (j=0; j<4; j++) {
      p_tmp = pt->v[j];
      vGlobNum[j] = mesh->point[p_tmp].tmp;
    }

    // // Do the actual split
    // switch (pt->flag) {
    // case 1: case 2: case 4: case 8: case 16: case 32: /* 1 edge split */
    //   if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    //       fprintf(stdout,"\n                  ---------------------------------------");
    //   if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    //     fprintf(stdout,"\n                  MyRank %d, Tetra k %d, MMG5_split1, Flag %d, vGlobNum %d-%d-%d-%d\n",
    //             parmesh->myrank,k,pt->flag,vGlobNum[0],vGlobNum[1],vGlobNum[2],vGlobNum[3]);
    //     //-- TODO
    //   ier = MMG5_split1_GlobNum(mesh,met,k,vx,vGlobNum,1,parmesh->myrank);
    //   ns++;
    //   break;

    // case 48: case 24: case 40: case 6: case 34: case 36:
    // case 20: case 5: case 17: case 9: case 3: case 10: /* 2 edges (same face) split */
    //   if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    //       fprintf(stdout,"\n                  ---------------------------------------");
    //   if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    //     fprintf(stdout,"\n                  MyRank %d, Tetra k %d, MMG5_split2sf, Flag %d, vGlobNum %d-%d-%d-%d\n",
    //             parmesh->myrank,k,pt->flag,vGlobNum[0],vGlobNum[1],vGlobNum[2],vGlobNum[3]);
    //   //-- TODO
    //   ier = MMG5_split2sf_GlobNum(mesh,met,k,vx,vGlobNum,1,parmesh->myrank);
    //   ns++;
    //   break;

    // case 7: case 25: case 42: case 52: /* 3 edges on conic configuration split */
    //   if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    //       fprintf(stdout,"\n                  ---------------------------------------");
    //   if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    //     fprintf(stdout,"\n                  MyRank %d, Tetra k %d, MMG5_split3cone_GlobNum, Flag %d, vGlobNum %d-%d-%d-%d\n",
    //             parmesh->myrank,k,pt->flag,vGlobNum[0],vGlobNum[1],vGlobNum[2],vGlobNum[3]);
    //   //-- TODO
    //   ier = MMG5_split3cone_GlobNum(mesh,met,k,vx,vGlobNum,1,parmesh->myrank);
    //   // ier = MMG5_split3cone(mesh,met,k,vx,1);
    //   ns++;
    //   break;

    // case 30: case 45: case 51:
    //   if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    //       fprintf(stdout,"\n                  ---------------------------------------");
    //   if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    //     fprintf(stdout,"\n                  MyRank %d, Tetra k %d, MMG5_split4op_GlobNum, Flag %d, vGlobNum %d-%d-%d-%d\n",
    //             parmesh->myrank,k,pt->flag,vGlobNum[0],vGlobNum[1],vGlobNum[2],vGlobNum[3]);
    //     //-- TODO
    //     ier = MMG5_split4op_GlobNum(mesh,met,k,vx,vGlobNum,1,parmesh->myrank);
    //   ns++;
    //   break;

    // default :
    //   assert(pt->flag == 0);
    //   if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    //       fprintf(stdout,"\n                  ---------------------------------------");
    //   if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    //     fprintf(stdout,"\n                  MyRank %d, Tetra k %d, Flag %d, NO SPLIT",parmesh->myrank,k,pt->flag);
    //   break;
    // }
    if ( !ier ) return 0;
  }

  // Print the number of split
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"     %7" MMG5_PRId " splitted\n",ns);
  // fprintf(stdout,"\n              ... PROC %d :: Nbr of split ns = %d. \n",parmesh->myrank,ns);

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n              ... STEP 6 :: Done. \n");

  // Delete the edges hash table
  MMG5_DEL_MEM(mesh,hash.item);

  // TODO : dealloc edge comm if not updated

  return ns;
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

  /* OK :: Create table of adjacency for tetra */
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  /* OK :: Check the compatibility of triangle orientation with tetra faces */
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

  /* OK :: Build hash table for initial edges */
  if ( !MMG5_hGeom(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem (0). Exit program.\n");
    return 0;
  }

  /* OK :: Set the triangles references to the tetrahedra faces and edges */
  if ( !MMG5_bdrySet(mesh) ) {
    fprintf(stderr,"\n  ## Problem in setting boundary. Exit program.\n");
    return 0;
  }

  /* OK :: Reset the mesh->info.isoref field everywhere */
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
  /* OK :: Initialize source point with input index */
  MMG5_int ip;
  for( ip = 1; ip <= mesh->np; ip++ )
    mesh->point[ip].src = ip;
#endif

  /* OK :: Compute vertices and triangles global numerotation */
  if ( !PMMG_Compute_verticesGloNum( parmesh,parmesh->comm ) ) {
    if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
      fprintf(stdout,"\n\n\n  -- WARNING: IMPOSSIBLE TO COMPUTE NODE GLOBAL NUMBERING\n\n\n");
      PMMG_RETURN_AND_FREE( parmesh, PMMG_LOWFAILURE );
    }
  }

  if ( !PMMG_Compute_trianglesGloNum( parmesh,parmesh->comm ) ) {
    if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
      fprintf(stdout,"\n\n\n  -- WARNING: IMPOSSIBLE TO COMPUTE TRIANGLE GLOBAL NUMBERING\n\n\n");
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

  /* OK :: Set ref to tetra according to the sign of the level-set */
  /* Comment for now as the level-set has not been performed yet it fails */
  // if ( !MMG3D_setref_ls(mesh,sol) ) {
  //   fprintf(stderr,"\n  ## Problem in setting references. Exit program.\n");
  //   return 0;
  // }

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

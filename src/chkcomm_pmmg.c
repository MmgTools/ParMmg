#include "parmmg.h"
#include "coorcell_pmmg.h"

/**
 * \param parmesh pointer toward a parmesh structure
 * \param min    minimal coordinates in each direction
 * \param max    maximal coordinates in each direction
 * \param delta  scaling value
 *
 * \return 1 if success, 0 if fail;
 *
 * Find the bounding box of the internal node communicator.
 *
 */
static
int PMMG_find_intNodeCommBoundingBox(PMMG_pParMesh parmesh,double min[3],
                                     double max[3],double *delta) {
  PMMG_coorCell *coor_list;
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  double         dd;
  int            ier,nitem,*intvalues,k,i,j,idx,ip;

  ier = 0;

  nitem = parmesh->int_node_comm->nitem;

  coor_list = NULL;
  PMMG_CALLOC(parmesh,coor_list,nitem,PMMG_coorCell,"node coordinates",
              return 0);
  PMMG_CALLOC(parmesh,parmesh->int_node_comm->intvalues,nitem,int,
              "intvalues array", goto end);
  intvalues = parmesh->int_node_comm->intvalues;

  /* Store the node coordinates of the nodes */
  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp  = &parmesh->listgrp[k];
    mesh = grp->mesh;
    assert ( mesh->info.delta &&  "missing scaling infos");

    dd = mesh->info.min[0]*mesh->info.min[0]
      + mesh->info.min[1]*mesh->info.min[1]
      + mesh->info.min[2]*mesh->info.min[2];

    assert ( fabs(mesh->info.delta-1.)<_MMG5_EPSD && dd<_MMG5_EPSD &&
             "scaled mesh... need to unscale it");

    for ( i=0; i<grp->nitem_int_node_comm; ++i ) {
      ip      = grp->node2int_node_comm_index1[i];
      idx     = grp->node2int_node_comm_index2[i];

      if ( !intvalues[idx] ) {
        for ( j=0; j<3; ++j )
          coor_list[idx].c[j] = mesh->point[ip].c[j];
      }
      else intvalues[idx] = 1;
    }
  }

  /* Find the bounding box of the internal comm */
  if ( !PMMG_find_coorCellListBoundingBox(coor_list,nitem,min,max,delta) )
    goto end;

  /* Success */
  ier = 1;

end:
  PMMG_DEL_MEM(parmesh,parmesh->int_node_comm->intvalues,nitem,int,
               "intvalues array");
  PMMG_DEL_MEM(parmesh,coor_list,nitem,PMMG_coorCell,"coor_list array");

  return ier;
}



/**
 * \param parmesh pointer toward a parmesh structure
 * \param min    minimal coordinates in each direction
 * \param max    maximal coordinates in each direction
 * \param delta  scaling value
 *
 * \return 1 if success, 0 if fail;
 *
 * Find the bounding box of the internal face communicator.
 *
 */
static
int PMMG_find_intFaceCommBoundingBox(PMMG_pParMesh parmesh,double min[3],
                                     double max[3],double *delta) {
  PMMG_coorCell *coor_list;
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_pTetra    pt;
  MMG5_pPoint    ppt;
  double         dd;
  int            ier,nitem,iel,ifac,*intvalues,k,i,j,l,idx,ip;

  ier = 0;

  nitem = parmesh->int_face_comm->nitem;

  PMMG_CALLOC(parmesh,coor_list,3*nitem,PMMG_coorCell,
              "node coordinates",return 0);
  PMMG_CALLOC(parmesh,parmesh->int_face_comm->intvalues,nitem,int,
              "face communicator",goto end);
  intvalues = parmesh->int_face_comm->intvalues;

  /* Store the node coordinates of the faces */
  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp  = &parmesh->listgrp[k];
    mesh = grp->mesh;
    assert ( mesh->info.delta &&  "missing scaling infos");

    dd = mesh->info.min[0]*mesh->info.min[0]
      + mesh->info.min[1]*mesh->info.min[1]
      + mesh->info.min[2]*mesh->info.min[2];

    assert ( fabs(mesh->info.delta-1.)<_MMG5_EPSD &&
             "scaled mesh... need to unscale it");

    for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
      iel   =  grp->face2int_face_comm_index1[i]/12;
      ifac  = (grp->face2int_face_comm_index1[i]%12)/3;

      assert ( iel && iel<=mesh->ne );
      assert ( 0<=ifac && ifac<4 );
      pt = &mesh->tetra[iel];

      idx = grp->face2int_face_comm_index2[i];

      if ( !intvalues[idx] ) {
        for ( l=0; l<3; ++l ) {
          ip = pt->v[_MMG5_idir[ifac][l]];
          assert ( ip && ip<=mesh->np );
          ppt = &mesh->point[ip];

          for ( j=0; j<3; ++j )
            coor_list[3*idx].c[j] = ppt->c[j];
        }
      }
      else intvalues[idx] = 1;
    }
  }

  /* Find the bounding box of the internal comm */
  if ( !PMMG_find_coorCellListBoundingBox(coor_list,nitem,min,max,
                                          delta) ) goto end;

  /* Success */
  ier = 1;

  end:
  PMMG_DEL_MEM(parmesh,parmesh->int_face_comm->intvalues,nitem,int,
               "face communicator");
  PMMG_DEL_MEM(parmesh,coor_list,3*nitem,PMMG_coorCell,"coor_list array");

  return ier;
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * \return 1 if success, 0 if fail
 *
 * Check the internal node communicator consitency by comparing the coordinates
 * of the vertices of the nodes that are in the same position in the
 * communicator. By the same time, check that we don't have more than 2 nodes
 * per position in the node communicator.
 *
 * For each group, check the unicity of the position in the internal
 * communicator.
 *
 */
int PMMG_check_intNodeComm( PMMG_pParMesh parmesh )
{
  PMMG_coorCell *coor_list;
  PMMG_pGrp     grp;
  MMG5_pMesh    mesh;
  double        dd,bb_min[3],bb_max[3],delta,dist[3],dist_norm;
  int ier;
  int ngrp = parmesh->ngrp;
  int nitem = parmesh->int_node_comm->nitem;
  int commIdx2 = 0;
  int commIdx1 = 0;
  int commSizeLoc = 0;
  int commIdx, k,j;

  if ( !parmesh->int_node_comm->nitem ) return 1;

  ier = 0;

  coor_list = NULL;

  /** Step 1: Find the internal communicator bounding box */
  if ( !PMMG_find_intNodeCommBoundingBox(parmesh,bb_min,bb_max,&delta) )
    return 0;

  /** Step 2: Fill coor_list (of size in_node_comm->nitem) with the (scaled)
   * coordinates of the boundary points. If the position is already used, check
   * the coordinates consistency and check that it is used only 2 times. */
  PMMG_CALLOC(parmesh,coor_list,nitem,PMMG_coorCell,"node coordinates",
              goto end);

  dd = 1./delta;
  for ( commIdx = 0; commIdx<nitem; commIdx++ )
    coor_list[commIdx].idx = PMMG_UNSET;

  for ( k = 0; k < ngrp; ++k ) {

    grp = parmesh->listgrp + k;
    mesh = grp->mesh;

    commSizeLoc = grp->nitem_int_node_comm;
    for ( commIdx = 0; commIdx < commSizeLoc; ++commIdx ) {

      commIdx2 = grp->node2int_node_comm_index2[ commIdx ];
      commIdx1 = grp->node2int_node_comm_index1[ commIdx ];

      if ( coor_list[ commIdx2 ].idx<0 ) {
        for ( j=0; j<3; ++j )
          coor_list[commIdx2].c[j] = dd*(mesh->point[commIdx1].c[j]-bb_min[j]);
        coor_list[ commIdx2 ].idx = commIdx;
      } else {

        dist_norm = 0.;
        for ( j = 0; j < 3; ++j ) {
          dist[j] = coor_list[ commIdx2 ].c[j]
            -dd*(mesh->point[ commIdx1 ].c[j]-bb_min[j]);
          dist_norm += dist[j]*dist[j];
        }

        if ( dist_norm>_MMG5_EPSD ) {
          fprintf(stderr,"  ## Error: %s: rank %d: group %d:\n"
                 "       2 different points (dist %e) in the same position (%d)"
                 " of the internal communicator:\n"
                 "       - point at position %d : %e %e %e\n"
                 "       - point at position %d : %e %e %e (point id=%d)\n",
                 __func__,parmesh->myrank,k,dist_norm,commIdx2,
                 coor_list[ commIdx2 ].idx,coor_list[ commIdx2 ].c[0],
                 coor_list[ commIdx2 ].c[1],coor_list[ commIdx2 ].c[2],
                 commIdx,dd*(mesh->point[ commIdx1 ].c[0]-bb_min[0]),
                 dd*(mesh->point[ commIdx1 ].c[1]-bb_min[1]),
                 dd*(mesh->point[ commIdx1 ].c[2]-bb_min[2]),commIdx1);
          goto end;
        }
      }
    }
  }

  /** Step 3: check that for a given group, each face has a unique position in
   * the internal node communicator */
  for ( k = 0; k < ngrp; ++k ) {
    grp = parmesh->listgrp + k;
    mesh = grp->mesh;
    commSizeLoc = grp->nitem_int_node_comm;

    for ( commIdx = 0; commIdx < nitem; ++commIdx )
      coor_list[ commIdx ].idx = -1;

    for ( commIdx = 0; commIdx < commSizeLoc; ++commIdx ) {
      commIdx2 = grp->node2int_node_comm_index2[ commIdx ];
      if ( coor_list[ commIdx2 ].idx>=0 ) {
        fprintf(stderr,"  ## Error: %s: rank %d: group %d:\n"
                "       2 points of the same node2int_node_comm array (%d %d)"
                " have the same position in the internal communicator (%d) ",
                __func__,parmesh->myrank,k,commIdx,coor_list[ commIdx2 ].idx,
                commIdx2);
        goto end;
      }
      assert ( coor_list[ commIdx2 ].idx<0 );

      coor_list[ commIdx2 ].idx = commIdx;
    }
  }
  /* Success */
  ier = 1;

end:
  PMMG_DEL_MEM(parmesh,coor_list,nitem,PMMG_coorCell,"coor_list array");

  return ier;
}


/**
 * \param parmesh pointer toward a parmesh structure
 *
 * \return 1 if success, 0 if fail
 *
 * Check the internal face communicator consitency by comparing the coordinates
 * of the vertices of the faces that are in the same position in the
 * communicator. Then check that we don't have more than 2 faces
 * per position in the face communicator.
 * Finally, for each group, check the unicity of the position in the internal
 * communicator.
 *
 */
int PMMG_check_intFaceComm( PMMG_pParMesh parmesh ) {
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_pTetra    pt;
  MMG5_pPoint    ppt;
  double         delta,dd,*doublevalues,dist[3],dist_norm,bb_min[3],bb_max[3];
  int            *intvalues;
  int            k,i,j,l,iel,ifac,iploc,ip,idx,idx_ori,nitem,ier;

  if ( !parmesh->int_face_comm->nitem ) return 1;

  ier = 0;

  /** Step 0: Find the internal communicator bounding box */
  if ( !PMMG_find_intFaceCommBoundingBox(parmesh,bb_min,bb_max,&delta) )
    return 0;

  /** Step 1: Fill int_face_comm->doublevalues with the (scaled) coordinates
   * of the boundary points. If the position is already used, check the
   * coordinates consistency and check that it is used only 2 times. */
  nitem = parmesh->int_face_comm->nitem;
  PMMG_CALLOC(parmesh,parmesh->int_face_comm->doublevalues,9*nitem,double,
              "face communicator",goto end);
  doublevalues = parmesh->int_face_comm->doublevalues;
  PMMG_MALLOC(parmesh,parmesh->int_face_comm->intvalues,nitem,int,
              "face communicator",goto end);
  intvalues = parmesh->int_face_comm->intvalues;

  for ( i=0; i<parmesh->int_face_comm->nitem; ++i )
    intvalues[i] = PMMG_UNSET;

  dd = 1./delta;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp  = &parmesh->listgrp[k];
    mesh = parmesh->listgrp[k].mesh;

    for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
      iel   =  grp->face2int_face_comm_index1[i]/12;
      ifac  = (grp->face2int_face_comm_index1[i]%12)/3;
      iploc = (grp->face2int_face_comm_index1[i]%12)%3;

      assert ( iel && iel<=mesh->ne );
      assert ( 0<=ifac && ifac<4 );
      pt = &mesh->tetra[iel];

      idx = grp->face2int_face_comm_index2[i];

      if ( intvalues[idx]>=0 ) {
        for ( l=0; l<3; ++l ) {
          /* List the face points from iploc to ensure that we start from the
           * same face point for the 2 mathing comminicators */
          ip = pt->v[_MMG5_idir[ifac][(l+iploc)%3]];

          assert ( ip && ip<=mesh->np );
          ppt = &mesh->point[ip];

          /* Travel the face in opposite direction for the two groups */
          switch(l)
          {
          case 0:
            idx_ori = 9*idx;
            break;
          case 1:
            idx_ori = 9*idx+6;
            break;
          case 2:
            idx_ori = 9*idx+3;
            break;
          }

          dist_norm = 0.;
          for ( j=0; j<3; ++j ) {
            dist[j] = doublevalues[idx_ori+j]-dd*(ppt->c[j]-bb_min[j]);
            dist_norm += dist[j]*dist[j];
          }

          if ( dist_norm > _MMG5_EPSD ) {
            fprintf(stderr,"  ## Error: %s: rank %d: group %d:\n"
                    "       2 different points (dist %e) in the same position (%d)"
                    " of the internal communicator:\n"
                    "       - face position %d : %e %e %e\n"
                    "       - face position %d : %e %e %e (elt %d, face %d)\n",
                    __func__,parmesh->myrank,k,dist_norm,idx,i,ppt->c[0],
                    ppt->c[1],ppt->c[2],intvalues[idx],doublevalues[idx_ori],
                    doublevalues[idx_ori+1],doublevalues[idx_ori+2],iel,ifac);
            goto end;
          }
        }
      }
      else {
        intvalues[idx] = i;
        for ( l=0; l<3; ++l ) {
          /* List the face points from iploc to ensure that we start from the
           * same face point for the 2 mathing comminicators */
          ip = pt->v[_MMG5_idir[ifac][(l+iploc)%3]];

          assert ( ip && ip<=mesh->np );
          ppt = &mesh->point[ip];

          for ( j=0; j<3; ++j )
            doublevalues[9*idx+3*l+j] = dd*(ppt->c[j]-bb_min[j]);
        }
      }
    }
  }

  /** Step 2: count the number of faces that point to the same position in the
   * internal communicator and check that we do not have more than 2 faces per
   * position. */
  for ( i=0; i<parmesh->int_face_comm->nitem; ++i )
    intvalues[i] = 0;

  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp  = &parmesh->listgrp[k];

    for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
      idx = grp->face2int_face_comm_index2[i];

      if ( ++intvalues[idx]>2 ) {
        fprintf(stderr,"  ## Error: %s: rank %d: more than 2 faces stored in"
                " the same position (%d) of the nodal communicator.\n",
                __func__,parmesh->myrank,idx );
        goto end;
      }
    }
  }

  /** Step 3: check that for a given group, each face has a unique position in
   * the internal face communicator */
  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp  = &parmesh->listgrp[k];
    mesh = parmesh->listgrp[k].mesh;

    for ( i=0; i<parmesh->int_face_comm->nitem; ++i )
      intvalues[i] = PMMG_UNSET;

    /* Reset the tetra flag, it will be used to store if a tetra face as been
     * listed in the communicator */
    for ( i=1; i<=mesh->ne; ++i )
      mesh->tetra[i].flag = 0;

    for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
      iel   =  grp->face2int_face_comm_index1[i]/12;
      ifac  = (grp->face2int_face_comm_index1[i]%12)/3;
      idx = grp->face2int_face_comm_index2[i];
      if ( intvalues[idx]>=0 ) {
        fprintf(stderr,"  ## Error: %s: rank %d: group %d:\n"
                "       2 different faces in the same position (%d)"
                " of the internal communicator:\n"
                "       - elt %d, face %d (position %d)\n"
                "       - elt %d, face %d \n",
                __func__,parmesh->myrank,k,idx,iel,ifac,i,intvalues[idx]/12,
                (intvalues[idx]%12)/3 );
        goto end;
      }
      intvalues[idx] = grp->face2int_face_comm_index1[i];

      pt = &mesh->tetra[iel];
      if ( MG_GET(pt->flag,ifac ) ) {
        printf("  ## Error: %s: grp %d: face %d of the tetra %d is stored twice"
               " in the internal communicator\n",__func__,k,ifac,iel);
        goto end;

      }
      MG_SET ( pt->flag,ifac );
    }
  }

  /* Success */
  ier = 1;

end:
  PMMG_DEL_MEM(parmesh,parmesh->int_face_comm->doublevalues,9*nitem,double,
               "face communicator");
  PMMG_DEL_MEM(parmesh,parmesh->int_face_comm->intvalues,nitem,int,
               "face communicator");

  return ier;
}


/**
 * \param parmesh pointer to current parmesh stucture
 *
 * \return 0 if fail, 1 otherwise
 *
 * Check the external node communicators consitency by comparing their size and
 * the coordinates of the listed points.
 *
 */
int PMMG_check_extNodeComm( PMMG_pParMesh parmesh )
{
  PMMG_pext_comm ext_node_comm;
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_pPoint    ppt;
  MPI_Request    *request;
  MPI_Status     *status;
  double         *rtosend,*rtorecv,*doublevalues,x,y,z,bb_min[3],bb_max[3];
  double         dd,delta,delta_all,bb_min_all[3];
  int            *r2send_size,*r2recv_size,color;
  int            k,i,j,ip,idx,ireq,nitem,nitem_color_out,ier;

  if ( !parmesh->next_node_comm ) return 1;

  ier = 0;
  r2send_size = NULL;
  r2recv_size = NULL;
  request     = NULL;
  status      = NULL;

  /** Step 1: Find the internal communicator bounding box */
  if ( !PMMG_find_intNodeCommBoundingBox(parmesh,bb_min,bb_max,&delta) )
    return 0;
  MPI_Allreduce ( &delta,&delta_all,1,MPI_DOUBLE,MPI_MAX,parmesh->comm);
  MPI_Allreduce ( bb_min,bb_min_all,3,MPI_DOUBLE,MPI_MIN,parmesh->comm);

  /** Step 2: Fill int_node_comm->doublevalues with the coordinates
   * of the boundary points */
  nitem = parmesh->int_node_comm->nitem;
  PMMG_CALLOC(parmesh,parmesh->int_node_comm->doublevalues,3*nitem,double,
              "node communicator",goto end);
  doublevalues = parmesh->int_node_comm->doublevalues;

  dd = 1./delta_all;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp  = &parmesh->listgrp[k];
    mesh = parmesh->listgrp[k].mesh;

    for ( i=0; i<grp->nitem_int_node_comm; ++i ) {
      ip  = grp->node2int_node_comm_index1[i];
      idx = grp->node2int_node_comm_index2[i];

      assert ( ip && ip<=mesh->np );

      ppt = &mesh->point[ip];
      for ( j=0; j<3; ++j )
        doublevalues[3*idx+j] = dd * (ppt->c[j] - bb_min_all[j]);
    }
  }

  /** Step 3: Send the values that need to be communicate to the suitable
   * processor */
  PMMG_MALLOC(parmesh,request,2*parmesh->next_node_comm,MPI_Request,
              "mpi request array",goto end);

  PMMG_MALLOC(parmesh,status,2*parmesh->next_node_comm,MPI_Status,
              "mpi status array",goto end);

  PMMG_CALLOC(parmesh,r2send_size,parmesh->next_node_comm,int,
              "size of the r2send array",goto end);

  ireq= 0;
  for ( k=0; k<parmesh->next_node_comm; ++k ) {

    ext_node_comm = &parmesh->ext_node_comm[k];

    if ( r2send_size[k] < 3*ext_node_comm->nitem ) {
      PMMG_REALLOC(parmesh,ext_node_comm->rtosend,3*ext_node_comm->nitem,
                   r2send_size[k],double,"rtosend",goto end);
      r2send_size[k] = 3*ext_node_comm->nitem;
    }

    /* Filling of the array to send */
    rtosend = ext_node_comm->rtosend;
    color   = ext_node_comm->color_out;
    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx = ext_node_comm->int_comm_index[i];
      for ( j=0; j<3; ++j )
        rtosend[3*i+j] = doublevalues[3*idx+j];
    }

    request[ireq]    = MPI_REQUEST_NULL;
    MPI_CHECK( MPI_Isend(&ext_node_comm->nitem,1,MPI_INT,color,
                         MPI_CHKCOMM_NODE_TAG,
                         parmesh->comm,&request[ireq++]),goto end );

    request[ireq]    = MPI_REQUEST_NULL;
    MPI_CHECK( MPI_Isend(rtosend,3*ext_node_comm->nitem,MPI_DOUBLE,color,
                         MPI_CHKCOMM_NODE_TAG+1,
                         parmesh->comm,&request[ireq++]),
               goto end );
  }

  /** Step 4: Recv the values from the senders and check:
   *
   * - that the communicators have the same size on both processors
   *
   * - that the coordinates of the points listed in the communicators are
         similar at epsilon machine
   */
  PMMG_CALLOC(parmesh,r2recv_size,parmesh->next_node_comm,int,
              "size of the r2recv array",goto end);
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    color         = ext_node_comm->color_out;

    MPI_CHECK( MPI_Recv(&nitem_color_out,1,MPI_INT,color,
                        MPI_CHKCOMM_NODE_TAG,parmesh->comm,
                        &status[0]), goto end );

    /* Check the size of the communicators */
    if ( nitem_color_out != ext_node_comm->nitem ) {
      fprintf(stderr,"  ## Error: %s: rank %d: the size of the external"
              " communicator %d->%d (%d) doesn't match with the size of the same"
              " external communicator on %d (%d)\n",__func__,parmesh->myrank,
              parmesh->myrank,color,ext_node_comm->nitem,color,nitem_color_out );
      goto end;
    }

    if ( r2recv_size[k] < 3*nitem_color_out ) {
      PMMG_REALLOC(parmesh,ext_node_comm->rtorecv,3*nitem_color_out,
                   r2recv_size[k],double,"rtorecv",goto end);
      r2recv_size[k] = 3*nitem_color_out;
    }
    rtorecv       = ext_node_comm->rtorecv;
    MPI_CHECK( MPI_Recv(rtorecv,3*nitem_color_out,MPI_DOUBLE,color,
                        MPI_CHKCOMM_NODE_TAG+1,parmesh->comm,
                        &status[0]), goto end );

    /* Check the values of the node in the communicator */
    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx = ext_node_comm->int_comm_index[i];

      x   = doublevalues[3*idx  ] - rtorecv[3*i  ];
      y   = doublevalues[3*idx+1] - rtorecv[3*i+1];
      z   = doublevalues[3*idx+2] - rtorecv[3*i+2];

      if ( x*x + y*y + z*z > _MMG5_EPSD ) {
        fprintf(stderr,"  ## Error: %s: rank %d:\n"
                "       2 different points (dist %e) in the same position (%d)"
                " of the external communicator %d %d (%d th item):\n"
                "       - point : %e %e %e\n"
                "       - point : %e %e %e\n",__func__,parmesh->myrank,
                x*x+y*y+z*z,idx,parmesh->myrank,color,i,doublevalues[3*idx],
                doublevalues[3*idx+1],doublevalues[3*idx+2],rtorecv[3*i],
                rtorecv[3*i+1],rtorecv[3*i+2]);
        goto end;
      }
    }
  }
  MPI_CHECK( MPI_Waitall(2*parmesh->next_node_comm,request,status), goto end );

  /* Success */
  ier = 1;

end:
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];

    // Change this and add to the external comm the possibility to not
    // unalloc/realloc every time, thus, here, we will be able to reset the
    // communicators without unallocated it
    PMMG_DEL_MEM(parmesh,ext_node_comm->rtosend,r2send_size[k],double,"r2send");
    PMMG_DEL_MEM(parmesh,ext_node_comm->rtorecv,r2recv_size[k],double,"r2recv");
  }

  PMMG_DEL_MEM(parmesh,r2send_size,parmesh->next_node_comm,int,
               "size of the r2send array");

  PMMG_DEL_MEM(parmesh,r2recv_size,parmesh->next_node_comm,int,
               "size of the r2send array");

  PMMG_DEL_MEM(parmesh,status,2*parmesh->next_node_comm,MPI_Status,
              "mpi status array");

  PMMG_DEL_MEM(parmesh,request,2*parmesh->next_node_comm,MPI_Request,
              "mpi request array");

  PMMG_DEL_MEM(parmesh,parmesh->int_node_comm->doublevalues,3*nitem,double,
               "node communicator");

  return ier;
}


/**
 * \param parmesh pointer to current parmesh stucture
 *
 * \return 0 if fail, 1 otherwise
 *
 * Check the external face communicators consitency by comparing their size and
 * the coordinates of the face points.
 *
 */
int PMMG_check_extFaceComm( PMMG_pParMesh parmesh )
{
  PMMG_pext_comm ext_face_comm;
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_pTetra    pt;
  MMG5_pPoint    ppt;
  MPI_Request    *request;
  MPI_Status     *status;
  double         *rtosend,*rtorecv,*doublevalues,dd,x,y,z;
  double         bb_min[3],bb_max[3],delta,delta_all,bb_min_all[3];
  int            *r2send_size,*r2recv_size,color;
  int            k,i,j,l,ireq,ip,iploc,iel,ifac,idx,nitem,nitem_color_out,ier;

  if ( !parmesh->next_face_comm ) return 1;

  ier = 0;

  r2send_size = NULL;
  r2recv_size = NULL;
  request     = NULL;
  status      = NULL;

  /** Step 0: Find the internal communicator bounding box */
  if ( !PMMG_find_intFaceCommBoundingBox(parmesh,bb_min,bb_max,&delta) )
    return 0;
  MPI_Allreduce ( &delta,&delta_all,1,MPI_DOUBLE,MPI_MAX,parmesh->comm);
  MPI_Allreduce ( bb_min,bb_min_all,3,MPI_DOUBLE,MPI_MIN,parmesh->comm);

  /** Step 1: Fill int_face_comm->doublevalues with the scaled coordinates
   * of the boundary points of the face */
  nitem = parmesh->int_face_comm->nitem;
  PMMG_CALLOC(parmesh,parmesh->int_face_comm->doublevalues,9*nitem,double,
              "face communicator",goto end);
  doublevalues = parmesh->int_face_comm->doublevalues;

  dd = 1./delta_all;
  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp  = &parmesh->listgrp[k];
    mesh = parmesh->listgrp[k].mesh;
    assert ( mesh->info.delta && "missing scaling infos");

    for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
      idx  = grp->face2int_face_comm_index2[i];

      iel   =  grp->face2int_face_comm_index1[i]/12;
      ifac  = (grp->face2int_face_comm_index1[i]%12)/3;
      iploc = (grp->face2int_face_comm_index1[i]%12)%3;

      assert ( iel && iel<=mesh->ne );
      assert ( 0<=ifac  && ifac <4 );
      assert ( 0<=iploc && iploc<3 );

      pt = &mesh->tetra[iel];

      for ( l=0; l<3; ++l ) {
        /* List the face points from iploc to ensure that we start from the same
         * face point for the 2 mathing comminicators */
        ip  = pt->v[_MMG5_idir[ifac][(l+iploc)%3]];

        assert ( ip && ip<=mesh->np );
        ppt = &mesh->point[ip];

        for ( j=0; j<3; ++j )
          doublevalues[9*idx+3*l+j]=dd*(ppt->c[j]-bb_min_all[j]);
      }
    }
  }

  /** Step 2: Send the values that need to be communicate to the suitable
   * processor */
  PMMG_MALLOC(parmesh,request,2*parmesh->next_face_comm,MPI_Request,
              "mpi request array",goto end);

  PMMG_MALLOC(parmesh,status,2*parmesh->next_face_comm,MPI_Status,
              "mpi status array",goto end);

  PMMG_CALLOC(parmesh,r2send_size,parmesh->next_face_comm,int,
              "size of the r2send array",goto end);

  ireq = 0;
  for ( k=0; k<parmesh->next_face_comm; ++k ) {

    ext_face_comm = &parmesh->ext_face_comm[k];

    if ( r2send_size[k] < 9*ext_face_comm->nitem ) {
      PMMG_REALLOC(parmesh,ext_face_comm->rtosend,9*ext_face_comm->nitem,
                   r2send_size[k],double,"rtosend",goto end);
      r2send_size[k] = 9*ext_face_comm->nitem;
    }

    /* Filling of the array to send */
    rtosend = ext_face_comm->rtosend;
    color   = ext_face_comm->color_out;
    for ( i=0; i<ext_face_comm->nitem; ++i ) {
      idx = ext_face_comm->int_comm_index[i];

      for ( j=0; j<9; ++j )
        rtosend[9*i+j] = doublevalues[9*idx+j];
    }

    request[ireq]    = MPI_REQUEST_NULL;
    MPI_CHECK( MPI_Isend(&ext_face_comm->nitem,1,MPI_INT,color,
                         MPI_CHKCOMM_FACE_TAG,
                         parmesh->comm,&request[ireq++]),goto end );

    request[ireq]    = MPI_REQUEST_NULL;
    MPI_CHECK( MPI_Isend(rtosend,9*ext_face_comm->nitem,MPI_DOUBLE,color,
                         MPI_CHKCOMM_FACE_TAG+1,
                         parmesh->comm,&request[ireq++]),
               goto end );
  }

  /** Step 3: Recv the values from the senders and check:
   *
   * - that the communicators have the same size on both processors
   *
   * - that the coordinates of the points listed in the communicators are
         similar at epsilon machine
   */
  PMMG_CALLOC(parmesh,r2recv_size,parmesh->next_face_comm,int,
              "size of the r2recv array",goto end);
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    color         = ext_face_comm->color_out;

    MPI_CHECK( MPI_Recv(&nitem_color_out,1,MPI_INT,color,
                        MPI_CHKCOMM_FACE_TAG,parmesh->comm,
                        &status[0]), goto end );

    /* Check the size of the communicators */
    if ( nitem_color_out != ext_face_comm->nitem ) {
      printf("  ## Error: %s: the size of the external communicator %d->%d"
             " doesn't match with the size of the same external communicator"
             " on %d (%d versus %d)\n",__func__,parmesh->myrank,color,color,
             ext_face_comm->nitem,nitem_color_out);
      goto end;
    }

    if ( r2recv_size[k] < 9*nitem_color_out ) {
      PMMG_REALLOC(parmesh,ext_face_comm->rtorecv,9*nitem_color_out,
                   r2recv_size[k],double,"rtorecv",goto end);
      r2recv_size[k] = 9*nitem_color_out;
    }
    rtorecv       = ext_face_comm->rtorecv;
    MPI_CHECK( MPI_Recv(rtorecv,9*nitem_color_out,MPI_DOUBLE,color,
                        MPI_CHKCOMM_FACE_TAG+1,parmesh->comm,
                        &status[0]), goto end );

    /* Check the values of the face in the communicator */
    for ( i=0; i<ext_face_comm->nitem; ++i ) {
      idx = ext_face_comm->int_comm_index[i];

      /* Travel the face in opposite direction for the two communicators */
      x   = doublevalues[9*idx  ] - rtorecv[9*i  ];
      y   = doublevalues[9*idx+1] - rtorecv[9*i+1];
      z   = doublevalues[9*idx+2] - rtorecv[9*i+2];
      if ( x*x + y*y + z*z > _MMG5_EPSD ) {
        printf("  ## Error: %s: item %d of the external communicator %d->%d:\n"
               "                     vertex %d: %e %e %e -- %e %e %e"
               " (dist = %e)\n",__func__,i,l,
               parmesh->myrank,color,doublevalues[9*idx],doublevalues[9*idx+1],
               doublevalues[9*idx+2],rtorecv[9*i],rtorecv[9*i+1],rtorecv[9*i+2],
               x*x + y*y + z*z);
        goto end;
      }

      x   = doublevalues[9*idx+3  ] - rtorecv[9*i+6  ];
      y   = doublevalues[9*idx+3+1] - rtorecv[9*i+6+1];
      z   = doublevalues[9*idx+3+2] - rtorecv[9*i+6+2];
      if ( x*x + y*y + z*z > _MMG5_EPSD ) {
        printf("  ## Error: %s: item %d of the external communicator %d->%d:\n"
               "                     vertex %d: %e %e %e -- %e %e %e"
               " (dist = %e)\n",__func__,i,l,
               parmesh->myrank,color,doublevalues[9*idx+3],doublevalues[9*idx+4],
               doublevalues[9*idx+5],rtorecv[9*i+6],rtorecv[9*i+7],rtorecv[9*i+8],
               x*x + y*y + z*z);
        goto end;
      }

      x   = doublevalues[9*idx+6  ] - rtorecv[9*i+3  ];
      y   = doublevalues[9*idx+6+1] - rtorecv[9*i+3+1];
      z   = doublevalues[9*idx+6+2] - rtorecv[9*i+3+2];
      if ( x*x + y*y + z*z > _MMG5_EPSD ) {
        printf("  ## Error: %s: item %d of the external communicator %d->%d:\n"
               "                     vertex %d: %e %e %e -- %e %e %e"
               " (dist = %e)\n",__func__,i,l,
               parmesh->myrank,color,doublevalues[9*idx+6],doublevalues[9*idx+7],
               doublevalues[9*idx+8],rtorecv[9*i+3],rtorecv[9*i+4],rtorecv[9*i+5],
               x*x + y*y + z*z);
        goto end;
      }
    }
  }
  MPI_CHECK( MPI_Waitall(2*parmesh->next_face_comm,request,status), goto end );

  /* Success */
  ier = 1;

end:
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];

    // Change this and add to the external comm the possibility to not
    // unalloc/realloc every time, thus, here, we will be able to reset the
    // communicators without unallocated it
    PMMG_DEL_MEM(parmesh,ext_face_comm->rtosend,r2send_size[k],double,"r2send");
    PMMG_DEL_MEM(parmesh,ext_face_comm->rtorecv,r2recv_size[k],double,"r2recv");
  }

  PMMG_DEL_MEM(parmesh,r2send_size,parmesh->next_face_comm,int,
               "size of the r2send array");

  PMMG_DEL_MEM(parmesh,r2recv_size,parmesh->next_face_comm,int,
               "size of the r2send array");

  PMMG_DEL_MEM(parmesh,status,2*parmesh->next_face_comm,MPI_Status,
              "mpi status array");

  PMMG_DEL_MEM(parmesh,request,2*parmesh->next_face_comm,MPI_Request,
              "mpi request array");

  PMMG_DEL_MEM(parmesh,parmesh->int_face_comm->doublevalues,9*nitem,double,
               "face communicator");

  return ier;
}

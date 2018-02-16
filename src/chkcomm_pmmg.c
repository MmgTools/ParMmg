#include "parmmg.h"
#include "coorcell_pmmg.h"

/**
 * \param parmesh pointer to current parmesh stucture
 *
 * \return 1 if tests succeeded
 *         0 if some test failed
 *
 *  This function does some checks on the internal communicator for each group
 *
 *    First test:
 *        allocate an array of int_node_comm->item_int_node struct point
 *            (ie three c[3]) elements (initialized to 0)
 *        for each grp communicator
 *            if array[indx2] is empty
 *                initialize to current group's point's coordinates
 *            else
 *                check that the existing point's coordinates match the current
 *                group's point's coordinates
 *      This test will detect mismatches between groups but will NOT detect:
 *          errors in elements that are not shared between the subgroups, ie
 *            elements interfacing groups on other processors
 *          if the same error happens in both groups, ie if the values in
 *            both groups match still they are are both wrong
 *
 *    Second test:
 *        elements should appear once and only once per mesh local communicator
 */

int PMMG_check_intNodeComm( PMMG_pParMesh parmesh )
{
  //NIKOS TODO: could try using std::vector<bool> and std::map in these tests to speed up things
  //NIKOS TODO: BOUND CHECK ARRAY ACCESSES: can't rely that the contents of the communicators are valid
  PMMG_coorCell   *coor_list;
  double dd,dd_mesh,bb_min[3],bb_max[3],delta,dist[3];
  int ret_val = 1;
  int ngrp = parmesh->ngrp;
  PMMG_pGrp grpCur = NULL;
  MMG5_pMesh meshCur = NULL;
  int commSizeGlo = parmesh->int_node_comm->nitem;
  int commIdx2 = 0;
  int commIdx1 = 0;
  int numFailed = 0;
  int commSizeLoc = 0;

  // Loop counter(s)
  int commIdx, grpId,i,j,ip,idx;

  // FIRST TEST:
  PMMG_CALLOC(parmesh,coor_list,commSizeGlo,PMMG_coorCell,"node coordinates",
              return 0);

  /* Store the node coordinates of the nodes */
  for ( grpId=0; grpId<parmesh->ngrp; ++grpId ) {
    grpCur  = &parmesh->listgrp[grpId];
    meshCur = grpCur->mesh;
    assert ( meshCur->info.delta &&  "missing scaling infos");

    dd = meshCur->info.min[0]*meshCur->info.min[0]
      + meshCur->info.min[1]*meshCur->info.min[1]
      + meshCur->info.min[2]*meshCur->info.min[2];

    assert ( fabs(meshCur->info.delta-1.)<_MMG5_EPSD && dd<_MMG5_EPSD &&
             "scaled mesh... need to unscale it");

    for ( i=0; i<grpCur->nitem_int_node_comm; ++i ) {
      ip      = grpCur->node2int_node_comm_index1[i];
      idx     = grpCur->node2int_node_comm_index2[i];
      for ( j=0; j<3; ++j )
        coor_list[idx].coor[j] = meshCur->point[ip].c[j];
    }
  }

  for ( commIdx = 0; commIdx<commSizeGlo; commIdx++ )
    coor_list[commIdx].pos = 0;

  /* Scale the coordinates depending to the bounding box ofthe internal comm */
  if ( !PMMG_scale_coorCellList(coor_list,commSizeGlo,bb_min,bb_max,&delta) )
    return 0;

  dd = 1./delta;

  for ( grpId = 0; grpId < ngrp; ++grpId ) {

    grpCur = parmesh->listgrp + grpId;
    meshCur = grpCur->mesh;
    numFailed = 0;

    assert ( meshCur->info.delta &&  "missing scaling infos");

    dd_mesh = meshCur->info.min[0]*meshCur->info.min[0]
      + meshCur->info.min[1]*meshCur->info.min[1]
      + meshCur->info.min[2]*meshCur->info.min[2];

    assert ( fabs(meshCur->info.delta-1.)<_MMG5_EPSD && dd_mesh<_MMG5_EPSD &&
             "scaled mesh... need to unscale it");



    commSizeLoc = grpCur->nitem_int_node_comm;
    for ( commIdx = 0; commIdx < commSizeLoc; ++commIdx ) {

      commIdx2 = grpCur->node2int_node_comm_index2[ commIdx ];
      commIdx1 = grpCur->node2int_node_comm_index1[ commIdx ];
//      printf( "+++++NIKOS[%d/%d]+++:: CIC: commSizeGlo = %d, i = %d, n2i_idx2(check position) = %d ",
//              grpId + 1, ngrp, commSizeGlo, commIdx, commIdx2);
      if ( !coor_list[ commIdx2 ].pos ) {
//        printf( "+++adding (%f,%f,%f) to (%f,%f,%f) \n",
//                meshCur->point[ commIdx1 ].c[0]-min[0]),
//                meshCur->point[ commIdx1 ].c[1]-min[1]),
//                meshCur->point[ commIdx1 ].c[2]-min[2]),
//                coor_list[ commIdx2 ].coor[0], coor_list[ commIdx2 ].coor[1],
//                coor_list[ commIdx2 ].coor[2] );
        coor_list[ commIdx2 ].coor[0] = dd*(meshCur->point[ commIdx1 ].c[0]-bb_min[0]);
        coor_list[ commIdx2 ].coor[1] = dd*(meshCur->point[ commIdx1 ].c[1]-bb_min[1]);
        coor_list[ commIdx2 ].coor[2] = dd*(meshCur->point[ commIdx1 ].c[2]-bb_min[2]);
        coor_list[ commIdx2 ].pos = 1;
      } else {
//        printf( "---checking (%f,%f,%f) to (%f,%f,%f) \n",
//                dd*(meshCur->point[ commIdx1 ].c[0]-min[0]),
//                dd*(meshCur->point[ commIdx1 ].c[1]-min[1]),
//                dd*(meshCur->point[ commIdx1 ].c[2]-min[2]),
//                coor_list[ commIdx2 ].coor[0], coor_list[ commIdx2 ].coor[1],
//                coor_list[ commIdx2 ].coor[2] );
        for ( int j = 0; j < 3; ++j )
          dist[j] = coor_list[ commIdx2 ].coor[j]
            -dd*(meshCur->point[ commIdx1 ].c[j]-bb_min[j]);

        if ( dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]>_MMG5_EPSD ) {
          ++numFailed;
          ret_val = 0;
          assert(0);
        }
      }
    }
    if ( numFailed )
      printf( "+++++NIKOS[%d/%d]:: INTERNAL COMMUNICATOR CHECKED: number of failed communicator elements: %d \n", ngrp, grpId+1, numFailed / 3 );
  }


  // SECOND TEST:
  for ( grpId = 0; grpId < ngrp; ++grpId ) {
    grpCur = parmesh->listgrp + grpId;
    meshCur = grpCur->mesh;
    commSizeLoc = grpCur->nitem_int_node_comm;

    for ( commIdx = 0; commIdx < commSizeGlo; ++commIdx )
      coor_list[ commIdx ].pos = -1;

    for ( commIdx = 0; commIdx < commSizeLoc; ++commIdx ) {
      commIdx2 = grpCur->node2int_node_comm_index2[ commIdx ];
      assert ( coor_list[ commIdx2 ].pos<0 );

      coor_list[ commIdx2 ].pos = commIdx;
    }
  }

  // THIRD TEST:
  //   all elements in commSizeGlo
  //     a) should appear in per mesh local communicators at least once
  //     b) should have no "gaps"
//THINK AGAIN ABOUT IT AND FINISH IT  for ( commIdx = 0; commIdx < commSizeGlo + 1; ++commIdx )
//THINK AGAIN ABOUT IT AND FINISH IT    check[ commIdx ].tmp = 0;
//THINK AGAIN ABOUT IT AND FINISH IT  for ( grpId = 0; grpId < ngrp; ++grpId ) {
//THINK AGAIN ABOUT IT AND FINISH IT    grpCur = parmesh->listgrp + grpId;
//THINK AGAIN ABOUT IT AND FINISH IT    meshCur = grpCur->mesh;
//THINK AGAIN ABOUT IT AND FINISH IT    commSizeLoc = grpCur->nitem_int_node_comm;
//THINK AGAIN ABOUT IT AND FINISH IT    for ( commIdx = 0; commIdx < commSizeLoc; ++commIdx ) {
//THINK AGAIN ABOUT IT AND FINISH IT      commIdx2 = grpCur->node2int_node_comm_index2[ commIdx ];
//THINK AGAIN ABOUT IT AND FINISH IT      if ( check[ commIdx2 ].tmp == 0. )
//THINK AGAIN ABOUT IT AND FINISH IT        check[ commIdx2 ].tmp = 1.;
//THINK AGAIN ABOUT IT AND FINISH IT    }
//THINK AGAIN ABOUT IT AND FINISH IT  }
//THINK AGAIN ABOUT IT AND FINISH IT  commIdx = 0;
//THINK AGAIN ABOUT IT AND FINISH IT  // Skip non zero elements
//THINK AGAIN ABOUT IT AND FINISH IT  while ( check[ ++commIdx ].tmp != 0 )
//THINK AGAIN ABOUT IT AND FINISH IT    ;
//THINK AGAIN ABOUT IT AND FINISH IT  if ( commIdx != parmesh->int_node_comm->nitem )
//THINK AGAIN ABOUT IT AND FINISH IT    ret_val = 0;
//THINK AGAIN ABOUT IT AND FINISH IT  printf( "+++++NIKOS[%d/%d]:: commIdx= %d, nitem= %d \n", ngrp, grpId+1, commIdx, parmesh->int_node_comm->nitem );

  // NIKOS TODO: CHECK THE idx1/idx2 pairs ?

  PMMG_DEL_MEM(parmesh,coor_list,commSizeGlo+1,PMMG_coorCell,
               "Deallocating check space:");

  return ret_val;
}


/**
 * \param parmesh pointer toward a parmesh structure
 *
 * \return 1 if success, 0 if fail
 *
 * Check the internal face communicator consitency by comparing the coordinates
 * of the vertices of the faces that are in the same position in the
 * communicator. By the same time, check that we don't have more than 2 faces
 * per position in the face communicator.
 *
 * For each group, check the unicity of the position in the internal
 * communicator.
 *
 */

int PMMG_check_intFaceComm( PMMG_pParMesh parmesh ) {
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_pTetra    pt;
  MMG5_pPoint    ppt;
  double         coor[3],*doublevalues,x,y,z;
  int            *intvalues;
  int            k,i,j,l,iel,ifac,iploc,ip,idx,idx_ori,nitem,ier;

  ier = 0;

  /** Step 1: Fill int_face_comm->doublevalues with the (unscaled) coordinates
   * of the boundary points. If the position is already used, check the
   * coordinates consistency and check that it is used only 2 times. */
  nitem = parmesh->int_face_comm->nitem;
  PMMG_CALLOC(parmesh,parmesh->int_face_comm->doublevalues,9*nitem,double,
              "face communicator",goto end);
  doublevalues = parmesh->int_face_comm->doublevalues;
  PMMG_CALLOC(parmesh,parmesh->int_face_comm->intvalues,nitem,int,
              "face communicator",goto end);
  intvalues = parmesh->int_face_comm->intvalues;


  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp  = &parmesh->listgrp[k];
    mesh = parmesh->listgrp[k].mesh;
    assert ( mesh->info.delta && "missing scaling infos");

    for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
      iel   =  grp->face2int_face_comm_index1[i]/12;
      ifac  = (grp->face2int_face_comm_index1[i]%12)/3;
      iploc = (grp->face2int_face_comm_index1[i]%12)%3;

      assert ( iel && iel<=mesh->ne );
      assert ( 0<=ifac && ifac<4 );
      pt = &mesh->tetra[iel];

      idx = grp->face2int_face_comm_index2[i];

      for ( l=0; l<3; ++l ) {
        /* List the face points from iploc to ensure that we start from the same
         * face point for the 2 mathing comminicators */
        ip = pt->v[_MMG5_idir[ifac][(l+iploc)%3]];

        assert ( ip && ip<=mesh->np );
        ppt = &mesh->point[ip];

        for ( j=0; j<3; ++j )
          coor[j] = ppt->c[j]*mesh->info.delta + mesh->info.min[j];

        if ( intvalues[idx] ) {
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
          x = doublevalues[idx_ori  ]-coor[0];
          y = doublevalues[idx_ori+1]-coor[1];
          z = doublevalues[idx_ori+2]-coor[2];

          if ( x*x + y*y + z*z > _MMG5_EPSD ) {
            printf("  ## Error: %s: rank %d: item %d of the internal"
                   " communicator:\n %e %e %e -- %e %e %e (dist = %e)\n",
                   __func__,parmesh->myrank,idx,doublevalues[idx_ori],
                   doublevalues[idx_ori+1],doublevalues[idx_ori+2],coor[0],
                   coor[1],coor[2],x*x + y*y + z*z);
            goto end;
          }
        }
        else {
          for ( j=0; j<3; ++j )
            doublevalues[9*idx+3*l+j] = coor[j];
        }
      }
      if ( ++intvalues[idx]>2 ) {
        printf("  ## Error: %s: rank %d: more than 2 faces stored in the"
               " same position (%d) of the nodal communicator.\n",
               __func__,parmesh->myrank,idx );
        assert(0);
        goto end;
      }
    }
  }

  /** Step 2: check that for a given group, each face has a unique position in
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
        printf("  ## Error: %s: grp %d: 2 faces of the group stored in the"
               " same position (%d) of the nodal communicator: tetra %d (face %d)"
               " -- %d (face %d).\n",__func__,k,idx,iel,ifac,
               intvalues[idx]/12,(intvalues[idx]%12)/3 );
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
  double         *rtosend,*rtorecv,*doublevalues,x,y,z;
  int            *r2send_size,*r2recv_size,color;
  int            k,i,j,ip,idx,ireq,nitem,nitem_color_out,ier;

  ier = 0;

  /** Step 1: Fill int_node_comm->doublevalues with the (unscaled) coordinates
   * of the boundary points */
  nitem = parmesh->int_node_comm->nitem;
  PMMG_CALLOC(parmesh,parmesh->int_node_comm->doublevalues,3*nitem,double,
              "node communicator",goto end);
  doublevalues = parmesh->int_node_comm->doublevalues;

  for ( k=0; k<parmesh->ngrp; ++k ) {
    grp  = &parmesh->listgrp[k];
    mesh = parmesh->listgrp[k].mesh;
    assert ( mesh->info.delta && "missing scaling infos");

    for ( i=0; i<grp->nitem_int_node_comm; ++i ) {
      ip  = grp->node2int_node_comm_index1[i];
      idx = grp->node2int_node_comm_index2[i];

      assert ( ip && ip<=mesh->np );

      ppt = &mesh->point[ip];
      for ( j=0; j<3; ++j )
        doublevalues[3*idx+j] = ppt->c[j]*mesh->info.delta + mesh->info.min[j];
    }
  }

  /** Step 2: Send the values that need to be communicate to the suitable
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
    MPI_CHECK( MPI_Isend(&ext_node_comm->nitem,1,MPI_INT,color,0,
                         parmesh->comm,&request[ireq++]),goto end );

    request[ireq]    = MPI_REQUEST_NULL;
    MPI_CHECK( MPI_Isend(rtosend,3*ext_node_comm->nitem,MPI_DOUBLE,color,1,
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

  PMMG_CALLOC(parmesh,r2recv_size,parmesh->next_node_comm,int,
              "size of the r2recv array",goto end);
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    color         = ext_node_comm->color_out;

    MPI_CHECK( MPI_Recv(&nitem_color_out,1,MPI_INT,color,0,parmesh->comm,
                        &status[0]), goto end );

    /* Check the size of the communicators */
    if ( nitem_color_out != ext_node_comm->nitem ) {
      printf("  ## Error: %s: the size of the external communicator %d->%d"
             " doesn't match with the size of the same external communicator"
             " on %d\n",__func__,parmesh->myrank,color,color );
      goto end;
    }

    if ( r2recv_size[k] < 3*nitem_color_out ) {
      PMMG_REALLOC(parmesh,ext_node_comm->rtorecv,3*nitem_color_out,
                   r2recv_size[k],double,"rtorecv",goto end);
      r2recv_size[k] = 3*nitem_color_out;
    }
    rtorecv       = ext_node_comm->rtorecv;
    MPI_CHECK( MPI_Recv(rtorecv,3*nitem_color_out,MPI_DOUBLE,color,1,parmesh->comm,
                        &status[0]), goto end );

    /* Check the values of the node in the communicator */
    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx = ext_node_comm->int_comm_index[i];

      x   = doublevalues[3*idx  ] - rtorecv[3*i  ];
      y   = doublevalues[3*idx+1] - rtorecv[3*i+1];
      z   = doublevalues[3*idx+2] - rtorecv[3*i+2];

      if ( x*x + y*y + z*z > _MMG5_EPSD ) {
        printf("  ## Error: %s: item %d of the external communicator %d->%d:\n"
               "                    %e %e %e -- %e %e %e (dist = %e)\n",
               __func__,i,
               parmesh->myrank,color,doublevalues[3*idx],doublevalues[3*idx+1],
               doublevalues[3*idx+2],rtorecv[3*i],rtorecv[3*i+1],rtorecv[3*i+2],
               x*x + y*y + z*z);
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
  double         *rtosend,*rtorecv,*doublevalues,x,y,z;
  int            *r2send_size,*r2recv_size,color;
  int            k,i,j,l,ireq,ip,iploc,iel,ifac,idx,nitem,nitem_color_out,ier;

  ier = 0;

  /** Step 1: Fill int_face_comm->doublevalues with the (unscaled) coordinates
   * of the boundary points of the face */
  nitem = parmesh->int_face_comm->nitem;
  PMMG_CALLOC(parmesh,parmesh->int_face_comm->doublevalues,9*nitem,double,
              "face communicator",goto end);
  doublevalues = parmesh->int_face_comm->doublevalues;

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
          doublevalues[9*idx+3*l+j]=ppt->c[j]*mesh->info.delta+mesh->info.min[j];
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
    MPI_CHECK( MPI_Isend(&ext_face_comm->nitem,1,MPI_INT,color,0,
                         parmesh->comm,&request[ireq++]),goto end );

    request[ireq]    = MPI_REQUEST_NULL;
    MPI_CHECK( MPI_Isend(rtosend,9*ext_face_comm->nitem,MPI_DOUBLE,color,1,
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

    MPI_CHECK( MPI_Recv(&nitem_color_out,1,MPI_INT,color,0,parmesh->comm,
                        &status[0]), goto end );

    /* Check the size of the communicators */
    if ( nitem_color_out != ext_face_comm->nitem ) {
      printf("  ## Error: %s: the size of the external communicator %d->%d"
             " doesn't match with the size of the same external communicator"
             " on %d\n",__func__,parmesh->myrank,color,color );
      goto end;
    }

    if ( r2recv_size[k] < 9*nitem_color_out ) {
      PMMG_REALLOC(parmesh,ext_face_comm->rtorecv,9*nitem_color_out,
                   r2recv_size[k],double,"rtorecv",goto end);
      r2recv_size[k] = 9*nitem_color_out;
    }
    rtorecv       = ext_face_comm->rtorecv;
    MPI_CHECK( MPI_Recv(rtorecv,9*nitem_color_out,MPI_DOUBLE,color,1,parmesh->comm,
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

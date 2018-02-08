/**
 * \file communicators_pmmg.c
 * \brief Functions related to the communicators construction
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria Soft)
 * \author Nikos Pattakos (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "linkedlist_pmmg.h"

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * \return 1 if success, 0 if fail.
 *
 * Build the node communicators (externals and internals) from the faces ones.
 *
 */
int PMMG_build_nodeCommFromFaces( PMMG_pParMesh parmesh ) {

  /** Build the external node communicators from the faces ones without looking
   * at the possibility to miss a processor that is not visible from the face
   * communicators: _0_\1/_2_ configuration ( from 0 it is not possible to say
   * that the node at edge intersction belongs to 2 */
  if ( !PMMG_build_simpleExtNodeComm(parmesh) ) {
    fprintf(stderr,"\n  ## Error: %s: unable to build the simple externals node"
            " communicators from the external faces communicators.\n",__func__);
    return 0;
  }

  /** Build the internal node communicator from the faces ones */
  if ( !PMMG_build_intNodeComm(parmesh) ) {
    fprintf(stderr,"\n  ## Error: %s: unable to build the internal node"
            " communicators from the internal faces communicators.\n",__func__);
    return 0;
  }

  if ( !PMMG_build_completeExtNodeComm(parmesh) ) {
    fprintf(stderr,"\n  ## Error: %s: unable to complete the external node"
            " communicators.\n",__func__);
    return 0;
  }

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * \return 1 if success, 0 if fail.
 *
 * Build the external node communicators from the faces ones without looking
 * at the possibility to miss a processor that is not visible from the face
 * communicators: _0_\1/_2_ configuration ( from 0 it is not possible to say
 * that the node at edge intersction belongs to 2.
 *
 */
int PMMG_build_simpleExtNodeComm( PMMG_pParMesh parmesh ) {
  PMMG_pGrp       grp;
  PMMG_pext_comm  ext_node_comm,ext_face_comm;
  PMMG_pint_comm  int_face_comm;
  MMG5_pMesh      mesh;
  MMG5_pTetra     pt;
  MMG5_pPoint     ppt;
  int             *intvalues,*meshId,idx1,idx2;
  int             *face2int_face_comm_index1,*face2int_face_comm_index2;
  int             next_face_comm,next_node_comm,nitem_ext_comm,nitem_int_comm;
  int             color_out,ier,i,j,k,iel,ifac,ip,iploc,grpid;

  ier = 0;
  assert ( !parmesh->ext_node_comm && "external comms must be deleted" );

  next_face_comm = parmesh->next_face_comm;
  next_node_comm = next_face_comm;

  PMMG_CALLOC(parmesh,parmesh->ext_node_comm,next_node_comm,PMMG_ext_comm,
              "ext_node_comm ",return 0);
  parmesh->next_node_comm = next_node_comm;

  /** Store in the internal communicator the index of the interface faces and
   * initialize the tmp field of the group points */
  int_face_comm = parmesh->int_face_comm;
  PMMG_CALLOC(parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,
              "face communicator",goto end);
  intvalues     = int_face_comm->intvalues;

  PMMG_CALLOC(parmesh,meshId,int_face_comm->nitem,int,
              "id of the mesh associated to the face stored in intvalues",
              goto end);

  for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
    grp  = &parmesh->listgrp[grpid];
    mesh = grp->mesh;

    for ( k=1; k<=mesh->np; ++k )
      mesh->point[k].tmp = PMMG_UNSET;

    face2int_face_comm_index1 = grp->face2int_face_comm_index1;
    face2int_face_comm_index2 = grp->face2int_face_comm_index2;
    for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
      idx1 = face2int_face_comm_index1[i];
      idx2 = face2int_face_comm_index2[i];
      if ( intvalues[idx2] ) continue;

      intvalues[idx2] = idx1;
      meshId[idx2]    = grpid;
    }
  }

  nitem_int_comm = 0;
  for ( k=0; k<next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    color_out = ext_face_comm->color_out;

    parmesh->ext_node_comm[k].color_in  = ext_face_comm->color_in;
    parmesh->ext_node_comm[k].color_out = color_out;

    /** Initialisation of the vertices of the interface faces */
    for ( i=0; i<ext_face_comm->nitem; ++i ) {
      idx2  = ext_face_comm->int_comm_index[i];
      iel   =  intvalues[idx2]/12;
      ifac  = (intvalues[idx2]%12)/3;
      grpid = meshId[idx2];

      assert ( iel && iel<=mesh->ne );
      assert ( 0<=ifac && ifac<4 );

      mesh  = parmesh->listgrp[grpid].mesh;
      pt    = &mesh->tetra[iel];
      for ( j=0; j<3; ++j ) {
        ip  = pt->v[_MMG5_idir[ifac][j]];

        assert ( ip && ip<=mesh->np );

        ppt = &mesh->point[ip];
        ppt->flag = PMMG_UNSET;
      }
    }

    /** Worst case allocation of the node communicator */
    parmesh->ext_node_comm[k].nitem = 3*ext_face_comm->nitem;
    PMMG_CALLOC(parmesh,parmesh->ext_node_comm[k].int_comm_index,
                parmesh->ext_node_comm[k].nitem,int,
                "external node communicator",goto end);
    ext_node_comm = &parmesh->ext_node_comm[k];


    /** Process the external face communicator and fill the node one */
    nitem_ext_comm = 0;
    for ( i=0; i<ext_face_comm->nitem; ++i ) {
      idx2  = ext_face_comm->int_comm_index[i];
      iel   =  intvalues[idx2]/12;
      ifac  = (intvalues[idx2]%12)/3;
      iploc = (intvalues[idx2]%12)%3;
      grpid = meshId[idx2];

      assert ( iel && iel<=mesh->ne );
      assert ( 0<=ifac && ifac<4 );

      mesh  = parmesh->listgrp[grpid].mesh;
      pt    = &mesh->tetra[iel];

      for  ( j=0; j<3; ++j ) {

        /* Ensure that we travel the face in opposite direction on the 2
         * processors */
        if ( color_out > parmesh->ext_node_comm[k].color_in )
          ip  = pt->v[_MMG5_idir[ifac][(iploc+j)%3]];
        else
          ip  = pt->v[_MMG5_idir[ifac][(iploc+3-j)%3]];

        assert ( ip && ip<=mesh->np );

        ppt = &mesh->point[ip];

        if ( ppt->tmp < 0 ) {
          /* Give a position to this point in the internal communicator */
          ppt->tmp = nitem_int_comm++;
        }
        if ( ppt->flag < 0 ) {
          /* Add this point to the current external communicator */
          ppt->flag = 1;
          ext_node_comm->int_comm_index[nitem_ext_comm++] = ppt->tmp;
        }
      }
    }
    assert ( nitem_ext_comm <= ext_node_comm->nitem );
    PMMG_REALLOC(parmesh,parmesh->ext_node_comm[k].int_comm_index,nitem_ext_comm,
                 parmesh->ext_node_comm[k].nitem,int,
                 "external node communicator",goto end);
    ext_node_comm = &parmesh->ext_node_comm[k];
    ext_node_comm->nitem = nitem_ext_comm;
  }
  parmesh->int_node_comm->nitem = nitem_int_comm;

  /* Success */
  ier = 1;

end:
  if ( !ier ) {
    /* Deallocation of the external comms because we have failed */
    if ( next_node_comm ) {
      for ( k=0; k<next_node_comm; ++k ) {
        PMMG_DEL_MEM(parmesh,parmesh->ext_node_comm[k].int_comm_index,
                     ext_node_comm->nitem,int,"external node communicator");
      }
    }
    PMMG_DEL_MEM(parmesh,parmesh->ext_node_comm,parmesh->next_node_comm
                 ,PMMG_ext_comm, "ext_node_comm ");
  }

  PMMG_DEL_MEM(parmesh,parmesh->int_face_comm->intvalues,
               parmesh->int_face_comm->nitem,int,
               "internal face communicator");

  PMMG_DEL_MEM(parmesh,meshId,parmesh->int_face_comm->nitem,int,
               "id of the mesh associated to the face stored in intvalues");
  return ier;
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * \return 1 if success, 0 if fail.
 *
 * Build the internal node communicators from the faces ones.
 *
 */
int PMMG_build_intNodeComm( PMMG_pParMesh parmesh ) {
  PMMG_pGrp       grp;
  MMG5_pMesh      mesh;
  MMG5_pTetra     pt;
  MMG5_pPoint     ppt;
  int             idx1;
  int             *face2int_face_comm_index1,*face2int_face_comm_index2;
  int             *node2int_node_comm_index1,*node2int_node_comm_index2;
  int             nitem_node,nitem_int_node_comm,nitem_tmp;
  int             ier,i,j,iel,ifac,ip,grpid;

  ier = 0;

  nitem_node = parmesh->int_node_comm->nitem;

  for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
    /** Reset the flag field of the nodes of the group */
    grp  = &parmesh->listgrp[grpid];
    mesh = grp->mesh;

    face2int_face_comm_index1 = grp->face2int_face_comm_index1;
    face2int_face_comm_index2 = grp->face2int_face_comm_index2;

    /** Worst case allocation for the node2int_node arrays */
    grp->nitem_int_node_comm = 3*grp->nitem_int_face_comm;
    PMMG_CALLOC(parmesh,grp->node2int_node_comm_index1,grp->nitem_int_node_comm,
                int,"node2int_node_comm_index1",goto end);

    PMMG_CALLOC(parmesh,grp->node2int_node_comm_index2,grp->nitem_int_node_comm,
                int,"node2int_node_comm_index2",goto end);

    node2int_node_comm_index1 = grp->node2int_node_comm_index1;
    node2int_node_comm_index2 = grp->node2int_node_comm_index2;

    nitem_int_node_comm = 0;

    for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
      idx1 = face2int_face_comm_index1[i];
      iel  =  idx1/12;
      ifac = (idx1%12)/3;

      assert ( iel && iel<=mesh->ne );
      assert ( 0<=ifac && ifac<4 );

      pt    = &mesh->tetra[iel];
      for ( j=0; j<3; ++j ) {
        ip  = pt->v[_MMG5_idir[ifac][j]];

        assert ( ip && ip<=mesh->np );

        ppt = &mesh->point[ip];
        ppt->flag = PMMG_UNSET;
      }
    }

    for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
      idx1 = face2int_face_comm_index1[i];
      iel  =  idx1/12;
      ifac = (idx1%12)/3;

      assert ( iel && iel<=mesh->ne );
      assert ( 0<=ifac && ifac<4 );

      pt    = &mesh->tetra[iel];
      for ( j=0; j<3; ++j ) {
        ip  = pt->v[_MMG5_idir[ifac][j]];

        assert ( ip && ip<=mesh->np );

        ppt = &mesh->point[ip];
        if ( ppt->tmp < 0 ) {
          /** Give a position int he internal communicator to this point */
          ppt->tmp = nitem_node++;
        }
        else {
          if ( ppt->flag!=PMMG_UNSET ) continue;
          ppt->flag = 1;
        }
        /** Update the nod2int_node_comm arrays */
        node2int_node_comm_index1[nitem_int_node_comm] = ip;
        node2int_node_comm_index2[nitem_int_node_comm] = ppt->tmp;
        ++nitem_int_node_comm;
      }
    }

    assert ( nitem_int_node_comm <= grp->nitem_int_node_comm );

    nitem_tmp = grp->nitem_int_node_comm;
    PMMG_REALLOC(parmesh,grp->node2int_node_comm_index1,nitem_int_node_comm,
                 nitem_tmp,int,"node2int_node_comm_index1",
                 goto end);
    grp->nitem_int_node_comm = nitem_int_node_comm;

    PMMG_REALLOC(parmesh,grp->node2int_node_comm_index2,nitem_int_node_comm,
                 nitem_tmp,int,"node2int_node_comm_index2",
                 PMMG_DEL_MEM(parmesh,node2int_node_comm_index2,nitem_tmp,
                              int,"node2int_node_comm_index2");
                 goto end);
  }
  parmesh->int_node_comm->nitem = nitem_node;

  /* Success */
  ier = 1;

end:
  if ( !ier ) {
    /* Dealloc of the communicators because we have failed */
    for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
      grp  = &parmesh->listgrp[grpid];

      PMMG_DEL_MEM(parmesh,grp->node2int_node_comm_index1,
                   grp->nitem_int_node_comm,int,"node2int_node_comm_index1");
      PMMG_DEL_MEM(parmesh,grp->node2int_node_comm_index2,
                   grp->nitem_int_node_comm,int,"node2int_node_comm_index2");
    }
  }

  return ier;
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * \return 1 if success, 0 if fail.
 *
 * Complete the external communicators by travelling through the processors and
 * faces to detect all the processors to which each node belongs.
 *
 */
int PMMG_build_completeExtNodeComm( PMMG_pParMesh parmesh ) {
  PMMG_pext_comm  ext_node_comm,*comm_ptr;
  PMMG_pint_comm  int_node_comm;
  PMMG_lnkdList   **proclists,list;
  int             *intvalues,nitem,nproclists,ier,ier2,k,i,j,idx,pos,rank,color;
  int             *itosend,*itorecv,*i2send_size,*i2recv_size,nitem2comm;
  int             *nitem_ext_comm,next_comm,val1_i,val2_i,val1_j,val2_j;;
  int8_t          glob_update,loc_update;
  MPI_Request     *request;
  MPI_Status      *status;

  ier = 0;

  rank          = parmesh->myrank;
  int_node_comm = parmesh->int_node_comm;
  nitem         = int_node_comm->nitem;

  PMMG_CALLOC(parmesh,int_node_comm->intvalues,nitem,int,"node communicator",
    return 0);
  intvalues     = int_node_comm->intvalues;

  /** Use intvalues to flag the already treated points of the communicator:
   * initialization to 0.  */
  for ( k=0; k<nitem; ++k ) intvalues[k] = 0;

  proclists = NULL;
  PMMG_CALLOC(parmesh,proclists,nitem,PMMG_lnkdList*,"array of linked lists",
              goto end);

  /* Reallocation of the list of external comms at maximal size (nprocs) to
   * avoid tricky treatment when filling it.*/
  PMMG_REALLOC(parmesh,parmesh->ext_node_comm,parmesh->nprocs,
               parmesh->next_node_comm,PMMG_ext_comm,
               "list of external communicators",goto end);
  next_comm = parmesh->next_node_comm;
  parmesh->next_node_comm = parmesh->nprocs;

  for ( k=next_comm; k<parmesh->nprocs; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    ext_node_comm->nitem          = 0;
    ext_node_comm->int_comm_index = NULL;
    ext_node_comm->color_in       = rank;
    ext_node_comm->color_out      = -1;
  }

  comm_ptr = NULL;
  PMMG_CALLOC(parmesh,comm_ptr,parmesh->nprocs,PMMG_pext_comm,
              "array of pointers toward the external communicators",
              goto end);

  /** Step 1: initialization of the list of the procs to which a point belongs
   * by the value of the current mpi rank */
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];

    if ( !ext_node_comm->nitem ) continue;

    comm_ptr[ext_node_comm->color_out] = ext_node_comm;

    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx = ext_node_comm->int_comm_index[i];

      if ( intvalues[idx] ) continue;

      PMMG_CALLOC(parmesh,proclists[idx],1,PMMG_lnkdList,"linked list pointer",
                  goto end);
      if ( !PMMG_lnkdListNew(parmesh,proclists[idx],PMMG_LISTSIZE) ) goto end;

      if ( !PMMG_add_cell2lnkdList( parmesh,proclists[idx],rank,idx ) )
        goto end;

      intvalues[idx] = 1;
    }
  }
  /* Fill the missing pointer toward the empty external communicators */
  for ( k=0; k<parmesh->nprocs; ++k ) {
    if ( !comm_ptr[k] ) {
      comm_ptr[k] = &parmesh->ext_node_comm[next_comm++];
      comm_ptr[k]->color_out = k;
    }
  }

  /** Step 2: While at least the proc list of 1 node is modified, send and
   * recieve the proc list of all the nodes to/from the other processors. At the
   * end of this loop, each node has the entire list of the proc to which it
   * belongs */
  PMMG_MALLOC(parmesh,request,parmesh->next_node_comm,MPI_Request,
              "mpi request array",goto end);

  PMMG_MALLOC(parmesh,status,parmesh->next_node_comm,MPI_Status,
              "mpi status array",goto end);

  PMMG_CALLOC(parmesh,i2send_size,parmesh->next_node_comm,int,
              "size of the i2send array",goto end);
  PMMG_CALLOC(parmesh,i2recv_size,parmesh->next_node_comm,int,
              "size of the i2recv array",goto end);

  if ( !PMMG_lnkdListNew(parmesh,&list,PMMG_LISTSIZE) ) goto end;

  do {
    glob_update = loc_update = 0;

    /** Send the list of procs to which belong each point of the communicator */
    for ( k=0; k<parmesh->next_node_comm; ++k ) {

      request[k] = MPI_REQUEST_NULL;
      ext_node_comm = &parmesh->ext_node_comm[k];

      /* Computation of the number of data to send to the other procs (we want
       * to send the the size of the linked list and the val1 and val2 fields of
       * each cell ) */
      nitem2comm = 0;
      for ( i=0; i<ext_node_comm->nitem; ++i ) {
        idx         = ext_node_comm->int_comm_index[i];
        nitem2comm += proclists[idx]->nitem*2+1;
      }

      if ( i2send_size[k] < nitem2comm ) {
        PMMG_REALLOC(parmesh,ext_node_comm->itosend,nitem2comm,i2send_size[k],int,
                     "itosend",goto end);
        i2send_size[k] = nitem2comm;
      }

      /* Filling of the array to send */
      pos     = 0;
      itosend = ext_node_comm->itosend;
      color   = ext_node_comm->color_out;
      for ( i=0; i<ext_node_comm->nitem; ++i ) {
        idx  = ext_node_comm->int_comm_index[i];
        pos += PMMG_packInArray_lnkdList(proclists[idx],&itosend[pos]);
      }
      MPI_CHECK( MPI_Isend(itosend,nitem2comm,MPI_INT,color,0,parmesh->comm,
                           &request[color]),goto end );
    }

    /** Recv the list of procs to which belong each point of the communicator */
    for ( k=0; k<parmesh->next_node_comm; ++k ) {
      ext_node_comm = &parmesh->ext_node_comm[k];

      if ( !ext_node_comm->nitem ) continue;

      color         = ext_node_comm->color_out;

      MPI_CHECK( MPI_Probe(color,0,parmesh->comm,&status[0] ),goto end);
      MPI_CHECK( MPI_Get_count(&status[0],MPI_INT,&nitem2comm),goto end);

      if ( i2recv_size[k] < nitem2comm ) {
        PMMG_REALLOC(parmesh,ext_node_comm->itorecv,nitem2comm,i2recv_size[k],
                     int,"itorecv",goto end);
        i2recv_size[k] = nitem2comm;
      }
      itorecv       = ext_node_comm->itorecv;
      MPI_CHECK( MPI_Recv(itorecv,nitem2comm,MPI_INT,color,0,parmesh->comm,
                          &status[0]), goto end );

      pos     = 0;
      for ( i=0; i<ext_node_comm->nitem; ++i ) {
        idx  = ext_node_comm->int_comm_index[i];
        assert ( idx>=0 );

        PMMG_reset_lnkdList( parmesh,&list );
        ier2 = PMMG_unpackArray_inLnkdList(parmesh,&list,&itorecv[pos]);

        if ( ier2 < 0 ) goto end;
        pos += ier2;

        ier2 = PMMG_merge_lnkdList(parmesh,proclists[idx],&list);
        if ( !ier2 ) goto end;

        loc_update |= (ier2%2);
      }
    }

    MPI_CHECK( MPI_Waitall(parmesh->next_node_comm,request,status), goto end );
    MPI_CHECK( MPI_Allreduce(&loc_update,&glob_update,1,MPI_INT8_T,MPI_LOR,
                             parmesh->comm),goto end);

  } while ( glob_update );

  /** Step 3: Cancel the old external communicator and build it again from the
   * list of proc of each node */
  PMMG_CALLOC(parmesh,nitem_ext_comm,parmesh->nprocs,int,
              "number of items in each external communicator",goto end);

  /* Remove the empty proc lists */
  nproclists = 0;
  for ( i=0; i<nitem; ++i ) {
    if ( intvalues[i] )
      proclists[nproclists++] = proclists[i];
  }

  /* Sort the list of procs to which each node belong to ensure that we will
   * build the external communicators in the same order on both processors
   * involved in an external comm */
  qsort(proclists,nproclists,sizeof(PMMG_lnkdList*),PMMG_compare_lnkdList);

  /* Remove the non unique paths */
  if ( nproclists ) {
    idx = 0;
    for ( i=1; i<nproclists; ++i ) {
      assert ( proclists[i]->nitem );
      if ( PMMG_compare_lnkdList(&proclists[i],&proclists[idx]) ) {
        ++idx;
        if ( idx != i ) {
          proclists[idx] = proclists[i];
        }
      }
    }
    nproclists = idx+1;
  }

  /* Double loop over the list of procs to which belong each point to add the
   * couples to the suitable external communicators */
  for ( k=0; k<nproclists; ++k ) {
    for ( i=0; i<proclists[k]->nitem; ++i ) {
      val1_i = proclists[k]->item[i].val1;
      val2_i = proclists[k]->item[i].val2;

      for ( j=i+1; j<proclists[k]->nitem; ++j ) {
        val1_j = proclists[k]->item[j].val1;
        val2_j = proclists[k]->item[j].val2;

        assert ( val1_i != val1_j );

        if ( val1_i == rank ) {
          ext_node_comm = comm_ptr[val1_j];
          assert ( ext_node_comm );

          if ( !ext_node_comm->nitem ) continue;

          if ( nitem_ext_comm[val1_j] == ext_node_comm->nitem ) {
            /* Reallocation */
            PMMG_REALLOC(parmesh,ext_node_comm->int_comm_index,
                         (int)((1.+PMMG_GAP)*ext_node_comm->nitem)+1,
                         ext_node_comm->nitem,int,
                         "external communicator",goto end);
            ext_node_comm->nitem = (int)((1.+PMMG_GAP)*ext_node_comm->nitem)+1;
            comm_ptr[val1_j] = ext_node_comm;
          }
          ext_node_comm->int_comm_index[nitem_ext_comm[val1_j]++] = val2_i;
        }
        else if ( val1_j == rank ) {
          ext_node_comm = comm_ptr[val1_i];
          assert ( ext_node_comm );
          if ( nitem_ext_comm[val1_j] == ext_node_comm->nitem ) {
            /* Reallocation */
            PMMG_REALLOC(parmesh,ext_node_comm->int_comm_index,
                          (int)((1.+PMMG_GAP)*ext_node_comm->nitem)+1,
                         ext_node_comm->nitem,int,
                         "external communicator",goto end);
            ext_node_comm->nitem = (int)((1.+PMMG_GAP)*ext_node_comm->nitem)+1;
            comm_ptr[val1_i] = ext_node_comm;
          }
          ext_node_comm->int_comm_index[nitem_ext_comm[val1_i]++] = val2_j;
        }
      }
    }
  }

  /* Pack and reallocate the external communicators */
  // The pack may be removed if we allow to have always nprocs external comms
  next_comm = 0;
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    if ( !ext_node_comm->nitem ) {
      /* Empty communicator */
      continue;
    }
    else {
      assert ( ext_node_comm->color_out>=0 && ext_node_comm->color_out!=rank );
      assert ( nitem_ext_comm[ext_node_comm->color_out] );

      PMMG_REALLOC(parmesh,ext_node_comm->int_comm_index,
                   nitem_ext_comm[ext_node_comm->color_out],
                   ext_node_comm->nitem,int,
                   "external communicator",goto end);
      ext_node_comm->nitem = nitem_ext_comm[ext_node_comm->color_out];

      if ( next_comm != k )
        parmesh->ext_node_comm[next_comm] = *ext_node_comm;
    }
    ++next_comm;
  }
  PMMG_REALLOC(parmesh,parmesh->ext_node_comm,
               next_comm,parmesh->next_node_comm,PMMG_ext_comm,
               "list of external communicator",goto end);
  parmesh->next_node_comm = next_comm;

  /* Success */
  ier = 1;

end:
  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,nitem,int,"node communicator");

  if ( proclists ) {
    for ( k=0; k<nitem; ++k ) {
      if ( proclists[k] ) {
        PMMG_DEL_MEM(parmesh,proclists[k]->item,proclists[k]->nitem_max,
                     PMMG_lnkdCell,"linked list array");

      }
    }
    PMMG_DEL_MEM(parmesh,proclists,nitem,PMMG_lnkdList*,"array of linked lists");
  }
  PMMG_DEL_MEM(parmesh,comm_ptr,parmesh->nprocs,PMMG_pext_comm,
              "array of pointers toward the external communicators");

  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];

    // Change this and add to the external comm the possibility to not
    // unalloc/realloc every time, thus, here, we will be able to reset the
    // communicators without unallocated it
    PMMG_DEL_MEM(parmesh,ext_node_comm->itosend,i2send_size[k],int,"i2send");
    PMMG_DEL_MEM(parmesh,ext_node_comm->itorecv,i2recv_size[k],int,"i2recv");
  }

  PMMG_DEL_MEM(parmesh,request,parmesh->next_node_comm,MPI_Request,
               "mpi request array");
  PMMG_DEL_MEM(parmesh,status,parmesh->next_node_comm,MPI_Status,
               "mpi status array");

  PMMG_DEL_MEM(parmesh,i2send_size,parmesh->next_node_comm,int,
               "size of the i2send array");

  PMMG_DEL_MEM(parmesh,i2recv_size,parmesh->next_node_comm,int,
               "size of the i2recv array");

  PMMG_DEL_MEM(parmesh,nitem_ext_comm,parmesh->nprocs,int,
               "number of items in each external communicator");

  PMMG_DEL_MEM(parmesh,list.item,list.nitem_max,PMMG_lnkdCell,
               "linked list array");

  return ier;
}

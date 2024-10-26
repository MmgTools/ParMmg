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
 * \file communicators_pmmg.c
 * \brief Functions related to the communicators construction
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria Soft)
 * \author Nikos Pattakos (Inria)
 * \author Luca Cirrottola (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "linkedlist_pmmg.h"
#include "coorcell_pmmg.h"

/**
 * \param parmesh pointer toward a parmesh structure
 * \param comm    external communicator to be freed
 *
 * deallocate all internal communicator's fields
 */
void PMMG_parmesh_int_comm_free( PMMG_pParMesh parmesh,PMMG_pInt_comm comm )
{
  if ( comm == NULL ) {
    return;
  }

  if ( NULL != comm->intvalues ) {
    assert ( comm->nitem != 0 && "incorrect parameters in internal communicator" );
    PMMG_DEL_MEM(parmesh,comm->intvalues,int,"int comm int array");
  }
  if ( NULL != comm->doublevalues ) {
    assert ( comm->nitem != 0 && "incorrect parameters in internal communicator" );
    PMMG_DEL_MEM(parmesh,
                 comm->doublevalues,double,"int comm double array");
  }
}

/**
 * \param parmesh  pointer toward a parmesh structure
 * \param listcomm external communicators to be freed
 * \param ncomm    parameter ncomm
 *
 * deallocate all external communicators's fields
 */
void PMMG_parmesh_ext_comm_free( PMMG_pParMesh parmesh,PMMG_pExt_comm listcomm,
                                 int ncomm )
{
  PMMG_pExt_comm comm;
  int i = 0;

  if ( listcomm == NULL )
    return;

  for( i = 0; i < ncomm; ++i ) {
    comm = &listcomm[i];
    if ( NULL != comm->int_comm_index ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(parmesh,comm->int_comm_index,int,"ext comm int array");
    }
    if ( NULL != comm->itosend ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(parmesh,comm->itosend,int,"ext comm itosend array");
    }
    if ( NULL != comm->itorecv ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(parmesh,comm->itorecv,int,"ext comm itorecv array");
    }
    if ( NULL != comm->rtosend ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(parmesh,comm->rtosend,double,"ext comm rtosend array");
    }
    if ( NULL != comm->rtorecv ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(parmesh,comm->rtorecv,double,"ext comm rtorecv array");
    }
  }
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param idx1    node2int_node_comm_index1 to be freed
 * \param idx2    node2int_node_comm_index2 to be freed
 * \param n       pointer to node2int_node_comm_nitem size
 *
 * Deallocate all the MMG3D meshes and their communicators and zero the size
 */
void PMMG_grp_comm_free( PMMG_pParMesh parmesh,int **idx1,int **idx2,
                                 int *n )
{
  PMMG_DEL_MEM(parmesh,*idx1,int,"group communicator");
  PMMG_DEL_MEM(parmesh,*idx2,int,"group communicator");
  *n = 0;
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * Deallocate the edge communicators of the parmesh
 *
 */
void PMMG_edge_comm_free( PMMG_pParMesh parmesh )
{
  PMMG_pGrp grp;
  int       k;

  assert( parmesh->ngrp == 1 );

  grp = &parmesh->listgrp[0];
  PMMG_grp_comm_free( parmesh,
                      &grp->edge2int_edge_comm_index1,
                      &grp->edge2int_edge_comm_index2,
                      &grp->nitem_int_edge_comm );

  PMMG_parmesh_int_comm_free( parmesh,parmesh->int_edge_comm );
  PMMG_parmesh_ext_comm_free( parmesh,parmesh->ext_edge_comm,parmesh->next_edge_comm);
  PMMG_DEL_MEM(parmesh, parmesh->ext_edge_comm,PMMG_Ext_comm,"ext edge comm");

  parmesh->next_edge_comm       = 0;

  if ( parmesh->int_edge_comm ) {
    parmesh->int_edge_comm->nitem = 0;
  }
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * Deallocate the nodal communicators of the parmesh
 *
 */
void PMMG_node_comm_free( PMMG_pParMesh parmesh )
{
  PMMG_pGrp grp;
  int       k;

  for ( k=0; k<parmesh->ngrp;++k ) {
    grp = &parmesh->listgrp[k];
    PMMG_grp_comm_free( parmesh,&grp->node2int_node_comm_index1,
                        &grp->node2int_node_comm_index2,
                        &grp->nitem_int_node_comm );
  }

  PMMG_parmesh_int_comm_free( parmesh,parmesh->int_node_comm);
  PMMG_parmesh_ext_comm_free( parmesh,parmesh->ext_node_comm,parmesh->next_node_comm);
  PMMG_DEL_MEM(parmesh, parmesh->ext_node_comm,PMMG_Ext_comm,"ext node comm");

  parmesh->next_node_comm       = 0;
  if ( parmesh->int_node_comm ) {
    parmesh->int_node_comm->nitem = 0;
  }
}

/**
 * \param parmesh pointer to parmesh structure
 * \param mesh pointer to the mesh structure
 * \param hpar hash table of parallel edges
 *
 * Build internal edge communicator.
 */
int PMMG_build_intEdgeComm( PMMG_pParMesh parmesh,MMG5_pMesh mesh,MMG5_HGeom *hpar ) {
  PMMG_pGrp      grp;
  PMMG_pInt_comm int_edge_comm;
  MMG5_pEdge     pa;
  MMG5_hgeom     *ph;
  int            k;

  assert( parmesh->ngrp == 1 );
  grp = &parmesh->listgrp[0];

  PMMG_CALLOC(parmesh,parmesh->int_edge_comm,1,PMMG_Int_comm,"int_edge_comm",return 0);
  int_edge_comm = parmesh->int_edge_comm;

  /* here edges should not be allocated, unless due to some debut I/O */
  if( mesh->na ){
    MMG5_DEL_MEM(mesh,mesh->edge);
    mesh->na = 0;
  }

  /* Count edges (the hash table only contains parallel edges) */
  for( k = 0; k <= hpar->max; k++ ) {
    ph = &hpar->geom[k];
    if( !(ph->a) ) continue;
    mesh->na++;
  }

  /* Create edge array */
  if ( mesh->na ) {
    MMG5_ADD_MEM(mesh,(mesh->na+1)*sizeof(MMG5_Edge),"edges",
                 return 0;
                 printf("  ## Warning: uncomplete mesh\n"));
    MMG5_SAFE_CALLOC(mesh->edge,mesh->na+1,MMG5_Edge,return 0);

    mesh->na = 0;
    for( k = 0; k <= hpar->max; k++ ) {
      ph = &hpar->geom[k];
      if ( !ph->a )  continue;
      /* Get edge */
      mesh->na++;
      mesh->edge[mesh->na].a    = ph->a;
      mesh->edge[mesh->na].b    = ph->b;
      /* Set links between hash table and array */
      ph->ref = mesh->na;
      mesh->edge[mesh->na].ref = k;
      /* Use base to keep track of the out rank */
      mesh->edge[mesh->na].base = PMMG_UNSET;
    }
  }

  /* Set nb. of items */
  int_edge_comm->nitem = mesh->na;
  grp->nitem_int_edge_comm = mesh->na;

  /* Set group indices to the edge array and the internal communicator */
  PMMG_MALLOC(parmesh,grp->edge2int_edge_comm_index1,grp->nitem_int_edge_comm,int,"edge2int_edge_comm_index1",return 0);
  PMMG_MALLOC(parmesh,grp->edge2int_edge_comm_index2,grp->nitem_int_edge_comm,int,"edge2int_edge_comm_index2",return 0);
  for( k = 0; k < grp->nitem_int_edge_comm; k++ ) {
    grp->edge2int_edge_comm_index1[k] = k+1;
    grp->edge2int_edge_comm_index2[k] = k;
  }

  return 1;
}

/**
 * \param parmesh pointer to parmesh structure
 * \param mesh pointer to the mesh structure
 * \param hpar hash table of parallel edges
 * \param ext_edge_comm pointer to the external edge communicator
 * \param pt pointer to the tetra
 * \param ifac face index
 * \param iloc local face vertex
 * \param j local edge to process
 * \param color used to mark already processed edges
 * \param pointer to the index of the next free position in the external comm
 * \return 0 if fail, 1 if success.
 *
 * Fill an item of the external edge communicator from a parallel face edge.
 */
int PMMG_fillExtEdgeComm_fromFace( PMMG_pParMesh parmesh,MMG5_pMesh mesh,MMG5_HGeom *hpar,
                                   PMMG_pExt_comm ext_edge_comm,MMG5_pTetra pt,int ifac,int iloc,int j,int color,int *item ) {
  MMG5_pEdge pa;
  int        edg;
  uint16_t   tag;
  int8_t     i1,i2;

  /* Take the edge opposite to vertex iloc+j on face ifac */
  i1 = MMG5_idir[ifac][(iloc+j+1)%3];
  i2 = MMG5_idir[ifac][(iloc+j+2)%3];
  if ( !MMG5_hGet( hpar, pt->v[i1], pt->v[i2], &edg, &tag ) ) return 0;
  pa = &mesh->edge[edg];
  /* Fill item and overwrite edge base with current color */
  if( pa->base != color ) {
    /* the position of the edge in the internal communicator is simply its
     * index-1 */
    ext_edge_comm->int_comm_index[(*item)++] = edg-1;
    pa->base = color;
  }
  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param comm pointer toward the MPI communicator to use: when called before
 * the first mesh balancing (at preprocessing stage) we have to use the
 * read_comm communicator (i.e. the communicator used to provide the inputs).
 * For all ather calls, comm has to be the communicator to use for computations.
 *
 * \return 1 if success, 0 if fail.
 *
 * Complete the external communicators by travelling through the processors and
 * faces to detect all the processors to which each edge belongs.
 *
 */
int PMMG_build_completeExtEdgeComm( PMMG_pParMesh parmesh, MPI_Comm comm ) {
  PMMG_pExt_comm    ext_edge_comm,*comm_ptr;
  PMMG_pInt_comm    int_edge_comm;
  PMMG_cellLnkdList **proclists,list;
  int               *intvalues,nitem,nproclists,ier,ier2,k,i,j,idx,pos,rank,color;
  int               *itosend,*itorecv,*i2send_size,*i2recv_size,nitem2comm;
  int               *nitem_ext_comm,next_comm,val1_i,val2_i,val1_j,val2_j;
  int               alloc_size;
  int8_t            glob_update,loc_update;
  MPI_Request       *request;
  MPI_Status        *status;

  ier = 0;

  rank          = parmesh->myrank;
  int_edge_comm = parmesh->int_edge_comm;
  nitem         = int_edge_comm->nitem;

  proclists       = NULL;
  comm_ptr        = NULL;
  request         = NULL;
  status          = NULL;
  i2send_size     = NULL;
  i2recv_size     = NULL;
  nitem_ext_comm  = NULL;
  list.item       = NULL;

  PMMG_CALLOC(parmesh,int_edge_comm->intvalues,nitem,int,"edge communicator",
    return 0);
  intvalues     = int_edge_comm->intvalues;

  /** Use intvalues to flag the already treated points of the communicator:
   * initialization to 0.  */
  for ( k=0; k<nitem; ++k ) intvalues[k] = 0;

  PMMG_CALLOC(parmesh,proclists,nitem,PMMG_cellLnkdList*,"array of linked lists",
              goto end);

  /* Reallocation of the list of external comms at maximal size (nprocs) to
   * avoid tricky treatment when filling it.*/
  PMMG_REALLOC(parmesh,parmesh->ext_edge_comm,parmesh->nprocs,
               parmesh->next_edge_comm,PMMG_Ext_comm,
               "list of external communicators",goto end);
  next_comm = parmesh->next_edge_comm;
  parmesh->next_edge_comm = parmesh->nprocs;

  for ( k=next_comm; k<parmesh->nprocs; ++k ) {
    ext_edge_comm = &parmesh->ext_edge_comm[k];
    ext_edge_comm->nitem          = 0;
    ext_edge_comm->nitem_to_share = 0;
    ext_edge_comm->int_comm_index = NULL;
    ext_edge_comm->itosend        = NULL;
    ext_edge_comm->itorecv        = NULL;
    ext_edge_comm->rtosend        = NULL;
    ext_edge_comm->rtorecv        = NULL;
    ext_edge_comm->color_in       = rank;
    ext_edge_comm->color_out      = PMMG_UNSET;
  }

  PMMG_CALLOC(parmesh,comm_ptr,parmesh->nprocs,PMMG_pExt_comm,
              "array of pointers toward the external communicators",
              goto end);

  /** Step 1: initialization of the list of the procs to which an edge belongs
   * by the value of the current mpi rank */
  for ( k=0; k<parmesh->next_edge_comm; ++k ) {
    ext_edge_comm = &parmesh->ext_edge_comm[k];

    if ( !ext_edge_comm->nitem ) continue;

    comm_ptr[ext_edge_comm->color_out] = ext_edge_comm;

    for ( i=0; i<ext_edge_comm->nitem; ++i ) {
      idx = ext_edge_comm->int_comm_index[i];

      if ( intvalues[idx] ) continue;

      PMMG_CALLOC(parmesh,proclists[idx],1,PMMG_cellLnkdList,"linked list pointer",
                  goto end);
      if ( !PMMG_cellLnkdListNew(parmesh,proclists[idx],idx,PMMG_LISTSIZE) ) goto end;

      if ( !PMMG_add_cell2lnkdList( parmesh,proclists[idx],rank,idx ) )
        goto end;

      intvalues[idx] = 1;
    }
  }
  /* Fill the missing pointer toward the empty external communicators */
  for ( k=0; k<parmesh->nprocs; ++k ) {
    if ( !comm_ptr[k] ) {

      assert ( next_comm < parmesh->nprocs );

      comm_ptr[k] = &parmesh->ext_edge_comm[next_comm++];
      comm_ptr[k]->color_out = k;
    }
  }

  /** Step 2: While at least the proc list of 1 edge is modified, send and
   * recieve the proc list of all the edges to/from the other processors. At the
   * end of this loop, each edge has the entire list of the proc to which it
   * belongs */
  alloc_size = parmesh->nprocs;
  PMMG_MALLOC(parmesh,request,    alloc_size,MPI_Request,"mpi request array",goto end);
  PMMG_MALLOC(parmesh,status,     alloc_size,MPI_Status,"mpi status array",goto end);
  PMMG_CALLOC(parmesh,i2send_size,alloc_size,int,"size of the i2send array",goto end);
  PMMG_CALLOC(parmesh,i2recv_size,alloc_size,int,"size of the i2recv array",goto end);

  if ( !PMMG_cellLnkdListNew(parmesh,&list,0,PMMG_LISTSIZE) ) goto end;

  do {
    glob_update = loc_update = 0;

    /** Send the list of procs to which belong each point of the communicator */
    for ( k=0; k<alloc_size; ++k ) {

      request[k] = MPI_REQUEST_NULL;
    }

    for ( k=0; k<parmesh->next_edge_comm; ++k ) {
      ext_edge_comm = &parmesh->ext_edge_comm[k];

      /* Computation of the number of data to send to the other procs (we want
       * to send the the size of the linked list and the val1 and val2 fields of
       * each cell ) */
      nitem2comm = 0;
      if ( !ext_edge_comm->nitem ) continue;

      for ( i=0; i<ext_edge_comm->nitem; ++i ) {
        idx         = ext_edge_comm->int_comm_index[i];
        nitem2comm += proclists[idx]->nitem*2+1;
      }

      if ( i2send_size[k] < nitem2comm ) {
        PMMG_REALLOC(parmesh,ext_edge_comm->itosend,nitem2comm,i2send_size[k],int,
                     "itosend",goto end);
        i2send_size[k] = nitem2comm;
      }

      /* Filling of the array to send */
      pos     = 0;
      itosend = ext_edge_comm->itosend;
      color   = ext_edge_comm->color_out;
      for ( i=0; i<ext_edge_comm->nitem; ++i ) {
        idx  = ext_edge_comm->int_comm_index[i];
        pos += PMMG_packInArray_cellLnkdList(proclists[idx],&itosend[pos]);
      }
      assert ( pos==nitem2comm );

      MPI_CHECK( MPI_Isend(itosend,nitem2comm,MPI_INT,color,
                           MPI_COMMUNICATORS_EDGE_TAG,comm,
                           &request[color]),goto end );
    }

    /** Recv the list of procs to which belong each point of the communicator */
    for ( k=0; k<parmesh->next_edge_comm; ++k ) {
      ext_edge_comm = &parmesh->ext_edge_comm[k];

      if ( !ext_edge_comm->nitem ) continue;

      color         = ext_edge_comm->color_out;

      MPI_CHECK( MPI_Probe(color,MPI_COMMUNICATORS_EDGE_TAG,comm,
                           &status[0] ),goto end);
      MPI_CHECK( MPI_Get_count(&status[0],MPI_INT,&nitem2comm),goto end);

      if ( i2recv_size[k] < nitem2comm ) {
        PMMG_REALLOC(parmesh,ext_edge_comm->itorecv,nitem2comm,i2recv_size[k],
                     int,"itorecv",goto end);
        i2recv_size[k] = nitem2comm;
      }

      if ( nitem2comm ) {

        itorecv       = ext_edge_comm->itorecv;
        MPI_CHECK( MPI_Recv(itorecv,nitem2comm,MPI_INT,color,
                            MPI_COMMUNICATORS_EDGE_TAG,comm,
                            &status[0]), goto end );

        pos     = 0;
        for ( i=0; i<ext_edge_comm->nitem; ++i ) {
          idx  = ext_edge_comm->int_comm_index[i];
          assert ( idx>=0 );

          PMMG_reset_cellLnkdList( parmesh,&list );
          ier2 = PMMG_unpackArray_inCellLnkdList(parmesh,&list,&itorecv[pos]);

          if ( ier2 < 0 ) goto end;
          pos += ier2;

          ier2 = PMMG_merge_cellLnkdList(parmesh,proclists[idx],&list);
          if ( !ier2 ) goto end;

          loc_update |= (ier2%2);
        }
      }
    }

    MPI_CHECK( MPI_Waitall(alloc_size,request,status), goto end );
    MPI_CHECK( MPI_Allreduce(&loc_update,&glob_update,1,MPI_INT8_T,MPI_LOR,
                             comm),goto end);

  } while ( glob_update );

  /** Step 3: Cancel the old external communicator and build it again from the
   * list of proc of each edge */
  PMMG_CALLOC(parmesh,nitem_ext_comm,parmesh->nprocs,int,
              "number of items in each external communicator",goto end);

  /* Remove the empty proc lists */
  nproclists = 0;
  for ( i=0; i<nitem; ++i ) {
    if ( intvalues[i] )
      proclists[nproclists++] = proclists[i];
  }

  /* Sort the list of procs to which each edge belong to ensure that we will
   * build the external communicators in the same order on both processors
   * involved in an external comm */
  qsort(proclists,nproclists,sizeof(PMMG_cellLnkdList*),PMMG_compare_cellLnkdList);

  /* Remove the non unique paths */
  if ( nproclists ) {
    idx = 0;
    for ( i=1; i<nproclists; ++i ) {
      assert ( proclists[i]->nitem );
      if ( PMMG_compare_cellLnkdList(&proclists[i],&proclists[idx]) ) {
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
          ext_edge_comm = comm_ptr[val1_j];
          assert ( ext_edge_comm );
          if ( nitem_ext_comm[val1_j] == ext_edge_comm->nitem || !ext_edge_comm->nitem ) {
            /* Reallocation */
            PMMG_REALLOC(parmesh,ext_edge_comm->int_comm_index,
                         (int)((1.+PMMG_GAP)*ext_edge_comm->nitem)+1,
                         ext_edge_comm->nitem,int,
                         "external communicator",goto end);
            ext_edge_comm->nitem = (int)((1.+PMMG_GAP)*ext_edge_comm->nitem)+1;
            comm_ptr[val1_j] = ext_edge_comm;
          }
          ext_edge_comm->int_comm_index[nitem_ext_comm[val1_j]++] = val2_i;
        }
        else if ( val1_j == rank ) {
          ext_edge_comm = comm_ptr[val1_i];
          assert ( ext_edge_comm );
          if ( nitem_ext_comm[val1_i] == ext_edge_comm->nitem || !ext_edge_comm->nitem ) {
            /* Reallocation */
            PMMG_REALLOC(parmesh,ext_edge_comm->int_comm_index,
                          (int)((1.+PMMG_GAP)*ext_edge_comm->nitem)+1,
                         ext_edge_comm->nitem,int,
                         "external communicator",goto end);
            ext_edge_comm->nitem = (int)((1.+PMMG_GAP)*ext_edge_comm->nitem)+1;
            comm_ptr[val1_i] = ext_edge_comm;
          }
          ext_edge_comm->int_comm_index[nitem_ext_comm[val1_i]++] = val2_j;
        }
      }
    }
  }

  /* Pack and reallocate the external communicators */
  // The pack may be removed if we allow to have always nprocs external comms
  next_comm = 0;
  for ( k=0; k<parmesh->next_edge_comm; ++k ) {
    ext_edge_comm = &parmesh->ext_edge_comm[k];
    if ( !ext_edge_comm->nitem ) {
      /* Empty communicator */
      continue;
    }
    else {
      assert ( ext_edge_comm->color_out>=0 && ext_edge_comm->color_out!=rank );

      PMMG_REALLOC(parmesh,ext_edge_comm->int_comm_index,
                   nitem_ext_comm[ext_edge_comm->color_out],
                   ext_edge_comm->nitem,int,
                   "external communicator",goto end);
      ext_edge_comm->nitem = nitem_ext_comm[ext_edge_comm->color_out];

      if ( next_comm != k )
        parmesh->ext_edge_comm[next_comm] = *ext_edge_comm;
    }
    ++next_comm;
  }
  PMMG_REALLOC(parmesh,parmesh->ext_edge_comm,
               next_comm,parmesh->next_edge_comm,PMMG_Ext_comm,
               "list of external communicator",goto end);
  parmesh->next_edge_comm = next_comm;

  /* Success */
  ier = 1;

end:
  PMMG_DEL_MEM(parmesh,int_edge_comm->intvalues,int,"edge communicator");

  if ( proclists ) {
    for ( k=0; k<nproclists; ++k ) {
      PMMG_DEL_MEM(parmesh,proclists[k]->item,PMMG_lnkdCell,"linked list array");
      PMMG_DEL_MEM(parmesh,proclists[k],PMMG_lnkdList,"linked list pointer");
    }
    PMMG_DEL_MEM(parmesh,proclists,PMMG_lnkdList*,"array of linked lists");
  }
  PMMG_DEL_MEM(parmesh,comm_ptr,PMMG_pExt_comm,
              "array of pointers toward the external communicators");

  for ( k=0; k<parmesh->next_edge_comm; ++k ) {
    ext_edge_comm = &parmesh->ext_edge_comm[k];

    // Change this and add to the external comm the possibility to not
    // unalloc/realloc every time, thus, here, we will be able to reset the
    // communicators without unallocated it
    PMMG_DEL_MEM(parmesh,ext_edge_comm->itosend,int,"i2send");
    PMMG_DEL_MEM(parmesh,ext_edge_comm->itorecv,int,"i2recv");
  }

  PMMG_DEL_MEM(parmesh,request,MPI_Request,"mpi request array");
  PMMG_DEL_MEM(parmesh,status,MPI_Status,"mpi status array");
  PMMG_DEL_MEM(parmesh,i2send_size,int,"size of the i2send array");
  PMMG_DEL_MEM(parmesh,i2recv_size,int,"size of the i2recv array");

  PMMG_DEL_MEM(parmesh,nitem_ext_comm,int,
               "number of items in each external communicator");

  PMMG_DEL_MEM(parmesh,list.item,PMMG_lnkdCell,"linked list array");

  return ier;
}

/**
 * \param parmesh pointer to parmesh structure
 * \param mesh pointer to the mesh structure
 * \param hpar hash table of parallel edges
 * \param comm pointer toward the MPI communicator to use: when called before
 * the first mesh balancing (at preprocessing stage) we have to use the
 * read_comm communicator (i.e. the communicator used to provide the inputs).
 * For all ather calls, comm has to be the communicator to use for computations.
 *
 * Build edge communicator.
 */
int PMMG_build_edgeComm( PMMG_pParMesh parmesh,MMG5_pMesh mesh,MMG5_HGeom *hpar,MPI_Comm comm) {
  PMMG_pGrp      grp;
  PMMG_pInt_comm int_face_comm,int_edge_comm;
  PMMG_pExt_comm ext_face_comm,ext_edge_comm;
  MMG5_pTetra    pt;
  MMG5_pEdge     pa;
  MMG5_hgeom     *ph;
  int            *nitems_ext_comm,color,k,i,idx,ie,ifac,iloc,j,item;
  int            edg;
  uint16_t       tag;
  int8_t         ia,i1,i2;

  assert( parmesh->ngrp == 1 );
  grp = &parmesh->listgrp[0];

  /** Build the internal edge communicator. It already contains ALL possible
   *  parallel edges (even for star configurations) every parallel edge
   *  necessarily belongs to a parallel face (unless the underlying global mesh
   *  is not connected).
   */
  if( !PMMG_build_intEdgeComm( parmesh,mesh,hpar ) ) return 0;

  int_face_comm = parmesh->int_face_comm;
  int_edge_comm = parmesh->int_edge_comm;

  /* Allocate internal communicator */
  PMMG_CALLOC(parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,"int_face_comm",return 0);


  /** Count edges in each external communicator seen from the face ones */
  PMMG_CALLOC(parmesh,nitems_ext_comm,parmesh->nprocs,int,"nitems_ext_comm",return 0);

  /* Expose face index to the external communicator */
  for( i = 0; i < grp->nitem_int_face_comm; i++ ) {
    k   = grp->face2int_face_comm_index1[i];
    idx = grp->face2int_face_comm_index2[i];
    int_face_comm->intvalues[idx] = k;
  }

  /* For each face communicator, get the edges */
  for( k = 0; k < parmesh->next_face_comm; k++ ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    color = ext_face_comm->color_out;
    for( i = 0; i < ext_face_comm->nitem; i++ ) {
      /* Get face */
      idx  =  ext_face_comm->int_comm_index[i];
      ie   =  int_face_comm->intvalues[idx]/12;
      ifac = (int_face_comm->intvalues[idx]%12)/3;
      iloc = (int_face_comm->intvalues[idx]%12)%3;
      /* Get face edges */
      pt = &mesh->tetra[ie];
      assert( MG_EOK(pt) );
      for( j = 0; j < 3; j++ ) {
        /* Take the edge opposite to vertex iloc on face ifac */
        i1 = MMG5_idir[ifac][(iloc+j+1)%3];
        i2 = MMG5_idir[ifac][(iloc+j+2)%3];
        if ( !MMG5_hGet( hpar, pt->v[i1], pt->v[i2], &edg, &tag ) ) return 0;
        pa = &mesh->edge[edg];
        /* Overwrite edge base with current color */
        if( pa->base != color ) {
          /* Count edge and mark it.
           * ext_face_comm are already ordered; use common face point to order
           * the edge communicator. */
          nitems_ext_comm[color]++;
          pa->base = color;
        }
      }
    }
  }

  /** First attempt to build the external communicator from the face ones */
  parmesh->next_edge_comm = parmesh->next_face_comm;
  PMMG_CALLOC(parmesh,parmesh->ext_edge_comm,parmesh->next_edge_comm,PMMG_Ext_comm,"ext_edge_comm",return 0);

  /* For each face communicator, fill the edge communicator */
  for( k = 0; k < parmesh->next_face_comm; k++ ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    color = ext_face_comm->color_out;

    ext_edge_comm = &parmesh->ext_edge_comm[k];
    ext_edge_comm->nitem = nitems_ext_comm[color];
    ext_edge_comm->color_in = parmesh->myrank;
    ext_edge_comm->color_out = color;
    PMMG_CALLOC(parmesh,ext_edge_comm->int_comm_index,ext_edge_comm->nitem,int,"int_comm_index",return 0);

    item = 0;
    for( i = 0; i < ext_face_comm->nitem; i++ ) {
      /* Get face */
      idx  =  ext_face_comm->int_comm_index[i];
      ie   =  int_face_comm->intvalues[idx]/12;
      ifac = (int_face_comm->intvalues[idx]%12)/3;
      iloc = (int_face_comm->intvalues[idx]%12)%3;
      /* Get face edges */
      pt = &mesh->tetra[ie];
      assert( MG_EOK(pt) );
      /* ext_face_comm are already ordered; use common face point to travel the
       * edges on the face in the same order. */
      if( parmesh->myrank < color ) {
        for( j = 0; j < 3; j++ ) {
          if( !PMMG_fillExtEdgeComm_fromFace( parmesh,mesh,hpar,ext_edge_comm,pt,ifac,iloc,j,parmesh->nprocs+color,&item ) ) return 0;
        }
      } else {
        for( j = 3; j > 0; j-- ) {
          if( !PMMG_fillExtEdgeComm_fromFace( parmesh,mesh,hpar,ext_edge_comm,pt,ifac,iloc,j%3,parmesh->nprocs+color,&item ) ) return 0;
        }
      }
    }
    assert( item == ext_edge_comm->nitem );
  }

  /** Complete the external edge communicator */
  if( !PMMG_build_completeExtEdgeComm( parmesh,comm ) ) return 0;

  /* Reorder edge nodes */
  if( !PMMG_color_commNodes( parmesh,comm ) ) return 0;
  MMG5_pPoint ppt0,ppt1;
  int swp;
  for( k = 1; k <= mesh->na; k++ ) {
    pa = &mesh->edge[k];
    ppt0 = &mesh->point[pa->a];
    ppt1 = &mesh->point[pa->b];
    /* Swap nodes so that the first one has the highest global label */
    if( ppt0->tmp < ppt1->tmp ){
      swp = pa->a;
      pa->a = pa->b;
      pa->b = swp;
    }
  }

  /** Check the external edge communicator */
  assert( PMMG_check_extEdgeComm( parmesh,comm ) );

  /* Free */
  if ( int_face_comm ) {
    PMMG_DEL_MEM(parmesh,int_face_comm->intvalues,int,"int_face_comm");
  }
  PMMG_DEL_MEM(parmesh,nitems_ext_comm,int,"nitem_int_face_comm");

  return 1;
}

/**
 * \param parmesh pointer to parmesh structure
 *
 * - Convert triangle index to tetra-face-node triplet, using the "cc" triangle
 *   field initialized by MMG5_bdryTria, and the global node indices temporarily
 *   stored in the points flag.
 * - Store the index triplet in group communicator index 1,
 * - Tag corresponding triangle edges and nodes as PARBDY.
 */
void PMMG_tria2elmFace_flags( PMMG_pParMesh parmesh ) {
  PMMG_pGrp    grp;
  MMG5_pMesh   mesh;
  MMG5_pTria   ptt;
  MMG5_pTetra  pt;
  int          kt,ie,ifac,iploc;
  int          i,imax,iloc,iglob;

  /* Only one group */
  grp = &parmesh->listgrp[0];
  mesh = grp->mesh;

  /* Process tria stored in index1 */
  for( i=0; i<grp->nitem_int_face_comm; i++ ) {
    kt    = grp->face2int_face_comm_index1[i];
    ptt   = &mesh->tria[kt];
    ie    = ptt->cc/4;
    ifac  = ptt->cc%4;

    /* Get triangle node with highest global index */
    iploc = imax  = 0;
    for( iloc=0; iloc<3; iloc++ ) {
      iglob = mesh->point[ptt->v[iloc]].flag;
      if( iglob > imax ) {
        imax   = iglob;
        iploc = iloc;
      }
    }

    /* Store ie-ifac-iploc in index1 */
    grp->face2int_face_comm_index1[i] = 12*ie+3*ifac+iploc;

    /* Set triangle and nodes as parallel */
    PMMG_tag_par_tria(ptt);
    for( iloc = 0; iloc < 3; iloc++ )
      PMMG_tag_par_node(&mesh->point[ptt->v[iloc]]);
  }
}

/**
 * \param parmesh pointer to parmesh structure
 *
 * - Convert triangle index to tetra-face-node triplet, using the "cc" triangle
 *   field initialized by MMG5_bdryTria, and the nodes coordinates.
 * - Store the index triplet in group communicator index 1,
 * - Tag corresponding triangle edges and nodes as PARBDY.
 */
void PMMG_tria2elmFace_coords( PMMG_pParMesh parmesh ) {
  PMMG_pGrp    grp;
  MMG5_pMesh   mesh;
  MMG5_pTria   ptt;
  MMG5_pTetra  pt;
  MMG5_pPoint  ppt;
  int          kt,ie,ifac,iploc;
  double       cmax[3];
  int          i,idim,iloc;

  /* Only one group */
  grp = &parmesh->listgrp[0];
  mesh = grp->mesh;

  /* Process tria stored in index1 */
  for( i=0; i<grp->nitem_int_face_comm; i++ ) {
    kt    = grp->face2int_face_comm_index1[i];
    ptt    = &mesh->tria[kt];
    ie     = ptt->cc/4;
    ifac   = ptt->cc%4;

    /* Get triangle node with highest coordinates */
    iploc = PMMG_tria_highestcoord(mesh,ptt->v);

    /* Store ie-ifac-iploc in index1 */
    grp->face2int_face_comm_index1[i] = 12*ie+3*ifac+iploc;

    /* Set triangle and nodes as parallel */
    PMMG_tag_par_tria(ptt);
    for( iloc = 0; iloc < 3; iloc++ )
      PMMG_tag_par_node(&mesh->point[ptt->v[iloc]]);
  }
}

/**
 * \param mesh  pointer toward the mesh structure
 * \param ptt_v indices of a triangle vertices
 * \return iploc (0, 1 or 2) local node index
 *
 * Get triangle node with highest coordinates
 *
 */
int PMMG_tria_highestcoord( MMG5_pMesh mesh, MMG5_int *ptt_v) {
  MMG5_pPoint  ppt;
  int          idim,iloc,iploc;
  double       cmax[3];

  /* Get triangle node with highest coordinates */
  iploc  = 0;
  ppt = &mesh->point[ptt_v[0]];
  cmax[0] = ppt->c[0];
  cmax[1] = ppt->c[1];
  cmax[2] = ppt->c[2];
  for( iloc=1; iloc<3; iloc++ ) {
    ppt = &mesh->point[ptt_v[iloc]];
    for( idim=0; idim<3; idim++ ) {
      if( ppt->c[idim] - cmax[idim] < -MMG5_EPSOK*20 ) break;
      if( ppt->c[idim] - cmax[idim] >  MMG5_EPSOK*20 ) {
        cmax[0] = ppt->c[0];
        cmax[1] = ppt->c[1];
        cmax[2] = ppt->c[2];
        iploc  = iloc;
        break;
      }
    }
  }

  return iploc;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \return 0 if fail, 1 if success
 *
 * Set group and internal communicator indexing, given mesh entities temporarily
 * stored in the external node communicator.
 *
 */
int PMMG_build_nodeCommIndex( PMMG_pParMesh parmesh ) {
  PMMG_pGrp      grp;
  PMMG_pInt_comm int_node_comm;
  PMMG_pExt_comm ext_node_comm;
  MMG5_pMesh     mesh;
  MMG5_pPoint    ppt;
  int            ip,nitem_int_node_comm,iext_comm,iext,iint;

  /* Only one group */
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;

  /* Use flag to mark visited points:
   * - points are initially marked with PMMG_NUL;
   * - when counting int comm items, points are marked with PMMG_UNSET;
   * - when filling communicators, points are marked with their index in the
   *   internal comm.
   */
  for( iext_comm = 0; iext_comm < parmesh->next_node_comm; iext_comm++ ) {
    ext_node_comm = &parmesh->ext_node_comm[iext_comm];
    for( iext = 0; iext < ext_node_comm->nitem; iext++ ) {
      ip = ext_node_comm->int_comm_index[iext];
      ppt = &mesh->point[ip];
      ppt->flag = PMMG_NUL;
    }
  }

  /** 1) Count internal communicator items (the same node can be in several
   *     external communicators) */
  nitem_int_node_comm = 0;
  for( iext_comm = 0; iext_comm < parmesh->next_node_comm; iext_comm++ ) {
    ext_node_comm = &parmesh->ext_node_comm[iext_comm];
    for( iext = 0; iext < ext_node_comm->nitem; iext++ ) {
      ip = ext_node_comm->int_comm_index[iext];
      ppt = &mesh->point[ip];
      if( ppt->flag == PMMG_UNSET ) continue;
      ppt->flag = PMMG_UNSET;
      nitem_int_node_comm++;
    }
  }

  /* Allocate group communicators */
  int_node_comm = parmesh->int_node_comm;
  int_node_comm->nitem = nitem_int_node_comm;
  PMMG_CALLOC(parmesh,grp->node2int_node_comm_index1,nitem_int_node_comm,int,"node2int_node_comm_index1",return 0);
  PMMG_CALLOC(parmesh,grp->node2int_node_comm_index2,nitem_int_node_comm,int,"node2int_node_comm_index2",return 0);
  grp->nitem_int_node_comm = nitem_int_node_comm;

  /** 2) Set communicators indexing */
  iint = 0;
  for( iext_comm = 0; iext_comm < parmesh->next_node_comm; iext_comm++ ) {
    ext_node_comm = &parmesh->ext_node_comm[iext_comm];
    for( iext = 0; iext < ext_node_comm->nitem; iext++ ) {
      ip = ext_node_comm->int_comm_index[iext];
      ppt = &mesh->point[ip];
      if( ppt->flag == PMMG_UNSET ) {
        /* New point */
        ppt->flag = iint;
        grp->node2int_node_comm_index1[iint] = ip;
        grp->node2int_node_comm_index2[iint] = iint;
        ext_node_comm->int_comm_index[iext] = iint++;
      } else {
        /* Point already in the int communicator: only the ext comm needs to be
         * updated */
        ext_node_comm->int_comm_index[iext] = ppt->flag;
      }
    }
  }
  assert( iint == nitem_int_node_comm );

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \return 0 if fail, 1 if success
 *
 * Set group and internal communicator indexing, given mesh entities temporarily
 * stored in the external face communicator.
 *
 */
int PMMG_build_faceCommIndex( PMMG_pParMesh parmesh, MMG5_int* permtria ) {
  PMMG_pGrp      grp;
  PMMG_pInt_comm int_face_comm;
  PMMG_pExt_comm ext_face_comm;
  int            nitem_int_face_comm,iext_comm,iext,iint;

  /* Only one group */
  grp = &parmesh->listgrp[0];

  /* Count internal communicator items (as each face can be only in one
   * ext_comm, simply ssum the nb of items in each ext_comm) */
  nitem_int_face_comm = 0;
  for( iext_comm = 0; iext_comm < parmesh->next_face_comm; iext_comm++ ) {
    ext_face_comm = &parmesh->ext_face_comm[iext_comm];
    nitem_int_face_comm += ext_face_comm->nitem;
  }

  /* Allocate group communicators */
  int_face_comm = parmesh->int_face_comm;
  int_face_comm->nitem = nitem_int_face_comm;
  PMMG_CALLOC(parmesh,grp->face2int_face_comm_index1,nitem_int_face_comm,int,"face2int_face_comm_index1",return 0);
  PMMG_CALLOC(parmesh,grp->face2int_face_comm_index2,nitem_int_face_comm,int,"face2int_face_comm_index2",return 0);
  grp->nitem_int_face_comm = nitem_int_face_comm;

  /* Set communicators indexing (faces in the ext comms are concatenated into
   * the internal comm) */
  iint = 0;
  for( iext_comm = 0; iext_comm < parmesh->next_face_comm; iext_comm++ ) {
    ext_face_comm = &parmesh->ext_face_comm[iext_comm];
    for( iext = 0; iext < ext_face_comm->nitem; iext++ ) {
      if (permtria) {
        grp->face2int_face_comm_index1[iint] = permtria[ext_face_comm->int_comm_index[iext]];
      }
      else {
        grp->face2int_face_comm_index1[iint] = ext_face_comm->int_comm_index[iext];
      }
      grp->face2int_face_comm_index2[iint] = iint;
      ext_face_comm->int_comm_index[iext] = iint++;
    }
  }
  assert( iint == nitem_int_face_comm );

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 *
 * \param comm pointer toward the MPI communicator to use: when called before
 * the first mesh balancing (at preprocessing stage) we have to use the
 * read_comm communicator (i.e. the communicator used to provide the inputs).
 * For all ather calls, comm has to be the communicator to use for computations.
 *
 * \return 0 if fail, 1 if succeed
 *
 * Construction of the face communicators from the node communicators
 *
 */
int PMMG_build_faceCommFromNodes( PMMG_pParMesh parmesh,MPI_Comm comm ) {
  PMMG_pExt_comm ext_node_comm;
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_pTria     ptt;
  MMG5_Hash      hash;
  int            **local_index,**global_index;
  int            *fNodes_loc,*fNodes_par,*fColors;
  int            nb_fNodes_loc,*nb_fNodes_par,*displs,*counter,*iproc2comm;
  int            kt,ia,ib,ic,i,icomm,iproc,iloc,iglob,myrank,next_face_comm,ier;

  myrank = parmesh->myrank;
  grp    = &parmesh->listgrp[0];
  mesh   = grp->mesh;

  /** 1) Store global node ids in point flags */
  /* Reset point flags */
  for( i=1; i<=mesh->np; i++ )
    mesh->point[i].flag = PMMG_NUL;

  /* Loop on ext node communicators to get global node IDs */
  for( icomm=0; icomm<parmesh->next_node_comm; icomm++ ) {
    ext_node_comm = &parmesh->ext_node_comm[icomm];
    if( !ext_node_comm->nitem_to_share ) continue;
    for( i=0; i<ext_node_comm->nitem_to_share; i++ ) {
      iloc  = ext_node_comm->itosend[i];
      iglob = ext_node_comm->itorecv[i];
      mesh->point[iloc].flag = iglob;
    }
    /* Deallocate arrays that have been used to store local/global interface
     * nodes enuleration */
    PMMG_DEL_MEM(parmesh,ext_node_comm->itosend,int,"ext comm itosend array");
    PMMG_DEL_MEM(parmesh,ext_node_comm->itorecv,int,"ext comm itorecv array");
    ext_node_comm->nitem_to_share = 0;
  }

  /** 2) Hash triangles with global node index: This avoids the occurrence of
   * non-boundary faces connected to three parallel nodes. */
  nb_fNodes_loc = 3*mesh->nt;
  PMMG_CALLOC(parmesh,fNodes_loc,nb_fNodes_loc,int,"fNodes_loc",return 0);
  PMMG_MALLOC(parmesh,fColors,2*mesh->nt,int,"fColors",return 0);
  for( i=0; i<2*mesh->nt; i++ )
    fColors[i] = PMMG_UNSET;

  if ( ! MMG5_hashNew(mesh,&hash,0.51*mesh->nt,1.51*mesh->nt) ) return 0;

  for (kt=1; kt<=mesh->nt; kt++) {
    ptt = &mesh->tria[kt];
    ia = mesh->point[ptt->v[0]].flag;
    ib = mesh->point[ptt->v[1]].flag;
    ic = mesh->point[ptt->v[2]].flag;
    if( ia && ib && ic ) {
      if ( !MMG5_hashFace(mesh,&hash,ia,ib,ic,kt) ) {
        MMG5_DEL_MEM(mesh,hash.item);
        return 0;
      }
      fNodes_loc[3*(kt-1)+0] = ia;
      fNodes_loc[3*(kt-1)+1] = ib;
      fNodes_loc[3*(kt-1)+2] = ic;
    }
  }


  /** 3) Communicate face nodes */
  PMMG_CALLOC(parmesh,nb_fNodes_par,parmesh->nprocs,int,"nb_fNodes_par",return 0);
  PMMG_CALLOC(parmesh,displs,parmesh->nprocs+1,int,"displs",return 0);

  MPI_CHECK( MPI_Allgather(&nb_fNodes_loc,1,MPI_INT,nb_fNodes_par,1,MPI_INT,comm),
             return 0 );
  
  displs[0] = 0;
  for( i=0; i<parmesh->nprocs; i++ )
    displs[i+1] = displs[i]+nb_fNodes_par[i];

  PMMG_CALLOC(parmesh,fNodes_par,displs[parmesh->nprocs],int,"fNodes_par",return 0);

  MPI_CHECK( MPI_Allgatherv(fNodes_loc,nb_fNodes_loc,MPI_INT,
                            fNodes_par,nb_fNodes_par,displs,MPI_INT,comm), return 0);


  /** 4) For each proc pair, get "other" tria in the local hash table,
   * count and store tria color. */
  PMMG_CALLOC(parmesh,counter,parmesh->nprocs,int,"counter",return 0);
  PMMG_CALLOC(parmesh,iproc2comm,parmesh->nprocs,int,"iproc2comm",return 0);
  for( iproc=0; iproc<parmesh->nprocs; iproc++ )
    iproc2comm[iproc] = PMMG_UNSET;

  /* Get face colors and count faces for each color */
  icomm = 0;
  for( iproc=0; iproc<parmesh->nprocs; iproc++ ) {
    if( iproc == myrank ) continue;
    for( i=displs[iproc]/3; i<displs[iproc+1]/3; i++ ) {
      ia = fNodes_par[3*i+0];
      ib = fNodes_par[3*i+1];
      ic = fNodes_par[3*i+2];
      if( ia && ib && ic ) {
        kt = MMG5_hashGetFace(&hash,ia,ib,ic);
        if( kt ) { /*it can be zero if (i) is an internal boundary on iproc  */
          if( iproc2comm[iproc] == PMMG_UNSET ) iproc2comm[iproc] = icomm++;
          /* Store face color and face global ID (starting from 1) on the other
           * proc */
          fColors[2*(kt-1)+0] = iproc;
          fColors[2*(kt-1)+1] = i+1;
          counter[iproc]++;
        }
      }
    }
  }
  assert( iproc2comm[myrank] == PMMG_UNSET );


  /** 5) Fill face communicators. */

  /* Set nb of communicators */
  next_face_comm = 0;
  for( iproc=0; iproc<parmesh->nprocs; iproc++ ) {
    if( iproc2comm[iproc] != PMMG_UNSET ) next_face_comm++;
  }
  ier = PMMG_Set_numberOfFaceCommunicators(parmesh,next_face_comm);

  PMMG_CALLOC(parmesh, local_index,next_face_comm,int*, "local_index pointer",return 0);
  PMMG_CALLOC(parmesh,global_index,next_face_comm,int*,"global_index pointer",return 0);
  for( iproc=0; iproc<parmesh->nprocs; iproc++ ) {
    if( iproc2comm[iproc] == PMMG_UNSET ) continue;
    /* Set communicator size and reset counter */
    icomm = iproc2comm[iproc];
    PMMG_CALLOC(parmesh, local_index[icomm],counter[iproc],int, "local_index array",return 0);
    PMMG_CALLOC(parmesh,global_index[icomm],counter[iproc],int,"global_index array",return 0);
    ier = PMMG_Set_ithFaceCommunicatorSize(parmesh,icomm,iproc,counter[iproc]);
    counter[iproc] = 0;
  }

  /* Create injective, non-surjective global face enumeration */
  for (kt=1; kt<=mesh->nt; kt++) {
    iproc = fColors[2*(kt-1)];
    if( iproc == PMMG_UNSET ) continue;
    iglob = fColors[2*(kt-1)+1];
    icomm = iproc2comm[iproc];
    i = counter[iproc]++;
    local_index[icomm][i] = kt;
    global_index[icomm][i] = MG_MIN(displs[myrank]/3+kt,iglob);
  }
 
 
  /* Fill and sort each communicator */
  for( iproc=0; iproc<parmesh->nprocs; iproc++ ) {
    if( iproc2comm[iproc] == PMMG_UNSET ) continue;
    icomm = iproc2comm[iproc];
    ier = PMMG_Set_ithFaceCommunicator_faces( parmesh, icomm, local_index[icomm],
                                              global_index[icomm], 1 );
  }

  /** 6) Set communicators indexing, convert tria index into iel face index */
  ier = PMMG_build_faceCommIndex( parmesh, NULL );
  PMMG_tria2elmFace_flags( parmesh );


  /* Free memory */
  MMG5_DEL_MEM(mesh,hash.item);
  PMMG_DEL_MEM(parmesh,fColors,int,"fColors");
  PMMG_DEL_MEM(parmesh,fNodes_loc,int,"fNodes_loc");
  PMMG_DEL_MEM(parmesh,fNodes_loc,int,"fNodes_par");
  PMMG_DEL_MEM(parmesh,nb_fNodes_par,int,"nb_fNodes_par");
  PMMG_DEL_MEM(parmesh,fNodes_par,int,"fNodes_par");
  PMMG_DEL_MEM(parmesh,displs,int,"displs");
  PMMG_DEL_MEM(parmesh,counter,int,"counter");
  PMMG_DEL_MEM(parmesh,iproc2comm,int,"iproc2comm");
  for( icomm=0; icomm<next_face_comm; icomm++ ) {
    PMMG_DEL_MEM(parmesh, local_index[icomm],int, "local_index array");
    PMMG_DEL_MEM(parmesh,global_index[icomm],int,"global_index array");
  }
  PMMG_DEL_MEM(parmesh, local_index,int*, "local_index pointer");
  PMMG_DEL_MEM(parmesh,global_index,int*,"global_index pointer");

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param comm pointer toward the MPI communicator to use: when called before
 * the first mesh balancing (at preprocessing stage) we have to use the
 * read_comm communicator (i.e. the communicator used to provide the inputs).
 * For all ather calls, comm has to be the communicator to use for computations.
 *
 * \return 1 if success, 0 if fail.
 *
 * Build the node communicators (externals and internals) from the faces ones.
 *
 */
int PMMG_build_nodeCommFromFaces( PMMG_pParMesh parmesh, MPI_Comm comm ) {
  int ier, ier_glob;

  assert ( PMMG_check_extFaceComm ( parmesh,comm ) );
  assert ( PMMG_check_intFaceComm ( parmesh ) );

  /** Build the internal node communicator from the faces ones */
  ier = PMMG_build_intNodeComm(parmesh);

  if ( !ier ) {
    fprintf(stderr,"\n  ## Error: %s: unable to build the internal node"
            " communicators from the internal faces communicators.\n",__func__);
  }
  else {
    /* The internal comm is filled, try to begin to fill the external one */
    assert ( PMMG_check_intNodeComm ( parmesh ) );

    /** Build the external node communicators from the faces ones without looking
     * at the possibility to miss a processor that is not visible from the face
     * communicators: _0_\1/_2_ configuration ( from 0 it is not possible to say
     * that the node at edge intersction belongs to 2 */
    ier = PMMG_build_simpleExtNodeComm(parmesh);
    if ( !ier ) {
      fprintf(stderr,"\n  ## Error: %s: unable to build the simple externals node"
              " communicators from the external faces communicators.\n",__func__);
    }
  }

  /* Check that all steps have successed until here (because the next function
   * involves MPI comms) */
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, comm);
  if ( !ier_glob ) return 0;

  /** Fill the external node communicator */
  ier = PMMG_build_completeExtNodeComm(parmesh,comm);
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, comm);
  if ( !ier ) {
    fprintf(stderr,"\n  ## Error: %s: unable to complete the external node"
            " communicators.\n",__func__);
  }
  if ( !ier_glob ) return 0;

  PMMG_MALLOC(parmesh,parmesh->int_node_comm->intvalues,
              parmesh->int_node_comm->nitem,int,"intvalues",return 0);

  if ( !PMMG_pack_nodeCommunicators(parmesh) ) return 0;

  if ( parmesh->int_node_comm->intvalues )
    PMMG_DEL_MEM(parmesh,parmesh->int_node_comm->intvalues,int,"intvalues");

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
  PMMG_pExt_comm  ext_node_comm,ext_face_comm;
  PMMG_pInt_comm  int_node_comm,int_face_comm;
  MMG5_pMesh      mesh;
  MMG5_pTetra     pt;
  MMG5_pPoint     ppt;
  int             *face_vertices,*flag,idx1,idx2;
  int             *node2int_node_comm_index1,*node2int_node_comm_index2;
  int             *face2int_face_comm_index1,*face2int_face_comm_index2;
  int             next_face_comm,next_node_comm,nitem_ext_comm;
  int             color_out,ier,i,j,k,iel,ifac,ip,iploc,grpid;

  ier = 0;
  assert ( !parmesh->ext_node_comm && "external comms must be deleted" );

  next_face_comm = parmesh->next_face_comm;
  next_node_comm = next_face_comm;

  PMMG_CALLOC(parmesh,parmesh->ext_node_comm,next_node_comm,PMMG_Ext_comm,
              "ext_node_comm ",return 0);
  parmesh->next_node_comm = next_node_comm;

  /** Step 1: For each face, store the position of its vertices in the internal
   * communicator */
  int_face_comm = parmesh->int_face_comm;
  PMMG_CALLOC(parmesh,face_vertices,3*int_face_comm->nitem,int,
              "position of face vertices in int_node_comm",goto end);

  for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
    grp  = &parmesh->listgrp[grpid];
    mesh = grp->mesh;

    for ( k=1; k<=mesh->np; ++k )
      mesh->point[k].tmp = PMMG_UNSET;

    node2int_node_comm_index1 = grp->node2int_node_comm_index1;
    node2int_node_comm_index2 = grp->node2int_node_comm_index2;
    for ( i=0; i<grp->nitem_int_node_comm; ++i ) {
      idx1 = node2int_node_comm_index1[i];
      idx2 = node2int_node_comm_index2[i];

      assert ( mesh->point[idx1].tmp < 0 && "node is not unique in the comm" );
      mesh->point[idx1].tmp = idx2;
    }

    face2int_face_comm_index1 = grp->face2int_face_comm_index1;
    face2int_face_comm_index2 = grp->face2int_face_comm_index2;
    for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
      iel   =  face2int_face_comm_index1[i]/12;
      ifac  = (face2int_face_comm_index1[i]%12)/3;
      iploc = (face2int_face_comm_index1[i]%12)%3;

      idx2 = face2int_face_comm_index2[i];

      pt = &mesh->tetra[iel];
      for ( j=0; j<3; ++j ) {
        ip  = pt->v[MMG5_idir[ifac][(j+iploc)%3]];
        ppt = &mesh->point[ip];
        assert ( ppt->tmp>=0 && "missing node in the communicator" );

        face_vertices[3*idx2+j] = ppt->tmp;
      }
    }
  }

  /** Step 2: fill the external communicators */
  int_node_comm = parmesh->int_node_comm;
  PMMG_CALLOC(parmesh,flag,int_node_comm->nitem,int,"node flag",goto end);

  for ( k=0; k<next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    color_out = ext_face_comm->color_out;

    parmesh->ext_node_comm[k].color_in  = ext_face_comm->color_in;
    parmesh->ext_node_comm[k].color_out = color_out;

    if ( ext_face_comm->nitem )
      /* Reset the node flags */
      for ( j=0; j<int_node_comm->nitem; ++j ) flag[j] = 0;

    /* Count the number of nodes in the communicator */
    nitem_ext_comm = 0;
    for ( i=0; i<ext_face_comm->nitem; ++i ) {
      idx2  = ext_face_comm->int_comm_index[i];

      for ( j=0; j<3; ++j ) {
        ip = face_vertices[3*idx2+j];
        assert ( ip >=0 );

        if ( flag[ip] > 0 ) continue;

        flag[ip] = 1;
        ++nitem_ext_comm;
      }
    }

    /* External communicator allocation */
    PMMG_CALLOC(parmesh,parmesh->ext_node_comm[k].int_comm_index,nitem_ext_comm,
                int,"external node communicator",goto end);
    parmesh->ext_node_comm[k].nitem = nitem_ext_comm;
    ext_node_comm = &parmesh->ext_node_comm[k];

    /* Process the external face communicator and fill the node one */
    if ( ext_face_comm->nitem )
      /* Reset the node flags */
      for ( j=0; j<int_node_comm->nitem; ++j ) flag[j] = 0;

    nitem_ext_comm = 0;

    if ( color_out > parmesh->ext_node_comm[k].color_in ) {
      /* Travel the face vertices in the natural direction */
      for ( i=0; i<ext_face_comm->nitem; ++i ) {
        idx2  = ext_face_comm->int_comm_index[i];

        for  ( j=0; j<3; ++j ) {
          ip = face_vertices[3*idx2+j];
          assert ( ip>=0 );

          if ( flag[ip] > 0 ) continue;
          flag[ip] = 1;
          ext_node_comm->int_comm_index[nitem_ext_comm++] = ip;
        }
      }
    }
    else {
      /* Travel the face vertices in the opposite direction */
      for ( i=0; i<ext_face_comm->nitem; ++i ) {
        idx2  = ext_face_comm->int_comm_index[i];

        for  ( j=0; j<3; ++j ) {
          ip = face_vertices[3*idx2+(3-j)%3];
          assert ( ip>=0 );

          if ( flag[ip] > 0 ) continue;
          flag[ip] = 1;
          ext_node_comm->int_comm_index[nitem_ext_comm++] = ip;
        }
      }
    }
    assert ( nitem_ext_comm == ext_node_comm->nitem );
  }

  /* Success */
  ier = 1;

end:
  if ( !ier ) {
    /* Deallocation of the external comms because we have failed */
    if ( next_node_comm ) {
      for ( k=0; k<next_node_comm; ++k ) {
        PMMG_DEL_MEM(parmesh,parmesh->ext_node_comm[k].int_comm_index,
                     int,"external node communicator");
      }
    }
    PMMG_DEL_MEM(parmesh,parmesh->ext_node_comm,PMMG_Ext_comm, "ext_node_comm ");
  }

  PMMG_DEL_MEM(parmesh,flag,int,"node flag");

  PMMG_DEL_MEM(parmesh,face_vertices,int,"position of face vertices in int_node_comm");

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
  PMMG_coorCell   *coor_list;
  MMG5_pMesh      mesh;
  MMG5_pTetra     pt;
  MMG5_pPoint     ppt;
  double          bb_min[3],bb_max[3],delta;
  int             *face2int_face_comm_index1;
  int             *node2int_node_comm_index1,*node2int_node_comm_index2;
  int             *shared_fac,*new_pos,nitem_node,first_nitem_node,pos;
  int             *face_vertices,ier,i,j,iel,ifac,ip,iploc,grpid,idx,fac_idx;
  int             nitem_node_init;
  int8_t          update;
#ifndef NDEBUG
  double dd,dist[3];
#endif

  ier = 0;

  /** Step 1: give a unique position in the internal communicator for each mesh
   * point but don't care about the unicity of the position for a point shared
   * by multiple groups */

  assert ( !parmesh->int_node_comm->nitem );
  nitem_node = 0;

  PMMG_MALLOC(parmesh,shared_fac,parmesh->int_face_comm->nitem,int,
              "Faces shared by 2 groups",goto end);
  for ( i=0; i<parmesh->int_face_comm->nitem; ++i )
    shared_fac[i] = PMMG_UNSET;

  for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
    /* Reset the flag field of the nodes of the group */
    grp  = &parmesh->listgrp[grpid];
    mesh = grp->mesh;

    for ( ip=1; ip<=mesh->np; ++ip )
      mesh->point[ip].tmp  = PMMG_UNSET;

    face2int_face_comm_index1 = grp->face2int_face_comm_index1;

    first_nitem_node = nitem_node;
    for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
      idx  = face2int_face_comm_index1[i];
      iel  =  idx/12;
      ifac = (idx%12)/3;

      assert ( iel && iel<=mesh->ne );
      assert ( 0<=ifac && ifac<4 );

      pt    = &mesh->tetra[iel];
      for ( j=0; j<3; ++j ) {
        ip  = pt->v[MMG5_idir[ifac][j]];

        assert ( ip && ip<=mesh->np );

        ppt = &mesh->point[ip];
        if ( ppt->tmp < 0 )
          /** Give a position in the internal communicator to this point */
          ppt->tmp = nitem_node++;
      }

      /* Mark the face as seen */
      ++shared_fac[grp->face2int_face_comm_index2[i]];
    }

    /* Allocations of the node2int_node arrays */
    PMMG_CALLOC(parmesh,grp->node2int_node_comm_index1,nitem_node-first_nitem_node,
                int,"node2int_node_comm_index1",goto end);

    grp->nitem_int_node_comm = nitem_node-first_nitem_node;
    PMMG_CALLOC(parmesh,grp->node2int_node_comm_index2,nitem_node-first_nitem_node,
                int,"node2int_node_comm_index2",goto end);

    /* Fill this arrays */
    node2int_node_comm_index1 = grp->node2int_node_comm_index1;
    node2int_node_comm_index2 = grp->node2int_node_comm_index2;
    idx = 0;
    for ( ip=1; ip<=mesh->np; ++ip ) {
      ppt = &mesh->point[ip];
      if ( ppt->tmp < 0 ) continue;

      node2int_node_comm_index1[idx] = ip;
      node2int_node_comm_index2[idx] = first_nitem_node+idx;
      ++idx;
    }
    assert ( idx == nitem_node-first_nitem_node && "missing item");
  }


  /** Step 2: remove some of the multiple positions and pack communicators */
  nitem_node_init = nitem_node;
  PMMG_MALLOC(parmesh,new_pos,nitem_node_init,int,"new pos in int_node_comm",goto end);
  PMMG_MALLOC(parmesh,face_vertices,3*parmesh->int_face_comm->nitem,int,
              "pos of face vertices in int_node_comm",goto end);

  /* Initialisation of the new_pos array with the current position */
 for ( i=0; i<nitem_node; ++i )
    new_pos[i] = i;

 /* Initialisation of the face_vertices array */
 for ( i=0; i<3*parmesh->int_face_comm->nitem; ++i )
   face_vertices[i]   = PMMG_UNSET;

 PMMG_CALLOC(parmesh,coor_list,nitem_node_init,PMMG_coorCell,"node coordinates",
             goto end);

#ifndef NDEBUG
 /* Store the node coordinates of the nodes */
 for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
   grp  = &parmesh->listgrp[grpid];
   mesh = grp->mesh;
   assert ( mesh->info.delta &&  "missing scaling infos");
   assert ( fabs(mesh->info.delta-1.)<MMG5_EPSD &&
            "scaled mesh... need to unscale it");

   for ( i=0; i<grp->nitem_int_node_comm; ++i ) {
     ip      = grp->node2int_node_comm_index1[i];
     idx     = grp->node2int_node_comm_index2[i];
     for ( j=0; j<3; ++j )
       coor_list[idx].c[j] = mesh->point[ip].c[j];
   }
 }
 /* Scale the coordinates depending to the bounding box ofthe internal comm */
 if ( !PMMG_scale_coorCellList(coor_list,nitem_node,bb_min,bb_max,&delta) )
   goto end;

#endif

  do {
    update = 0;
    /* Mark all the interface faces as not seen */
    for ( i=0; i<parmesh->int_face_comm->nitem; ++i )
      shared_fac[i] = -abs(shared_fac[i]);

    for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
      grp  = &parmesh->listgrp[grpid];
      mesh = grp->mesh;

      /* Reset tmp to the old position in the communicator */
      for ( i=0; i<grp->nitem_int_node_comm; ++i ) {
        ip = grp->node2int_node_comm_index1[i];
        mesh->point[ip].tmp  = grp->node2int_node_comm_index2[i];;
      }

      for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
        idx   = grp->face2int_face_comm_index1[i];
        iel   =  idx/12;
        ifac  = (idx%12)/3;
        iploc = (idx%12)%3;

        fac_idx = grp->face2int_face_comm_index2[i];

        pt = &mesh->tetra[iel];
        if ( shared_fac[fac_idx]<=0 ) {
          /* First time that we see this face: store it in the natural direction
           * and mark it as seen */
          shared_fac[fac_idx] = abs(shared_fac[fac_idx]);
          for ( j=0; j<3; ++j ) {
            ip  = pt->v[MMG5_idir[ifac][(j+iploc)%3]];
            ppt = &mesh->point[ip];
            assert ( ppt->tmp>=0 );

            pos = 3*fac_idx+j;
            if ( face_vertices[pos] < 0 ) {
              face_vertices[pos] = new_pos[mesh->point[ip].tmp];
              update = 1;
            }
            else if ( new_pos[face_vertices[pos]]>new_pos[mesh->point[ip].tmp] ) {
              new_pos[face_vertices[pos]] = new_pos[mesh->point[ip].tmp];
              update = 1;
            }
            else if ( new_pos[face_vertices[pos]]<new_pos[mesh->point[ip].tmp] ) {
              new_pos[mesh->point[ip].tmp] = new_pos[face_vertices[pos]];
              update = 1;
            }

#ifndef NDEBUG
            double scaled_coor[3];

            assert ( mesh->info.delta && "missing scaling infos" );
            assert ( fabs(mesh->info.delta-1.)<MMG5_EPSD &&
                     "scaled mesh... need to unscale it");

            /* Scale the point coor using the internal node comm scaling data */
            assert ( delta );
            dd = 1./delta;
            for ( j=0; j<3; ++j )
              scaled_coor[j] = dd*(mesh->point[ip].c[j]-bb_min[j]);

            /* Compute the distance between the points */
            dd = 0;
            for ( j=0; j<3; ++j ) {
              dist[j] = coor_list[new_pos[mesh->point[ip].tmp]].c[j]
                - scaled_coor[j];
              dd += dist[j]*dist[j];
            }
            assert ( dd < PMMG_EPSCOOR2 );
#endif
          }
        }
        else {
          /* Second time we see this face: store it in the opposite direction */
          for ( j=0; j<3; ++j ) {
            ip  = pt->v[MMG5_idir[ifac][(iploc+3-j)%3]];
            ppt = &mesh->point[ip];
            assert ( ppt->tmp>=0 );

            pos = 3*fac_idx+j;
            if ( face_vertices[pos] < 0 ) {
              face_vertices[pos] = new_pos[mesh->point[ip].tmp];
              update = 1;
            }
            else if ( new_pos[face_vertices[pos]]>new_pos[mesh->point[ip].tmp] ) {
              new_pos[face_vertices[pos]] = new_pos[mesh->point[ip].tmp];
              update = 1;
            }
            else if ( new_pos[face_vertices[pos]]<new_pos[mesh->point[ip].tmp] ) {
              new_pos[mesh->point[ip].tmp] = new_pos[face_vertices[pos]];
              update = 1;
            }
#ifndef NDEBUG
            double scaled_coor[3];

            /* Unscale the point coordinates using the internal node comm
             * scaling data */
            assert ( delta );
            dd = 1./delta;
            for ( j=0; j<3; ++j )
              scaled_coor[j] = dd*(mesh->point[ip].c[j]-bb_min[j]);

            /* Compute the distance between the points */
            dd = 0.;
            for ( j=0; j<3; ++j ) {
              dist[j] = coor_list[new_pos[mesh->point[ip].tmp]].c[j]
                -scaled_coor[j];
              dd += dist[j]*dist[j];
            }
            assert ( dd < PMMG_EPSCOOR2 );
#endif
          }
        }
      }
    }
  } while ( update );

  /* Update node2int_node_comm arrays */
  for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
    grp  = &parmesh->listgrp[grpid];

    for ( i=0; i<grp->nitem_int_node_comm; ++i ) {
      idx = grp->node2int_node_comm_index2[i];
      grp->node2int_node_comm_index2[i] = new_pos[idx];
    }
  }

  /* Repove the empty positions in int_node_comm */
  for ( i=0; i<nitem_node; ++i )
    new_pos[i] = PMMG_UNSET;

  j = 0;
  for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
    grp  = &parmesh->listgrp[grpid];

    for ( i=0; i<grp->nitem_int_node_comm; ++i ) {

      idx = grp->node2int_node_comm_index2[i];
      if ( new_pos[idx]<0 ) new_pos[idx] = j++;

      grp->node2int_node_comm_index2[i] = new_pos[idx];
    }
  }
  nitem_node = j;

  /** Step 3: Remove the last doublon by comparing the listed coordinates */
  /* Store the node coordinates of the nodes */
  for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
    grp  = &parmesh->listgrp[grpid];
    mesh = grp->mesh;

    for ( i=0; i<grp->nitem_int_node_comm; ++i ) {
      ip      = grp->node2int_node_comm_index1[i];
      idx     = grp->node2int_node_comm_index2[i];
      for ( j=0; j<3; ++j )
        coor_list[idx].c[j] = mesh->point[ip].c[j];
    }
  }
  /* Scale the coordinates depending to the bounding box of the internal comm */
  if ( !PMMG_scale_coorCellList(coor_list,nitem_node,bb_min,bb_max,&delta) )
    goto end;

  /* Store the point position in the internal communicator */
  for ( i=0; i<nitem_node; ++i ) {
    new_pos[i]       = i;
    coor_list[i].idx = i;
  }

  /* Travel the list and remove the identic nodes (use naive algorithm after
   * issues using point sorting). */
  if ( nitem_node ) {

    /* Detection of duplicated valies and assignation of a unique position in
     * the internal communicator. For now these positions are not contiguous. */
    for ( i = 0; i < nitem_node-1; ++i ) {

      if ( new_pos[coor_list[i].idx] < i ) {
        /* Point is a duplication so it has already been compared to everyone */
        continue;
      }

      for ( j=i+1; j<nitem_node; ++j ) {
        if ( !PMMG_compare_coorCell(&coor_list[i],&coor_list[j]) ) {
          new_pos[coor_list[j].idx] = new_pos[coor_list[i].idx];
        }
      }
    }
  }

  /* Update node2int_node_comm arrays */
  for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
    grp  = &parmesh->listgrp[grpid];

    for ( i=0; i<grp->nitem_int_node_comm; ++i ) {

      idx = grp->node2int_node_comm_index2[i];
      grp->node2int_node_comm_index2[i] = new_pos[idx];
    }
  }

  /* Repove the empty positions in int_node_comm */
  for ( i=0; i<nitem_node; ++i )
    new_pos[i] = PMMG_UNSET;

  j = 0;
  for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
    grp  = &parmesh->listgrp[grpid];

    for ( i=0; i<grp->nitem_int_node_comm; ++i ) {

      idx = grp->node2int_node_comm_index2[i];
      if ( new_pos[idx]<0 ) new_pos[idx] = j++;

      grp->node2int_node_comm_index2[i] = new_pos[idx];
    }
  }
  nitem_node = j;

  /** Step 4: Update the number of items in the internal node communicator */
  parmesh->int_node_comm->nitem = nitem_node;

  /* Success */
  ier = 1;

end:
  if ( !ier ) {
    /* Dealloc of the communicators because we have failed */
    for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
      grp  = &parmesh->listgrp[grpid];

      PMMG_DEL_MEM(parmesh,grp->node2int_node_comm_index1,
                   int,"node2int_node_comm_index1");
      PMMG_DEL_MEM(parmesh,grp->node2int_node_comm_index2,
                   int,"node2int_node_comm_index2");
    }
  }

  PMMG_DEL_MEM(parmesh,new_pos,int,"new pos in int_node_comm");
  PMMG_DEL_MEM(parmesh,face_vertices,int,"pos of face vertices in int_node_comm");
  PMMG_DEL_MEM(parmesh,shared_fac,int,"Faces shared by 2 groups");
  PMMG_DEL_MEM(parmesh,coor_list,PMMG_coorCell,"node coordinates");

  return ier;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param comm pointer toward the MPI communicator to use: when called before
 * the first mesh balancing (at preprocessing stage) we have to use the
 * read_comm communicator (i.e. the communicator used to provide the inputs).
 * For all ather calls, comm has to be the communicator to use for computations.
 *
 * \return 1 if success, 0 if fail.
 *
 * Complete the external communicators by travelling through the processors and
 * faces to detect all the processors to which each node belongs.
 *
 */
int PMMG_build_completeExtNodeComm( PMMG_pParMesh parmesh, MPI_Comm comm ) {
  PMMG_pExt_comm    ext_node_comm,*comm_ptr;
  PMMG_pInt_comm    int_node_comm;
  PMMG_cellLnkdList **proclists,list;
  int               *intvalues,nitem,nproclists,ier,ier2,k,i,j,idx,pos,rank,color;
  int               *itosend,*itorecv,*i2send_size,*i2recv_size,nitem2comm;
  int               *nitem_ext_comm,next_comm,val1_i,val2_i,val1_j,val2_j;
  int               alloc_size;
  int8_t            glob_update,loc_update;
  MPI_Request       *request;
  MPI_Status        *status;

  ier = 0;

  rank          = parmesh->myrank;
  int_node_comm = parmesh->int_node_comm;
  nitem         = int_node_comm->nitem;

  proclists       = NULL;
  comm_ptr        = NULL;
  request         = NULL;
  status          = NULL;
  i2send_size     = NULL;
  i2recv_size     = NULL;
  nitem_ext_comm  = NULL;
  list.item       = NULL;

  PMMG_CALLOC(parmesh,int_node_comm->intvalues,nitem,int,"node communicator",
    return 0);
  intvalues     = int_node_comm->intvalues;

  /** Use intvalues to flag the already treated points of the communicator:
   * initialization to 0.  */
  for ( k=0; k<nitem; ++k ) intvalues[k] = 0;

  PMMG_CALLOC(parmesh,proclists,nitem,PMMG_cellLnkdList*,"array of linked lists",
              goto end);

  /* Reallocation of the list of external comms at maximal size (nprocs) to
   * avoid tricky treatment when filling it.*/
  PMMG_REALLOC(parmesh,parmesh->ext_node_comm,parmesh->nprocs,
               parmesh->next_node_comm,PMMG_Ext_comm,
               "list of external communicators",goto end);

  next_comm = parmesh->next_node_comm;
  parmesh->next_node_comm = parmesh->nprocs;

  for ( k=next_comm; k<parmesh->nprocs; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    ext_node_comm->nitem          = 0;
    ext_node_comm->nitem_to_share = 0;
    ext_node_comm->int_comm_index = NULL;
    ext_node_comm->itosend        = NULL;
    ext_node_comm->itorecv        = NULL;
    ext_node_comm->rtosend        = NULL;
    ext_node_comm->rtorecv        = NULL;
    ext_node_comm->color_in       = rank;
    ext_node_comm->color_out      = PMMG_UNSET;
  }

  PMMG_CALLOC(parmesh,comm_ptr,parmesh->nprocs,PMMG_pExt_comm,
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

      PMMG_CALLOC(parmesh,proclists[idx],1,PMMG_cellLnkdList,"linked list pointer",
                  goto end);
      if ( !PMMG_cellLnkdListNew(parmesh,proclists[idx],idx,PMMG_LISTSIZE) ) goto end;

      if ( !PMMG_add_cell2lnkdList( parmesh,proclists[idx],rank,idx ) )
        goto end;

      intvalues[idx] = 1;
    }
  }
  /* Fill the missing pointer toward the empty external communicators */
  for ( k=0; k<parmesh->nprocs; ++k ) {
    if ( !comm_ptr[k] ) {

      assert ( next_comm < parmesh->nprocs );

      comm_ptr[k] = &parmesh->ext_node_comm[next_comm++];
      comm_ptr[k]->color_out = k;
    }
  }

  /** Step 2: While at least the proc list of 1 node is modified, send and
   * recieve the proc list of all the nodes to/from the other processors. At the
   * end of this loop, each node has the entire list of the proc to which it
   * belongs */
  alloc_size = parmesh->nprocs;
  PMMG_MALLOC(parmesh,request,    alloc_size,MPI_Request,"mpi request array",goto end);
  PMMG_MALLOC(parmesh,status,     alloc_size,MPI_Status,"mpi status array",goto end);
  PMMG_CALLOC(parmesh,i2send_size,alloc_size,int,"size of the i2send array",goto end);
  PMMG_CALLOC(parmesh,i2recv_size,alloc_size,int,"size of the i2recv array",goto end);

  if ( !PMMG_cellLnkdListNew(parmesh,&list,0,PMMG_LISTSIZE) ) goto end;

  do {
    glob_update = loc_update = 0;

    /** Send the list of procs to which belong each point of the communicator */
    for ( k=0; k<alloc_size; ++k ) {

      request[k] = MPI_REQUEST_NULL;
    }

    for ( k=0; k<parmesh->next_node_comm; ++k ) {
      ext_node_comm = &parmesh->ext_node_comm[k];

      /* Computation of the number of data to send to the other procs (we want
       * to send the the size of the linked list and the val1 and val2 fields of
       * each cell ) */
      nitem2comm = 0;
      if ( !ext_node_comm->nitem ) continue;

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
        pos += PMMG_packInArray_cellLnkdList(proclists[idx],&itosend[pos]);
      }
      assert ( pos==nitem2comm );

      MPI_CHECK( MPI_Isend(itosend,nitem2comm,MPI_INT,color,
                           MPI_COMMUNICATORS_NODE_TAG,comm,
                           &request[color]),goto end );
    }

    /** Recv the list of procs to which belong each point of the communicator */
    for ( k=0; k<parmesh->next_node_comm; ++k ) {
      ext_node_comm = &parmesh->ext_node_comm[k];

      if ( !ext_node_comm->nitem ) continue;

      color         = ext_node_comm->color_out;

      MPI_CHECK( MPI_Probe(color,MPI_COMMUNICATORS_NODE_TAG,comm,
                           &status[0] ),goto end);
      MPI_CHECK( MPI_Get_count(&status[0],MPI_INT,&nitem2comm),goto end);

      if ( i2recv_size[k] < nitem2comm ) {
        PMMG_REALLOC(parmesh,ext_node_comm->itorecv,nitem2comm,i2recv_size[k],
                     int,"itorecv",goto end);
        i2recv_size[k] = nitem2comm;
      }

      if ( nitem2comm ) {

        itorecv       = ext_node_comm->itorecv;
        MPI_CHECK( MPI_Recv(itorecv,nitem2comm,MPI_INT,color,
                            MPI_COMMUNICATORS_NODE_TAG,comm,
                            &status[0]), goto end );

        pos     = 0;
        for ( i=0; i<ext_node_comm->nitem; ++i ) {
          idx  = ext_node_comm->int_comm_index[i];
          assert ( idx>=0 );

          PMMG_reset_cellLnkdList( parmesh,&list );
          ier2 = PMMG_unpackArray_inCellLnkdList(parmesh,&list,&itorecv[pos]);

          if ( ier2 < 0 ) goto end;
          pos += ier2;

          ier2 = PMMG_merge_cellLnkdList(parmesh,proclists[idx],&list);
          if ( !ier2 ) goto end;

          loc_update |= (ier2%2);
        }
      }
    }

    MPI_CHECK( MPI_Waitall(alloc_size,request,status), goto end );
    MPI_CHECK( MPI_Allreduce(&loc_update,&glob_update,1,MPI_INT8_T,MPI_LOR,
                             comm),goto end);

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
  qsort(proclists,nproclists,sizeof(PMMG_cellLnkdList*),PMMG_compare_cellLnkdList);

  /* Remove the non unique paths */
  if ( nproclists ) {
    idx = 0;
    for ( i=1; i<nproclists; ++i ) {
      assert ( proclists[i]->nitem );
      if ( PMMG_compare_cellLnkdList(&proclists[i],&proclists[idx]) ) {
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
          if ( nitem_ext_comm[val1_j] == ext_node_comm->nitem || !ext_node_comm->nitem ) {
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
          if ( nitem_ext_comm[val1_i] == ext_node_comm->nitem || !ext_node_comm->nitem ) {
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
               next_comm,parmesh->next_node_comm,PMMG_Ext_comm,
               "list of external communicator",goto end);
  parmesh->next_node_comm = next_comm;

  /* Success */
  ier = 1;

end:
  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,int,"node communicator");

  if ( proclists ) {
    for ( k=0; k<nproclists; ++k ) {
      PMMG_DEL_MEM(parmesh,proclists[k]->item,PMMG_lnkdCell,"linked list array");
      PMMG_DEL_MEM(parmesh,proclists[k],PMMG_lnkdList,"linked list pointer");
    }
    PMMG_DEL_MEM(parmesh,proclists,PMMG_lnkdList*,"array of linked lists");
  }
  PMMG_DEL_MEM(parmesh,comm_ptr,PMMG_pExt_comm,
              "array of pointers toward the external communicators");

  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];

    // Change this and add to the external comm the possibility to not
    // unalloc/realloc every time, thus, here, we will be able to reset the
    // communicators without unallocated it
    PMMG_DEL_MEM(parmesh,ext_node_comm->itosend,int,"i2send");
    PMMG_DEL_MEM(parmesh,ext_node_comm->itorecv,int,"i2recv");
  }

  PMMG_DEL_MEM(parmesh,request,MPI_Request,"mpi request array");
  PMMG_DEL_MEM(parmesh,status,MPI_Status,"mpi status array");
  PMMG_DEL_MEM(parmesh,i2send_size,int,"size of the i2send array");
  PMMG_DEL_MEM(parmesh,i2recv_size,int,"size of the i2recv array");

  PMMG_DEL_MEM(parmesh,nitem_ext_comm,int,
               "number of items in each external communicator");

  PMMG_DEL_MEM(parmesh,list.item,PMMG_lnkdCell,"linked list array");

  return ier;
}

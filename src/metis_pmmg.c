/**
 * \file metis_pmmg.c
 * \brief Partition mesh using metis
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "metis_pmmg.h"
#include "linkedlist_pmmg.h"


/**
 * \param parmesh pointer toward the parmesh structure.
 * \param hgrp pointer toward the hatable.
 * \param hsiz initial size of hash table.
 * \param hmax maximal size of the hash table.
 *
 * \return 1 if success, 0 if fail.
 *
 * Initialisation of the hash table of group adjacency.
 *
 */
static inline
int PMMG_hashNew( PMMG_pParMesh parmesh,PMMG_HGrp *hash,int hsiz,int hmax ) {
  int k;

  /* adjust hash table params */
  hash->siz  = hsiz+1;
  hash->max  = hmax + 2;
  hash->nxt  = hash->siz;

  PMMG_CALLOC(parmesh,hash->item,hash->max+1,PMMG_hgrp,"group hash table",return 0);

  for (k=0; k<hash->siz; ++k )
    hash->item[k].adj = PMMG_UNSET;

  for (k=hash->siz; k<hash->max; ++k) {
    hash->item[k].adj = PMMG_UNSET;
    hash->item[k].nxt = k+1;
  }

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param hash pointer toward the hash table of group adja
 * \param k    local id (in the parmesh) of the group to which we want to add
 *             an adjacent.
 * \param adj  global id of the adjacent group
 *
 * \return 2 if adj is already listed, 1 if we success to add \a adj to the hash
 * table, 0 if fail
 *
 * Add the group id \a adj to the hash table at key \a k+1. The group is added
 * in such way that the key \a k+1 gives a list of groups sorted in ascending
 * order (except for the first item of the list that is not ordered).
 *
 */
static inline
int PMMG_hashGrp( PMMG_pParMesh parmesh,PMMG_HGrp *hash, int k, idx_t adj ) {
  PMMG_hgrp  *ph;
  int        tmp_nxt,j,newsize;

  ph  = &hash->item[k+1];

  if ( ph->adj == adj ) return 2;

  tmp_nxt = 0;
  if ( PMMG_UNSET != ph->adj ) {
    while ( ph->nxt && (ph->nxt<hash->max) && (hash->item[ph->nxt].adj<adj) )
      ph = &hash->item[ph->nxt];

    if ( hash->item[ph->nxt].adj == adj ) return 2;

    tmp_nxt   = ph->nxt;
    ph->nxt   = hash->nxt;
    ph        = &hash->item[hash->nxt];

    if ( hash->nxt >= hash->max-1 ) {
      newsize = MG_MAX((int)((1+PMMG_GAP)*hash->max),hash->max+1);
      PMMG_RECALLOC(parmesh,hash->item,newsize,hash->max,PMMG_hgrp,
                    "grp hash table: if we pass too much time here"
                    " (on true cases (more than 20 procs)), it means"
                    " that the hmax value is not a good estimation of the"
                    " interprocessor adjacency... TO CHANGE",return 0);
      hash->max = newsize;

      /* ph pointer may be false after realloc */
      ph        = &hash->item[hash->nxt];

      for (j=ph->nxt; j<hash->max; j++) {
        hash->item[j].adj = PMMG_UNSET;
        hash->item[j].nxt = j+1;
      }
    }
    hash->nxt = ph->nxt;
  }

  /* insert new group */
  ph->adj = adj;
  ph->nxt = tmp_nxt;

  return 1;
}



int PMMG_correct_meshElts2metis( PMMG_pParMesh parmesh,idx_t* part,idx_t ne,idx_t nproc ) {
  PMMG_lnkdList **partlist;
  idx_t iproc,ie,dummy;
  int nempt,iempt;

  /* Initialize lists */
  PMMG_CALLOC(parmesh,partlist,nproc,PMMG_lnkdList*,"array of list pointers",return 0);
  for( iproc=0; iproc<nproc; iproc++ ) {
    PMMG_CALLOC(parmesh,partlist[iproc],1,PMMG_lnkdList,"linked list pointer",return 0);
    if( !PMMG_lnkdListNew(parmesh,partlist[iproc],iproc,PMMG_LISTSIZE) ) return 0;
  }

  /* Fill the lists */
  for( ie=0; ie<ne; ie++ ) {
    iproc = part[ie];
    if( !PMMG_add_cell2lnkdList(parmesh,partlist[iproc],ie,iproc) ) return 0;
  }

  /* Sort lists based on nb. of entities, in ascending order */
  qsort(partlist,nproc,sizeof(PMMG_lnkdList*),PMMG_compare_lnkdList);

  /* Count empty partitions */
  nempt = 0;
  for( iproc=0; iproc<nproc; iproc++ )
    if( !partlist[iproc]->nitem ) nempt++;
  assert( nempt < nproc );
  if( !nempt ) return 1;

  /** Correct partitioning */
  iempt = 0;
  iproc = nproc-1;
  while( nempt ) {
    /* Get next "reservoir" proc */
    while( partlist[iproc]->nitem < partlist[iproc-1]->nitem )
      iproc--;
    /* Pop entity ie from iproc, add to iempt */
    if( !PMMG_pop_cell_lnkdList(parmesh,partlist[iproc],&ie,&dummy) ) return 0;
    if( !PMMG_add_cell2lnkdList(parmesh,partlist[iempt],ie,iempt) ) return 0;
    /* Update partition table and go on to next empty proc */
    part[ie] = partlist[iempt]->id;
    iempt++;
    nempt--;
  }

  /* Deallocate lists */
  for( iproc=0; iproc<nproc; iproc++ )
    PMMG_DEL_MEM(parmesh,partlist[iproc]->item,PMMG_lnkdCell,"linked list");
  PMMG_DEL_MEM(parmesh,partlist,PMMG_lnkdList*,"array of linked lists");

  return 1;
}

/**
 * \param parmesh pointer toward the PMMG parmesh structure
 * \param mesh pointer toward a MMG5 mesh structure
 * \param xadj pointer toward the position of the elt adjacents in adjncy
 * \param adjncy pointer toward the list of the adjacent of each elt
 * \param nadjncy number of data in adjncy array
 * \param memAv pointer toward the available memory (to update)
 *
 * \return  1 if success, 0 if fail
 *
 * Build the metis graph with the mesh elements as metis nodes.
 *
 * \warning the mesh must be packed
 *
 */
int PMMG_graph_meshElts2metis( PMMG_pParMesh parmesh,MMG5_pMesh mesh,
                               idx_t **xadj,idx_t **adjncy,
                               idx_t *nadjncy,size_t *memAv) {
  size_t     memMaxOld;
  int        *adja;
  int        j,k,iadr,jel,count,nbAdj,ier;

  /** Step 1: mesh adjacency creation */

  /* Give the available memory to the mesh */
  memMaxOld     = mesh->memMax;
  mesh->memMax += *memAv;
  if ( (!mesh->adja) && (1 != MMG3D_hashTetra( mesh, 1 )) ) {
    fprintf( stderr,"  ## PMMG Hashing problem (1).\n" );
    return 0;
  }
  /* Update the available memory */
  mesh->memMax = mesh->memCur;
  *memAv      -= (mesh->memMax - memMaxOld);

  /** Step 2: build the metis graph */

  /* Give the available memory to the parmesh */
  memMaxOld        = parmesh->memMax;
  parmesh->memMax += *memAv;

  PMMG_CALLOC(parmesh, (*xadj), mesh->ne+1, idx_t, "allocate xadj",
              return 0);

  /** 1) Count the number of adjacent of each elements and fill xadj */
  (*xadj)[0] = 0;
  (*nadjncy) = 0;
  for( k = 1; k <= mesh->ne; k++ ) {
    nbAdj = 0;
    iadr = 4*(k-1) + 1;
    adja = &mesh->adja[iadr];
    for( j = 0; j < 4; j++ )
      if( adja[j] )
        nbAdj++;

    (*nadjncy)+= nbAdj;
    (*xadj)[k] = (*nadjncy);
  }

  /** 2) List the adjacent of each elts in adjncy */
  ier = 1;
  ++(*nadjncy);
  PMMG_CALLOC(parmesh, (*adjncy), (*nadjncy), idx_t, "allocate adjncy",
              ier = 0;);

  if ( ier ) {
    count = 0;
    for( k = 1; k <= mesh->ne; k++ ) {
      iadr = 4*(k-1) + 1;
      adja = &mesh->adja[iadr];
      for ( j = 0; j < 4; j++ ) {
        jel = adja[j] / 4;

        if ( !jel )
          continue;

        (*adjncy)[count++] = jel-1;
      }
      assert( count == ( (*xadj)[k] ) );
    }
  }
  else {
    PMMG_DEL_MEM(parmesh, xadj, idx_t, "deallocate xadj" );
  }

  parmesh->memMax = parmesh->memCur;
  *memAv -= (parmesh->memMax - memMaxOld);

  return ier;
}

/**
 * \param parmesh pointer toward the PMMG parmesh structure
 * \param vtxdist pointer toward the description of the node distribution
 * \param xadj pointer toward the position of the elt adjacents in adjncy
 * \param adjncy pointer toward the list of the adjacent of each elt
 * \param nadjncy number of data in adjncy array
 * \param vwgt pointer toward the metis node weights
 * \param wgtflag how to apply the metis weights
 * \param numflag numbering style (C versus frotran)
 * \param ncon number of of weights per metis node
 * \param nproc number of partitions asked
 * \param tpwgt pointer toward the fraction of weight to send to each domain
 * \param ubvec imbalance tolerance for each vertex weight
 *
 * \return  1 if success, 0 if fail
 *
 * Build the metis graph with the mesh elements as metis nodes.
 *
 */
int PMMG_graph_parmeshGrps2parmetis( PMMG_pParMesh parmesh,idx_t **vtxdist,
                                     idx_t **xadj,idx_t **adjncy,idx_t *nadjncy,
                                     idx_t **vwgt,idx_t *wgtflag,idx_t *numflag,
                                     idx_t *ncon,idx_t nproc,real_t **tpwgts,
                                     real_t **ubvec) {
  PMMG_pGrp      grp;
  PMMG_pExt_comm ext_face_comm;
  PMMG_pInt_comm int_face_comm;
  PMMG_HGrp      hash;
  PMMG_hgrp      *ph;
  MMG5_pMesh     mesh;
  MMG5_pTetra    pt;
  MPI_Comm       comm;
  MPI_Status     status;
  int            *face2int_face_comm_index2;
  int            *intvalues,*itosend,*itorecv;
  int            found,color;
  int            ngrp,myrank,nitem,k,igrp,igrp_adj,i,idx;

  *wgtflag = 2; /* Weights applies on vertices */
  *numflag = 0; /* C-style numbering */
  *ncon    = 1; /* number of weight per metis node */

  comm   = parmesh->comm;
  grp    = parmesh->listgrp;
  myrank = parmesh->myrank;
  ngrp   = parmesh->ngrp;

  /** Step 1: Fill vtxdist array with the range of groups local to each
   * processor */
  PMMG_CALLOC(parmesh,*vtxdist,nproc+1,idx_t,"parmetis vtxdist", return 0);

  MPI_CHECK( MPI_Allgather(&ngrp,1,MPI_INT,&(*vtxdist)[1],1,MPI_INT,comm),
             goto fail_1 );

  for ( k=1; k<=nproc; ++k )
    (*vtxdist)[k] += (*vtxdist)[k-1];

  /** Step 2: Fill weights array with the number of MG_PARBDY face per group */
  PMMG_CALLOC(parmesh,*vwgt,ngrp,idx_t,"parmetis vwgt", goto fail_1);

  for ( igrp=0; igrp<ngrp; ++igrp ) {
    mesh = parmesh->listgrp[igrp].mesh;

    if ( !mesh ) {
      (*vwgt)[igrp] = 1;
      continue;
    }

    for ( k=1; k<=mesh->ne; ++k ) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) ) continue;

      (*vwgt)[igrp] += pt->mark;
    }
  }

  /* Fill tpwgts */
  PMMG_CALLOC(parmesh,*tpwgts,(*ncon)*nproc,real_t,"parmetis tpwgts", goto fail_2);
  for ( k=0; k < (*ncon)*nproc ; ++k )
    (*tpwgts)[k] = 1./(double)nproc;

  /* Fill ubvec */
  PMMG_CALLOC(parmesh,*ubvec,(*ncon),real_t,"parmetis ubvec", goto fail_3);
  for ( k=0; k < (*ncon); ++k )
    (*ubvec)[k] = PMMG_UBVEC_DEF;

  /** Step 3: Fill the internal communicator with the greater index of the 2
   * groups to which the face belong */
  PMMG_CALLOC(parmesh,*xadj,ngrp+1,idx_t,"parmetis xadj", goto fail_4);

  int_face_comm = parmesh->int_face_comm;

  PMMG_MALLOC(parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,
              "face communicator",goto fail_5);

  /* Face communicator initialization */
  intvalues = parmesh->int_face_comm->intvalues;
  for ( k=0; k < int_face_comm->nitem; ++k )
    intvalues[k] = PMMG_UNSET;

  for ( igrp=ngrp-1; igrp>=0; --igrp ) {
    grp                       = &parmesh->listgrp[igrp];
    face2int_face_comm_index2 = grp->face2int_face_comm_index2;

    for ( k=0; k<grp->nitem_int_face_comm; ++k )
      if ( PMMG_UNSET == intvalues[face2int_face_comm_index2[k] ] )
        intvalues[face2int_face_comm_index2[k]]= igrp;
  }

  /** Step 4: Send and receive external communicators filled by the group id of
   * the neighbours (through the faces) */
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    nitem         = ext_face_comm->nitem;
    color         = ext_face_comm->color_out;

    PMMG_CALLOC(parmesh,ext_face_comm->itosend,nitem,int,"itosend array",
                goto fail_6);
    itosend = ext_face_comm->itosend;

    PMMG_CALLOC(parmesh,ext_face_comm->itorecv,nitem,int,"itorecv array",
                goto fail_6);
    itorecv       = ext_face_comm->itorecv;

    for ( i=0; i<nitem; ++i ) {
      idx            = ext_face_comm->int_comm_index[i];
      itosend[i]     = intvalues[idx] ;
      /* Mark the face as boundary in the intvalues array */
      intvalues[idx] = PMMG_UNSET;
    }

    MPI_CHECK(
      MPI_Sendrecv(itosend,nitem,MPI_INT,color,MPI_PARMESHGRPS2PARMETIS_TAG,
                   itorecv,nitem,MPI_INT,color,MPI_PARMESHGRPS2PARMETIS_TAG,
                   comm,&status),goto fail_6 );
  }

   /** Step 5: Process the external communicators to count for each group the
    * adjacent groups located on another processor and fill the sorted linked
    * list of adjacency */

  /* hash is used to store the sorted list of adjacent groups to a group */
  if ( !PMMG_hashNew(parmesh,&hash,ngrp+1,PMMG_NBADJA_GRPS*ngrp+1) )
    goto fail_6;

  for (  k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    itosend       = ext_face_comm->itosend;
    itorecv       = ext_face_comm->itorecv;
    nitem         = ext_face_comm->nitem;

    /* i2send array contains the group id of the boundary faces and i2recv the
     * group id of the same face in the other proc */
    for ( i=0; i<nitem; ++i ) {
      /* Get the group id of the face in our proc and the group id of the face
         in the adjacent proc */
      igrp     = itosend[i];
      igrp_adj = itorecv[i];

      assert ( igrp != PMMG_UNSET );
      assert ( igrp_adj != PMMG_UNSET );

      /* Search the neighbour in the hash table and insert it if not found */
      found = PMMG_hashGrp(parmesh,&hash,igrp,igrp_adj
                           +(*vtxdist)[ext_face_comm->color_out]);

      if ( !found ) {
        fprintf(stderr,"  ## Error: %s: unable to add a new group in adjacency"
                " hash table.\n",__func__);
        goto fail_7;
      }
      if ( found==2 ) continue; // The group is already in the adja list

      ++(*xadj)[ igrp+1 ];
    }
  }

  /** Step 6: Foreach group, process the internal communicator arrays and count
   * the number of neighbours of the group. Update at the same time the linked
   * list of adjacency */
  for ( igrp=0; igrp<ngrp; ++igrp ) {
    grp                       = &parmesh->listgrp[igrp];
    face2int_face_comm_index2 = grp->face2int_face_comm_index2;

    for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
      igrp_adj = intvalues[face2int_face_comm_index2[i]];

      if ( igrp_adj==-1 || igrp_adj<=igrp ) continue;

      /* Search the neighbour in the hash table and insert it if not found */
       found = PMMG_hashGrp(parmesh,&hash,igrp,igrp_adj+(*vtxdist)[myrank] );

      if ( !found ) {
        fprintf(stderr,"  ## Error: %s: unable to add a new group in adjacency"
                " hash table.\n",__func__);
        goto fail_7;
      }

      if ( found==2 ) continue; // The group is already in the adja list

      ++(*xadj)[ igrp+1 ];

      /* add igrp to the sorted list of adjacnts to igrp_adj */
      found = PMMG_hashGrp(parmesh,&hash,igrp_adj,igrp+(*vtxdist)[myrank] );
      if ( 1!=found ) {
        fprintf(stderr,"  ## Error: %s: unable to add a new group in adjacency"
                " hash table.\n",__func__);
        goto fail_7;
      }
      ++(*xadj)[ igrp_adj+1 ];
    }
  }

  /** Step 7: xadj array contains the number of adja per group, fill it for
   * Metis (it must contains the index of the first adja of the group in the
   * array adjcncy that list all the group adja in a 1D array) */
  for ( igrp=1; igrp <= parmesh->ngrp; ++igrp )
    (*xadj)[igrp] += (*xadj)[igrp-1];

  /** Step 8: Fill adjncy array at metis format */
  PMMG_CALLOC(parmesh,*adjncy,(*xadj)[ngrp],idx_t,"adjcncy parmetis array",
              goto fail_7);

  (*nadjncy) = 0;
  for ( igrp=0; igrp<=ngrp; ++igrp ) {

    ph = &hash.item[igrp+1];
    if ( PMMG_UNSET==ph->adj ) continue;

    (*adjncy)[(*nadjncy)++] = ph->adj;

    while ( ph->nxt ) {
      ph                = &hash.item[ph->nxt];
      (*adjncy)[(*nadjncy)++] = ph->adj;
    }
  }
  assert ( (*nadjncy)==(*xadj)[ngrp] );

  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    nitem = ext_face_comm->nitem;

    if ( ext_face_comm->itorecv )
      PMMG_DEL_MEM(parmesh,ext_face_comm->itorecv,int,"itorecv array");
    if ( ext_face_comm->itosend )
      PMMG_DEL_MEM(parmesh,ext_face_comm->itosend,int,"itosend array");
  }
  PMMG_DEL_MEM(parmesh,hash.item,PMMG_hgrp,"group hash table");
  PMMG_DEL_MEM(parmesh,int_face_comm->intvalues,int,"face communicator");

  return 1;

fail_7:
  if ( hash.item )
    PMMG_DEL_MEM(parmesh,hash.item,PMMG_hgrp,"group hash table");
fail_6:
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    if ( ext_face_comm->itorecv )
      PMMG_DEL_MEM(parmesh,ext_face_comm->itorecv,int,"itorecv array");
  }
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    if ( ext_face_comm->itosend )
      PMMG_DEL_MEM(parmesh,ext_face_comm->itosend,int,"itosend array");
  }
  if ( int_face_comm->intvalues )
    PMMG_DEL_MEM(parmesh,int_face_comm->intvalues,int,"face communicator");
fail_5:
  PMMG_DEL_MEM(parmesh,*xadj,idx_t,"parmetis xadj");
fail_4:
  PMMG_DEL_MEM(parmesh,*ubvec,real_t,"parmetis ubvec");
fail_3:
  PMMG_DEL_MEM(parmesh,*tpwgts,real_t,"parmetis tpwgts");
fail_2:
  PMMG_DEL_MEM(parmesh,*vwgt,idx_t,"parmetis vwgt");
fail_1:
  PMMG_DEL_MEM(parmesh,*vtxdist,idx_t,"parmetis vtxdist");

  return 0;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param part pointer of an array containing the partitions (at the end)
 * \param nproc number of partitions asked
 *
 * \return  1 if success, 0 if fail
 *
 * Use metis to partition the first mesh in the list of meshes into nproc groups
 *
 */
int PMMG_part_meshElts2metis( PMMG_pParMesh parmesh, idx_t* part, idx_t nproc )
{
  PMMG_pGrp  grp = parmesh->listgrp;
  MMG5_pMesh mesh = grp[0].mesh;
  size_t     memAv;
  idx_t      *xadj = NULL;
  idx_t      *adjncy = NULL;
  idx_t      nelt = mesh->ne;
  idx_t      ncon = 1; // number of balancing constraint
  idx_t      objval = 0;
  int        adjsize;
  int        ier = 0;
  int        status = 1;

  /* Fit the parmesh and the meshes in memory and compute the available memory */
  parmesh->memMax = parmesh->memCur;
  parmesh->listgrp[0].mesh->memMax = parmesh->listgrp[0].mesh->memCur;
  memAv = parmesh->memGloMax-parmesh->memMax-parmesh->listgrp[0].mesh->memMax;

  /** Build the graph */
  if ( !PMMG_graph_meshElts2metis(parmesh,mesh,&xadj,&adjncy,&adjsize,&memAv) )
    return 0;

  /* Give the memory to the parmesh */
  parmesh->memMax += memAv;

  /** Call metis and get the partition array */
  ier = METIS_PartGraphKway( &nelt,&ncon,xadj,adjncy,NULL,NULL,NULL,&nproc,
                             NULL,NULL,NULL,&objval, part );
  if ( ier != METIS_OK ) {
    switch ( ier ) {
      case METIS_ERROR_INPUT:
        fprintf(stderr, "METIS_ERROR_INPUT: input data error\n" );
        break;
      case METIS_ERROR_MEMORY:
        fprintf(stderr, "METIS_ERROR_MEMORY: could not allocate memory error\n" );
        break;
      case METIS_ERROR:
        fprintf(stderr, "METIS_ERROR: generic error\n" );
        break;
      default:
        fprintf(stderr, "METIS_ERROR: update your METIS error handling\n" );
        break;
    }
    status = 0;
  }

  /** Correct partitioning to avoid empty partitions */
  if( !PMMG_correct_meshElts2metis( parmesh,part,nelt,nproc ) ) return 0;

  PMMG_DEL_MEM(parmesh, adjncy, idx_t, "deallocate adjncy" );
  PMMG_DEL_MEM(parmesh, xadj, idx_t, "deallocate xadj" );

  return status;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param part pointer of an array containing the partitions (at the end)
 * \param nproc number of partitions asked
 *
 * \return  1 if success, 0 if fail
 *
 * Use metis to partition the first mesh in the list of meshes into nproc groups
 *
 */
int PMMG_part_parmeshGrps2parmetis( PMMG_pParMesh parmesh,idx_t* part,idx_t nproc )
{
  real_t     *tpwgts,*ubvec;
  idx_t      *xadj,*adjncy,*vwgt,*vtxdist,adjsize,edgecut;
  idx_t      wgtflag,numflag,ncon,options[3];
  int        ngrp,nprocs,ier;

  ngrp   = parmesh->ngrp;
  nprocs = parmesh->nprocs;
  ier    = 1;

  /** Build the parmetis graph */
  xadj   = adjncy = vwgt = vtxdist = NULL;
  tpwgts = ubvec  =  NULL;
  options[0] = 0;

  if ( !PMMG_graph_parmeshGrps2parmetis(parmesh,&vtxdist,&xadj,&adjncy,&adjsize,
                                        &vwgt,&wgtflag,&numflag,&ncon,
                                        nproc,&tpwgts,&ubvec) ) {
    fprintf(stderr,"\n  ## Error: Unable to build parmetis graph.\n");
    return 0;
  }

  /** Call parmetis and get the partition array */
  if ( 2 < nprocs + ngrp ) {
    if ( ParMETIS_V3_PartKway( vtxdist,xadj,adjncy,vwgt,NULL,&wgtflag,&numflag,
                               &ncon,&nproc,tpwgts,ubvec,options,&edgecut,part,
                               &parmesh->comm) != METIS_OK ) {
        fprintf(stderr,"\n  ## Error: Parmetis fails.\n" );
        ier = 0;
    }
  }

  PMMG_DEL_MEM(parmesh, adjncy, idx_t, "deallocate adjncy" );
  PMMG_DEL_MEM(parmesh, xadj, idx_t, "deallocate xadj" );
  PMMG_DEL_MEM(parmesh, ubvec, real_t,"parmetis ubvec");
  PMMG_DEL_MEM(parmesh, vwgt, idx_t, "deallocate vwgt" );
  PMMG_DEL_MEM(parmesh, tpwgts, real_t, "deallocate tpwgts" );
  PMMG_DEL_MEM(parmesh, vtxdist, idx_t, "deallocate vtxdist" );

  return ier;
}

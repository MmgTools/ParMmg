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
 * \file distributegrps_pmmg.c
 * \brief Group distribution on the processors
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \author Luca Cirrotola (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 */
#include <mpi.h>
#include "parmmg.h"
#include "linkedlist_pmmg.h"

/**
 * \param parmesh pointer toward a parmesh structure
 * \param interactions pointer toward the interaction list (to fill)
 * \param interaction_map map of interactions (to fill)
 *
 * \return -1 if fail, the number of interactions otherwise (size of \a
 * interactions array)
 *
 * Compute the map of interactions between processors and count the number of
 * interactions. \a interaction_map[k] is 1 if I have interactions with proc \a
 * k. \a interactions[2k:2k+1] contains the index of the \f$k^th \f$ pair of
 * procs that interacts. It is shuffled to maximize comm overlapping.
 *
 */
int PMMG_interactionMap(PMMG_pParMesh parmesh,int **interactions,int **interaction_map) {

  PMMG_pGrp      grp;
  PMMG_pExt_comm ext_comm;
  MPI_Comm       comm;
  MPI_Status     status;
  MPI_Request    *request;

  int           *interact_list;
  int           *recv_count,*loc_recv_list,*ninter_all;
  int           *dist_list,*interaction_shuffling,*interaction_idx;
  int           buf[2],mybuf[2],bufs[2];
  int           *neighbours,*receivers;
  int           nneighs_max,nrecv_max,idx,mysndr;
  int           nprocs,myrank,i,j,k,l,proc,ninteractions;
  int           ier;

  if ( parmesh->info.imprim >= PMMG_VERB_DETQUAL ) {
    fprintf(stdout,"\n  ** MAP OF INTERACTIONS COMPUTATION\n");
  }

  ier = 1;
  myrank = parmesh->myrank;
  comm   = parmesh->comm;
  nprocs = parmesh->nprocs;

  /** Step 1: list of receivers of myrank */
  loc_recv_list = NULL;
  PMMG_CALLOC ( parmesh,loc_recv_list,nprocs,int,"loc_recv_list",ier=0 );

  if ( ier ) {
    for ( k=0; k<parmesh->ngrp; ++k ) {
      grp = parmesh->listgrp + k;
      if ( grp->flag != myrank ) {
        ++loc_recv_list[grp->flag];
      }
    }
  }

  /** Step 2: compute number of neighbours and receivers per proc */
  mybuf[0] = mybuf[1] = buf[0] = buf[1] = 0;

  /* number of neighbours */
  for ( i=0; i< parmesh->next_face_comm; ++i ) {
    ext_comm = &parmesh->ext_face_comm[i];
    if ( ext_comm->nitem ) {
      assert ( ext_comm->color_out != myrank );
      ++mybuf[0];
    }
  }

  /* number of receivers per proc */
  if ( ier ) {
    for ( k=0; k<nprocs; ++k ) {
      if ( loc_recv_list[k] ) {
        ++mybuf[1];
      }
    }
  }

  MPI_CHECK( MPI_Allreduce(mybuf,buf ,2,MPI_INT,MPI_MAX,parmesh->comm),ier=0 );
  MPI_CHECK( MPI_Allreduce(mybuf,bufs,2,MPI_INT,MPI_SUM,parmesh->comm),ier=0 );

  nneighs_max   = buf[0];
  nrecv_max     = buf[1];
  ninteractions = bufs[1];

  /** Step 3: compute interactions */
  PMMG_CALLOC( parmesh,(*interactions),2*ninteractions,int,"interactions",ier = 0 );

  idx = 0;
  for ( k=0; k<nprocs; ++k ) {
    if ( loc_recv_list[k] ) {
      ++idx;
    }
  }

  /* Compute the number of interactions from rank 0 to \a myrank (because each rank
   * will write in a different area of the \a interactions array) */
  ninter_all = NULL;
  PMMG_CALLOC ( parmesh,ninter_all,nprocs,int,"ninter_all",ier=0 );
  MPI_CHECK( MPI_Allgather(&idx,1,MPI_INT,ninter_all,1,MPI_INT,comm),ier=0 );

  idx = 0;
  for ( k=0; k<myrank; ++k ) {
    idx += ninter_all[k];
  }
  PMMG_DEL_MEM ( parmesh,ninter_all,int,"ninter_all" );

  for ( k=0; k<nprocs; ++k ) {
    if ( loc_recv_list[k] ) {
      (*interactions)[   idx*2] = myrank;
      (*interactions)[1+ idx*2] = k;
      ++idx;
    }
  }

  /* Now everone knows the entire list of interactions of the others */
  MPI_CHECK(MPI_Allreduce(MPI_IN_PLACE,(*interactions),2*ninteractions,MPI_INT,MPI_MAX,comm),ier=0);

  /** Step 4: compute list of proc neighbors */
  neighbours    = NULL;
  interact_list = NULL;
  PMMG_MALLOC ( parmesh,neighbours,   nneighs_max*nprocs,int,"neighbours",ier=0 );
  PMMG_CALLOC ( parmesh,interact_list,nneighs_max       ,int,"interact_list",ier=0 );

  idx = 0;
  for ( i=0; i<parmesh->next_face_comm; ++i ) {
    ext_comm = &parmesh->ext_face_comm[i];
    if ( ext_comm->nitem ) {
      assert ( ext_comm->color_out != myrank );
      /* Store color_out + 1 to distinguish uninitialized values from 0 rank */
      interact_list[idx++] = ext_comm->color_out+1 ;
    }
  }

  /* Now everone knows the entire list of neighbours of the others (a
   * allgatherv is possible too here) */
  MPI_CHECK(
    MPI_Allgather(interact_list,nneighs_max,MPI_INT,neighbours,nneighs_max,MPI_INT,comm)
    ,ier=0);

  /** Step 5: get list of recv for each proc */
  receivers = NULL;
  PMMG_MALLOC   ( parmesh,receivers,    nrecv_max*nprocs,int,"receivers"   ,ier=0 );
  PMMG_REALLOC  ( parmesh,interact_list,nrecv_max,nneighs_max,int,"interact_list" ,ier=0 );
  memset ( interact_list,0, nrecv_max * sizeof(int) );

  idx = 0;
  for ( k=0; k<nprocs; ++k ) {
    if (loc_recv_list[k] ) {
      /* Store proc index + 1 to distinguish uninitialized values from 0 rank */
      interact_list[idx++] = k+1;
    }
  }

  /* Now everone knows the entire list of receivers of the others (a
   * allgatherv is possible too here) */
  MPI_CHECK(MPI_Allgather(interact_list,nrecv_max,MPI_INT,receivers,nrecv_max,MPI_INT,comm),ier=0);
  PMMG_DEL_MEM ( parmesh, interact_list,int,"interact_list");

  /** Step 6: flag all the interacting procs */
  PMMG_CALLOC ( parmesh,(*interaction_map),nprocs,int,"interaction_map" ,ier=0 );

  recv_count = NULL;
  PMMG_MALLOC ( parmesh,recv_count,nprocs,int,"recv_count" ,ier=0 );

  for ( k=0; k<nprocs; ++k ) {
    if ( loc_recv_list[k] ) {
      /* flag the receiver */
      (*interaction_map)[k] = 1;
    }
  }

  /* flag my neighbours and the procs where they may be transferred */
  for ( k=0; k<nneighs_max; ++k ) {
    idx = k + myrank*nneighs_max;
    if ( neighbours[idx] ) {
      proc = neighbours[idx]-1;

      /* flag the neighbour */
      (*interaction_map)[proc] = 1;
      for ( i=0; i<nrecv_max; ++i ) {
        j = i + proc*nrecv_max;
        if ( receivers[j] ) {
          (*interaction_map)[receivers[j]-1] = 1;
        }
      }
    }
  }

  /* flag the other processors that will receive data from the receiver */
  for ( k=0; k<nprocs; ++k ) {
    for ( i=0; i<nrecv_max; ++i ) {
      j = i + k*nrecv_max;
      if ( receivers[j] ) {
        for ( proc=0; proc<nprocs; ++proc ) {
          for ( l=0; l<nrecv_max; ++l ) {
            idx = l + proc*nrecv_max;
            if ( receivers[idx] == myrank+1 ) {
              (*interaction_map)[proc] = 1;
            }
          }
        }
      }
    }
  }

  /** Step 7: count the procs exchanging data with me */
  MPI_Allreduce((*interaction_map),recv_count,nprocs,MPI_INTEGER,MPI_SUM,comm);

  /** Step 8: build interaction map based on group exchange */
  request = NULL;
  PMMG_MALLOC ( parmesh, request, nprocs, MPI_Request, "request_list", ier = 0);
  for ( k=0; k<nprocs; ++k ) {
    request[k] = MPI_REQUEST_NULL;
    if ( (*interaction_map)[k] ) {
      MPI_CHECK(
        MPI_Isend(&myrank,1,MPI_INT,k,MPI_TRANSFER_GRP_TAG,comm,&request[k]),
        ier=0);
    }
  }

  for ( k=0; k<recv_count[myrank]; ++k ) {
    MPI_CHECK (
      MPI_Recv ( &mysndr, 1, MPI_INT,MPI_ANY_SOURCE,MPI_TRANSFER_GRP_TAG,comm,&status),
      ier=0);
    (*interaction_map)[mysndr] = 1;
  }
  MPI_CHECK ( MPI_Waitall(nprocs,request,MPI_STATUSES_IGNORE), ier=0 );
  PMMG_DEL_MEM ( parmesh, request,MPI_Request,"request_list");

  (*interaction_map)[myrank] = 1;

  /** Step 9: compute ratio of active versus total procs */
  if  ( parmesh->info.imprim0 > PMMG_VERB_DETQUAL ) {

    mybuf[0] = 0;
    for ( k=0; k<nprocs; ++k ) {
      mybuf[0] += (*interaction_map)[k];
    }
    mybuf[1] = nprocs;

    MPI_CHECK ( MPI_Allreduce(mybuf,buf,2,MPI_INT,MPI_SUM,comm), ier=0 );

    if  ( parmesh->info.imprim > PMMG_VERB_DETQUAL ) {
      fprintf(stdout,"      ratio of the interaction map: %12.4f\n",
              (float)buf[0]/(float)buf[1]);
    }
  }

  PMMG_DEL_MEM ( parmesh, loc_recv_list,int,"loc_recv_list");
  PMMG_DEL_MEM ( parmesh, neighbours,int,"neighbours");
  PMMG_DEL_MEM ( parmesh, receivers ,int,"receivers");

  PMMG_DEL_MEM ( parmesh, recv_count,int,"rcv_count");

  /** Step 10: Shuffling to maximize overlapping */
  interaction_shuffling = interaction_idx = NULL;
  PMMG_MALLOC  ( parmesh, interaction_shuffling, ninteractions,int,"interaction_shuffling" ,ier=0 );
  PMMG_MALLOC  ( parmesh, interaction_idx,       ninteractions,int,"interaction_idx" ,ier=0 );

  /* Generation of a random vector */
  if ( myrank == parmesh->info.root ) {
    for ( k=0; k<ninteractions; ++k) {
      interaction_shuffling[k]=rand()%(2*nprocs);
    }
  }
  MPI_CHECK(MPI_Bcast(interaction_shuffling,ninteractions,MPI_INT,parmesh->info.root,comm),ier=0);

  /* Sort the interaction */
  PMMG_sort_iarray(parmesh,NULL,interaction_shuffling,interaction_idx,ninteractions);

  /* change the order of the interactions */
  for ( k=0; k<ninteractions; ++k ) {
    interaction_shuffling[k] = (*interactions)[2*interaction_idx[k]];
  }
  for ( k=0; k<ninteractions; ++k) {
    (*interactions)[2*k] = interaction_shuffling[k];
  }

  for ( k=0; k<ninteractions; ++k) {
    interaction_shuffling[k] = (*interactions)[2*interaction_idx[k]+1];
  }
  for ( k=0; k<ninteractions; ++k) {
    (*interactions)[2*k+1] = interaction_shuffling[k];
  }
  PMMG_DEL_MEM ( parmesh, interaction_shuffling, int, "interaction_shuffling");
  PMMG_DEL_MEM ( parmesh, interaction_idx      , int, "interaction_idx"      );

  /* Compute the distance between the interactions */
  if ( parmesh->ddebug ) {
    dist_list = NULL;
    PMMG_MALLOC  ( parmesh, dist_list, nprocs,int,"dist_list" ,ier=0 );

    int last_interaction = 0;
    int dist = 2*nprocs;

    for ( k=0; k<ninteractions; ++k ) {
      j = (*interactions)[2*k];
      l = (*interactions)[2*k+1];
      int interacts = ((*interaction_map)[j] || (*interaction_map)[l] ) ? 1 : 0;

      for ( i=0; i<nprocs; ++i ) {
        dist_list[i] = 2*nprocs;
      }
      if ( interacts ) {
        if ( last_interaction )  {
          dist_list[myrank] = dist;
        }
        last_interaction = k;
        dist = 0;
      }
      else {
        ++dist;
      }

      MPI_Allreduce(MPI_IN_PLACE,dist_list,nprocs,MPI_INT,MPI_MIN,comm);

      if ( myrank==parmesh->info.root ) {
        int min = INT_MAX;
        for ( i=0; i<nprocs; ++i ) {
          min = MG_MIN ( dist_list[i],min );
        }
        fprintf ( stdout,"     interaction %8d, ranks %8d %8d, distance %8d.\n",
                  k,j,l,min );
      }
    }
    PMMG_DEL_MEM ( parmesh, dist_list, int,"dist_list" );
  }

  if ( !ier ) {
    return -1;
  }
  else {
    return ninteractions;
  }
}

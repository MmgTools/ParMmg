/**
 * \file loadbalancing_pmmg.c
 * \brief Load balancing after a remeshing step
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \author Nikos Pattakos (Inria)
 * \author Luca Cirrottola (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 */
#include "parmmg.h"
#include "metis_pmmg.h"

#warning Luca: to remove when nold_grp will be updated
/**
 * \param parmesh pointer toward the parmesh structure.
 * \return The number of groups on the current process once that elements have
 * been marked for splitting on a merged mesh on the current proc.
 *
 */
int PMMG_get_ngrp( PMMG_pParMesh parmesh ) {
  PMMG_pGrp   grp;
  MMG5_pMesh  mesh;
  MMG5_pTetra pt;
  int         ie,igrp,ngrp;

  /* It has to be called on a merged partition */
  assert( parmesh->ngrp == 1 );

  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;
 
  /* Retrieve the grp ID from the tetra mark field */
  ngrp = 0;
  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    if( !MG_EOK(pt) ) continue;
    igrp  = pt->mark / parmesh->nprocs;
    if( igrp > ngrp ) ngrp = igrp;
  }
  ngrp++;

  return ngrp;
}


/**
 * \param parmesh pointer toward the parmesh structure.
 * \param mesh pointer toward the mesh structure.
 * \param start index of the starting tetrahedra.
 * \param ip local index of the point in the tetrahedra \a start.
 * \param list pointer toward the list of the tetra in the volumic ball of
 * \a ip.
 * \return 0 if fail and the number of the tetra in the ball otherwise.
 *
 * Fill the volumic ball (i.e. filled with tetrahedra) of point \a ip in tetra
 * \a start. Results are stored under the form \f$4*kel + jel\f$, kel = number
 * of the tetra, jel = local index of p within kel.
 * Mark each tetrahedron in the ball with the maximum between its own value and
 * the value brought by the point.
 *
 */
int PMMG_mark_boulevolp( PMMG_pParMesh parmesh, MMG5_pMesh mesh,int ngrp,int base_front, int ip, int * list){
  MMG5_pTetra  pt,pt1;
  MMG5_pPoint  ppt1;
  int    *adja,nump,ilist,base,cur,k,k1;
  int     start,color,iloc;
  char    j,l,i;

  /* Get point color */
  start  = mesh->point[ip].s / 4;
  iloc   = mesh->point[ip].s % 4;
  color  = mesh->point[ip].tmp;

  base = ++mesh->base;
  pt   = &mesh->tetra[start];
  nump = pt->v[iloc];
  assert( nump == ip );

  /* Store initial tetrahedron */
  pt->flag = base;
  list[0] = 4*start + iloc;
  ilist=1;

  /* Explore list and travel by adjacency through elements sharing p */
  cur = 0;
  while ( cur < ilist ) {
    k = list[cur] / 4;
    i = list[cur] % 4; // index of point p in tetra k
    adja = &mesh->adja[4*(k-1)+1];

    for (l=0; l<3; l++) {
      i  = MMG5_inxt3[i];
      k1 = adja[i];
      if ( !k1 )  continue;
      k1 /= 4;
      pt1 = &mesh->tetra[k1];
      if ( pt1->flag == base )  continue;
      if ( pt1->mark < color ) {
        pt1->mark = color;
        for (j=0; j<4; j++) {
          ppt1 = &mesh->point[pt1->v[j]];
          /* Mark and flag new interface points */
          if ( ppt1->flag < base_front ) {
            ppt1->tmp  = color;
            ppt1->s    = 4*k1+j;
            ppt1->flag = base_front+1;
          }
        }
      }
      pt1->flag = base;
      for (j=0; j<4; j++)
        if ( pt1->v[j] == nump )  break;
      assert(j<4);
      /* overflow */
      if ( ilist > MMG3D_LMAX-3 )  return 0;
      list[ilist] = 4*k1+j;
      ilist++;
    }
    cur++;
  }
  return ilist;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param mesh pointer toward the mesh structure.
 *
 * \return The flag base value used to track the front.
 *
 * Mark interface points from the maximum tetra mark field in their ball, and
 * flag them as mesh->base.
 *
 */
int PMMG_mark_interfacePoints( PMMG_pParMesh parmesh, MMG5_pMesh mesh,int ngrp ) {
  MMG5_pTetra pt;
  MMG5_pPoint ppt;
  int         ip,ie,iloc;

  /* New base flag */
  mesh->base++;

  /* Reset point flag field */
  for( ip = 1; ip <= mesh->np; ip++ )
    mesh->point[ip].tmp = PMMG_UNSET;

  /* Mark interface points with the maximum color */
  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    if( !MG_EOK(pt) ) continue;

    for( iloc = 0; iloc < 4; iloc++ ) {
      ppt = &mesh->point[pt->v[iloc]];

      /* Mark new point */
      if( ppt->tmp == PMMG_UNSET ) {
        ppt->tmp = pt->mark;
        ppt->s   = 4*ie+iloc;
        continue;
      }

      /* Mark and flag interface point */
      if( ppt->tmp < pt->mark ) {
        ppt->tmp  = pt->mark;
        ppt->s    = 4*ie+iloc;
        ppt->flag = mesh->base;
      }
    }
  }

  return mesh->base;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param part groups partitions array
 *
 * \return 1 if success.
 *
 * Fill the groups partitions array by retrieving the proc ID from the mark
 * field of their first valid tetra.
 *
 */
int PMMG_part_getProcs( PMMG_pParMesh parmesh,int *part ) {
  PMMG_pGrp   grp;
  MMG5_pMesh  mesh;
  MMG5_pTetra pt;
  idx_t       *vtxdist;
  int         igrp,ie,ier;

  for( igrp = 0; igrp < parmesh->ngrp; igrp++ ) {
    grp  = &parmesh->listgrp[igrp];
    mesh = grp->mesh;
    for( ie = 1; ie <= mesh->ne; ie++ ) {
      pt = &mesh->tetra[ie];
      if( !MG_EOK(pt) ) continue;
      part[igrp] = pt->mark % parmesh->nprocs;
      break;
    }
  }
  ier = 1;

  PMMG_CALLOC(parmesh,vtxdist,parmesh->nprocs+1,idx_t,"parmetis vtxdist", return 0);

  MPI_CHECK( MPI_Allgather(&parmesh->ngrp,1,MPI_INT,&vtxdist[1],1,MPI_INT,parmesh->comm),
             PMMG_DEL_MEM(parmesh,vtxdist,idx_t,"parmetis vtxdist"); return 0 );
  for( int iproc = 0; iproc < parmesh->nprocs; iproc++ )
    vtxdist[iproc+1] += vtxdist[iproc];

  ier = PMMG_correct_parmeshGrps2parmetis( parmesh, vtxdist, part, parmesh->nprocs );

  PMMG_DEL_MEM(parmesh,vtxdist,idx_t,"parmetis vtxdist");

  return ier;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param part groups partitions array
 * \param ngrps array of nb of old grps on each proc
 *
 * \return the number of groups required on the current proc.
 *
 * Fill the groups partitions array by retrieving the grp ID from the mark field
 * of each tetra.
 *
 */
int PMMG_part_getInterfaces( PMMG_pParMesh parmesh,int *part,int *ngrps ) {
  PMMG_pGrp   grp;
  MMG5_pMesh  mesh;
  MMG5_pTetra pt;
  int *map_grps; 
  int         igrp,iproc,color;
  int         sumngrps[parmesh->nprocs+1];
  int         ie,i,count;

  /* It has to be called on a merged partition */
  assert( parmesh->ngrp == 1 );

  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;

  sumngrps[0] = 0;
  for( iproc = 0; iproc < parmesh->nprocs; iproc++ )
    sumngrps[iproc+1] = sumngrps[iproc]+ngrps[iproc];

  PMMG_CALLOC(parmesh,map_grps,sumngrps[parmesh->nprocs],int,"map_grps",
              return 0);
 
  /* Retrieve the grp ID from the tetra mark field */
  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    if( !MG_EOK(pt) ) continue;
    igrp  = pt->mark / parmesh->nprocs;
    iproc = pt->mark % parmesh->nprocs;
    color = igrp+sumngrps[iproc];
    map_grps[color] = 1;
    part[ie-1] = color;
  }

  /* Pack the map */
  count = 0;
  for( i = 0; i < parmesh->nprocs; i++ ) {
    iproc = ( i+parmesh->myrank ) % parmesh->nprocs;
    for( igrp = 0; igrp < ngrps[iproc]; igrp++ ) {
    color = igrp+sumngrps[iproc];
    if( map_grps[color] ) map_grps[color] = count++;
    }
  }

  /* Permute the partition array */
  for( ie = 1; ie <= mesh->ne; ie++ )
    part[ie-1] = map_grps[part[ie-1]];

  PMMG_DEL_MEM(parmesh,map_grps,int,"map_grps");

  /* Return the nb of groups on the current proc */
  return count;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param part groups partitions array
 *
 * Move old groups interfaces through an advancing-front method.
 *
 */
int PMMG_part_moveInterfaces( PMMG_pParMesh parmesh ) {
  PMMG_pGrp    grp;
  MMG5_pMesh   mesh;
  MMG5_pTetra  pt,pt1;
  MMG5_pxTetra pxt;
  MMG5_pPoint  ppt;
  PMMG_pInt_comm int_node_comm;
  PMMG_pExt_comm ext_node_comm;
  MPI_Comm       comm;
  MPI_Status     status;
  int          *node2int_node_comm_index1,*node2int_node_comm_index2;
  int          *intvalues,*itosend,*itorecv;
  int          nlayers;
  int          nprocs,ngrp,base_front;
  int          igrp,k,i,idx,ip,ie,ifac,je,ne,nitem,color,color_out;
  int          list[MMG3D_LMAX+2];
  int          ier=1;

  comm   = parmesh->comm;
  assert( parmesh->ngrp == 1 );
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;
  node2int_node_comm_index1 = grp->node2int_node_comm_index1;
  node2int_node_comm_index2 = grp->node2int_node_comm_index2;

  nprocs = parmesh->nprocs;
  ngrp   = PMMG_get_ngrp( parmesh );

  /* Mark interface points with the maximum color */
  base_front = PMMG_mark_interfacePoints( parmesh, mesh, ngrp );

  /* Reset internal communicator */
  int_node_comm = parmesh->int_node_comm;
  PMMG_CALLOC( parmesh,int_node_comm->intvalues,int_node_comm->nitem,int,"intvalues",return 0);
  intvalues = int_node_comm->intvalues;
  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    idx = node2int_node_comm_index2[i];
    intvalues[idx] = PMMG_UNSET;
  }

  /* Save grp index and proc in the internal communicator */
  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    idx = node2int_node_comm_index2[i];
    ip  = node2int_node_comm_index1[i];
    ppt = &mesh->point[ip];
    assert( MG_VOK(ppt) );
    intvalues[idx] = ppt->tmp;  // contains nprocs*igrp+iproc
  }

  /* Exchange values on the interfaces among procs */
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    nitem         = ext_node_comm->nitem;
    color         = ext_node_comm->color_out;

    PMMG_CALLOC(parmesh,ext_node_comm->itosend,nitem,int,"itosend array",
                return 0);
    itosend = ext_node_comm->itosend;

    PMMG_CALLOC(parmesh,ext_node_comm->itorecv,nitem,int,"itorecv array",
                return 0);
    itorecv       = ext_node_comm->itorecv;

    for ( i=0; i<nitem; ++i ) {
      idx            = ext_node_comm->int_comm_index[i];
      itosend[i]     = intvalues[idx] ;
    }

    MPI_CHECK(
      MPI_Sendrecv(itosend,nitem,MPI_INT,color,MPI_PARMESHGRPS2PARMETIS_TAG,
                   itorecv,nitem,MPI_INT,color,MPI_PARMESHGRPS2PARMETIS_TAG,
                   comm,&status),return 0 );

    for ( i=0; i<nitem; ++i ) {
      idx            = ext_node_comm->int_comm_index[i];
      intvalues[idx] = itorecv[i];
    }

  }

  /* Update grp index and proc after communication */
  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    idx = node2int_node_comm_index2[i];
    ip  = node2int_node_comm_index1[i];
    ppt = &mesh->point[ip];
    assert( MG_VOK(ppt) );
    if( intvalues[idx] > ppt->tmp ) {
      ppt->tmp = intvalues[idx];
      ppt->flag = mesh->base;
    }
  }

  /* Move interfaces */
  nlayers = 2;
  for( i = 0; i < nlayers; i++ ) {

    /* Mark tetra in the ball of interface points */
    for( ip = 1; ip <= mesh->np; ip++ ) {
      ppt = &mesh->point[ip];
      if( !MG_VOK(ppt) ) continue;

      /* Skip not-interface points */
      if( ppt->flag != base_front ) continue;

      /* Advance the front: New interface points will be flagged as
       * base_front+1 */
      ier = PMMG_mark_boulevolp( parmesh, mesh, ngrp, base_front, ip, list);
      if( !ier ) break;

    }
    if( !ier ) break;

    /* Update flag base for next wave */
    base_front++;
  }

  PMMG_DEL_MEM( parmesh,int_node_comm->intvalues,int,"intvalues" );
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    PMMG_DEL_MEM(parmesh,ext_node_comm->itosend,int,"itosend array");
    PMMG_DEL_MEM(parmesh,ext_node_comm->itorecv,int,"itorecv array");
  }
  
  return ier;
}

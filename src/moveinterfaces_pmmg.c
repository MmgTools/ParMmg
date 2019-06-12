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
int PMMG_mark_boulevolp( PMMG_pParMesh parmesh, MMG5_pMesh mesh, int base_front, int ip, int * list){
  MMG5_pTetra  pt,pt1;
  MMG5_pPoint  ppt1;
  int    *adja,nump,ilist,base,cur,k,k1;
  int     nprocs,ngrp,shift,start,color;
  char    j,l,i,j1;

  base = ++mesh->base;
  pt   = &mesh->tetra[start];
  nump = pt->v[ip];

  /* Get point color */
  nprocs = parmesh->nprocs;
  ngrp   = parmesh->ngrp;
  shift  = nprocs*ngrp;

  start  = mesh->point[ip].tmp / shift;
  color  = mesh->point[ip].tmp % shift;

  /* Store initial tetrahedron */
  pt->flag = base;
  list[0] = 4*start + ip;
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
      if ( pt1->mark <= color ) continue;
      pt1->flag = base;
      pt1->mark = color;
      for (j=0; j<4; j++) {
        if ( pt1->v[j] == nump )  j1 = j;
        else {
          ppt1 = &mesh->point[pt1->v[j]];
          /* Mark and flag new interface points */
          if ( ppt1->flag < base_front ) {
            ppt1->tmp = color + shift*k1;
            ppt1->flag = base_front+1;
          }
        }
      }
      assert(j1<4);
      /* overflow */
      if ( ilist > MMG3D_LMAX-3 )  return 0;
      list[ilist] = 4*k1+j1;
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
int PMMG_mark_interfacePoints( PMMG_pParMesh parmesh, MMG5_pMesh mesh ) {
  MMG5_pTetra pt;
  MMG5_pPoint ppt;
  int         nprocs,ngrp,shift;
  int         ip,ie,iloc;

  nprocs = parmesh->nprocs;
  ngrp   = parmesh->ngrp;
  shift  = nprocs*ngrp;

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
        ppt->tmp = pt->mark + shift*ie;
        continue;
      }

      /* Mark and flag interface point */
      if( (ppt->tmp % shift) < pt->mark ) {
        ppt->tmp = pt->mark+shift*ie;
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
 * Fill the groups partitions array by retrieving the grp ID from the mark field
 * of each tetra.
 *
 */
void PMMG_part_getInterfaces( PMMG_pParMesh parmesh,int *part ) {
  PMMG_pGrp   grp;
  MMG5_pMesh  mesh;
  MMG5_pTetra pt;
  int         ie;

  /* It has to be called on a merged partition */
  assert( parmesh->ngrp == 1 );

  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;

  /* Retrieve the grp ID from the tetra mark field */
  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    part[ie-1] = pt->mark/parmesh->nprocs;
  }
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param part groups partitions array
 *
 * Move old groups interfaces by retrieving the grp ID from the mark field
 * of each tetra.
 *
 */
int PMMG_part_moveInterfaces( PMMG_pParMesh parmesh ) {
  PMMG_pGrp    grp;
  MMG5_pMesh   mesh;
  MMG5_pTetra  pt,pt1;
  MMG5_pxTetra pxt;
  MMG5_pPoint  ppt;
  PMMG_pExt_comm ext_node_comm;
  MPI_Comm       comm;
  MPI_Status     status;
  int          *node2int_node_comm_index1,*node2int_node_comm_index2;
  int          *intvalues,*itosend,*itorecv;
  int          *adja;
  int          nlayers;
  int          nprocs,ngrp,shift,base_front;
  int          igrp,k,i,idx,ip,ie,ifac,je,ne,ne_min,nitem,color,color_out;
  int          list[MMG3D_LMAX+2];
  int          ier;

  ne_min = 6; //FIXME

  intvalues = parmesh->int_node_comm->intvalues;

  comm   = parmesh->comm;
  assert( parmesh->ngrp == 0 );
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;
  node2int_node_comm_index1 = grp->node2int_node_comm_index1;
  node2int_node_comm_index2 = grp->node2int_node_comm_index2;

  nprocs = parmesh->nprocs;
  ngrp   = parmesh->ngrp;
  shift  = nprocs*ngrp;

  /* Mark interface points with the maximum color */
  base_front = PMMG_mark_interfacePoints( parmesh, mesh );

  /* Reset internal communicator */
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
    intvalues[idx] = ppt->tmp % shift;  // contains nprocs*igrp+iproc
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
    ie = ppt->tmp / shift;
    if( intvalues[idx] > (ppt->tmp % shift) ) {
      ppt->tmp = intvalues[idx] + shift*ie;
      ppt->flag = mesh->base;
    }
  }

  /* Move interfaces */
  nlayers = 2;
  for( i = 0; i < nlayers; i++ ) {

    /* Mark interface points with the maximum color (only in the first wave */
    if( i ) base_front = PMMG_mark_interfacePoints( parmesh, mesh );

    /* Mark tetra in the ball of interface points */
    for( ip = 1; ip <= mesh->np; ip++ ) {
      ppt = &mesh->point[ip];
      if( !MG_VOK(ppt) ) continue;

      /* Skip not-interface points */
      if( ppt->flag != base_front ) continue;

      /* Advance the front: New interface points will be flagged as
       * base_front+1 */
      ier = PMMG_mark_boulevolp( parmesh, mesh, base_front, ip, list);
    }

    /* Update flag base for next wave */
    base_front++;
  }

  return 1;
}

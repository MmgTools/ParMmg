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



/**
 * \param parmesh pointer toward the parmesh structure.
 * \param color color of the group.
 * \return iproc the proc index.
 *
 * Get rank from a global group ID, according to the format
 * parmesh->nprocs*igrp+parmesh->myrank.
 *
 */
int PMMG_get_proc( PMMG_pParMesh parmesh,int color ) {
  return color % parmesh->nprocs;
}


/**
 * \param parmesh pointer toward the parmesh structure.
 * \param color color of the group.
 * \return igrp the group index.
 *
 * Get group from a global group ID, according to the format
 * parmesh->nprocs*igrp+parmesh->myrank.
 *
 */
int PMMG_get_grp( PMMG_pParMesh parmesh,int color ) {
  return color / parmesh->nprocs;
}

int PMMG_compare_grps( PMMG_pParMesh parmesh,int *vtxdist,int *map,
    int color0, int color1 ) {
  int igrp0,iproc0,idx0;
  int igrp1,iproc1,idx1;

  igrp0 = PMMG_get_grp( parmesh, color0 );
  igrp1 = PMMG_get_grp( parmesh, color1 );

  iproc0 = PMMG_get_proc( parmesh, color0 );
  iproc1 = PMMG_get_proc( parmesh, color1 );

  idx0 = igrp0 + vtxdist[iproc0];
  idx1 = igrp1 + vtxdist[iproc1];

  assert( map[idx0] && map[idx1] );

  if( map[idx0] > map[idx1] )
    return 1;
  else
    return 0;

}
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param igrp index of the group.
 *
 * Set a global group ID, according to the format
 * parmesh->nprocs*igrp+parmesh->myrank.
 *
 */
int PMMG_set_color( PMMG_pParMesh parmesh,int igrp ) {
  return parmesh->nprocs*igrp+parmesh->myrank;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param igrp index of the group.
 *
 * Store a global group ID into the tetra mark field.
 *
 */
void PMMG_set_color_tetra( PMMG_pParMesh parmesh,int igrp ) {
  PMMG_pGrp   grp;
  MMG5_pMesh  mesh;
  MMG5_pTetra pt;
  int         color,ie;

  grp   = &parmesh->listgrp[igrp];
  mesh  = grp->mesh;
  color = PMMG_set_color(parmesh,igrp);

  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    if( !MG_EOK(pt) ) continue;
    pt->mark = color;
  }
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param mesh pointer toward the mesh structure.
 * \param color color of the group to scan.
 * \param len length of the list to be merged.
 * \param otetra tetra of different color to be used for merging.
 * \return 0 if fail, 1 if success.
 *
 * Merge the subgroup into a neighbour subgroup with different color.
 *
 */
int PMMG_merge_subgroup( PMMG_pParMesh parmesh,MMG5_pMesh mesh,int color,
                         int *list,int len,int otetra ) {
  MMG5_pTetra      pt,pto;
  int              cur,k;

  pto = &mesh->tetra[otetra];

  for( cur = 0; cur < len; cur++ ) {
    k = list[cur];
    pt = &mesh->tetra[k];
    /** Merge it */
    pt->mark = pto->mark;
    pt->flag = pto->flag;
  }

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param mesh pointer toward the mesh structure.
 * \param color color of the group to scan.
 * \return 0 if fail, 1 if success.
 *
 */
int PMMG_list_contiguous( PMMG_pParMesh parmesh,MMG5_pMesh mesh,
                           int start,int *list,int *list_head,int *list_len,
                           int *list_base,int *list_otetra ) {
  MMG5_pTetra      pt,pt1;
  int              *adja,ilist,cur,k,k1,l;
  int              base,color;

  /* New flavour */
  base = ++mesh->base;

  /* Get color */
  color = mesh->tetra[start].mark;

  /* Store initial tetrahedron */
  list[0] = start;
  ilist = 1;

  /* Flag initial tetra as treated */
  mesh->tetra[start].flag = base;

  /** Explore list and fill it by adjacency */
  cur = 0;
  while ( cur < ilist ) {
    k = list[cur];
    pt = &mesh->tetra[k];
    adja = &mesh->adja[4*(k-1)+1];

    /* Add neighbours to the list */
    for (l=0; l<4; l++) {
      k1 = adja[l];
      if ( !k1 )  continue;
      k1 /= 4;
      pt1 = &mesh->tetra[k1];
      /* Skip already visited tetra (by this or another list */
      if ( pt1->flag )  continue;
      /* Skip tetra with different color */
      if ( pt1->mark != color ) continue;
      /* Flag tetra as treated */
      pt1->flag = base;
      /* Add tetra to the list */
      assert( ilist <= mesh->ne );
      list[ilist] = k1;
      ilist++;
    }
    cur++;
  }

  /** Return list head, length, flavour, and a neighbour color */
  *list_head   = start;
  *list_len    = ilist;
  *list_base   = base;
  *list_otetra = PMMG_UNSET;
  for( cur = 0; cur < ilist; cur++ ) {
    k = list[cur];
    pt = &mesh->tetra[k];
    adja = &mesh->adja[4*(k-1)+1];
    /* Add neighbours to the list */
    for (l=0; l<4; l++) {
      k1 = adja[l];
      if ( !k1 )  continue;
      k1 /= 4;
      pt1 = &mesh->tetra[k1];
      if ( pt1->flag != base ) {
        *list_otetra = k1;
        break;
      }
    }
    if( *list_otetra != PMMG_UNSET ) break;
  }

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param igrp index of the group to check.
 * \return 0 if fail, 1 if success.
 *
 * Check contiguity of the group mesh.
 *
 */
int PMMG_check_contiguity( PMMG_pParMesh parmesh,int igrp ) {
  MMG5_pMesh const mesh = parmesh->listgrp[igrp].mesh;
  MMG5_pTetra      pt;
  int              *list;
  int              next_head,next_len,next_base,next_otetra;
  int              start,k,counter;

  PMMG_MALLOC(parmesh,list,mesh->ne,int,"tetra list",return 0);

  /* Initialize counter */
  counter = 0;

  /* Reset tetra flag */
  for( k = 1; k <= mesh->ne; k++ )
    mesh->tetra[k].flag = 0;

  /** 1) Find the first subgroup */
  start = 1;
  if( !PMMG_list_contiguous( parmesh, mesh, start, list, &next_head,
        &next_len, &next_base, &next_otetra ) ) return 0;
  counter += next_len;

  /** Find the next list head */
  start++;
  while( start <= mesh->ne ) {
    pt = &mesh->tetra[start];
    if( !pt->flag ) break;
    start++;
  }

  /** 2) Look for new subgroups until all the mesh is scanned */
  while( start <= mesh->ne ) {

    if( !PMMG_list_contiguous( parmesh, mesh, start, list, &next_head,
          &next_len, &next_base, &next_otetra ) ) return 0;
    counter += next_len;

    /* Find the next list head */
    start++;
    while( start <= mesh->ne ) {
      pt = &mesh->tetra[start];
      if( !pt->flag ) break;
      start++;
    }
  }

  /* Check that all the elements have been listed */
  for( start = 1; start <= mesh->ne; start++ )
    assert( mesh->tetra[start].flag );
  assert( counter == mesh->ne );

  PMMG_DEL_MEM(parmesh,list,int,"tetra list");

  return 1;
}


/**
 * \param parmesh pointer toward the parmesh structure.
 * \param igrp index of the group to check.
 * \param color color of the group to make contiguous.
 * \param main_list pointer to the first tetra list.
 * \param next_list pointer to the second tetra list.
 * \param counter pointer to the remaining number of tetra of the given color
 * \return 0 if fail, 1 if success.
 *
 */
int PMMG_fix_subgrp_contiguity( PMMG_pParMesh parmesh,int color,int *main_list,
                                int *next_list,int *counter ) {
  MMG5_pMesh const mesh = parmesh->listgrp[0].mesh;
  MMG5_pTetra      pt;
  int              main_head,main_len,main_base,main_otetra;
  int              next_head,next_len,next_base,next_otetra;
  int              start,k;

  /* Only works on a merged group */
  assert( parmesh->ngrp == 1 );

  /** 1) Find the first subgroup with the given color */
  start = 1;
  while( start <= mesh->ne ) {
    pt = &mesh->tetra[start];
    if( pt->mark == color ) break; /* break when a new subgroup is found */
    start++;
  }

  if( start <= mesh->ne ) {
    if( !PMMG_list_contiguous( parmesh, mesh, start, main_list, &main_head,
          &main_len, &main_base, &main_otetra ) ) return 0;
    *counter += main_len;
  }

  /** Find the next list head */
  start++;
  while( start <= mesh->ne ) {
    pt = &mesh->tetra[start];
    if( !pt->flag && (pt->mark == color) ) break;
    start++;
  }

  /** 2) Look for new subgroups until all the mesh is scanned */
  while( start <= mesh->ne ) {

    if( !PMMG_list_contiguous( parmesh, mesh, start, next_list, &next_head,
          &next_len, &next_base, &next_otetra ) ) return 0;
    *counter += next_len;


    /* Compare the next subgroup with the main one */
    if( next_len > main_len ) {
      /* Merge main */
      if( main_otetra == PMMG_UNSET ) {
        fprintf(stderr,"\n### Error: Cannot merge main subgroup on proc %d\n",parmesh->myrank);
        return 0;
      } else {
        if( !PMMG_merge_subgroup( parmesh, mesh, color, main_list, main_len, main_otetra ) )
          return 0;
        /* Decrease counter if the merged list will be scanned again */
        if( !mesh->tetra[main_otetra].flag ) *counter -= main_len;
        /* Swap */
        main_head   = next_head;
        main_len    = next_len;
        main_base   = next_base;
        main_otetra = next_otetra;
      }
    } else {
      /* Merge next */
      if( next_otetra == PMMG_UNSET ) {
        fprintf(stderr,"\n### Error: Cannot merge next subgroup on proc %d\n",parmesh->myrank);
        return 0;
      } else {
        if( !PMMG_merge_subgroup( parmesh, mesh, color, next_list, next_len, next_otetra ) )
          return 0;
        /* Decrease counter if the merged list will be scanned again */
        if( !mesh->tetra[next_otetra].flag ) *counter -= next_len;
      }
    }

    /* Find the next list head */
    start++;
    while( start <= mesh->ne ) {
      pt = &mesh->tetra[start];
      if( !pt->flag && (pt->mark == color) ) break;
      start++;
    }
  }

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param counter pointer to the number of tetra of counted/OK tetra.
 * \return 0 if fail, 1 if success.
 *
 */
int PMMG_fix_contiguity( PMMG_pParMesh parmesh,int *counter ) {
  MMG5_pMesh const mesh = parmesh->listgrp[0].mesh;
  int              *main_list,*next_list;
  int              k,color,igrp;

  /* Only works on a merged group */
  assert( parmesh->ngrp == 1 );

  /* Allocate tetra lists */
  PMMG_MALLOC(parmesh,main_list,mesh->ne,int,"tetra main list",return 0);
  PMMG_MALLOC(parmesh,next_list,mesh->ne,int,"tetra next list",return 0);

  /* Initialize list counter */
  *counter = 0;

  /* Reset the tetra flag */
  for( k = 1; k <= mesh->ne; k++ )
    mesh->tetra[k].flag = 0;

  /* Loop on the old groups */
  for( igrp = 0; igrp < parmesh->nold_grp; igrp++ ) {

    /* Get group color */
    color = PMMG_set_color( parmesh, igrp );

    /* Check and fix contiguity of the old group */
    if( !PMMG_fix_subgrp_contiguity( parmesh, color, main_list, next_list,
          counter ) ) {
      PMMG_DEL_MEM(parmesh,main_list,int,"tetra main list");
      PMMG_DEL_MEM(parmesh,next_list,int,"tetra next list");
      return 0;
    };

  }

  /* Check that all the contiguous tetra have been visited */
#ifndef NDEBUG
  int count = 0;
  for( k = 1; k <= mesh->ne; k++ )
    if( mesh->tetra[k].flag ) count++;
  assert( count == *counter );
#endif

  /* Deallocate lists and return */
  PMMG_DEL_MEM(parmesh,main_list,int,"tetra main list");
  PMMG_DEL_MEM(parmesh,next_list,int,"tetra next list");

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param counter pointer to the number of tetra of the grp color
 * \return 0 if fail, 1 if success.
 *
 * Check that subgroups created by interface migration are reachable from
 * the target groups.
 *
 */
int PMMG_check_reachability( PMMG_pParMesh parmesh,int *counter ) {
  PMMG_pGrp    grp;
  MMG5_pMesh   mesh;
  MMG5_pTetra  pt;
  PMMG_pInt_comm int_face_comm;
  PMMG_pExt_comm ext_face_comm;
  MPI_Comm       comm;
  MPI_Status     status;
  int          *face2int_face_comm_index1,*face2int_face_comm_index2;
  int          *intvalues,*itosend,*itorecv,rank_out;
  int          *list;
  int          next_head,next_len,next_base,next_otetra;
  int          nitem,color;
  int          ie,i,idx,k;

  comm   = parmesh->comm;
  assert( parmesh->ngrp == 1 );
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;
  face2int_face_comm_index1 = grp->face2int_face_comm_index1;
  face2int_face_comm_index2 = grp->face2int_face_comm_index2;


  /* Reset internal communicator */
  int_face_comm = parmesh->int_face_comm;
  PMMG_CALLOC( parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,"intvalues",return 0);
  intvalues = int_face_comm->intvalues;
  for( i = 0; i < grp->nitem_int_face_comm; i++ ) {
    idx = face2int_face_comm_index2[i];
    intvalues[idx] = PMMG_UNSET;
  }

  /* Save grp index and proc in the internal communicator */
  for( i = 0; i < grp->nitem_int_face_comm; i++ ) {
    idx = face2int_face_comm_index2[i];
    ie  = face2int_face_comm_index1[i]/12;
    pt = &mesh->tetra[ie];
    assert( MG_EOK(pt) );
    intvalues[idx] = pt->mark;
  }

  /* Exchange values on the interfaces among procs */
  for ( k = 0; k < parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    nitem         = ext_face_comm->nitem;
    color         = ext_face_comm->color_out;

    PMMG_CALLOC(parmesh,ext_face_comm->itosend,nitem,int,"itosend array",
                return 0);
    itosend = ext_face_comm->itosend;

    PMMG_CALLOC(parmesh,ext_face_comm->itorecv,nitem,int,"itorecv array",
                return 0);
    itorecv       = ext_face_comm->itorecv;

    for ( i=0; i<nitem; ++i ) {
      idx            = ext_face_comm->int_comm_index[i];
      itosend[i]     = intvalues[idx] ;
    }

#warning Luca: change this tag
    MPI_CHECK(
      MPI_Sendrecv(itosend,nitem,MPI_INT,color,MPI_PARMESHGRPS2PARMETIS_TAG,
                   itorecv,nitem,MPI_INT,color,MPI_PARMESHGRPS2PARMETIS_TAG,
                   comm,&status),return 0 );

    for ( i=0; i<nitem; ++i ) {
      idx            = ext_face_comm->int_comm_index[i];
      intvalues[idx] = itorecv[i];
    }

  }

  PMMG_MALLOC(parmesh,list,mesh->ne,int,"tetra list",return 0);

  /* Ignore values if coming from a proc different than color_out */
  for( k = 0; k < parmesh->next_face_comm; k++ ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    rank_out = ext_face_comm->color_out;
    for( i = 0; i < ext_face_comm->nitem; i++ ) {
      idx = ext_face_comm->int_comm_index[i];
      if( PMMG_get_proc( parmesh,intvalues[idx] ) != rank_out )
        intvalues[idx] = PMMG_UNSET;
    }
  }

  /* Reach as many tetra as possible from the interface through walk search */
  for( i = 0; i < grp->nitem_int_face_comm; i++ ) {
    idx = face2int_face_comm_index2[i];
    ie  = face2int_face_comm_index1[i]/12;
    pt  = &mesh->tetra[ie];
    color = intvalues[idx];
    /* Skip faces not on the external comm */
    if( color == PMMG_UNSET ) continue;
    /* Skip tetra with color different from the interface */
    if( pt->mark != color ) continue;
    /* Skip already seen tetra */
    if( pt->flag ) continue;
    /* Flag the reachable adjacents */
    if( !PMMG_list_contiguous( parmesh, mesh, ie, list, &next_head,
          &next_len, &next_base, &next_otetra ) ) return 0;
    *counter += next_len;
    assert( *counter <= mesh->ne );
  }

#ifndef NDEBUG
  int count = 0;
  for( k = 1; k <= mesh->ne; k++ )
    if( mesh->tetra[k].flag ) count++;
  assert( count == *counter );
#endif

  /* Merge the unreachable subgroups */
  ie = 1;
  while( ie <= mesh->ne ) {
    pt = &mesh->tetra[ie];
    if( !pt->flag ) break; /* break when an unseen element is found */
    ie++;
  }

  while( ie <= mesh->ne ) {
    color = pt->mark;

    if( !PMMG_list_contiguous( parmesh, mesh, ie, list, &next_head,
          &next_len, &next_base, &next_otetra ) ) return 0;
    *counter += next_len;

    if( next_otetra == PMMG_UNSET ) {
      fprintf(stderr,"\n### Error: Cannot merge unreachable subgroup on proc %d\n",parmesh->myrank);
      return 0;
    } else {
      printf("Merging unseen %d into %d \n",color,mesh->tetra[next_otetra].mark);
      if( !PMMG_merge_subgroup( parmesh, mesh, color, list, next_len, next_otetra ) )
        return 0;
      /* Decrease counter if the merged list will be scanned again */
      if( !mesh->tetra[next_otetra].flag ) *counter -= next_len;
    }
    assert( *counter <= mesh->ne );

    /* Find the next list head */
    ie++;
    while( ie <= mesh->ne ) {
      pt = &mesh->tetra[ie];
      if( !pt->flag ) break;
      ie++;
    }
  }


  /* Check that all tetra have been visited, and taken into account in a list */
#ifndef NDEBUG
  for( ie = 1; ie <= mesh->ne; ie++ )
    assert( mesh->tetra[ie].flag );
#endif
  assert( *counter == mesh->ne );


  PMMG_DEL_MEM( parmesh,int_face_comm->intvalues,int,"intvalues" );
  for ( k = 0; k < parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    PMMG_DEL_MEM(parmesh,ext_face_comm->itosend,int,"itosend array");
    PMMG_DEL_MEM(parmesh,ext_face_comm->itorecv,int,"itorecv array");
  }

  PMMG_DEL_MEM(parmesh,list,int,"tetra list");

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param ngrps pointer to the number of groups on each proc.
 * \return The number of groups on the current process after front advancement..
 *
 */
int PMMG_count_grpsPerProc( PMMG_pParMesh parmesh,int *ngrps ) {
  MMG5_pMesh const mesh = parmesh->listgrp[0].mesh;
  MMG5_pTetra pt;
  int sumngrps[parmesh->nprocs+1];
  int *map_grps;
  int ie,igrp,iproc,color;
  int ngrp;

  sumngrps[0] = 0;
  for( iproc = 0; iproc < parmesh->nprocs; iproc++ )
    sumngrps[iproc+1] = sumngrps[iproc]+ngrps[iproc];

  PMMG_CALLOC(parmesh,map_grps,sumngrps[parmesh->nprocs],int,"map_grps",
              return 0);
#warning Luca: Largely inefficient, loop on communicators instead
  /* Retrieve the grp ID from the tetra mark field */
  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    if( !MG_EOK(pt) ) continue;
    igrp = PMMG_get_grp( parmesh, pt->mark );
    iproc = PMMG_get_proc( parmesh, pt->mark );
    color = igrp+sumngrps[iproc];
    map_grps[color] = 1;
  }

  /* Count grps on local proc */
  ngrp = 0;
  for( igrp = 0; igrp < sumngrps[parmesh->nprocs]; igrp++ ) {
    if( map_grps[igrp] ) ngrp++;
  }

  PMMG_DEL_MEM(parmesh,map_grps,int,"map_grps");
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
int PMMG_mark_boulevolp( PMMG_pParMesh parmesh,MMG5_pMesh mesh,int *vtxdist,
    int *map,int *nelem,int ngrp,int base_front, int ip, int * list){
  MMG5_pTetra  pt,pt1;
  MMG5_pPoint  ppt1;
  int    *adja,nump,ilist,base,cur,k,k1,j1;
  int    igrp;
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
  /** Analyse initial tetra: It could be on the other side of the front */
  if ( PMMG_compare_grps( parmesh, vtxdist, map, pt->mark, color ) ) {
    if( PMMG_get_proc( parmesh, pt->mark ) == parmesh->myrank ) {
      igrp = PMMG_get_grp( parmesh, pt->mark );
      if( nelem[igrp] == PMMG_METIS_NELEM_MIN ) return 1;
      nelem[igrp]--;
    }
    pt->mark = color;
    for (j=0; j<4; j++) {
      ppt1 = &mesh->point[pt->v[j]];
      /* Mark and flag points not on the current front */
      if ( ppt1->flag != base_front ) {
        ppt1->tmp  = color;
        ppt1->s    = 4*start+j;
        ppt1->flag = base_front+1;
      }
    }
  }

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
      /* Skip already visited tetra */
      if ( pt1->flag == base )  continue;
      /* Get local node index */
      for (j=0; j<4; j++)
        if ( pt1->v[j] == nump ) { j1 = j; break; }
      assert(j1<4);
      /** Flag not-owned tetra but don't put it in the list */
      if ( PMMG_compare_grps( parmesh, vtxdist, map, color, pt1->mark ) ) {
        pt1->flag = base;
        continue;
      }
      /** Mark owned tetra and its vertices */
      if ( PMMG_compare_grps( parmesh, vtxdist, map, pt1->mark, color ) ) {
        if( PMMG_get_proc( parmesh, pt1->mark ) == parmesh->myrank ) {
          igrp = PMMG_get_grp( parmesh, pt1->mark );
          if( nelem[igrp] == PMMG_METIS_NELEM_MIN ) return 1;
          nelem[igrp]--;
        }
        pt1->mark = color;
        for (j=0; j<4; j++) {
          ppt1 = &mesh->point[pt1->v[j]];
          /* Mark and flag points not on the current front */
          if ( ppt1->flag != base_front ) {
            ppt1->tmp  = color;
            ppt1->s    = 4*k1+j;
            ppt1->flag = base_front+1;
          }
        }
      }
      /** Flag tetra and put it in the list */
      pt1->flag = base;
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
int PMMG_mark_interfacePoints( PMMG_pParMesh parmesh,MMG5_pMesh mesh,int ngrp,
                               int *vtxdist,int *map ) {
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
      if( PMMG_compare_grps( parmesh,vtxdist,map,ppt->tmp,pt->mark ) ) {
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
      part[igrp] = PMMG_get_proc( parmesh, pt->mark );
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
    igrp  = PMMG_get_grp( parmesh, pt->mark );
    iproc = PMMG_get_proc( parmesh, pt->mark );
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
int PMMG_init_map( PMMG_pParMesh parmesh,int **vtxdist,int **map ) {
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_pTetra    pt;
  MPI_Comm       comm;
  int            igrp;
  int            ngrp,nproc,myrank,k,ie;

  assert( parmesh->ngrp == 1 );

  comm   = parmesh->comm;
  myrank = parmesh->myrank;
  ngrp   = parmesh->nold_grp;
  nproc  = parmesh->nprocs;

  grp                       = &parmesh->listgrp[0];
  mesh                      = grp->mesh;

  /** Step 1: Fill vtxdist array with the range of groups local to each
   * processor */
  PMMG_CALLOC(parmesh,*vtxdist,nproc+1,int,"vtxdist", return 0);

  MPI_CHECK( MPI_Allgather(&ngrp,1,MPI_INT,&(*vtxdist)[1],1,MPI_INT,comm),
             return 0 );

  for ( k=1; k<=nproc; ++k )
    (*vtxdist)[k] += (*vtxdist)[k-1];

  PMMG_CALLOC(parmesh,*map,(*vtxdist)[nproc],int,"map", return 0);

  /* Count the nb of tetra for each old group */
  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    if( !MG_EOK(pt) ) continue;
    igrp = PMMG_get_grp( parmesh, pt->mark );
    ++(*map)[ igrp + (*vtxdist)[parmesh->myrank] ];
  }

  return 1;
}

int PMMG_sort_procs( PMMG_pParMesh parmesh,int **vtxdist,int **map ) {
  PMMG_pGrp      grp;
  PMMG_pExt_comm ext_node_comm;
  PMMG_pInt_comm int_node_comm;
  MMG5_pMesh     mesh;
  MMG5_pPoint    ppt;
  MPI_Comm       comm;
  MPI_Status     status;
  int            *node2int_node_comm_index1,*node2int_node_comm_index2;
  int            *intvalues,*itosend,*itorecv;
  int            color,igrp,iproc;
  int            ngrp,nproc,myrank,nitem,k,i,idx,ip;

  assert( parmesh->ngrp == 1 );

  comm   = parmesh->comm;
  myrank = parmesh->myrank;
  ngrp   = parmesh->nold_grp;
  nproc  = parmesh->nprocs;

  grp                       = &parmesh->listgrp[0];
  mesh                      = grp->mesh;
  node2int_node_comm_index1 = grp->node2int_node_comm_index1;
  node2int_node_comm_index2 = grp->node2int_node_comm_index2;

  /** Step 2: Fill the internal communicator.*/
  int_node_comm = parmesh->int_node_comm;
  PMMG_MALLOC(parmesh,int_node_comm->intvalues,int_node_comm->nitem,int,
              "node communicator",return 0);

  /* Node communicator initialization */
  intvalues = parmesh->int_node_comm->intvalues;

  /*Fill the internal communicator */
  for ( k=0; k<grp->nitem_int_node_comm; ++k ) {
    ip   =  node2int_node_comm_index1[k];
    ppt = &mesh->point[ip];
    assert( MG_VOK(ppt) );
    intvalues[node2int_node_comm_index2[k]] = ppt->tmp;
  }


  /** Step 3: Send and receive external communicators filled by the (group id +
   * ishift) of the neighbours (through the nodes) */
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    nitem         = ext_node_comm->nitem;
    color         = ext_node_comm->color_out;

    PMMG_CALLOC(parmesh,ext_node_comm->itosend,nitem,int,"itosend array",
                goto fail_3);
    itosend = ext_node_comm->itosend;

    PMMG_CALLOC(parmesh,ext_node_comm->itorecv,nitem,int,"itorecv array",
                goto fail_4);
    itorecv       = ext_node_comm->itorecv;

    for ( i=0; i<nitem; ++i ) {
      idx            = ext_node_comm->int_comm_index[i];
      itosend[i]     = intvalues[idx] ;
    }

    MPI_CHECK(
      MPI_Sendrecv(itosend,nitem,MPI_INT,color,MPI_PARMESHGRPS2PARMETIS_TAG,
                   itorecv,nitem,MPI_INT,color,MPI_PARMESHGRPS2PARMETIS_TAG,
                   comm,&status),goto fail_5 );
  }

  /** Step 4: Process the external communicators to count for each group the
   * adjacent groups located on another processor and fill the sorted linked
   * list of adjacency */

  for (  k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    itosend       = ext_node_comm->itosend;
    itorecv       = ext_node_comm->itorecv;
    nitem         = ext_node_comm->nitem;

    /* i2send array contains the group id of the boundary nodes and i2recv the
     * group id of the same node in the other proc */
    for ( i=0; i<nitem; ++i ) {
      /* Get the group id (+ishift) of the node in our proc and the group id
       * (+ishift) of the node in the adjacent proc */
      color = itorecv[i];
      igrp  = PMMG_get_grp( parmesh, color );
      iproc = PMMG_get_proc( parmesh, color );
      (*map)[ igrp + (*vtxdist)[iproc] ] = 1;


    }
  }

  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    if ( ext_node_comm->itorecv )
      PMMG_DEL_MEM(parmesh,ext_node_comm->itorecv,int,"itorecv array");
    if ( ext_node_comm->itosend )
      PMMG_DEL_MEM(parmesh,ext_node_comm->itosend,int,"itosend array");
  }
 
  /** Step 5: Send and receive external communicators filled by the (group id +
   * ishift) of the neighbours (through the nodes) */
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    color         = ext_node_comm->color_out;
    nitem         = (*vtxdist)[myrank+1]-(*vtxdist)[myrank];

    PMMG_CALLOC(parmesh,ext_node_comm->itosend,nitem,int,"itosend array",
                goto fail_3);
    itosend = ext_node_comm->itosend;

    PMMG_CALLOC(parmesh,ext_node_comm->itorecv,nitem,int,"itorecv array",
                goto fail_4);
    itorecv       = ext_node_comm->itorecv;

    for ( i=0; i<nitem; ++i ) {
      idx = i + (*vtxdist)[myrank];
      itosend[i] = (*map)[idx];
    }

    MPI_CHECK(
      MPI_Sendrecv(itosend,nitem,MPI_INT,color,MPI_PARMESHGRPS2PARMETIS_TAG,
                   itorecv,nitem,MPI_INT,color,MPI_PARMESHGRPS2PARMETIS_TAG,
                   comm,&status),goto fail_5 );
  }

  /** Step 6: Process the external communicators to count for each group the
   * adjacent groups located on another processor and fill the sorted linked
   * list of adjacency */

  for (  k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    itosend       = ext_node_comm->itosend;
    itorecv       = ext_node_comm->itorecv;
    color         = ext_node_comm->color_out;
    nitem         = (*vtxdist)[color+1]-(*vtxdist)[color];

    /* i2send array contains the group id of the boundary nodes and i2recv the
     * group id of the same node in the other proc */
    for ( i=0; i<nitem; ++i ) {
      idx = i + (*vtxdist)[color];
      if( (*map)[idx] ) (*map)[idx] = itorecv[i];
    }
  }

  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    if ( ext_node_comm->itorecv )
      PMMG_DEL_MEM(parmesh,ext_node_comm->itorecv,int,"itorecv array");
    if ( ext_node_comm->itosend )
      PMMG_DEL_MEM(parmesh,ext_node_comm->itosend,int,"itosend array");
  }
  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,int,"node communicator");
  return 1;

fail_5:
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    if ( ext_node_comm->itorecv )
      PMMG_DEL_MEM(parmesh,ext_node_comm->itorecv,int,"itorecv array");
  }
fail_4:
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    if ( ext_node_comm->itosend )
      PMMG_DEL_MEM(parmesh,ext_node_comm->itosend,int,"itosend array");
  }
fail_3:
  if ( int_node_comm->intvalues )
    PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,int,"node communicator");
  return 0;
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
  int          *vtxdist,*map,*nelem;
  int          nlayers;
  int          nprocs,ngrp,base_front;
  int          igrp,k,i,idx,ip,ie,ifac,je,ne,nitem,color,color_out;
  int          list[MMG3D_LMAX+2];
  int          ier=1,ier_glob;

  comm   = parmesh->comm;
  assert( parmesh->ngrp == 1 );
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;
  node2int_node_comm_index1 = grp->node2int_node_comm_index1;
  node2int_node_comm_index2 = grp->node2int_node_comm_index2;

  nprocs = parmesh->nprocs;
  ngrp   = parmesh->nold_grp;

  if( !PMMG_init_map( parmesh, &vtxdist, &map ) ) return 0;

  /* Mark each point with the maximum color among the tetras in the ball,
   * flag interface points */
  base_front = PMMG_mark_interfacePoints( parmesh, mesh, ngrp, vtxdist, map );

  if( !PMMG_sort_procs( parmesh, &vtxdist, &map ) ) return 0;

  PMMG_CALLOC( parmesh,nelem,parmesh->nold_grp,int,"nelem",return 0);
  for( igrp = 0; igrp < parmesh->nold_grp; igrp++ )
    nelem[igrp] = map[igrp+vtxdist[parmesh->myrank]];

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
    intvalues[idx] = ppt->tmp;  // contains the point color
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

#warning Luca: change this tag
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
    if( PMMG_compare_grps( parmesh, vtxdist, map, ppt->tmp, intvalues[idx] ) ) {
      ppt->tmp = intvalues[idx];
      ppt->flag = base_front;
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
      ier = PMMG_mark_boulevolp( parmesh, mesh, vtxdist, map, nelem, ngrp,
                                 base_front, ip, list);
      if( !ier ) break;

    }
    if( !ier ) break;

    /* Update flag base for next wave */
    base_front++;
  }

#ifndef NDEBUG
  PMMG_check_contiguity( parmesh,0 );
#endif
  int counter;

  ier = PMMG_fix_contiguity( parmesh, &counter );
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);
  if( !ier_glob ) return 0;

  ier = PMMG_check_reachability( parmesh, &counter );
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);
  if( !ier_glob ) return 0;

  PMMG_DEL_MEM( parmesh,nelem,int,"nelem" );
  PMMG_DEL_MEM( parmesh,int_node_comm->intvalues,int,"intvalues" );
  for ( k = 0; k < parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    PMMG_DEL_MEM(parmesh,ext_node_comm->itosend,int,"itosend array");
    PMMG_DEL_MEM(parmesh,ext_node_comm->itorecv,int,"itorecv array");
  }

  PMMG_DEL_MEM(parmesh,map,int,"map");
  PMMG_DEL_MEM(parmesh,vtxdist,int,"vtxdist");

  return ier;
}

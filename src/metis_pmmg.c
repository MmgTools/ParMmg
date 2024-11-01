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
 * \file metis_pmmg.c
 * \brief Partition mesh using metis
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "metis_pmmg.h"
#include "linkedlist_pmmg.h"
#include "mmg3dexterns_private.h"

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param xadj array of shifts (CSR).
 * \param adjncy array of adjacents (CSR).
 * \param adjwgt array of weights (CSR).
 * \param filename filename prefix.
 *
 * \return 1 if success, 0 if fail.
 *
 * Save mesh graph and weight on file in Medit format.
 *
 */
int PMMG_saveGraph( PMMG_pParMesh parmesh,idx_t *xadj,idx_t *adjncy,
    idx_t *adjwgt,const char *filename ) {
  PMMG_pGrp  grp  = &parmesh->listgrp[0];
  MMG5_pMesh mesh = grp->mesh;
  MMG5_pTetra pt;
  MMG5_pPoint ppt;
  char *sname,*smesh,*ssol;
  FILE *fmesh,*fsol;
  double wgt = 0.0;
  idx_t istart,istop,iadj,nadj;
  int ip,k,j,jel;

  size_t snamelen = strlen(filename)+1+4*sizeof(int);
  PMMG_CALLOC(parmesh,sname,snamelen,char,"file name prefix",return 0);
  PMMG_CALLOC(parmesh,smesh,snamelen+5,char,"mesh file name",return 0);
  PMMG_CALLOC(parmesh,ssol,snamelen+4,char,"sol file name",return 0);

  snprintf(sname,snamelen,"%s-P%02d-I%02d",filename,parmesh->myrank,parmesh->iter);
  strcpy(smesh,sname);
  strcat(smesh,".mesh");
  strcpy(ssol,sname);
  strcat(ssol,".sol");

  /* Open files and write headers */
  fmesh = fopen(smesh,"w");
  fprintf(fmesh,"MeshVersionFormatted 2\n");
  fprintf(fmesh,"\nDimension 3\n");

  fsol  = fopen(ssol,"w");
  fprintf(fsol,"MeshVersionFormatted 2\n");
  fprintf(fsol,"\nDimension 3\n");

  /* Write vertices */
  fprintf(fmesh,"\nVertices\n%d\n",mesh->np);
  for( ip = 1; ip <= mesh->np; ip++ ) {
    ppt = &mesh->point[ip];
    fprintf(fmesh,"%f %f %f %d\n",
        ppt->c[0],
        ppt->c[1],
        ppt->c[2],
        ppt->ref);
  }

  /* Write triangles and solution on triangles */
  fprintf(fmesh,"\nTriangles\n%d\n",xadj[mesh->ne]/2);
  fprintf(fsol,"\nSolAtTriangles\n%d\n %d %d\n",xadj[mesh->ne]/2,1,MMG5_Scalar);
  for( k = 1; k <= mesh->ne; k++ ) {
    pt   = &mesh->tetra[k];
    istart = xadj[k-1];
    istop  = xadj[k];
    nadj   = istop-istart;
    for ( j = 0; j < nadj; j++ ) {
      iadj = istart+j;
      jel = adjncy[iadj]+1;
      if ( jel < k ) continue;

      /* Save triangle */
      fprintf(fmesh,"%d %d %d %d\n",
          pt->v[MMG5_idir[j][0]],
          pt->v[MMG5_idir[j][1]],
          pt->v[MMG5_idir[j][2]],
          0);

      /* Save weight */
      if( adjwgt ) wgt = (double)adjwgt[iadj];
      fprintf(fsol,"%f\n",wgt);
    }
  }

  /* Close files */
  fprintf(fmesh,"\n\nEnd");
  fprintf(fsol,"\n\nEnd");
  fclose(fmesh);
  fclose(fsol);

  /* Free memory */
  PMMG_DEL_MEM(parmesh,sname,char,"file name prefix");
  PMMG_DEL_MEM(parmesh,smesh,char,"mesh file name");
  PMMG_DEL_MEM(parmesh,ssol,char,"sol file name");
  return 1;
}

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

  for (k=0; k<hash->siz; ++k ) {
    hash->item[k].adj = PMMG_UNSET;
    hash->item[k].wgt = PMMG_NUL;
  }

  for (k=hash->siz; k<hash->max; ++k) {
    hash->item[k].adj = PMMG_UNSET;
    hash->item[k].wgt = PMMG_NUL;
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
int PMMG_hashGrp( PMMG_pParMesh parmesh,PMMG_HGrp *hash, int k, idx_t adj,
                  idx_t wgt ) {
  PMMG_hgrp  *ph;
  int        tmp_nxt,j,newsize;

  ph  = &hash->item[k+1];

  if ( ph->adj == adj ) {
    ph->wgt += wgt;
    return 2;
  }

  tmp_nxt = 0;
  if ( PMMG_UNSET != ph->adj ) {
    while ( ph->nxt && (ph->nxt<hash->max) && (hash->item[ph->nxt].adj<adj) )
      ph = &hash->item[ph->nxt];

    if ( hash->item[ph->nxt].adj == adj ) {
      hash->item[ph->nxt].wgt += wgt;
      return 2;
    }

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
        hash->item[j].wgt = PMMG_NUL;
        hash->item[j].nxt = j+1;
      }
    }
    hash->nxt = ph->nxt;
  }

  /* insert new group */
  ph->adj  = adj;
  ph->wgt += wgt;
  ph->nxt  = tmp_nxt;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met  pointer toward the met structure.
 * \param tag  face tag to be checked
 *
 * Compute and store metis weight in the tetra qual field.
 *
 */
void PMMG_computeWgt_mesh( MMG5_pMesh mesh,MMG5_pSol met,int tag ) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  int          ie,ifac;

  /* Reset quality field */
  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    if( !MG_EOK(pt) ) continue;
    if( !pt->xt ) continue;
    pt->qual = 0.0;
  }

  /* INcrement weight for a given face tag */
  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    if( !MG_EOK(pt) ) continue;
    if( !pt->xt ) continue;
    pxt = &mesh->xtetra[pt->xt];
    for( ifac = 0; ifac < 4; ifac++ )
      if( pxt->ftag[ifac] & tag )
        pt->qual += PMMG_computeWgt( mesh, met, pt, ifac );
  }

}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met  pointer toward the met structure.
 * \param pt   pointer toward the tetrahedron structure.
 * \param ifac face index of the tetrahedron.
 *
 * \return The weight value
 *
 * Compute an element weight to be used for the metis weight on parallel
 * interfaces.
 *
 */
double PMMG_computeWgt( MMG5_pMesh mesh,MMG5_pSol met,MMG5_pTetra pt,int ifac ) {
  double       len,res,alpha=28.0;
  int          i,ia;

  if ( met && met->m ) {
    res = 0.0;
    for( i=0; i<3; i++ ) {
      ia = MMG5_iarf[ifac][i];
      len = MMG5_lenedg(mesh,met,ia,pt);
      if( len <= 1.0 )
        res += len-1.0;
      else
        res += 1.0/len-1.0;
    }
    res = MG_MIN(1.0/exp(alpha*res/3.0),PMMG_WGTVAL_HUGEINT);
  }
  else {
    res = PMMG_WGTVAL_HUGEINT;
  }

  return res;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return The number of contiguous subgroups.
 *
 * Check group mesh contiguity by counting the number of adjacent element
 * subgroups.
 *
 */
int PMMG_check_part_contiguity( PMMG_pParMesh parmesh,idx_t *xadj,idx_t *adjncy,
    idx_t *part, idx_t nnodes, idx_t nparts ) {
  idx_t       *flag,ipart,inode,jnode,iadj,istart;
  int         nnodes_part,ncolors,maxcolors,nb_seen,done;

  /** Flags */
  PMMG_CALLOC( parmesh, flag, nnodes, idx_t, "graph node flags", return 0 );

  /** Count the nb. of graph subgroups for each part */
  maxcolors = 0;
  for( ipart = 0; ipart < nparts; ipart++ ) {

    /** Reset node flag of the current partition and count nodes on partition */
    nnodes_part = 0;
    for( inode = 0; inode < nnodes; inode++ )
      if( part[inode] == ipart ) {
        flag[inode] = 0;
        ++nnodes_part;
      }

    /** Loop on all nodes to count the subgroups */
    ncolors = 0;
    nb_seen = 0;
    while( nb_seen < nnodes_part ) {

      /** Find first node not seen */
      for( inode = 0; inode < nnodes; inode++ ) {
        /* Remain on current partition */
        if( part[inode] != ipart ) continue;

        if( flag[inode] == 0 ) {
          /* New color */
          ++ncolors;
          istart = inode;
          /* Mark element */
          flag[inode] = ncolors;
          ++nb_seen;
          done = 0;
          /* Exit */
          break;
        }
      }

      /** Loop on colored elements until the subgroup is full */
      while( !done ) {
        done = 1;
        for( inode = istart; inode < nnodes; inode++ ) {
          /* Remain on current partition */
          if( part[inode] != ipart ) continue;

          /* Skip unseen nodes or different colors */
          if( flag[inode] != ncolors ) continue;

          /** Loop on adjacents */
          for( iadj = xadj[inode]; iadj < xadj[inode+1]; iadj++ ) {
            jnode = adjncy[iadj];
            if( part[jnode] != ipart ) continue;
            /** Mark with current color (if not already seen) */
            if( flag[jnode] == 0 ) {
              flag[jnode] = ncolors;
              ++nb_seen;
              done = 0;
            }
          }
        }
      }
    }

    if( ncolors > 1 && parmesh->ddebug ) {
      fprintf(stderr,"\n  ## Warning: %d contiguous subgroups found on part %d, proc %d.\n",
              ncolors,ipart,parmesh->myrank);
    }

    /** Update the max nb of subgroups found */
    if( ncolors > maxcolors ) maxcolors = ncolors;
  }

  PMMG_DEL_MEM( parmesh,flag,idx_t,"graph node flags");

  return ncolors;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 if success
 *
 * Check group mesh contiguity and reset global option for forcing contiguity
 * in Metis if discontiguous groups are found.
 *
 */
int PMMG_checkAndReset_grps_contiguity( PMMG_pParMesh parmesh ) {
  int contigresult, ier;

  /** Check grps contiguity */
  if( (parmesh->info.loadbalancing_mode == PMMG_LOADBALANCING_metis) &&
      (parmesh->info.contiguous_mode) ) {

    ier = PMMG_check_grps_contiguity( parmesh );
    if( !ier ) {
      if ( parmesh->ddebug ) {
        fprintf(stderr,"\n  ## Error %s: Unable to count mesh contiguous subgroups.\n",
                __func__);
      }
    } else if( ier>1 ) {
      if ( parmesh->ddebug ) {
        fprintf(stderr,"\n  ## Warning %s: Group meshes are not contiguous. Reverting to discontiguous mode.\n",
                __func__);
      }
      parmesh->info.contiguous_mode = PMMG_NUL;
      ier = 1;
    }

    /* Check that the same option is applied on all procs */
    MPI_Allreduce( &parmesh->info.contiguous_mode, &contigresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
    parmesh->info.contiguous_mode = contigresult;

  } else {
    ier = 1;
  }

  return ier;
}


/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return The number of contiguous subgroups.
 *
 * Check group mesh contiguity by counting the number of adjacent element
 * subgroups.
 *
 */
int PMMG_check_grps_contiguity( PMMG_pParMesh parmesh ) {
  PMMG_pGrp   grp;
  MMG5_pMesh  mesh;
  MMG5_pTetra pt,pt1;
  int         *adja,igrp,ie,je,ifac;
  int         ncolors,maxcolors,nb_seen,mark_notseen,done,istart;

  /** Labels */
  mark_notseen = 0;

  /** Count the nb. of mesh subgroups for each group */
  maxcolors = 0;
  for( igrp = 0; igrp < parmesh->ngrp; igrp++ ) {
    grp  = &parmesh->listgrp[igrp];
    mesh = grp->mesh;

    if ( !mesh->adja ) {
      if ( !MMG3D_hashTetra(mesh,0) ) {
        fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
        return 0;
      }
    }

    /** Reset tetra flag */
    for( ie = 1; ie < mesh->ne+1; ie++ )
      mesh->tetra[ie].flag = mark_notseen;

    /** Loop on all elements to count the subgroups */
    ncolors = 0;
    nb_seen = 0;
    while( nb_seen < mesh->ne ) {

      /** Find first element not seen */
      for( ie = 1; ie < mesh->ne+1; ie++ ) {
        pt = &mesh->tetra[ie];
        if( pt->flag == mark_notseen ) {
          /* New color */
          ++ncolors;
          istart = ie;
          /* Mark element */
          pt->flag = ncolors;
          ++nb_seen;
          done = 0;
          /* Exit */
          break;
        }
      }

      /** Loop on colored elements until the subgroup is full */
      while( !done ) {
        done = 1;
        for( ie = istart; ie < mesh->ne+1; ie++ ) {
          pt   = &mesh->tetra[ie];
          adja = &mesh->adja[4*(ie-1)+1];

          /* Skip unseen elts or different colors */
          if( pt->flag != ncolors ) continue;

          /** Loop on adjacents */
          for( ifac = 0; ifac < 4; ifac++ ) {
            je = adja[ifac]/4;
            pt1 = &mesh->tetra[je];
            /** Mark with current color (if not already seen) */
            if( je && pt1->flag == mark_notseen ) {
              pt1->flag = ncolors;
              ++nb_seen;
              done = 0;
            }
          }
        }
      }
    }

    if ( ncolors > 1 && parmesh->ddebug ) {
      fprintf(stderr,"\n  ## Warning: %d contiguous subgroups found on grp %d, proc %d.\n",
              ncolors,igrp,parmesh->myrank);
    }

    /** Update the max nb of subgroups found */
    if( ncolors > maxcolors ) maxcolors = ncolors;
  }

  return ncolors;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param part    elements partition array
 * \param ne      nb of elements
 * \param nproc   nb of groups for partitioning
 *
 * \return 1 if no empty partitions or successfully corrected, 0 if fail
 *
 * Check if metis has returned empty partitions, correct partitioning if so.
 *
 */
int PMMG_correct_meshElts2metis( PMMG_pParMesh parmesh,idx_t* part,idx_t ne,idx_t nproc ) {
  PMMG_valLnkdList **partlist;
  idx_t            iproc,ie;
  int              nempt,iempt;


  /* Initialize lists */
  PMMG_CALLOC(parmesh,partlist,nproc,PMMG_valLnkdList*,"array of list pointers",return 0);
  for( iproc=0; iproc<nproc; iproc++ ) {
    PMMG_CALLOC(parmesh,partlist[iproc],1,PMMG_valLnkdList,"linked list pointer",return 0);
    if( !PMMG_valLnkdListNew(parmesh,partlist[iproc],iproc,PMMG_LISTSIZE) ) return 0;
  }

  /* Fill the lists */
  assert ( ne >= nproc &&  "not enough elements for the number of partitions" );

  for( ie=0; ie<ne; ie++ ) {
    iproc = part[ie];
    if( !PMMG_add_val2lnkdList(parmesh,partlist[iproc],ie) ) return 0;
  }

  /* Sort lists based on nb. of entities, in ascending order */
  qsort(partlist,nproc,sizeof(PMMG_valLnkdList*),PMMG_compare_valLnkdListLen);

  /* Count empty partitions */
  nempt = 0;
  for( iproc=0; iproc<nproc; iproc++ ) {
    if( partlist[iproc]->nitem ) break;
    nempt++;
  }

  assert( nempt < nproc );
  if( !nempt ) {
    /* Deallocate lists and return */
    for( iproc=0; iproc<nproc; iproc++ ) {
      PMMG_DEL_MEM(parmesh,partlist[iproc]->item,PMMG_lnkdVal,"linked list array");
      PMMG_DEL_MEM(parmesh,partlist[iproc],PMMG_valLnkdList,"linked list pointer");
    }
    PMMG_DEL_MEM(parmesh,partlist,PMMG_valLnkdList*,"array of linked lists");

    return 1;
  }

  fprintf(stdout,"   ### Warning: Empty partitions on proc %d, nelts %d\n",parmesh->myrank,ne);
  /** Correct partitioning */
  assert ( nproc > 1 );
  iproc = nproc-1;
  iempt = 0;
  while( nempt ) {
    /* Get next "reservoir" proc */
    if ( iproc == nproc-1 ) {
      while( partlist[iproc]->nitem <= partlist[iproc-1]->nitem ) {

        /* list are sorted depending to their number of items so iproc has more
         * items than iproc-1 */
        assert ( partlist[iproc]->nitem == partlist[iproc-1]->nitem );
      iproc--;
      }
      iproc--;
    }
    ++iproc;

    /* if ne > nproc, normally, we can fill the empty procs without emptying a
     * proc with only 1 item */
    assert ( partlist[iproc]->nitem > 1 && "not enough elements for the"
             " number of partitions");

    /* Pop entity ie from iproc, add to iempt */
    if( !PMMG_pop_val_lnkdList(parmesh,partlist[iproc],&ie) ) return 0;
    if( !PMMG_add_val2lnkdList(parmesh,partlist[iempt],ie) ) return 0;
    /* Update partition table and go on to next empty proc */
    part[ie] = partlist[iempt]->id;
    iempt++;
    nempt--;
  }

  /* Deallocate lists */
  for( iproc=0; iproc<nproc; iproc++ ) {
    PMMG_DEL_MEM(parmesh,partlist[iproc]->item,PMMG_lnkdVal,"linked list array");
    PMMG_DEL_MEM(parmesh,partlist[iproc],PMMG_valLnkdList,"linked list pointer");
  }
  PMMG_DEL_MEM(parmesh,partlist,PMMG_valLnkdList*,"array of linked lists");

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param vtxdist parmetis structure for nb of groups on each proc
 * \param mypart  local groups partition array
 * \param nproc   nb of procss for partitioning
 *
 * \return 1 if no empty partitions or successfully corrected, 0 if fail
 *
 * Check if parmetis has returned empty partitions, correct partitioning if so.
 *
 */
int PMMG_correct_parmeshGrps2parmetis( PMMG_pParMesh parmesh,idx_t *vtxdist,
                                       idx_t* mypart,idx_t nproc ) {
  PMMG_valLnkdList **partlist;
  idx_t             *part;
  idx_t              iproc,ie,ne,*recvcounts;
  int                myrank,ngrp,nempt,iempt;
  MPI_Comm           comm;

  myrank   = parmesh->myrank;
  ngrp     = parmesh->ngrp;
  ne       = vtxdist[nproc];
  comm     = parmesh->comm;

  /** Step 1: Fill part array with the partitions local to each processor */
  PMMG_CALLOC(parmesh,part,ne,idx_t,"parmetis part", return 0);
  PMMG_CALLOC(parmesh,recvcounts,nproc,idx_t,"recvcounts", return 0);

  for( iproc = 0; iproc<nproc; iproc++ )
    recvcounts[iproc] = vtxdist[iproc+1]-vtxdist[iproc];

  MPI_CHECK( MPI_Allgatherv(mypart,ngrp,MPI_INT,
                            part,recvcounts,vtxdist,MPI_INT,comm), return 0);

  PMMG_DEL_MEM(parmesh,recvcounts,idx_t,"recvcounts");

  /* Initialize lists */
  PMMG_CALLOC(parmesh,partlist,nproc,PMMG_valLnkdList*,"array of list pointers",return 0);
  for( iproc=0; iproc<nproc; iproc++ ) {
    PMMG_CALLOC(parmesh,partlist[iproc],1,PMMG_valLnkdList,"linked list pointer",return 0);
    if( !PMMG_valLnkdListNew(parmesh,partlist[iproc],iproc,PMMG_LISTSIZE) ) return 0;
  }

  /* Fill the lists */
  for( ie=0; ie<ne; ie++ ) {
    iproc = part[ie];
    if( !PMMG_add_val2lnkdList(parmesh,partlist[iproc],ie) ) return 0;
  }

  /* Sort lists based on nb. of entities, in ascending order */
  qsort(partlist,nproc,sizeof(PMMG_valLnkdList*),PMMG_compare_valLnkdListLen);

  /* Count empty partitions */
  nempt = 0;
  for( iproc=0; iproc<nproc; iproc++ )
    if( !partlist[iproc]->nitem ) nempt++;
  assert( nempt < nproc );
  if( !nempt ) {
    /* Deallocate and return */
    PMMG_DEL_MEM(parmesh,part,idx_t,"parmetis part");

    for( iproc=0; iproc<nproc; iproc++ ) {
      PMMG_DEL_MEM(parmesh,partlist[iproc]->item,PMMG_lnkdVal,"linked list array");
      PMMG_DEL_MEM(parmesh,partlist[iproc],PMMG_valLnkdList,"linked list pointer");
    }
    PMMG_DEL_MEM(parmesh,partlist,PMMG_valLnkdList*,"array of linked lists");

   return 1;
  }


  /** Correct partitioning */
  iempt = 0;
  while( nempt ) {
    /* Get next "reservoir" proc */
    iproc = nproc-1;
    while( partlist[iproc]->nitem <= partlist[iproc-1]->nitem )
      iproc--;
    /* Pop entity ie from iproc, add to iempt */
    if( !PMMG_pop_val_lnkdList(parmesh,partlist[iproc],&ie) ) return 0;
    if( !PMMG_add_val2lnkdList(parmesh,partlist[iempt],ie) ) return 0;
    /* Update partition table and go on to next empty proc */
    part[ie] = partlist[iempt]->id;
    iempt++;
    nempt--;
  }

  /** Update the local part */
  for( ie=0; ie<ngrp; ie++ )
    mypart[ie] = part[vtxdist[myrank]+ie];

  /* Deallocations */
  PMMG_DEL_MEM(parmesh,part,idx_t,"parmetis part");

  for( iproc=0; iproc<nproc; iproc++ ) {
    PMMG_DEL_MEM(parmesh,partlist[iproc]->item,PMMG_lnkdVal,"linked list array");
    PMMG_DEL_MEM(parmesh,partlist[iproc],PMMG_valLnkdList,"linked list pointer");
  }
  PMMG_DEL_MEM(parmesh,partlist,PMMG_valLnkdList*,"array of linked lists");

  return 1;
}


/**
 * \param parmesh pointer toward the PMMG parmesh structure
 * \param mesh pointer toward a MMG5 mesh structure
 * \param xadj pointer toward the position of the elt adjacents in adjncy
 * \param adjncy pointer toward the list of the adjacent of each elt
 * \param nadjncy number of data in adjncy array
 *
 * \return  1 if success, 0 if fail
 *
 * Build the metis graph with the mesh elements as metis nodes.
 *
 * \warning the mesh must be packed
 *
 */
int PMMG_graph_meshElts2metis( PMMG_pParMesh parmesh,MMG5_pMesh mesh,MMG5_pSol met,
                               idx_t **xadj,idx_t **adjncy,idx_t **adjwgt,
                               idx_t *nadjncy ) {
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  int          *adja;
  int          j,k,iadr,jel,count,nbAdj,wgt,ier;

  /** Step 1: mesh adjacency creation */

  if ( (!mesh->adja) && (1 != MMG3D_hashTetra( mesh, 1 )) ) {
    fprintf( stderr,"  ## PMMG Hashing problem (1).\n" );
    return 0;
  }

  /** Step 2: build the metis graph */

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
  PMMG_CALLOC(parmesh, (*adjncy), (*nadjncy), idx_t, "allocate adjncy", ier=0;);
  if( !ier ) {
    PMMG_DEL_MEM(parmesh, (*xadj), idx_t, "deallocate xadj" );
    return ier;
  }
  /* Don't compute weights at mesh distribution, or if output load balancing is required at last iter */
  if( (parmesh->iter != PMMG_UNSET) &&
      ((parmesh->iter < parmesh->niter-1) || parmesh->info.nobalancing) ) {
    PMMG_CALLOC(parmesh, (*adjwgt), (*nadjncy), idx_t, "allocate adjwgt", ier=0;);
    if( !ier ) {
      PMMG_DEL_MEM(parmesh, (*xadj), idx_t, "deallocate xadj" );
      PMMG_DEL_MEM(parmesh, (*adjncy), idx_t, "deallocate adjncy" );
      return ier;
    }
  }

  count = 0;
  for( k = 1; k <= mesh->ne; k++ ) {
    iadr = 4*(k-1) + 1;
    adja = &mesh->adja[iadr];
    pt   = &mesh->tetra[k];
    for ( j = 0; j < 4; j++ ) {
      jel = adja[j] / 4;
      if ( !jel ) continue;

      /* Assign graph edge weights */
      if( *adjwgt ) {
        /* Compute weight using face edge size */
        wgt = (int)PMMG_computeWgt(mesh,met,pt,j);

        (*adjwgt)[count] = MG_MAX(wgt,1);
      }

      (*adjncy)[count++]   = jel-1;
    }
    assert( count == ( (*xadj)[k] ) );
  }

  return ier;
}

/**
 * \param parmesh pointer toward the PMMG parmesh structure
 * \param vtxdist pointer toward the description of the node distribution
 * \param xadj pointer toward the position of the elt adjacents in adjncy
 * \param adjncy pointer toward the list of the adjacent of each elt
 * \param nadjncy number of data in adjncy array
 * \param vwgt pointer toward the metis node weights
 * \param adjvwgt pointer toward the metis edge weights
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
                                     idx_t **vwgt,idx_t **adjwgt,idx_t *wgtflag,
                                     idx_t *numflag,idx_t *ncon,idx_t nproc,
                                     real_t **tpwgts,real_t **ubvec) {
  PMMG_pGrp      grp;
  PMMG_pExt_comm ext_face_comm;
  PMMG_pInt_comm int_face_comm;
  PMMG_HGrp      hash;
  PMMG_hgrp      *ph;
  MMG5_pMesh     mesh;
  MMG5_pSol      met;
  MMG5_pTetra    pt;
  MMG5_pxTetra   pxt;
  MPI_Comm       comm;
  MPI_Status     status;
  int            *face2int_face_comm_index1,*face2int_face_comm_index2;
  int            *intvalues,*itosend,*itorecv;
  double         *doublevalues,*rtosend,*rtorecv;
  int            found,color;
  int            ngrp,myrank,nitem,k,igrp,igrp_adj,i,idx,ie,ifac,ishift,wgt;

  if( (parmesh->iter == parmesh->niter-1) && !parmesh->info.nobalancing ) {
    /* Switch off weights for output load balancing */
    *wgtflag = PMMG_WGTFLAG_NONE;
  } else {
    /* Default weight choice for parmetis */
    *wgtflag = PMMG_WGTFLAG_DEF;
  }
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
   * groups to which the face belong. Use a minus sign to mark old parallel
   * faces.*/
  PMMG_CALLOC(parmesh,*xadj,ngrp+1,idx_t,"parmetis xadj", goto fail_4);

  int_face_comm = parmesh->int_face_comm;

  PMMG_MALLOC(parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,
              "face communicator",goto fail_5);
  PMMG_CALLOC(parmesh,int_face_comm->doublevalues,int_face_comm->nitem,double,
              "face communicator",goto fail_5);

  /* Face communicator initialization */
  intvalues = parmesh->int_face_comm->intvalues;
  doublevalues = parmesh->int_face_comm->doublevalues;
  for ( k=0; k < int_face_comm->nitem; ++k )
    intvalues[k] = PMMG_UNSET;

  /*Fill the internal communicator with the greater index of the 2 groups to
   * which the face belong (igrp+ishift, to avoid ambiguity on grp 0 and with
   * PMMG_UNSET) */
  ishift = abs(PMMG_UNSET)+1;
  for ( igrp=ngrp-1; igrp>=0; --igrp ) {
    grp                       = &parmesh->listgrp[igrp];
    mesh                      = grp->mesh;
    met                       = grp->met;
    face2int_face_comm_index1 = grp->face2int_face_comm_index1;
    face2int_face_comm_index2 = grp->face2int_face_comm_index2;

    for ( k=0; k<grp->nitem_int_face_comm; ++k ) {
      ie   =  face2int_face_comm_index1[k]/12;
      ifac = (face2int_face_comm_index1[k]%12)/3;
      pt = &mesh->tetra[ie];
      assert( MG_EOK(pt) && pt->xt );
      pxt = &mesh->xtetra[pt->xt];

      /* Increase grp weight by the inverse of the interface element quality */
      if( pxt->ftag[ifac] & MG_OLDPARBDY )
        doublevalues[face2int_face_comm_index2[k]] += PMMG_computeWgt(mesh,met,pt,ifac);

      if ( PMMG_UNSET == intvalues[face2int_face_comm_index2[k] ] ) {

        /* Save group ID with a minus sign if the face was parallel in the
         * previous adaptation iteration */
        if( pxt->ftag[ifac] & MG_OLDPARBDY )
          intvalues[face2int_face_comm_index2[k]]= -(igrp+ishift);
        else
          intvalues[face2int_face_comm_index2[k]]= igrp+ishift;

      }
    }
  }

  /** Step 4: Send and receive external communicators filled by the (group id +
   * ishift) of the neighbours (through the faces) */
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    nitem         = ext_face_comm->nitem;
    color         = ext_face_comm->color_out;

    /* If the communicator is empty, dont try to communicate. This case happens when
       the mesh was loaded from an HDF5 file with more procs than partitions. */
    if ( !nitem ) continue;

    PMMG_CALLOC(parmesh,ext_face_comm->itosend,nitem,int,"itosend array",
                goto fail_6);
    itosend = ext_face_comm->itosend;

    PMMG_CALLOC(parmesh,ext_face_comm->itorecv,nitem,int,"itorecv array",
                goto fail_6);
    itorecv       = ext_face_comm->itorecv;

    PMMG_CALLOC(parmesh,ext_face_comm->rtosend,nitem,double,"rtosend array",
                goto fail_6);
    rtosend = ext_face_comm->rtosend;

    PMMG_CALLOC(parmesh,ext_face_comm->rtorecv,nitem,double,"rtorecv array",
                goto fail_6);
    rtorecv       = ext_face_comm->rtorecv;

    for ( i=0; i<nitem; ++i ) {
      idx            = ext_face_comm->int_comm_index[i];
      rtosend[i]     = doublevalues[idx] ;
      itosend[i]     = intvalues[idx] ;
      /* Mark the face as boundary in the intvalues array */
      intvalues[idx] = PMMG_UNSET;
    }

    MPI_CHECK(
      MPI_Sendrecv(itosend,nitem,MPI_INT,color,MPI_PARMESHGRPS2PARMETIS_TAG,
                   itorecv,nitem,MPI_INT,color,MPI_PARMESHGRPS2PARMETIS_TAG,
                   comm,&status),goto fail_6 );
    MPI_CHECK(
      MPI_Sendrecv(rtosend,nitem,MPI_DOUBLE,color,MPI_PARMESHGRPS2PARMETIS_TAG,
                   rtorecv,nitem,MPI_DOUBLE,color,MPI_PARMESHGRPS2PARMETIS_TAG,
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
    rtosend       = ext_face_comm->rtosend;
    rtorecv       = ext_face_comm->rtorecv;
    nitem         = ext_face_comm->nitem;

    /* i2send array contains the group id of the boundary faces and i2recv the
     * group id of the same face in the other proc */
    for ( i=0; i<nitem; ++i ) {
      /* Get the group id (+ishift) of the face in our proc and the group id
       * (+ishift) of the face in the adjacent proc */
      igrp     = itosend[i];
      igrp_adj = itorecv[i];

      /* Put high weight on old parallel faces (marked by a minus sign)
       * and reset signs; then, get grp indices starting from 0 */
      if( igrp_adj <= -ishift ) {
        assert( igrp <= -ishift );
        wgt = (int)( rtosend[i]+rtorecv[i] );
        igrp     *= -1;
        igrp_adj *= -1;
      } else {
        assert( igrp >=  ishift );
        wgt = 0;
      }
      igrp     -= ishift;
      igrp_adj -= ishift;

      assert ( igrp != PMMG_UNSET );
      assert ( igrp_adj != PMMG_UNSET );

      /* Search the neighbour in the hash table and insert it if not found */
      found = PMMG_hashGrp(parmesh,&hash,igrp,igrp_adj
                           +(*vtxdist)[ext_face_comm->color_out],wgt);

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

      /* Put high weight on old parallel faces (marked by a minus sign) and
       * reset signs; then, get grp index starting from 0 */
      if ( igrp_adj <= -ishift ) {
        wgt = (int)doublevalues[face2int_face_comm_index2[i]];
        igrp_adj *= -1;
      } else {
        wgt = 0;
      }
      igrp_adj -= ishift;

      if ( igrp_adj==PMMG_UNSET || igrp_adj<=igrp ) continue;

      /* Insert igrp_adj in the sorted list of igrp if not already found,
       * increment weight if found */
      found = PMMG_hashGrp(parmesh,&hash,igrp,igrp_adj+(*vtxdist)[myrank],wgt);
      if ( !found ) {
        fprintf(stderr,"  ## Error: %s: unable to add a new group in adjacency"
                " hash table.\n",__func__);
        goto fail_7;
      }
      /* Count group only if not already in the adja list */
      if ( found!=2 ) ++(*xadj)[ igrp+1 ];

      /* Insert igrp in the sorted list of igrp_adj if not already found,
       * increment weight if found */
      found = PMMG_hashGrp(parmesh,&hash,igrp_adj,igrp+(*vtxdist)[myrank],wgt);
      if ( !found ) {
        fprintf(stderr,"  ## Error: %s: unable to add a new group in adjacency"
                " hash table.\n",__func__);
        goto fail_7;
      }
      /* Count group only if not already in the adja list */
      if ( found !=2 ) ++(*xadj)[ igrp_adj+1 ];
    }
  }
#ifndef NDEBUG
  /* Check that the sum of shared faces for each grp is equal to the nb. of its
   * faces in the communicator */
  for ( igrp=0; igrp<ngrp; ++igrp ) {
    grp  = &parmesh->listgrp[igrp];

    found = 0;
    ph = &hash.item[igrp+1];
    if( PMMG_UNSET==ph->adj ) continue;


    found += ph->wgt;
    while ( ph->nxt ) {
      ph = &hash.item[ph->nxt];
      found += ph->wgt;
    }
  }
#endif

  /** Step 7: xadj array contains the number of adja per group, fill it for
   * Metis (it must contains the index of the first adja of the group in the
   * array adjcncy that list all the group adja in a 1D array) */
  for ( igrp=1; igrp <= parmesh->ngrp; ++igrp )
    (*xadj)[igrp] += (*xadj)[igrp-1];

  /** Step 8: Fill adjncy array at metis format */
  PMMG_CALLOC(parmesh,*adjncy,(*xadj)[ngrp],idx_t,"adjcncy parmetis array",
              goto fail_7);
  PMMG_CALLOC(parmesh,*adjwgt,(*xadj)[ngrp],idx_t,"parmetis adjwgt",
              goto fail_7);

  (*nadjncy) = 0;
  for ( igrp=0; igrp<=ngrp; ++igrp ) {

    ph = &hash.item[igrp+1];
    if ( PMMG_UNSET==ph->adj ) continue;

    (*adjncy)[(*nadjncy)]   = ph->adj;
    (*adjwgt)[(*nadjncy)++] = MG_MAX(ph->wgt,1);

    while ( ph->nxt ) {
      ph                = &hash.item[ph->nxt];
      (*adjncy)[(*nadjncy)]   = ph->adj;
      (*adjwgt)[(*nadjncy)++] = MG_MAX(ph->wgt,1);
    }
  }
  assert ( (*nadjncy)==(*xadj)[ngrp] );

#ifndef NDEBUG
  /* Print graph to file */
/*  FILE* fid;
  char filename[48];
  int len = PMMG_count_digits(parmesh->myrank);
  snprintf(filename,11+len*sizeof(int),"graph_proc%d",parmesh->myrank);
  fid = fopen(filename,"w");
  for( k = 0; k < parmesh->nprocs+1; k++ )
    fprintf(fid,"%d\n",(*vtxdist)[k]);
  fprintf(fid,"\n");
  for( k = 0; k < parmesh->ngrp+1; k++ )
    fprintf(fid,"%d\n",(*xadj)[k]);
  fprintf(fid,"\n");
  for( k = 0; k < (*xadj)[parmesh->ngrp]; k++ )
    fprintf(fid,"%d %d\n",(*adjncy)[k],(*adjwgt)[k]);  
  fclose(fid);*/
#endif

  /* Nullify unnecessary weights */
  switch (*wgtflag) {
    case PMMG_WGTFLAG_NONE:
      PMMG_DEL_MEM(parmesh,*vwgt,idx_t,"parmetis vwgt");
      PMMG_DEL_MEM(parmesh,*adjwgt,idx_t,"parmetis adjwgt");
      *vwgt = *adjwgt = NULL;
      break;
    case PMMG_WGTFLAG_ADJ:
      PMMG_DEL_MEM(parmesh,*vwgt,idx_t,"parmetis vwgt");
      *vwgt = NULL;
      break;
    case PMMG_WGTFLAG_VTX:
      PMMG_DEL_MEM(parmesh,*adjwgt,idx_t,"parmetis adjwgt");
      *adjwgt = NULL;
      break;
    default:
      break;
  }

  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    nitem = ext_face_comm->nitem;

    if ( ext_face_comm->itorecv )
      PMMG_DEL_MEM(parmesh,ext_face_comm->itorecv,int,"itorecv array");
    if ( ext_face_comm->itosend )
      PMMG_DEL_MEM(parmesh,ext_face_comm->itosend,int,"itosend array");
    if ( ext_face_comm->rtorecv )
      PMMG_DEL_MEM(parmesh,ext_face_comm->rtorecv,double,"rtorecv array");
    if ( ext_face_comm->rtosend )
      PMMG_DEL_MEM(parmesh,ext_face_comm->rtosend,double,"rtosend array");
  }
  PMMG_DEL_MEM(parmesh,hash.item,PMMG_hgrp,"group hash table");
  PMMG_DEL_MEM(parmesh,int_face_comm->intvalues,int,"face communicator");
  PMMG_DEL_MEM(parmesh,int_face_comm->doublevalues,double,"face communicator");

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
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    if ( ext_face_comm->rtorecv )
      PMMG_DEL_MEM(parmesh,ext_face_comm->rtorecv,double,"rtorecv array");
  }
  for ( k=0; k<parmesh->next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    if ( ext_face_comm->rtosend )
      PMMG_DEL_MEM(parmesh,ext_face_comm->rtosend,double,"rtosend array");
  }
  if ( int_face_comm->intvalues )
    PMMG_DEL_MEM(parmesh,int_face_comm->intvalues,int,"face communicator");
  if ( int_face_comm->doublevalues )
    PMMG_DEL_MEM(parmesh,int_face_comm->doublevalues,double,"face communicator");
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
  MMG5_pSol  met  = grp[0].met;
  idx_t      *xadj,*adjncy,*vwgt,*adjwgt;
  idx_t      adjsize;
  idx_t      nelt = mesh->ne;
  idx_t      ncon = 1; // number of balancing constraint
  idx_t      options[METIS_NOPTIONS];
  idx_t      objval = 0;
  int        ier = 0;
  int        status = 1;

  xadj = adjncy = vwgt = adjwgt = NULL;

  if (!nelt) return 1;

  /* Set contiguity of partitions if using Metis also for graph partitioning */
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_CONTIG] = ( parmesh->info.contiguous_mode &&
    (parmesh->info.loadbalancing_mode & PMMG_LOADBALANCING_metis) );

  /** Build the graph */
  if ( !PMMG_graph_meshElts2metis(parmesh,mesh,met,&xadj,&adjncy,&adjwgt,&adjsize) )
    return 0;


  /** Call metis and get the partition array */
  if( nproc >= 8 ) {
    ier = METIS_PartGraphKway( &nelt,&ncon,xadj,adjncy,vwgt,NULL,adjwgt,&nproc,
                               NULL,NULL,options,&objval, part );
  }
  else
    ier = METIS_PartGraphRecursive( &nelt,&ncon,xadj,adjncy,vwgt,NULL,adjwgt,&nproc,
                               NULL,NULL,options,&objval, part );
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

  PMMG_DEL_MEM(parmesh, adjwgt, idx_t, "deallocate adjwgt" );
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
int PMMG_part_parmeshGrps2metis( PMMG_pParMesh parmesh,idx_t* part,idx_t nproc )
{
  real_t     *tpwgts,*ubvec;
  idx_t      *xadj,*adjncy,*vwgt,*adjwgt,*vtxdist,adjsize;
  idx_t      *xadj_seq,*adjncy_seq,*vwgt_seq,*adjwgt_seq,*part_seq;
  idx_t      sendcounts,*recvcounts,*displs;
  idx_t      wgtflag,numflag;
  idx_t      ncon = 1; // number of balancing constraint
  idx_t      options[METIS_NOPTIONS];
  idx_t      objval = 0;
  int        ngrp,nprocs,ier;
  int        iproc,root,ip,status;

  ngrp   = parmesh->ngrp;
  nprocs = parmesh->nprocs;
  ier    = 1;

  /** Build the parmetis graph */
  xadj   = adjncy = vwgt = adjwgt = vtxdist = NULL;
  tpwgts = ubvec  =  NULL;

  if ( !PMMG_graph_parmeshGrps2parmetis(parmesh,&vtxdist,&xadj,&adjncy,&adjsize,
                                        &vwgt,&adjwgt,&wgtflag,&numflag,&ncon,
                                        nproc,&tpwgts,&ubvec) ) {
    fprintf(stderr,"\n  ## Error: Unable to build parmetis graph.\n");
    return 0;
  }

  /** Gather the graph on proc 0 */
  root = 0;
  vwgt_seq = adjwgt_seq = NULL;
  PMMG_CALLOC(parmesh,recvcounts,nproc,idx_t,"recvcounts", return 0);
  PMMG_CALLOC(parmesh,displs,nproc,idx_t,"displs", return 0);

  /** xadj, vwgt */
  for( iproc = 0; iproc<nproc; iproc++ ) {
    recvcounts[iproc] = vtxdist[iproc+1]-vtxdist[iproc];
    displs[iproc] = vtxdist[iproc];
  }

  xadj_seq = NULL;
  if(parmesh->myrank == root)
    PMMG_CALLOC(parmesh,xadj_seq,vtxdist[nproc]+1,idx_t,"xadj_seq", return 0);

  MPI_CHECK( MPI_Gatherv(&xadj[1],recvcounts[parmesh->myrank],MPI_INT,
                         &xadj_seq[1],recvcounts,displs,MPI_INT,
                         root,parmesh->comm), return 0);

  if(parmesh->myrank == root)
    for( iproc = 0; iproc < nproc; iproc++ )
      for( ip = 1; ip <= vtxdist[iproc+1]-vtxdist[iproc]; ip++ )
          xadj_seq[vtxdist[iproc]+ip] += xadj_seq[vtxdist[iproc]];

  if(wgtflag == PMMG_WGTFLAG_VTX || wgtflag == PMMG_WGTFLAG_BOTH ) {
    if(parmesh->myrank == root)
      PMMG_CALLOC(parmesh,vwgt_seq,vtxdist[nproc]+1,idx_t,"vwgt_seq", return 0);

    MPI_CHECK( MPI_Gatherv(vwgt,recvcounts[parmesh->myrank],MPI_INT,
                           vwgt_seq,recvcounts,displs,MPI_INT,
                           root,parmesh->comm), return 0);
  }

  /** adjncy, adjwgt */
  sendcounts = xadj[recvcounts[parmesh->myrank]];
  MPI_CHECK( MPI_Allgather(&sendcounts,1,MPI_INT,
                           recvcounts,1,MPI_INT,parmesh->comm), return 0);

  displs[0] = 0;
  for( iproc = 0; iproc<nproc-1; iproc++ ) {
    displs[iproc+1] = displs[iproc]+recvcounts[iproc];
  }

  adjncy_seq = NULL;
  if ( parmesh->myrank == root )
    PMMG_CALLOC(parmesh,adjncy_seq,xadj_seq[vtxdist[nproc]],idx_t,"xadj_seq", return 0);

  MPI_CHECK( MPI_Gatherv(adjncy,recvcounts[parmesh->myrank],MPI_INT,
                         adjncy_seq,recvcounts,displs,MPI_INT,
                         root,parmesh->comm), return 0);

  if(wgtflag == PMMG_WGTFLAG_ADJ || wgtflag == PMMG_WGTFLAG_BOTH ) {
    if(parmesh->myrank == root)
      PMMG_CALLOC(parmesh,adjwgt_seq,xadj_seq[vtxdist[nproc]],idx_t,"xadj_seq", return 0);

    MPI_CHECK( MPI_Gatherv(adjwgt,recvcounts[parmesh->myrank],MPI_INT,
                           adjwgt_seq,recvcounts,displs,MPI_INT,
                           root,parmesh->comm), return 0);
  }


  PMMG_DEL_MEM(parmesh,recvcounts,idx_t,"recvcounts");
  PMMG_DEL_MEM(parmesh,displs,idx_t,"displs");


  /** Call metis and get the partition array */
  if ( nprocs > 1 ) {

    if(parmesh->myrank == root) {
      PMMG_CALLOC(parmesh,part_seq,vtxdist[nproc],idx_t,"part_seq", return 0);


      /* Set contiguity of partitions */
      METIS_SetDefaultOptions(options);
      options[METIS_OPTION_CONTIG] = parmesh->info.contiguous_mode;

      /** Call metis and get the partition array */
      if( nprocs >= 8 )
        status = METIS_PartGraphKway( &vtxdist[nproc],&ncon,xadj_seq,adjncy_seq,
                                      vwgt_seq,NULL,adjwgt_seq,&nproc,
                                      NULL,NULL,options,&objval, part_seq );
      else
        status = METIS_PartGraphRecursive( &vtxdist[nproc],&ncon,xadj_seq,adjncy_seq,
                                      vwgt_seq,NULL,adjwgt_seq,&nproc,
                                      NULL,NULL,options,&objval, part_seq );
 

      if ( status != METIS_OK ) {
        switch ( status ) {
          case METIS_ERROR_INPUT:
            fprintf(stderr, "Group redistribution --- METIS_ERROR_INPUT: input data error\n" );
            break;
          case METIS_ERROR_MEMORY:
            fprintf(stderr, "Group redistribution --- METIS_ERROR_MEMORY: could not allocate memory error\n" );
            break;
          case METIS_ERROR:
            fprintf(stderr, "Group redistribution --- METIS_ERROR: generic error\n" );
            break;
          default:
            fprintf(stderr, "Group redistribution --- METIS_ERROR: update your METIS error handling\n" );
            break;
        }
        return 0;
      }
#ifndef NDEBUG
      /* Print graph to file */
/*      FILE* fid;
      char filename[48];
      snprintf(filename,17,"part_centralized");
      fid = fopen(filename,"w");
      for( iproc = 0; iproc < vtxdist[nproc]; iproc++ )
        fprintf(fid,"%d\n",part_seq[iproc]);
      fclose(fid);*/
#endif
    }

    /** Scatter the partition array */
    PMMG_CALLOC(parmesh,recvcounts,nproc,idx_t,"recvcounts", return 0);
    for( iproc = 0; iproc<nproc; iproc++ )
      recvcounts[iproc] = vtxdist[iproc+1]-vtxdist[iproc];
    assert(recvcounts[parmesh->myrank] == parmesh->ngrp);

    MPI_CHECK( MPI_Scatterv(part_seq,recvcounts,vtxdist,MPI_INT,
                            part,recvcounts[parmesh->myrank],MPI_INT,
                            root,parmesh->comm), return 0);
    PMMG_DEL_MEM(parmesh,recvcounts,idx_t,"recvcounts");

    /** Correct partitioning to avoid empty procs */
    if( !PMMG_correct_parmeshGrps2parmetis(parmesh,vtxdist,part,nproc) ) return 0;

    if(parmesh->myrank == root) PMMG_DEL_MEM(parmesh,part_seq,idx_t,"part_seq");

  }

  PMMG_DEL_MEM(parmesh, adjncy, idx_t, "deallocate adjncy" );
  PMMG_DEL_MEM(parmesh, xadj, idx_t, "deallocate xadj" );
  PMMG_DEL_MEM(parmesh, ubvec, real_t,"parmetis ubvec");
  PMMG_DEL_MEM(parmesh, tpwgts, real_t, "deallocate tpwgts" );
  PMMG_DEL_MEM(parmesh, vtxdist, idx_t, "deallocate vtxdist" );
  switch (wgtflag) {
    case PMMG_WGTFLAG_ADJ:
      PMMG_DEL_MEM(parmesh, adjwgt, idx_t, "deallocate adjwgt" );
      break;
    case PMMG_WGTFLAG_VTX:
      PMMG_DEL_MEM(parmesh, vwgt, idx_t, "deallocate vwgt" );
      break;
    case PMMG_WGTFLAG_BOTH:
      PMMG_DEL_MEM(parmesh, vwgt, idx_t, "deallocate vwgt" );
      PMMG_DEL_MEM(parmesh, adjwgt, idx_t, "deallocate adjwgt" );
      break;
    default:
      break;
  }

  if(parmesh->myrank == root) {
    PMMG_DEL_MEM(parmesh,xadj_seq,idx_t,"xadj_seq");
    PMMG_DEL_MEM(parmesh,adjncy_seq,idx_t,"adjcncy_seq");
    switch (wgtflag) {
     case PMMG_WGTFLAG_ADJ:
        PMMG_DEL_MEM(parmesh,adjwgt_seq,idx_t,"adjwgt_seq");
        break;
     case PMMG_WGTFLAG_VTX:
        PMMG_DEL_MEM(parmesh,vwgt_seq,idx_t,"vwgt_seq");
        break;
     case PMMG_WGTFLAG_BOTH:
        PMMG_DEL_MEM(parmesh,vwgt_seq,idx_t,"vwgt_seq");
        PMMG_DEL_MEM(parmesh,adjwgt_seq,idx_t,"adjwgt_seq");
        break;
      default:
        break;
    }
  }


  return ier;
}

#ifdef USE_PARMETIS
/**
 * \param parmesh pointer toward the parmesh structure
 * \param part pointer of an array containing the partitions (at the end)
 * \param nproc number of partitions asked
 *
 * \return  1 if success, 0 if fail
 *
 * Use parmetis to partition the first mesh in the list of meshes into nproc
 * groups
 *
 */
int PMMG_part_parmeshGrps2parmetis( PMMG_pParMesh parmesh,idx_t* part,idx_t nproc )
{
  real_t     *tpwgts,*ubvec;
  idx_t      *xadj,*adjncy,*vwgt,*adjwgt,*vtxdist,adjsize,edgecut;
  idx_t      wgtflag,numflag,ncon,options[3];
  int        ngrp,nprocs,ier;

  ngrp   = parmesh->ngrp;
  nprocs = parmesh->nprocs;
  ier    = 1;

  /** Build the parmetis graph */
  xadj   = adjncy = vwgt = adjwgt = vtxdist = NULL;
  tpwgts = ubvec  =  NULL;
  options[0] = 0;

  if ( !PMMG_graph_parmeshGrps2parmetis(parmesh,&vtxdist,&xadj,&adjncy,&adjsize,
                                        &vwgt,&adjwgt,&wgtflag,&numflag,&ncon,
                                        nproc,&tpwgts,&ubvec) ) {
    fprintf(stderr,"\n  ## Error: Unable to build parmetis graph.\n");
    return 0;
  }

  /** Call parmetis and get the partition array */
  if ( 2 < nprocs + ngrp ) {
    if ( ParMETIS_V3_PartKway( vtxdist,xadj,adjncy,vwgt,adjwgt,&wgtflag,&numflag,
                               &ncon,&nproc,tpwgts,ubvec,options,&edgecut,part,
                               &parmesh->comm) != METIS_OK ) {
        fprintf(stderr,"\n  ## Error: Parmetis fails.\n" );
        ier = 0;
    }
  }

  /** Correct partitioning to avoid empty procs */
  if( !PMMG_correct_parmeshGrps2parmetis(parmesh,vtxdist,part,nproc) ) return 0;

  PMMG_DEL_MEM(parmesh, adjncy, idx_t, "deallocate adjncy" );
  PMMG_DEL_MEM(parmesh, xadj, idx_t, "deallocate xadj" );
  PMMG_DEL_MEM(parmesh, ubvec, real_t,"parmetis ubvec");
  PMMG_DEL_MEM(parmesh, tpwgts, real_t, "deallocate tpwgts" );
  PMMG_DEL_MEM(parmesh, vtxdist, idx_t, "deallocate vtxdist" );
  switch (wgtflag) {
    case PMMG_WGTFLAG_ADJ:
      PMMG_DEL_MEM(parmesh, adjwgt, idx_t, "deallocate adjwgt" );
      break;
    case PMMG_WGTFLAG_VTX:
      PMMG_DEL_MEM(parmesh, vwgt, idx_t, "deallocate vwgt" );
      break;
    case PMMG_WGTFLAG_BOTH:
      PMMG_DEL_MEM(parmesh, vwgt, idx_t, "deallocate vwgt" );
      PMMG_DEL_MEM(parmesh, adjwgt, idx_t, "deallocate adjwgt" );
      break;
    default:
      break;
  }

  return ier;
}

#endif

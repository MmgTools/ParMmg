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

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param graph graph structure.
 *
 * Nullify graph arrays and variables.
 *
 */
void PMMG_graph_init( PMMG_pParMesh parmesh,PMMG_pGraph graph ) {
  graph->nvtxs   = 0;
  graph->nadjncy = 0;
  graph->npart   = 0;
  graph->vtxdist = NULL;
  graph->xadj    = NULL;
  graph->adjncy  = NULL;
  graph->vwgt    = NULL;
  graph->adjwgt  = NULL;
  graph->map     = NULL;
  graph->wgtflag = 0;
  graph->numflag = 0;
  graph->ncon    = 0;
  graph->tpwgts  = NULL;
  graph->ubvec   = NULL;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param graph_dst destination graph structure.
 * \param graph_src source graph structure.
 *
 * Copy graph variables (except for arrays).
 *
 */
void PMMG_graph_copy( PMMG_pParMesh parmesh,PMMG_pGraph graph_dst,
    PMMG_pGraph graph_src ) {
  graph_dst->nvtxs   = graph_src->nvtxs;
  graph_dst->nadjncy = graph_src->nadjncy;
  graph_dst->npart   = graph_src->npart;
  graph_dst->wgtflag = graph_src->wgtflag;
  graph_dst->numflag = graph_src->numflag;
  graph_dst->ncon    = graph_src->ncon;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param graph graph structure.
 *
 * Deallocate graph arrays.
 *
 */
void PMMG_graph_free( PMMG_pParMesh parmesh,PMMG_pGraph graph ) {
  if( graph->vtxdist )
    PMMG_DEL_MEM( parmesh, graph->vtxdist, idx_t,  "deallocate vtxdist" );
  if( graph->xadj )
    PMMG_DEL_MEM( parmesh, graph->xadj,    idx_t,  "deallocate xadj" );
  if( graph->adjncy )
    PMMG_DEL_MEM( parmesh, graph->adjncy,  idx_t,  "deallocate adjncy" );
  if( graph->vwgt )
    PMMG_DEL_MEM( parmesh, graph->vwgt,   idx_t, "deallocate vwgt" );
  if( graph->adjwgt )
    PMMG_DEL_MEM( parmesh, graph->adjwgt,   idx_t, "deallocate adjwgt" );
  if( graph->map )
    PMMG_DEL_MEM( parmesh, graph->map,      idx_t,  "deallocate map" );
  if( graph->tpwgts )
    PMMG_DEL_MEM( parmesh, graph->tpwgts,  real_t, "deallocate tpwgts" );
  if( graph->ubvec )
    PMMG_DEL_MEM( parmesh, graph->ubvec,   real_t, "deallocate ubvec");
}

void PMMG_subgraph_free( PMMG_pParMesh parmesh,PMMG_pGraph graph ) {
  PMMG_graph_free( parmesh,graph );
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param subgraph pointer to the subgraph structure.
 * \param nvtxs number of nodes in the subgraphs.
 * \param nadjncy number of arcs in the subgraph.
 * \param nvtxsmax maximum number of nodes in the subgraph.
 * \param nadjncymax maximum number of arcs in the subgraph.
 *
 * Allocate graph arrays and initialize variables.
 *
 */
int PMMG_subgraph_init( PMMG_pParMesh parmesh,PMMG_pGraph subgraph,
                        int nvtxsmax, int nadjncymax ) {
  subgraph->nvtxs   = 0;
  subgraph->nadjncy = 0;
  subgraph->npart   = 1;
  subgraph->wgtflag = 0;
  subgraph->numflag = PMMG_WGTFLAG_BOTH;
  subgraph->ncon    = 1;
  subgraph->vtxdist = NULL;
  PMMG_CALLOC( parmesh, subgraph->xadj, nvtxsmax+1, idx_t,
               "allocate xadj", return 0);
  PMMG_CALLOC( parmesh, subgraph->adjncy, nadjncymax+1, idx_t,
               "allocate adjncy", return 0);
  PMMG_CALLOC( parmesh, subgraph->vwgt, nvtxsmax+1, idx_t,
               "allocate vwgt", return 0);
  PMMG_CALLOC( parmesh, subgraph->adjwgt, nadjncymax+1, idx_t,
               "allocate adjwgt", return 0);
  PMMG_CALLOC( parmesh, subgraph->map, nvtxsmax, idx_t,
               "allocate map", return 0);
  subgraph->tpwgts  = NULL;
  subgraph->ubvec   = NULL;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param graph pointer toward the graph structure.
 * \param map pointer toward an integer map for graph nodes.
 * \param dim dimension for visualization (2D or 3D).
 * \param coordinates pointer toward a table of coordinates (for visualization).
 *
 * \return 1 if success, 0 if fail.
 *
 * Save graph in Medit format.
 *
 */
int PMMG_graph_save( PMMG_pParMesh parmesh,PMMG_pGraph graph,idx_t *map,int ndim,
                     double coords[][ndim], int *mask ) {
  FILE *fid;
  char name[48];
  int d;
  idx_t ivtx,jvtx,iadj;

  /* Open file and write headers */
  sprintf(name,"mygraph%d.mesh",parmesh->myrank);
  fid = fopen(name,"w");
  fprintf(fid,"MeshVersionFormatted 2\n");
  fprintf(fid,"\nDimension %d\n",ndim);

  /* Write vertices */
  fprintf(fid,"\nVertices\n%d\n",graph->nvtxs);
  for( ivtx = 0; ivtx < graph->nvtxs; ivtx++ ) {
    if( mask )
      jvtx = mask[ivtx];
    else
      jvtx = ivtx;
    for( d = 0; d < ndim; d++ )
      fprintf(fid,"%f ",coords[jvtx][d]);
    fprintf(fid,"0 \n");
  }

  /* Write edges */
  fprintf(fid,"\nEdges\n%d\n",graph->nadjncy);
  for( ivtx = 0; ivtx < graph->nvtxs; ivtx++ ) {
    for( iadj = graph->xadj[ivtx]; iadj < graph->xadj[ivtx+1]; iadj++ ) {
      jvtx = graph->adjncy[iadj];
      fprintf(fid,"%d %d %d\n",ivtx+1,jvtx+1,0);
    }
  }

  /* Close file */
  fprintf(fid,"\n\nEnd");
  fclose(fid);


  /* Open files and write headers */
  sprintf(name,"mygraph%d.sol",parmesh->myrank);
  fid = fopen(name,"w");
  fprintf(fid,"MeshVersionFormatted 2\n");
  fprintf(fid,"\nDimension 3\n");

  /* Write solution */
  fprintf(fid,"\nSolAtVertices\n%d\n %d %d\n",graph->nvtxs,1,MMG5_Scalar);
  for( ivtx = 0; ivtx < graph->nvtxs; ivtx++ ) {
    fprintf(fid,"%f\n",(double)map[ivtx]);
  }

  /* Close file */
  fprintf(fid,"\n\nEnd");
  fclose(fid);

  return 1;
}

int PMMG_graph_set( PMMG_pParMesh parmesh,PMMG_pGraph graph,
                    int nvtxs,int nadjncy,int *xadj,int *adjncy,
                    int *vwgt,int *adjwgt,int *vtxdist,int *map ){
  idx_t i;

  graph->nvtxs = nvtxs;
  graph->nadjncy = nadjncy;

  PMMG_CALLOC( parmesh,graph->xadj,graph->nvtxs+1,idx_t,"xadj",return 0 );
  for( i = 0; i < graph->nvtxs+1; i++ ) {
    graph->xadj[i] = xadj[i];
  }
  assert( graph->xadj[graph->nvtxs] == graph->nadjncy );

  PMMG_CALLOC( parmesh,graph->adjncy,graph->nadjncy+1,idx_t,"adjncy",return 0 );
  for( i = 0; i < graph->nadjncy+1; i++ ) {
    graph->adjncy[i] = adjncy[i];
  }
  PMMG_CALLOC( parmesh,graph->vtxdist,parmesh->nprocs+1,idx_t,"vtxdist",return 0 );
  for( i = 0; i < parmesh->nprocs+1; i++ ) {
    graph->vtxdist[i] = vtxdist[i];
  }
  if( vwgt ) {
    PMMG_CALLOC( parmesh,graph->vwgt,graph->nvtxs,idx_t,"vwgt",return 0 );
    for( i = 0; i < graph->nvtxs; i++ ) {
      graph->vwgt[i] = vwgt[i];
    }
  }
  if( adjwgt ) {
    PMMG_CALLOC( parmesh,graph->adjwgt,graph->nadjncy,idx_t,"adjwgt",return 0 );
    for( i = 0; i < graph->nadjncy; i++ ) {
      graph->adjwgt[i] = adjwgt[i];
    }
  }
  if( map ) {
    PMMG_CALLOC( parmesh,graph->map,graph->nvtxs,idx_t,"map",return 0 );
    for( i = 0; i < graph->nvtxs; i++ ) {
      graph->map[i] = map[i];
    }
  }

  return 1;
}

int PMMG_graph_test( PMMG_pParMesh parmesh,int nvtxs,int nadjncy,
                     int *xadj,int *adjncy,int *vtxdist,
                     int *color,int ndim, double coords[][ndim] ) {
  PMMG_graph graph;

  PMMG_graph_init( parmesh, &graph );
//  if( !PMMG_graph_set( parmesh, &graph, nvtxs, nadjncy, xadj, adjncy, vtxdist ) )
//    return 0;
  PMMG_graph_save( parmesh, &graph, color, ndim, coords, NULL );
  PMMG_graph_free( parmesh, &graph );
  return 1;
}

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

  PMMG_CALLOC(parmesh,sname,strlen(filename)+9,char,"file name prefix",return 0);
  PMMG_CALLOC(parmesh,smesh,strlen(filename)+15,char,"mesh file name",return 0);
  PMMG_CALLOC(parmesh,ssol,strlen(filename)+15,char,"sol file name",return 0);
  sprintf(sname,"%s-P%02d-I%02d",filename,parmesh->myrank,parmesh->iter);
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

int PMMG_bandGraph_grps( PMMG_pParMesh parmesh,
                         idx_t *vtxdist_in,idx_t *xadj_in,idx_t *adjncy_in,
                         idx_t **vtxdist_b,idx_t **xadj_b,idx_t **adjncy_b) {
  idx_t *vtx_tmp;
  int myrank = parmesh->myrank;
  int nvtx_in = vtxdist_in[myrank+1]-vtxdist_in[myrank];
  int nvtx_b,nadj_b;
  int i,k,iadj,kadj;

  PMMG_CALLOC(parmesh,vtx_tmp,nvtx_in,idx_t,"vtx_tmp",return 0);

  for( i = 0; i < nvtx_in; i++ )
    vtx_tmp[i] = PMMG_UNSET;

  /* Mark and count vertices adjacent to a parallel interface */
  nvtx_b = 0;
  for( i = 0; i < nvtx_in; i++ ) {
    for( k = 0; k < xadj_in[i+1]-xadj_in[i]; k++ ) {
      if( (adjncy_in[k] <  vtxdist_in[myrank]) ||
          (adjncy_in[k] >= vtxdist_in[myrank+1]) ) {
        vtx_tmp[i] = nvtx_b++;
        break;
      }
    }
  }

  PMMG_CALLOC(parmesh,*xadj_b,nvtx_b,idx_t,"xadj_b",return 0);

  /* Count adjacents in the band */
  for( i = 0; i < nvtx_in; i++ ) {
    if( vtx_tmp[i] == PMMG_UNSET ) continue;
    for( k = 0; k < xadj_in[i+1]-xadj_in[i]; k++ ) {
      if( (adjncy_in[k] <  vtxdist_in[myrank]) ||
          (adjncy_in[k] >= vtxdist_in[myrank+1]) ) {
        iadj = adjncy_in[k]-vtxdist_in[myrank];
        if( vtx_tmp[iadj] == PMMG_UNSET ) continue;
        (*xadj_b)[vtx_tmp[i]]++;
      }
    }
  }

  nadj_b = 0;
  for( i = 0; i < nvtx_b; i++ )
    nadj_b += (*xadj_b)[i];

  PMMG_CALLOC(parmesh,*adjncy_b,nadj_b,idx_t,"adjncy_b",return 0);

  /* Fill adjacency in the band */
  for( i = 0; i < nvtx_in; i++ ) {
    if( vtx_tmp[i] == PMMG_UNSET ) continue;
    kadj = 0;
    for( k = 0; k < xadj_in[i+1]-xadj_in[i]; k++ ) {
      if( (adjncy_in[k] <  vtxdist_in[myrank]) ||
          (adjncy_in[k] >= vtxdist_in[myrank+1]) ) {
        iadj = adjncy_in[k]-vtxdist_in[myrank];
        if( vtx_tmp[iadj] == PMMG_UNSET ) continue;
        (*adjncy_b)[kadj++] = iadj;
      }
    }
  }

  PMMG_DEL_MEM(parmesh,vtx_tmp,idx_t,"vtx_tmp");
  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param hash pointer toward the hatable.
 *
 * Reset the values in the hash table of group adjacency.
 *
 */
void PMMG_hashReset( PMMG_pParMesh parmesh,PMMG_HGrp *hash ) {
  int k;

  for (k=0; k<hash->siz; ++k ) {
    hash->item[k].adj = PMMG_UNSET;
    hash->item[k].wgt = PMMG_NUL;
  }

  for (k=hash->siz; k<hash->max; ++k) {
    hash->item[k].adj = PMMG_UNSET;
    hash->item[k].wgt = PMMG_NUL;
    hash->item[k].nxt = k+1;
  }
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param hash pointer toward the hatable.
 * \param hsiz initial size of hash table.
 * \param hmax maximal size of the hash table.
 *
 * \return 1 if success, 0 if fail.
 *
 * Initialisation of the hash table of group adjacency.
 *
 */
int PMMG_hashNew( PMMG_pParMesh parmesh,PMMG_HGrp *hash,int hsiz,int hmax ) {

  /* adjust hash table params */
  hash->siz  = hsiz+1;
  hash->max  = hmax + 2;
  hash->nxt  = hash->siz;

  PMMG_CALLOC(parmesh,hash->item,hash->max+1,PMMG_hgrp,"group hash table",return 0);

  PMMG_hashReset( parmesh,hash );

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
                               PMMG_pGraph graph ) {
  MMG5_pTetra  pt;
  int          *adja;
  int          j,k,iadr,jel,count,nbAdj,wgt,ier;

  /** Step 1: mesh adjacency creation */

  if ( (!mesh->adja) && (1 != MMG3D_hashTetra( mesh, 1 )) ) {
    fprintf( stderr,"  ## PMMG Hashing problem (1).\n" );
    return 0;
  }

  /** Step 2: build the metis graph */
  graph->nvtxs = mesh->ne;
  graph->ncon = 1;
  PMMG_CALLOC(parmesh, graph->xadj, graph->nvtxs+1, idx_t, "allocate xadj",
              return 0);

  /** 1) Count the number of adjacent of each elements and fill xadj */
  graph->xadj[0] = 0;
  graph->nadjncy = 0;
  for( k = 1; k <= mesh->ne; k++ ) {
    nbAdj = 0;
    iadr = 4*(k-1) + 1;
    adja = &mesh->adja[iadr];
    for( j = 0; j < 4; j++ )
      if( adja[j] )
        nbAdj++;

    graph->nadjncy += nbAdj;
    graph->xadj[k]  = graph->nadjncy;
  }

  /** 2) List the adjacent of each elts in adjncy */
  ier = 1;
  ++graph->nadjncy;
  PMMG_CALLOC(parmesh, graph->adjncy, graph->nadjncy, idx_t, "allocate adjncy", ier=0;);
  if( !ier ) {
    PMMG_DEL_MEM(parmesh, graph->xadj, idx_t, "deallocate xadj" );
    return ier;
  }
  /* Don't compute weights at mesh distribution, or if output load balancing is required at last iter */
  if( (parmesh->iter != PMMG_UNSET) &&
      ((parmesh->iter < parmesh->niter-1) || parmesh->info.nobalancing) ) {
    PMMG_CALLOC(parmesh, graph->adjwgt, graph->nadjncy, idx_t, "allocate adjwgt", ier=0;);
    if( !ier ) {
      PMMG_DEL_MEM(parmesh, graph->xadj,   idx_t, "deallocate xadj" );
      PMMG_DEL_MEM(parmesh, graph->adjncy, idx_t, "deallocate adjncy" );
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
      if( graph->adjwgt ) {
        /* Compute weight using face edge size */
        wgt = (int)PMMG_computeWgt(mesh,met,pt,j);

        graph->adjwgt[count] = MG_MAX(wgt,1);
      }

      graph->adjncy[count++] = jel-1;
    }
    assert( count == ( graph->xadj[k] ) );
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
int PMMG_graph_parmeshGrps2parmetis( PMMG_pParMesh parmesh,PMMG_pGraph graph ) {
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

  graph->nvtxs = parmesh->ngrp;
  graph->npart = parmesh->nprocs;
  if( (parmesh->iter == parmesh->niter-1) && !parmesh->info.nobalancing ) {
    /* Switch off weights for output load balancing */
    graph->wgtflag = PMMG_WGTFLAG_NONE;
  } else {
    /* Default weight choice for parmetis */
    graph->wgtflag = PMMG_WGTFLAG_DEF;
  }
  graph->numflag = 0; /* C-style numbering */
  graph->ncon    = 1; /* number of weight per metis node */

  comm   = parmesh->comm;
  grp    = parmesh->listgrp;
  myrank = parmesh->myrank;
  ngrp   = parmesh->ngrp;

  /** Step 1: Fill vtxdist array with the range of groups local to each
   * processor */
  PMMG_CALLOC(parmesh,graph->vtxdist,parmesh->nprocs+1,idx_t,"parmetis vtxdist", return 0);

  MPI_CHECK( MPI_Allgather(&graph->nvtxs,1,MPI_INT,&(graph->vtxdist)[1],1,MPI_INT,
             comm),goto fail_1 );

  for ( k=1; k<=parmesh->nprocs; ++k )
    (graph->vtxdist)[k] += (graph->vtxdist)[k-1];

  /** Step 2: Fill weights array with the number of MG_PARBDY face per group */
  PMMG_CALLOC(parmesh,graph->vwgt,graph->nvtxs,idx_t,"parmetis vwgt", goto fail_1);

  for ( igrp=0; igrp<ngrp; ++igrp ) {
    mesh = parmesh->listgrp[igrp].mesh;

    if ( !mesh ) {
      (graph->vwgt)[igrp] = 1;
      continue;
    }

    for ( k=1; k<=mesh->ne; ++k ) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) ) continue;

      (graph->vwgt)[igrp] += pt->mark;
    }
  }

  /* Fill tpwgts */
  PMMG_CALLOC(parmesh,graph->tpwgts,graph->ncon*graph->npart,real_t,"parmetis tpwgts", goto fail_2);
  for ( k=0; k < graph->ncon*graph->npart ; ++k )
    graph->tpwgts[k] = 1./(double)graph->npart;

  /* Fill ubvec */
  PMMG_CALLOC(parmesh,graph->ubvec,graph->ncon,real_t,"parmetis ubvec", goto fail_3);
  for ( k=0; k < graph->ncon; ++k )
    graph->ubvec[k] = PMMG_UBVEC_DEF;

  /** Step 3: Fill the internal communicator with the greater index of the 2
   * groups to which the face belong. Use a minus sign to mark old parallel
   * faces.*/
  PMMG_CALLOC(parmesh,graph->xadj,graph->nvtxs+1,idx_t,"parmetis xadj", goto fail_4);
  PMMG_CALLOC(parmesh,graph->map,graph->nvtxs,idx_t,"parmetis map", goto fail_5);

  int_face_comm = parmesh->int_face_comm;

  PMMG_MALLOC(parmesh,int_face_comm->intvalues,int_face_comm->nitem,int,
              "face communicator",goto fail_6);
  PMMG_CALLOC(parmesh,int_face_comm->doublevalues,int_face_comm->nitem,double,
              "face communicator",goto fail_6);

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

    PMMG_CALLOC(parmesh,ext_face_comm->itosend,nitem,int,"itosend array",
                goto fail_7);
    itosend = ext_face_comm->itosend;

    PMMG_CALLOC(parmesh,ext_face_comm->itorecv,nitem,int,"itorecv array",
                goto fail_7);
    itorecv       = ext_face_comm->itorecv;

    PMMG_CALLOC(parmesh,ext_face_comm->rtosend,nitem,double,"rtosend array",
                goto fail_7);
    rtosend = ext_face_comm->rtosend;

    PMMG_CALLOC(parmesh,ext_face_comm->rtorecv,nitem,double,"rtorecv array",
                goto fail_7);
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
                   comm,&status),goto fail_7 );
    MPI_CHECK(
      MPI_Sendrecv(rtosend,nitem,MPI_DOUBLE,color,MPI_PARMESHGRPS2PARMETIS_TAG,
                   rtorecv,nitem,MPI_DOUBLE,color,MPI_PARMESHGRPS2PARMETIS_TAG,
                   comm,&status),goto fail_7 );
  }

   /** Step 5: Process the external communicators to count for each group the
    * adjacent groups located on another processor and fill the sorted linked
    * list of adjacency */

  /* hash is used to store the sorted list of adjacent groups to a group */
  if ( !PMMG_hashNew(parmesh,&hash,ngrp+1,PMMG_NBADJA_GRPS*ngrp+1) )
    goto fail_7;

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
                           +graph->vtxdist[ext_face_comm->color_out],wgt);

      if ( !found ) {
        fprintf(stderr,"  ## Error: %s: unable to add a new group in adjacency"
                " hash table.\n",__func__);
        goto fail_8;
      }
      if ( found==2 ) continue; // The group is already in the adja list

      ++graph->xadj[ igrp+1 ];
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
      found = PMMG_hashGrp(parmesh,&hash,igrp,igrp_adj+graph->vtxdist[myrank],wgt);
      if ( !found ) {
        fprintf(stderr,"  ## Error: %s: unable to add a new group in adjacency"
                " hash table.\n",__func__);
        goto fail_8;
      }
      /* Count group only if not already in the adja list */
      if ( found!=2 ) ++graph->xadj[ igrp+1 ];

      /* Insert igrp in the sorted list of igrp_adj if not already found,
       * increment weight if found */
      found = PMMG_hashGrp(parmesh,&hash,igrp_adj,igrp+graph->vtxdist[myrank],wgt);
      if ( !found ) {
        fprintf(stderr,"  ## Error: %s: unable to add a new group in adjacency"
                " hash table.\n",__func__);
        goto fail_8;
      }
      /* Count group only if not already in the adja list */
      if ( found !=2 ) ++graph->xadj[ igrp_adj+1 ];
    }

    /* Store the active/inactive flag in the graph map */
    graph->map[igrp] = !grp->isNotActive;
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
    graph->xadj[igrp] += graph->xadj[igrp-1];

  /** Step 8: Fill adjncy array at metis format */
  PMMG_CALLOC(parmesh,graph->adjncy,graph->xadj[graph->nvtxs],idx_t,"adjcncy parmetis array",
              goto fail_8);
  PMMG_CALLOC(parmesh,graph->adjwgt,graph->xadj[graph->nvtxs],idx_t,"parmetis adjwgt",
              goto fail_8);

  graph->nadjncy = 0;
  for ( igrp=0; igrp<=ngrp; ++igrp ) {

    ph = &hash.item[igrp+1];
    if ( PMMG_UNSET==ph->adj ) continue;

    graph->adjncy[graph->nadjncy]   = ph->adj;
    graph->adjwgt[graph->nadjncy++] = MG_MAX(ph->wgt,1);

    while ( ph->nxt ) {
      ph                = &hash.item[ph->nxt];
      graph->adjncy[graph->nadjncy]   = ph->adj;
      graph->adjwgt[graph->nadjncy++] = MG_MAX(ph->wgt,1);
    }
  }
  assert ( graph->nadjncy==graph->xadj[graph->nvtxs] );

#ifndef NDEBUG
  /* Print graph to file */
/*  FILE* fid;
  char filename[48];
  sprintf(filename,"graph_proc%d",parmesh->myrank);
  fid = fopen(filename,"w");
  for( k = 0; k < parmesh->nprocs+1; k++ )
    fprintf(fid,"%d\n",(graph->vtxdist)[k]);
  fprintf(fid,"\n");
  for( k = 0; k < parmesh->ngrp+1; k++ )
    fprintf(fid,"%d\n",(graph->xadj)[k]);
  fprintf(fid,"\n");
  for( k = 0; k < (graph->xadj)[parmesh->ngrp]; k++ )
    fprintf(fid,"%d %d\n",(graph->adjncy)[k],(graph->adjwgt)[k]);
  fclose(fid);*/
#endif

  /* Nullify unnecessary weights */
  switch( graph->wgtflag ) {
    case PMMG_WGTFLAG_NONE:
      PMMG_DEL_MEM(parmesh,graph->vwgt,idx_t,"parmetis vwgt");
      PMMG_DEL_MEM(parmesh,graph->adjwgt,idx_t,"parmetis adjwgt");
      graph->vwgt = graph->adjwgt = NULL;
      break;
    case PMMG_WGTFLAG_ADJ:
      PMMG_DEL_MEM(parmesh,graph->vwgt,idx_t,"parmetis vwgt");
      graph->vwgt = NULL;
      break;
    case PMMG_WGTFLAG_VTX:
      PMMG_DEL_MEM(parmesh,graph->adjwgt,idx_t,"parmetis adjwgt");
      graph->adjwgt = NULL;
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

fail_8:
  if ( hash.item )
    PMMG_DEL_MEM(parmesh,hash.item,PMMG_hgrp,"group hash table");
fail_7:
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
fail_6:
  PMMG_DEL_MEM(parmesh,graph->map,idx_t,"parmetis map");
fail_5:
  PMMG_DEL_MEM(parmesh,graph->xadj,idx_t,"parmetis xadj");
fail_4:
  PMMG_DEL_MEM(parmesh,graph->ubvec,real_t,"parmetis ubvec");
fail_3:
  PMMG_DEL_MEM(parmesh,graph->tpwgts,real_t,"parmetis tpwgts");
fail_2:
  PMMG_DEL_MEM(parmesh,graph->vwgt,idx_t,"parmetis vwgt");
fail_1:
  PMMG_DEL_MEM(parmesh,graph->vtxdist,idx_t,"parmetis vtxdist");

  return 0;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param graph pointer toward the graph structure (centralized)
 * \param subgraph pointer toward the subgraph structure (centralized)
 * \param root rank of the process centralizing the graphs
 *
 * \return  1 if success, 0 if fail
 *
 * Create the map from nodes in the centralized graph to nodes in the active
 * subgraph, count nodes in the active subgraph and compute the number of
 * partitions.
 *
 */
int PMMG_subgraph_map_active( PMMG_pParMesh parmesh,PMMG_pGraph graph,
                              PMMG_pGraph subgraph,int root ) {
  idx_t  ivtx,sumwgts;
  int    iproc,igrp;

  /* return if not on the correct rank */
  if( parmesh->myrank != root ) return 1;

   /* Create map, count subgraph nodes and sum weights */
  sumwgts = 0;
  for( ivtx = 0; ivtx < graph->nvtxs; ivtx++ ) {
    if( graph->map[ivtx] ) {
      /* Count node */
      graph->map[ivtx] = subgraph->nvtxs++;
      /* Sum weight */
      sumwgts += graph->vwgt ? graph->vwgt[ivtx] : 1;
    } else {
      graph->map[ivtx] = PMMG_UNSET;
    }
  }

  /* Compute the number of parts for the active subgraph in order to fit on
   * the minimum possible number of processes */
  subgraph->npart = PMMG_howManyGroups( sumwgts,
                                        abs(parmesh->info.target_mesh_size) );

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param graph pointer toward the graph structure (centralized)
 * \param subgraph pointer toward the subgraph structure (centralized)
 * \param part partitioning array for the subgraph
 * \param root rank of the process centralizing the graphs
 *
 * \return  1 if success, 0 if fail
 *
 * Create the map from nodes in the centralized graph to nodes in the collapsed
 * subgraph, count nodes in the collapsed subgraph and set the number of
 * partitions.
 *
 */
int PMMG_subgraph_map_collapse( PMMG_pParMesh parmesh,PMMG_pGraph graph,
                                PMMG_pGraph subgraph,idx_t *part,
                                int root ) {
  idx_t ivtx;
  int   iinactive,ninactive;

  /* return if not on the correct rank */
  if( parmesh->myrank != root ) return 1;

  ninactive = 0;
  for( ivtx = 0; ivtx < graph->nvtxs; ivtx++ ) {
    if( graph->map[ivtx] != PMMG_UNSET )
      graph->map[ivtx] = part[graph->map[ivtx]];
    else
      ninactive++;
  }
  iinactive = 0;
  for( ivtx = 0; ivtx < graph->nvtxs; ivtx++ ) {
    if( (graph->map[ivtx] != PMMG_UNSET) && (graph->map[ivtx] < ninactive) )
      graph->map[ivtx] += ninactive;
    else
      graph->map[ivtx] = iinactive++;
  }

  /* Number of unique nodes in the collapsed graph */
  subgraph->nvtxs = ninactive+subgraph->npart;

  /* Number of partitions for the collapsed graph */
  subgraph->npart = parmesh->nprocs;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param graph pointer toward the graph structure (centralized)
 * \param part partitioning array for the centralized graph
 * \param part_sub partitioning array for the subgraph
 *
 * Remap the partitioning array from the collapsed subgraph to the centralized
 * graph.
 *
 */
void PMMG_graph_map_subpart( PMMG_pParMesh parmesh, PMMG_pGraph graph,
    idx_t *part, idx_t *part_sub ) {
  idx_t ivtx;

  for( ivtx = 0; ivtx < graph->nvtxs; ivtx++ ) {
      part[ivtx] = part_sub[graph->map[ivtx]];
    }

}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param graph pointer toward the graph structure (centralized)
 * \param subgraph pointer toward the subgraph structure (centralized)
 * \param hash pointer toward the array of lists structure for nodes adjacency
 *
 * \return  1 if success, 0 if fail
 *
 * Extract and build a subgraph from a graph (both centralized) given a map from
 * graph nodes to subgraph nodes. Use a temporary array of lists to fill nodes
 * adjacencies (already allocated and reset).
 *
 */
int PMMG_subgraph_build( PMMG_pParMesh parmesh,PMMG_pGraph graph,
                         PMMG_pGraph subgraph,PMMG_HGrp *hash ) {
  PMMG_hgrp *ph;
  idx_t      ivtx,jvtx,iadj,wgt;
  int        k,ier = 1;

  /* Reset number of arcs in the graph (the number of nodes has already been
   * reset) */
  subgraph->nadjncy = 0;

  /* Loop on nodes that have an image in the subgraph */
  for( ivtx = 0; ivtx < graph->nvtxs; ivtx++ ) {
    if( graph->map[ivtx] != PMMG_UNSET ) {
      /* Sum weights */
      wgt = graph->vwgt ? graph->vwgt[ivtx] : 1;
      subgraph->vwgt[graph->map[ivtx]] += wgt;
      /* Loop on active adjacents */
      for( iadj = graph->xadj[ivtx]; iadj < graph->xadj[ivtx+1]; iadj++ ) {
        jvtx = graph->adjncy[iadj];
        /* If the graph is collapsed, avoid edges onto the same node */
        if( (graph->map[jvtx] != PMMG_UNSET) &&
            (graph->map[jvtx] != graph->map[ivtx]) ) {
          /* Hash pair and sum dual weight */
          wgt = graph->adjwgt ? graph->adjwgt[iadj] : 1;
          ier = PMMG_hashGrp(parmesh,hash,graph->map[ivtx],graph->map[jvtx],wgt);
          if( !ier ) {
            return 0;
          } else if( ier == 1 ) {
            /* Count adjacent */
            subgraph->nadjncy++;
          }
        }
      }
      /* Fill inverse map frpm subgraph to graph */
      subgraph->map[graph->map[ivtx]] = ivtx;
    }
  }


  /* Loop on the adjacency of active nodes */
  for( k = 0; k < subgraph->nvtxs; k++ ) {

    /* Initialize the next free adjacency position as the last free position */
    subgraph->xadj[k+1] = subgraph->xadj[k];

    /* Store adjacents and dual weights */
    ph = &hash->item[k+1];
    if( ph->adj != PMMG_UNSET ) {
      subgraph->adjncy[subgraph->xadj[k+1]]   = ph->adj;
      subgraph->adjwgt[subgraph->xadj[k+1]++] = ph->wgt;
      while( ph->nxt ) {
        ph = &hash->item[ph->nxt];
        subgraph->adjncy[subgraph->xadj[k+1]]   = ph->adj;
        subgraph->adjwgt[subgraph->xadj[k+1]++] = ph->wgt;
      }
    }
  }
  assert( subgraph->nadjncy == subgraph->xadj[subgraph->nvtxs] );

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param part pointer of an array containing the tetrahedra partitions
 * \param npart pointer to the number of partitions asked
 * \param activelist array of flags for active/inactive parts
 *
 * \return  1 if success, 0 if fail
 *
 * Further partition the partitioned tetrahedra into active and inactive
 * subparts. Update the total number \var npart of parts.
 *
 */
int PMMG_part_meshElts_graded( PMMG_pParMesh parmesh, idx_t* part, idx_t *npart,
    int8_t *activelist ) {
  typedef struct {
    int head;
    idx_t id;
  } PMMG_listpart;
  PMMG_listpart *listpart;
  MMG5_pMesh     mesh;
  MMG5_pTetra    pt;
  MMG5_pPoint    ppt;
  idx_t          *mypart;
  int            ie,i,ip,ipart,inew;

  assert( parmesh->ngrp == 1 );
  mesh = parmesh->listgrp[0].mesh;

#ifdef USE_POINTMAP
  PMMG_CALLOC(parmesh,listpart,2*(*npart),PMMG_listpart,"listpart",return 0);

  /* Loop on tetra, update partitioning if a graded vertex is found */
  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    if( !MG_VOK(pt) ) continue;
    mypart = &part[ie-1];
    if( (*mypart) < (*npart) ) {
      /* Look for graded vertices */
      for( i = 0; i < 4; i++ ) {
        ip = pt->v[i];
        ppt = &mesh->point[ip];
        if( ppt->src < 0 ) {
          /* Update partitioning */
          (*mypart) += (*npart);
          /* Exit loop if a graded vertex is found */
          break;
        }
      }
    }
    /* Store tetra as list head if empty */
    if( !listpart[*mypart].head ) {
      listpart[*mypart].head = ie;
    }
  }

  /* Pack partitioning */
  inew = 0;
  for( ipart = 0; ipart < 2*(*npart); ipart++ ) {
    if( listpart[ipart].head ) {
      listpart[ipart].id = inew++;
      /* Flag active parts */
      if( ipart >= (*npart) ) {
        activelist[listpart[ipart].id] = 1;
      }
    } else {
      listpart[ipart].id = PMMG_UNSET;
    }
  }
  assert( inew <= 2*(*npart) );

  /* Loop on tetra and permute partitioning for all parts */
  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    if( !MG_VOK(pt) ) continue;
    mypart   = &part[ie-1];
    /* Permute tetra partitioning */
    assert( listpart[*mypart].id > PMMG_UNSET );
    (*mypart) = listpart[*mypart].id;
  }

  /* Update number of groups */
  (*npart) = inew;

  /* Free memory and return */
  PMMG_DEL_MEM(parmesh,listpart,PMMG_listpart,"listpart");
#endif

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param part pointer of an array containing the partitions (at the end)
 * \param npart number of partitions asked
 *
 * \return  1 if success, 0 if fail
 *
 * Use metis to partition the first mesh in the list of meshes into nproc groups
 *
 */
int PMMG_part_meshElts2metis( PMMG_pParMesh parmesh, idx_t* part, idx_t npart )
{
  PMMG_pGrp  grp = parmesh->listgrp;
  MMG5_pMesh mesh = grp[0].mesh;
  MMG5_pSol  met  = grp[0].met;
  PMMG_graph graph;
  idx_t      options[METIS_NOPTIONS];
  idx_t      objval = 0;
  int        ier = 0;
  int        status = 1;

  /* Set contiguity of partitions if using Metis also for graph partitioning */
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_CONTIG] = ( parmesh->info.contiguous_mode &&
    (parmesh->info.loadbalancing_mode & PMMG_LOADBALANCING_metis) );

  /** Build the graph */
  PMMG_graph_init( parmesh,&graph );
  graph.npart = npart;
  if ( !PMMG_graph_meshElts2metis(parmesh,mesh,met,&graph) )
    return 0;


  /** Call metis and get the partition array */
  if( npart >= 8 ) {
    ier = METIS_PartGraphKway( &graph.nvtxs,&graph.ncon,graph.xadj,graph.adjncy,
                               graph.vwgt,NULL,graph.adjwgt,&graph.npart,
                               NULL,NULL,options,&objval,part );
  }
  else
    ier = METIS_PartGraphRecursive( &graph.nvtxs,&graph.ncon,graph.xadj,
                                    graph.adjncy,graph.vwgt,NULL,graph.adjwgt,
                                    &graph.npart,NULL,NULL,options,&objval,part );
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
  if( !PMMG_correct_meshElts2metis( parmesh,part,graph.nvtxs,graph.npart ) )
    return 0;

  PMMG_graph_free( parmesh,&graph );

  return status;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param graph pointer toward the distributed graph structure
 * \param graph_seq pointer toward the centralized graph structure
 * \param root rank of the process to gather the graph onto.
 *
 * \return  1 if success, 0 if fail
 *
 * Gather a distributed graph on a root process.
 *
 */
int PMMG_graph_gather( PMMG_pParMesh parmesh, PMMG_pGraph graph,
                       PMMG_pGraph graph_seq, int root ) {
  idx_t sendcounts,*recvcounts,*displs;
  int   iproc,ip;

  PMMG_graph_init( parmesh, graph_seq );
  PMMG_graph_copy( parmesh, graph_seq, graph );

  /** vtxdist */
  if(parmesh->myrank == root) {
    PMMG_CALLOC(parmesh,graph_seq->vtxdist,parmesh->nprocs+1,idx_t,"vtxdist", return 0 );
    for( iproc = 0; iproc < parmesh->nprocs+1; iproc++ ) {
      graph_seq->vtxdist[iproc] = graph->vtxdist[iproc];
    }
    graph_seq->nvtxs = graph_seq->vtxdist[parmesh->nprocs];
  }

  PMMG_CALLOC( parmesh,recvcounts,parmesh->nprocs,idx_t,"recvcounts", return 0 );
  PMMG_CALLOC( parmesh,displs,    parmesh->nprocs,idx_t,"displs",     return 0 );

  /** xadj, vwgt */
  for( iproc = 0; iproc<parmesh->nprocs; iproc++ ) {
    recvcounts[iproc] = graph->vtxdist[iproc+1]-graph->vtxdist[iproc];
    displs[iproc]     = graph->vtxdist[iproc];
  }

  if(parmesh->myrank == root)
    PMMG_CALLOC(parmesh,graph_seq->xadj,graph->vtxdist[parmesh->nprocs]+1,idx_t,"xadj_seq", return 0 );

  MPI_CHECK( MPI_Gatherv(&graph->xadj[1],recvcounts[parmesh->myrank],MPI_INT,
                         &graph_seq->xadj[1],recvcounts,displs,MPI_INT,
                         root,parmesh->comm), return 0);

  if(parmesh->myrank == root)
    PMMG_CALLOC(parmesh,graph_seq->map,graph->vtxdist[parmesh->nprocs],idx_t,"map_seq", return 0 );

  MPI_CHECK( MPI_Gatherv(&graph->map[0],recvcounts[parmesh->myrank],MPI_INT,
                         &graph_seq->map[0],recvcounts,displs,MPI_INT,
                         root,parmesh->comm), return 0);

  if(parmesh->myrank == root) {
    for( iproc = 0; iproc < parmesh->nprocs; iproc++ )
      for( ip = 1; ip <= graph->vtxdist[iproc+1]-graph->vtxdist[iproc]; ip++ )
          graph_seq->xadj[graph->vtxdist[iproc]+ip] += graph_seq->xadj[graph->vtxdist[iproc]];
    graph_seq->nadjncy = graph_seq->xadj[graph->vtxdist[parmesh->nprocs]];
  }

  if(graph->wgtflag == PMMG_WGTFLAG_VTX || graph->wgtflag == PMMG_WGTFLAG_BOTH ) {
    if(parmesh->myrank == root)
      PMMG_CALLOC(parmesh,graph_seq->vwgt,graph->vtxdist[parmesh->nprocs]+1,idx_t,"vwgt_seq", return 0);

    MPI_CHECK( MPI_Gatherv(graph->vwgt,recvcounts[parmesh->myrank],MPI_INT,
                           graph_seq->vwgt,recvcounts,displs,MPI_INT,
                           root,parmesh->comm), return 0);
  }

  /** adjncy, adjwgt */
  sendcounts = graph->xadj[recvcounts[parmesh->myrank]];
  MPI_CHECK( MPI_Allgather(&sendcounts,1,MPI_INT,
                           recvcounts,1,MPI_INT,parmesh->comm), return 0);

  displs[0] = 0;
  for( iproc = 0; iproc<parmesh->nprocs-1; iproc++ ) {
    displs[iproc+1] = displs[iproc]+recvcounts[iproc];
  }

  if ( parmesh->myrank == root )
    PMMG_CALLOC(parmesh,graph_seq->adjncy,graph_seq->xadj[graph->vtxdist[parmesh->nprocs]],idx_t,"adjncy_seq", return 0);

  MPI_CHECK( MPI_Gatherv(graph->adjncy,recvcounts[parmesh->myrank],MPI_INT,
                         graph_seq->adjncy,recvcounts,displs,MPI_INT,
                         root,parmesh->comm), return 0);

  if(graph->wgtflag == PMMG_WGTFLAG_ADJ || graph->wgtflag == PMMG_WGTFLAG_BOTH ) {
    if(parmesh->myrank == root)
      PMMG_CALLOC(parmesh,graph_seq->adjwgt,graph_seq->xadj[graph->vtxdist[parmesh->nprocs]],idx_t,"adjwgt_seq", return 0);

    MPI_CHECK( MPI_Gatherv(graph->adjwgt,recvcounts[parmesh->myrank],MPI_INT,
                           graph_seq->adjwgt,recvcounts,displs,MPI_INT,
                           root,parmesh->comm), return 0);
  }


  PMMG_DEL_MEM(parmesh,recvcounts,idx_t,"recvcounts");
  PMMG_DEL_MEM(parmesh,displs,idx_t,"displs");

  return 1;
}

int PMMG_subgraph_part( PMMG_pParMesh parmesh, PMMG_pGraph graph, idx_t *part ){
  idx_t      options[METIS_NOPTIONS];
  idx_t      objval = 0;
  int        status;

  /** Call metis and get the partition array */
  if( graph->npart > 1 ) {

    /* Set contiguity of partitions */
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_CONTIG] = parmesh->info.contiguous_mode;

    /** Call metis and get the partition array */
    if( graph->npart >= 8 )
      status = METIS_PartGraphKway( &graph->nvtxs,&graph->ncon,
                                    graph->xadj,graph->adjncy,
                                    graph->vwgt,NULL,graph->adjwgt,
                                    &graph->npart,
                                    NULL,NULL,options,&objval, part );
    else
      status = METIS_PartGraphRecursive( &graph->nvtxs,&graph->ncon,
                                         graph->xadj,graph->adjncy,
                                         graph->vwgt,NULL,graph->adjwgt,
                                         &graph->npart,
                                         NULL,NULL,options,&objval, part );


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
  }

  return 1;
}

int PMMG_part_scatter( PMMG_pParMesh parmesh,PMMG_pGraph graph,int *part,
                       int *part_seq,int root ) {
  idx_t *recvcounts;
  int   iproc;

  /** Scatter the partition array */
  PMMG_CALLOC(parmesh,recvcounts,parmesh->nprocs,idx_t,"recvcounts", return 0);
  for( iproc = 0; iproc<parmesh->nprocs; iproc++ )
    recvcounts[iproc] = graph->vtxdist[iproc+1]-graph->vtxdist[iproc];
  assert(recvcounts[parmesh->myrank] == parmesh->ngrp);

  MPI_CHECK( MPI_Scatterv(part_seq,recvcounts,graph->vtxdist,MPI_INT,
                          part,recvcounts[parmesh->myrank],MPI_INT,
                          root,parmesh->comm), return 0);
  PMMG_DEL_MEM(parmesh,recvcounts,idx_t,"recvcounts");

  return 1;
}

int PMMG_part_active( PMMG_pParMesh parmesh, idx_t *part ) {
  PMMG_graph graph,graph_seq,subgraph;
  PMMG_HGrp  hash;
  idx_t      *map,*part_seq,*part_sub;
  idx_t      options[METIS_NOPTIONS];
  idx_t      objval = 0;
  int        ier;
  int        iproc,root,status,iinactive,ninactive,ivtx;

  ier    = 1;

  /** Step 0: Build the distributed graph */
  PMMG_graph_init( parmesh, &graph );

  if ( !PMMG_graph_parmeshGrps2parmetis( parmesh,&graph ) ) {
    fprintf(stderr,"\n  ## Error: Unable to build parmetis graph.\n");
    return 0;
  }

  /** Gather the graph on proc 0 */
  root = 0;

  if( !PMMG_graph_gather( parmesh, &graph, &graph_seq, root ) )
    return 0;

  if(parmesh->myrank == root) {

    /* Initialize subgraph and hash table */
    if( !PMMG_subgraph_init( parmesh, &subgraph, graph_seq.nvtxs, graph_seq.nadjncy ) )
      return 0;
    /* hash is used to store the sorted list of adjacent groups to a group */
    if ( !PMMG_hashNew( parmesh,&hash,graph_seq.nvtxs+1,
                        PMMG_NBADJA_GRPS*graph_seq.nvtxs+1) )
      return 0;

    /* Create map from centralized to reduced graph, compute number of nodes in
     * the reduced graph and the target number of parts */
    if( !PMMG_subgraph_map_active( parmesh, &graph_seq, &subgraph, root ) )
      return 0;

    /* Extract the reduced subgraph */
    if( !PMMG_subgraph_build( parmesh, &graph_seq, &subgraph, &hash ) )
      return 0;

    /* Allocate the partition array, sized with the maximum number
     * of active nodes (reuse it for the active subgraph and the collapsed graph */
    PMMG_CALLOC(parmesh,part_seq,graph_seq.nvtxs,idx_t,"part_seq", return 0);
    PMMG_CALLOC(parmesh,part_sub,graph_seq.nvtxs,idx_t,"part_sub", return 0);

    /* Partition the active subgraph */
    if( !PMMG_subgraph_part( parmesh, &subgraph, part_sub ) )
      return 0;

    /* Collapse map using the partition array, count inactive nodes */
    if( !PMMG_subgraph_map_collapse( parmesh,&graph_seq,&subgraph,part_sub,
                                     root ) )
      return 0;

    /* reset hash */
    PMMG_hashReset( parmesh,&hash );

    /* extract collapsed graph */
    if( !PMMG_subgraph_build( parmesh, &graph_seq, &subgraph, &hash ) )
      return 0;

    /* Partition the collapsed graph */
    if( !PMMG_subgraph_part( parmesh, &subgraph, part_sub ) )
      return 0;

    /* Update the original partition array */
    PMMG_graph_map_subpart( parmesh, &graph_seq, part_seq, part_sub );

  } /* end of sequential part */

  /** Scatter the partition array */
  if( !PMMG_part_scatter( parmesh, &graph, part, part_seq, root ) )
    return 0;

//  /** Correct partitioning to avoid empty procs */
//  if( !PMMG_correct_parmeshGrps2parmetis(parmesh,graph.vtxdist,part,
//      graph.npart) ) return 0;

  if( parmesh->myrank == root ) {
    PMMG_DEL_MEM( parmesh,part_seq,idx_t,"part_seq" );
    PMMG_subgraph_free( parmesh,&subgraph );
    PMMG_graph_free( parmesh,&graph_seq );
  }
  PMMG_graph_free( parmesh,&graph );


  return 1;
}


/**
 * \param parmesh pointer toward the parmesh structure
 * \param part pointer of an array containing the partitions (at the end)
 * \param npart number of partitions asked
 *
 * \return  1 if success, 0 if fail
 *
 * Use metis to partition the first mesh in the list of meshes into nproc groups
 *
 */
int PMMG_part_parmeshGrps2metis( PMMG_pParMesh parmesh,idx_t* part,idx_t npart )
{
  PMMG_graph graph,graph_seq;
  idx_t      *part_seq;
  idx_t      *recvcounts;
  idx_t      options[METIS_NOPTIONS];
  idx_t      objval = 0;
  int        ier;
  int        iproc,root,status;

  ier    = 1;

  /** Build the parmetis graph */
  PMMG_graph_init( parmesh, &graph );

  if ( !PMMG_graph_parmeshGrps2parmetis( parmesh,&graph ) ) {
    fprintf(stderr,"\n  ## Error: Unable to build parmetis graph.\n");
    return 0;
  }

  /* Store the asked number of partitions */
  graph.npart = npart;

  /** Gather the graph on proc 0 */
  root = 0;

  if( !PMMG_graph_gather( parmesh, &graph, &graph_seq, root ) )
    return 0;


  /** Call metis and get the partition array */
  if ( graph.npart > 1 ) {

    if(parmesh->myrank == root) {
      PMMG_CALLOC(parmesh,part_seq,graph.vtxdist[parmesh->nprocs],idx_t,"part_seq", return 0);


      /* Set contiguity of partitions */
      METIS_SetDefaultOptions(options);
      options[METIS_OPTION_CONTIG] = parmesh->info.contiguous_mode;

      /** Call metis and get the partition array */
      if( graph.npart >= 8 )
        status = METIS_PartGraphKway( &graph.vtxdist[parmesh->nprocs],&graph.ncon,
                                      graph_seq.xadj,graph_seq.adjncy,
                                      graph_seq.vwgt,NULL,graph_seq.adjwgt,
                                      &graph_seq.npart,
                                      NULL,NULL,options,&objval, part_seq );
      else
        status = METIS_PartGraphRecursive( &graph.vtxdist[parmesh->nprocs],&graph.ncon,
                                           graph_seq.xadj,graph_seq.adjncy,
                                           graph_seq.vwgt,NULL,graph_seq.adjwgt,
                                           &graph_seq.npart,
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
      sprintf(filename,"part_centralized");
      fid = fopen(filename,"w");
      for( iproc = 0; iproc < vtxdist[nproc]; iproc++ )
        fprintf(fid,"%d\n",part_seq[iproc]);
      fclose(fid);*/
#endif
    }

    /** Scatter the partition array */
    PMMG_CALLOC(parmesh,recvcounts,parmesh->nprocs,idx_t,"recvcounts", return 0);
    for( iproc = 0; iproc<parmesh->nprocs; iproc++ )
      recvcounts[iproc] = graph.vtxdist[iproc+1]-graph.vtxdist[iproc];
    assert(recvcounts[parmesh->myrank] == parmesh->ngrp);

    MPI_CHECK( MPI_Scatterv(part_seq,recvcounts,graph.vtxdist,MPI_INT,
                            part,recvcounts[parmesh->myrank],MPI_INT,
                            root,parmesh->comm), return 0);
    PMMG_DEL_MEM(parmesh,recvcounts,idx_t,"recvcounts");

    /** Correct partitioning to avoid empty procs */
    if( !PMMG_correct_parmeshGrps2parmetis(parmesh,graph.vtxdist,part,
          graph.npart) ) return 0;

    if(parmesh->myrank == root) PMMG_DEL_MEM(parmesh,part_seq,idx_t,"part_seq");

  }

  PMMG_graph_free( parmesh,&graph );
  PMMG_graph_free( parmesh,&graph_seq );

  return ier;
}

#ifdef USE_PARMETIS
/**
 * \param parmesh pointer toward the parmesh structure
 * \param part pointer of an array containing the partitions (at the end)
 * \param npart number of partitions asked
 *
 * \return  1 if success, 0 if fail
 *
 * Use parmetis to partition the first mesh in the list of meshes into nproc
 * groups
 * \warning Refactored, but not tested
 *
 */
int PMMG_part_parmeshGrps2parmetis( PMMG_pParMesh parmesh,idx_t* part,idx_t npart )
{
  PMMG_graph graph;
  idx_t      edgecut,options[3];
  int        ngrp,nprocs,ier;

  ngrp   = parmesh->ngrp;
  nprocs = parmesh->nprocs;
  ier    = 1;

  /** Build the parmetis graph */
  PMMG_graph_init( parmesh,&graph );
  graph.npart = npart;
  options[0] = 0;

  if ( !PMMG_graph_parmeshGrps2parmetis(parmesh,&graph) ) {
    fprintf(stderr,"\n  ## Error: Unable to build parmetis graph.\n");
    return 0;
  }

  /** Call parmetis and get the partition array */
  if ( 2 < nprocs + ngrp ) {
    if ( ParMETIS_V3_PartKway( graph.vtxdist,graph.xadj,graph.adjncy,graph.vwgt,
                               graph.adjwgt,&graph.wgtflag,&graph.numflag,
                               &graph.ncon,&graph.npart,graph.tpwgts,graph.ubvec,options,&edgecut,part,
                               &parmesh->comm) != METIS_OK ) {
        fprintf(stderr,"\n  ## Error: Parmetis fails.\n" );
        ier = 0;
    }
  }

  /** Correct partitioning to avoid empty procs */
  if( !PMMG_correct_parmeshGrps2parmetis(parmesh,graph.vtxdist,part,graph.npart) ) return 0;

  PMMG_graph_free( parmesh,&graph );

  return ier;
}

#endif

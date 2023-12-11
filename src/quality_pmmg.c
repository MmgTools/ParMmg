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

#include "parmmg.h"
#include <stddef.h>
#include "inlined_functions_3d_private.h"
#include "mmgexterns_private.h"

typedef struct {
  double min;
  int iel, iel_grp, cpu;
} min_iel_t;

static int PMMG_count_nodes_par(PMMG_pParMesh parmesh,PMMG_pGrp grp){
  MMG5_pTetra    pt;
  MMG5_pPoint    ppt;
  PMMG_pInt_comm int_node_comm;
  int            *intvalues,base;
  int            i,ip,idx,np,iel;

  grp->mesh->base++;
  base = grp->mesh->base;
  int_node_comm = parmesh->int_node_comm;
  intvalues = int_node_comm->intvalues;

  /** Reset flags */
  for( ip = 1; ip <= grp->mesh->np; ip++ ) {
    ppt = &grp->mesh->point[ip];
    ppt->flag = 0;
  }

  /** Initialize counter */
  np = 0;

  /* 1) Count points if not marked in the internal communicator,
   *    then mqrk them in the internal communicator and flag them. */
  for( i = 0; i < grp->nitem_int_node_comm; i++ ) {
    ip  = grp->node2int_node_comm_index1[i];
    idx = grp->node2int_node_comm_index2[i];
    if( !intvalues[idx] ) {
      intvalues[idx] = base;
      np++;
    }
    grp->mesh->point[ip].flag = base;
  }

  /* 2) Count all other points (touched by a tetra) */
  for( iel = 1; iel <= grp->mesh->ne; iel++ ) {
    pt = &grp->mesh->tetra[iel];
    if( !MG_EOK(pt) ) continue;
    for( i = 0; i < 4; i++ ) {
      ppt = &grp->mesh->point[pt->v[i]];
      if( !ppt->flag ) {
        ppt->flag = base;
        np++;
      }
    }
  }

  return np;
}

static void PMMG_min_iel_compute(void *in1, void* out1, int *len, MPI_Datatype *dptr )
{
  min_iel_t *in, *out;

  in  = (min_iel_t *)in1;
  out = (min_iel_t *)out1;

  int i;
  for ( i=0; i<*len; i++) {
    if ( in[ i ].min < out[ i ]. min ) {
      out[ i ].min = in[ i ].min;
      out[ i ].iel = in[ i ].iel;
      out[ i ].iel_grp = in[ i ].iel_grp;
      out[ i ].cpu = in[ i ].cpu;
    }
  }
}

typedef struct {
  double avlen,lmin,lmax;
  int    ned,amin,bmin,amax,bmax,nullEdge,hl[9];
  int    cpu_min,cpu_max;
} PMMG_lenStats;

static void PMMG_compute_lenStats( void* in1,void* out1,int *len, MPI_Datatype *dptr )
{
  PMMG_lenStats *in,*out;
  int i,j;

  in  = (PMMG_lenStats *)in1;
  out = (PMMG_lenStats *)out1;

  for ( i=0; i<*len; i++ ) {

    out[i].avlen    += in[i].avlen;
    out[i].ned      += in[i].ned;
    out[i].nullEdge += in[i].nullEdge;

    for ( j=0; j<9; ++j ) {
      out[i].hl[j] += in[i].hl[j];
    }


    if ( in[i].lmin < out[i].lmin ) {
      out[i].lmin    = in[i].lmin;
      out[i].amin    = in[i].amin;
      out[i].bmin    = in[i].bmin;
      out[i].amax    = in[i].amax;
      out[i].bmax    = in[i].bmax;
      out[i].cpu_min = in[i].cpu_min;
    }

    if ( in[i].lmax > out[i].lmax ) {
      out[i].lmax    = in[i].lmax;
      out[i].amax    = in[i].amax;
      out[i].bmax    = in[i].bmax;
      out[i].amax    = in[i].amax;
      out[i].bmax    = in[i].bmax;
      out[i].cpu_max = in[i].cpu_max;
    }

  }
}

/**
 * \param parmesh pointer to parmesh structure
 * \param opt PMMG_INQUA if called before the Mmg call, PMMG_OUTQUA otherwise
 * \param isCentral 1 for centralized mesh (no parallel communication), 0 for
 * distributed mesh
 * \param comm pointer toward the MPI communicator to use: when called before
 * the first mesh balancing (at preprocessing stage) we have to use the
 * read_comm communicator (i.e. the communicator used to provide the inputs).
 * For all ather calls, comm has to be the communicator to use for computations.
 *
 * \return 1 if success, 0 if fail;
 *
 * Print quality histogram among all group meshes and all processors
 */
int PMMG_qualhisto( PMMG_pParMesh parmesh, int opt, int isCentral, MPI_Comm comm )
{
  PMMG_pGrp    grp;
  PMMG_pInt_comm int_node_comm;
  PMMG_pExt_comm ext_node_comm;
  int          *intvalues;
  int          i, j, k, iel_grp;
  int          np_cur,ne_cur;
  int64_t      np, np_result, ne, ne_result;
  double       max, max_cur, max_result;
  double       avg, avg_cur, avg_result;
  double       min, min_cur;
  int          iel, iel_cur;
  int          good, good_cur,good_result,optimLES,optimLES_result;
  int          med, med_cur, med_result;
  int          his[PMMG_QUAL_HISSIZE],his_cur[PMMG_QUAL_HISSIZE],his_result[PMMG_QUAL_HISSIZE];
  int          nrid, nrid_cur, nrid_result,ier;
  MPI_Op       iel_min_op;
  MPI_Datatype mpi_iel_min_t;
  MPI_Datatype types[ PMMG_QUAL_MPISIZE ] = { MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT };
  min_iel_t    min_iel, min_iel_result = { DBL_MAX, 0, 0, 0 };
  MPI_Aint     disps[ PMMG_QUAL_MPISIZE ] = { offsetof( min_iel_t, min ),
                                              offsetof( min_iel_t, iel ),
                                              offsetof( min_iel_t, iel_grp ),
                                              offsetof( min_iel_t, cpu ) };
  int lens[ PMMG_QUAL_MPISIZE ]           = { 1, 1, 1, 1 };

  /* Calculate the quality values for local process */
  iel_grp = 0;
  np = 0;
  ne = 0;
  max = max_cur = DBL_MIN;
  avg = 0.;
  min = min_cur = DBL_MAX;
  iel = 0;
  good = 0;
  med = 0;
  grp = &parmesh->listgrp[0];
  optimLES = ( grp && grp->mesh ) ? grp->mesh->info.optimLES : 0;

  /* Reset node intvalues (in order to avoid counting parallel nodes twice) */
  int_node_comm = parmesh->int_node_comm;
  if( int_node_comm ) {
    PMMG_CALLOC( parmesh,int_node_comm->intvalues,int_node_comm->nitem,int,"intvalues",return 0);
    intvalues = int_node_comm->intvalues;

    /* Mark nodes not to be counted if the outer rank is lower than myrank */
    for( k = 0; k < parmesh->next_node_comm; k++ ) {
      ext_node_comm = &parmesh->ext_node_comm[k];
      if( parmesh->myrank > ext_node_comm->color_out ) continue;
      for( i = 0; i < ext_node_comm->nitem; i++ )
        intvalues[ext_node_comm->int_comm_index[i]] = 1;
    }
  }


  for ( i = 0; i < PMMG_QUAL_HISSIZE; ++i )
    his[ i ] = 0;

  nrid = 0;
  for ( i = 0; i < parmesh->ngrp; ++i ) {
    grp  = &parmesh->listgrp[ i ];

    nrid_cur = 0;

    if ( grp->mesh->info.optimLES ) {
      MMG3D_computeLESqua(grp->mesh,grp->met,&ne_cur,&max_cur,&avg_cur,&min_cur,
                          &iel_cur,&good_cur,&med_cur,his_cur,parmesh->info.imprim0);
    }
    else {
      if ( opt == PMMG_INQUA ) {
        MMG3D_computeInqua( grp->mesh, grp->met, &ne_cur, &max_cur, &avg_cur, &min_cur,
                             &iel_cur, &good_cur, &med_cur, his_cur,parmesh->info.imprim0 );
      }
      else {
        assert ( opt == PMMG_OUTQUA );
        MMG3D_computeOutqua( grp->mesh, grp->met, &ne_cur, &max_cur, &avg_cur, &min_cur,
                             &iel_cur, &good_cur, &med_cur, his_cur, &nrid_cur,parmesh->info.imprim0 );

      }
    }

    if( !int_node_comm )
      np_cur = grp->mesh->np;
    else
      np_cur = PMMG_count_nodes_par( parmesh,grp );
    np   += (int64_t)np_cur;
    ne   += (int64_t)ne_cur;
    avg  += avg_cur;
    med  += med_cur;
    good += good_cur;

    if ( max_cur > max )
      max = max_cur;

    if ( min_cur < min ) {
      min = min_cur;
      iel = iel_cur;
      iel_grp = i;
    }

    for ( j = 0; j < PMMG_QUAL_HISSIZE; ++j )
      his[ j ] += his_cur[ j ];

    nrid += nrid_cur;
  }

  if( int_node_comm )
    PMMG_DEL_MEM( parmesh,int_node_comm->intvalues,int,"intvalues" );

  if ( parmesh->info.imprim0 <= PMMG_VERB_VERSION )
    return 1;

  /* Calculate the quality values for all processes */
  if( isCentral ) {
    np_result = np;
    ne_result = ne;
    avg_result = avg;
    med_result = med;
    good_result = good;
    max_result = max;
    optimLES_result = optimLES;
  } else {
    MPI_Reduce( &np, &np_result, 1, MPI_INT64_T, MPI_SUM, 0, comm );
    MPI_Reduce( &ne, &ne_result, 1, MPI_INT64_T, MPI_SUM, 0, comm );
    MPI_Reduce( &avg, &avg_result, 1, MPI_DOUBLE, MPI_SUM, 0, comm );
    MPI_Reduce( &med, &med_result, 1, MPI_INT, MPI_SUM, 0, comm );
    MPI_Reduce( &good, &good_result, 1, MPI_INT, MPI_SUM, 0, comm );
    MPI_Reduce( &max, &max_result, 1, MPI_DOUBLE, MPI_MAX, 0, comm );
    MPI_Reduce( &optimLES,&optimLES_result,1,MPI_INT,MPI_MAX,0,comm );
  }

  min_iel.min = min;
  min_iel.iel = iel;
  min_iel.iel_grp = iel_grp;
  min_iel.cpu = parmesh->myrank;

  if( isCentral ) {
    min_iel_result.min = min_iel.min;
    min_iel_result.iel = min_iel.iel;
    min_iel_result.iel_grp = min_iel.iel_grp;
    min_iel_result.cpu = min_iel.cpu;
    for ( i = 0; i < PMMG_QUAL_HISSIZE; ++i )
      his_result[i] = his[i];
    nrid_result = nrid;
  } else {
    MPI_Type_create_struct( PMMG_QUAL_MPISIZE, lens, disps, types, &mpi_iel_min_t );
    MPI_Type_commit( &mpi_iel_min_t );

    MPI_Op_create( PMMG_min_iel_compute, 1, &iel_min_op );
    MPI_Reduce( &min_iel, &min_iel_result, 1, mpi_iel_min_t, iel_min_op, 0, comm );
    MPI_Reduce( his, his_result, PMMG_QUAL_HISSIZE, MPI_INT, MPI_SUM, 0, comm );
    MPI_Reduce( &nrid, &nrid_result, 1, MPI_INT, MPI_SUM, 0, comm );
    MPI_Type_free( &mpi_iel_min_t );
    MPI_Op_free( &iel_min_op );
  }


  if ( parmesh->myrank == 0 ) {

    if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
      fprintf(stdout,"\n  -- PARALLEL MESH QUALITY");

      if ( optimLES_result ) {
        fprintf( stdout," (LES)" );
      }

      fprintf( stdout, "  %"PRId64"   %"PRId64"\n", np_result, ne_result );

      fprintf( stdout, "     BEST   %8.6f  AVRG.   %8.6f  WRST.   %8.6f (",
               max_result, avg_result / ne_result, min_iel_result.min);

      if ( parmesh->ngrp>1 )
        fprintf( stdout, "GROUP %d - ",min_iel_result.iel_grp);

      if ( parmesh->nprocs>1 )
        fprintf( stdout, "PROC %d - ",min_iel_result.cpu);

      fprintf( stdout,"ELT %d)\n", min_iel_result.iel );
    }

    ier =
      MMG3D_displayQualHisto_internal( ne_result, max_result, avg_result,
                                       min_iel_result.min, min_iel_result.iel,
                                       good_result, med_result, his_result,
                                       nrid_result,optimLES_result,
                                       parmesh->info.imprim );
    if ( !ier ) return 0;
  }

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric structure.
 * \param avlen average length (to fill).
 * \param lmin minimal length (to fill).
 * \param lmax max length (to fill).
 * \param ned number of edges (to fill).
 * \param amin (to fill).
 * \param bmin (to fill).
 * \param amax (to fill).
 * \param bmax (to fill).
 * \param nullEdge (to fill).
 * \param metRidTyp (to fill).
 * \param bd_in (to fill).
 * \param hl (to fill).
 * \param comm pointer toward the MPI communicator to use: when called before
 * the first mesh balancing (at preprocessing stage) we have to use the
 * read_comm communicator (i.e. the communicator used to provide the inputs).
 * For all ather calls, comm has to be the communicator to use for computations.
 *
 * \return 2 without metric, 0 if fail, 1 otherwise.
 *
 * Compute the required information to print the length histogram
 *
 */
int PMMG_computePrilen( PMMG_pParMesh parmesh,MMG5_pMesh mesh, MMG5_pSol met, double* avlen,
                        double* lmin, double* lmax, int* ned, int* amin, int* bmin, int* amax,
                        int* bmax, int* nullEdge,int8_t metRidTyp, double** bd_in,
                        int hl[9],MPI_Comm comm )
{
  PMMG_pGrp       grp;
  PMMG_pInt_comm  int_edge_comm;
  PMMG_pExt_comm  ext_edge_comm;
  MMG5_pTetra     pt;
  MMG5_pPoint     ppt;
  MMG5_HGeom      hpar;
  MMG5_Hash       hash;
  int             *intvalues,idx;
  double          len;
  int             i,k,ia,np,nq,n;
  int             ref;
  int16_t         tag;
  int8_t          i0,i1,ier;
  static double   bd[9]= {0.0, 0.3, 0.6, 0.7071, 0.9, 1.3, 1.4142, 2.0, 5.0};

  *bd_in = bd;
  memset(hl,0,9*sizeof(int));
  *ned = 0;
  *avlen = 0.0;
  *lmax = 0.0;
  *lmin = 1.e30;
  *amin = *amax = *bmin = *bmax = 0;
  *nullEdge = 0;

  if ( (!met) || (!met->m) ) {
    /* the functions that computes the edge length cannot be called without an
     * allocated metric */
    return 2;
  }

  /* Hash parallel edges in the mesh */
  if ( PMMG_hashPar(mesh,&hpar) != PMMG_SUCCESS ) return 0;

  /* Build parallel edge communicator */
  if( !PMMG_build_edgeComm( parmesh,mesh,&hpar,comm ) ) return 0;

  /* Initialize internal communicator with current rank */
  int_edge_comm = parmesh->int_edge_comm;
  PMMG_MALLOC(parmesh,int_edge_comm->intvalues,int_edge_comm->nitem,int,"intvalues",return 0);
  intvalues = int_edge_comm->intvalues;
  for( i = 0; i < int_edge_comm->nitem; i++ )
    intvalues[i] = parmesh->myrank;

  /* Assign parallel edge to the lowest rank */
  for( k = 0; k < parmesh->next_edge_comm; k++ ) {
    ext_edge_comm = &parmesh->ext_edge_comm[k];
    for( i = 0; i < ext_edge_comm->nitem; i++ ) {
      idx = ext_edge_comm->int_comm_index[i];
      if( ext_edge_comm->color_out < intvalues[idx] )
        intvalues[idx] = ext_edge_comm->color_out;
    }
  }

  /* Hash all edges in the mesh */
  if ( !MMG5_hashNew(mesh,&hash,mesh->np,7*mesh->np) )  return 0;

  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    for(ia=0; ia<6; ia++) {
      i0 = MMG5_iare[ia][0];
      i1 = MMG5_iare[ia][1];
      np = pt->v[i0];
      nq = pt->v[i1];

      if(!MMG5_hashEdge(mesh,&hash,np,nq,0)){
        fprintf(stderr,"  ## Error: %s: function MMG5_hashEdge return 0\n",
                __func__);
        return 0;
      }
    }
  }


  /** Pop edges from hash table, and analyze their length */

  /* 1) On parallel edges */
  grp = &parmesh->listgrp[0];
  for( i = 0; i < grp->nitem_int_edge_comm; i++ ) {
    ia  = grp->edge2int_edge_comm_index1[i];
    idx = grp->edge2int_edge_comm_index2[i];

    /* Only analyse owned edges */
    if( intvalues[idx] != parmesh->myrank ) continue;

    np = mesh->edge[ia].a;
    nq = mesh->edge[ia].b;
    if( MMG5_hGet(&hpar,np,nq,&ref,&tag) ) {
      assert( tag & MG_BDY );

      /* Remove edge from hash ; ier = 1 if edge has been found */
      ier = MMG5_hashPop(&hash,np,nq);
      if( ier ) {
        assert ( met->m );
        if ( (!metRidTyp) && met->size==6 ) {
          assert ( met->m );
          /* We pass here if metric is aniso without metRidTyp */
          len = MMG5_lenSurfEdg33_ani(mesh,met,np,nq,(tag & MG_GEO));
        }
        else {
          /* We pass here if metric is aniso with metRidTyp or iso with
           * allocated metric. Note that the lenSurfEdg function segfault if called
           * with met==NULL or met->m==NULL */
          len = MMG5_lenSurfEdg(mesh,met,np,nq,0);
        }


        if ( !len ) {
          ++(*nullEdge);
        }
        else {
          *avlen += len;
          (*ned)++;

          if( len < (*lmin) ) {
            *lmin = len;
            *amin = np;
            *bmin = nq;
          }

          if ( len > (*lmax) ) {
            *lmax = len;
            *amax = np;
            *bmax = nq;
          }

          /* Locate size of edge among given table */
          for(i=0; i<8; i++) {
            if ( bd[i] <= len && len < bd[i+1] ) {
              hl[i]++;
              break;
            }
          }
          if( i == 8 ) hl[8]++;
        }
      }
    }

    /* Mark edge as analysed */
    intvalues[idx] = parmesh->nprocs;
  }


  /* 2) On internal edges */
  for(k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;
    n = 0;
    for(i=0 ; i<4 ; i++) {
      ppt = &mesh->point[pt->v[i]];
      if(!(MG_SIN(ppt->tag) || MG_NOM & ppt->tag) && (ppt->tag & MG_GEO)) continue;
      n++;
    }
    if(!n) {
      continue;
    }
    for(ia=0; ia<6; ia++) {
      i0 = MMG5_iare[ia][0];
      i1 = MMG5_iare[ia][1];
      np = pt->v[i0];
      nq = pt->v[i1];

      /* Remove edge from hash ; ier = 1 if edge has been found */
      ier = MMG5_hashPop(&hash,np,nq);
      if( ier ) {
        assert ( met->m );
        if ( (!metRidTyp) && met->size==6 && met->m ) {
          len = MMG5_lenedg33_ani(mesh,met,ia,pt);
        }
        else
          len = MMG5_lenedg(mesh,met,ia,pt);


        if ( !len ) {
          ++(*nullEdge);
        }
        else {
          *avlen += len;
          (*ned)++;

          if( len < (*lmin) ) {
            *lmin = len;
            *amin = np;
            *bmin = nq;
          }

          if ( len > (*lmax) ) {
            *lmax = len;
            *amax = np;
            *bmax = nq;
          }

          /* Locate size of edge among given table */
          for(i=0; i<8; i++) {
            if ( bd[i] <= len && len < bd[i+1] ) {
              hl[i]++;
              break;
            }
          }
          if( i == 8 ) hl[8]++;
        }
      }
    }
  }

  PMMG_DEL_MEM(parmesh,int_edge_comm->intvalues,int,"intvalues");
  PMMG_edge_comm_free( parmesh );
  MMG5_DEL_MEM(mesh,hpar.geom);
  MMG5_DEL_MEM(mesh,hash.item);
  MMG5_DEL_MEM(mesh,mesh->edge);
  mesh->na = 0;

  return 1;
}

/**
 * \param parmesh pointer to parmesh structure
 * \param metRidTyp Type of storage of ridges metrics: 0 for classic storage,
 * \param isCentral 1 for centralized mesh, 0 for distributed mesh.
 * \param comm pointer toward the MPI communicator to use: when called before
 * the first mesh balancing (at preprocessing stage) we have to use the
 * read_comm communicator (i.e. the communicator used to provide the inputs).
 * For all ather calls, comm has to be the communicator to use for computations.
 *
 * \return 1 if success, 0 if fail;
 *
 * Resume edge length histo computed on each procs on the root processor
 *
 * \warning to rewrite : the interface edges are counted multiple times in the
 * final histo
 *
 * \warning for now, only callable on "merged" parmeshes (=1 group per parmesh)
 *
 */
int PMMG_prilen( PMMG_pParMesh parmesh, int8_t metRidTyp, int isCentral, MPI_Comm comm )
{
  MMG5_pMesh    mesh;
  MMG5_pSol     met;
  double        dned,*bd;
  int           ier,ieresult;
  PMMG_lenStats lenStats,lenStats_result;
  MPI_Op        mpi_lenStats_op;
  MPI_Datatype  mpi_lenStats_t;
  MPI_Datatype  types[ PMMG_LENSTATS_MPISIZE ] = { MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,
                                                   MPI_INT, MPI_INT, MPI_INT,
                                                   MPI_INT, MPI_INT, MPI_INT,
                                                   MPI_INT, MPI_INT, MPI_INT };
  MPI_Aint      disps[ PMMG_LENSTATS_MPISIZE ] = { offsetof( PMMG_lenStats, avlen ),
                                                   offsetof( PMMG_lenStats, lmin  ),
                                                   offsetof( PMMG_lenStats, lmax  ),
                                                   offsetof( PMMG_lenStats, ned   ),
                                                   offsetof( PMMG_lenStats, amin  ),
                                                   offsetof( PMMG_lenStats, bmin  ),
                                                   offsetof( PMMG_lenStats, amax  ),
                                                   offsetof( PMMG_lenStats, bmax  ),
                                                   offsetof( PMMG_lenStats, nullEdge),
                                                   offsetof( PMMG_lenStats, hl    ),
                                                   offsetof( PMMG_lenStats, cpu_min),
                                                   offsetof( PMMG_lenStats, cpu_max) };
  int lens[ PMMG_LENSTATS_MPISIZE ]            = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 9, 1, 1 };



  ier  = ieresult = 1;
  mesh = NULL;

  if ( parmesh->ngrp > 1 ) {
    printf("  ## Warning:%s: this function must be called with at most 1"
           "group per processor. Exit function.\n",__func__);
    ier = 0;
  }

  MPI_Type_create_struct( PMMG_LENSTATS_MPISIZE, lens, disps, types, &mpi_lenStats_t );
  MPI_Type_commit( &mpi_lenStats_t );
  MPI_Op_create( PMMG_compute_lenStats, 1, &mpi_lenStats_op );

  lenStats.avlen = 0.;
  lenStats.lmin = DBL_MAX;
  lenStats.lmax = 0.;
  lenStats.ned = 0;
  lenStats.amin = lenStats.amax = lenStats.bmin = lenStats.bmax = 0;
  lenStats.nullEdge = 0;
  memset(lenStats.hl,0,9*sizeof(int));
  lenStats.cpu_min = lenStats.cpu_max = parmesh->myrank;

  if ( parmesh->ngrp==1 ) {
    mesh = parmesh->listgrp[0].mesh;
    met = parmesh->listgrp[0].met;
    if( isCentral ) {
      /* If metric is not allocated or if hash table alloc fails, the next
       * function returns 0, which allows to detect that we cannot print the
       * edge length histo. */
      ier = MMG3D_computePrilen( mesh, met,
                                 &lenStats.avlen, &lenStats.lmin,
                                 &lenStats.lmax, &lenStats.ned, &lenStats.amin,
                                 &lenStats.bmin, &lenStats.amax, &lenStats.bmax,
                                 &lenStats.nullEdge, metRidTyp, &bd,
                                 lenStats.hl );
    }
    else {
      /* The next function returns 0 if the hash table alloc fails and 2 if
       * called without metric (in this case we are not able to compute the edge
       * lengths). It allows to detect:
       *   - if we can't print the histo due to an alloc error (ier = 0 on 1 MPI
       *     process at least)
       *   - if we can't print the histo because the metric is not allocated
       *     (ier=2 on all the MPI process)
       *   - if we can print the histo (metric is allocated on at least 1 MPI
       *     process and no MPI process fail, thus ier is at least 1 on all the
       *     MPI proc but may be 2 on some of them)
       */
      ier = PMMG_computePrilen( parmesh, mesh, met,
                                &lenStats.avlen, &lenStats.lmin,
                                &lenStats.lmax, &lenStats.ned, &lenStats.amin,
                                &lenStats.bmin, &lenStats.amax, &lenStats.bmax,
                                &lenStats.nullEdge, metRidTyp, &bd,
                                lenStats.hl,comm );
    }
  }

  if( isCentral )
    ieresult = ier;
  else
    MPI_Reduce( &ier, &ieresult,1, MPI_INT, MPI_MIN, parmesh->info.root, comm );

  if ( (ieresult==0) || ieresult==2 ) {
    /* We are not able to print the histogram */
    MPI_Type_free( &mpi_lenStats_t );
    MPI_Op_free( &mpi_lenStats_op );
    return 0;
  }

  if( isCentral )
    memcpy(&lenStats_result,&lenStats,sizeof(PMMG_lenStats));
  else
    MPI_Reduce( &lenStats, &lenStats_result, 1, mpi_lenStats_t, mpi_lenStats_op,
                0, comm );

  MPI_Type_free( &mpi_lenStats_t );
  MPI_Op_free( &mpi_lenStats_op );

  if ( parmesh->myrank == parmesh->info.root ) {
    dned                  = (double)lenStats_result.ned;
    lenStats_result.avlen = lenStats_result.avlen / dned;

    fprintf(stdout,"\n  -- RESULTING EDGE LENGTHS (ROUGH EVAL.) %d \n",lenStats_result.ned);
    fprintf(stdout,"     AVERAGE LENGTH         %12.4f\n",lenStats_result.avlen);
    fprintf(stdout,"     SMALLEST EDGE LENGTH   %12.4f   %6d %6d",
            lenStats_result.lmin,lenStats_result.amin,lenStats_result.bmin);
    if ( parmesh->nprocs>1 ) {
      fprintf(stdout," (PROC %d)\n",lenStats_result.cpu_min);
    }
    else { fprintf(stdout,"\n"); }

    fprintf(stdout,"     LARGEST  EDGE LENGTH   %12.4f   %6d %6d",
            lenStats_result.lmax,lenStats_result.amax,lenStats_result.bmax);
    if ( parmesh->nprocs>1 ) {
      fprintf(stdout," (PROC %d)\n",lenStats_result.cpu_max);
    }
    else { fprintf(stdout,"\n"); }


    MMG5_displayLengthHisto_internal ( lenStats_result.ned,lenStats_result.amin,lenStats_result.bmin,
                                       lenStats_result.lmin,lenStats_result.amax,
                                       lenStats_result.bmax,lenStats_result.lmax,
                                       lenStats_result.nullEdge,bd,
                                       lenStats_result.hl,1,parmesh->info.imprim);
  }

  return 1;
}

/**
 * \param parmesh pointer to parmesh structure
 * \param metRidTyp Type of storage of ridges metrics: 0 for classic storage,
 *
 * \return 1 if success, 0 if fail;
 *
 * Compute elements quality in a given metrics.
 *
 */
int PMMG_tetraQual( PMMG_pParMesh parmesh,int8_t metRidTyp ) {
  PMMG_pGrp grp;
  int       igrp;

  for( igrp = 0; igrp < parmesh->ngrp; igrp++ ){
    grp  = &parmesh->listgrp[igrp];
    if( !MMG3D_tetraQual( grp->mesh, grp->met, metRidTyp ) ) {
      fprintf(stderr,"\n  ## Quality computation problem.\n");
      return 0;
    }
  }

  return 1;
}

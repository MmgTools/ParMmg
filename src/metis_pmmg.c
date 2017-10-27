/**
 * \file metis_pmmg.c
 * \brief Partition mesh using metis
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "metis_pmmg.h"
#include "mpitypes_pmmg.h"

/**
 * \param parmesh pointer toward the PMMG parmesh structure
 * \param mesh pointer toward a MMG5 mesh structure
 * \param xadj pointer toward the position of the elt adjacents in adjncy
 * \param adjncy pointer toward the list of the adjacent of each elt
 * \param nadjncy number of data in adjncy array
 * \param nproc number of partitions asked
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
                               idx_t *nadjncy,idx_t nproc ) {
  int        *adja;
  int        j,k,iadr,jel,count,nbAdj;

  /** Step 1: mesh adjacency creation */
  if ( (!mesh->adja) && (1 != MMG3D_hashTetra( mesh, 1 )) ) {
    fprintf( stderr,"  ## PMMG Hashing problem (1).\n" );
    return 0;
  }

  /** Step 2: build the metis graph */
  PMMG_CALLOC(parmesh, (*xadj), mesh->ne+1, idx_t, "allocate xadj",
              goto end_mesh_adja);

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
  ++(*nadjncy);
  PMMG_CALLOC(parmesh, (*adjncy), (*nadjncy), idx_t, "allocate adjncy",
              goto end_xadj);

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

  return 1;

#warning IS MESH->ADJA DEALLOCATION CORRECT? IF I CAN REMOVE IT, THEN MULTIPLE RETURN PATHS ARE NOT REQUIRED
end_mesh_adja:
  PMMG_DEL_MEM(mesh,mesh->adja,4*mesh->nemax+5,int,"deallocate mesh->adja");
end_xadj:
  PMMG_DEL_MEM(parmesh, xadj, mesh->ne + 1, idx_t, "deallocate xadj" );

  return 0;
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

  *wgtflag = 2; /* Weights applies on vertices */
  *numflag = 0; /* C-style numbering */
  *ncon    = 1; /* number of weight per metis node */

  return 1;
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
  idx_t      *xadj = NULL;
  idx_t      *adjncy = NULL;
  idx_t      nelt = mesh->ne;
  idx_t      ncon = 1; // number of balancing constraint
  idx_t      objval = 0;
  int        adjsize;
  int        ier = 0;

  if ( !PMMG_graph_meshElts2metis(parmesh,mesh,&xadj,&adjncy,&adjsize,nproc) )
    return 0;

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
    goto end_adjncy;
  }

  PMMG_DEL_MEM(parmesh, adjncy, adjsize, idx_t, "deallocate adjncy" );
  PMMG_DEL_MEM(parmesh, xadj, mesh->ne + 1, idx_t, "deallocate xadj" );

  return 1;

end_adjncy:
  PMMG_DEL_MEM(parmesh, adjncy, adjsize, idx_t, "deallocate adjncy" );
#warning IS MESH->ADJA DEALLOCATION CORRECT? IF I CAN REMOVE IT, THEN MULTIPLE RETURN PATHS ARE NOT REQUIRED
  PMMG_DEL_MEM(mesh, mesh->adja, 4*mesh->nemax + 5, int, "deallocate mesh->adja");
  PMMG_DEL_MEM(parmesh, xadj, mesh->ne + 1, idx_t, "deallocate xadj" );

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
int PMMG_part_parmeshGrps2parmetis( PMMG_pParMesh parmesh,idx_t* part,idx_t nproc )
{
  real_t     *tpwgts,*ubvec;
  idx_t      *xadj,*adjncy,*vwgt,*vtxdist,adjsize;
  idx_t      wgtflag,numflag,ncon;
  int        ier;

  /** Build the parmetis graph */
  xadj   = adjncy = vwgt = vtxdist = NULL;
  tpwgts = ubvec  =  NULL;
  if ( !PMMG_graph_parmeshGrps2parmetis(parmesh,&vtxdist,&xadj,&adjncy,&adjsize,
                                        &vwgt,&wgtflag,&numflag,&ncon,
                                        nproc,&tpwgts,&ubvec) )
    return 0;

  /** Call parmetis and get the partition array */
  /* ier = ParMETIS_V3_PartKway( vtxdist,xadj,adjncy,vwgt,NULL,&wgtflag,&numflag,&ncon, */
  /*                             &nproc,tpwgts,ubvec,NULL,NULL,part,parmesh->comm); */
  /* if ( ier != METIS_OK ) { */
  /*   switch ( ier ) { */
  /*     case METIS_ERROR_INPUT: */
  /*       fprintf(stderr, "METIS_ERROR_INPUT: input data error\n" ); */
  /*       break; */
  /*     case METIS_ERROR_MEMORY: */
  /*       fprintf(stderr, "METIS_ERROR_MEMORY: could not allocate memory error\n" ); */
  /*       break; */
  /*     case METIS_ERROR: */
  /*       fprintf(stderr, "METIS_ERROR: generic error\n" ); */
  /*       break; */
  /*     default: */
  /*       fprintf(stderr, "METIS_ERROR: update your METIS error handling\n" ); */
  /*       break; */
  /*   } */
  /*   goto end_adjncy; */
  /* } */

  /* PMMG_DEL_MEM(parmesh, adjncy, adjsize, idx_t, "deallocate adjncy" ); */
  /* PMMG_DEL_MEM(parmesh, xadj, ngrps+1, idx_t, "deallocate xadj" ); */
  /* PMMG_DEL_MEM(parmesh, vwgt, ngrps+1, idx_t, "deallocate vwgt" ); */
  /* PMMG_DEL_MEM(parmesh, adjwgt, mesh->ne + 1, idx_t, "deallocate adjwgt" ); */
  /* PMMG_DEL_MEM(parmesh, vtxdist, nproc+1, idx_t, "deallocate adjwgt" ); */

  return 1;

/* end_adjncy: */
/*   PMMG_DEL_MEM(parmesh, adjncy, adjsize, idx_t, "deallocate adjncy" ); */
/*   PMMG_DEL_MEM(parmesh, xadj, mesh->ne + 1, idx_t, "deallocate xadj" ); */

/*   return 0; */
}

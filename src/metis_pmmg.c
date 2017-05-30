/**
 * \file metisfunctions.c
 * \brief Partionning and repartitionning functions (with metis calls).
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */

#include "metis_pmmg.h"
#include "mpitypes_pmmg.h"

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param xadj idx_t pointer.
 * \param adjncy idx_t pointer.
 *
 * \return PMMG_FAILURE
 *         PMMG_SUCCESS
 *
 * Allocate and fill the xadj and adjncy arrays to be used in metis calls
 */
int
PMMG_mesh2metis( PMMG_pParMesh parmesh, idx_t** xadj, idx_t** adjncy, int *adjsize )
{
  PMMG_pGrp  grp = NULL;
  MMG5_pMesh mesh = NULL;
  int        *adja = NULL;
  int        j,k,iadr,jel,count,totCount,nbAdj;

  grp = parmesh->listgrp;
  mesh = grp[0].mesh;

  /* create tetra adjacency */
  if ( (!mesh->adja) && !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"  ## PMMG Hashing problem (1).\n");
    return PMMG_FAILURE;
  }

  /*mesh -> graph*/
  PMMG_CALLOC(parmesh, *xadj, mesh->ne + 1, idx_t,
              "allocate xadj", return PMMG_FAILURE);

  /*count neighboors*/
  (*xadj)[0]  = 0;
  totCount = 0;
  for( k = 1; k <= mesh->ne; k++ ) {
    nbAdj = 0;
    iadr = 4*(k-1) + 1;
    adja = &mesh->adja[iadr];
    for( j = 0; j < 4; j++ )
      if( adja[j] )
        nbAdj++;
    assert( nbAdj <= 4 );
    totCount += nbAdj;
    (*xadj)[k] = totCount;
  }


  PMMG_CALLOC(parmesh, *adjncy, totCount + 1, idx_t, "allocate adjncy",
              PMMG_FREE(parmesh,*xadj,mesh->ne + 1, idx_t, "deallocate xadj");
              return PMMG_FAILURE);
  *adjsize = totCount + 1;

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

    if ( count != (*xadj)[k] )
      printf( "count %d %d %d\n", k, count, (*xadj)[k] );
    assert( count == ( (*xadj)[k] ) );
  }

  return PMMG_SUCCESS;
}

int PMMG_metis_wrapper( idx_t* nelt, idx_t* ncon, idx_t* xadj, idx_t* adjncy,
                        idx_t* nproc, idx_t *objval, idx_t* part )
{
  idx_t ier = 0;
  ier = METIS_PartGraphKway( nelt, ncon, xadj, adjncy, NULL/*vwgt*/,
                             NULL/*vsize*/, NULL/*adjwgt*/, nproc,
                             NULL/*tpwgts*/, NULL/*ubvec*/, NULL/*options*/,
                             objval, part );
  /* ier =  METIS_PartGraphRecursive(&nelt,&ncon,xadj,adjncy,NULL/\*vwgt*\/,NULL/\*vsize*\/, */
  /*                         NULL/\*adjwgt*\/,&nproc,NULL/\*tpwgts*\/, */
  /*                         NULL/\*ubvec*\/,NULL/\*options*\/,&objval,part); */
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
    return PMMG_FAILURE;
  }
  return PMMG_SUCCESS;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param part pointer of an array containing the partitions (at the end).
 * \return  PMMG_SUCCESS
 *          PMMG_FAILURE
 *
 * Proc 0 call metis to part the mesh. Then broadcasts the partition.
 */
int PMMG_metispartitioning( PMMG_pParMesh parmesh, idx_t* part )
{
  PMMG_pGrp  grp;
  MMG5_pMesh mesh;
  idx_t      *xadj,*adjncy;
  idx_t      nelt,ncon,nproc,objval;
  int adjsize = 0;

  grp = parmesh->listgrp;
  mesh = grp[0].mesh;
  if ( !parmesh->myrank ) {

    if ( PMMG_SUCCESS != PMMG_mesh2metis( parmesh, &xadj, &adjncy, &adjsize ) )
      return PMMG_FAILURE;

    /*call metis*/
    ncon       = 1;/*number of balancing constraint*/
    nelt       = mesh->ne;
    nproc      = parmesh->nprocs;
    if (   PMMG_metis_wrapper( &nelt, &ncon, xadj, adjncy, &nproc, &objval, part )
        != PMMG_SUCCESS )
      return PMMG_FAILURE;

    PMMG_FREE(parmesh, adjncy, adjsize, idx_t, "deallocate adjncy" );
    PMMG_FREE(parmesh, xadj, mesh->ne + 1, idx_t, "deallocate xadj" );
  }

  /* Send at each proc the graph */
  MPI_Bcast( &part[0], mesh->ne, MPI_INT, 0, parmesh->comm );

  return PMMG_SUCCESS;
}

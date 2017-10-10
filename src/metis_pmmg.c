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
 * \param parmesh pointer toward the parmesh structure
 * \param part pointer of an array containing the partitions (at the end)
 * \return  PMMG_SUCCESS
 *          PMMG_FAILURE
 *
 * Use metis to partition the first mesh in the list of meshes into nproc groups
 */
int PMMG_partition_metis( PMMG_pParMesh parmesh, idx_t* part, idx_t nproc )
{
  PMMG_pGrp  grp = parmesh->listgrp;
  MMG5_pMesh mesh = grp[0].mesh;
  idx_t      *xadj = NULL;
  idx_t      *adjncy = NULL;
  int        *adja = NULL;
  idx_t      nelt = mesh->ne;
#warning NIKOS TODO: experiment with number of balancing constraints
  idx_t      ncon = 1; // number of balancing constraint
  idx_t      objval = 0;
  int        adjsize = 0;
  int        ier = 0;
  int        j, k, iadr, jel, count, totCount, nbAdj;


  // Create the adjacency vectors that the call to metis requires

  if ( (!mesh->adja) && (1 != MMG3D_hashTetra( mesh, 1 )) ) {
    fprintf( stderr,"  ## PMMG Hashing problem (1).\n" );
    return PMMG_FAILURE;
  }

  // mesh to graph
  PMMG_CALLOC(parmesh, xadj, mesh->ne + 1, idx_t, "allocate xadj",
              goto end_mesh_adja);

  // count neighboors
  xadj[0]  = 0;
  totCount = 0;
  for( k = 1; k <= mesh->ne; k++ ) {
    nbAdj = 0;
    iadr = 4*(k-1) + 1;
    adja = &mesh->adja[iadr];
    for( j = 0; j < 4; j++ )
      if( adja[j] )
        nbAdj++;

    totCount += nbAdj;
    xadj[k] = totCount;
  }

  PMMG_CALLOC(parmesh, adjncy, totCount + 1, idx_t, "allocate adjncy",
              goto end_xadj);
  adjsize = totCount + 1;

  count = 0;
  for( k = 1; k <= mesh->ne; k++ ) {
    iadr = 4*(k-1) + 1;
    adja = &mesh->adja[iadr];
    for ( j = 0; j < 4; j++ ) {
      jel = adja[j] / 4;

      if ( !jel )
        continue;

      adjncy[count++] = jel-1;
    }

    if ( count != xadj[k] )
      printf( "count %d %d %d\n", k, count, xadj[k] );
    assert( count == ( xadj[k] ) );
  }


  // Adjacency vectors are ready, call metis
  ier = METIS_PartGraphKway( &nelt, &ncon, xadj, adjncy, NULL/*vwgt*/,
                             NULL/*vsize*/, NULL/*adjwgt*/, &nproc,
                             NULL/*tpwgts*/, NULL/*ubvec*/, NULL/*options*/,
                             &objval, part );
  //ier = METIS_PartGraphRecursive( &nelt, &ncon, xadj, adjncy, NULL/\*vwgt*\/,
  //                                NULL/\*vsize*\/, NULL/\*adjwgt*\/, &nproc,
  //                                NULL/\*tpwgts*\/, NULL/\*ubvec*\/,
  //                                NULL/\*options*\/, &objval, part );
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
  return PMMG_SUCCESS;

end_adjncy:
  PMMG_DEL_MEM(parmesh, adjncy, adjsize, idx_t, "deallocate adjncy" );
#warning IS MESH->ADJA DEALLOCATION CORRECT? IF I CAN REMOVE IT, THEN MULTIPLE RETURN PATHS ARE NOT REQUIRED
end_mesh_adja:
  PMMG_DEL_MEM(mesh, mesh->adja, 4*mesh->nemax + 5, int, "deallocate mesh->adja"); 
end_xadj:
  PMMG_DEL_MEM(parmesh, xadj, mesh->ne + 1, idx_t, "deallocate xadj" );
  return PMMG_FAILURE;
}

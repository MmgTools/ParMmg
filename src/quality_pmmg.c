#include "parmmg.h"

/**
 * \param parmesh pointer to parmesh structure
 *
 * \return 1 if success, 0 if fail;
 *
 * Print quality histogram among all group meshes and all processors
 */
int PMMG_outqua( PMMG_pParMesh parmesh )
{
  PMMG_pGrp grp;
  int i, j;
  int ne, ne_cur, ne_result;
  double max, max_cur, max_result;
  double avg, avg_cur, avg_result;
  double min, min_cur, min_result;
  int iel, iel_cur;
  int good, good_cur, good_result;
  int med, med_cur, med_result;
  const int HIS_SIZE = 5;
  int his[ HIS_SIZE ], his_cur[ HIS_SIZE ], his_result[ HIS_SIZE ];
  int nrid, nrid_cur, nrid_result;

  // Calculate the quality values for local process
  ne = 0;
  max = FLT_MIN;
  avg = 0.;
  min = FLT_MAX;
  iel = 0;
  good = 0;
  med = 0;
  for ( i = 0; i < HIS_SIZE; ++i )
    his[ i ] = 0;
  nrid = 0;
  for ( i = 0; i < parmesh->ngrp; ++i ) {
    grp  = &parmesh->listgrp[ i ];
    MMG3D_computeOutqua( grp->mesh, grp->met, &ne_cur, &max_cur, &avg_cur, &min_cur,
                         &iel_cur, &good_cur, &med_cur, his_cur, &nrid_cur );
    ne += ne_cur;
    avg += avg_cur;
    med += med_cur;
    good += good_cur;
    if ( max_cur > max )
      max = max_cur;
    if ( min_cur < min ) {
      min = min_cur;
      iel = iel_cur;
    }
    for ( j = 0; j < HIS_SIZE; ++j )
      his[ j ] += his_cur[ j ];
    nrid += nrid_cur;
  }

  // Calculate the quality values for all processes
  MPI_Reduce( &ne, &ne_result, 1, MPI_INT, MPI_SUM, 0, parmesh->comm );
  MPI_Reduce( &avg, &avg_result, 1, MPI_DOUBLE, MPI_SUM, 0, parmesh->comm );
  MPI_Reduce( &med, &med_result, 1, MPI_INT, MPI_SUM, 0, parmesh->comm );
  MPI_Reduce( &good, &good_result, 1, MPI_INT, MPI_SUM, 0, parmesh->comm );
  MPI_Reduce( &max, &max_result, 1, MPI_DOUBLE, MPI_MAX, 0, parmesh->comm );
  MPI_Reduce( &min, &min_result, 1, MPI_DOUBLE, MPI_MIN, 0, parmesh->comm );
  MPI_Reduce( his, his_result, HIS_SIZE, MPI_INT, MPI_SUM, 0, parmesh->comm );
  MPI_Reduce( &nrid, &nrid_result, 1, MPI_INT, MPI_SUM, 0, parmesh->comm );
  if ( parmesh->myrank == 0 ) {
    grp = &parmesh->listgrp[0];
    if ( MMG3D_displayQualHisto( ne_result, max_result, avg_result, min_result,
                                 iel_cur, good_result, med_result, his_result,
                                 nrid_result,grp->mesh->info.optimLES,
                                 grp->mesh->info.imprim) )
      return 0;
  }
  return 1;
}


#include "parmmg.h"

/**
 * \param parmesh pointer toward a parmesh structure
 * \param 
 * \param
 * \param
 *
 * \return 1 if success, 0 if fail;
 *
 * 
 *
 */
//NIKOS DEL ME: int PMMG_print_lengths_histogram( PMMG_pParMesh parmesh )
int PMMG_prilen( PMMG_pParMesh parmesh )
{
  //PMMG_pGrp grp = parmesh->listgrp;
  //MMG5_pMesh mesh = grp->mesh;

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param rapmin the minimal (normalized) element quality (over the multiple parmeshes)
 * \param rapmax the maximal one
 * \param rapavg the average one
 * \param his    the number of element per quality span (0-0.2,0.2-0.4...0.8-1)
 *
 * \return 1 if success, 0 if fail;
 *
 * 
 *
 */
int PMMG_computeOutqua( PMMG_pParMesh parmesh, int* rapmin, int* rapmax, int* rapavg, int his[5] )
{
}

//NIKOS DEL ME int PMMG_print_quality_histogram( PMMG_pParMesh parmesh )
int PMMG_outqua( PMMG_pParMesh parmesh )
{
  PMMG_pGrp grp;
  int i, j;
  int ne, ne_cur;
  double max, max_cur;
  double avg, avg_cur;
  double min, min_cur;
  int iel, iel_cur;
  int good, good_cur;
  int med, med_cur;
  int his[ 5 ], his_cur[ 5 ];
  int nrid, nrid_cur;

  // Calculate the quality values for local process
  ne = 0;
  max = FLT_MIN;
  avg = 0.;
  min = FLT_MAX;
  iel = 0;
  good = 0;
  med = 0;
  for ( i = 0; i < 5; ++i )
    his[ i ] = 0;
  nrid = 0;
  for ( i = 0; i < parmesh->ngrp; ++i ) {
    grp  = &parmesh->listgrp[ i ];
 //   PMMG_computeOutqua( parmesh, mesh, &rapmin_cur, &rapmax_cur, &rapavg_cur, his_cur );
    MMG3D_computeOutqua( grp->mesh, grp->met, &ne_cur, &max_cur, &avg_cur, &min_cur,
                         &iel_cur, &good_cur, &med_cur, his_cur, &nrid_cur );
 //   printf(" %d - grp %d: ne %d, max %f, avg %f, min %f, iel %8d, good %8d, med %8d,"
   //        "his %5d-%5d-%5d-%5d-%5d, nrid %d\n",
     //      parmesh->myrank, i, ne_cur, max_cur, avg_cur, min_cur, iel_cur, good_cur, 
       //    med_cur, his_cur[0], his_cur[1], his_cur[2], his_cur[3], his_cur[4], nrid_cur );
    ne += ne_cur;
    avg += avg_cur;
    med += med_cur;
    good += good_cur;
    if ( max_cur > max )
      max = max_cur;
    if ( min_cur < min )
      min = min_cur;
    for ( j = 0; j < 5; ++j )
      his[ j ] += his_cur[ j ];
    nrid += nrid_cur;
  }

  //MPI_Bcast( MPI_MIN, rapmin );
  //MPI_Bcast( MPI_MAX, rapmax );
if ( parmesh->myrank == 1 ){
  printf(" %d : max %f, min %f, avg: %f, good(> 0.12): %f \n", parmesh->myrank, max, min, avg/ne, 100.*((float)good/(float)ne) );
   Αυτό εδώ δεν μπορεί να μείνει ως έχει γιατί τα references στο mesh/met δεν είναι σωστά. ειτε μπορώ να δοκιμάσω σε ποιο mesh/met αναφέρεται το iel (που είναι προβληματικό με απλο MPI_Bcast(MPI_MAX/MIN/klp) είτε να αλλάξει η συνάρτηση ε
ίτε να κάνω copy paste σε κάτι που να μην τα έχει μέσα.
  MMG3D_displayQualHisto( parmesh->listgrp[ 0 ].mesh, parmesh->listgrp[ 0 ].met, ne, max, avg, min, iel, good, med, his, nrid );
}
  return 1;
}



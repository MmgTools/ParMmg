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
int PMMG_mark_boulevolp( PMMG_pParMesh parmesh, MMG5_pMesh mesh, int start, int ip, int * list){
  MMG5_pTetra  pt,pt1;
  int    *adja,nump,ilist,base,cur,k,k1;
  int     nprocs,ngrp,colorp;
  char    j,l,i;

  base = ++mesh->base;
  pt   = &mesh->tetra[start];
  nump = pt->v[ip];

  /* Get point color */
  nprocs = parmesh->nprocs;
  ngrp   = parmesh->ngrp;
  colorp = mesh->point[ip].tmp / (nprocs*ngrp+1) - 1;

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
      if ( pt1->mark <= colorp ) continue;
      pt1->flag = base;
      pt1->mark = colorp;
      for (j=0; j<4; j++)
        if ( pt1->v[j] == nump )  break;
      assert(j<4);
      /* overflow */
      if ( ilist > MMG3D_LMAX-3 )  return 0;
      list[ilist] = 4*k1+j;
      ilist++;
    }
    cur++;
  }

  mesh->point[ip].tmp *= -1;
  return ilist;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param mesh pointer toward the mesh structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * Mark interface points from the tetra mark field.
 *
 * */
int PMMG_mark_interfacePoints( PMMG_pParMesh parmesh, MMG5_pMesh mesh ) {
  MMG5_pTetra pt;
  MMG5_pPoint ppt;
  int         nprocs,ngrp,shift;
  int         colort,colorp,ie,iloc;

  nprocs = parmesh->nprocs;
  ngrp   = parmesh->ngrp;
  shift  = nprocs*ngrp+1;

  /* Mark interface points with the maximum color */
  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    if( !MG_EOK(pt) ) continue;
    colort = pt->mark+1;

    for( iloc = 0; iloc < 4; iloc++ ) {
      ppt = &mesh->point[pt->v[iloc]];
      colorp = abs( ppt->tmp/shift );

      if( !colorp ) {
        ppt->tmp = -( colort+shift*ie );
        continue;
      }
      else if( colorp < colort )
        ppt->tmp = colort+shift*ie;
    }
  }

  return 1;
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

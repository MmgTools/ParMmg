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
 * \file coorcell_pmmg.c
 * \brief functions to manage cells of coordinates and index
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#include "coorcell_pmmg.h"

/**
 * \param a  pointer toward a PMMG_coorCell structure.
 * \param b  pointer toward a PMMG_coorCell structure.
 *
 * \return 1 if a is different from b, 0 if they are
 * equals to within a given tolerance.
 *
 * Compare 2 coor cells (can be used inside the qsort C fnuction), first on
 * their x-coordinates, second ond their y-coordinates then on their
 * z-coordinates.
 *
 */
int PMMG_compare_coorCell (const void * a, const void * b) {
  PMMG_coorCell *cell1,*cell2;
  double        dist[3],tol;
  int           k;

  cell1 = (PMMG_coorCell*)a;
  cell2 = (PMMG_coorCell*)b;

  /* Before this step, each group has been scaled and unscaled 3 times. Each
   * scaling operator can generate an error of the order of 3x10e-15, thus,
   * the error of approximation of the point coordinates is of order
   * 18*10e-15 < 20*MMG5_EPSOK.
   * In a rerun, this would be done twice, so take 50*MMG5_EPSOK. */
  tol  = 50.0;

  for ( k=0; k<3; ++k ) {
    dist[k] = cell1->c[k]-cell2->c[k];

    if ( fabs(dist[k]) >  MMG5_EPSOK*tol ) return 1;

  }

  assert ( dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]<MMG5_EPSD*tol*tol
    && "Wrong epsilon machine consistency: function to improve");

  return 0;
}

/**
 * \param list   array of PMMG_coorCell
 * \param nitem  number of items in the list
 * \param min    minimal coordinates in each direction
 * \param max    maximal coordinates in each direction
 * \param delta  scaling value
 *
 * \return 1 if success, 0 if fail;
 *
 * Find the bounding box and fill the scaling data (\a min \a max \a delta) for
 * the points listed in the \a list array.
 *
 */
int PMMG_find_coorCellListBoundingBox (PMMG_coorCell *list,int nitem,
                                       double min[3],double max[3],
                                       double *delta) {
  int    i,j;

  /* Bounding box computation */
  for (i=0; i<3; i++) {
    min[i] =  DBL_MAX;
    max[i] = -DBL_MAX;
  }
  for ( i=0; i<nitem; ++i ) {
    for (j=0; j<3; j++) {
      if ( list[i].c[j] > max[j] ) max[j] = list[i].c[j];
      if ( list[i].c[j] < min[j] ) min[j] = list[i].c[j];
    }
  }
  (*delta) = 0.0;
  for (i=0; i<3; i++)
    if ( max[i]-min[i] > (*delta) )  (*delta) = max[i]-min[i];

  return 1;
}

/**
 * \param list   array of PMMG_coorCell
 * \param nitem  number of items in the list
 * \param min    minimal coordinates in each direction
 * \param max    maximal coordinates in each direction
 * \param delta  scaling value
 *
 * \return 1 if success, 0 if fail;
 *
 * Scale the coordinates listed in the \a list array and fill the scaling data
 * (\a min \a max \a delta).
 *
 */
int PMMG_scale_coorCellList (PMMG_coorCell *list,int nitem,double min[3],
                             double max[3],double *delta) {
  double dd;
  int    i,j;

  /* Bounding box computation */
  if ( !PMMG_find_coorCellListBoundingBox(list,nitem,min,max,delta) )
    return 0;

  /* Coordinate scaling */
  dd = 1./(*delta);
  for ( i=0; i<nitem; ++i ) {
    for (j=0; j<3; j++) list[i].c[j] = dd*(list[i].c[j]-min[j]);
  }

  return 1;
}

/**
 * \param list   array of PMMG_coorCell
 * \param nitem  number of items in the list
 * \param min    minimal coordinates in each direction
 * \param max    maximal coordinates in each direction
 * \param delta  scaling value
 *
 * \return 1 if success, 0 if fail;
 *
 * Unscale the coordinates listed in the \a list array.
 *
 */
int PMMG_unscale_coorCellList (PMMG_coorCell *list,int nitem,double min[3],
                               double max[3],double delta) {
  int    i,j;

  /* Coordinate scaling */
  for ( i=0; i<nitem; ++i )
    for (j=0; j<3; j++) list[i].c[j] = delta*list[i].c[j]+min[j];

  return 1;
}

/**
 * \file coorcell_pmmg.h
 * \brief coorcell_pmmg.c header file
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#ifndef COORCELL_PMMG_H

#define COORCELL_PMMG_H

#include "parmmg.h"

/**
 * \struct PMMG_coorCell
 *
 * \brief Cell containing a point coordinate and its position in an array
 *
 */
typedef struct {
  double  coor[3]; /*!< point coordinates */
  int     pos;     /*!< position of the coordinate in array */
} PMMG_coorCell;


int PMMG_compare_coorCell (const void * a, const void * b);
int PMMG_scale_coorCellList (PMMG_coorCell*,int,double*,double*,double*);
int PMMG_unscale_coorCellList (PMMG_coorCell*,int,double*,double*,double);

#endif

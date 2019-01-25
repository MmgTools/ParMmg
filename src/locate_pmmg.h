/**
 * \file locate_pmmg.h
 * \brief Point localization for interpolation on a new mesh.
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \author Luca Cirrottola (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#ifndef LOCATE_PMMG_H

#define LOCATE_PMMG_H

#include "parmmg.h"

/** \struct PMMG_baryCoord
 *
 * \brief Struct containing the index and value of a barycentric coordinate
 *
 */
typedef struct {
  int    idx; /*!< direction */
  double val; /*!< coordinate value */
} PMMG_baryCoord;

#endif

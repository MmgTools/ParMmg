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
 * \file metis_pmmg.h
 * \brief metis_pmmg.c header file
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */

#ifndef METIS_PMMG_H

#define METIS_PMMG_H

#include "parmmg.h"
#include <metis.h>

#ifdef USE_PARMETIS
#include <parmetis.h>
#endif

/* Available choices for the wgtflag parameter of ParMetis */
#define PMMG_WGTFLAG_NONE  0
#define PMMG_WGTFLAG_ADJ   1
#define PMMG_WGTFLAG_VTX   2
#define PMMG_WGTFLAG_BOTH  3

/**
 * \def PMMG_WGTFLAG_DEF
 *
 * value for the wgtflag parameter of ParMetis
 *
 */
#define PMMG_WGTFLAG_DEF   (PMMG_WGTFLAG_ADJ)

/**
 * \def PMMG_WGTVAL_HUGEINT
 *
 * Huge integer weight for parallel faces in load balancing
 *
 */
#define PMMG_WGTVAL_HUGEINT   1000000

/**
 * \def PMMG_UBVEC_DEF
 *
 * value for the ubvec vector of ParMetis
 *
 */
#define PMMG_UBVEC_DEF     1.05

/**
 * \def PMMG_CONTIG_DEF
 *
 * value for option[METIS_OPTION_CONTIG] in Metis
 *
 */
#define PMMG_CONTIG_DEF     0


/**
 * \def PMMG_NBADJA_GRPS
 *
 * Supposed maximum number of groups adjacent to one given group
 *
 */
#define PMMG_NBADJA_GRPS   30

/**
 * \struct PMMG_hgrp
 *
 * \brief Cell of the hastable HGrp
 *
 */
typedef struct {
  idx_t   adj; /*!< item to add */
  idx_t   wgt; /*!< Edge weight */
  int     nxt; /*!< Next element of the linked list */
} PMMG_hgrp;

/**
 * \struct PMMG_HGrp
 *
 * \brief Hash table to store the list of adjacents groups to a given group.
 *
 */
typedef struct {
  int          siz,max,nxt;
  PMMG_hgrp    *item;
} PMMG_HGrp;

int PMMG_checkAndReset_grps_contiguity( PMMG_pParMesh parmesh );
int PMMG_check_grps_contiguity( PMMG_pParMesh parmesh );
int PMMG_graph_meshElts2metis(PMMG_pParMesh,MMG5_pMesh,MMG5_pSol,idx_t**,idx_t**,idx_t**,idx_t*);
int PMMG_part_meshElts2metis( PMMG_pParMesh,idx_t*,idx_t);
int PMMG_graph_parmeshGrps2parmetis(PMMG_pParMesh,idx_t**,idx_t**,idx_t**,idx_t*,
                                    idx_t**,idx_t**,idx_t*,idx_t*,idx_t*,idx_t,
                                    real_t**,real_t**);
int PMMG_part_parmeshGrps2parmetis(PMMG_pParMesh,idx_t*,idx_t);
int PMMG_part_parmeshGrps2metis(PMMG_pParMesh,idx_t*,idx_t);
int PMMG_correct_parmeshGrps2parmetis( PMMG_pParMesh parmesh, idx_t *vtxdist,
                                       idx_t* mypart,idx_t nproc );

#endif

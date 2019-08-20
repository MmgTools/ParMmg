/**
 * \file moveinterfaces_pmmg.h
 * \brief moveinterfaces_pmmg.c header file
 * \author Luca Cirrottola (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#ifndef MOVEINTERFACES_PMMG_H
#define MOVEINTERFACES_PMMG_H

#include "parmmg.h"

#ifdef __cplusplus
extern "C" {
#endif

void PMMG_set_color_tetra( PMMG_pParMesh parmesh,int igrp );

#ifdef __cplusplus
}
#endif

#endif

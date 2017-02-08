/**
 * \file libparmmg.h
 * \brief API headers for the parmmg library
 * \author
 * \version
 * \date 11 2016
 * \copyright
 */

#ifndef _PMMGLIB_H
#define _PMMGLIB_H

#ifdef __cplusplus
extern "C" {
#endif

#include "metis.h"

#include "libparmmgtypes.h"
#include "mmg3d.h"

#define PMMG_VER   "1.0.0"
#define PMMG_REL   "2016"
#define PMMG_CPY   "Copyright (c) Bx INP/INRIA, 2016-"
#define PMMG_STR   "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"

/* API_functions_pmmg.c */
int PMMG_Init_parMesh(const int starter,...);

/* inout_3d.c */
int PMMG_loadMesh(PMMG_pParMesh ,const char *);
int PMMG_saveMesh(PMMG_pParMesh ,const char *);
int PMMG_loadSol(PMMG_pParMesh ,const char *);
int PMMG_saveSol(PMMG_pParMesh ,const char *);

/* libparmmg.c */
/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return \ref PMMG_SUCCESS if success, \ref PMMG_LOWFAILURE if fail but a
 * conform mesh is saved or \ref PMMG_STRONGFAILURE if fail and we can't save
 * the mesh.
 *
 * Main program for the parallel remesh library.
 *
 **/
int PMMG_parmmglib(PMMG_pParMesh parmesh);

/* metisfunctions.c */
int PMMG_metispartitioning(PMMG_pParMesh ,idx_t *);

/*distributemesh*/
int PMMG_distributeMesh(PMMG_pParMesh ,idx_t *);

/*mergeMesh*/
int PMMG_mergeMesh(PMMG_pParMesh);

#ifdef __cplusplus
}
#endif

#endif

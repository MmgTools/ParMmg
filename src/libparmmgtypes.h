#include "mmg/libmmgtypes.h"
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
/**
 * \def PMMG_SUCCESS
 *
 * Return value for success.
 *
 */
#define PMMG_SUCCESS       0
/**
 * \def PMMG_LOWFAILURE
 *
 * Return value if the remesh process failed but we can save a conform
 * mesh.
 *
 */
#define PMMG_LOWFAILURE    1
/**
 * \def PMMG_STRONGFAILURE
 *
 * Return value if the remesh process failed and the mesh is
 * non-conform.
 *
 */
#define PMMG_STRONGFAILURE 2


/**
 * Types
 */
/**
 * \struct PARMMG_int_comm
 * \brief internal communicator structure.
 */
typedef struct {
  
} PMMG_int_comm;
typedef PMMG_int_comm  * PMMG_pint_comm;

/**
 * \struct PARMMG_ext_comm
 * \brief external communicator structure.
 */
typedef struct {
  
} PMMG_ext_comm;
typedef PMMG_ext_comm  * PMMG_pext_comm;

/**
 * \struct PARMMG_Grp
 * \brief Grp mesh structure.
 */
typedef struct {
  MMG5_pMesh   mesh;  /*mesh definition : coordinates, tetra etc..*/
  MMG5_pSol    sol;   /*physical solutions defined on each point of the mesh*/
  MMG5_pSol    met;   /*metric*/

  
} PMMG_Grp;
typedef PMMG_Grp  * PMMG_pGrp;


/**
 * \struct PMMG_ParMesh
 * \brief ParMmg mesh structure.
 */
typedef struct {
  
  /* mpi info*/
  int      nprocs;
  int      myrank;
  
  /* grp */
  int       neltpergrp; /*!< Average number of elements per group */
  int       ngrp;       /*!< Number of grp */
  PMMG_pGrp listgrp;    /*!< List of grp */

  /* internal communicators */
  PMMG_pint_comm  int_face_comm; /*!< Internal face communicator */

  /* external communicators */
  PMMG_pext_comm send_ext_face_comm; 

} PMMG_ParMesh;
typedef PMMG_ParMesh  * PMMG_pParMesh;



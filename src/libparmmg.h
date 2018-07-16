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

#include "libparmmgtypes.h"
#include "mmg3d.h"
#include "metis.h"


#define PMMG_VER   "1.0.0"
#define PMMG_REL   "2016"
#define PMMG_CPY   "Copyright (c) Bx INP/INRIA, 2016-"
#define PMMG_STR   "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"

/**
 * \enum PMMG_Param
 * \brief Input parameters for mmg library.
 *
 * Input parameters for mmg library. Options prefixed by \a
 * PMMG_IPARAM asked for integers values ans options prefixed by \a
 * PMMG_DPARAM asked for real values.
 *
 */
enum PMMG_Param {
  PMMG_IPARAM_verbose,           /*!< [-10..10], Tune level of verbosity */
  PMMG_IPARAM_mem,               /*!< [n/-1], Set memory size to n Mbytes or keep the default value */
  PMMG_IPARAM_debug,             /*!< [1/0], Turn on/off debug mode */
  PMMG_IPARAM_angle,             /*!< [1/0], Turn on/off angle detection */
  PMMG_IPARAM_iso,               /*!< [1/0], Level-set meshing */
  PMMG_IPARAM_lag,               /*!< [-1/0/1/2], Lagrangian option */
  PMMG_IPARAM_optim,             /*!< [1/0], Optimize mesh keeping its initial edge sizes */
  PMMG_IPARAM_optimLES,          /*!< [1/0], Strong mesh optimization for Les computations */
  PMMG_IPARAM_noinsert,          /*!< [1/0], Avoid/allow point insertion */
  PMMG_IPARAM_noswap,            /*!< [1/0], Avoid/allow edge or face flipping */
  PMMG_IPARAM_nomove,            /*!< [1/0], Avoid/allow point relocation */
  PMMG_IPARAM_nosurf,            /*!< [1/0], Avoid/allow surface modifications */
  PMMG_IPARAM_numberOfLocalParam,/*!< [n], Number of local parameters */
  PMMG_IPARAM_anisosize,         /*!< [1/0], Turn on/off anisotropic metric creation when no metric is provided */
  PMMG_IPARAM_octree,            /*!< [n], Specify the max number of points per octree cell (DELAUNAY) */
  PMMG_DPARAM_angleDetection,    /*!< [val], Value for angle detection */
  PMMG_DPARAM_hmin,              /*!< [val], Minimal mesh size */
  PMMG_DPARAM_hmax,              /*!< [val], Maximal mesh size */
  PMMG_DPARAM_hausd,             /*!< [val], Control global Hausdorff distance (on all the boundary surfaces of the mesh) */
  PMMG_DPARAM_hgrad,             /*!< [val], Control gradation */
  PMMG_DPARAM_ls,                /*!< [val], Value of level-set */
  PMMG_PARAM_size,               /*!< [n], Number of parameters */
};


/**< Subgroups target size for a fast remeshing step */
static const int REMESHER_TARGET_MESH_SIZE = 16 * 1024;
/**< Subgroups target size for parmetis */
static const int METIS_TARGET_MESH_SIZE = 8000;


/* API_functions_pmmg.c */
/**
 * \param parmesh pointer toward the parmesh structure
 * \param comm MPI communicator for ParMmg
 *
 * \return \ref PMMG_SUCCESS if success, \ref PMMG_LOWFAILURE if fail but a
 * conform mesh is saved or \ref PMMG_STRONGFAILURE if fail and we can't save
 * the mesh.
 *
 * Allocation of the parmesh structure and initialisation of the parameters
 *
 **/
int PMMG_Init_parMesh( PMMG_pParMesh *parmesh,MPI_Comm comm );

/* libparmmg.c */
/**
 * \param parmesh pointer toward the parmesh structure (boundary entities are
 * stored into MMG5_Tria, MMG5_Edge... structures)
 *
 * \return \ref PMMG_SUCCESS if success, \ref PMMG_LOWFAILURE if fail but a
 * conform mesh is saved or \ref PMMG_STRONGFAILURE if fail and we can't save
 * the mesh.
 *
 * Main program for the parallel remesh library.
 *
 **/
int PMMG_parmmglib(PMMG_pParMesh parmesh);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param comm MPI communicator for ParMmg
 *
 * Initialization of the input parameters (stored in the Info structure).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_INIT_PARAMETERS(parmesh,comm)\n
 * >     PMMG_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     MPI_Comm,INTENT(IN)           :: comm\n
 * >   END SUBROUTINE\n
 *
 */
void  PMMG_Init_parameters(PMMG_pParMesh parmesh,MPI_Comm comm);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param iparam integer parameter to set (see \a MMG3D_Param structure).
 * \param val value for the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set integer parameter \a iparam at value \a val.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_IPARAMETER(parmesh,iparam,val,retval)\n
 * >     PMMG_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: iparam,val\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Set_iparameter(PMMG_pParMesh parmesh, int iparam, int val);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Set integer parameter \a iparam at value \a val.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_DISTRIBUTE_MESH(parmesh,retval)\n
 * >     PMMG_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_distribute_mesh(PMMG_pParMesh);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Set integer parameter \a iparam at value \a val.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_DISTRIBUTE_MESH(parmesh,retval)\n
 * >     PMMG_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_merge_parmesh( PMMG_pParMesh parmesh );

#ifdef __cplusplus
}
#endif

#endif

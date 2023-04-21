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
 * \file libparmmg.h
 * \brief API headers for the parmmg library
 * \author
 * \version
 * \date 11 2016
 * \copyright
 */

#ifndef _PMMGLIB_H
#define _PMMGLIB_H

#include "libparmmgtypes.h"

#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

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
  PMMG_IPARAM_mmgVerbose,        /*!< [-10..10], Tune level of verbosity of Mmg */
  PMMG_IPARAM_mem,               /*!< [n/-1], Set memory size to n Mbytes or keep the default value */
  PMMG_IPARAM_debug,             /*!< [1/0], Turn on/off debug mode */
  PMMG_IPARAM_distributedOutput, /*!< [0/1], Turn off/on distributed output */
  PMMG_IPARAM_mmgDebug,          /*!< [1/0], Turn on/off debug mode */
  PMMG_IPARAM_angle,             /*!< [1/0], Turn on/off angle detection */
  PMMG_IPARAM_iso,               /*!< [1/0], Level-set meshing */
  PMMG_IPARAM_lag,               /*!< [-1/0/1/2], Lagrangian option */
  PMMG_IPARAM_opnbdy,            /*!< [0/1], Enable preservation of open boundaries */
  PMMG_IPARAM_optim,             /*!< [1/0], Optimize mesh keeping its initial edge sizes */
  PMMG_IPARAM_optimLES,          /*!< [1/0], Strong mesh optimization for Les computations */
  PMMG_IPARAM_nofem,             /*!< [1/0], Generate a non finite element mesh */
  PMMG_IPARAM_noinsert,          /*!< [1/0], Avoid/allow point insertion */
  PMMG_IPARAM_noswap,            /*!< [1/0], Avoid/allow edge or face flipping */
  PMMG_IPARAM_nomove,            /*!< [1/0], Avoid/allow point relocation */
  PMMG_IPARAM_nosurf,            /*!< [1/0], Avoid/allow surface modifications */
  PMMG_IPARAM_numberOfLocalParam,/*!< [n], Number of local parameters */
  PMMG_IPARAM_anisosize,         /*!< [1/0], Turn on/off anisotropic metric creation when no metric is provided */
  PMMG_IPARAM_octree,            /*!< [n], Specify the max number of points per octree cell (DELAUNAY) */
  PMMG_IPARAM_meshSize,          /*!< [n], Target mesh size of Mmg (advanced use) */
  PMMG_IPARAM_nobalancing,       /*!< [1/0], Deactivate load balancing of the output mesh */
  PMMG_IPARAM_metisRatio,        /*!< [n], wanted ratio # mesh / # metis super nodes (advanced use) */
  PMMG_IPARAM_ifcLayers,         /*!< [n], Number of layers of interface displacement */
  PMMG_DPARAM_groupsRatio,       /*!< [val], Allowed imbalance between current and desired groups size */
  PMMG_IPARAM_APImode,           /*!< [0/1], Initialize parallel library through interface faces or nodes */
  PMMG_IPARAM_globalNum,         /*!< [1,0], Compute nodes and triangles global numbering in output */
  PMMG_IPARAM_niter,             /*!< [n], Set the number of remeshing iterations */
  PMMG_DPARAM_angleDetection,    /*!< [val], Value for angle detection */
  PMMG_DPARAM_hmin,              /*!< [val], Minimal mesh size */
  PMMG_DPARAM_hmax,              /*!< [val], Maximal mesh size */
  PMMG_DPARAM_hsiz,              /*!< [val], Constant mesh size */
  PMMG_DPARAM_hausd,             /*!< [val], Control global Hausdorff distance (on all the boundary surfaces of the mesh) */
  PMMG_DPARAM_hgrad,             /*!< [val], Control gradation */
  PMMG_DPARAM_hgradreq,          /*!< [val], Control gradation from required entities */
  PMMG_DPARAM_ls,                /*!< [val], Value of level-set */
  PMMG_PARAM_size,               /*!< [n], Number of parameters */
};


/* API_functions_pmmg.c */
/* init structures */
/**
 * \param starter dummy argument used to initialize the variadic argument list
 * \param ... variadic arguments that depend to the parmesh fields that you want
 * to init
 *
 * You need to provide at least the following arguments:
 * the \a PMMG_ARG_start keyword to start the list of variadic arguments
 * the \a PMMG_ARG_ppParMesh keyword to say that the next argument is a pointer
 *  toward a pointer toward a parmesh
 * a pointer toward a pointer toward a parmesh
 * the \a PMMG_ARG_pMesh keyword to initialize a \a mesh pointer inside your \a parmesh
 * the \a PMMG_ARG_pMet keyword to initialize a \a metric pointer inside your \a parmesh
 * the \a PMMG_ARG_dim keyword to set the mesh dimension
 * the \a PMMG_ARG_MPIComm keyword to set the MPI Communicator in which parmmg will work
 * the \a PMMG_ARG_end keyword to end the list of variadic args.
 *
 * Example :
 * PMMG_Init_parmesh(PMMG_ARG_start,PMMG_ARG_ppParMesh,your_parmesh_address,
 *   PMMG_ARG_pMesh,PMMG_ARG_pMet,PMMG_ARG_dim,mesh_dimension,
 *   PMMG_ARG_MPIComm,mpi_communicator,PMMG_ARG_end)
 *
 * \return 1 if success, 0 if fail
 *
 * ParMMG structures allocation and initialization.
 *
 * \remark No fortran interface to allow variadic arguments.
 *
 */
 int PMMG_Init_parMesh(const int starter,...);

/* libparmmg.c */
/**
 * \param parmesh pointer toward the parmesh structure (boundary entities are
 * stored into MMG5_Tria, MMG5_Edge... structures)
 *
 * \return \ref PMMG_SUCCESS if success, \ref PMMG_LOWFAILURE if fail but we can
 * one unscaled mesh per proc or \ref PMMG_STRONGFAILURE if fail and
 * we can't return one unscaled mesh per proc.
 *
 * Main program for the parallel remesh library for distributed meshes.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_parmmglib_distributed(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 **/
int PMMG_parmmglib_distributed(PMMG_pParMesh parmesh);

/**
 * \param parmesh pointer toward the parmesh structure (boundary entities are
 * stored into MMG5_Tria, MMG5_Edge... structures)
 *
 * \return \ref PMMG_SUCCESS if success, \ref PMMG_LOWFAILURE if fail but we can
 * return a centralized and unscaled mesh or \ref PMMG_STRONGFAILURE if fail and
 * we can't return a centralized and unscaled mesh.
 *
 * Main program for the parallel remesh library for centralized meshes
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_parmmglib_centralized(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 **/
int PMMG_parmmglib_centralized(PMMG_pParMesh parmesh);

/**
 * \param parmesh pointer toward the parmesh structure (boundary entities are
 * stored into MMG5_Tria, MMG5_Edge... structures)
 *
 * \return \ref PMMG_SUCCESS if success, \ref PMMG_LOWFAILURE if fail but we can
 * return a centralized and unscaled mesh or \ref PMMG_STRONGFAILURE if fail and
 * we can't return a centralized and unscaled mesh.
 *
 * Main program for the parallel isovalue discretisation library for centralized
 * meshes
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_parmmgls_centralized(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 **/
int PMMG_parmmgls_centralized(PMMG_pParMesh parmesh);

/**
 * \param parmesh pointer toward the parmesh structure (boundary entities are
 * stored into MMG5_Tria, MMG5_Edge... structures)
 *
 * \return \ref PMMG_SUCCESS if success, \ref PMMG_LOWFAILURE if fail but we can
 * return a centralized and unscaled mesh or \ref PMMG_STRONGFAILURE if fail and
 * we can't return a centralized and unscaled mesh.
 *
 * Main program for the parallel isovalue discretisation library and remesh
 * library for centralized meshes
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_parmmg_centralized(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 **/
int PMMG_parmmg_centralized(PMMG_pParMesh parmesh);

/* init file names */
/**
 * \param parmesh pointer toward a parmesh structure.
 * \param meshin input mesh name.
 * \return 1.
 *
 * Set the name of input mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_INPUTMESHNAME(parmesh,meshin,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: meshin\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Set_inputMeshName(PMMG_pParMesh parmesh,const char* meshin);
/**
 * \param parmesh pointer toward a parmesh structure.
 * \param meshout name of the output mesh file.
 * \return 1.
 *
 * Set the name of output mesh file.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_OUTPUTMESHNAME(parmesh,meshout,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: meshout\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Set_outputMeshName(PMMG_pParMesh parmesh, const char* meshout);
/**
 * \param parmesh pointer toward a parmesh structure.
 * \param solin name of the input solution file.
 * \return 1.
 *
 * Set the name of input solution file. The call to this function is mandatory
 * to load a solution field.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_INPUTSOLSNAME(parmesh,solin,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: solin\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Set_inputSolsName(PMMG_pParMesh parmesh,const char* solin);
/**
 * \param parmesh pointer toward a parmesh structure.
 * \param metin name of the input metric file.
 * \return 1.
 *
 * Set the name of input metric file. The parmesh must have been initialized
 * with a metric field (otherwise the metric structure is not allocated).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_INPUTMETNAME(parmesh,metin,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: metin\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Set_inputMetName(PMMG_pParMesh parmesh,const char* metin);
/**
 * \param parmesh pointer toward a parmesh structure.
 * \param lsin name of the input level-set file.
 * \return 1.
 *
 * Set the name of input level-set file. The parmesh must have been initialized
 * with a level-set field (otherwise the ls structure is not allocated).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_INPUTLSNAME(parmesh,lsin,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: lsin\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Set_inputLsName(PMMG_pParMesh parmesh,const char* lsin);
/**
 * \param parmesh pointer toward a parmesh structure.
 * \param dispin name of the input displacement file.
 * \return 1.
 *
 * Set the name of input displacement file. The parmesh must have been
 * initialized with a displacement field (otherwise the displacement structure
 * is not allocated).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_INPUTDISPNAME(parmesh,dispin,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: dispin\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Set_inputDispName(PMMG_pParMesh parmesh,const char* dispin);
/**
 * \param parmesh pointer toward a parmesh structure.
 * \param solout name of the output solution file.
 * \return 0 if failed, 1 otherwise.
 *
 *  Set the name of output solution file. If not called, an automatic output
 * name is computed from the path of the output mesh and the input solutions
 * name: <inputfield>.sol -> <outputmeshpath>/<inputfield>.o.sol
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_OUTPUTSOLSNAME(parmesh,solout,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: solout\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Set_outputSolsName(PMMG_pParMesh parmesh, const char* solout);
/**
 * \param parmesh pointer toward a parmesh structure.
 * \param metout name of the output metric file.
 * \return 0 if failed, 1 otherwise.
 *
 *  Set the name of output metric file.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_OUTPUTMETNAME(parmesh,metout,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: metout\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Set_outputMetName(PMMG_pParMesh parmesh, const char* metout);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param comm MPI communicator for ParMmg
 *
 * Initialization of the input parameters (stored in the Info structure).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_INIT_PARAMETERS(parmesh,comm)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER,INTENT(IN)            :: comm\n
 * >   END SUBROUTINE\n
 *
 */
void  PMMG_Init_parameters(PMMG_pParMesh parmesh,MPI_Comm comm);

/**
 * \param parmesh   Pointer towards the parmesh structure.
 * \param np        Number of vertices.
 * \param ne        Number of tetrahedra.
 * \param nprism    Number of prisms.
 * \param nt        Number of triangles.
 * \param nquad     Number of quads.
 * \param na        Number of edges.
 * \return          0 if failed, 1 otherwise.
 *
 * Set the number of vertices, tetrahedra, prisms, triangles, quadrilaterals and
 * edges of the mesh and allocate the associated tables.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_MESHSIZE(parmesh, np, ne, nprism, nt, nquad, na, &\n
 * >                                retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: np, ne, nprism, nt, nquad, na\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_meshSize(PMMG_pParMesh parmesh, int np, int ne, int nprism, int nt,
                      int nquad, int na);
/**
 * \param parmesh   Pointer towards the parmesh structure.
 * \param nsols     Number of solutions per entity.
 * \param nentities Number of entities.
 * \param typSol    Array of size nsol listing the type of the solutions
 *                  (scalar, vectorial...).
 * \return          0 if failed, 1 otherwise.
 *
 * Set the number of solutions per entity, the type of entity to which applies
 * each solution and the type of each solution.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_SOLSATVERTICESSIZE(parmesh, nsols,nentities,typSol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: nsols,nentities\n
 * >     INTEGER, INTENT(IN)           :: typSol(*)\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_solsAtVerticesSize(PMMG_pParMesh parmesh,int nsols,
                                   int nentities,int *typSol);

/**
 * \param parmesh   Pointer towards the parmesh structure.
 * \param typEntity Type of entity to which applies the metric (vertices, triangles...).
 * \param np        number of vertices
 * \param typMet    Type of the metric (scalar, vectorial...).
 * \return          0 if failed, 1 otherwise.
 *
 * Set the type of entity to which applies the metric and the type of the metric.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_METSIZE(parmesh,typEntity,np,typMet,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: np\n
 * >     INTEGER, INTENT(IN)           :: typEntity, typMet\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_metSize(PMMG_pParMesh parmesh,int typEntity,int np,int typMet);


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
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: iparam,val\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Set_iparameter(PMMG_pParMesh parmesh, int iparam, int val);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param dparam double parameter to set (see \a MMG3D_Param structure).
 * \param val value for the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set double parameter \a dparam at value \a val.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_DPARAMETER(parmesh,dparam,val,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: dparam\n
 * >     REAL(KIND=8), INTENT(IN)      :: val\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Set_dparameter(PMMG_pParMesh parmesh, int iparam, double val);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param typ type of entity (triangle, edge,...).
 * \param ref reference of the entity.
 * \param hmin minimal edge size.
 * \param hmax maximal edge size.
 * \param hausd value of the Hausdorff number.
 * \return 0 if failed, 1 otherwise.
 *
 * Set local parameters: set the hausdorff value at \a val for all
 * elements of type \a typ and reference \a ref.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_LOCALPARAMETER(parmesh,typ,ref,& \n
 * >                                      hmin,hmax,hausd,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: typ,ref\n
 * >     REAL(KIND=8), INTENT(IN)      :: hmin,hmax,hausd\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Set_localParameter(PMMG_pParMesh parmesh, int typ,
                             int ref,double hmin,double hmax,double hausd);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Free names stored in the parmesh
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_FREE_NAMES(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Free_names(PMMG_pParMesh parmesh);

/**
 * \param starter dummy argument used to initialize the variadic argument list.
 * \param ... variadic arguments to list the structure that must be deallocated.
 * For now, you must provide the PMMG_ARG_ppParMesh keyword and your parmesh
 * address as arguments.
 *
 * \return 1 if success, 0 if fail
 *
 * Deallocations before return.
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 * \remark no Fortran interface to allow variadic args.
 *
 */
  int PMMG_Free_all(const int starter,...);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Set integer parameter \a iparam at value \a val.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_DISTRIBUTEMESH_CENTRALIZED(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_distributeMesh_centralized(PMMG_pParMesh);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Set integer parameter \a iparam at value \a val.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_DISTRIBUTE_MESH(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
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
 * >   SUBROUTINE PMMG_MERGE_PARMESH(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_merge_parmesh( PMMG_pParMesh parmesh );


/**
 * \param parmesh  pointer toward the group structure.
 * \param c0   coordinate of the point along the first dimension.
 * \param c1   coordinate of the point along the second dimension.
 * \param c2   coordinate of the point along the third dimension.
 * \param ref  point reference.
 * \param pos  position of the point in the mesh.
 * \return 1.
 *
 * Set vertex of coordinates \a c0, \a c1,\a c2 and reference \a ref
 * at position \a pos in mesh structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_VERTEX(parmesh,c0,c1,c2,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     REAL(KIND=8), INTENT(IN)      :: c0,c1,c2\n
 * >     INTEGER, INTENT(IN)           :: ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_vertex(PMMG_pParMesh parmesh, double c0, double c1, double c2,
                    int ref, int pos);

/**
 * \param parmesh      pointer toward the group structure.
 * \param vertices table of the points coor.
 * The coordinates of the \f$i^{th}\f$ point are stored in
 * vertices[(i-1)*3]\@3.
 * \param refs     table of points references.
 * The ref of the \f$i^th\f$ point is stored in refs[i-1].
 * \return 1.
 *
 * Set vertices coordinates and references in mesh structure
 * (wrapper for MMG3D function).
 *
 * \remark Fortran interface: (commentated in order to allow to pass
 * \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE PMMG_SET_VERTICES(parmesh,vertices,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * > !    REAL(KIND=8), INTENT(IN)      :: vertices(*)\n
 * > !    INTEGER,INTENT(IN)            :: refs(*)\n
 * > !    INTEGER, INTENT(OUT)          :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int PMMG_Set_vertices(PMMG_pParMesh parmesh, double *vertices, int *refs);

/**
 * \param parmesh pointer toward the group structure.
 * \param v0  first vertex of tetrahedron.
 * \param v1  second vertex of tetrahedron.
 * \param v2  third vertex of tetrahedron.
 * \param v3  fourth vertex of tetrahedron.
 * \param ref tetrahedron reference.
 * \param pos tetrahedron position in the mesh.
 * \return 0  if failed, 1 otherwise.
 *
 * Set tetrahedra of vertices \a v0, \a v1,\a v2,\a v3 and reference
 * \a ref at position \a pos in mesh structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_TETRAHEDRON(parmesh,v0,v1,v2,v3,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: v0,v1,v2,v3,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_tetrahedron(PMMG_pParMesh parmesh, int v0, int v1, int v2, int v3,
                         int ref, int pos);

/**
 * \param parmesh   pointer toward the group structure.
 * \param tetra vertices of the tetras of the mesh
 * Vertices of the \f$i^{th}\f$ tetra are stored in tetra[(i-1)*4]\@4.
 * \param refs  table of the tetrahedra references.
 * References of the \f$i^{th}\f$ tetra is stored in refs[i-1].
 * \return 0    if failed, 1 otherwise.
 *
 * Set vertices and references of the mesh tetrahedra.
 *
 * \remark Fortran interface: (commentated in
 * order to allow to pass \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE PMMG_SET_TETRAHEDRA(parmesh,tetra,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)     :: parmesh\n
 * > !    INTEGER, DIMENSION(*), INTENT(IN) :: tetra,refs\n
 * > !    INTEGER, INTENT(OUT)              :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int PMMG_Set_tetrahedra(PMMG_pParMesh parmesh, int *tetra, int *refs);

/**
 * \param parmesh pointer toward the group structure.
 * \param v0  first vertex of prism.
 * \param v1  second vertex of prism.
 * \param v2  third vertex of prism.
 * \param v3  fourth vertex of prism.
 * \param v4  fifth vertex of prism.
 * \param v5  sixth vertex of prism.
 * \param ref prism reference.
 * \param pos prism position in the mesh.
 * \return 0  if failed, 1 otherwise.
 *
 * Set prisms of vertices \a v0, \a v1,\a v2,\a v3,\a v4,\a v5 and reference
 * \a ref at position \a pos in mesh structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_PRISM(parmesh,v0,v1,v2,v3,v4,v5,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: v0,v1,v2,v3,v4,v5,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_prism(PMMG_pParMesh parmesh, int v0, int v1, int v2,
                   int v3, int v4, int v5, int ref, int pos);

/**
 * \param parmesh    pointer toward the mesh structure.
 * \param prisms vertices of the prisms of the mesh
 * Vertices of the \f$i^{th}\f$ prism are stored in prism[(i-1)*6]\@6.
 * \param refs   table of the prisms references.
 * References of the \f$i^{th}\f$ prisms is stored in refs[i-1].
 * \return 0     if failed, 1 otherwise.
 *
 * Set vertices and references of the mesh prisms (wrapper for MMG3D
 * function).
 *
 * \remark Fortran interface: (commentated in
 * order to allow to pass \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE PMMG_SET_PRISMS(parmesh,prisms,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)     :: parmesh\n
 * > !    INTEGER, DIMENSION(*), INTENT(IN) :: prisms,refs\n
 * > !    INTEGER, INTENT(OUT)              :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int PMMG_Set_prisms(PMMG_pParMesh parmesh, int *prisms, int *refs);

/**
 * \param parmesh pointer toward the group structure.
 * \param v0  first vertex of triangle.
 * \param v1  second vertex of triangle.
 * \param v2  third vertex of triangle.
 * \param ref triangle reference.
 * \param pos triangle position in the mesh.
 * \return 0  if failed, 1 otherwise.
 *
 * Set triangle of vertices \a v0, \a v1, \a v2 and reference \a ref
 * at position \a pos in mesh structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_TRIANGLE(parmesh,v0,v1,v2,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: v0,v1,v2,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_triangle(PMMG_pParMesh parmesh, int v0, int v1, int v2,
                      int ref,int pos);

/**
 * \param parmesh  pointer toward the group structure.
 * \param tria pointer toward the table of the tria vertices
 * Vertices of the \f$i^{th}\f$ tria are stored in tria[(i-1)*3]\@3.
 * \param refs pointer toward the table of the triangle references.
 * refs[i-1] is the ref of the \f$i^{th}\f$ tria.
 * \return 0   if failed, 1 otherwise.
 *
 * Set vertices and references of the mesh triangles (wrapper for MMG3D
 * function).
 *
 * \remark Fortran interface: (commentated in
 * order to allow to pass \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE PMMG_SET_TRIANGLES(parmesh,tria,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)    :: parmesh\n
 * > !    INTEGER,DIMENSION(*), INTENT(IN) :: tria,refs\n
 * > !    INTEGER, INTENT(OUT)             :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int PMMG_Set_triangles(PMMG_pParMesh parmesh, int *tria, int *refs);

/**
 * \param parmesh pointer toward the group structure.
 * \param v0  first vertex of quadrilateral.
 * \param v1  second vertex of quadrilateral.
 * \param v2  third vertex of quadrilateral.
 * \param v3  fourth vertex of quadrilateral.
 * \param ref quadrilateral reference.
 * \param pos quadrilateral position in the mesh.
 * \return 0  if failed, 1 otherwise.
 *
 * Set quadrilateral of vertices \a v0, \a v1, \a v2, \a v3 and reference \a ref
 * at position \a pos in mesh structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_QUADRILATERAL(parmesh,v0,v1,v2,v3,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: v0,v1,v2,v3,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_quadrilateral(PMMG_pParMesh parmesh, int v0, int v1, int v2, int v3,
                           int ref,int pos);

/**
 * \param parmesh   pointer toward the group structure.
 * \param quads pointer toward the table of the quads vertices
 * Vertices of the \f$i^{th}\f$ quad are stored in quads[(i-1)*3]\@3.
 * \param refs  pointer toward the table of the quadrilateral references.
 * refs[i-1] is the ref of the \f$i^{th}\f$ quad.
 * \return 0    if failed, 1 otherwise.
 *
 * Set vertices and references of the mesh quadrilaterals
 * (wrapper for MMG3D function).
 *
 * \remark Fortran interface: (commentated in
 * order to allow to pass \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE PMMG_SET_QUADRILATERALS(parmesh,quads,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)    :: parmesh\n
 * > !    INTEGER,DIMENSION(*), INTENT(IN) :: quads,refs\n
 * > !    INTEGER, INTENT(OUT)             :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int PMMG_Set_quadrilaterals(PMMG_pParMesh parmesh, int *quads, int *refs);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param v0  first extremity of the edge.
 * \param v1  second extremity of the edge.
 * \param ref edge reference.
 * \param pos edge position in the mesh.
 * \return 0  if failed, 1 otherwise.
 *
 * Set edges of extremities \a v0, \a v1 and reference \a ref at
 * position \a pos in mesh structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_EDGE(parmesh,v0,v1,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: v0,v1,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_edge(PMMG_pParMesh parmesh, int v0, int v1, int ref, int pos);

/**
 * \param parmesh pointer toward the group structure.
 * \param k   vertex index.
 * \return 1.
 *
 * Set corner at point \a pos.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_CORNER(parmesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_corner(PMMG_pParMesh parmesh, int k);

/**
 * \param parmesh pointer toward the group structure.
 * \param k   vertex index.
 * \return 1.
 *
 * Set point \a k as required (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_REQUIREDVERTEX(parmesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_requiredVertex(PMMG_pParMesh parmesh, int k);

/**
 * \param parmesh pointer toward the group structure.
 * \param k   element index.
 * \return 1.
 *
 * Set element \a k as required (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_REQUIREDTETRAHEDRON(parmesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_requiredTetrahedron(PMMG_pParMesh parmesh, int k);

/**
 * \param parmesh    pointer toward the group structure.
 * \param reqIdx table of the indices of the required elements.
 * \param nreq   number of required elements
 * \return 1.
 *
 * Set the required Tetra (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_REQUIREDTETRAHEDRA(parmesh,reqIdx,nreq,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)    :: parmesh\n
 * >     INTEGER, DIMENSION(*),INTENT(IN) :: reqIdx\n
 * >     INTEGER, INTENT(IN)              :: nreq\n
 * >     INTEGER, INTENT(OUT)             :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_requiredTetrahedra(PMMG_pParMesh parmesh, int *reqIdx, int nreq);

/**
 * \param parmesh pointer toward the group structure.
 * \param k   triangle index.
 * \return 1.
 *
 * Set triangle \a k as required (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_REQUIREDTRIANGLE(parmesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_requiredTriangle(PMMG_pParMesh parmesh, int k);

/**
 * \param parmesh    pointer toward the group structure.
 * \param reqIdx table of the indices of the required trias.
 * \param nreq   number of required trias
 * \return 1.
 *
 * Set the required triangles (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_REQUIREDTRIANGLES(parmesh,reqIdx,nreq,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)    :: parmesh\n
 * >     INTEGER, DIMENSION(*),INTENT(IN) :: reqIdx\n
 * >     INTEGER, INTENT(IN)              :: nreq\n
 * >     INTEGER, INTENT(OUT)             :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_requiredTriangles(PMMG_pParMesh parmesh, int *reqIdx, int nreq);

/**
 * \param parmesh pointer toward the group structure.
 * \param k   edge index.
 * \return 1.
 *
 * Set ridge at edge \a k (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_RIDGE(parmesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_ridge(PMMG_pParMesh parmesh, int k);

/**
 * \param parmesh pointer toward the group structure.
 * \param k   edge index.
 * \return 1.
 *
 * Set edge \a k as required (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_REQUIREDEDGE(parmesh,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_requiredEdge(PMMG_pParMesh parmesh, int k);

/**
 * \param parmesh pointer toward the group structure.
 * \param k   point index
 * \param n0  x componant of the normal at point \a k.
 * \param n1  y componant of the normal at point \a k.
 * \param n2  z componant of the normal at point \a k.
 *
 * \return 1 if success.
 *
 * Set normals (n0,n1,n2) at point \a k (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_NORMALATVERTEX(parmesh,k,n0,n1,n2,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     REAL(KIND=8), INTENT(IN)      :: n0,n1,n2\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_normalAtVertex(PMMG_pParMesh parmesh, int k, double n0, double n1,
                              double n2);

/**
 * \param parmesh pointer toward a parmesh structure
 * \param i position of the solution field that we want to set.
 * \param s table of the solutions at mesh vertices. The solution at vertex \a k
 * is given by s[k-1] for a scalar sol, s[3*(k-1)]\@3 for a vectorial solution
 * and s[6*(k-1)]\@6 for a tensor solution.
 * \param pos position of the vertex to which applies the solution.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Set values of the ith solution array at the vertex \a pos.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_ITHSOL_INSOLSATVERTICES(parmesh,i,s,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: i,pos\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Set_ithSol_inSolsAtVertices(PMMG_pParMesh parmesh,int i, double* s,int pos);

/**
 * \param parmesh pointer toward a parmesh structure
 * \param i position of the solution field that we want to set.
 * \param s table of the solutions at mesh vertices. The solution at vertex \a k
 * is given by s[k-1] for a scalar sol, s[3*(k-1)]\@3 for a vectorial solution
 * and s[6*(k-1)]\@6 for a tensor solution.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Set values of the solution at the ith field of the solution array.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_ITHSOLS_INSOLSATVERTICES(parmesh,i,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: i\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Set_ithSols_inSolsAtVertices(PMMG_pParMesh parmesh,int i, double* s);

/**
 * \param parmesh pointer toward the group structure.
 * \param m   metrics scalar value.
 * \param pos position of the solution in the mesh.
 * \return 0  if failed, 1 otherwise.
 *
 * Set scalar value \a m at position \a pos in metrics structure (wrapper
 * for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_SCALARMET(parmesh,m,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     REAL(KIND=8), INTENT(IN)      :: m\n
 * >     INTEGER, INTENT(IN)           :: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_scalarMet(PMMG_pParMesh parmesh, double m, int pos);

/**
 * \param parmesh pointer toward the group structure.
 * \param m   table of the metrics solutions values.
 * m[i-1] is the solution at vertex i.
 * \return 0  if failed, 1 otherwise.
 *
 * Set scalar metrics at mesh vertices (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_SCALARMETS(parmesh,m,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)         :: parmesh\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: m\n
 * >     INTEGER, INTENT(OUT)                  :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_scalarMets(PMMG_pParMesh parmesh, double *m);

/**
 * \param parmesh pointer toward the group structure.
 * \param vx  x value of the vectorial solution.
 * \param vy  y value of the vectorial solution.
 * \param vz  z value of the vectorial solution.
 * \param pos position of the solution in the mesh (begin to 1).
 * \return 0  if failed, 1 otherwise.
 *
 * Set vectorial value \f$(v_x,v_y,v_z)\f$ at position \a pos in metrics
 * structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_VECTORMET(parmesh,vx,vy,vz,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     REAL(KIND=8), INTENT(IN)      :: vx,vy,vz\n
 * >     INTEGER, INTENT(IN)           :: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_vectorMet(PMMG_pParMesh parmesh, double vx,double vy, double vz,
                        int pos);

/**
 * \param parmesh  pointer toward the group structure.
 * \param mets table of the vectorial metrics
 * mets[3*(i-1)]\@3 is the solution at vertex i
 * \return 0   if failed, 1 otherwise.
 *
 * Set vectorial solutions at mesh vertices (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_VECTORMETS(parmesh,mets,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)         :: parmesh\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: mets\n
 * >     INTEGER, INTENT(OUT)                  :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_vectorMets(PMMG_pParMesh parmesh, double *mets);

/**
 * \param parmesh pointer toward the group structure.
 * \param m11 value of the tensorial metrics at position (1,1) in the tensor
 * \param m12 value of the tensorial metrics at position (1,2) in the tensor
 * \param m13 value of the tensorial metrics at position (1,3) in the tensor
 * \param m22 value of the tensorial metrics at position (2,2) in the tensor
 * \param m23 value of the tensorial metrics at position (2,3) in the tensor
 * \param m33 value of the tensorial metrics at position (3,3) in the tensor
 * \param pos position of the metrics in the mesh (begin to 1).
 * \return 0 if failed, 1 otherwise.
 *
 * Set tensorial values at position \a pos in metrics
 * structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_TENSORMET(parmesh,m11,m12,m13,m22,m23,m33,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     REAL(KIND=8), INTENT(IN)      :: m11,m12,m13,m22,m23,m33\n
 * >     INTEGER, INTENT(IN)           :: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_tensorMet(PMMG_pParMesh parmesh, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos);

/**
 * \param parmesh  pointer toward the group structure.
 * \param mets table of the tensorial solutions.
 * mets[6*(i-1)]\@6 is the solution at vertex i
 * \return 0   if failed, 1 otherwise.
 *
 * Set tensorial values at position \a pos in solution
 * structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_TENSORMETS(parmesh,mets,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)         :: parmesh\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: mets\n
 * >     INTEGER, INTENT(OUT)                  :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_tensorMets(PMMG_pParMesh parmesh, double *mets);

/**
 * \param parmesh   Pointer towards the parmesh structure.
 * \param np        Number of vertices.
 * \param ne        Number of tetrahedra.
 * \param nprism    Number of prisms.
 * \param nt        Number of triangles.
 * \param nquad     Number of quads.
 * \param na        Number of edges.
 * \return          0 if failed, 1 otherwise.
 *
 * Get the number of vertices, tetrahedra, prisms, triangles, quadrilaterals and
 * edges of the mesh and allocate the associated tables.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_MESHSIZE(parmesh, np, ne, nprism, nt, nquad, na, &\n
 * >                                retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER                       :: np, ne, nprism, nt, nquad, na\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Get_meshSize(PMMG_pParMesh parmesh,int *np,int *ne,int *nprism,int *nt,
                      int *nquad, int *na);
/**
 * \param parmesh   Pointer towards the parmesh structure.
 * \param nsols     Number of solutions per entity.
 * \param nentities Number of entities.
 * \param typSol    Array of size nsol listing the type of the solutions
 *                  (scalar, vectorial...).
 * \return          0 if failed, 1 otherwise.
 *
 * Get the number of solutions per entity, the type of entity to which applies
 * each solution and the type of each solution.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_SOLSATVERTICESSIZE(parmesh, nsols,nentities,typSol,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER                       :: nsols,nentities\n
 * >     INTEGER                       :: typSol(*)\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Get_solsAtVerticesSize(PMMG_pParMesh parmesh,int *nsols,
                                int *nentities,int *typSol);

/**
 * \param parmesh   Pointer towards the parmesh structure.
 * \param typEntity Type of entity to which applies the metric (vertices, triangles...).
 * \param np        number of vertices
 * \param typMet    Type of the metric (scalar, vectorial...).
 * \return          0 if failed, 1 otherwise.
 *
 * Get the type of entity to which applies the metric and the type of the metric.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_METSIZE(parmesh,typEntity,np,typMet,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER                       :: np\n
 * >     INTEGER                       :: typEntity, typMet\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Get_metSize(PMMG_pParMesh parmesh,int *typEntity,int *np,int *typMet);

/**
 * \param parmesh pointer toward the group structure.
 * \param c0  pointer toward the coordinate of the point along the first
 * dimension.
 * \param c1  pointer toward the coordinate of the point along the second
 * dimension.
 * \param c2  pointer toward the coordinate of the point along the third
 * dimension.
 * \param ref pointer to the point reference.
 * \param isCorner   pointer toward the flag saying if point is corner.
 * \param isRequired pointer toward the flag saying if point is required.
 * \return 1.
 *
 * Get coordinates \a c0, \a c1,\a c2 and reference \a ref of next
 * vertex of mesh (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_VERTEX(parmesh,c0,c1,c2,ref,isCorner,isRequired, &\n
 * >                               retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     REAL(KIND=8), INTENT(OUT)     :: c0,c1,c2\n
 * >     INTEGER                       :: ref,isCorner,isRequired\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Get_vertex(PMMG_pParMesh parmesh, double* c0, double* c1, double* c2,
                     int* ref,int* isCorner, int* isRequired);

/**
 * \param parmesh      pointer toward the group structure.
 * \param vertices pointer toward the table of the points coordinates.
 * The coordinates of the \f$i^{th}\f$ point are stored in
 * vertices[(i-1)*3]\@3.
 * \param refs     pointer to the table of the point references.
 * The ref of the \f$i^th\f$ point is stored in refs[i-1].
 * \param areCorners pointer toward the table of the flags saying if
 * points are corners.
 * areCorners[i-1]=1 if the \f$i^{th}\f$ point is corner.
 * \param areRequired pointer toward the table of flags saying if points
 * are required. areRequired[i-1]=1 if the \f$i^{th}\f$ point is required.
 * \return 1.
 *
 * Get the coordinates and references of the mesh vertices (wrapper for
 * MMG3D function).
 *
 * \remark Fortran interface: (commentated in order to allow to pass \%val(0)
 * instead of the refs, areCorners or areRequired arrays)
 *
 * > ! SUBROUTINE PMMG_GET_VERTICES(parmesh,vertices,refs,areCorners,&\n
 * > !                               areRequired,retval)\n
 * > !   MMG5_DATA_PTR_T,INTENT(INOUT)          :: parmesh\n
 * > !   REAL(KIND=8),DIMENSION(*), INTENT(OUT) :: vertices\n
 * > !   INTEGER, DIMENSION(*)                  :: refs,areCorners,areRequired\n
 * > !   INTEGER, INTENT(OUT)                   :: retval\n
 * > ! END SUBROUTINE\n
 *
 */
int  PMMG_Get_vertices(PMMG_pParMesh parmesh, double* vertices, int* refs,
                       int* areCorners, int* areRequired);

/**
 * \param parmesh pointer toward the group structure.
 * \param v0  pointer toward the first vertex of tetrahedron.
 * \param v1  pointer toward the second vertex of tetrahedron.
 * \param v2  pointer toward the third vertex of tetrahedron.
 * \param v3  pointer toward the fourth vertex of tetrahedron.
 * \param ref pointer toward the tetrahedron reference.
 * \param isRequired pointer toward the flag saying if tetrahedron is
 *  required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices \a v0, \a v1, \a v2, \a v3 and reference \a ref of
 * next tetra of mesh (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_TETRAHEDRON(parmesh,v0,v1,v2,v3,ref,isRequired,&\n
 * >                                    retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: v0,v1,v2,v3\n
 * >     INTEGER                       :: ref,isRequired\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Get_tetrahedron(PMMG_pParMesh parmesh, int* v0, int* v1, int* v2,
                           int* v3,int* ref, int* isRequired);

/**
 * \param parmesh   pointer toward the group structure.
 * \param tetra pointer toward the table of the tetrahedra vertices.
 * Vertices of the \f$i^{th}\f$ tetra are stored in tetra[(i-1)*4]\@4.
 * \param refs  pointer toward the table of the tetrahedron references.
 * References of the \f$i^{th}\f$ tetra is stored in refs[i-1].
 * \param areRequired pointer toward the table of the flags saying if the
 *  tetrahedra are required. areRequired[i-1]=1 if the \f$i^{th}\f$ tetra
 * is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices and references of the mesh tetrahedra (wrapper for MMG3D
 * function).
 *
 * \remark Fortran interface: (commentated in order to allow to pass \%val(0)
 * instead of the refs, areCorners or areRequired arrays)
 *
 * > !  SUBROUTINE PMMG_GET_TETRAHEDRA(parmesh,tetra,refs,areRequired,&\n
 * > !                                   retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * > !    INTEGER, DIMENSION(*),INTENT(OUT) :: tetra\n
 * > !    INTEGER, DIMENSION(*)         :: refs,areRequired\n
 * > !    INTEGER, INTENT(OUT)          :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int  PMMG_Get_tetrahedra(PMMG_pParMesh parmesh, int* tetra,int* refs,
                          int* areRequired);

/**
 * \param parmesh pointer toward the mesh structure.
 * \param v0  pointer toward the first vertex of prism.
 * \param v1  pointer toward the second vertex of prism.
 * \param v2  pointer toward the third vertex of prism.
 * \param v3  pointer toward the fourth vertex of prism.
 * \param v4  pointer toward the fifth vertex of prism.
 * \param v5  pointer toward the sixth vertex of prism.
 * \param ref pointer toward the prism reference.
 * \param isRequired pointer toward the flag saying if prism is
 *  required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices \a v0, \a v1, \a v2, \a v3, \a v4, \a v5 and reference \a ref of
 * next prism of mesh (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_PRISM(parmesh,v0,v1,v2,v3,v4,v5,ref,isRequired,&\n
 * >                                    retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: v0,v1,v2,v3,v4,v5\n
 * >     INTEGER                       :: ref,isRequired\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Get_prism(PMMG_pParMesh parmesh, int* v0, int* v1, int* v2,
                    int* v3, int* v4, int* v5, int* ref, int* isRequired);

/**
 * \param parmesh    pointer toward the group structure.
 * \param prisms pointer toward the table of the prisms vertices.
 * Vertices of the \f$i^{th}\f$ prism are stored in prisms[(i-1)*6]\@6.
 * \param refs   pointer toward the table of the prism references.
 * References of the \f$i^{th}\f$ prism is stored in refs[i-1].
 * \param areRequired pointer toward the table of the flags saying if the
 *  prisms are required. areRequired[i-1]=1 if the \f$i^{th}\f$ prism
 * is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices and references of the mesh prisms (wrapper for MMG3D function).
 *
 * \remark Fortran interface: (commentated in order to allow to pass \%val(0)
 * instead of the refs, areCorners or areRequired arrays)
 *
 * > !  SUBROUTINE PMMG_GET_PRISMS(parmesh,prisms,refs,areRequired,&\n
 * > !                              retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)     :: parmesh\n
 * > !    INTEGER, DIMENSION(*),INTENT(OUT) :: prisms\n
 * > !    INTEGER, DIMENSION(*)             :: refs,areRequired\n
 * > !    INTEGER, INTENT(OUT)              :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int  PMMG_Get_prisms(PMMG_pParMesh parmesh, int* prisms,int* refs,
                          int* areRequired);

/**
 * \param parmesh pointer toward the group structure.
 * \param v0  pointer toward the first vertex of triangle.
 * \param v1  pointer toward the second vertex of triangle.
 * \param v2  pointer toward the third vertex of triangle.
 * \param ref pointer toward the triangle reference.
 * \param isRequired pointer toward the flag saying if triangle is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices \a v0,\a v1,\a v2 and reference \a ref of next
 * triangle of mesh (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_TRIANGLE(parmesh,v0,v1,v2,ref,isRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: v0,v1,v2\n
 * >     INTEGER                       :: ref,isRequired\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Get_triangle(PMMG_pParMesh parmesh, int* v0, int* v1, int* v2, int* ref,
                       int* isRequired);

/**
 * \param parmesh  pointer toward the group structure.
 * \param tria pointer toward the table of the triangles vertices
 * Vertices of the \f$i^{th}\f$ tria are stored in tria[(i-1)*3]\@3.
 * \param refs pointer toward the table of the triangles references.
 * refs[i-1] is the ref of the \f$i^{th}\f$ tria.
 * \param areRequired pointer toward table of the flags saying if triangles
 * are required. areRequired[i-1]=1 if the \f$i^{th}\f$ tria
 * is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices and references of the mesh triangles (wrapper for MMG3D function).
 *
 * \remark Fortran interface: (Commentated in order to allow to pass \%val(0)
 * instead of the refs or areRequired arrays)
 *
 * > !  SUBROUTINE PMMG_GET_TRIANGLES(parmesh,tria,refs,areRequired,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)     :: parmesh\n
 * > !    INTEGER, DIMENSION(*),INTENT(OUT) :: tria\n
 * > !    INTEGER, DIMENSION(*)             :: refs,areRequired\n
 * > !    INTEGER, INTENT(OUT)              :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int  PMMG_Get_triangles(PMMG_pParMesh parmesh, int* tria, int* refs,
                        int* areRequired);

/**
 * \param parmesh pointer toward the group structure.
 * \param v0  pointer toward the first vertex of quadrilateral.
 * \param v1  pointer toward the second vertex of quadrilateral.
 * \param v2  pointer toward the third vertex of quadrilateral.
 * \param v3  pointer toward the fourth vertex of quadrilateral.
 * \param ref pointer toward the quadrilateral reference.
 * \param isRequired pointer toward the flag saying if quadrilateral is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices \a v0,\a v1,\a v2,\a v3 and reference \a ref of next
 * quadrilateral of mesh (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_QUADRILATERAL(parmesh,v0,v1,v2,v3,ref,isRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: v0,v1,v2,v3\n
 * >     INTEGER                       :: ref,isRequired\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Get_quadrilateral(PMMG_pParMesh parmesh, int* v0, int* v1, int* v2,int* v3,
                            int* ref, int* isRequired);

/**
 * \param parmesh   pointer toward the group structure.
 * \param quads pointer toward the table of the quadrilaterals vertices
 * Vertices of the \f$i^{th}\f$ quad are stored in tria[(i-1)*4]\@4.
 * \param refs  pointer toward the table of the quadrilaterals references.
 * refs[i-1] is the ref of the \f$i^{th}\f$ quad.
 * \param areRequired pointer toward table of the flags saying if quadrilaterals
 * are required. areRequired[i-1]=1 if the \f$i^{th}\f$ quad
 * is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices and references of the mesh quadrilaterals (wrapper for MMG3D
 * function).
 *
 * \remark Fortran interface: (Commentated in order to allow to pass \%val(0)
 * instead of the refs or areRequired arrays)
 *
 * > !  SUBROUTINE PMMG_GET_QUADRILATERALS(parmesh,quads,refs,areRequired,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)     :: parmesh\n
 * > !    INTEGER, DIMENSION(*),INTENT(OUT) :: quads\n
 * > !    INTEGER, DIMENSION(*)             :: refs,areRequired\n
 * > !    INTEGER, INTENT(OUT)              :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int  PMMG_Get_quadrilaterals(PMMG_pParMesh parmesh, int* quads, int* refs,
                             int* areRequired);

/**
 * \param parmesh pointer toward the group structure.
 * \param e0  pointer toward the first extremity of the edge.
 * \param e1  pointer toward the second  extremity of the edge.
 * \param ref pointer toward the edge reference.
 * \param isRidge pointer toward the flag saying if the edge is ridge.
 * \param isRequired pointer toward the flag saying if the edge is required.
 * \return 0 if failed, 1 otherwise.
 *
 * Get extremities \a e0, \a e1 and reference \a ref of next edge of mesh
 * (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_EDGE(parmesh,e0,e1,ref,isRidge,isRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: e0,e1\n
 * >     INTEGER                       :: ref,isRidge,isRequired\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Get_edge(PMMG_pParMesh parmesh, int* e0, int* e1, int* ref,
                   int* isRidge, int* isRequired);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param edges pointer toward the array of edges.
 * Vertices of the \f$i^{th}\f$ edge are stored in edge[(i-1)*2]\@2.
 * \param refs edges references. refs[i-1] is the ref of the \f$i^{th}\f$ edge.
 * \return 0 if failed, 1 otherwise.
 *
 * Set vertices and references of the mesh edges.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_EDGES(parmesh,edges,refs,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: edges(*),refs(*)\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Set_edges(PMMG_pParMesh parmesh, int *edges, int* refs);
/**
 * \param parmesh pointer toward the mesh structure.
 * \param edges pointer toward the array of edges.
 * Vertices of the \f$i^{th}\f$ edge are stored in edge[(i-1)*2]\@2.
 * \param refs edges references. refs[i-1] is the ref of the \f$i^{th}\f$ edge.
 * \param areRidges 1 if the edge is a ridge, 0 otherwise.
 * \param areRequired 1 if the edge is required, 0 otherwise.
 * \return 0 if failed, 1 otherwise.
 *
 * Get vertices and references of the mesh edges.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_EDGES(parmesh,edges,refs,areRidges,areRequired,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: edges(*)\n
 * >     INTEGER, INTENT(OUT)          :: refs(*),areRequired(*),areRidges(*)\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Get_edges(PMMG_pParMesh parmesh,int *edges,int* refs,
                     int *areRidges,int *areRequired);
/**
 * \param parmesh pointer toward the mesh structure.
 * \param k   point index
 * \param n0  x componant of the normal at point \a k.
 * \param n1  y componant of the normal at point \a k.
 * \param n2  z componant of the normal at point \a k.
 *
 * \return 1 if success.
 *
 * Get normals (n0,n1,n2) at point \a k (wrapper for PMMG function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_NORMALATVERTEX(parmesh,k,n0,n1,n2,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     REAL(KIND=8)                  :: n0,n1,n2\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Get_normalAtVertex(PMMG_pParMesh parmesh, int k, double *n0, double *n1,
                             double *n2) ;

/**
 * \param parmesh pointer toward the array of solutions
 * \param i position of the solution field that we want to get.
 * \param s table of the solutions at mesh vertices. The solution at vertex \a k
 * is given by s[k-1] for a scalar sol, s[3*(k-1)]\@3 for a vectorial solution
 * and s[6*(k-1)]\@6 for a tensor solution.
 * \param pos position of the vertew towhich applies the solution.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Get values of the ith field of the solution array at vertex \a pos
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_ITHSOL_INSOLSATVERTICES(parmesh,i,s,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: i,pos\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int  PMMG_Get_ithSol_inSolsAtVertices(PMMG_pParMesh parmesh,int i, double* s,int pos);


/**
 * \param parmesh pointer toward the array of solutions
 * \param i position of the solution field that we want to get.
 * \param s table of the solutions at mesh vertices. The solution at vertex \a k
 * is given by s[k-1] for a scalar sol, s[3*(k-1)]\@3 for a vectorial solution
 * and s[6*(k-1)]\@6 for a tensor solution.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Get values of the solution at the ith field of the solution array.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_ITHSOLS_INSOLSATVERTICES(parmesh,i,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: i\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: s\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int  PMMG_Get_ithSols_inSolsAtVertices(PMMG_pParMesh parmesh,int i, double* s);

/**
 * \param parmesh pointer toward the group structure.
 * \param m   pointer toward the scalar metrics value.
 * \return 0  if failed, 1 otherwise.
 *
 * Get solution \a m of next vertex of mesh (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_SCALARMET(parmesh,m,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     REAL(KIND=8), INTENT(OUT)     :: m\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Get_scalarMet(PMMG_pParMesh parmesh, double* m);

/**
 * \param parmesh pointer toward the group structure.
 * \param m   table of the scalar solutions at mesh vertices. s[i-1] is
 * the solution at vertex i.
 * \return 0  if failed, 1 otherwise.
 *
 * Get metrics at mesh vertices (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_SCALARMETS(parmesh,m,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: parmesh\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: m\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG_Get_scalarMets(PMMG_pParMesh parmesh, double* m);

/**
 * \param parmesh pointer toward the group structure.
 * \param vx  x value of the vectorial solution.
 * \param vy  y value of the vectorial solution.
 * \param vz  z value of the vectorial solution.
 * \return 0  if failed, 1 otherwise.
 *
 * Get vectorial metrics \f$(v_x,v_y,vz)\f$ of next vertex of mesh
 * (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_VECTORMET(parmesh,vx,vy,vz,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     REAL(KIND=8), INTENT(OUT)     :: vx,vy,vz\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Get_vectorMet(PMMG_pParMesh parmesh, double* vx, double* vy, double* vz);

/**
 * \param parmesh  pointer toward the group structure.
 * \param mets table of the solutions at mesh vertices. mets[3*(i-1)]\@3 is
 * the solution at vertex i.
 * \return 0   if failed, 1 otherwise.
 *
 * Get vectorial metrics at mesh vertices (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_VECTORMETS(parmesh,mets,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)          :: parmesh\n
 * >     REAL(KIND=8), DIMENSION(*),INTENT(OUT) :: mets\n
 * >     INTEGER, INTENT(OUT)                   :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Get_vectorMets(PMMG_pParMesh parmesh, double* mets);

/**
 * \param parmesh pointer toward the group structure.
 * \param m11 pointer toward the position (1,1) in the metrics tensor.
 * \param m12 pointer toward the position (1,2) in the metrics tensor.
 * \param m13 pointer toward the position (1,3) in the metrics tensor.
 * \param m22 pointer toward the position (2,2) in the metrics tensor.
 * \param m23 pointer toward the position (2,3) in the metrics tensor.
 * \param m33 pointer toward the position (3,3) in the metrics tensor.
 * \return 0  if failed, 1 otherwise.
 *
 * Get tensorial metrics of next vertex of mesh (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_TENSORMET(parmesh,m11,m12,m13,m22,m23,m33,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     REAL(KIND=8), INTENT(OUT)     :: m11,m12,m13,m22,m23,m33\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Get_tensorMet(PMMG_pParMesh parmesh, double *m11,double *m12, double *m13,
                       double *m22,double *m23, double *m33);

/**
 * \param parmesh  pointer toward the group structure.
 * \param mets table of the metrics at mesh vertices.
 * mets[6*(i-1)]\@6 is the solution at vertex i.
 * \return 0   if failed, 1 otherwise.
 *
 * Get tensorial metrics at mesh vertices (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_TENSORMETS(parmesh,mets,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)           :: parmesh\n
 * >     REAL(KIND=8), DIMENSION(*), INTENT(OUT) :: mets\n
 * >     INTEGER, INTENT(OUT)                    :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Get_tensorMets(PMMG_pParMesh parmesh, double *mets);

/* libparmmg_tools.c: Tools for the library */
/**
 * \param parmesh pointer to pmmg structure
 * \return 0 if fail, 1 if success.
 *
 * Print the default parameters values
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_DEFAULTVALUES(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_defaultValues( PMMG_pParMesh parmesh );

/**
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 1 on success
 *         0 on failure
 *
 * Parse command line arguments.
 *
 * \remark no matching fortran function.
 * \remark each proc read the parameters.
 *
 */
int PMMG_parsar( int argc, char *argv[], PMMG_pParMesh parmesh );

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 1 on success
 *         0 on failure
 *
 * Parse parameter file.
 *
 * \remark no matching fortran function.
 * \remark each proc read the file of parameters.
 *
 */
int PMMG_parsop( PMMG_pParMesh parmesh );

/**
 * \param parmesh pointer toward the parmesh structure
 * \param prog pointer toward the program name.
 * \param return 1 if success, 0 if fail.
 *
 * Print help for parmmg options.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_USAGE(parmesh,prog,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: prog\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_usage( PMMG_pParMesh parmesh, char * const prog);

/* input/output functions */
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return 1 if success, 0 or -1 otherwise
 *
 * Read distributed mesh data. Insert rank index to the mesh name.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_LOADMESH_DISTRIBUTED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_loadMesh_distributed(PMMG_pParMesh parmesh,const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return 1 if success, 0 or -1 otherwise
 *
 * Read centralized mesh data.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_LOADMESH_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_loadMesh_centralized(PMMG_pParMesh parmesh,const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return -1 data invalid, 0 no file, 1 ok.
 *
 * Load metric field. The solution file must contains only 1 solution: the
 * metric
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_LOADMET_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_loadMet_centralized(PMMG_pParMesh parmesh,const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return -1 data invalid, 0 no file, 1 ok.
 *
 * Load metric field. The solution file must contains only 1 solution: the
 * metric. Insert rank index in the file name.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_LOADMET_DISTRIBUTED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_loadMet_distributed(PMMG_pParMesh parmesh,const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return -1 data invalid, 0 no file, 1 ok.
 *
 * Load displacement field. The solution file must contains only 1 solution.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_LOADDISP_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_loadDisp_centralized(PMMG_pParMesh parmesh,const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return -1 data invalid, 0 no file, 1 ok.
 *
 * Load level-set field. The solution file must contains only 1 solution.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_LOADLS_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_loadLs_centralized(PMMG_pParMesh parmesh,const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return -1 data invalid, 0 no file, 1 ok.
 *
 * Load displacement, level-set or metric field depending on the
 * option setted. The solution file must contains only 1 solution.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_LOADSOL_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_loadSol_centralized(PMMG_pParMesh parmesh,const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return -1 data invalid, 0 no file, 1 ok.
 *
 * Load distributed displacement, level-set or metric field depending on the
 * option setted. The solution file must contains only 1 solution.
 * Insert rank index to the mesh name.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_LOADSOL_DISTRIBUTED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_loadSol_distributed(PMMG_pParMesh parmesh,const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return -1 data invalid, 0 no file, 1 ok.
 *
 * Load 1 or more solutions in a solution file at medit file format.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_LOADALLSOLS_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_loadAllSols_centralized(PMMG_pParMesh parmesh,const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return -1 data invalid, 0 no file, 1 ok.
 *
 * Load 1 or more distributed solutions in a solution file at medit file format.
 * Insert rank index in the file name.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_LOADALLSOLS_DSITRIBUTED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_loadAllSols_distributed(PMMG_pParMesh parmesh,const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename pointer toward the name of file.

 * \return 0 if failed, 1 otherwise.
 *
 * Save mesh data for a centralized mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SAVEMESH_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_saveMesh_centralized(PMMG_pParMesh parmesh, const char *filename);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename pointer toward the name of file.

 * \return 0 if failed, 1 otherwise.
 *
 * Save mesh data for a distributed mesh (the MPI rank index is added to the
 * filename)
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SAVEMESH_DISTRIBUTED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_saveMesh_distributed(PMMG_pParMesh parmesh, const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write isotropic or anisotropic metric.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SAVEMET_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_saveMet_centralized(PMMG_pParMesh parmesh, const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write isotropic or anisotropic metric of a distributed mesh (insert rank
 * index to filename).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SAVEMET_DISTRIBUTED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_saveMet_distributed(PMMG_pParMesh parmesh, const char *filename);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write 1 or more than 1 solution in a file at medit format.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SAVEALLSOLS_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_saveAllSols_centralized(PMMG_pParMesh parmesh, const char *filename);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write 1 or more than 1 solution in a file at medit format for a distributed
 * mesh (insert rank index to filename).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SAVEALLSOLS_DISTRIBUTED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_saveAllSols_distributed(PMMG_pParMesh parmesh, const char *filename);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and 0 or 1 data at pvtu Vtk file format (.pvtu extension).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SAVEPVTUMESH(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_savePvtuMesh(PMMG_pParMesh parmesh, const char * filename);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and a list of data fields at pvtu Vtk file format (.pvtu extension).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SAVEPVTUMESH_AND_ALLDATA(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_savePvtuMesh_and_allData(PMMG_pParMesh parmesh, const char * filename);

/**
 * \param parmesh pointer toward the parmesh structure
 * \param next_comm number of communicators
 * \return 0 if failed, 1 otherwise.
 *
 * Set the number of node communicators (i.e. parallel node interfaces) seen
 * by the current process.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_NUMBEROFNODECOMMUNICATORS(parmesh,next_comm,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)            :: next_comm\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Set_numberOfNodeCommunicators(PMMG_pParMesh parmesh, int next_comm);
/**
 * \param parmesh pointer toward the parmesh structure
 * \param next_comm number of communicators
 * \return 0 if failed, 1 otherwise.
 *
 * Set the number of face communicators (i.e. parallel face interfaces) seen
 * by the current process.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_NUMBEROFFACECOMMUNICATORS(parmesh,next_comm,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)            :: next_comm\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Set_numberOfFaceCommunicators(PMMG_pParMesh parmesh, int next_comm);
/**
 * \param parmesh pointer toward the parmesh structure
 * \param ext_comm_index index of the communicator
 * \param color_out rank of the outward process
 * \param nitem number of entities in the communicator
 * \return 0 if failed, 1 otherwise.
 *
 * Set size (number of entities) and color (rank of the outward process) of a
 * given node communicator.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_ITHNODECOMMUNICATORSIZE(parmesh,ext_comm_index,color_out,nitem,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)            :: ext_comm_index\n
 * >     INTEGER, INTENT(IN)            :: color_out\n
 * >     INTEGER, INTENT(IN)            :: nitem\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Set_ithNodeCommunicatorSize(PMMG_pParMesh parmesh, int ext_comm_index, int color_out, int nitem);
/**
 * \param parmesh pointer toward the parmesh structure
 * \param ext_comm_index index of the communicator
 * \param color_out rank of the outward process
 * \param nitem number of entities in the communicator
 * \return 0 if failed, 1 otherwise.
 *
 * Set size (number of entities) and color (rank of the outward process) of a
 * given face communicator.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_ITHFACECOMMUNICATORSIZE(parmesh,ext_comm_index,color_out,nitem,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)            :: ext_comm_index\n
 * >     INTEGER, INTENT(IN)            :: color_out\n
 * >     INTEGER, INTENT(IN)            :: nitem\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Set_ithFaceCommunicatorSize(PMMG_pParMesh parmesh, int ext_comm_index, int color_out, int nitem);
/**
 * \param parmesh pointer toward the parmesh structure
 * \param ext_comm_index index of the communicator
 * \param local_index array of local mesh IDs of interface entities
 * \param global_index array of global mesh IDs of interface entities
 * \param isNotOrdered flag for reordering interface entities if not already done
 * \return 0 if failed, 1 otherwise.
 *
 * Set the nodes on a parallel interface. Global numbering is used to reorder
 * interface entities if isNotOrdered is equal to 1; otherwise, entities need to
 * be listed in the same order on the two sides of the interface (isNotOrdered
 * equal to 0) and ParMmg assumes this ordering is valid, but global indices
 * are still needed by ParMmg to internally match interface faces.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_ITHNODECOMMUNICATOR_NODES(parmesh,ext_comm_index,&\n
 * >                                                 local_index,global_index,&\n
 * >                                                 isNotOrdered,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)    :: parmesh\n
 * >     INTEGER, INTENT(IN)               :: ext_comm_index\n
 * >     INTEGER, DIMENSION(*), INTENT(IN) :: local_index\n
 * >     INTEGER, DIMENSION(*), INTENT(IN) :: global_index\n
 * >     INTEGER, INTENT(IN)               :: isNotOrdered\n
 * >     INTEGER, INTENT(OUT)              :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Set_ithNodeCommunicator_nodes(PMMG_pParMesh parmesh, int ext_comm_index, int* local_index, int* global_index, int isNotOrdered);
/**
 * \param parmesh pointer toward the parmesh structure
 * \param ext_comm_index index of the communicator
 * \param local_index array of local mesh IDs of interface entities
 * \param global_index array of global mesh IDs of interface entities
 * \param isNotOrdered flag for reordering interface entities if not already done
 * \return 0 if failed, 1 otherwise.
 *
 * Set the faces on a parallel interface. Global numbering is used to reorder
 * interface entities if isNotOrdered is equal to 1; otherwise, entities need to
 * be listed in the same order on the two sides of the interface (isNotOrdered
 * equal to 0) and global indices are not read.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_ITHFACECOMMUNICATOR_FACES(parmesh,ext_comm_index,&\n
 * >                                                 local_index,global_index,&\n
 * >                                                 isNotOrdered,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)    :: parmesh\n
 * >     INTEGER, INTENT(IN)               :: ext_comm_index\n
 * >     INTEGER, DIMENSION(*), INTENT(IN) :: local_index\n
 * >     INTEGER, DIMENSION(*), INTENT(IN) :: global_index\n
 * >     INTEGER, INTENT(IN)               :: isNotOrdered\n
 * >     INTEGER, INTENT(OUT)              :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Set_ithFaceCommunicator_faces(PMMG_pParMesh parmesh, int ext_comm_index, int* local_index, int* global_index, int isNotOrdered);
/**
 * \param parmesh pointer toward the parmesh structure
 * \param next_comm number of communicators
 * \return 0 if failed, 1 otherwise.
 *
 * Get the number of node communicators (i.e. parallel node interfaces) seen
 * by the current process.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_NUMBEROFNODECOMMUNICATORS(parmesh,next_comm,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)           :: next_comm\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Get_numberOfNodeCommunicators(PMMG_pParMesh parmesh, int *next_comm);
/**
 * \param parmesh pointer toward the parmesh structure
 * \param next_comm number of communicators
 * \return 0 if failed, 1 otherwise.
 *
 * Get the number of face communicators (i.e. parallel face interfaces) seen
 * by the current process.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_NUMBEROFFACECOMMUNICATORS(parmesh,next_comm,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)           :: next_comm\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Get_numberOfFaceCommunicators(PMMG_pParMesh parmesh, int *next_comm);
/**
 * \param parmesh pointer toward the parmesh structure
 * \param ext_comm_index index of the communicator
 * \param color_out rank of the outward process
 * \param nitem number of entities in the communicator
 * \return 0 if failed, 1 otherwise.
 *
 * Get size (number of entities) and color (rank of the outward process) of a
 * given node communicator.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_ITHNODECOMMUNICATORSIZE(parmesh,ext_comm_index,color_out,nitem,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)            :: ext_comm_index\n
 * >     INTEGER, INTENT(OUT)           :: color_out\n
 * >     INTEGER, INTENT(OUT)           :: nitem\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Get_ithNodeCommunicatorSize(PMMG_pParMesh parmesh, int ext_comm_index, int *color_out, int *nitem);
/**
 * \param parmesh pointer toward the parmesh structure
 * \param ext_comm_index index of the communicator
 * \param color_out rank of the outward process
 * \param nitem number of entities in the communicator
 * \return 0 if failed, 1 otherwise.
 *
 * Get size (number of entities) and color (rank of the outward process) of a
 * given face communicator.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_ITHFACECOMMUNICATORSIZE(parmesh,ext_comm_index,color_out,nitem,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)            :: ext_comm_index\n
 * >     INTEGER, INTENT(OUT)           :: color_out\n
 * >     INTEGER, INTENT(OUT)           :: nitem\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 */
  int PMMG_Get_ithFaceCommunicatorSize(PMMG_pParMesh parmesh, int ext_comm_index, int *color_out, int *nitem);
/**
 * \param parmesh pointer toward the parmesh structure
 * \param local_index array of local mesh IDs of interface entities
 * \return 0 if failed, 1 otherwise.
 *
 * \warning Non callable from a fortran code as Fortran cannot assign a **int
 * with differing allocations on each index.
 * \ref PMMG_Get_ithNodeCommunicator_nodes should be used instead.
 *
 * Get the nodes on a parallel interface.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_NODECOMMUNICATOR_NODES(parmesh,local_index,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)       :: parmesh\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)   :: local_index\n
 * >     INTEGER, INTENT(OUT)                 :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Get_NodeCommunicator_nodes(PMMG_pParMesh parmesh, int** local_index);
/**
 * \param parmesh pointer toward the parmesh structure
 * \param ext_comm_index index of the communicator
 * \param local_index array of local mesh IDs of specified interface entities
 * \return 0 if failed, 1 otherwise.
 *
 * Get the nodes on a parallel interface for a given node communicator.
 * To be used for Fortran users in place of \ref PMMG_Get_NodeCommunicator_nodes.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_ITHNODECOMMUNICATOR_NODES(parmesh,ext_comm_index,local_index,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)       :: parmesh\n
 * >     INTEGER, INTENT(IN)                  :: ext_comm_index\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)   :: local_index\n
 * >     INTEGER, INTENT(OUT)                 :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Get_ithNodeCommunicator_nodes(PMMG_pParMesh parmesh, int ext_comm_index, int* local_index);
/**
 * \param parmesh pointer toward the parmesh structure
 * \param ext_comm_index index of the communicator
 * \param local_index array of local mesh IDs of interface entities
 * \return 0 if failed, 1 otherwise.
 *
 * Get the faces on a parallel interface.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_FACECOMMUNICATOR_FACES(parmesh,local_index,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)       :: parmesh\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)   :: local_index\n
 * >     INTEGER, INTENT(OUT)                 :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Get_FaceCommunicator_faces(PMMG_pParMesh parmesh, int** local_index);
/**
 * \param parmesh pointer toward the parmesh structure
 * \param ncomm number of input communicators
 * \param nitem array of number of items in each communicator
 * \param color array of ranks of the outward processes
 * \param local_index array of input local mesh IDs of interface entities
 * (for each interface)
 * \return 0 if failed, 1 otherwise.
 *
 * Check the node communicators built through the "Set" functions against input
 * data. It can be called only after a "dry run" of parmmg (i.e. a parallel
 * library call performing 0 remeshing iterations).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_CHECK_SET_NODECOMMUNICATORS(parmesh,ncomm,nitem,&\n
 * >                                               color,local_index,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)      :: parmesh\n
 * >     INTEGER, INTENT(IN)                 :: ncomm\n
 * >     INTEGER, DIMENSION(*), INTENT(IN)   :: nitem\n
 * >     INTEGER, DIMENSION(*), INTENT(IN)   :: color\n
 * >     INTEGER, DIMENSION(*), INTENT(IN)   :: local_index\n
 * >     INTEGER, INTENT(OUT)                :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Check_Set_NodeCommunicators(PMMG_pParMesh parmesh,int ncomm,int* nitem,
                                   int* color, int** local_index);
/**
 * \param parmesh pointer toward the parmesh structure
 * \param ncomm number of input communicators
 * \param nitem array of number of items in each input communicator
 * \param color array of ranks of the outward processes
 * \param trianodes array (nt*3) of vertices of input interface triangles
 * (for each interface)
 * \return 0 if failed, 1 otherwise.
 *
 * Check the face communicators built through the "Set" functions against input
 * data. It can be called only after a "dry run" of parmmg (i.e. a parallel
 * library call performing 0 remeshing iterations).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_CHECK_SET_FACECOMMUNICATORS(parmesh,ncomm,nitem,&\n
 * >                                               color,trianodes,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)      :: parmesh\n
 * >     INTEGER, INTENT(IN)                 :: ncomm\n
 * >     INTEGER, DIMENSION(*), INTENT(IN)   :: nitem\n
 * >     INTEGER, DIMENSION(*), INTENT(IN)   :: color\n
 * >     INTEGER, DIMENSION(*), INTENT(IN)   :: trianodes\n
 * >     INTEGER, INTENT(OUT)                :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Check_Set_FaceCommunicators(PMMG_pParMesh parmesh,int ncomm,int* nitem,
                                   int* color, int** trianodes);
/**
 * \param parmesh pointer toward the parmesh structure
 * \param ncomm_in number of input communicators
 * \param nitem_in array of number of items in each input communicator
 * \param color_in array of ranks of the outward processes
 * \param local_index_in array of input local mesh IDs of interface entities
 * (for each interface)
 * \param ncomm_out number of output communicators
 * \param nitem_out array of number of items in each output communicator
 * \param color_out array of ranks of the outward processes
 * \param local_index_out array of output local mesh IDs of interface entities
 * (for each interface)
 * \return 0 if failed, 1 otherwise.
 *
 * Check the node communicators collected through the "Get" functions against
 * input data. It can be called only after a "dry run" of parmmg (i.e. a
 * parallel library call performing 0 remeshing iterations). Entities and
 * interfaces need not be ordered.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_CHECK_GET_NODECOMMUNICATORS(parmesh,ncomm_in,nitem_in,&\n
 * >                                               color_in,local_index_in,&\n
 * >                                               ncomm_out,nitem_out,&\n
 * >                                               color_out,local_index_out,&\n
 * >                                               retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)       :: parmesh\n
 * >     INTEGER, INTENT(IN)                  :: ncomm_in\n
 * >     INTEGER, DIMENSION(*), INTENT(IN)    :: nitem_in\n
 * >     INTEGER, DIMENSION(*), INTENT(IN)    :: color_in\n
 * >     INTEGER, DIMENSION(*), INTENT(IN)    :: local_index_in\n
 * >     INTEGER, INTENT(OUT)                 :: ncomm_out\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)   :: nitem_out\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)   :: color_out\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)   :: local_index_out\n
 * >     INTEGER, INTENT(OUT)                 :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Check_Get_NodeCommunicators(PMMG_pParMesh parmesh,
                                       int ncomm_in,int* nitem_in,
                                       int* color_in, int** local_index_in,
                                       int ncomm_out,int* nitem_out,
                                       int* color_out, int** local_index_out);
/**
 * \param parmesh pointer toward the parmesh structure
 * \param ncomm_in number of input communicators
 * \param nitem_in array of number of items in each input communicator
 * \param color_in array of ranks of the outward processes
 * \param trianodes_in array (nt*3) of vertices of input interface triangles
 * (for each interface)
 * \param ncomm_out number of output communicators
 * \param nitem_out array of number of items in each output communicator
 * \param color_out array of ranks of the outward processes
 * \param trianodes_out array (nt*3)of vertices of output interface triangles
 * (for each interface)
 * \return 0 if failed, 1 otherwise.
 *
 * Check the face communicators collected through the "Get" functions against
 * input data. It can be called only after a "dry run" of parmmg (i.e. a
 * parallel library call performing 0 remeshing iterations). Entities and
 * interfaces need not be ordered.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_CHECK_GET_FACECOMMUNICATORS(parmesh,ncomm_in,nitem_in,&\n
 * >                                               color_in,trianodes_in,&\n
 * >                                               ncomm_out,nitem_out,&\n
 * >                                               color_out,trianodes_out,&\n
 * >                                               retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)       :: parmesh\n
 * >     INTEGER, INTENT(IN)                  :: ncomm_in\n
 * >     INTEGER, DIMENSION(*), INTENT(IN)    :: nitem_in\n
 * >     INTEGER, DIMENSION(*), INTENT(IN)    :: color_in\n
 * >     INTEGER, DIMENSION(*), INTENT(IN)    :: trianodes_in\n
 * >     INTEGER, INTENT(OUT)                 :: ncomm_out\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)   :: nitem_out\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)   :: color_out\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)   :: trianodes_out\n
 * >     INTEGER, INTENT(OUT)                 :: retval\n
 * >   END SUBROUTINE\n
 *
 */
 int PMMG_Check_Get_FaceCommunicators(PMMG_pParMesh parmesh,
                                       int ncomm_in,int* nitem_in,
                                       int* color_in, int** trianodes_in,
                                       int ncomm_out,int* nitem_out,
                                       int* color_out, int** trianodes_out);

/**
 * \param parmesh pointer toward parmesh structure.
 * \param idx_glob pointer to the global node numbering.
 * \param owner pointer to the rank of the process owning the node.
 * \return 1 if success, 0 if fail.
 *
 * Get global node numbering (starting from 1) and rank of the process owning
 * the node.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_VERTEXGLONUM(parmesh,idx_glob,owner,&\n
 * >                                    retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)       :: parmesh\n
 * >     INTEGER, INTENT(OUT)                 :: idx_glob\n
 * >     INTEGER, INTENT(OUT)                 :: owner\n
 * >     INTEGER, INTENT(OUT)                 :: retval\n
 * >   END SUBROUTINE\n
 */
int PMMG_Get_vertexGloNum( PMMG_pParMesh parmesh, int *idx_glob, int *owner );

/**
 * \param parmesh pointer toward parmesh structure.
 * \param idx_glob array of global nodes numbering.
 * \param owner array of ranks of processes owning each node.
 * \return 1 if success, 0 if fail.
 *
 * Get global nodes numbering (starting from 1) and ranks of processes owning
 * each node.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_VERTICESGLONUM(parmesh,idx_glob,owner,&\n
 * >                                      retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)       :: parmesh\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)   :: idx_glob\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)   :: owner\n
 * >     INTEGER, INTENT(OUT)                 :: retval\n
 * >   END SUBROUTINE\n
 */
int PMMG_Get_verticesGloNum( PMMG_pParMesh parmesh, int *idx_glob, int *owner );

/**
 * \param parmesh pointer toward parmesh structure.
 * \param idx_glob pointer to the global triangle numbering.
 * \param owner pointer to the rank of the process owning the node.
 * \return 1 if success, 0 if fail.
 *
 * Get global node numbering (starting from 1) and rank of the process owning
 * the triangle.
 * If of the triangle is simply a parallel face (but not a boundary), its owner
 * will be negative.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_TRIANGLEGLONUM(parmesh,idx_glob,owner,&\n
 * >                                    retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)       :: parmesh\n
 * >     INTEGER, INTENT(OUT)                 :: idx_glob\n
 * >     INTEGER, INTENT(OUT)                 :: owner\n
 * >     INTEGER, INTENT(OUT)                 :: retval\n
 * >   END SUBROUTINE\n
 */
int PMMG_Get_triangleGloNum( PMMG_pParMesh parmesh, int *idx_glob, int *owner );

/**
 * \param parmesh pointer toward parmesh structure.
 * \param idx_glob array of global nodes numbering.
 * \param owner array of ranks of processes owning each node.
 * \return 1 if success, 0 if fail.
 *
 * Get global nodes numbering (starting from 1) and ranks of processes owning
 * each node.
 * If of the triangle is simply a parallel face (but not a boundary), its owner
 * will be negative.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_TRIANGLESGLONUM(parmesh,idx_glob,owner,&\n
 * >                                      retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)       :: parmesh\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)   :: idx_glob\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)   :: owner\n
 * >     INTEGER, INTENT(OUT)                 :: retval\n
 * >   END SUBROUTINE\n
 */
int PMMG_Get_trianglesGloNum( PMMG_pParMesh parmesh, int *idx_glob, int *owner );

/**
 * \param parmesh pointer toward parmesh structure
 * \param owner IDs of the process owning each interface node
 * \param idx_glob global IDs of interface nodes
 * \param nunique nb of non-redundant interface nodes on current rank
 * \param ntot totat nb of non-redundant interface nodes
 *
 * Create global IDs for nodes on parallel interfaces.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_NODECOMMUNICATOR_OWNERS(parmesh,owner,idx_glob,&\n
 * >                                               nunique,ntot,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)       :: parmesh\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)   :: owner\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)   :: idx_glob\n
 * >     INTEGER, INTENT(OUT)                 :: nunique\n
 * >     INTEGER, INTENT(OUT)                 :: ntot\n
 * >     INTEGER, INTENT(OUT)                 :: retval\n
 * >   END SUBROUTINE\n
 */
int PMMG_Get_NodeCommunicator_owners(PMMG_pParMesh parmesh,int **owner,int **idx_glob,int *nunique,int *ntot);

/**
 * \param parmesh pointer toward parmesh structure
 * \param owner IDs of the process owning each interface triangle
 * \param idx_glob global IDs of interface triangles
 * \param nunique nb of non-redundant interface triangles on current rank
 * \param ntot totat nb of non-redundant interface triangles
 *
 * Create global IDs for triangles on parallel interfaces.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_FACECOMMUNICATOR_OWNERS(parmesh,owner,idx_glob,&\n
 * >                                               nunique,ntot,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)       :: parmesh\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)   :: owner\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)   :: idx_glob\n
 * >     INTEGER, INTENT(OUT)                 :: nunique\n
 * >     INTEGER, INTENT(OUT)                 :: ntot\n
 * >     INTEGER, INTENT(OUT)                 :: retval\n
 * >   END SUBROUTINE\n
 */
int PMMG_Get_FaceCommunicator_owners(PMMG_pParMesh parmesh,int **owner,int **idx_glob,int *nunique,int *ntot);


/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * Set function pointers.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SETFUNC(parmesh)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)     :: parmesh\n
 * >   END SUBROUTINE\n
 *
 */
void PMMG_setfunc( PMMG_pParMesh parmesh );

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename file name (if null, print on stdout).
 *
 * \return 0 if fail, 1 otherwise
 *
 * Print parallel communicator in ASCII format.
 * \remark Mostly for debug purposes.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_PRINTCOMMUNICATOR(parmesh,API_mode,idx_loc,idx_glob,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)      :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)        :: filename\n
 * >     INTEGER, INTENT(IN)                 :: strlen\n
 * >     INTEGER, INTENT(OUT)                :: retval\n
 * >   END SUBROUTINE\n
 */
int PMMG_printCommunicator( PMMG_pParMesh parmesh,const char *filename );

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param ktri index of the boundary triangle.
 * \param ktet pointer toward an integer that will contains the tetra index.
 * \param iface pointer toward the triangle in \a ktet.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Fill \a ktet by the indice of a tetra to which belong a boundary triangle
 * and \a iface by the indice of the triangle in the tetra.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_TETFROMTRIA(parmesh,ktri,ktet,iface,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(IN)              :: parmesh\n
 * >     INTEGER, INTENT(IN)                      :: ktri\n
 * >     INTEGER, INTENT(OUT)                     :: ktet,iface\n
 * >     INTEGER, INTENT(OUT)                     :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Get_tetFromTria(PMMG_pParMesh parmesh, int ktri, int *ktet, int *iface);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param ktri index of the boundary triangle.
 * \param ktet array of size 2 that will contain the indices of the tetra
 * (filled by the function).
 * \param iface pointer toward an array of size 2 that will contains the indices
 * of the faces of the tetras \a ktet[i] that corresponds to the boundary tria
 * \a ktri.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Fill \a ktet by the indices of the tetra to which belong a boundary triangle
 * and \a iface by the indices of the faces of the tetras that correspond to the
 * triangle. Fill ktet[1] and iface[1] by 0 if the triangle belongs to 1 tetra only.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_TETSFROMTRIA(parmesh,ktri,ktet,iface,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(IN)              :: parmesh\n
 * >     INTEGER, INTENT(IN)                      :: ktri\n
 * >     INTEGER, DIMENSION(2), INTENT(OUT)       :: ktet,iface\n
 * >     INTEGER, INTENT(OUT)                     :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG_Get_tetsFromTria(PMMG_pParMesh parmesh, int ktri, int ktet[2], int iface[2]);


#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

#endif

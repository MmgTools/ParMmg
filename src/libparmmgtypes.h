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
 * \file libparmmgtypes.h
 * \brief parmmg types and functions that must be accessible to the library users
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#ifndef _LIBPARMMGTYPES_H
#define _LIBPARMMGTYPES_H

#include "mmg/common/libmmgtypes.h"
#include "pmmgversion.h"
#include <mpi.h>


/**
 * \def PMMG_SUCCESS
 *
 * Return value for success.
 */
#define PMMG_SUCCESS       0
/**
 * \def PMMG_LOWFAILURE
 *
 * Return value if the remesh process failed but we can save a conform
 * mesh.
 */
#define PMMG_LOWFAILURE    1
/**
 * \def PMMG_STRONGFAILURE
 *
 * Return value if the remesh process failed and the mesh is
 * non-conform.
 */
#define PMMG_STRONGFAILURE 2
/**
 * \def PMMG_FAILURE
 *
 * Return value of failure, caller is to decide how to proceed further,
 * ie whether remesh process failed and/or the mesh is non-conformant
 */
#define PMMG_FAILURE  4

/**
 * \def PMMG_TMPFAILURE
 *
 * Return failure not yet handled by developers: MPI_abort is called but should
 * be removed with a clean error handling with no deadlocks.
 */
#define PMMG_TMPFAILURE PMMG_STRONGFAILURE

/**
 * \def PMMG_ARG_start
 *
 * To begin a list of variadic arguments (mandatory first arg for all our
 * variadic functions)
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG_ARG_start  1
/**
 * \def PMMG_ARG_ppParMesh
 *
 * pointer toward a pointer toward a parMesh structure
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG_ARG_ppParMesh  2
/**
 * \def PMMG_ARG_pMesh
 *
 * PMMG_pMesh structure
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG_ARG_pMesh  3
/**
 * \def PMMG_ARG_pMet
 *
 * PMMG_pSol structure storing a metric field
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG_ARG_pMet   4
/**
 * \def PMMG_ARG_pSols
 *
 * PMMG_pSol structure storing an array of solutions
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG_ARG_pSols  5
/**
 * \def PMMG_ARG_pDisp
 *
 * PMMG_pSol structure storing a displacement field
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG_ARG_pDisp  6
/**
 * \def PMMG_ARG_pLs
 *
 * PMMG_pSol structure storing level-set function
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG_ARG_pLs  7
/**
 * \def PMMG_ARG_ngroups
 *
 * Number of groups per processor.
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */

#define PMMG_ARG_ngroups  8
/**
 * \def PMMG_ARG_MPIComm
 *
 * MPI Communicator
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG_ARG_MPIComm  9
/**
 * \def PMMG_ARG_dim
 *
 * mesh dimension
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG_ARG_dim  10
/**
 * \def PMMG_ARG_end
 *
 * To end a list of variadic argument (mandatory last argument for all our
 * variadic functions)
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG_ARG_end    11

/**
 * \def PMMG_REDISTRIBUTION_graph_balancing
 *
 * Use graph redistribution to load balance the parallel mesh
 *
 */
#define PMMG_REDISTRIBUTION_graph_balancing 0

/**
 * \def PMMG_REDISTRIBUTION_ifc_displacement
 *
 * Use interface displacement to redistribute the parallel mesh
 *
 */
#define PMMG_REDISTRIBUTION_ifc_displacement 1

/**
 * \def PMMG_REDISTRIBUTION_mode
 *
 * Choose how to redistribute the parallel mesh
 *
 */
#define PMMG_REDISTRIBUTION_mode PMMG_REDISTRIBUTION_ifc_displacement

/**
 * \def PMMG_LOADBALANCING_metis
 *
 * Use metis to compute and balance the graph during the loadbalancing step
 *
 */
#define PMMG_LOADBALANCING_metis 1

/**
 * \def PMMG_LOADBALANCING_parmetis
 *
 * Use parmetis to compute and balance the graph during the loadbalancing step
 *
 */
#define PMMG_LOADBALANCING_parmetis 2

/**
 * \def PMMG_APIDISTRIB_faces
 *
 * Use parallel faces information to build communicators with distributed API
 * functions
 *
 */
#define PMMG_APIDISTRIB_faces 0

/**
 * \def PMMG_APIDISTRIB_nodes
 *
 * Use parallel nodes information to build communicators with distributed API
 * functions
 *
 */
#define PMMG_APIDISTRIB_nodes 1

/**
 * \def PMMG_UNSET
 *
 * Initialization value
 *
 */
#define PMMG_UNSET     -1

/**
 * \def PMMG_ON
 */
#define PMMG_ON     1

/**
 * \def PMMG_OFF
 */
#define PMMG_OFF     0

/**
 * \def PMMG_GAP
 *
 * Gap value for reallocation
 *
 */
#define PMMG_GAP     0.2

/**
 * \enum PMMG_IO_entities
 * \brief Type of mesh entities that are saved in/loaded from HDF5 files.
 */
enum PMMG_IO_entities {
  PMMG_IO_Vertex,
  PMMG_IO_Edge,
  PMMG_IO_Tria,
  PMMG_IO_Quad,
  PMMG_IO_Tetra,
  PMMG_IO_Prism,
  PMMG_IO_Corner,
  PMMG_IO_RequiredVertex,
  PMMG_IO_ParallelVertex,
  PMMG_IO_Ridge,
  PMMG_IO_RequiredEdge,
  PMMG_IO_ParallelEdge,
  PMMG_IO_RequiredTria,
  PMMG_IO_ParallelTria,
  PMMG_IO_RequiredQuad,
  PMMG_IO_ParallelQuad,
  PMMG_IO_RequiredTetra,
  PMMG_IO_ParallelTetra,
  PMMG_IO_Normal,
  PMMG_IO_Tangent,
  PMMG_IO_ENTITIES_size, // Number of type of entities that can be saved
  PMMG_IO_Required, // To enable/disable saving of all type of required entites
  PMMG_IO_Parallel // To enable/disable saving of all type of parallel entites
};

/**
 * Types
 */
/**
 * \struct PMMG_int_comm
 * \brief internal communicator structure.
 */
typedef struct {
  int     nitem; /*!< Nb items in the communicator */
  int*    intvalues;  /*!< Array of integer */
  double* doublevalues;  /*!< Array of double */

} PMMG_Int_comm;
typedef PMMG_Int_comm  * PMMG_pInt_comm;

/**
 * \struct PMMG_ext_comm
 * \brief external communicator structure.
 */
typedef struct {
  int        color_in;  /*!< Color of the hosting processor */
  int        color_out; /*!< Color of the remote processor */

  int        nitem; /*!< Nb items in the communicator */
  int        nitem_to_share; /*!< Nb items in the *send/recv arrays */

  int*       int_comm_index; /*!< Index of the items in the internal communicators */

  int*       itosend; /*!< Array to send the data to the remote processor */
  int*       itorecv; /*!< Array to receive the data to the remote processor */
  double*    rtosend; /*!< Array to send the data to the remote processor */
  double*    rtorecv; /*!< Array to receive the data to the remote processor */

} PMMG_Ext_comm;
typedef PMMG_Ext_comm  * PMMG_pExt_comm;

/**
 * \struct PMMG_Grp
 * \brief Grp mesh structure.
 */
typedef struct {
  MMG5_pMesh   mesh;  /*!< mesh definition : coordinates, tetra etc.. */
  MMG5_pSol    field; /*!< physical solutions defined on each point of the mesh and interpolated from init to final mesh */
  MMG5_pSol    met;   /*!< metric */
  MMG5_pSol    disp;  /*!< displacement */
  MMG5_pSol    ls;    /*!< level-set */

  /* communicators */
  int          nitem_int_node_comm;       /*!< Nb nodes of this grp in internal communicator*/
  int*         node2int_node_comm_index1; /*!< List of interface nodes (local index)*/
  int*         node2int_node_comm_index2; /*!< List of index in internal communicator (where put the interface nodes)*/

  int          nitem_int_edge_comm;/*!< Nb edges of this grp in internal communicator*/
  int*         edge2int_edge_comm_index1; /*!< List of interface edges (local index)*/
  int*         edge2int_edge_comm_index2; /*!< List of index in internal communicator (where put the interface edges)*/

  int          nitem_int_face_comm;/*!< Nb faces of this grp in internal communicator*/
  int*         face2int_face_comm_index1; /*!< List of interface faces (local index)*/
  int*         face2int_face_comm_index2; /*!< List of index in internal communicator (where put the interface faces)*/
  int          flag;
} PMMG_Grp;
typedef PMMG_Grp  * PMMG_pGrp;

/**
 * \struct PMMG_Info
 * \brief Store input parameters of the run.
 */
typedef struct {

  int imprim;  /*!< ParMmg verbosity (may be non-null only on zero rank) */
  int imprim0; /*!< ParMmg verbosity of the zero rank */
  int mem;     /*!< memory asked by user */
  int iso;     /*!< ls mode (not yet available) */
  int isosurf; /*!< ls mode on boundaries only (not yet available) */
  int lag;     /*!< lagrangian motion (not yet available) */
  int root;    /*!< MPI root rank */
  int setfem;  /*!< fem mesh (no elt with more than 1 bdy face */
  int mmg_imprim; /*!< 1 if the user has manually setted the mmg verbosity */
  int repartitioning; /*!< way to perform mesh repartitioning */
  int pure_partitioning; /*!< enable/disable pure mesh partitioning mode  */
  int ifc_layers;  /*!< nb of layers for interface displacement */
  double grps_ratio;  /*!< allowed imbalance ratio between current and demanded groups size */
  int nobalancing; /*!< switch off final load balancing */
  int loadbalancing_mode; /*!< way to perform the loadbalanding (see LOADBALANCING) */
  int contiguous_mode; /*!< force/don't force partitions contiguity */
  int metis_ratio; /*!< wanted ratio between the number of meshes and the number of metis super nodes */
  int target_mesh_size; /*!< target mesh size for Mmg */
  int API_mode; /*!< use faces or nodes information to build communicators */
  int globalNum; /*!< compute nodes and triangles global numbering in output */
  int globalVNumGot; /*!< have global nodes actually been calculated */
  int globalTNumGot; /*!< have global triangles actually been calculated */
  int fmtout; /*!< store the output format asked */
  int io_entities[PMMG_IO_ENTITIES_size]; /* Array to store entities to save in some I/O formats */
  int8_t sethmin; /*!< 1 if user set hmin, 0 otherwise (needed for multiple library calls) */
  int8_t sethmax; /*!< 1 if user set hmin, 0 otherwise (needed for multiple library calls) */
  uint8_t inputMet; /*!< 1 if User prescribe a metric or a size law */
  int npartin; /*!< nb of partitions of the mesh in the input HDF5 file */
  MPI_Comm read_comm; /*!< MPI comm containing the procs that read the mesh (HDF5 input) */
} PMMG_Info;

/**
 * \struct PMMG_overlap
 * \brief Overlap structure.
 */
typedef struct {
  int color_in;     /*!< Color of the hosting processor */
  int color_out;    /*!< Color of the remote processor */
  int np_in2out;    /*!< Nbr of points sends    from color_in to   color_out */
  int np_out2in;    /*!< Nbr of points receives on   color_in from color_out */
  int nt_in2out;    /*!< Nbr of tetra  sends    from color_in to   color_out */
  int nt_out2in;    /*!< Nbr of tetra  receives on   color_in from color_out */
  int *hash_in2out; /*!< Hash table to find pts index on color_out from pts index on color_in */
  int *hash_out2in; /*!< Hash table to find pts index on color_in  from pts index on color_out */

} PMMG_Overlap;
typedef PMMG_Overlap * PMMG_pOverlap;

/**
 * \struct PMMG_ParMesh
 * \brief ParMmg mesh structure.
 */
typedef struct {

  /* mpi info */
  MPI_Comm    comm;   /*!< Global communicator of all parmmg processes */
  int         nprocs; /*!< Number of processes in global communicator */
  int         myrank; /*!< Rank in global communicator */
  int         size_shm; /*!< Number or MPI process per Node */

  /* mem info */
  size_t    memGloMax; /*!< Maximum memory available to all structs */
  size_t    memMax; /*!< Maximum memory parmesh is allowed to allocate */
  size_t    memCur; /*!< Currently allocated memory */

  /* file names */
  char     *meshin,*meshout;
  char     *metin,*metout;
  char     *lsin,*lsout;
  char     *dispin;
  char     *fieldin,*fieldout;

  /* grp */
  int       ngrp;       /*!< Number of grp used inside the parmesh. It can
                         * differ from listgrp size (for example if inputs have
                         * been provided on a different number of processes than
                         * the ones used for computation) */
  PMMG_pGrp listgrp;    /*!< List of grp */
  int       nold_grp;       /*!< Number of old grp */
  PMMG_pGrp old_listgrp;    /*!< List of old grp */

  /* internal communicators */
  PMMG_pInt_comm  int_node_comm; /*!< Internal node communicator (only one PMMG_Int_comm, it is not an array) */
  PMMG_pInt_comm  int_edge_comm; /*!< Internal edge communicator */
  PMMG_pInt_comm  int_face_comm; /*!< Internal face communicator */

  /* external communicators */
  int            next_node_comm; /*!< Number of external node communicator */
  PMMG_pExt_comm ext_node_comm;  /*!< External communicators (in increasing order w.r. to the remote proc index) */
  int            next_edge_comm; /*!< Number of external edge communicator */
  PMMG_pExt_comm ext_edge_comm;  /*!< External communicators (in increasing order w.r. to the remote proc index) */
  int            next_face_comm; /*!< Number of external face communicator */
  PMMG_pExt_comm ext_face_comm;  /*!< External communicators (in increasing order w.r. to the remote proc index) */

  /* overlap variables */
  PMMG_pOverlap overlap; /*!<  Overlap variables */

  /* global variables */
  int            ddebug; //! Debug level
  int            iter;   //! Current adaptation iteration
  int            niter;  //! Number of adaptation iterations

  /* parameters of the run */
  PMMG_Info      info; /*!< \ref PMMG_Info structure */

} PMMG_ParMesh;
typedef PMMG_ParMesh  * PMMG_pParMesh;

#endif

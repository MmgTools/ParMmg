/**
 * \file libparmmgtypes.h
 * \brief parmmg types and functions that must be accessible to the library users
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */


#include "mmg/mmg3d/libmmgtypes.h"
#include <mpi.h>

#ifndef _LIBPARMMGTYPES_H
#define _LIBPARMMGTYPES_H


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
 * \def PMMG_UNSET
 *
 * Initialization value
 *
 */
#define PMMG_UNSET     -1

/**
 * \def PMMG_GAP
 *
 * Gap value for reallocation
 *
 */
#define PMMG_GAP     0.2

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
  MMG5_pSol    sol;  /*!< physical solutions defined on each point of the mesh */
  MMG5_pSol    met;   /*!< metric */
  MMG5_pSol    disp;  /*!< displacement */

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
  int root;    /*!< MPI root rank */
  int fem;     /*!< fem mesh (no elt with more than 1 bdy face */
  int mmg_imprim; /*!< 1 if the user has manually setted the mmg verbosity */
} PMMG_Info;


/**
 * \struct PMMG_ParMesh
 * \brief ParMmg mesh structure.
 */
typedef struct {

  /* mpi info */
  MPI_Comm    comm;   /*!< Global communicator of all parmmg processes */
  int         nprocs; /*!< Number of processes in global communicator */
  int         myrank; /*!< Rank in global communicator */

  /* verbosity */
  int         imprim;

  /* mem info */
  size_t    memGloMax; /*!< Maximum memory available to all structs */
  size_t    memMax; /*!< Maximum memory parmesh is allowed to allocate */
  size_t    memCur; /*!< Currently allocated memory */

  /* grp */
  int       ngrp;       /*!< Number of grp */
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

  /* global variables */
  int            ddebug; //! Debug level
  int            niter;  //! Number of adaptation iterations

  /* parameters of the run */
  PMMG_Info      info; /*!< \ref PMMG_Info structure */

} PMMG_ParMesh;
typedef PMMG_ParMesh  * PMMG_pParMesh;
#endif

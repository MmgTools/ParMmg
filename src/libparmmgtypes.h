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
 * \def PMMG_ARG_ngroups
 *
 * Number of groups per processor.
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG_ARG_ngroups  3
/**
 * \def PMMG_ARG_comm
 *
 * MPI Communicator
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG_ARG_comm  4
/**
 * \def PMMG_ARG_end
 *
 * To end a list of variadic argument (mandatory last argument for all our
 * variadic functions)
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG_ARG_end    5

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
 * \struct PARMMG_int_comm
 * \brief internal communicator structure.
 */
typedef struct {
  int     nitem; /*!< Nb items in the communicator */
  int*    intvalues;  /*!< Array of integer */
  double* doublevalues;  /*!< Array of double */

} PMMG_int_comm;
typedef PMMG_int_comm  * PMMG_pint_comm;

/**
 * \struct PARMMG_ext_comm
 * \brief external communicator structure.
 */
typedef struct {
  int        color_in;  /*!< Color of the hosting processor */
  int        color_out; /*!< Color of the remote processor */

  int        nitem; /*!< Nb items in the communicator */

  int*       int_comm_index; /*!< Index of the items in the internal communicators */

  int*       itosend; /*!< Array to send the data to the remote processor */
  int*       itorecv; /*!< Array to receive the data to the remote processor */
  double*    rtosend; /*!< Array to send the data to the remote processor */
  double*    rtorecv; /*!< Array to receive the data to the remote processor */
} PMMG_ext_comm;
typedef PMMG_ext_comm  * PMMG_pext_comm;

/**
 * \struct PARMMG_Grp
 * \brief Grp mesh structure.
 */
typedef struct {
  MMG5_pMesh   mesh;  /*!< mesh definition : coordinates, tetra etc.. */
  MMG5_pSol    sol;   /*!< physical solutions defined on each point of the mesh */
  MMG5_pSol    met;   /*!< metric */
  MMG5_pSol    disp;  /*!< displacement */

  /*communicators*/
  int          nitem_int_node_comm;       /*!< Nb nodes of this grp in internal communicator*/
  int*         node2int_node_comm_index1; /*!< List of interface nodes (local index)*/
  int*         node2int_node_comm_index2; /*!< List of index in internal communicator (where put the interface nodes)*/

  int          nitem_int_edge_comm;/*!< Nb edges of this grp in internal communicator*/
  int*         edge2int_edge_comm_index1; /*!< List of interface edges (local index)*/
  int*         edge2int_edge_comm_index2; /*!< List of index in internal communicator (where put the interface edges)*/

  int          nitem_int_face_comm;/*!< Nb faces of this grp in internal communicator*/
  int*         face2int_face_comm_index1; /*!< List of interface faces (local index)*/
  int*         face2int_face_comm_index2; /*!< List of index in internal communicator (where put the interface faces)*/

} PMMG_Grp;
typedef PMMG_Grp  * PMMG_pGrp;


/**
 * \struct PMMG_ParMesh
 * \brief ParMmg mesh structure.
 */
typedef struct {

  /* mpi info */
  MPI_Comm comm;   /*!< Global communicator of all parmmg processes */
  int      nprocs; /*!< Number of processes in global communicator */
  int      myrank; /*!< Rank in global communicator */

  /* mem info */
  long long memGloMax; /*!< Maximum memory available to all structs */
  long long memMax; /*!< Maximum memory parmesh is allowed to allocate */
  long long memCur; /*!< Currently allocated memory */

  /* grp */
  int       ngrp;       /*!< Number of grp */
  PMMG_pGrp listgrp;    /*!< List of grp */

  /* internal communicators */
  PMMG_pint_comm  int_node_comm; /*!< Internal node communicator (only one PMMG_int_comm, it is not an array) */
  PMMG_pint_comm  int_edge_comm; /*!< Internal edge communicator */
  PMMG_pint_comm  int_face_comm; /*!< Internal face communicator */

  /* external communicators */
  int            next_node_comm; /*!< Number of external node communicator */
  PMMG_pext_comm ext_node_comm;  /*!< External communicators (in increasing order w.r. to the remote proc index) */
  int            next_edge_comm; /*!< Number of external edge communicator */
  PMMG_pext_comm ext_edge_comm;  /*!< External communicators (in increasing order w.r. to the remote proc index) */
  int            next_face_comm; /*!< Number of external face communicator */
  PMMG_pext_comm ext_face_comm;  /*!< External communicators (in increasing order w.r. to the remote proc index) */

  /* global variables */
  int            ddebug; //! Debug level
  int            niter;  //! Number of adaptation iterations
} PMMG_ParMesh;
typedef PMMG_ParMesh  * PMMG_pParMesh;
#endif

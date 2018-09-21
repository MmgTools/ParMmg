#ifndef MPI_PMMG_H
#define MPI_PMMG_H
/**
 * \file mpi_pmmg.h
 * \brief Mpi tools that are used in different part of parMmg
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \author Nikos Pattakos (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 */
#include <mpi.h>

#define MPI_COMMUNICATORS_NODE_TAG      1000
#define MPI_CHKCOMM_NODE_TAG            2000
#define MPI_CHKCOMM_FACE_TAG            3000
#define MPI_DISTRIBUTEGRPS_MESHSIZE_TAG 4000
#define MPI_PARMESHGRPS2PARMETIS_TAG    5000
#define MPI_SENDGRP_TAG                 6000
#define MPI_SENDEXTFACECOMM_TAG         7000
#define MPI_TRANSFER_GRP_TAG            8000

#define MPI_CHECK(func_call,on_failure) do {                            \
    int mpi_ret_val;                                                    \
                                                                        \
    mpi_ret_val = func_call;                                            \
                                                                        \
    switch ( mpi_ret_val ) {                                            \
    case ( MPI_SUCCESS ):                                               \
      break;                                                            \
    case ( MPI_ERR_COMM ):                                              \
      fprintf(stderr," ## Error: %s:%d: Invalid MPI communicator\n.",   \
              __func__,__LINE__ );                                      \
      on_failure;                                                       \
      break;                                                            \
    case ( MPI_ERR_TYPE ):                                              \
      fprintf(stderr," ## Error: %s:%d: Invalid MPI datatype\n.",       \
              __func__,__LINE__ );                                      \
      on_failure;                                                       \
      break;                                                            \
    case ( MPI_ERR_COUNT ):                                             \
      fprintf(stderr," ## Error: %s:%d: Invalid MPI count argument\n.", \
              __func__,__LINE__ );                                      \
      on_failure;                                                       \
      break;                                                            \
    case ( MPI_ERR_TAG ):                                               \
      fprintf(stderr," ## Error: %s:%d: Invalid MPI tag\n.",            \
              __func__,__LINE__ );                                      \
      on_failure;                                                       \
      break;                                                            \
    case ( MPI_ERR_RANK ):                                              \
      fprintf(stderr," ## Error: %s:%d: Invalid MPI rank\n.",           \
              __func__,__LINE__ );                                      \
      on_failure;                                                       \
      break;                                                            \
    case ( MPI_ERR_OTHER ):                                             \
      fprintf(stderr," ## Error: %s:%d: MPI error.\n.",                 \
              __func__,__LINE__ );                                      \
      on_failure;                                                       \
      break;                                                            \
    default:                                                            \
      fprintf(stderr," ## Error: %s:%d: Unexpected MPI error\n.",       \
              __func__,__LINE__ );                                      \
      on_failure;                                                       \
      break;                                                            \
    }                                                                   \
  } while(0)
#endif

#define RUN_ON_ROOT_AND_BCAST(func_call,root,myrank,on_failure) do {    \
    int ierloc;                                                         \
                                                                        \
    if ( myrank == root ) {                                             \
      ierloc = func_call;                                               \
    }                                                                   \
    MPI_CHECK( MPI_Bcast(&ierloc,1,MPI_INT,root,parmesh->comm),on_failure); \
    if ( !ierloc ) {                                                    \
      on_failure;                                                       \
    }                                                                   \
                                                                        \
  } while(0)

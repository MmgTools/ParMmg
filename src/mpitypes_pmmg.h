#ifndef MPITYPES_PMMG_H
#define MPITYPES_PMMG_H
/**
 * \file mpitypes.h
 * \brief Mpi types management header file.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 */
#include <mpi_pmmg.h>
#include "libmmgtypes.h"

int PMMG_create_MPI_lightPoint(MPI_Datatype *mpi_light_point);

int PMMG_create_MPI_Point(MPI_Datatype *mpi_point);

int PMMG_create_MPI_xPoint(MPI_Datatype *mpi_xPoint);

int PMMG_create_MPI_lightTetra(MPI_Datatype *mpi_light_tetra);

int PMMG_create_MPI_Tetra(MPI_Datatype *mpi_tetra);

int PMMG_create_MPI_Edge(MPI_Datatype *mpi_edge);

int PMMG_create_MPI_Tria(MPI_Datatype *mpi_tria);

int PMMG_create_MPI_xTetra(MPI_Datatype *mpi_xtetra);

int PMMG_Free_MPI_meshDatatype( MPI_Datatype*,MPI_Datatype*,
                                MPI_Datatype*,MPI_Datatype*);

#endif

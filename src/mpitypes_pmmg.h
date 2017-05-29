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
#include <mpi.h>
#include "libmmgtypes.h"
void PMMG_create_MPI_lightPoint(MMG5_pPoint point, MPI_Datatype *mpi_light_point);

void PMMG_create_MPI_Point(MMG5_pPoint point, MPI_Datatype *mpi_point);

void PMMG_create_MPI_xPoint(MMG5_pxPoint xPoint, MPI_Datatype *mpi_xPoint);

void PMMG_create_MPI_lightTetra(MMG5_pTetra tetra, MPI_Datatype *mpi_light_tetra);

void PMMG_create_MPI_Tetra(MMG5_pTetra tetra, MPI_Datatype *mpi_tetra);

void PMMG_create_MPI_Edge(MMG5_pEdge edge, MPI_Datatype *mpi_edge);

void PMMG_create_MPI_Tria(MMG5_pTria tria, MPI_Datatype *mpi_tria);

void PMMG_create_MPI_xTetra(MMG5_pxTetra xTetra, MPI_Datatype *mpi_xtetra);
#endif

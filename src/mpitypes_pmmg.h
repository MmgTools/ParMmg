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

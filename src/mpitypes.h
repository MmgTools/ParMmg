/**
 * \file mpitypes.h
 * \brief Mpi types management header file.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 */


int PMMG_create_MPI_Point(MMG5_pPoint point, MPI_Datatype *mpi_type_point);

int PMMG_create_MPI_xPoint(MMG5_pxPoint xPoint, MPI_Datatype *mpi_type_xPoint);

int PMMG_create_MPI_Tetra(MMG5_pTetra tetra, MPI_Datatype *mpi_type_tetra);

int PMMG_create_MPI_xTetra(MMG5_pxTetra xTetra, MPI_Datatype *mpi_type_xTetra);

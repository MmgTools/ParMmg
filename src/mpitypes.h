/**
 * \file mpitypes.h
 * \brief Mpi types management header file.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 */

int PMMG_create_MPI_ParMesh(PMMG_pParMesh parmesh, MPI_Datatype *mpi_type_parmesh,
                            MPI_Datatype *mpi_type_grp, MPI_Datatype *mpi_type_mesh,
                            MPI_Datatype *mpi_type_point,MPI_Datatype *mpi_type_xpoint,
                            MPI_Datatype *mpi_type_tetra,MPI_Datatype *mpi_type_xtetra);

void PMMG_free_MPI_ParMesh(PMMG_pParMesh parmesh, MPI_Datatype *mpi_type_parmesh,
                           MPI_Datatype *mpi_type_grp, MPI_Datatype *mpi_type_mesh,
                           MPI_Datatype *mpi_type_point,MPI_Datatype *mpi_type_xpoint,
                           MPI_Datatype *mpi_type_tetra,MPI_Datatype *mpi_type_xtetra );

int PMMG_create_MPI_Int_Comm(PMMG_pint_comm int_comm, MPI_Datatype *mpi_type_int_comm);

void PMMG_free_MPI_Int_Comm(MPI_Datatype *mpi_type_int_comm);

int PMMG_create_MPI_Ext_Comm(PMMG_pext_comm ext_comm, MPI_Datatype *mpi_type_ext_comm);

void PMMG_free_MPI_Ext_Comm(MPI_Datatype *mpi_type_ext_comm);

int PMMG_create_MPI_Grp(PMMG_pGrp grp, MPI_Datatype *mpi_type_grp,
                        MPI_Datatype *mpi_type_mesh,
                        MPI_Datatype *mpi_type_point,MPI_Datatype *mpi_type_xpoint,
                        MPI_Datatype *mpi_type_tetra,MPI_Datatype *mpi_type_xtetra);

void PMMG_free_MPI_Grp(MPI_Datatype *mpi_type_grp,
                       MPI_Datatype *mpi_type_mesh,
                       MPI_Datatype *mpi_type_point,MPI_Datatype *mpi_type_xpoint,
                       MPI_Datatype *mpi_type_tetra,MPI_Datatype *mpi_type_xtetra);

int PMMG_create_MPI_Mesh(MMG5_pMesh mesh, MPI_Datatype *mpi_type_mesh,
                         MPI_Datatype *mpi_type_point,MPI_Datatype *mpi_type_xpoint,
                         MPI_Datatype *mpi_type_tetra,MPI_Datatype *mpi_type_xtetra);

void PMMG_free_MPI_Mesh(MPI_Datatype *mpi_type_mesh,
                        MPI_Datatype *mpi_type_point,MPI_Datatype *mpi_type_xpoint,
                        MPI_Datatype *mpi_type_tetra,MPI_Datatype *mpi_type_xtetra);

int PMMG_create_MPI_Point(MMG5_pPoint point, MPI_Datatype *mpi_type_point);

void PMMG_free_MPI_Point(MPI_Datatype *mpi_type_point);

int PMMG_create_MPI_xPoint(MMG5_pxPoint xPoint, MPI_Datatype *mpi_type_xPoint);

void PMMG_free_MPI_xPoint(MPI_Datatype *mpi_type_xPoint);

int PMMG_create_MPI_Tetra(MMG5_pTetra tetra, MPI_Datatype *mpi_type_tetra);

void PMMG_free_MPI_Tetra(MPI_Datatype *mpi_type_tetra);

int PMMG_create_MPI_xTetra(MMG5_pxTetra xTetra, MPI_Datatype *mpi_type_xTetra);

void PMMG_free_MPI_xTetra(MPI_Datatype *mpi_type_xTetra);

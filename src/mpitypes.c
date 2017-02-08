/**
 * \file mpitypes.c
 * \brief Mpi types management.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 */

#include "parmmg.h"


/**
 * \param point pointer toward the point we want to communicate
 * \param mpi_type_point new MPI data type
 * \return 0 if fail, 1 otherwise.
 *
 * Create an MPI data type named mpi_type_point to allow the comminication of
 * choosen fields of a point (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */

int PMMG_create_MPI_Point(MMG5_pPoint point, MPI_Datatype *mpi_type_point) {
  int          i,blck_lengths[5] = {3, 3, 1, 1, 1};
  MPI_Aint     displs[5];
  MPI_Datatype types[5] = {MPI_DOUBLE,MPI_DOUBLE,MPI_INT,MPI_INT,MPI_UB};


  MPI_Get_address(&point->c[0],       &displs[0]);
  MPI_Get_address(&point->n,          &displs[1]);
  MPI_Get_address(&point->ref,        &displs[2]);
  MPI_Get_address(&point->xp,         &displs[3]);

  /* Relative displacement from field 0 to field i */
  for ( i=3 ; i>= 0; --i ) {
    displs[i] -= displs[0];
  }
  displs[4] = sizeof(MMG5_Point);

  MPI_Type_create_struct(5, blck_lengths, displs, types, mpi_type_point);

  MPI_Type_commit(mpi_type_point);

  return(1);
}

/**
 * \param xPoint pointer toward the xPoint we want to communicate
 * \param mpi_type_xpoint new MPI data type
 * \return 0 if fail, 1 otherwise.
 *
 * Create an MPI data type named mpi_type_xpoint to allow the comminication of
 * choosen fields of a xPoint (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */

int PMMG_create_MPI_xPoint(MMG5_pxPoint xPoint, MPI_Datatype *mpi_type_xpoint) {
  int          i,blck_lengths[3] = {3,3,1};
  MPI_Aint     displs[3];
  MPI_Datatype types[3] = {MPI_DOUBLE,MPI_DOUBLE,MPI_UB};


  MPI_Get_address(&xPoint->n1[0], &displs[0]);
  MPI_Get_address(&xPoint->n2[0], &displs[1]);

  /* Relative displacement from field 0 to field i */
  for ( i=1 ; i>= 0; --i ) {
    displs[i] -= displs[0];
  }
  displs[2] = sizeof(MMG5_xPoint);

  MPI_Type_create_struct(3, blck_lengths, displs, types, mpi_type_xpoint);

  MPI_Type_commit(mpi_type_xpoint);

  return(1);
}

/**
 * \param tetra pointer toward the tetra we want to communicate
 * \param mpi_type_tetra new MPI data type
 * \return 0 if fail, 1 otherwise.
 *
 * Create an MPI data type named mpi_type_tetra to allow the comminication of
 * choosen fields of a tetra (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */

int PMMG_create_MPI_Tetra(MMG5_pTetra tetra, MPI_Datatype *mpi_type_tetra) {
  int          i,blck_lengths[5] = {1, 4, 1, 1, 1};
  MPI_Aint     displs[5];
  MPI_Datatype types[5] = {MPI_LB,MPI_INT,MPI_INT,MPI_INT,MPI_UB};

  MPI_Get_address(&tetra->qual, &displs[0]);
  MPI_Get_address(&tetra->v[0], &displs[1]);
  MPI_Get_address(&tetra->ref,  &displs[2]);
  MPI_Get_address(&tetra->xt,   &displs[3]);

  /* Relative displacement from field 0 to field i */
  for ( i=3 ; i>= 0; --i ) {
    displs[i] -= displs[0];
  }
  displs[4] = sizeof(MMG5_Tetra);

  MPI_Type_create_struct(5, blck_lengths, displs, types, mpi_type_tetra);

  MPI_Type_commit(mpi_type_tetra);

  return(1);
}

/**
 * \param xTetra pointer toward the xTetra we want to communicate
 * \param mpi_type_xtetra new MPI data type
 * \return 0 if fail, 1 otherwise.
 *
 * Create an MPI data type named mpi_type_xtetra to allow the comminication of
 * choosen fields of a xTetra (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */

int PMMG_create_MPI_xTetra(MMG5_pxTetra xTetra, MPI_Datatype *mpi_type_xtetra) {
  MPI_Aint     displs[3];
  MPI_Datatype types[3] = {MPI_INT,MPI_INT,MPI_UB};
  int          i,blck_lengths[3] = {4,6,1};

  MPI_Get_address(&xTetra->ref[0], &displs[0]);
  MPI_Get_address(&xTetra->edg[0], &displs[1]);

 /* Relative displacement from field 0 to field i */
  for ( i=1 ; i>= 0; --i ) {
    displs[i] -= displs[0];
  }
  displs[2] = sizeof(MMG5_pxTetra);

  MPI_Type_create_struct(3, blck_lengths, displs, types, mpi_type_xtetra);

  MPI_Type_commit(mpi_type_xtetra);

  return(1);
}

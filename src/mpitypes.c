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
  MPI_Aint     displs[5],lb,ub;
  MPI_Datatype mpi_type_noextent;
  MPI_Datatype types[5] = {MPI_DOUBLE,MPI_DOUBLE,MPI_INT,MPI_INT,MPI_INT16_T};

  MPI_Get_address(&(point[0]),       &lb);
  MPI_Get_address(&(point[0].c[0]),  &displs[0]);
  MPI_Get_address(&(point[0].n),     &displs[1]);
  MPI_Get_address(&(point[0].ref),   &displs[2]);
  MPI_Get_address(&(point[0].xp),    &displs[3]);
  MPI_Get_address(&(point[0].tag),   &displs[4]);
  MPI_Get_address(&(point[1]),       &ub);

  /* Relative displacement from field 0 to field i */
  for ( i=4 ; i>= 0; --i ) {
    displs[i] -= lb;
  }

  MPI_Type_create_struct(5, blck_lengths, displs, types, &mpi_type_noextent);

  MPI_Type_create_resized(mpi_type_noextent,lb,ub-lb,mpi_type_point);

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
  int          i,blck_lengths[2] = {3,3};
  MPI_Aint     displs[2],lb,ub;
  MPI_Datatype mpi_type_noextent;
  MPI_Datatype types[2] = {MPI_DOUBLE,MPI_DOUBLE};

  MPI_Get_address(&(xPoint[0]),       &lb);
  MPI_Get_address(&(xPoint[0].n1[0]), &displs[0]);
  MPI_Get_address(&(xPoint[0].n2[0]), &displs[1]);
  MPI_Get_address(&(xPoint[1]),       &ub);

  /* Relative displacement from field 0 to field i */
  for ( i=1 ; i>= 0; --i ) {
    displs[i] -= lb;
  }

  MPI_Type_create_struct(2, blck_lengths, displs, types, &mpi_type_noextent);

  MPI_Type_create_resized(mpi_type_noextent,lb,ub-lb,mpi_type_xpoint);

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
  int          i,blck_lengths[3] = {4, 1, 1};
  MPI_Aint     displs[3],lb,ub;
  MPI_Datatype mpi_type_noextent;
  MPI_Datatype types[3] = {MPI_INT,MPI_INT,MPI_INT};

  MPI_Get_address(&(tetra[0]),      &lb);
  MPI_Get_address(&(tetra[0].v[0]), &displs[0]);
  MPI_Get_address(&(tetra[0].ref),  &displs[1]);
  MPI_Get_address(&(tetra[0].xt),   &displs[2]);
  MPI_Get_address(&(tetra[1]),      &ub);

  /* Relative displacement from field 0 to field i */
  for ( i=2 ; i>= 0; --i ) {
    displs[i] -= lb;
  }

  MPI_Type_create_struct(3, blck_lengths, displs, types, &mpi_type_noextent);

  MPI_Type_create_resized(mpi_type_noextent,lb,ub-lb,mpi_type_tetra);

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
  MPI_Aint     displs[2],lb,ub;
  MPI_Datatype mpi_type_noextent;
  MPI_Datatype types[2] = {MPI_INT,MPI_INT};
  int          i,blck_lengths[2] = {4,6};

  MPI_Get_address(&(xTetra[0])       , &lb);
  MPI_Get_address(&(xTetra[0].ref[0]), &displs[0]);
  MPI_Get_address(&(xTetra[0].edg[0]), &displs[1]);
  MPI_Get_address(&(xTetra[0])       , &ub);

 /* Relative displacement from field 0 to field i */
  for ( i=1 ; i>= 0; --i ) {
    displs[i] -= lb;
  }

  MPI_Type_create_struct(2, blck_lengths, displs, types, &mpi_type_noextent);

  MPI_Type_create_resized(mpi_type_noextent,lb,ub-lb,mpi_type_xtetra);

  MPI_Type_commit(mpi_type_xtetra);

  return(1);
}

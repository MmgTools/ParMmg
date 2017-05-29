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
#include "mpitypes_pmmg.h"

/**
 * \param point pointer toward the point we want to communicate
 * \param mpi_light_point new MPI data type
 *
 * Create an MPI data type to allow the communication of
 * choosen fields of a point (only those who needs to be communicated) when the
 * mesh hast not yet been analyzed (xp field must not be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */
void PMMG_create_MPI_lightPoint(MMG5_pPoint point, MPI_Datatype *mpi_light_point)
{
  int          i,blck_lengths[4] = {3, 3, 1, 1};
  MPI_Aint     displs[4],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[4] = {MPI_DOUBLE,MPI_DOUBLE,MPI_INT,MPI_INT16_T};

  MPI_Get_address(&(point[0]),       &lb);
  MPI_Get_address(&(point[0].c[0]),  &displs[0]);
  MPI_Get_address(&(point[0].n),     &displs[1]);
  MPI_Get_address(&(point[0].ref),   &displs[2]);
  MPI_Get_address(&(point[0].tag),   &displs[3]);
  MPI_Get_address(&(point[1]),       &ub);

  /* Relative displacement from field 0 to field i */
  for ( i=3 ; i>= 0; --i )
    displs[i] -= lb;

  MPI_Type_create_struct(4, blck_lengths, displs, types, &mpi_noextent);

  MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_light_point);

  MPI_Type_commit(mpi_light_point);
}

/**
 * \param point pointer toward the point we want to communicate
 * \param mpi_point new MPI data type
 *
 * Create an MPI data type named mpi_point to allow the communication of
 * choosen fields of a point (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */
void PMMG_create_MPI_Point(MMG5_pPoint point, MPI_Datatype *mpi_point)
{
  int          i,blck_lengths[5] = {3, 3, 1, 1, 1};
  MPI_Aint     displs[5],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[5] = {MPI_DOUBLE,MPI_DOUBLE,MPI_INT,MPI_INT,MPI_INT16_T};

  MPI_Get_address(&(point[0]),       &lb);
  MPI_Get_address(&(point[0].c[0]),  &displs[0]);
  MPI_Get_address(&(point[0].n),     &displs[1]);
  MPI_Get_address(&(point[0].ref),   &displs[2]);
  MPI_Get_address(&(point[0].xp),    &displs[3]);
  MPI_Get_address(&(point[0].tag),   &displs[4]);
  MPI_Get_address(&(point[1]),       &ub);

  /* Relative displacement from field 0 to field i */
  for ( i=4 ; i>= 0; --i )
    displs[i] -= lb;

  MPI_Type_create_struct(5, blck_lengths, displs, types, &mpi_noextent);

  MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_point);

  MPI_Type_commit(mpi_point);
}

/**
 * \param xPoint pointer toward the xPoint we want to communicate
 * \param mpi_xpoint new MPI data type
 *
 * Create an MPI data type named mpi_xpoint to allow the communication of
 * choosen fields of a xPoint (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */
void PMMG_create_MPI_xPoint(MMG5_pxPoint xPoint, MPI_Datatype *mpi_xpoint)
{
  int          i,blck_lengths[2] = {3,3};
  MPI_Aint     displs[2],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[2] = {MPI_DOUBLE,MPI_DOUBLE};

  MPI_Get_address(&(xPoint[0]),       &lb);
  MPI_Get_address(&(xPoint[0].n1[0]), &displs[0]);
  MPI_Get_address(&(xPoint[0].n2[0]), &displs[1]);
  MPI_Get_address(&(xPoint[1]),       &ub);

  /* Relative displacement from field 0 to field i */
  for ( i=1 ; i>= 0; --i )
    displs[i] -= lb;

  MPI_Type_create_struct(2, blck_lengths, displs, types, &mpi_noextent);

  MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_xpoint);

  MPI_Type_commit(mpi_xpoint);
}

/**
 * \param tetra pointer toward the tetra we want to communicate
 * \param mpi_light_tetra new MPI data type
 *
 * Create an MPI data type to allow the communication of
 * choosen fields of a tetra (only those who needs to be communicated) when the
 * mesh has not be analyzed (no need to communicate the xt field).
 *
 * \warning to fill when we need to communicate additional things
 */
void PMMG_create_MPI_lightTetra(MMG5_pTetra tetra, MPI_Datatype *mpi_light_tetra)
{
  int          i,blck_lengths[3] = {4, 1, 1};
  MPI_Aint     displs[3],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[3] = {MPI_INT,MPI_INT,MPI_INT16_T};

  MPI_Get_address(&(tetra[0]),      &lb);
  MPI_Get_address(&(tetra[0].v[0]), &displs[0]);
  MPI_Get_address(&(tetra[0].ref),  &displs[1]);
  MPI_Get_address(&(tetra[0].tag),  &displs[2]);
  MPI_Get_address(&(tetra[1]),      &ub);

  /* Relative displacement from field 0 to field i */
  for ( i=2 ; i>= 0; --i )
    displs[i] -= lb;

  MPI_Type_create_struct(3, blck_lengths, displs, types, &mpi_noextent);

  MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_light_tetra);

  MPI_Type_commit(mpi_light_tetra);
}

/**
 * \param tetra pointer toward the tetra we want to communicate
 * \param mpi_tetra new MPI data type
 *
 * Create an MPI data type to allow the communication of
 * choosen fields of a tetra (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */
void PMMG_create_MPI_Tetra(MMG5_pTetra tetra, MPI_Datatype *mpi_tetra)
{
  int          i,blck_lengths[4] = {4, 1, 1, 1};
  MPI_Aint     displs[4],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[4] = {MPI_INT,MPI_INT,MPI_INT,MPI_INT16_T};

  MPI_Get_address(&(tetra[0]),      &lb);
  MPI_Get_address(&(tetra[0].v[0]), &displs[0]);
  MPI_Get_address(&(tetra[0].ref),  &displs[1]);
  MPI_Get_address(&(tetra[0].xt),   &displs[2]);
  MPI_Get_address(&(tetra[0].tag),  &displs[3]);
  MPI_Get_address(&(tetra[1]),      &ub);

  /* Relative displacement from field 0 to field i */
  for ( i=3 ; i>= 0; --i )
    displs[i] -= lb;

  MPI_Type_create_struct(4, blck_lengths, displs, types, &mpi_noextent);

  MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_tetra);

  MPI_Type_commit(mpi_tetra);
}

/**
 * \param edge pointer toward the edge we want to communicate
 * \param mpi_edge new MPI data type
 *
 * Create an MPI data type to allow the communication of
 * choosen fields of a edge (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */
void PMMG_create_MPI_Edge(MMG5_pEdge edge, MPI_Datatype *mpi_edge)
{
  int          i,blck_lengths[4] = {1, 1, 1, 1};
  MPI_Aint     displs[4],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[4] = {MPI_INT,MPI_INT,MPI_INT,MPI_INT16_T};

  MPI_Get_address(&(edge[0]),      &lb);
  MPI_Get_address(&(edge[0].a),    &displs[0]);
  MPI_Get_address(&(edge[0].b),    &displs[1]);
  MPI_Get_address(&(edge[0].ref),  &displs[2]);
  MPI_Get_address(&(edge[0].tag),  &displs[3]);
  MPI_Get_address(&(edge[1]),      &ub);

  /* Relative displacement from field 0 to field i */
  for ( i=3 ; i>= 0; --i ) {
    displs[i] -= lb;
  }

  MPI_Type_create_struct(4, blck_lengths, displs, types, &mpi_noextent);

  MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_edge);

  MPI_Type_commit(mpi_edge);
}

/**
 * \param tria pointer toward the tria we want to communicate
 * \param mpi_tria new MPI data type
 *
 * Create an MPI data type to allow the communication of
 * choosen fields of a tria (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */
void PMMG_create_MPI_Tria(MMG5_pTria tria, MPI_Datatype *mpi_tria)
{
  int          i,blck_lengths[3] = {3, 1, 3};
  MPI_Aint     displs[3],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[3] = {MPI_INT,MPI_INT,MPI_INT16_T};

  MPI_Get_address(&(tria[0]),       &lb);
  MPI_Get_address(&(tria[0].v[0]),  &displs[0]);
  MPI_Get_address(&(tria[0].ref),   &displs[1]);
  MPI_Get_address(&(tria[0].tag[0]),&displs[2]);
  MPI_Get_address(&(tria[1]),       &ub);

  /* Relative displacement from field 0 to field i */
  for ( i=2 ; i>= 0; --i )
    displs[i] -= lb;

  MPI_Type_create_struct(3, blck_lengths, displs, types, &mpi_noextent);

  MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_tria);

  MPI_Type_commit(mpi_tria);
}

/**
 * \param xTetra pointer toward the xTetra we want to communicate
 * \param mpi_xtetra new MPI data type
 *
 * Create an MPI data type to allow the communication of
 * choosen fields of a xTetra (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */

void PMMG_create_MPI_xTetra(MMG5_pxTetra xTetra, MPI_Datatype *mpi_xtetra)
{
  MPI_Aint     displs[4],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[4] = {MPI_INT,MPI_INT,MPI_INT16_T,MPI_INT16_T};
  int          i,blck_lengths[4] = {4,6,4,6};

  MPI_Get_address(&(xTetra[0])       ,  &lb);
  MPI_Get_address(&(xTetra[0].ref[0]),  &displs[0]);
  MPI_Get_address(&(xTetra[0].edg[0]),  &displs[1]);
  MPI_Get_address(&(xTetra[0].ftag[0]), &displs[2]);
  MPI_Get_address(&(xTetra[0].tag[0]),  &displs[3]);
  MPI_Get_address(&(xTetra[1])       ,  &ub);

 /* Relative displacement from field 0 to field i */
  for ( i=3 ; i>= 0; --i )
    displs[i] -= lb;

  MPI_Type_create_struct(4, blck_lengths, displs, types, &mpi_noextent);

  MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_xtetra);

  MPI_Type_commit(mpi_xtetra);
}

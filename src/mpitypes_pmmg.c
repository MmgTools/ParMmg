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
 * \param mpi_light_point new MPI data type
 *
 * Create an MPI data type to allow the communication of
 * choosen fields of a point (only those who needs to be communicated) when the
 * mesh hast not yet been analyzed (xp field must not be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */
int PMMG_create_MPI_lightPoint(MPI_Datatype *mpi_light_point)
{
  MMG5_Point   point[2];
  int          i,blck_lengths[4] = {3, 3, 1, 1};
  MPI_Aint     displs[4],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[4] = {MPI_DOUBLE,MPI_DOUBLE,MPI_INT,MPI_INT16_T};

  MPI_CHECK( MPI_Get_address(&(point[0]),       &lb),return 0 );
  MPI_CHECK( MPI_Get_address(&(point[0].c[0]),  &displs[0]),return 0 );
  MPI_CHECK( MPI_Get_address(&(point[0].n),     &displs[1]),return 0 );
  MPI_CHECK( MPI_Get_address(&(point[0].ref),   &displs[2]),return 0 );
  MPI_CHECK( MPI_Get_address(&(point[0].tag),   &displs[3]),return 0 );
  MPI_CHECK( MPI_Get_address(&(point[1]),       &ub),return 0 );

  /* Relative displacement from field 0 to field i */
  for ( i=3 ; i>= 0; --i )
    displs[i] -= lb;

  MPI_CHECK( MPI_Type_create_struct(4, blck_lengths, displs, types, &mpi_noextent),
             return 0);

  MPI_CHECK( MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_light_point),
             return 0);

  MPI_CHECK( MPI_Type_commit(mpi_light_point),return 0);

  return 1;
}

/**
 * \param mpi_point new MPI data type
 *
 * Create an MPI data type named mpi_point to allow the communication of
 * choosen fields of a point (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */
int PMMG_create_MPI_Point(MPI_Datatype *mpi_point)
{
  MMG5_Point   point[2];
  int          i,blck_lengths[6] = {3, 3, 1, 1, 1, 1};
  MPI_Aint     displs[6],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[6] = {MPI_DOUBLE,MPI_DOUBLE,MPI_INT,MPI_INT,MPI_INT,MPI_INT16_T};

  MPI_CHECK( MPI_Get_address(&(point[0]),       &lb),return 0 );
  MPI_CHECK( MPI_Get_address(&(point[0].c[0]),  &displs[0]),return 0 );
  MPI_CHECK( MPI_Get_address(&(point[0].n),     &displs[1]),return 0 );
  MPI_CHECK( MPI_Get_address(&(point[0].src),   &displs[2]),return 0 );
  MPI_CHECK( MPI_Get_address(&(point[0].ref),   &displs[3]),return 0 );
  MPI_CHECK( MPI_Get_address(&(point[0].xp),    &displs[4]),return 0 );
  MPI_CHECK( MPI_Get_address(&(point[0].tag),   &displs[5]),return 0 );
  MPI_CHECK( MPI_Get_address(&(point[1]),       &ub),return 0 );

  /* Relative displacement from field 0 to field i */
  for ( i=5 ; i>= 0; --i )
    displs[i] -= lb;

  MPI_CHECK( MPI_Type_create_struct(6, blck_lengths, displs, types, &mpi_noextent),
             return 0 );

  MPI_CHECK( MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_point),return 0);

  MPI_CHECK( MPI_Type_commit(mpi_point),return 0);

  return 1;
}

/**
 * \param mpi_xpoint new MPI data type
 *
 * Create an MPI data type named mpi_xpoint to allow the communication of
 * choosen fields of a xPoint (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */
int PMMG_create_MPI_xPoint(MPI_Datatype *mpi_xpoint)
{
  MMG5_xPoint  xPoint[2];
  int          i,blck_lengths[2] = {3,3};
  MPI_Aint     displs[2],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[2] = {MPI_DOUBLE,MPI_DOUBLE};

  MPI_CHECK( MPI_Get_address(&(xPoint[0]),       &lb),return 0);
  MPI_CHECK( MPI_Get_address(&(xPoint[0].n1[0]), &displs[0]),return 0);
  MPI_CHECK( MPI_Get_address(&(xPoint[0].n2[0]), &displs[1]),return 0);
  MPI_CHECK( MPI_Get_address(&(xPoint[1]),       &ub),return 0);

  /* Relative displacement from field 0 to field i */
  for ( i=1 ; i>= 0; --i )
    displs[i] -= lb;

  MPI_CHECK( MPI_Type_create_struct(2, blck_lengths, displs, types, &mpi_noextent),
             return 0 );

  MPI_CHECK( MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_xpoint),return 0);

  MPI_CHECK( MPI_Type_commit(mpi_xpoint),return 0);

  return 1;
}

/**
 * \param mpi_light_tetra new MPI data type
 *
 * Create an MPI data type to allow the communication of
 * choosen fields of a tetra (only those who needs to be communicated) when the
 * mesh has not be analyzed (no need to communicate the xt field).
 *
 * \warning to fill when we need to communicate additional things
 */
int PMMG_create_MPI_lightTetra(MPI_Datatype *mpi_light_tetra)
{
  MMG5_Tetra   tetra[2];
  int          i,blck_lengths[3] = {4, 1, 1};
  MPI_Aint     displs[3],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[3] = {MPI_INT,MPI_INT,MPI_INT16_T};

   MPI_CHECK( MPI_Get_address(&(tetra[0]),      &lb),return 0);
   MPI_CHECK( MPI_Get_address(&(tetra[0].v[0]), &displs[0]),return 0);
   MPI_CHECK( MPI_Get_address(&(tetra[0].ref),  &displs[1]),return 0);
   MPI_CHECK( MPI_Get_address(&(tetra[0].tag),  &displs[2]),return 0);
   MPI_CHECK( MPI_Get_address(&(tetra[1]),      &ub),return 0);

  /* Relative displacement from field 0 to field i */
  for ( i=2 ; i>= 0; --i )
    displs[i] -= lb;

  MPI_CHECK( MPI_Type_create_struct(3, blck_lengths, displs, types, &mpi_noextent),
             return 0);

  MPI_CHECK( MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_light_tetra),
             return 0);

  MPI_CHECK( MPI_Type_commit(mpi_light_tetra),return 0);;

  return 1;
}

/**
 * \param mpi_tetra new MPI data type
 *
 * Create an MPI data type to allow the communication of
 * choosen fields of a tetra (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */
int PMMG_create_MPI_Tetra(MPI_Datatype *mpi_tetra)
{
  MMG5_Tetra   tetra[2];
  int          i,blck_lengths[6] = {4, 1, 1, 1, 1, 1};
  MPI_Aint     displs[6],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[6] = {MPI_INT,MPI_INT,MPI_INT,MPI_INT16_T,MPI_INT,MPI_DOUBLE};

   MPI_CHECK( MPI_Get_address(&(tetra[0]),      &lb),return 0);
   MPI_CHECK( MPI_Get_address(&(tetra[0].v[0]), &displs[0]),return 0);
   MPI_CHECK( MPI_Get_address(&(tetra[0].ref),  &displs[1]),return 0);
   MPI_CHECK( MPI_Get_address(&(tetra[0].xt),   &displs[2]),return 0);
   MPI_CHECK( MPI_Get_address(&(tetra[0].tag),  &displs[3]),return 0);
   MPI_CHECK( MPI_Get_address(&(tetra[0].mark), &displs[4]),return 0);
   MPI_CHECK( MPI_Get_address(&(tetra[0].qual), &displs[5]),return 0);
   MPI_CHECK( MPI_Get_address(&(tetra[1]),      &ub),return 0);

  /* Relative displacement from field 0 to field i */
  for ( i=5 ; i>= 0; --i )
    displs[i] -= lb;

  MPI_CHECK( MPI_Type_create_struct(6, blck_lengths, displs, types, &mpi_noextent),
             return 0);

  MPI_CHECK( MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_tetra),return 0);

  MPI_CHECK( MPI_Type_commit(mpi_tetra),return 0);

  return 1;
}

/**
 * \param mpi_edge new MPI data type
 *
 * Create an MPI data type to allow the communication of
 * choosen fields of a edge (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */
int PMMG_create_MPI_Edge(MPI_Datatype *mpi_edge)
{
  MMG5_Edge    edge[2];
  int          i,blck_lengths[4] = {1, 1, 1, 1};
  MPI_Aint     displs[4],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[4] = {MPI_INT,MPI_INT,MPI_INT,MPI_INT16_T};

   MPI_CHECK( MPI_Get_address(&(edge[0]),      &lb),return 0);
   MPI_CHECK( MPI_Get_address(&(edge[0].a),    &displs[0]),return 0);
   MPI_CHECK( MPI_Get_address(&(edge[0].b),    &displs[1]),return 0);
   MPI_CHECK( MPI_Get_address(&(edge[0].ref),  &displs[2]),return 0);
   MPI_CHECK( MPI_Get_address(&(edge[0].tag),  &displs[3]),return 0);
   MPI_CHECK( MPI_Get_address(&(edge[1]),      &ub),return 0);

  /* Relative displacement from field 0 to field i */
  for ( i=3 ; i>= 0; --i ) {
    displs[i] -= lb;
  }

  MPI_CHECK( MPI_Type_create_struct(4, blck_lengths, displs, types, &mpi_noextent),
             return 0);

  MPI_CHECK( MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_edge),return 0);

  MPI_CHECK( MPI_Type_commit(mpi_edge),return 0);

  return 1;
}

/**
 * \param mpi_tria new MPI data type
 *
 * Create an MPI data type to allow the communication of
 * choosen fields of a tria (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */
int PMMG_create_MPI_Tria(MPI_Datatype *mpi_tria)
{
  MMG5_Tria    tria[2];
  int          i,blck_lengths[3] = {3, 1, 3};
  MPI_Aint     displs[3],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[3] = {MPI_INT,MPI_INT,MPI_INT16_T};

   MPI_CHECK( MPI_Get_address(&(tria[0]),       &lb),return 0);
   MPI_CHECK( MPI_Get_address(&(tria[0].v[0]),  &displs[0]),return 0);
   MPI_CHECK( MPI_Get_address(&(tria[0].ref),   &displs[1]),return 0);
   MPI_CHECK( MPI_Get_address(&(tria[0].tag[0]),&displs[2]),return 0);
   MPI_CHECK( MPI_Get_address(&(tria[1]),       &ub),return 0);

  /* Relative displacement from field 0 to field i */
  for ( i=2 ; i>= 0; --i )
    displs[i] -= lb;

  MPI_CHECK( MPI_Type_create_struct(3, blck_lengths, displs, types, &mpi_noextent),
             return 0);

  MPI_CHECK( MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_tria),return 0);

  MPI_CHECK( MPI_Type_commit(mpi_tria),return 0);

  return 1;
}

/**
 * \param mpi_xtetra new MPI data type
 *
 * Create an MPI data type to allow the communication of
 * choosen fields of a xTetra (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */
int PMMG_create_MPI_xTetra(MPI_Datatype *mpi_xtetra)
{
  MMG5_xTetra  xTetra[2];
  MPI_Aint     displs[4],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[4] = {MPI_INT,MPI_INT,MPI_INT16_T,MPI_INT16_T};
  int          i,blck_lengths[4] = {4,6,4,6};

  MPI_CHECK( MPI_Get_address(&(xTetra[0])       ,  &lb),return 0);
  MPI_CHECK( MPI_Get_address(&(xTetra[0].ref[0]),  &displs[0]),return 0);
  MPI_CHECK( MPI_Get_address(&(xTetra[0].edg[0]),  &displs[1]),return 0);
  MPI_CHECK( MPI_Get_address(&(xTetra[0].ftag[0]), &displs[2]),return 0);
  MPI_CHECK( MPI_Get_address(&(xTetra[0].tag[0]),  &displs[3]),return 0);
  MPI_CHECK( MPI_Get_address(&(xTetra[1])       ,  &ub),return 0);

 /* Relative displacement from field 0 to field i */
  for ( i=3 ; i>= 0; --i )
    displs[i] -= lb;

  MPI_CHECK( MPI_Type_create_struct(4, blck_lengths, displs, types, &mpi_noextent),
             return 0);

  MPI_CHECK( MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_xtetra),return 0);

  MPI_CHECK( MPI_Type_commit(mpi_xtetra),return 0);

  return 1;
}

/**
 * \param mpi_point  pointer toward an MPI_Datatype
 * \param mpi_xpoint pointer toward an MPI_Datatype
 * \param mpi_tetra  pointer toward an MPI_Datatype
 * \param mpi_xtetra pointer toward an MPI_Datatype
 *
 * If used, free the \a mpi_point, \a mpi_xpoint, ... MPI_Datatype
 *
 */
int PMMG_Free_MPI_meshDatatype( MPI_Datatype *mpi_point,
                                MPI_Datatype *mpi_xpoint,
                                MPI_Datatype *mpi_tetra,
                                MPI_Datatype *mpi_xtetra ) {

  if ( *mpi_xtetra ) {
    MPI_Type_free( mpi_xtetra );
    *mpi_xtetra = MPI_DATATYPE_NULL;
  }

  if ( *mpi_tetra ) {
    MPI_Type_free( mpi_tetra );
    *mpi_tetra = MPI_DATATYPE_NULL;
  }

  if ( *mpi_xpoint ) {
    MPI_Type_free( mpi_xpoint );
    *mpi_xpoint = MPI_DATATYPE_NULL;
  }

  if ( *mpi_point ) {
    MPI_Type_free( mpi_point );
    *mpi_point = MPI_DATATYPE_NULL;
  }

  return 1;
}

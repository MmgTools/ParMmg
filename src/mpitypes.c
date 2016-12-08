/**
 * \file mpitypes.c
 * \brief Mpi types management.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 */

#include "libparmmg.h"

/**
 * \param parmesh pointer toward the parmesh we want to communicate
 * \param mpi_type_parmesh new MPI data type for parmesh
 * \param mpi_type_int_comm new MPI data type for internal communicator
 * \param mpi_type_ext_comm new MPI data type for external communicator
 * \param mpi_type_grp new MPI data type for groupe
 * \param mpi_type_mesh new MPI data type for mesh
 * \param mpi_type_point new MPI data type for point
 * \param mpi_type_xpoint new MPI data type for xpoint
 * \param mpi_type_tetra new MPI data type for tetra
 * \param mpi_type_xtetra new MPI data type for xtetra
 * \return 0 if fail, 1 otherwise.
 *
 * Create an MPI data type named mpi_type_parmesh to allow the comminication of
 * choosen fields of a parmesh (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 *
 * \remark we suppose that we have merged groups over each proc, thus we only
 * have 1 group per proc.
 */
int PMMG_create_MPI_ParMesh(PMMG_pParMesh parmesh, MPI_Datatype *mpi_type_parmesh,
                            MPI_Datatype *mpi_type_int_comm,MPI_Datatype *mpi_type_ext_comm,
                            MPI_Datatype *mpi_type_grp, MPI_Datatype *mpi_type_mesh,
                            MPI_Datatype *mpi_type_point,MPI_Datatype *mpi_type_xpoint,
                            MPI_Datatype *mpi_type_tetra,MPI_Datatype *mpi_type_xtetra) {

  int          blck_lengths[4],i;
  MPI_Datatype types[4];
  MPI_Aint     displs[4];

  if ( !PMMG_create_MPI_Int_Comm(parmesh->int_node_comm,mpi_type_int_comm) )
    return(0);

  if ( !PMMG_create_MPI_Ext_Comm(&parmesh->ext_node_comm[0],mpi_type_ext_comm) )
    return(0);

  if ( !PMMG_create_MPI_Grp(&parmesh->listgrp[0],mpi_type_grp,mpi_type_mesh,
                            mpi_type_point,mpi_type_xpoint,
                            mpi_type_tetra,mpi_type_xtetra) )  return(0);

  blck_lengths[0] = 1; // grp
  blck_lengths[1] = 1; // int_comm
  blck_lengths[2] = 1; // next_node_comm
  blck_lengths[3] = parmesh->next_node_comm; // ext_node_comm

  types[0] = *mpi_type_grp;
  types[1] = *mpi_type_int_comm;
  types[2] = MPI_INT;
  types[3] = *mpi_type_ext_comm;

  MPI_Get_address(&parmesh->listgrp[0],       &displs[0]);
  MPI_Get_address(&parmesh->int_node_comm,    &displs[1]);
  MPI_Get_address(&parmesh->next_node_comm,   &displs[2]);
  MPI_Get_address(&parmesh->ext_node_comm,    &displs[3]);

  /* Relative displacement from field 0 to field i */
  for ( i=3 ; i>= 0; --i ) {
    displs[i] -= displs[0];
  }

  MPI_Type_create_struct(4, blck_lengths,displs,types,mpi_type_parmesh);

  MPI_Type_commit(mpi_type_parmesh);

  return(1);
}

/**
 * \param mpi_type_parmesh MPI data type to delete
 * \param mpi_type_int_comm MPI data type to delete
 * \param mpi_type_ext_comm MPI data type to delete
 * \param mpi_type_grp MPI data type to delete
 * \param mpi_type_mesh MPI data type to delete
 * \param mpi_type_point MPI data type to delete
 * \param mpi_type_xpoint MPI data type to delete
 * \param mpi_type_tetra MPI data type to delete
 * \param mpi_type_xtetra MPI data type to delete
 *
 * Free the MPI data type named mpi_type_parmesh
 *
 */
void PMMG_free_MPI_ParMesh(MPI_Datatype *mpi_type_parmesh,
                           MPI_Datatype *mpi_type_int_comm,MPI_Datatype *mpi_type_ext_comm,
                           MPI_Datatype *mpi_type_grp, MPI_Datatype *mpi_type_mesh,
                           MPI_Datatype *mpi_type_point,MPI_Datatype *mpi_type_xpoint,
                           MPI_Datatype *mpi_type_tetra,MPI_Datatype *mpi_type_xtetra ) {

   MPI_Type_free(mpi_type_parmesh);

   PMMG_free_MPI_Int_Comm( mpi_type_int_comm );

   PMMG_free_MPI_Ext_Comm( mpi_type_ext_comm );

   PMMG_free_MPI_Grp(mpi_type_grp,mpi_type_mesh,mpi_type_point,mpi_type_xpoint,
                     mpi_type_tetra,mpi_type_xtetra);

}

/**
 * \param int_comm pointer toward the PMMG_int_comm we want to communicate
 * \param mpi_type_int_comm new MPI data type
 * \return 0 if fail, 1 otherwise.
 *
 * Create an MPI data type named mpi_type_int_comm to allow the comminication of
 * choosen fields of a int_comm (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */

int PMMG_create_MPI_Int_Comm(PMMG_pint_comm int_comm, MPI_Datatype *mpi_type_int_comm) {
  int          blck_lengths[2],i;
  MPI_Datatype types[2];
  MPI_Aint     displs[2];

  blck_lengths[0] = 1; // nitem
  blck_lengths[1] = int_comm->nitem; // intvalues

  types[0] = MPI_INT;
  types[1] = MPI_INT;

  MPI_Get_address(&int_comm->nitem,       &displs[0]);
  MPI_Get_address(&int_comm->intvalues,   &displs[1]);

  /* Relative displacement from field 0 to field i */
  for ( i=1 ; i>= 0; --i ) {
    displs[i] -= displs[0];
  }

  MPI_Type_create_struct(2, blck_lengths, displs, types, mpi_type_int_comm);

  MPI_Type_commit(mpi_type_int_comm);

  return(1);
}

/**
 * \param mpi_type_int_comm new MPI data type
 *
 * Free the MPI data type named mpi_type_int_comm
 *
 */

void PMMG_free_MPI_Int_Comm(MPI_Datatype *mpi_type_int_comm) {

  MPI_Type_free ( mpi_type_int_comm );

}

/**
 * \param ext_comm pointer toward the PMMG_ext_comm we want to communicate
 * \param mpi_type_ext_comm new MPI data type
 * \return 0 if fail, 1 otherwise.
 *
 * Create an MPI data type named mpi_type_ext_comm to allow the comminication of
 * choosen fields of a ext_comm (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */

int PMMG_create_MPI_Ext_Comm(PMMG_pext_comm ext_comm, MPI_Datatype *mpi_type_ext_comm) {
  int          blck_lengths[4],i;
  MPI_Datatype types[4];
  MPI_Aint     displs[4];

  blck_lengths[0] = 1; // color_in
  blck_lengths[1] = 1; // color_out
  blck_lengths[2] = 1; // nitem
  blck_lengths[3] = ext_comm->nitem; // intvalues

  types[0] = MPI_INT;
  types[1] = MPI_INT;
  types[2] = MPI_INT;
  types[3] = MPI_INT;

  MPI_Get_address(&ext_comm->color_in,       &displs[0]);
  MPI_Get_address(&ext_comm->color_out,      &displs[1]);
  MPI_Get_address(&ext_comm->nitem,          &displs[2]);
  MPI_Get_address(&ext_comm->int_comm_index, &displs[3]);

  /* Relative displacement from field 0 to field i */
  for ( i=3 ; i>= 0; --i ) {
    displs[i] -= displs[0];
  }

  MPI_Type_create_struct(4, blck_lengths, displs, types, mpi_type_ext_comm);

  MPI_Type_commit(mpi_type_ext_comm);

  return(1);
}

/**
 * \param mpi_type_ext_comm new MPI data type
 *
 * Free the MPI data type named mpi_type_ext_comm
 *
 */

void PMMG_free_MPI_Ext_Comm(MPI_Datatype *mpi_type_ext_comm) {

  MPI_Type_free ( mpi_type_ext_comm );

}


/**
 * \param grp pointer toward the group we want to communicate
 * \param mpi_type_grp new MPI data type for groupe
 * \param mpi_type_mesh new MPI data type for mesh
 * \param mpi_type_point new MPI data type for point
 * \param mpi_type_xpoint new MPI data type for xpoint
 * \param mpi_type_tetra new MPI data type for tetra
 * \param mpi_type_xtetra new MPI data type for xtetra
 * \return 0 if fail, 1 otherwise.
 *
 * Create an MPI data type named mpi_type_grp to allow the comminication of
 * choosen fields of a grp (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */

int PMMG_create_MPI_Grp(PMMG_pGrp grp, MPI_Datatype *mpi_type_grp,
                        MPI_Datatype *mpi_type_mesh,
                        MPI_Datatype *mpi_type_point,MPI_Datatype *mpi_type_xpoint,
                        MPI_Datatype *mpi_type_tetra,MPI_Datatype *mpi_type_xtetra) {
  int          blck_lengths[4],i;
  MPI_Datatype types[4];
  MPI_Aint     displs[4];

  if( !PMMG_create_MPI_Mesh( grp->mesh,mpi_type_mesh,
                             mpi_type_point,mpi_type_xpoint,
                             mpi_type_tetra,mpi_type_xtetra ) )   return(0);

  blck_lengths[0] = 1; // mesh
  blck_lengths[1] = 1; // nitem_int_node_comm
  blck_lengths[2] = grp->nitem_int_node_comm; // node2int_node_comm_index1
  blck_lengths[3] = grp->nitem_int_node_comm; // node2int_node_comm_index2

  types[0] = *mpi_type_mesh;
  types[1] = types[2] = types[3] = MPI_INT;

  MPI_Get_address( grp->mesh,                      &displs[0]);
  MPI_Get_address(&grp->nitem_int_node_comm,       &displs[1]);
  MPI_Get_address(&grp->node2int_node_comm_index1, &displs[2]);
  MPI_Get_address(&grp->node2int_node_comm_index2, &displs[3]);

  /* Relative displacement from field 0 to field i */
  for ( i=3 ; i>= 0; --i ) {
    displs[i] -= displs[0];
  }

  MPI_Type_create_struct(4, blck_lengths, displs, types, mpi_type_grp);

  MPI_Type_commit(mpi_type_grp);

  return(1);
}

/**
 * \param mpi_type_grp grp MPI data type
 * \param mpi_type_mesh mesh MPI data type for mesh
 * \param mpi_type_point point MPI data type for point
 * \param mpi_type_xpoint xpoint MPI data type for xpoint
 * \param mpi_type_tetra tetra MPI data type for tetra
 * \param mpi_type_xtetra xtetra MPI data type for xtetra
 *
 * Free the MPI data type named mpi_type_grp
 *
 */

void PMMG_free_MPI_Grp(MPI_Datatype *mpi_type_grp, MPI_Datatype *mpi_type_mesh,
                       MPI_Datatype *mpi_type_point,MPI_Datatype *mpi_type_xpoint,
                       MPI_Datatype *mpi_type_tetra,MPI_Datatype *mpi_type_xtetra) {

  PMMG_free_MPI_Mesh ( mpi_type_mesh,mpi_type_point,mpi_type_xpoint,
                       mpi_type_tetra,mpi_type_xtetra);

  MPI_Type_free ( mpi_type_grp );

}


/**
 * \param mesh pointer toward the mesh we want to communicate
 * \param mpi_type_mesh new MPI data type for mesh
 * \param mpi_type_point new MPI data type for point
 * \param mpi_type_xpoint new MPI data type for xpoint
 * \param mpi_type_tetra new MPI data type for tetra
 * \param mpi_type_xtetra new MPI data type for xtetra
 * \return 0 if fail, 1 otherwise.
 *
 * Create an MPI data type named mpi_type_mesh to allow the comminication of
 * choosen fields of a mesh (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */

int PMMG_create_MPI_Mesh(MMG5_pMesh mesh, MPI_Datatype *mpi_type_mesh,
                         MPI_Datatype *mpi_type_point,MPI_Datatype *mpi_type_xpoint,
                         MPI_Datatype *mpi_type_tetra,MPI_Datatype *mpi_type_xtetra) {
  int          blck_lengths[8],i;
  MPI_Datatype types[8];
  MPI_Aint     displs[8];

  if( !PMMG_create_MPI_Point( &mesh->point[1] , mpi_type_point) )   return(0);
  if( !PMMG_create_MPI_xPoint(&mesh->xpoint[1], mpi_type_xpoint) )  return(0);

  if( !PMMG_create_MPI_Tetra( &mesh->tetra[1],  mpi_type_tetra) )  return(0);
  if( !PMMG_create_MPI_xTetra(&mesh->xtetra[1], mpi_type_xtetra) )  return(0);

  blck_lengths[0] = 1; // np
  blck_lengths[1] = 1; // ne
  blck_lengths[2] = 1; // xp
  blck_lengths[3] = 1; // xt
  blck_lengths[4] = mesh->np; // point
  blck_lengths[5] = mesh->xp; // xpoint
  blck_lengths[6] = mesh->ne; // tetra
  blck_lengths[7] = mesh->xt; // xtetra

  types[0] = MPI_INT;
  types[1] = MPI_INT;
  types[2] = MPI_INT;
  types[3] = MPI_INT;
  types[4] = *mpi_type_point;
  types[5] = *mpi_type_xpoint;
  types[6] = *mpi_type_tetra;
  types[7] = *mpi_type_xtetra;

  MPI_Get_address(&mesh->np,     &displs[0]);
  MPI_Get_address(&mesh->ne,     &displs[1]);
  MPI_Get_address(&mesh->xp,     &displs[2]);
  MPI_Get_address(&mesh->xt,     &displs[3]);
  MPI_Get_address(&mesh->point,  &displs[4]);
  MPI_Get_address(&mesh->xpoint, &displs[5]);
  MPI_Get_address(&mesh->tetra,  &displs[6]);
  MPI_Get_address(&mesh->xtetra, &displs[7]);

  /* Relative displacement from field 0 to field i */
  for ( i=7 ; i>= 0; --i ) {
    displs[i] -= displs[0];
  }

  MPI_Type_create_struct(8, blck_lengths, displs, types, mpi_type_mesh);

  MPI_Type_commit(mpi_type_mesh);

  return(1);
}

/**
 * \param mpi_type_mesh mesh MPI data type
 * \param mpi_type_point point MPI data type
 * \param mpi_type_xpoint xPoint MPI data type
 * \param mpi_type_tetra tetra MPI data type
 * \param mpi_type_xtetra xTetra MPI data type
 *
 * Free the MPI data type named mpi_type_mesh
 *
 */

void PMMG_free_MPI_Mesh(MPI_Datatype *mpi_type_mesh,
                        MPI_Datatype *mpi_type_point,MPI_Datatype *mpi_type_xpoint,
                        MPI_Datatype *mpi_type_tetra,MPI_Datatype *mpi_type_xtetra) {

  PMMG_free_MPI_Point( mpi_type_point);
  PMMG_free_MPI_xPoint( mpi_type_xpoint);

  PMMG_free_MPI_Tetra( mpi_type_tetra);
  PMMG_free_MPI_xTetra( mpi_type_xtetra);

  MPI_Type_free ( mpi_type_mesh );

}


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
  int          i,blck_lengths[4] = {3, 3, 1, 1};
  MPI_Aint     displs[4];
  MPI_Datatype types[4] = {MPI_DOUBLE,MPI_DOUBLE,MPI_INT,MPI_INT};


  MPI_Get_address(&point->c[0], &displs[0]);
  MPI_Get_address(&point->n,    &displs[1]);
  MPI_Get_address(&point->ref,  &displs[2]);
  MPI_Get_address(&point->xp,   &displs[3]);

  /* Relative displacement from field 0 to field i */
  for ( i=3 ; i>= 0; --i ) {
    displs[i] -= displs[0];
  }

  MPI_Type_create_struct(4, blck_lengths, displs, types, mpi_type_point);

  MPI_Type_commit(mpi_type_point);

  return(1);
}

/**
 * \param mpi_type_point point MPI data type
 *
 * Free the MPI data type named mpi_type_point
 *
 */

void PMMG_free_MPI_Point(MPI_Datatype *mpi_type_point) {

   MPI_Type_free ( mpi_type_point );

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
  MPI_Aint     displs[2];
  MPI_Datatype types[2] = {MPI_DOUBLE,MPI_DOUBLE};


  MPI_Get_address(&xPoint->n1[0], &displs[0]);
  MPI_Get_address(&xPoint->n2[0], &displs[1]);

  /* Relative displacement from field 0 to field i */
  for ( i=1 ; i>= 0; --i ) {
    displs[i] -= displs[0];
  }

  MPI_Type_create_struct(2, blck_lengths, displs, types, mpi_type_xpoint);

  MPI_Type_commit(mpi_type_xpoint);

  return(1);
}

/**
 * \param mpi_type_xpoint xPoint MPI data type
 *
 * Free the MPI data type named mpi_type_xpoint
 *
 */

void PMMG_free_MPI_xPoint(MPI_Datatype *mpi_type_xpoint) {

   MPI_Type_free ( mpi_type_xpoint );

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
  MPI_Aint     displs[3];
  MPI_Datatype types[3] = {MPI_INT,MPI_INT,MPI_INT};


  MPI_Get_address(&tetra->v[0], &displs[0]);
  MPI_Get_address(&tetra->ref,  &displs[1]);
  MPI_Get_address(&tetra->xt,   &displs[2]);

  /* Relative displacement from field 0 to field i */
  for ( i=2 ; i>= 0; --i ) {
    displs[i] -= displs[0];
  }

  MPI_Type_create_struct(3, blck_lengths, displs, types, mpi_type_tetra);

  MPI_Type_commit(mpi_type_tetra);

  return(1);
}

/**
 * \param mpi_type_tetra tetra MPI data type
 *
 * Free the MPI data type named mpi_type_tetra
 *
 */

void PMMG_free_MPI_Tetra(MPI_Datatype *mpi_type_tetra) {

   MPI_Type_free ( mpi_type_tetra );

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
  MPI_Aint     displs[2];
  MPI_Datatype types[2] = {MPI_INT,MPI_INT};
  int          i,blck_lengths[2] = {4,6};

  MPI_Get_address(&xTetra->ref[0], &displs[0]);
  MPI_Get_address(&xTetra->edg[0], &displs[1]);

 /* Relative displacement from field 0 to field i */
  for ( i=1 ; i>= 0; --i ) {
    displs[i] -= displs[0];
  }

  MPI_Type_create_struct(2, blck_lengths, displs, types, mpi_type_xtetra);

  MPI_Type_commit(mpi_type_xtetra);

  return(1);
}

/**
 * \param mpi_type_xtetra xTetra MPI data type
 *
 * Free the MPI data type named mpi_type_xtetra
 *
 */

void PMMG_free_MPI_xTetra(MPI_Datatype *mpi_type_xtetra) {

   MPI_Type_free ( mpi_type_xtetra );

}

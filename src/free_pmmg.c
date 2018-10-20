/**
 * \file free_pmmg.c
 * \brief Free parmmg structures
 * \copyright GNU Lesser General Public License.
 *
 * Functions to free the parmesh structure
 *
 */

#include "parmmg.h"

/**
 * \param parmesh pointer toward a parmesh structure
 * \param idx1    node2int_node_comm_index1 to be freed
 * \param idx2    node2int_node_comm_index2 to be freed
 * \param n       pointer to node2int_node_comm_nitem size
 *
 * Deallocate all the MMG3D meshes and their communicators and zero the size
 */
static void PMMG_parmesh_grp_comm_free( PMMG_pParMesh parmesh,
                                        int **idx1, int **idx2, int *n )
{
  PMMG_DEL_MEM(parmesh,*idx1,int,"group communicator");
  PMMG_DEL_MEM(parmesh,*idx2,int,"group communicator");
  *n = 0;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param listgrp group of MMG3D meshes in parmesh
 * \param ngrp    number of mmg meshes in listgrp
 *
 * Deallocate all the MMG3D meshes and their communicators
 */
void PMMG_listgrp_free( PMMG_pParMesh parmesh, PMMG_pGrp *listgrp, int ngrp )
{
  int k;

  for ( k = 0; k < ngrp; ++k )
    PMMG_grp_free( parmesh, listgrp[0] + k );

  PMMG_DEL_MEM(parmesh,*listgrp,PMMG_Grp,"Deallocating listgrp container");
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param grp     group to free
 *
 * Deallocate all the MMG3D meshes and their communicators
 */
void PMMG_grp_free( PMMG_pParMesh parmesh, PMMG_pGrp grp )
{
  PMMG_parmesh_grp_comm_free( parmesh,
                              &grp->node2int_node_comm_index1,
                              &grp->node2int_node_comm_index2,
                              &grp->nitem_int_node_comm);
  PMMG_parmesh_grp_comm_free( parmesh,
                              &grp->edge2int_edge_comm_index1,
                              &grp->edge2int_edge_comm_index2,
                              &grp->nitem_int_edge_comm);
  PMMG_parmesh_grp_comm_free( parmesh,
                              &grp->face2int_face_comm_index1,
                              &grp->face2int_face_comm_index2,
                              &grp->nitem_int_face_comm);
  MMG3D_Free_all( MMG5_ARG_start,
                  MMG5_ARG_ppMesh, &grp->mesh,
                  MMG5_ARG_ppMet, &grp->met,
                  MMG5_ARG_end );
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * Free parmesh communicators that are allocated
 */
void PMMG_parmesh_Free_Comm( PMMG_pParMesh parmesh )
{
  PMMG_parmesh_int_comm_free( parmesh, parmesh->int_node_comm );
  PMMG_DEL_MEM(parmesh, parmesh->int_node_comm, PMMG_Int_comm, "int node comm");
  PMMG_parmesh_int_comm_free( parmesh, parmesh->int_edge_comm );
  PMMG_DEL_MEM(parmesh, parmesh->int_edge_comm, PMMG_Int_comm, "int edge comm");
  PMMG_parmesh_int_comm_free( parmesh, parmesh->int_face_comm );
  PMMG_DEL_MEM(parmesh, parmesh->int_face_comm, PMMG_Int_comm, "int face comm");
 
  PMMG_parmesh_ext_comm_free( parmesh, parmesh->ext_node_comm, parmesh->next_node_comm );
  PMMG_DEL_MEM(parmesh, parmesh->ext_node_comm,PMMG_Ext_comm, "ext node comm");
  PMMG_parmesh_ext_comm_free( parmesh, parmesh->ext_edge_comm, parmesh->next_edge_comm );
  PMMG_DEL_MEM(parmesh, parmesh->ext_edge_comm,PMMG_Ext_comm, "ext edge comm");
  PMMG_parmesh_ext_comm_free( parmesh, parmesh->ext_face_comm, parmesh->next_face_comm );
  PMMG_DEL_MEM(parmesh, parmesh->ext_face_comm,PMMG_Ext_comm, "ext face comm");
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * Free parmesh listgrp that are allocated
 */
void PMMG_parmesh_Free_Listgrp( PMMG_pParMesh parmesh )
{
  PMMG_listgrp_free( parmesh, &parmesh->listgrp, parmesh->ngrp );
}

/**
 * \param parmesh pointer toward a parmesh structure.
 * \param listgrp pointer toward an array of groups.
 * \param ngrp number of groups in the array.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Delete empty mesh and associated metrics and solutions in list of groups
 *
 */
int PMMG_clean_emptyMesh( PMMG_pParMesh parmesh, PMMG_pGrp listgrp, int ngrp ) {
  MMG5_pMesh mesh;
  int        k;

  for ( k=0; k<ngrp; ++k ) {
    mesh = listgrp[k].mesh;

    if ( !mesh ) continue;

    if ( (!mesh->ne) && (!mesh->np) ) {
      PMMG_grp_free( parmesh, &listgrp[k] );
    }
  }

  return 1;
}

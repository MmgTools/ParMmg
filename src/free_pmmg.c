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
 * \param comm    external communicator to be freed
 *
 * deallocate all internal communicator's fields
 */
static void PMMG_parmesh_int_comm_free( PMMG_pParMesh parmesh,
                                        PMMG_pInt_comm comm )
{
  if ( comm == NULL )
    return;

  if ( NULL != comm->intvalues ) {
    assert ( comm->nitem != 0 && "incorrect parameters in internal communicator" );
    PMMG_DEL_MEM(parmesh,comm->intvalues,comm->nitem,int,"int comm int array");
  }
  if ( NULL != comm->doublevalues ) {
    assert ( comm->nitem != 0 && "incorrect parameters in internal communicator" );
    PMMG_DEL_MEM(parmesh,
                 comm->doublevalues,comm->nitem,double,"int comm double array");
  }
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param comm    external communicator to be freed
 * \param ncomm   parameter ncomm
 *
 * deallocate all external communicators's fields
 */
static void PMMG_parmesh_ext_comm_free( PMMG_pParMesh parmesh,
                                        PMMG_pExt_comm comm, int ncomm )
{
  int i = 0;

  if ( comm == NULL )
    return;

  for( i = 0; i < ncomm; ++i ) {
    if ( NULL != comm->int_comm_index ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(parmesh,comm->int_comm_index,comm->nitem,int,"ext comm int array");
    }
    if ( NULL != comm->itosend ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(parmesh,comm->itosend,comm->nitem,int,"ext comm itosend array");
    }
    if ( NULL != comm->itorecv ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(parmesh,comm->itorecv,comm->nitem,int,"ext comm itorecv array");
    }
    if ( NULL != comm->rtosend ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(parmesh,comm->rtosend,comm->nitem,int,"ext comm rtosend array");
    }
    if ( NULL != comm->rtorecv ) {
      assert ( comm->nitem != 0 && "incorrect parameters in external communicator" );
      PMMG_DEL_MEM(parmesh,comm->rtorecv,comm->nitem,int,"ext comm rtorecv array");
    }
  }
}

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
  PMMG_DEL_MEM(parmesh,*idx1,*n,int,"group communicator");
  PMMG_DEL_MEM(parmesh,*idx2,*n,int,"group communicator");
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

  PMMG_DEL_MEM(parmesh,*listgrp,ngrp,PMMG_Grp,"Deallocating listgrp container");
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
  PMMG_parmesh_int_comm_free( parmesh, parmesh->int_edge_comm );
  PMMG_parmesh_int_comm_free( parmesh, parmesh->int_face_comm );

  PMMG_parmesh_ext_comm_free( parmesh, parmesh->ext_node_comm, parmesh->next_node_comm );
  PMMG_DEL_MEM(parmesh, parmesh->ext_node_comm, parmesh->next_node_comm,
            PMMG_Ext_comm, "ext node comm");
  PMMG_parmesh_ext_comm_free( parmesh, parmesh->ext_edge_comm, parmesh->next_edge_comm );
  PMMG_DEL_MEM(parmesh, parmesh->ext_edge_comm, parmesh->next_edge_comm,
            PMMG_Ext_comm, "ext edge comm");
  PMMG_parmesh_ext_comm_free( parmesh, parmesh->ext_face_comm, parmesh->next_face_comm );
  PMMG_DEL_MEM(parmesh, parmesh->ext_face_comm, parmesh->next_face_comm,
            PMMG_Ext_comm, "ext face comm");
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * Free parmesh listgrp that are allocated
 */
void PMMG_parmesh_Free_Listgrp( PMMG_pParMesh parmesh )
{
  PMMG_listgrp_free( parmesh, &parmesh->listgrp, parmesh->ngrp );

  PMMG_DEL_MEM(parmesh,parmesh->listgrp,1,PMMG_Grp,"deallocating groups container");
}

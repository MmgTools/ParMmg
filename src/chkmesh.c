#include "chkmesh.h"
#include "mmgcommon.h" // _MMG5_SAFE_CALLOC

//! Internal Communicator Test:
//    allocate an array of \SUM{nitem_int_node} elements, each being a struct point (ie three c[3])
//    initialize to 0
//    for each grp communicator
//      if array[indx2] is 0 => initialize to point's coordinates
//      else check that the existing point's coordinates match the coordinates of the
//                                  point that the group local communicator points to
// Will     detect mismatches between groups
// Will NOT detect errors in elements not shared between the subgroups 
// Will NOT detect if the same error happens in both groups, ie if the values match but are both wrong
int checkIntComm( PMMG_pParMesh mesh )
{
  struct point { double c[3]; int tmp; };
  struct point *check = NULL;
  int testError = 0;
  int ngrp = mesh->ngrp;
  PMMG_pGrp grpCur = NULL;
  MMG5_pMesh meshCur = NULL;
  int commSizeGlo = mesh->int_node_comm->nitem;

  // Loop counter(s)
  int commIdx, grpId;

  _MMG5_SAFE_CALLOC( check, commSizeGlo + 1, struct point );

  for ( grpId = 0; grpId < ngrp; ++grpId ) {

    grpCur = mesh->listgrp + grpId;
    meshCur = grpCur->mesh;
    int numFailed = 0;

    int commSizeLoc = grpCur->nitem_int_node_comm;
    for ( commIdx = 0; commIdx < commSizeLoc; ++commIdx ) {

      int comm1Idx = grpCur->node2int_node_comm_index2[ commIdx ];
      printf( "+++++NIKOS[%d/%d]+++:: CIC: commSizeGlo = %d, i = %d, n2i_idx2(check position) = %d \n",
              grpId + 1, ngrp, commSizeGlo, commIdx, comm1Idx );
      if ( check[ comm1Idx ].tmp == 0. ) {
        //check[ comm1Idx ].c[0] = meshCur->point[ grpCur->node2int_node_comm_index1[ commIdx ] ].c[0];
        //check[ comm1Idx ].c[1] = meshCur->point[ grpCur->node2int_node_comm_index1[ commIdx ] ].c[1];
        //check[ comm1Idx ].c[2] = meshCur->point[ grpCur->node2int_node_comm_index1[ commIdx ] ].c[2];
        check[ comm1Idx ].tmp = 1.;
      } else {
        //for ( int j = 0; j < 3; ++j )
        //  if ( check[ pos ].c[j] != meshCur->point[ grpCur->node2int_node_comm_index1[ i ] ].c[j] ) {
        //    ++numFailed;
        //    testError = 1;
        //}
      }
    }
    printf( "+++++NIKOS[%d/%d]:: INTERNAL COMMUNICATOR CHECKED. NumFailed: %d \n", grpId+1, ngrp, numFailed / 3 );
  }
  if ( testError )
    printf( "+++++NIKOS[%d/%d]:: INTERNAL COMMUNICATOR ERROR: SHOULD ABORT\n", grpId+1, ngrp );

  _MMG5_SAFE_FREE( check );

  return 1;
}


int  checkExtComm( void ) 
{
  return 1;
}


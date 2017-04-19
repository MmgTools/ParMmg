#include "chkmesh.h"
#include "mmgcommon.h" // _MMG5_SAFE_CALLOC

//! Check Internal Communicator for each group:
//
//    First test:
//      allocate an array of int_node_comm->item_int_node struct point (ie three c[3]) elements (initialized to 0)
//      for each grp communicator
//        if array[indx2] is empty => initialize to current group's point's coordinates
//        else check that the existing point's coordinates match the  current group's point's coordinates
//    Will detect mismatches between groups but will NOT detect:
//      errors in elements not shared between the subgroups
//      if the same error happens in both groups, ie if the values match but are both wrong
//
//    Second test:
//NIKOS TODO: ADD CHECK: elements should appear one and only one per mesh local communicator
//
//    Third test:
//NIKOS TODO: ADD CHECK: all elements in commSizeGlo
//                         a) should appear in per mesh local communicators at least once
//                         b) should have no "gaps"
int checkIntComm( PMMG_pParMesh mesh )
{
  struct point { double c[3]; int tmp; };
  struct point *check = NULL;
  int testError = 0;
  int ngrp = mesh->ngrp;
  PMMG_pGrp grpCur = NULL;
  MMG5_pMesh meshCur = NULL;
  int commSizeGlo = mesh->int_node_comm->nitem;
  int commIdx2 = 0;
  int commIdx1 = 0;
  int numFailed = 0;
  int commSizeLoc = 0;

  // Loop counter(s)
  int commIdx, grpId;

  _MMG5_SAFE_CALLOC( check, commSizeGlo + 1, struct point );

  for ( grpId = 0; grpId < ngrp; ++grpId ) {

    grpCur = mesh->listgrp + grpId;
    meshCur = grpCur->mesh;
    numFailed = 0;

    commSizeLoc = grpCur->nitem_int_node_comm;
    for ( commIdx = 0; commIdx < commSizeLoc; ++commIdx ) {

      commIdx2 = grpCur->node2int_node_comm_index2[ commIdx ];
      commIdx1 = grpCur->node2int_node_comm_index1[ commIdx ];
//      printf( "+++++NIKOS[%d/%d]+++:: CIC: commSizeGlo = %d, i = %d, n2i_idx2(check position) = %d ",
//              grpId + 1, ngrp, commSizeGlo, commIdx, commIdx2);
      if ( check[ commIdx2 ].tmp == 0. ) {
//        printf( "+++adding (%f,%f,%f) to (%f,%f,%f) \n",
//                meshCur->point[ commIdx1 ].c[0],
//                meshCur->point[ commIdx1 ].c[1],
//                meshCur->point[ commIdx1 ].c[2],
//                check[ commIdx2 ].c[0], check[ commIdx2 ].c[1], check[ commIdx2 ].c[2] );
        check[ commIdx2 ].c[0] = meshCur->point[ commIdx1 ].c[0];
        check[ commIdx2 ].c[1] = meshCur->point[ commIdx1 ].c[1];
        check[ commIdx2 ].c[2] = meshCur->point[ commIdx1 ].c[2];
        check[ commIdx2 ].tmp = 1.;
      } else {
//        printf( "---checking (%f,%f,%f) to (%f,%f,%f) \n",
//                meshCur->point[ commIdx1 ].c[0],
//                meshCur->point[ commIdx1 ].c[1],
//                meshCur->point[ commIdx1 ].c[2],
//                check[ commIdx2 ].c[0], check[ commIdx2 ].c[1], check[ commIdx2 ].c[2] );
        for ( int j = 0; j < 3; ++j )
          if ( check[ commIdx2 ].c[j] != meshCur->point[ commIdx1 ].c[j] ) {
            ++numFailed;
            testError = 1;
        }
      }
    }
    printf( "+++++NIKOS[%d/%d]:: INTERNAL COMMUNICATOR CHECKED. NumFailed: %d \n", ngrp, grpId+1, numFailed / 3 );
  }

  _MMG5_SAFE_FREE( check );

  return testError;
}


int  checkExtComm( void )
{
  return 1;
}

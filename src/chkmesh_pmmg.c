#include "parmmg.h"
#include "chkmesh_pmmg.h"

//! Check Internal Communicator for each group.
//    3 checks are implemented, see inside the function for specifics
int PMMG_checkIntComm( PMMG_pParMesh mesh )
{
  //NIKOS TODO: could try using std::vector<bool> and std::map in these tests to speed up things
  //NIKOS TODO: BOUND CHECK ARRAY ACCESSES: can't rely that the contents of the communicators are valid
  struct point { double c[3]; int tmp; };
  struct point *check = NULL;
  int testFailed = 0;
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

  // FIRST TEST:
  //   allocate an array of int_node_comm->item_int_node struct point (ie three c[3]) elements (initialized to 0)
  //   for each grp communicator
  //     if array[indx2] is empty => initialize to current group's point's coordinates
  //     else check that the existing point's coordinates match the  current group's point's coordinates
  // Will detect mismatches between groups but will NOT detect:
  //   errors in elements not shared between the subgroups
  //   if the same error happens in both groups, ie if the values match but are both wrong
  PMMG_CALLOC( mesh, check, commSizeGlo + 1, struct point, "Allocating check space: " );

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
            testFailed = 1;
        }
      }
    }
    printf( "+++++NIKOS[%d/%d]:: INTERNAL COMMUNICATOR CHECKED: number of failed communicator elements: %d \n", ngrp, grpId+1, numFailed / 3 );
  }



  // SECOND TEST:
  //   elements should appear once and only once per mesh local communicator
  for ( grpId = 0; grpId < ngrp; ++grpId ) {
    grpCur = mesh->listgrp + grpId;
    meshCur = grpCur->mesh;
    commSizeLoc = grpCur->nitem_int_node_comm;

    for ( commIdx = 0; commIdx < commSizeGlo + 1; ++commIdx )
      check[ commIdx ].tmp = 0;

    for ( commIdx = 0; commIdx < commSizeLoc; ++commIdx ) {
      commIdx2 = grpCur->node2int_node_comm_index2[ commIdx ];
      if ( check[ commIdx2 ].tmp == 0. )
        check[ commIdx2 ].tmp = 1.;
      else
        testFailed = 1;
    }
  }



  // THIRD TEST:
  //   all elements in commSizeGlo
  //     a) should appear in per mesh local communicators at least once
  //     b) should have no "gaps"
//THINK AGAIN ABOUT IT AND FINISH IT  for ( commIdx = 0; commIdx < commSizeGlo + 1; ++commIdx )
//THINK AGAIN ABOUT IT AND FINISH IT    check[ commIdx ].tmp = 0;
//THINK AGAIN ABOUT IT AND FINISH IT  for ( grpId = 0; grpId < ngrp; ++grpId ) {
//THINK AGAIN ABOUT IT AND FINISH IT    grpCur = mesh->listgrp + grpId;
//THINK AGAIN ABOUT IT AND FINISH IT    meshCur = grpCur->mesh;
//THINK AGAIN ABOUT IT AND FINISH IT    commSizeLoc = grpCur->nitem_int_node_comm;
//THINK AGAIN ABOUT IT AND FINISH IT    for ( commIdx = 0; commIdx < commSizeLoc; ++commIdx ) {
//THINK AGAIN ABOUT IT AND FINISH IT      commIdx2 = grpCur->node2int_node_comm_index2[ commIdx ];
//THINK AGAIN ABOUT IT AND FINISH IT      if ( check[ commIdx2 ].tmp == 0. )
//THINK AGAIN ABOUT IT AND FINISH IT        check[ commIdx2 ].tmp = 1.;
//THINK AGAIN ABOUT IT AND FINISH IT    }
//THINK AGAIN ABOUT IT AND FINISH IT  }
//THINK AGAIN ABOUT IT AND FINISH IT  commIdx = 0;
//THINK AGAIN ABOUT IT AND FINISH IT  // Skip non zero elements
//THINK AGAIN ABOUT IT AND FINISH IT  while ( check[ ++commIdx ].tmp != 0 )
//THINK AGAIN ABOUT IT AND FINISH IT    ;
//THINK AGAIN ABOUT IT AND FINISH IT  if ( commIdx != mesh->int_node_comm->nitem )
//THINK AGAIN ABOUT IT AND FINISH IT    testFailed = 1;
//THINK AGAIN ABOUT IT AND FINISH IT  printf( "+++++NIKOS[%d/%d]:: commIdx= %d, nitem= %d \n", ngrp, grpId+1, commIdx, mesh->int_node_comm->nitem );

  // NIKOS TODO: CHECK THE idx1/idx2 pairs ?

  PMMG_FREE( mesh, check, (commSizeGlo + 1) * sizeof(struct point), "Deallocating check space:" );
  return testFailed;
}


int PMMG_checkExtComm( void )
{
  return 1;
}

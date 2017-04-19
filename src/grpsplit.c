/**
 * \file grpsplit.c
 * \brief Split groups into sub groups.
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \author Nikos Pattakos (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */
#include "libparmmgtypes.h" // PMMG_pGrp
#include "metis.h" // idx_t
#include "mmgcommon.h" // _MMG5_SAFE_CALLOC
#include "mmg3d.h" //_MMG5_idir[4][3]
#include "libparmmg.h" // PMMG_mesh2metis
#include "grpsplit.h"

// Subgroups target size. It is chosen arbitrarily to help assist the remesher work faster
static const int REMESHER_TARGET_MESH_SIZE = 1024; // * 128;

static int HowManyGroups ( const int nelem )
{
  int ngrp = nelem / REMESHER_TARGET_MESH_SIZE;

  if ( ngrp == 0 )
    return ( 1 );
  else if ( ngrp * REMESHER_TARGET_MESH_SIZE < nelem )
    return ( ngrp + 1 );

  return ( ngrp );
}


//NIKOS TODO: there is some code duplication. eg this metis calling is also in the metisfunction.c

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Split one group into several groups.
 *
 */
int PMMG_splitGrps( PMMG_pParMesh parmesh )
{
  PMMG_pGrp grpOld = parmesh->listgrp;
  PMMG_pGrp grpsNew = NULL;
  PMMG_pGrp grpCur = NULL;
  MMG5_pMesh meshOld = parmesh->listgrp->mesh;
  MMG5_pMesh meshCur = NULL;
  MMG5_pTetra tetraCur = NULL; // pointer to the tetra being processed in the subgroup mesh struct
  MMG5_pxTetra pxt;
  int *countPerGrp = NULL;

  idx_t ngrp = 1;
  idx_t *part = NULL;
  idx_t *xadj = NULL;
  idx_t *adjncy = NULL;
  //NIKOS TODO: experiment with the number of balancing constraints
  idx_t ncon = 1;/*number of balancing constraint*/
  idx_t nelt = meshOld->ne;
  idx_t objval;
  int ier,j;

  /* counters for tetra, point, while constructing a subgroup */
  int tetPerGrp = 0;
  int poiPerGrp = 0;
  int *adja = NULL;
  int vindex = 0;
  int adjidx = 0;

  // Loop counter vars
  int i, grpId, poi, tet, fac;

  printf( "+++++NIKOS+++++[%d/%d]: mesh has: %d(%d) #points, %d(%d) #edges, %d(%d) #tria and %d(%d) tetras(elements)\n",
          parmesh->myrank+1, parmesh->nprocs, meshOld->np, meshOld->npi, meshOld->na, meshOld->nai, meshOld->nt, meshOld->nti, meshOld->ne, meshOld->nei );

//NIKOS TODO _MMG3D_bdryBuild( meshOld ); //note: no error checking
//NIKOS TODO MMG3D_saveMesh( meshOld, "mesh-ORIG.mesh" );
  ngrp = HowManyGroups( meshOld->ne );
  /* Does the group need to be further subdivided to subgroups or not? */
  if ( ngrp == 1 )  {
    printf( "[%d-%d]: %d group is enough, no need to create sub groups.\n", parmesh->myrank+1, parmesh->nprocs, ngrp );
    return 1;
  } else {
    printf( "[%d-%d]: %d groups required, splitting into sub groups...\n", parmesh->myrank+1, parmesh->nprocs, ngrp );
  }



  // use metis to partition the mesh into the computed number of groups needed
  // part array contains the groupID computed by metis for each tetra
  _MMG5_SAFE_CALLOC( part, meshOld->ne, idx_t );

  PMMG_mesh2metis( parmesh, &xadj, &adjncy );

  ier =  METIS_PartGraphKway( &nelt          , &ncon          , xadj , adjncy        , NULL/*vwgt*/ ,
                              NULL/*vsize*/  , NULL/*adjwgt*/ , &ngrp, NULL/*tpwgts*/, NULL/*ubvec*/,
                              NULL/*options*/, &objval, part );
  //ier =  METIS_PartGraphRecursive( &nelt, &ncon, xadj, adjncy, NULL/\*vwgt*\/, NULL/\*vsize*\/,
  //                                 NULL/\*adjwgt*\/, &nproc, NULL/\*tpwgts*\/,
  //                                 NULL/\*ubvec*\/, NULL/\*options*\/, &objval, part );
  if ( ier != METIS_OK ) {
    fprintf(stderr, "Metis returned error value: " );
    switch ( ier ) {
      case METIS_ERROR_INPUT:
        fprintf(stderr, "METIS_ERROR_INPUT: input data error\n" );
        break;
      case METIS_ERROR_MEMORY:
        fprintf(stderr, "METIS_ERROR_MEMORY: could not allocate memory error\n" );
        break;
      case METIS_ERROR:
        fprintf(stderr, "METIS_ERROR: generic error\n" );
        break;
    }
    _MMG5_SAFE_FREE( xadj );
    _MMG5_SAFE_FREE( adjncy );
    _MMG5_SAFE_FREE( part );
    return (0);
  }



  /* count_per_grp: count elements per group */
  _MMG5_SAFE_CALLOC( countPerGrp, ngrp, int );
  for ( tet = 0; tet < meshOld->ne ; ++tet )
    ++countPerGrp[ part[ tet ] ];
  for ( i = 0; i < ngrp ; i++ )
    printf( "+++++NIKOS+++++[%d/%d]: group[%d] has %d elements\n", parmesh->myrank+1, parmesh->nprocs, i, countPerGrp[i] );



  /* Allocate list of subgroups struct and allocate memory*/
  _MMG5_SAFE_CALLOC( grpsNew, ngrp, PMMG_Grp );
  //NIKOS TODO: do MMG3D_Init_mesh/MMG3D_Set_meshSize offer sth? otherwise allocate on your own
  for ( grpId = 0; grpId < ngrp; ++grpId ) {
    grpCur = &grpsNew[grpId];
    grpCur->mesh = NULL;
    grpCur->sol  = NULL;
    grpCur->disp = NULL;
    //_MMG5_SAFE_CALLOC( ptr, size, type );
    //_MMG5_SAFE_REALLOC( ptr, size, type, message );
    MMG3D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &grpCur->mesh,
                     MMG5_ARG_ppMet, &grpCur->sol, MMG5_ARG_end );
    //grpCur->mesh->info.mem = 300;
    MMG3D_Set_meshSize( grpCur->mesh, countPerGrp[grpId], countPerGrp[grpId], 0, 0, 0, 0 );

    meshCur = grpCur->mesh;
    meshCur->xtmax = meshOld->xtmax;
    _MMG5_ADD_MEM( meshCur, ( meshCur->xtmax + 1 ) * sizeof(MMG5_xTetra), "boundary tetrahedra",
                   fprintf( stderr,"  Exit program.\n" );
                   exit( EXIT_FAILURE ) );
    _MMG5_SAFE_CALLOC( meshCur->xtetra, meshCur->xtmax + 1, MMG5_xTetra );

    /* memory to store normals for boundary points */
    meshCur->xpmax  = meshCur->npmax;
    _MMG5_ADD_MEM(meshCur, ( meshCur->xpmax + 1 ) * sizeof(MMG5_xPoint), "boundary points", return(0) );
    _MMG5_SAFE_CALLOC( meshCur->xpoint, meshCur->xpmax + 1, MMG5_xPoint );
    _MMG5_ADD_MEM( meshCur, ( 4 * meshCur->nemax + 5 ) * sizeof(int), "adjacency table",
                   fprintf( stderr,"  Exit program.\n" );
                   exit( EXIT_FAILURE ) );
    _MMG5_SAFE_CALLOC( meshCur->adja, 4 * meshCur->nemax + 5, int );

    _MMG5_SAFE_CALLOC( grpCur->node2int_node_comm_index1, meshCur->ne / 3, int );
    _MMG5_SAFE_CALLOC( grpCur->node2int_node_comm_index2, meshCur->ne / 3, int );
  }



/////////////////// algean algo /////////////////////////////
//  // Initialize flag/tmp/etc vars
//  for ( poi = 1; poi < meshOld->np + 1; ++poi )
//    meshOld->point[poi].flag = 0;
//  // Will use the MMG5_Point.tmp field to store the position of each vertex in
//  // the node2index2 communicator buffer or -1.
//  for ( poi = 1; poi < meshOld->np + 1; ++poi )
//    meshOld->point[poi].tmp = -1;
//  for ( i = 0; i < grpCur->nitem_int_node_comm; ++i )
//    meshOld->point[ grpOld->node2int_node_comm_index1 [i ] ].tmp = grpOld->node2int_node_comm_index2[ i ];
//
//  // Process tetrahedrons
//  for ( tet = 1; tet < meshOld->ne + 1; ++tet ) {
//    int gr_idx = part[ tet - 1 ];
//    MMG5_pTetra pt = &meshOld->tetra[ tet ];
//
//    for ( int fac = 0 ; fac < 4 ; ++fac ) {
//
//      int adja_gr_idx = part[ meshOld->adja[ 4 * ( tet - 1 ) + 1 + fac] - 1 ];
//      if ( adja_gr_idx != gr_idx ) {
//        for ( int ver = 0; ver < 3; ++ver ) {
//          ip_loc =
//          ip =
//          ppt = &meshOld->point[ip];
//        }
//      }
//    }
//  }
/////////////////// algean algo /////////////////////////////

  //NIKOS TODO: LOOP OVER part ngrp TIMES or USE A tmp[NGROUPS][NP] ARRAY AND LOOP ONLY ONCE?
  //            it wastes memory (eg 10 groups x 100k tetra = 4Mb of ints) but only loops over part once
  // use point[].tmp field to "remember" index in internal communicator of vertices. specifically:
  //   place a copy of vertices' node2index2 position at point[].tmp field or -1 if they are not in the comm
  for ( poi = 1; poi < meshOld->np + 1; ++poi )
    meshOld->point[poi].tmp = -1;
  for ( i = 0; i < grpOld->nitem_int_node_comm; i++ )
    meshOld->point[ grpOld->node2int_node_comm_index1[ i ] ].tmp = grpOld->node2int_node_comm_index2[ i ];
  for ( grpId = 0 ; grpId < ngrp ; ++grpId ) {
    grpCur = grpsNew + grpId;
    meshCur = grpCur->mesh;

    //NIKOS TODO: this could be replaced by a vector<bool> flag
    // use point[].flag field to "remember" assigned local(ie in subgroup) numbering
    for ( poi = 1; poi < meshOld->np + 1; ++poi )
      meshOld->point[poi].flag = 0;

    tetPerGrp = 0;
    poiPerGrp = 0;
    for ( tet = 1; tet < meshOld->ne + 1; ++tet ) {
      // MMG3D_Tetra.flag is used to update adjacency vector: 
      //   if the tetra belongs to the group we store the local tetrahedron id in tet.flag
      meshOld->tetra[tet].flag = 0; 

      // Skip elements that do not belong in the group processed in this iteration
      if ( grpId != part[ tet - 1 ] )
        continue;

      ++tetPerGrp;
      tetraCur = meshCur->tetra + tetPerGrp;
      meshOld->tetra[tet].flag = tetPerGrp;

      // add tetrahedron to subgroup (copy from original group)
      memcpy( tetraCur, &meshOld->tetra[tet], sizeof(MMG5_Tetra) );
      tetraCur->base = 0;
      tetraCur->mark = 0;
      tetraCur->flag = 0;
      //NIKOS TODO: there are also cases where we have to update existing boundaries,
      //            eg when boundary tetras are also interfaced with new groups
      // Handle xTetras
      if ( tetraCur->xt != 0 ) {
        ++meshCur->xt;
        if ( meshCur->xt > meshCur->xtmax ) {
          /* realloc of xtetra table */
          _MMG5_TAB_RECALLOC(meshCur,meshCur->xtetra,meshCur->xtmax,0.2,MMG5_xTetra,
                             "larger xtetra table",
                             meshCur->xt--;
                             return(0));
        }
        memcpy( meshCur->xtetra + meshCur->xt, &meshOld->xtetra[ meshOld->tetra[tet].xt], sizeof(MMG5_xTetra) );
        tetraCur->xt = meshCur->xt;
        //    } else {
        //      printf ("there may be more: new boundary between subgroups\n");
      }

      //printf ( "+++++NIKOS[%d/%d]:: tetra %4d(%4d), mesh xt: %3d\n", grpId+1, ngrp, tet, tetPerGrp, meshOld->tetra[tet].xt);
      /* Add tetrahedron vertices in points struct and adjust tetrahedron vertices indices */
      for ( poi = 0; poi < 4 ; ++poi ) {
        // if it is the first time this point is seen in this subgroup
        if ( 0 == meshOld->point[ meshOld->tetra[tet].v[poi] ].flag ) {
  
          // Add point in subgroup point array
          ++poiPerGrp;
          memcpy( meshCur->point + poiPerGrp, &meshOld->point[ meshOld->tetra[tet].v[poi] ], sizeof(MMG5_Point) );
  
          // update tetra vertex reference
          tetraCur->v[poi] = poiPerGrp;
  
          // "Remember" the assigned subgroup point id
          meshOld->point[ meshOld->tetra[tet].v[poi] ].flag = poiPerGrp;
  
          // Add point in subgroup's communicator if it already was in group's communicator
          if ( meshCur->point[ poiPerGrp ].tmp != -1 ) {
            grpCur->node2int_node_comm_index1[ grpCur->nitem_int_node_comm ] = poiPerGrp;
            grpCur->node2int_node_comm_index2[ grpCur->nitem_int_node_comm ] = meshCur->point[ poiPerGrp ].tmp;
            ++grpCur->nitem_int_node_comm;
          }

          // Handle xPoints
          if ( meshCur->point[poiPerGrp].xp != 0 ) {
            ++meshCur->xp;
            memcpy( meshCur->xpoint + meshCur->xp, &meshOld->xpoint[ meshOld->point[ meshOld->tetra[tet].v[poi] ].xp ], sizeof(MMG5_xPoint) );
            meshCur->point[poiPerGrp].xp = meshCur->xp;
            //        } else {
            //        printf( " some other more border cases, think about it \n" );
          }
        // point is already included in this subgroup
        } else {
          // update tetra vertex reference
          tetraCur->v[poi] = meshOld->point[ meshOld->tetra[tet].v[poi] ].flag;
        }
        //printf ( "adja: %5d, point[%4d-%4d-%d].xt: %d \n", meshOld->adja[ 4 * (tet - 1) + 1 + poi ], meshOld->tetra[tet].v[poi], tetPerGrp, poi, meshOld->point[ meshOld->tetra[tet].v[poi] ].xp);
      }

      adja = &meshCur->adja[ 4 * ( tetPerGrp - 1 ) + 1 ];
      memcpy( adja, &meshOld->adja[4 * ( tet - 1 ) + 1 ], 4 * sizeof(int) );
      for ( fac = 0; fac < 4; ++fac ) {
        if ( adja[ fac ] != 0 ) {
          adjidx = adja[ fac ] / 4;
          vindex = adja[ fac ] % 4;
          // new boundary face, set to 0
          if ( part[ adjidx - 1 ] != grpId ) {
            adja[ fac ] = 0;

            /*creation of the interface faces : ref 0 and tag MG_PARBDY*/
            if( !meshCur->tetra[tet].xt ) {
              /*create a new xt*/
              ++meshCur->xt;
              if ( meshCur->xt > meshCur->xtmax ) {
                /* realloc of xtetra table */
                _MMG5_TAB_RECALLOC(meshCur,meshCur->xtetra,meshCur->xtmax,0.2,MMG5_xTetra,
                                   "larger xtetra table",
                                   meshCur->xt--;
                                   return(0));
              }
              memcpy( meshCur->xtetra + meshCur->xt, &meshOld->xtetra[ meshOld->tetra[tet].xt], sizeof(MMG5_xTetra) );
              tetraCur->xt = meshCur->xt;
            }
            pxt = &meshCur->xtetra[tetraCur->xt];
            pxt->ref[fac] = 0;
            pxt->ftag[fac] |= (MG_PARBDY + MG_BDY + MG_REQ);

            for ( j=0; j<3; ++j )
              pxt->tag[_MMG5_iarf[fac][j]] |= (MG_PARBDY + MG_BDY + MG_REQ);

            // if the adjacent number is already processed
          } else if ( adjidx < tet ) { //NIKOS TODO: ie if meshOld->tetra[tet].flag != 0
            adja[i] =  4 * meshOld->tetra[adjidx].flag  + vindex;
            // Now that we know the local gr_idx for both tetra, also update the adjacency entry for the adjacent face
            meshCur->adja[ 4 * ( meshOld->tetra[adjidx].flag - 1) + 1 + vindex ] = 4 * tetPerGrp + fac;
          }
        }
      }
      adja = &meshOld->adja[ 4 * ( tet - 1 ) + 1 ];
      for ( fac = 0; fac < 4; ++fac ) {
        adjidx = adja[ fac ] / 4;

        if ( adjidx && grpId != part[ adjidx - 1 ] ) {
          for ( poi = 0; poi < 3; ++poi ) {
            if ( meshCur->point[ tetraCur->v[ _MMG5_idir[fac][poi] ] ].tmp == -1 ) {
              grpCur->node2int_node_comm_index1[ grpCur->nitem_int_node_comm ] = tetraCur->v[ _MMG5_idir[fac][poi] ];
              grpCur->node2int_node_comm_index2[ grpCur->nitem_int_node_comm ] = grpCur->nitem_int_node_comm + 1;
              meshOld->point[ meshOld->tetra[tet].v[ _MMG5_idir[fac][poi]  ] ].tmp = grpCur->nitem_int_node_comm + 1;
              meshCur->point[ tetraCur->v[ _MMG5_idir[fac][poi] ] ].tmp = grpCur->nitem_int_node_comm + 1;
              ++grpCur->nitem_int_node_comm;
            }
          }
        }
      }
    }
    assert( (meshCur->ne == tetPerGrp) && "Error in PMMG_splitGrps" );
    printf( "+++++NIKOS[%d/%d]:: %d points in group, %d tetra (expected: %d)ed.%d nitem in int communicator\n",
      grpId+1, ngrp, poiPerGrp, tetPerGrp, meshCur->ne, grpCur->nitem_int_node_comm );
  }





//NIKOS TODO  for ( grpId = 0 ; grpId < ngrp ; grpId++ ) {
//NIKOS TODO    meshCur = grpsNew[grpId].mesh;
//NIKOS TODO    //MMG5_saveMshMesh( meshCur, mesMMG5_pSol sol,const char *filename)
//NIKOS TODO    char name[20];
//NIKOS TODO    sprintf( name, "mesh-p%d-%02d.mesh", parmesh->myrank+1, grpId );
//NIKOS TODO    _MMG3D_bdryBuild(meshCur); //note: no error checking
//NIKOS TODO    MMG3D_saveMesh( meshCur, name );
//NIKOS TODO  }

#ifndef NDEBUG
#warning NIKOS: move this to a unit test and find a nicer way to integrate this here
  //! Communicator Test:
  //    allocate an array of \SUM{nitem_int_node} elements, each being a struct point (ie three c[3])
  //    initialize to 0
  //    for each grp communicator
  //      if array[indx2] is 0 => initialize to point's coordinates
  //      else check that the existing point's coordinates match the coordinates of the
  //                                  point that the group local communicator points to
  // Will     detect mismatches between groups
  // Will NOT detect errors in elements not shared between the subgroups 
  // Will NOT detect if the same error happens in both groups, ie if the values match but are both wrong
  struct point { double c[3]; int tmp; };
  struct point *check = NULL;
  int tot_index = 0;
  int testError = 0;
  for ( grpId = 0; grpId < ngrp; ++grpId )
    tot_index += grpsNew[grpId].nitem_int_node_comm;
  //  printf( "+++++NIKOS[%d/%d]:: allocating to check up to %d int node comm elements \n", grpId+1, ngrp, tot_index );
  _MMG5_SAFE_CALLOC( check, tot_index + 1, struct point );
  for ( grpId = 0; grpId < ngrp; ++grpId ) {
    grpCur = grpsNew + grpId;
    meshCur = grpCur->mesh;
    int numFailed = 0;
    for ( i = 0; i < grpCur->nitem_int_node_comm; ++i ) {
      int pos = grpCur->node2int_node_comm_index2[ i ];
      //printf( "+++++NIKOS[%d/%d]:: CIC: nitem_iterator=%d, n2i_n_idx2=%d, check[%d].c[0]=%f against element %d in local group with value %f\n", 
      //        grpId + 1, ngrp,    i,    pos,   pos, check[ pos ].c[0],
      //        grpCur->node2int_node_comm_index1[ i ], meshCur->point[ grpCur->node2int_node_comm_index1[ i ] ].c[0] );
      //printf( "+++++NIKOS[%d/%d]:: CHECKING INTERNAL COMMUNICATOR: [%d/%d] :: idx2 = %d .",
      //        grpId+1, ngrp, 
      //        i, grpCur->nitem_int_node_comm, 
      //        pos );
      if ( check[ pos ].tmp == 0. ) {
        check[ pos ].c[0] = meshCur->point[ grpCur->node2int_node_comm_index1[ i ] ].c[0];
        check[ pos ].c[1] = meshCur->point[ grpCur->node2int_node_comm_index1[ i ] ].c[1];
        check[ pos ].c[2] = meshCur->point[ grpCur->node2int_node_comm_index1[ i ] ].c[2];
        check[ pos ].tmp = 1.;
        //printf( "\t adding %d: (%f,%f,%f) \n", pos, check[ pos ].c[0], check[ pos ].c[1], check[ pos ].c[2] );
      } else {
        // printf( "\t comparing %d: (%f,%f,%f) to (%f,%f,%f) \n", pos, check[ pos ].c[0], check[ pos ].c[1], check[ pos ].c[2], meshCur->point[ grpCur->node2int_node_comm_index1[ i ] ].c[0] , meshCur->point[ grpCur->node2int_node_comm_index1[ i ] ].c[1] , meshCur->point[ grpCur->node2int_node_comm_index1[ i ] ].c[2]  );
        for ( int j = 0; j < 3; ++j )
          if ( check[ pos ].c[j] != meshCur->point[ grpCur->node2int_node_comm_index1[ i ] ].c[j] ) {
            ++numFailed;
            testError = 1;
            printf( "+++++NIKOS[%d/%d]:: CHECKING INTERNAL COMMUNICATOR FAILED: check[%d].c[%d](=%f) != meshCur->point[ grpCur->node2int_node_comm_index1[ %d ](=%d) ].c[%d](=%f) \n",
                    grpId + 1, ngrp,
                    pos, j, check[ pos ].c[j],
                    i, j, grpCur->node2int_node_comm_index1[ i ], meshCur->point[ grpCur->node2int_node_comm_index1[ i ] ].c[j] );
          }
      }
    }
       printf( "+++++NIKOS[%d/%d]:: INTERNAL COMMUNICATOR CHECKED. NumFailed: %d \n", grpId+1, ngrp, numFailed / 3 );
  }
  if ( testError )
    printf( "+++++NIKOS[%d/%d]:: INTERNAL COMMUNICATOR ERROR: SHOULD ABORT\n", grpId+1, ngrp );
  _MMG5_SAFE_FREE( xadj );
#endif

  //#error NIKOS: CHANGE THE MEMORY ALLOCATIONS WITH PROPER ALLOCATION+REALLOCATION ??
  //MMG3D_Free_all( MMG5_ARG_start,
  //                MMG5_ARG_ppMesh, &(parmesh->listgrp->mesh), MMG5_ARG_ppMet, &(parmesh->listgrp->sol),
  //                MMG5_ARG_end);
  //parmesh->listgrp = grpsNew;
  //parmesh->ngrp = ngrp;

  _MMG5_SAFE_FREE( xadj );
  _MMG5_SAFE_FREE( adjncy );
  _MMG5_SAFE_FREE( part );
  _MMG5_SAFE_FREE( countPerGrp );
  return( 1 );
}

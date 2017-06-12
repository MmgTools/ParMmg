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
#include "parmmg.h"
#include "metis_pmmg.h"
#include "libparmmgtypes.h" // PMMG_pGrp
#include "mmg3d.h" //_MMG5_idir[4][3]
#include "grpsplit_pmmg.h"
#include "chkmesh_pmmg.h"

// Subgroups target size. Chosen arbitrarily to help assist the remesher work faster
static const int REMESHER_TARGET_MESH_SIZE = 128 * 1024;

static int howManyGroups ( const int nelem )
{
  int ngrp = nelem / REMESHER_TARGET_MESH_SIZE;

  if ( ngrp == 0 )
    return ( 1 );
  else if ( ngrp * REMESHER_TARGET_MESH_SIZE < nelem )
    return ( ngrp + 1 );

  return ( ngrp );
}


/**
 * \param to       mesh to copy xtetra item to
 * \param from     mesh to copy xtetra item from
 * \param location location of xtetra to copy in original mesh
 *
 *  Append a tetraedron to the mesh, increasing the array if the allocated
 *  array is not big enough.
 */
static int xtetraAppend( MMG5_pMesh to, MMG5_pMesh from, int location )
{
  // Adjust the value of scale to reallocate memory more or less agressively
  const float scale = 2.f;
  if ( (to->xt + 1) > to->xtmax ) {
    PMMG_RECALLOC(to, to->xtetra, scale*to->xtmax+1, to->xtmax+1,MMG5_xTetra,
                  "larger xtetra table",return PMMG_FAILURE);
    to->xtmax = scale * to->xtmax;
  }
  ++to->xt;
  memcpy( &to->xtetra[ to->xt ],
          &from->xtetra[ from->tetra[ location ].xt ], sizeof(MMG5_xTetra) );
  return PMMG_SUCCESS;
}
/**
 * \param to       mesh to copy xpoint item to
 * \param from     mesh to copy xpoint item from
 * \param location location of xpoint to copy in original mesh
 * \param point    xpoint's number in tetrahedron
 *
 *  Append a point in the point list, increasing the array if the allocated
 *  array is not big enough.
 */
static int xpointAppend( MMG5_pMesh to, MMG5_pMesh from, int tetrahedron, int point )
{
  // Adjust the value of scale to reallocate memory more or less agressively
  const float scale = 2.f;
  if ( (to->xp + 1) > to->xpmax + 1 ) {
    PMMG_RECALLOC(to, to->xpoint, scale*to->xpmax+1, to->xpmax+1,MMG5_xPoint,
                  "larger xpoint table",return PMMG_FAILURE);
    to->xpmax = scale * to->xpmax;
  }
  ++to->xp;
  memcpy( &to->xpoint[ to->xp ],
          &from->xpoint[ from->point[ from->tetra[tetrahedron].v[point] ].xp ],
          sizeof(MMG5_xPoint) );
  return PMMG_SUCCESS;
}
/**
 * \param parmesh parmmg struct pointer
 * \param grp     working group pointer
 * \param idx1    index1 value
 * \param idx2    index2 value
 *
 *  Append new values in  group's internal communitor, resizing the buffers if
 *  required
 */
static int n2incAppend( PMMG_pParMesh parmesh, PMMG_pGrp grp, int idx1, int idx2 )
{
  // Adjust the value of scale to reallocate memory more or less agressively
  const float scale = 2.f;
  if ( grp->nitem_int_node_comm + 1 < grp->n2inc_max ) {
    PMMG_RECALLOC(parmesh, grp->node2int_node_comm_index1,
                  scale * grp->n2inc_max, grp->n2inc_max, int,
                  "increasing n2inc_idx1",return PMMG_FAILURE);
    PMMG_RECALLOC(parmesh, grp->node2int_node_comm_index2,
                  scale * grp->n2inc_max, grp->n2inc_max, int,
                  "increasing n2inc_idx1",return PMMG_FAILURE);
    grp->n2inc_max  = grp->n2inc_max * scale;
  }
  grp->node2int_node_comm_index1[ grp->nitem_int_node_comm ] = idx1;
  grp->node2int_node_comm_index2[ grp->nitem_int_node_comm ] = idx2;
  ++grp->nitem_int_node_comm;
  return PMMG_SUCCESS;
}

#warning NIKOS TODO: there is at least code duplication here and in distributegrp, eg the xpoint/xtetra reallocs

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return PMMG_FAILURE
 *         PMMG_SUCCESS
 *
 * if the existing group of only one mesh is too big, split it into into several
 * meshes.
 */
int PMMG_splitGrps( PMMG_pParMesh parmesh )
{
  PMMG_pGrp const grpOld = parmesh->listgrp;
  PMMG_pGrp grpsNew = NULL;
  PMMG_pGrp grpCur = NULL;
  MMG5_pMesh const meshOld = parmesh->listgrp->mesh;
  MMG5_pMesh meshCur = NULL;
  MMG5_pTetra tetraCur = NULL;
  MMG5_pxTetra pxt;
  int *countPerGrp = NULL;
  int ret_val = PMMG_SUCCESS; // returned value (unless set otherwise)
  long long int memMeshTotal = 0;

  idx_t ngrp = 1;
  idx_t *part = NULL;
  int ier, j;

  // counters for tetra, point, while constructing a subgroup
  int tetPerGrp = 0;
  int poiPerGrp = 0;
  int *adja = NULL;
  int vindex = 0;
  int adjidx = 0;

  // Loop counter vars
  int i, grpId, poi, tet, fac;

  assert ( (parmesh->ngrp == 1) && " splitGrps can not split m to n groups");
  printf( "+++++NIKOS+++++[%d/%d]: mesh has: %d(%d) #points, %d(%d) #edges, %d(%d) #tria and %d(%d) tetras(elements)\n",
          parmesh->myrank+1, parmesh->nprocs, meshOld->np, meshOld->npi, meshOld->na, meshOld->nai, meshOld->nt, meshOld->nti, meshOld->ne, meshOld->nei );

  ngrp = howManyGroups( meshOld->ne );
  // Does the group need to be further subdivided to subgroups or not?
  if ( ngrp == 1 )  {
    fprintf( stdout,
            "[%d-%d]: %d group is enough, no need to create sub groups.\n",
            parmesh->myrank+1, parmesh->nprocs, ngrp );
    return PMMG_SUCCESS;
  } else {
    fprintf( stdout,
             "[%d-%d]: %d groups required, splitting into sub groups...\n",
             parmesh->myrank+1, parmesh->nprocs, ngrp );
  }


  // Check if there is enough free memory to allocate the new group
#warning NIKOS TODO: it lacks an estimate of the internal communicator sizes
  if (  parmesh->memCur + 2 * parmesh->listgrp[0].mesh->memCur
      > parmesh->memMax ) {
    fprintf( stderr, "Not enough memory to create listgrp struct\n" );
    return PMMG_FAILURE;
  }

  // use metis to partition the mesh into the computed number of groups needed
  // part array contains the groupID computed by metis for each tetra
  PMMG_CALLOC(parmesh,part,meshOld->ne,idx_t,"metis buffer ", return PMMG_FAILURE);

  if (   PMMG_partition_metis( parmesh, part, ngrp )
      != PMMG_SUCCESS ) {
    ret_val = PMMG_FAILURE;
    goto fail_part;
  }


  /* count_per_grp: count elements per group */
  PMMG_CALLOC(parmesh,countPerGrp,ngrp,int,"counter buffer ",
              ret_val = PMMG_FAILURE;goto fail_part);
  for ( tet = 0; tet < meshOld->ne ; ++tet )
    ++countPerGrp[ part[ tet ] ];
  for ( i = 0; i < ngrp ; i++ )
    printf( "+++++NIKOS+++++[%d/%d]: group[%d] has %d elements\n", parmesh->myrank+1, parmesh->nprocs, i, countPerGrp[i] );



  /* Allocate list of subgroups struct and allocate memory*/
  PMMG_CALLOC(parmesh,grpsNew,ngrp,PMMG_Grp,"subgourp list ",
              ret_val = PMMG_FAILURE; goto fail_counters);
  for ( grpId = 0; grpId < ngrp; ++grpId ) {
    grpCur = &grpsNew[grpId];
    grpCur->mesh = NULL;
    grpCur->met  = NULL;
    grpCur->disp = NULL;

    MMG3D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &grpCur->mesh,
                     MMG5_ARG_ppMet, &grpCur->met, MMG5_ARG_end );

    meshCur = grpCur->mesh;

    /* Copy the mesh filenames */
    if ( !MMG5_Set_inputMeshName( meshCur, meshOld->namein ) ) {
      ret_val = PMMG_FAILURE;
      goto fail_sgrp;
    }
#warning NIKOS: Why are these commented out? Enable or delete?
    //if ( !MMG5_Set_inputSolName( meshCur, grpCur->met, grpOld->met->namein ) ) {
    //  ret_val = PMMG_FAILURE;
    //  goto fail_sgrp;
    //}
    if ( !MMG5_Set_outputMeshName( meshCur, meshOld->nameout ) ) {
      ret_val = PMMG_FAILURE;
      goto fail_sgrp;
    }
    //if ( !MMG5_Set_outputSolName( meshCur, grpCur->met, grpOld->met->nameout ) ) {
    //  ret_val = PMMG_FAILURE;
    //  goto fail_sgrp;
    //}

    MMG3D_Set_meshSize( grpCur->mesh, countPerGrp[grpId], countPerGrp[grpId],
                        0, 0, 0, 0 );

    /* Copy the info structure of the initial mesh: it contains the remeshing
     * options */
    memcpy(&(grpCur->mesh->info),&(meshOld->info),sizeof(MMG5_Info) );

    meshCur->xtmax = 128;
    PMMG_CALLOC(meshCur,meshCur->xtetra,meshCur->xtmax+1,MMG5_xTetra,
                "msh boundary xtetra", ret_val = PMMG_FAILURE;goto fail_sgrp);

    /* memory to store normals for boundary points */
    meshCur->xpmax  = 128;
    PMMG_CALLOC(meshCur,meshCur->xpoint,meshCur->xpmax+1,MMG5_xPoint,
                "boundary points", ret_val = PMMG_FAILURE;goto fail_sgrp);

    PMMG_CALLOC(meshCur,meshCur->adja,4*meshCur->nemax+5,int,"adjacency table",
                ret_val = PMMG_FAILURE;goto fail_sgrp);

    grpCur->n2inc_max = 128;
    PMMG_CALLOC(parmesh,grpCur->node2int_node_comm_index1,grpCur->n2inc_max,int,
                "subgroup internal1 communicator ", ret_val = PMMG_FAILURE;goto fail_sgrp);
    PMMG_CALLOC(parmesh,grpCur->node2int_node_comm_index2,grpCur->n2inc_max,int,
               "subgroup internal2 communicator ", ret_val = PMMG_FAILURE;goto fail_sgrp);
    printf( "+++++NIKOS+++++[%d/%d]:\t meshCur %p,\t xtmax %d - xtetra:%p,\t xpmax %d - xpoint %p,\t nemax %d - adja %p,\t ne %d- index1 %p index2 %p \n",
            parmesh->myrank + 1, parmesh->nprocs,
            meshCur, meshCur->xtmax, meshCur->xtetra,
            meshCur->xpmax, meshCur->xpoint,
            meshCur->nemax, meshCur->adja,
            grpCur->n2inc_max, grpCur->node2int_node_comm_index1, grpCur->node2int_node_comm_index2);
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
  // use point[].tmp field to "remember" index in internal communicator of
  // vertices. specifically:
  //   place a copy of vertices' node2index2 position at point[].tmp field
  //   or -1 if they are not in the comm
  for ( poi = 1; poi < meshOld->np + 1; ++poi )
    meshOld->point[poi].tmp = -1;
  for ( i = 0; i < grpOld->nitem_int_node_comm; i++ )
    meshOld->point[ grpOld->node2int_node_comm_index1[ i ] ].tmp = grpOld->node2int_node_comm_index2[ i ];
  for ( grpId = 0 ; grpId < ngrp ; ++grpId ) {
    grpCur = &grpsNew[grpId];
    meshCur = grpCur->mesh;

    // use point[].flag field to "remember" assigned local(in subgroup) numbering
    for ( poi = 1; poi < meshOld->np + 1; ++poi )
      meshOld->point[poi].flag = 0;

    // Loop over tetras and choose the ones to add in the submesh being constructed
    tetPerGrp = 0;
    poiPerGrp = 0;
    for ( tet = 1; tet < meshOld->ne + 1; ++tet ) {
      // MMG3D_Tetra.flag is used to update adjacency vector:
      //   if the tetra belongs to the group we store the local tetrahedron id
      //   in the tet.flag
      meshOld->tetra[tet].flag = 0;

      // Skip elements that do not belong in the group processed in this iteration
      if ( grpId != part[ tet - 1 ] )
        continue;

      ++tetPerGrp;
      assert( ( tetPerGrp < meshCur->nemax ) && "overflowing tetra array?" );
      tetraCur = meshCur->tetra + tetPerGrp;
      meshOld->tetra[tet].flag = tetPerGrp;

      // add tetrahedron to subgroup (copy from original group)
      memcpy( tetraCur, &meshOld->tetra[tet], sizeof(MMG5_Tetra) );
      tetraCur->base = 0;
      tetraCur->mark = 0;
      tetraCur->flag = 0;
#warning NIKOS TODO: there are also cases where we have to update existing boundaries, eg when boundary tetras are also interfaced with new groups
      // xTetra: this element was already an xtetra (in meshOld)
      if ( tetraCur->xt != 0 ) {
        if ( PMMG_SUCCESS != xtetraAppend( meshCur, meshOld, tet ) ) {
          ret_val = PMMG_FAILURE;
          goto fail_sgrp;
        }
        tetraCur->xt = meshCur->xt;
      }
#warning NIKOS TODO: are there or should these be just deleted?
      //else
      //  printf ("there may be more: new boundary between subgroups\n");

      // Add tetrahedron vertices in points struct and
      // adjust tetrahedron vertices indices
      for ( poi = 0; poi < 4 ; ++poi ) {
        if ( 0 == meshOld->point[ meshOld->tetra[tet].v[poi] ].flag ) { // 1st time that this point is seen in this subgroup

          // Add point in subgroup point array
          ++poiPerGrp;
          assert( (poiPerGrp < meshCur->npmax) && "overflowing mesh points" );
          memcpy( meshCur->point + poiPerGrp,
                  &meshOld->point[ meshOld->tetra[tet].v[poi] ],
                  sizeof(MMG5_Point) );

          // update tetra vertex reference
          tetraCur->v[poi] = poiPerGrp;

          // "Remember" the assigned subgroup point id
          meshOld->point[ meshOld->tetra[tet].v[poi] ].flag = poiPerGrp;

          // Add point in subgroup's communicator if it already was in group's
          // communicator
          if ( meshCur->point[ poiPerGrp ].tmp != -1 ) {
            if (  n2incAppend( parmesh, grpCur, poiPerGrp,
                               meshCur->point[ poiPerGrp ].tmp )
                != PMMG_SUCCESS ) {
              ret_val = PMMG_FAILURE;
              goto fail_sgrp;
            }
          }
          // xPoints: this was already a boundary point
          if ( meshCur->point[poiPerGrp].xp != 0 ) {
            if ( PMMG_SUCCESS != xpointAppend( meshCur, meshOld, tet, poi ) ) {
              ret_val = PMMG_FAILURE;
              goto fail_sgrp;
            }
            meshCur->point[poiPerGrp].xp = meshCur->xp;
          }
#warning NIKOS: are there other border cases or simply delete these?
          //else
          //  printf( " some other more border cases, think about it \n" );
        } else { // point is already included in this subgroup, update current tetra vertex reference
          tetraCur->v[poi] = meshOld->point[ meshOld->tetra[tet].v[poi] ].flag;
        }
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

            /* creation of the interface faces : ref 0 and tag MG_PARBDY */
            if( !meshCur->tetra[tet].xt ) {
              if ( PMMG_SUCCESS != xtetraAppend( meshCur, meshOld, tet ) ) {
                ret_val = PMMG_FAILURE;
                goto fail_sgrp;
              }
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
            // Now that we know the local gr_idx for both tetra we should also
            // update the adjacency entry for the adjacent face
            assert(  (4 * (meshOld->tetra[adjidx].flag-1) + 1 + vindex ) < (4 * meshCur->nemax+5)
                &&"adja overflow" );
            meshCur->adja[ 4 * (meshOld->tetra[adjidx].flag-1) + 1 + vindex ] =
              4 * tetPerGrp + fac;
          }
        }
      }
      adja = &meshOld->adja[ 4 * ( tet - 1 ) + 1 ];
      for ( fac = 0; fac < 4; ++fac ) {
        adjidx = adja[ fac ] / 4;

        if ( adjidx && grpId != part[ adjidx - 1 ] ) {
          for ( poi = 0; poi < 3; ++poi ) {
            if ( meshCur->point[ tetraCur->v[ _MMG5_idir[fac][poi] ] ].tmp == -1 ) {
              assert ( grpCur->nitem_int_node_comm < grpCur->n2inc_max && "if increasing the internal comm size is correct, add call to n2incAppend here");
              grpCur->node2int_node_comm_index1[ grpCur->nitem_int_node_comm ] =
                tetraCur->v[ _MMG5_idir[fac][poi] ];
              grpCur->node2int_node_comm_index2[ grpCur->nitem_int_node_comm ] =
                parmesh->int_node_comm->nitem + 1;
              meshOld->point[ meshOld->tetra[tet].v[ _MMG5_idir[fac][poi]  ] ].tmp =
                parmesh->int_node_comm->nitem + 1;
              meshCur->point[ tetraCur->v[ _MMG5_idir[fac][poi] ] ].tmp =
                parmesh->int_node_comm->nitem + 1;
              ++grpCur->nitem_int_node_comm;
              ++parmesh->int_node_comm->nitem;
            }
          }
        }
      }
    }
    assert( (meshCur->ne == tetPerGrp) && "Error in PMMG_splitGrps" );
    printf( "+++++NIKOS[%d/%d]:: %d points in group, %d tetra (expected: %d)ed.%d nitem in int communicator\n",
            ngrp, grpId+1, poiPerGrp, tetPerGrp, meshCur->ne, grpCur->nitem_int_node_comm );
  }



//NIKOS TODO  for ( grpId = 0 ; grpId < ngrp ; grpId++ ) {
//NIKOS TODO    meshCur = grpsNew[grpId].mesh;
//NIKOS TODO    Is there an additional call to MMG missing (to compact meshes or sth?)
//NIKOS TODO    //MMG5_saveMshMesh( meshCur, mesMMG5_pSol met,const char *filename)
//NIKOS TODO    char name[20];
//NIKOS TODO    sprintf( name, "mesh-p%d-%02d.mesh", parmesh->myrank+1, grpId );
//NIKOS TODO    _MMG3D_bdryBuild(meshCur); //note: no error checking
//NIKOS TODO    MMG3D_saveMesh( meshCur, name );
//NIKOS TODO  }

#warning NIKOS: DEALLOCATE OLD STRUCTS!!!!! MEM LEAKING HERE
  parmesh->listgrp = grpsNew;
  parmesh->ngrp = ngrp;

  PMMG_PMesh_SetMemMax(parmesh, 50);

#ifndef NDEBUG
  if ( PMMG_checkIntComm ( parmesh ) ) {
    fprintf ( stderr, " INTERNAL COMMUNICATOR CHECK FAILED \n" );
    ret_val = PMMG_FAILURE;
    goto fail_sgrp;
  }
#endif

  //MMG3D_Free_all( MMG5_ARG_start,
  //                MMG5_ARG_ppMesh, &(parmesh->listgrp->mesh), MMG5_ARG_ppMet, &(parmesh->listgrp->met),
  //                MMG5_ARG_end);
  //parmesh->listgrp = grpsNew;
  //parmesh->ngrp = ngrp;

  // No error so far, skip deallocation of lstgrps
  goto fail_counters;

  // fail_sgrp deallocates any mesh that has been allocated in listgroup.
  // Should be executed only if an error has occured
fail_sgrp:
  for ( grpId = 0; grpId < ngrp; ++grpId ) {
    meshCur = grpsNew[grpId].mesh;
    if ( grpCur->node2int_node_comm_index2 != NULL )
      PMMG_DEL_MEM(parmesh,grpCur->node2int_node_comm_index2,meshCur->ne/3,int,"subgroup internal2 communicator ");
    if ( grpCur->node2int_node_comm_index1 != NULL )
      PMMG_DEL_MEM(parmesh,grpCur->node2int_node_comm_index1,meshCur->ne/3,int,"subgroup internal1 communicator ");
    if ( meshCur->adja != NULL )
      PMMG_DEL_MEM(meshCur,meshCur->adja,4*meshCur->nemax+5,int,"adjacency table");
    if ( meshCur->xpoint != NULL )
      PMMG_DEL_MEM(meshCur,meshCur->xpoint,meshCur->xpmax+1,MMG5_xPoint,"boundary points");
    if ( meshCur->xtetra != NULL )
      PMMG_DEL_MEM(meshCur,meshCur->xtetra,meshCur->xtmax+1,MMG5_xTetra,"msh boundary tetra");
#warning NIKOS: ADD DEALLOC/WHATEVER FOR EACH MESH:    MMG3D_DeInit_mesh() or STH
  }
  PMMG_DEL_MEM(parmesh,countPerGrp,ngrp,int,"counter buffer ");
  // these labels should be executed as part of normal code execution before
  // returning as well as error handling
fail_counters:
  PMMG_DEL_MEM(parmesh,countPerGrp,ngrp,int,"counter buffer ");
fail_part:
  PMMG_DEL_MEM(parmesh,part,meshOld->ne,idx_t,"free metis buffer ");
  return ret_val;
}

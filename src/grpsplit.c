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
#include "libparmmg.h" // PMMG_mesh2metis
#include "grpsplit.h"

// Subgroups target size. It is chosen arbitrarily to help assist the remesher work faster
static const int REMESHER_TARGET_MESH_SIZE = 2000;

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
  PMMG_pGrp grpsOld = parmesh->listgrp;
  MMG5_pMesh mesh = grpsOld->mesh;
  PMMG_pGrp grpsNew = NULL;
  int *countPerGrp = NULL;

  idx_t ngrp = 1;
  idx_t *part = NULL;
  idx_t *xadj = NULL;
  idx_t *adjncy = NULL;
//NIKOS TODO: experiment with the number of balancing constraints
  idx_t ncon = 1;/*number of balancing constraint*/
  idx_t nelt = mesh->ne;
  idx_t objval;
  int ier;

  int tetPerGrp = 0;
  int poiPerGrp = 0;
int removeMe = 0; //!!!!!NIKOS this is used only to keep the last iteration number to print in debugging message.2be deleted
int removeMe2 = 0; //!!!!!NIKOS this is used only to keep the number of already existing xtetras to print in debugging message.2be deleted

  // Loop counter vars
  int i, grpId, poi, tet;

printf( "+++++NIKOS+++++[%d/%d]: mesh has: %d(%d) #points, %d(%d) #edges, %d(%d) #tria and %d(%d) tetras(elements).",
        parmesh->myrank+1, parmesh->nprocs, mesh->np, mesh->npi, mesh->na, mesh->nai, mesh->nt, mesh->nti, mesh->ne, mesh->nei );

  ngrp = HowManyGroups( mesh->ne );
  /* Check whether we need to further subdivide the groups or not */
  if ( ngrp == 1 )  {
    printf( "[%d-%d]: %d subgroup is needed, no sub groups to create.\n", parmesh->myrank+1, parmesh->nprocs, ngrp );
    return 1;
  } else {
    printf( "[%d-%d]: %d subgroups are needed, creating sub groups...\n", parmesh->myrank+1, parmesh->nprocs, ngrp );
  }

  // use metis to partition the mesh into the computed number of groups needed

  // part array contains each tetras groupID
  _MMG5_SAFE_CALLOC( part, mesh->ne, idx_t );

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
  for ( i = 0; i < mesh->ne ; i++ )
    countPerGrp[part[i]]++;
for ( i = 0; i < ngrp ; i++ )
  printf( "+++++NIKOS+++++[%d/%d]: group[%d] has %d elements\n", parmesh->myrank+1, parmesh->nprocs, i, countPerGrp[i] );

  /* Create list of subgroups struct */

  _MMG5_SAFE_CALLOC( grpsNew, ngrp, PMMG_Grp );

//NIKOS TODO: do MMG3D_Init_mesh/MMG3D_Set_meshSize offer sth? otherwise allocate on your own
  for ( i = 0; i < ngrp; ++i ) {
    PMMG_pGrp grp = &grpsNew[i];
    grp->mesh = NULL;
    grp->sol  = NULL;
    grp->disp = NULL;
    //_MMG5_SAFE_CALLOC( ptr, size, type );
    //_MMG5_SAFE_REALLOC( ptr, size, type, message );
    MMG3D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &grp->mesh,
                     MMG5_ARG_ppMet, &grp->sol, MMG5_ARG_end );
    MMG3D_Set_meshSize( grp->mesh, countPerGrp[i], countPerGrp[i], 0, 0, 0, 0 );

    grp->mesh->ne = countPerGrp[i];
printf( "+++++NIKOS+++++[%d/%d]: lstgrp[%d] @ %p, ", parmesh->myrank+1, parmesh->nprocs, i, grp );
printf( "np = %d/%d, na = %d/%d, nt = %d/%d, ne = %d/%d \n",
        grp->mesh->np, grp->mesh->npmax,
        grp->mesh->na, grp->mesh->namax,
        grp->mesh->nt, grp->mesh->ntmax,
        grp->mesh->ne, grp->mesh->nemax );
printf( "\n" );
  }

//NIKOS TODO: LOOP OVER part ngrp TIMES or USE A tmp[NGROUPS][NP] ARRAY AND LOOP ONLY ONCE? it wastes memory (eg 10 groups x 100k tetra = 4Mb of ints) but only loops over part once
  for( grpId = 0 ; grpId < ngrp ; grpId++ ) {
    MMG5_pMesh curMesh = grpsNew[grpId].mesh;

    // Will use the MMG5_Point.flag field to assign local numbering in the newly created subgroups.
    for ( poi = 1; poi < mesh->np + 1; poi++ )
      mesh->point[poi].flag = 0;

    for ( tet = 1; tet < mesh->ne + 1; tet++ ) {
      MMG5_pTetra curGrpTetra = NULL;

      /* Skip elements that are not included in the group being processed */
      if ( grpId != part[tet-1] )
        continue;

      removeMe = tet; //!!!!!NIKOS TODO delete

      ++tetPerGrp;
      curGrpTetra = curMesh->tetra + tetPerGrp;

      /* Add tetrahedron to group */
if ( tetPerGrp == 1 )
  printf( "+++++NIKOS[%d/%d]:: Adding tetra: from %d to %d\n", grpId+1, ngrp, tet, tetPerGrp );
  //printf( "+++++NIKOS[%d/%d]:: copied %zd bytes from %p to %p\n", grpId+1, ngrp, sizeof(MMG5_Tetra), &mesh->tetra[tet], curGrpTetra );
      memcpy( curGrpTetra, &mesh->tetra[tet], sizeof(MMG5_Tetra) );

      /* Add tetrahedron vertices in points and adjust tetrahedron vertices indices */
      for ( poi = 0; poi < 4 ; poi++ ) {
        // point hasnt been assigned a group local point id
        if ( 0 == mesh->point[ mesh->tetra[tet].v[poi] ].flag ) {

          // Add point in group mesh point array
          ++poiPerGrp;
//printf( "+++++NIKOS[%d/%d]:: copied %zd point bytes from %p to %p\n", grpId+1, ngrp, sizeof(MMG5_Point), &mesh->point[mesh->tetra[tet].v[poi]], curMesh->point + poiPerGrp );
          memcpy( curMesh->point + poiPerGrp, &mesh->point[ mesh->tetra[tet].v[poi] ], sizeof(MMG5_Point) );

          // "Remember" group local point id
          mesh->point[ mesh->tetra[tet].v[poi] ].flag = poiPerGrp;

          // update group mesh tetra point reference
          curGrpTetra->v[poi] = poiPerGrp;
//printf( "+++++NIKOS[%d/%d]:: point %d has xp %d\n", grpId+1, ngrp, mesh->tetra[tet].v[poi], mesh->point[mesh->tetra[tet].v[poi]].xp );
if ( mesh->point[mesh->tetra[tet].v[poi]].xp ) //NEW POINT WAS ALREADY POINTING TO AN xp.
  ++removeMe2;
       // point has already been assigned a group local point id
       } else {
         // update group mesh tetra point reference
         curGrpTetra->v[poi] = mesh->point[mesh->tetra[tet].v[poi]].flag;
       }
    }

    // Take care of the xtetra/xp now
    //  xPoints/xTetras:
    //  => when mesh is split there are new boundary tetras on top of the already existing ones
    //  => xPoints/xTetras that were on the boundary before will continue to be so
    //  => Plus there will be new ones as well. which is the criterion? perhaps we can look for them after the groups are created?
    for ( i = 0 ; i < 4; i++ ) { LOOP OVER THE adja VALUES
      if ( globalmesh->adja[4*( current global element number - 1 ) + 1 + i] ) 
        ;
      else
    //printf( "+++++NIKOS[%d/%d]:: Current tetra id: %4d, xt value: %4d ", grpId+1, ngrp, tetPerGrp, curGrpTetra->xt );
    //for ( int poi = 0; poi < 4 ; poi++ )
    //  printf( "poi %1d(%3d) - xp %3d ", poi, curGrpTetra->v[poi], curMesh->point[curGrpTetra->v[poi]].xp );
    //printf( "\n" );
    //      if ( mesh->tetra[tet].xt ) {
    //++removeMe2;
    //printf( "+++++NIKOS[%d/%d]:: tetra %d points to non zero xt %d\n", grpId+1, ngrp, tet, mesh->tetra[tet].xt );
    //        printf( "+++++NIKOS[%d/%d]:: tetra %d(%d) points to xt: %d\n", grpId+1, ngrp, tet, tetPerGrp,mesh->tetra[tet].xt );
    printf( "+++++NIKOS[%d/%d]:: Adding tetra: from %d to %d\n", grpId+1, ngrp, removeMe, tetPerGrp );
    printf( "+++++NIKOS[%d/%d]:: added %d points(%d xt ref found) in group with %d tetra out of %d tetras expected.\n", grpId+1, ngrp, poiPerGrp, removeMe2, tetPerGrp, curMesh->ne );
    //printf( "+++++NIKOS[%d/%d]:: added %d points in group with %d tetra out of %d tetras expected.inherited %d xtetras\n", grpId+1, ngrp, poiPerGrp, tetPerGrp, curMesh->ne, removeMe2 );
    //}
    }
  }
//#error NIKOS: CHANGE THE MEMORY ALLOCATIONS WITH PROPER ALLOCATION+REALLOCATION
//NIKOS TODO: FREE ME PROPERLY!
  //parmesh->ngrp = ngrp;
  //free(grpsOld);
  //grpsOld = grpsNew;

  _MMG5_SAFE_FREE( xadj );
  _MMG5_SAFE_FREE( adjncy );
  _MMG5_SAFE_FREE( part );
  _MMG5_SAFE_FREE( countPerGrp );
  return( 1 );
}

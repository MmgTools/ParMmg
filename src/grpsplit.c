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

// Subgroups target size. It is chosen arbitrarily to help assist the remesher
// work faster
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

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Split one group into several groups.
 *
 */
int PMMG_splitGrps( PMMG_pParMesh parmesh )
{
	PMMG_pGrp grp = parmesh->listgrp;
	MMG5_pMesh mesh = grp->mesh;

	printf( "+++++NIKOS+++++[%d/%d]: I have %d group(s) and the first lstgroup is at %p.",
			parmesh->myrank+1, parmesh->nprocs, parmesh->ngrp, parmesh->listgrp );
	printf( "\n\tMesh has: %d(%d) #points, %d(%d) #edges, %d(%d) #tria and %d(%d) tetras(elements).",
			mesh->np, mesh->npi, mesh->na, mesh->nai, mesh->nt, mesh->nti, mesh->ne, mesh->nei );

	idx_t ngrp = HowManyGroups( mesh->ne );
	printf( " Calculated that %d groups are needed.", ngrp );
	/* Check whether we need to further subdivide the groups or not */
	if ( ngrp == 1 )  {
		printf( " Nothing to do, returning.\n" );
		return 1;
	} else {
		printf( " Things to do, continueing: \n" );
	}

	/* Call metis*/
	idx_t *part = NULL;
	/** Call metis for partionning*/
//NIKOS TODO: there is some code duplication. eg this metis calling is also in the metisfunction.c
	_MMG5_SAFE_CALLOC( part, mesh->ne, idx_t );
{
	idx_t *xadj = NULL;
	idx_t *adjncy = NULL;
//NIKOS TODO: experiment with the number of balancing constraints
	idx_t ncon = 1;/*number of balancing constraint*/
	idx_t nelt = mesh->ne;
	idx_t objval;

	PMMG_mesh2metis(parmesh, &xadj, &adjncy);

	int ier =  METIS_PartGraphKway( &nelt          , &ncon          , xadj , adjncy        , NULL/*vwgt*/ ,
	                                NULL/*vsize*/  , NULL/*adjwgt*/ , &ngrp, NULL/*tpwgts*/, NULL/*ubvec*/,
	                                NULL/*options*/, &objval, part );
	//int ier =  METIS_PartGraphRecursive( &nelt, &ncon, xadj, adjncy, NULL/\*vwgt*\/, NULL/\*vsize*\/,
	//                                     NULL/\*adjwgt*\/, &nproc, NULL/\*tpwgts*\/,
	//                                     NULL/\*ubvec*\/, NULL/\*options*\/, &objval, part );
	_MMG5_SAFE_FREE(xadj);
	_MMG5_SAFE_FREE(adjncy);
//NIKOS TODO: better error handling, esp memory deallocation
	if ( ier != METIS_OK )
		return (0);
}

	/* count_per_grp: count elements per group */
	int *countPerGrp = NULL;
	_MMG5_SAFE_CALLOC( countPerGrp, ngrp, int );
	for ( int i = 0; i < ngrp ; i++ )
		countPerGrp[i] = 0;
	for ( int i = 0; i < mesh->ne ; i++ )
		countPerGrp[part[i]]++;
	for ( int i = 0; i < ngrp ; i++ )
		printf( "+++++NIKOS+++++[%d/%d]: group[%d] has %d elements\n", parmesh->myrank+1, parmesh->nprocs, i, countPerGrp[i] );

	/* Create subgroups: the mesh contained in each group and the internal communicators */
	PMMG_pGrp listgrp = NULL;
	_MMG5_SAFE_CALLOC( listgrp, ngrp, PMMG_Grp );

	for ( int i = 0; i < ngrp; ++i ) {
		PMMG_pGrp grp = &listgrp[i];
		grp->mesh = NULL;
		grp->sol  = NULL;
		grp->disp = NULL;
		MMG3D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &grp->mesh,
		                 MMG5_ARG_ppMet, &grp->sol, MMG5_ARG_end );
//NIKOS TODO: allcate on your own, no need to use MMG3D_Set_meshSize
		//_MMG5_SAFE_CALLOC( ptr, size, type );
		//_MMG5_SAFE_REALLOC( ptr, size, type, message );
		MMG3D_Set_meshSize( grp->mesh, countPerGrp[i], countPerGrp[i], 0, 0, 0, 0 );

		grp->mesh->ne = countPerGrp[i];
	}

	for ( int i = 0; i < ngrp ; i++ ) {
		printf( " lstgrp[%d] @ %p, ", i, listgrp + i );
		printf( "np = %d/%d, na = %d/%d, nt = %d/%d, ne = %d/%d \n",
		        listgrp[i].mesh->np, listgrp[i].mesh->npmax,
		        listgrp[i].mesh->na, listgrp[i].mesh->namax,
		        listgrp[i].mesh->nt, listgrp[i].mesh->ntmax,
		        listgrp[i].mesh->ne, listgrp[i].mesh->nemax );
		printf( "\n" );
	}

	/* Loop over tetras and copy them to the groups metis assigned them */
//NIKOS TODO: LOOP OVER part ngrp TIMES or USE A tmp[NGROUPS][NP] ARRAY AND LOOP ONLY ONCE? it wastes memory (eg 10 groups x 100k tetra = 4Mb of ints) but only loops over part once
	for( int grpId = 0 ; grpId < ngrp ; grpId++ ) {

//NIKOS TODO: is this redundant given the mesh struct has just been allocated?
		/* Will use the field flag to assign local numbering in the newly created subgroups. Initialize to 0 */
		//NIKOS PREPEI NA MIDENISEIS GIA NA TO XRISIMOPOIHSEIS OS flag FIELD
		for ( int poi = 0; poi < mesh->np ; poi++ )
			mesh->point[poi].flag = 0;

		MMG5_pMesh curMesh = listgrp[grpId].mesh;
		int tetPerGrp = 0;
		int poiPerGrp = 0;
		for ( int tet = 0; tet < mesh->ne ; tet++ ) {

			/* Skip elements that are not included in the group being processed */
			if ( grpId != part[tet] )
				continue;

			++tetPerGrp;
			MMG5_pTetra curGrpTetra = curMesh->tetra + tetPerGrp;

			/* Add tetrahedron to group */
			//printf( "+++++NIKOS[%d/%d]:: copied %zd bytes from %p to %p\n", grpId+1, ngrp, sizeof(MMG5_Tetra), &mesh->tetra[tet], curGrpTetra );
			printf( "+++++NIKOS[%d/%d]:: Adding tetra: from %d to %d\n", grpId+1, ngrp, tet, tetPerGrp );
			memcpy( curGrpTetra, &mesh->tetra[tet], sizeof(MMG5_Tetra) );

			/* Add tetrahedron vertices in points and adjust tetrahedron vertices indices */
			for ( int poi = 0; poi < 4 ; poi++ ) {
				// point hasnt been assigned a group local point id
				if ( 0 == mesh->point[mesh->tetra[tet].v[poi]].flag ) {
					// Add point in group mesh point array
					++poiPerGrp;

					//printf( "+++++NIKOS[%d/%d]:: copied %zd point bytes from %p to %p\n", grpId+1, ngrp, sizeof(MMG5_Point), &mesh->point[mesh->tetra[tet].v[poi]], curMesh->point + poiPerGrp );
					memcpy( curMesh->point + poiPerGrp, &mesh->point[mesh->tetra[tet].v[poi]], sizeof(MMG5_Point) );

					// "Remember" group local point id
					mesh->point[mesh->tetra[tet].v[poi]].flag = poiPerGrp;

					// update group mesh tetra point reference
					curGrpTetra->v[poi] = poiPerGrp;
				// point has already been assigned a group local point id
				} else {
					// update group mesh tetra point reference
					curGrpTetra->v[poi] = mesh->point[mesh->tetra[tet].v[poi]].flag;
				}
			}

			// Take care of the xtetra/xp now
			//printf( "+++++NIKOS[%d/%d]:: Current tetra id: %4d, xt value: %4d ", grpId+1, ngrp, tetPerGrp, curGrpTetra->xt );
			//for ( int poi = 0; poi < 4 ; poi++ )
			//	printf( "poi %1d(%3d) - xp %3d ", poi, curGrpTetra->v[poi], curMesh->point[curGrpTetra->v[poi]].xp );
			//printf( "\n" );
		}
		printf( "+++++NIKOS[%d/%d]:: added %d points in group with %d tetra out of %d tetras expected\n", grpId+1, ngrp, poiPerGrp, tetPerGrp, curMesh->ne );
	}

//NIKOS TODO: FREE ME PROPERLY!
	//parmesh->ngrp = ngrp;
	//free(parmesh->listgrp);
	//parmesh->listgrp = listgrp;

	_MMG5_SAFE_FREE( countPerGrp );
	_MMG5_SAFE_FREE( part );
	return( 1 );
}

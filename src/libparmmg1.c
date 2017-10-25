/**
 * \file libparmmg1.c
 * \brief Parallel remeshing library
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version
 * \copyright
 *
 */
#include "parmmg.h"

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return PMMG_FAILURE
 *         PMMG_SUCCESS
 *
 * Pack the sparse meshes of each group and create triangles and edges before
 * getting out of library
 *
 */
int PMMG_packParMesh( PMMG_pParMesh parmesh )
{
  PMMG_pGrp   grp;
  MMG5_pMesh  mesh;
  MMG5_pSol   met;
  MMG5_pSol   disp;
  MMG5_pPoint ppt,pptnew;
  int         *node2int_node_comm_index1;
  int         np,nbl,nc,k,igrp;
  int         iadr;

  for ( igrp=0; igrp<parmesh->ngrp; ++igrp ) {
    grp                       = &parmesh->listgrp[igrp];
    mesh                      = grp->mesh;
    met                       = grp->met;
    node2int_node_comm_index1 = grp->node2int_node_comm_index1;
    disp                      = grp->disp;

    /* compact vertices */
    if ( !MMG3D_count_packedPoints(mesh,&np,&nc) ) return 0;

    /* node index update in internal communicator */
    for (k=0; k<grp->nitem_int_node_comm; ++k) {
      iadr = node2int_node_comm_index1[k];

      ppt = &mesh->point[iadr];
      assert ( MG_VOK(ppt) );

      node2int_node_comm_index1[k] = ppt->tmp;
    }

    /* compact tetrahedra */
    if ( mesh->adja ) {
      if ( !MMG3D_pack_tetraAndAdja(mesh) ) return 0;
    }
    else {
      if ( !MMG3D_pack_tetra(mesh) ) return 0;
    }

    /* update prisms and quads vertex indices */
    if ( !MMG3D_pack_prismsAndQuads(mesh) ) return 0;

    /* compact metric */
    if ( met && met->m )
      if ( !MMG3D_pack_sol(mesh,met) ) return 0;

    /* compact displacement */
    if ( disp && disp->m )
      if ( !MMG3D_pack_sol(mesh,disp) ) return 0;

    /* compact vertices */
    // pack_points fail because we miss some xp after merging
    // if ( !MMG3D_pack_points(mesh) ) return 0;
    np  = 0;
    nbl = 1;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )
        continue;
      np++;
      if ( k!=nbl ) {
        pptnew = &mesh->point[nbl];
        memmove(pptnew,ppt,sizeof(MMG5_Point));
        memset(ppt,0,sizeof(MMG5_Point));
        ppt->tag = MG_NUL;
      }
      nbl++;
    }
    mesh->np = np;

    /* create prism adjacency */
    if ( !MMG3D_hashPrism(mesh) ) {
      fprintf(stderr,"\n  ## Error: %s: prism hashing problem. Exit program.\n",
              __func__);
      return(0);
    }

    /* Remove the MG_REQ tags added by the nosurf option */
    MMG3D_unset_reqBoundaries(mesh);

    if ( mesh->info.imprim ) {
      fprintf(stdout,"     NUMBER OF VERTICES   %8d   CORNERS %8d\n",mesh->np,nc);
      fprintf(stdout,"     NUMBER OF ELEMENTS   %8d\n",mesh->ne);
    }

    /* to could save the mesh, the adjacency have to be correct */
    if ( mesh->info.ddebug && (!_MMG5_chkmsh(mesh,1,1) ) ) {
      fprintf(stderr,"  ##  Problem. Invalid mesh.\n");
      return PMMG_FAILURE;
    }
  }

  return PMMG_SUCCESS;
}


/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return PMMG_SUCCESS
 *         PMMG_FAILURE
 *
 * Check the validity of the input mesh data (tetra orientation, solution
 * compatibility with respect to the provided mesh, Mmg options).
 *
 */
#warning NIKOS: more appropriate name: PMMG_CHECK_MESHES
int PMMG_check_inputData(PMMG_pParMesh parmesh)
{
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  int        k;

  if ( !parmesh->myrank && parmesh->listgrp[0].mesh->info.imprim )
    fprintf(stdout,"\n  -- PMMG: CHECK INPUT DATA\n");

  for ( k=0; k<parmesh->ngrp; ++k ) {
    mesh = parmesh->listgrp[k].mesh;
    met  = parmesh->listgrp[k].met;

    /* Check options */
    if ( mesh->info.lag > -1 ) {
      fprintf(stderr,
              "  ## Error: lagrangian mode unavailable (MMG3D_IPARAM_lag):\n");
      return PMMG_FAILURE;
    } else if ( mesh->info.iso ) {
      fprintf(stderr,"  ## Error: level-set discretisation unavailable"
              " (MMG3D_IPARAM_iso):\n");
      return PMMG_FAILURE;
    } else if ( mesh->info.optimLES && met->size==6 ) {
      fprintf(stdout,"  ## Error: strong mesh optimization for LES methods"
              " unavailable (MMG3D_IPARAM_optimLES) with an anisotropic metric.\n");
      return PMMG_FAILURE;
    }
    /* specific meshing */
    if ( met->np ) {
      if ( mesh->info.optim ) {
        printf("\n  ## ERROR: MISMATCH OPTIONS: OPTIM OPTION CAN NOT BE USED"
               " WITH AN INPUT METRIC.\n");
        return PMMG_FAILURE;
      }

      if ( mesh->info.hsiz>0. ) {
        printf("\n  ## ERROR: MISMATCH OPTIONS: HSIZ OPTION CAN NOT BE USED"
               " WITH AN INPUT METRIC.\n");
        return PMMG_FAILURE;
      }
    }

    if ( mesh->info.optim &&  mesh->info.hsiz>0. ) {
      printf("\n  ## ERROR: MISMATCH OPTIONS: HSIZ AND OPTIM OPTIONS CAN NOT BE USED"
             " TOGETHER.\n");
      return PMMG_FAILURE;
    }

    /* load data */
    _MMG5_warnOrientation(mesh);

    if ( met->np && (met->np != mesh->np) ) {
      fprintf(stdout,"  ## WARNING: WRONG METRIC NUMBER. IGNORED\n");
      _MMG5_DEL_MEM(mesh,met->m,(met->size*(met->npmax+1))*sizeof(double));
      met->np = 0;
    } else if ( met->size!=1 && met->size!=6 ) {
      fprintf(stderr,"  ## ERROR: WRONG DATA TYPE.\n");
      return PMMG_FAILURE;
    }
  }
  if ( !parmesh->myrank && parmesh->listgrp[0].mesh->info.imprim )
    fprintf(stdout,"  -- CHECK INPUT DATA COMPLETED\n");

  return PMMG_SUCCESS;
}

/**
 * \param parmesh pointer toward a parmesh structure where the boundary entities
 * are stored into xtetra and xpoint strucutres
 *
 * Main program of the parallel remeshing library: split the meshes over each
 * proc into groups, then perform niter of sequential remeshing of each group
 * (with moving of the proc boundaries between two iterations) and last, merge
 * the groups over each proc.
 *
 * \return PMMG_STRONGFAILURE if  we can't save the mesh (non-conform),
 *         PMMG_LOWFAILURE    if  we can save the mesh
 *         PMMG_SUCCESS
 */
int PMMG_parmmglib1( PMMG_pParMesh parmesh )
{
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  int        it, i;

  /** Groups creation */
  if ( PMMG_SUCCESS != PMMG_split_grps( parmesh,REMESHER_TARGET_MESH_SIZE ) )
    return PMMG_STRONGFAILURE;

//DEBUGGING: grplst_meshes_to_saveMesh(parmesh->listgrp, 1, parmesh->myrank, "Begin_libparmmg1_proc");

  /** Reset the boundary fields between the old mesh size and the new one (Mmg
   * uses this fields assiming they are setted to 0)/ */
  for ( i=0; i<parmesh->ngrp; ++i ) {
    mesh         = parmesh->listgrp[i].mesh;
    memset(&mesh->xtetra[mesh->xt+1],0,(mesh->xtmax-mesh->xt)*sizeof(MMG5_xTetra));
    memset(&mesh->xpoint[mesh->xp+1],0,(mesh->xpmax-mesh->xp)*sizeof(MMG5_xPoint));

    /* if ( !_MMG3D_analys(mesh) ) return PMMG_STRONGFAILURE; */
  }

  /** Mesh adaptation */
  for ( it = 0; it < parmesh->niter; ++it ) {
    for ( i=0; i<parmesh->ngrp; ++i ) {
      mesh         = parmesh->listgrp[i].mesh;
      met          = parmesh->listgrp[i].met;

#ifdef PATTERN
      if ( 1 != _MMG5_mmg3d1_pattern( mesh, met ) )
        goto failed;
#else
      if ( 1 != _MMG5_mmg3d1_delone( mesh, met ) )
        goto failed;
#endif
      /** load Balancing at group scale and communicators reconstruction */
      if ( !PMMG_loadBalancing(parmesh) ) goto failed;
    }
  }

#warning add adjacendy update in merge_grp and merge_parmesh function and remove this
#warning NIKOS: I am not sure what is happening here. Finish error handling
  for ( i=0; i<parmesh->ngrp; ++i ) {
    mesh = parmesh->listgrp[i].mesh;
    met  = parmesh->listgrp[i].met;
    PMMG_DEL_MEM(mesh, mesh->adja, 4 * mesh->nemax + 5, int, "adjacency table" );
  }
  if ( PMMG_SUCCESS != PMMG_packParMesh(parmesh) )
    return PMMG_STRONGFAILURE;

//DEBUGGING:  saveGrpsToMeshes( parmesh->listgrp, parmesh->ngrp, parmesh->myrank, "2-packParMesh" );
  if ( !PMMG_merge_grps(parmesh) ) return PMMG_STRONGFAILURE;

  /* if ( !MMG3D_hashTetra(parmesh->listgrp[0].mesh,0) ) return PMMG_STRONGFAILURE; */
  /* char filename[200]; */
  /* sprintf(filename,"End_libparmmg1_proc%d.mesh",parmesh->myrank); */
  /* MMG3D_bdryBuild(parmesh->listgrp[0].mesh); */
  /* PMMG_saveMesh(parmesh,filename); */

  return PMMG_SUCCESS;

  /** mmg3d1_delone failure */
#warning NIKOS: These lines are exactly the same as in 334-344, handle ret_val correctly and remove code duplication
failed:
#warning add adjacendy update in merge grp and merge_parmesh function and remove this
  for ( i=0; i<parmesh->ngrp; ++i ) {
    mesh = parmesh->listgrp[i].mesh;
    met  = parmesh->listgrp[i].met;
    PMMG_DEL_MEM(mesh, mesh->adja, 4 * mesh->nemax + 5, int, "adjacency table" );
  }

  if ( PMMG_SUCCESS != PMMG_packParMesh(parmesh) )
    return PMMG_STRONGFAILURE;

  if ( !PMMG_merge_grps(parmesh) ) return PMMG_STRONGFAILURE;

  return PMMG_LOWFAILURE;
}

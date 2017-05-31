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
#include "grpsplit_pmmg.h" //PMMG_splitGrps

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
  MMG5_pTetra pt,ptnew;
  MMG5_pPrism pp;
  MMG5_pQuad  pq;
  MMG5_pPoint ppt,pptnew;
  int         *node2int_node_comm_index1;
  int         np,nc,k,ne,nbl,imet,imetnew,i,igrp;
  int         iadr;//,iadrnew,iadrv,*adjav,*adja,*adjanew,voy;

  for ( igrp=0; igrp<parmesh->ngrp; ++igrp ) {
    grp                       = &parmesh->listgrp[igrp];
    mesh                      = grp->mesh;
    met                       = grp->met;
    node2int_node_comm_index1 = grp->node2int_node_comm_index1;
    disp                      = grp->disp;
    mesh->memMax = PMMG_PMesh_SetMemMax(parmesh, mesh->memCur);

    /* compact vertices */
    np = nc = 0;
    for (k=1; k<=mesh->np; k++) {
      ppt = &mesh->point[k];
      if ( !MG_VOK(ppt) )
        continue;
      ppt->tmp = ++np;

      if ( mesh->info.nosurf && (ppt->tag & MG_NOSURF) )
        ppt->tag &= ~MG_REQ;

      if ( ppt->tag & MG_CRN )
        nc++;

      ppt->ref = abs(ppt->ref);
    }

    /* node index update in internal communicator */
    for (k=0; k<grp->nitem_int_node_comm; ++k) {
      iadr = node2int_node_comm_index1[k];

      ppt = &mesh->point[iadr];
      assert ( MG_VOK(ppt) );

      node2int_node_comm_index1[k] = ppt->tmp;
    }

    /* compact tetrahedra */
    ne  = 0;
    nbl = 1;
    for (k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( !MG_EOK(pt) )
        continue;

      pt->v[0] = mesh->point[pt->v[0]].tmp;
      pt->v[1] = mesh->point[pt->v[1]].tmp;
      pt->v[2] = mesh->point[pt->v[2]].tmp;
      pt->v[3] = mesh->point[pt->v[3]].tmp;
      ne++;
      if ( k!=nbl ) {
        ptnew = &mesh->tetra[nbl];
        memcpy(ptnew,pt,sizeof(MMG5_Tetra));
#warning update adjacency table: to keep for when the adjacency update will be implemented in mergeGrp and mergeParmesh functions
      /*   iadr = 4*(k-1) + 1; */
      /*   adja = &mesh->adja[iadr]; */
      /*   iadrnew = 4*(nbl-1) + 1; */
      /*   adjanew = &mesh->adja[iadrnew]; */
      /*   for(i=0 ; i<4 ; i++) { */
      /*     adjanew[i] = adja[i]; */
      /*     if(!adja[i]) continue; */
      /*     iadrv = 4*(adja[i]/4-1) +1; */
      /*     adjav = &mesh->adja[iadrv]; */
      /*     voy = i; */
      /*     adjav[adja[i]%4] = 4*nbl + voy; */
      /*   } */
      }
      nbl++;
    }
    mesh->ne = ne;

    /* update prisms and quads vertex indices */
    for (k=1; k<=mesh->nprism; k++) {
      pp = &mesh->prism[k];
      if ( !MG_EOK(pp) )
        continue;

      pp->v[0] = mesh->point[pp->v[0]].tmp;
      pp->v[1] = mesh->point[pp->v[1]].tmp;
      pp->v[2] = mesh->point[pp->v[2]].tmp;
      pp->v[3] = mesh->point[pp->v[3]].tmp;
      pp->v[4] = mesh->point[pp->v[4]].tmp;
      pp->v[5] = mesh->point[pp->v[5]].tmp;
    }
    for (k=1; k<=mesh->nquad; k++) {
      pq = &mesh->quadra[k];
      if ( !MG_EOK(pq) )
        continue;

      pq->v[0] = mesh->point[pq->v[0]].tmp;
      pq->v[1] = mesh->point[pq->v[1]].tmp;
      pq->v[2] = mesh->point[pq->v[2]].tmp;
      pq->v[3] = mesh->point[pq->v[3]].tmp;
    }

    /* compact metric */
    nbl = 1;
    if ( met && met->m ) {
      for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) )
          continue;
        imet    = k   * met->size;
        imetnew = nbl * met->size;

        for (i=0; i<met->size; i++)
          met->m[imetnew + i] = met->m[imet + i];
        ++nbl;
      }
    }

    /* compact displacement */
    nbl = 1;
    if ( disp && disp->m ) {
      for (k=1; k<=mesh->np; k++) {
        ppt = &mesh->point[k];
        if ( !MG_VOK(ppt) )
          continue;
        imet    = k   * disp->size;
        imetnew = nbl * disp->size;

        for (i=0; i<disp->size; i++)
          disp->m[imetnew + i] = disp->m[imet + i];
        ++nbl;
      }
    }

    /* compact vertices */
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
    if ( met && met->m )
      met->np  = np;
    if ( disp && disp->m )
      disp->np = np;

    /* create prism adjacency */
    if ( !MMG3D_hashPrism(mesh) ) {
      fprintf(stderr,"  ## Prism hashing problem. Exit program.\n");
      return PMMG_FAILURE;
    }

    /* Remove the MG_REQ tags added by the nosurf option */
    if ( mesh->info.nosurf ) {
      for (k=1; k<=mesh->ne; k++) {
        pt = &mesh->tetra[k];
        if ( MG_EOK(pt) &&  pt->xt ) {

          for (i=0; i<6; i++) {
            if ( mesh->xtetra[pt->xt].tag[i] & MG_NOSURF ) {
              mesh->xtetra[pt->xt].tag[i] &= ~MG_REQ;
              mesh->xtetra[pt->xt].tag[i] &= ~MG_NOSURF;
            }
          }
        }
      }
    }

    if ( mesh->info.imprim ) {
      fprintf(stdout,"     NUMBER OF VERTICES   %8d   CORNERS %8d\n",mesh->np,nc);
      fprintf(stdout,"     NUMBER OF ELEMENTS   %8d\n",mesh->ne);
    }

    for(k=1 ; k<=mesh->np ; k++)
      mesh->point[k].tmp = 0;

    mesh->npnil = mesh->np + 1;
    for(k=mesh->npnil; k<mesh->npmax-1; k++)
      mesh->point[k].tmp  = k+1;

    mesh->nenil = mesh->ne + 1;
    for(k=mesh->nenil; k<mesh->nemax-1; k++)
      mesh->tetra[k].v[3] = k+1;

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
  if ( PMMG_SUCCESS != PMMG_splitGrps( parmesh ) )
    return PMMG_STRONGFAILURE;

  /* sprintf(filename,"Begin_libparmmg1_proc%d.mesh",parmesh->myrank); */
  /* _MMG3D_bdryBuild(parmesh->listgrp[0].mesh); */
  /* PMMG_saveMesh(parmesh,filename); */
  /* if ( !_MMG3D_analys(mesh) ) return PMMG_STRONGFAILURE; */

  /** Mesh adaptation */
  for ( it = 0; it < parmesh->niter; ++it ) {
    for ( i=0; i<parmesh->ngrp; ++i ) {
      mesh         = parmesh->listgrp[i].mesh;
      met          = parmesh->listgrp[i].met;
      mesh->memMax = PMMG_PMesh_SetMemMax( parmesh, mesh->memCur );

#ifdef PATTERN
      if ( 1 != _MMG5_mmg3d1_pattern( mesh, met ) )
        goto failed;
#else
      if ( 1 != _MMG5_mmg3d1_delone( mesh, met ) )
        goto failed;
#endif

#warning Do we need to update the communicators? Does Mmg renum the boundary nodes with -nosurf option?
      /** load Balancing at group scale and communicators reconstruction */

    }
  }

#warning add adjacendy update in mergeGrp and mergeParMesh function and remove this
#warning NIKOS: I am not sure what is happening here. Finish error handling
  for ( i=0; i<parmesh->ngrp; ++i ) {
    mesh = parmesh->listgrp[i].mesh;
    met  = parmesh->listgrp[i].met;
    PMMG_FREE(mesh, mesh->adja, 4 * mesh->nemax + 5, int, "adjacency table" );
  }
  if ( PMMG_SUCCESS != PMMG_packParMesh(parmesh) )
    return PMMG_STRONGFAILURE;

  if ( PMMG_SUCCESS != PMMG_mergeGrps(parmesh) )
    return PMMG_STRONGFAILURE;

  /* if ( !MMG3D_hashTetra(parmesh->listgrp[0].mesh,0) ) return PMMG_STRONGFAILURE; */
  /* char filename[200]; */
  /* sprintf(filename,"End_libparmmg1_proc%d.mesh",parmesh->myrank); */
  /* _MMG3D_bdryBuild(parmesh->listgrp[0].mesh); */
  /* PMMG_saveMesh(parmesh,filename); */

  return PMMG_SUCCESS;

  /** mmg3d1_delone failure */
#warning NIKOS: These lines are exactly the same as in 334-344, handle ret_val correctly and remove code duplication
failed:
#warning add adjacendy update in mergeGrp and mergeParMesh function and remove this
  for ( i=0; i<parmesh->ngrp; ++i ) {
    mesh = parmesh->listgrp[i].mesh;
    met  = parmesh->listgrp[i].met;
    PMMG_FREE(mesh, mesh->adja, 4 * mesh->nemax + 5, int, "adjacency table" );
  }

  if ( PMMG_SUCCESS != PMMG_packParMesh(parmesh) )
    return PMMG_STRONGFAILURE;

  if ( PMMG_SUCCESS != PMMG_mergeGrps(parmesh) )
    return PMMG_STRONGFAILURE;

  return PMMG_LOWFAILURE;
}

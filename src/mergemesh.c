/**
 * \file mergemesh.c
 * \brief Merge the mesh.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 */


#include "libparmmg.h"


/**
 * \param parmesh pointer toward the mesh structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Merge the groups (mesh, internal communica) on the processor.
 *
 * \remark not yet tested
 */
int PMMG_mergeGrps(PMMG_pParMesh parmesh) {
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh0,mesh;
  MMG5_pSol      sol0,sol;
  MMG5_pTetra    pt0,pt;
  MMG5_pxTetra   pxt0,pxt;
  MMG5_pPoint    ppt0,ppt;
  MMG5_pxPoint   pxp0,pxp;
  PMMG_pint_comm int_node_comm;
  int            nitem_int_node_comm,nitem_int_node_comm0,*intvalues;
  int            *node2int_node_comm_index1,*node2int_node_comm0_index1;
  int            *node2int_node_comm_index2,*node2int_node_comm0_index2;
  int            idx,idx1,idx2,np,imsh,k,i,ie,ip;

  if ( parmesh->ngrp == 1 )  return(1);

  grp  = parmesh->listgrp;

  /** First step: Add all the interfaces points in grp 0 and update its internal
   * communicator */
  int_node_comm              = parmesh->int_node_comm;
  intvalues                  = int_node_comm->intvalues;

  mesh0                      = grp[0].mesh;
  nitem_int_node_comm0       = grp[0].nitem_int_node_comm;
  node2int_node_comm0_index1 = grp[0].node2int_node_comm_index1;
  node2int_node_comm0_index2 = grp[0].node2int_node_comm_index2;

  np = mesh0->np;

  /* Reset the tmp field of points (grp 0) */
  for ( k=1; k<=np; k++ ) {
    mesh0->point[k].tmp = 0;

    /* Travel through internal communicator to recover a unique numbering of
     * points (grp 0) */
    for ( k=0; k<nitem_int_node_comm0; ++k ) {
      intvalues[node2int_node_comm0_index2[k]] = node2int_node_comm0_index1[k];
    }
  }

  idx = nitem_int_node_comm;
  for ( imsh=1; imsh<parmesh->ngrp; ++imsh ) {
    mesh                      =  grp[imsh].mesh;
    nitem_int_node_comm       =  grp[imsh].nitem_int_node_comm;
    node2int_node_comm_index1 =  grp[imsh].node2int_node_comm_index1;
    node2int_node_comm_index2 =  grp[imsh].node2int_node_comm_index2;

    /* Reset the tmp field of points */
    for ( k=1; k<=mesh->np; k++ )
      mesh->point[k].tmp = 0;

    /* Travel through internal communicators to recover a unique numbering of
     * points: add the new interface points in grp 0, update the internal
     * communicator and save the new point idx in the tmp field of point */
    for ( k=0; k<nitem_int_node_comm; ++k ) {
      idx2 = node2int_node_comm_index2[k];

      if ( !intvalues[idx2] ) {
        idx1 = node2int_node_comm_index1[k];
        ppt = &mesh->point[idx1];

        /* Add the point ppt to the group 0 */
        ip = _MMG3D_newPt(mesh0,ppt->c,ppt->tag);
        if ( !ip ) {
          /* reallocation of point table */
          _MMG5_POINT_REALLOC(mesh0,sol0,ip,mesh0->gap,
                              printf("  ## Error: unable to merge group points\n");
                              _MMG5_INCREASE_MEM_MESSAGE();
                              return(0);
                              ,ppt->c,ppt->tag);
        }
        ppt->tmp = ip;
        /* Add xpoint if needed */
        if ( ppt->xp ) {
          pxp  = &mesh->xpoint[ppt->xp];
          ppt0 = &mesh0->point[ip];
          pxp0 = &mesh0->xpoint[ppt0->xp];
          memcpy(pxp0,pxp,sizeof(MMG5_xPoint));
        }
        if ( sol0->m ) {
          assert(sol->m);
          memcpy(&sol0->m[sol0->size*ip],&sol->m[sol->size*k],sol->size*sizeof(double));
        }
        node2int_node_comm0_index1[++idx] = ip;
        node2int_node_comm0_index2[  idx] = idx2;
        intvalues[idx2]                   = ip;
      }
      else {
        /* The point exists: store its index in grp 0 */
        ppt      = &mesh0->point[intvalues[idx2]];
        ppt->tmp = intvalues[idx2];
      }
    }
  }


  /** Second step: merge the meshes inside grp[0]->mesh */
  mesh0 = grp[0].mesh;
  sol0  = grp[0].sol;
#warning Do we need a smart fit of the mesh size?
  np = mesh0->np;
  for ( imsh=0; imsh<parmesh->ngrp; ++imsh ) {
    mesh = grp[imsh].mesh;
    sol  = grp[imsh].sol;

    /* Add the new points to our mesh */
    for ( k=1; k<=mesh->np; k++ ) {
      ppt = &mesh->point[k];

      if ( !ppt->tmp ) {
        ip = _MMG3D_newPt(mesh0,ppt->c,ppt->tag);
        if ( !ip ) {
          /* reallocation of point table */
          _MMG5_POINT_REALLOC(mesh0,sol0,ip,mesh0->gap,
                              printf("  ## Error: unable to merge group points\n");
                              _MMG5_INCREASE_MEM_MESSAGE();
                              return(0);
                              ,ppt->c,ppt->tag);
        }
        ppt->tmp = ip;
        /* Add xpoint if needed */
        if ( ppt->xp ) {
#warning add only the "true" xtetras (not those linked to the interfaces)
          pxp  = &mesh->xpoint[ppt->xp];
          ppt0 = &mesh0->point[ip];
          pxp0 = &mesh0->xpoint[ppt0->xp];
          memcpy(pxp0,ppt0,sizeof(MMG5_xPoint));
        }
        if ( sol0->m ) {
          assert(sol->m);
          memcpy(&sol0->m[sol0->size*ip],&sol->m[sol->size*k],sol->size*sizeof(double));
        }
      }
    }

    /* Add the new tetra to our mesh */
    for ( k=1; k<=mesh->ne; k++ ) {
      pt  = &mesh->tetra[k];
      ie  = _MMG3D_newElt(mesh0);
      if ( !ie ) {
         /* reallocation of tetra table */
        _MMG5_TETRA_REALLOC(mesh0,ie,mesh0->gap,
                            fprintf(stderr,"  ## Error: unable to merge group elts.\n");
                            _MMG5_INCREASE_MEM_MESSAGE();
                            return(0););
      }
      pt0 = &mesh0->tetra[ie];

      for ( i=0; i<4; ++i )
        pt0->v[i] = mesh->point[pt->v[i]].tmp;

      /* Add xtetra if needed */
      if ( pt->xt ) {
#warning add only the "true" xtetras (not those linked to the interfaces)
        pxt = &mesh->xtetra[pt->xt];

        mesh0->xt++;
        if ( mesh0->xt > mesh0->xtmax ) {
          /* realloc of xtetra table */
          _MMG5_TAB_RECALLOC(mesh0,mesh0->xtetra,mesh0->xtmax,0.2,MMG5_xTetra,
                             "larger xtetra table",
                             mesh0->xt--;
                             return(0));
        }
        pt0->xt = mesh0->xt;
        pxt0 = &mesh0->xtetra[pt0->xt];
        memcpy(pxt0,pxt,sizeof(MMG5_xTetra));
      }
    }

    /* Free the mesh */
    MMG3D_Free_all(MMG5_ARG_start,
                   MMG5_ARG_ppMesh,mesh,MMG5_ARG_ppMet,sol,
                   MMG5_ARG_end);
  }
  _MMG5_SAFE_REALLOC(parmesh->listgrp,1,PMMG_Grp,"(mergeGrps) listgrp");
  parmesh->ngrp = 1;

  return(1);
}


/**
 * \param parmesh pointer toward the mesh structure.
 * \return 0 if fail, 1 otherwise.
 *
 * Merge the mesh through the processors: The processors send their parmesh to
 * proc 0 that merge all the parmeshes.
 *
 */
int PMMG_mergeMesh(PMMG_pParMesh parmesh) {
  PMMG_pParMesh  *allParmeshes;
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
#warning add solution/metric communication
  MMG5_pSol      sol;
  PMMG_pint_comm int_node_comm;
  MPI_Comm       comm;
  MPI_Datatype   mpi_type_parmesh,mpi_type_grp,mpi_type_mesh;
  MPI_Datatype   mpi_type_int_node_comm,mpi_type_ext_node_comm;
  MPI_Datatype   mpi_type_point,mpi_type_xpoint,mpi_type_tetra,mpi_type_xtetra;
  int            nprocs,rank,imsh,k,idx;
  char           filename[13];

#warning to trash
  printf( " je suis le proc %d\n", parmesh->myrank);

  rank   = parmesh->myrank;
  nprocs = parmesh->nprocs;
  comm   = parmesh->comm;

  /** Step 1: merge the groups over each procs and return 1 group per proc.
   * This group contains a packed mesh where the triangle and edges are not
   * reconstructed (the mesh contains tetra and xtetra). */
  if ( !PMMG_mergeGrps(parmesh) ) return(0);

  /** Step 2: allocate an array of PMMG_parmesh that will store the PMMG_parmesh of
   *  each proc. Allocate internal and external communicators and fill its. */
  if ( !rank ) {
    _MMG5_SAFE_CALLOC(allParmeshes,nprocs,PMMG_pParMesh);
  }

  int_node_comm = parmesh->int_node_comm;
  _MMG5_SAFE_CALLOC(int_node_comm->intvalues,int_node_comm->nitem,int);

  grp = &parmesh->listgrp[0];
  for ( k=0; k<grp->nitem_int_node_comm; ++k ) {
    idx = grp->node2int_node_comm_index2[k];
    int_node_comm->intvalues[idx] = grp->node2int_node_comm_index1[k];
  }

  /** Step 3: Procs send their parmeshes to Proc 0 and Proc 0 recieve the data */
  if ( !PMMG_create_MPI_ParMesh(parmesh,&mpi_type_parmesh,
                                &mpi_type_int_node_comm,&mpi_type_ext_node_comm,
                                &mpi_type_grp,&mpi_type_mesh,
                                &mpi_type_point,&mpi_type_xpoint,
                                &mpi_type_tetra,&mpi_type_xtetra) )  return(0);

#warning MPI_gather does not work (maybe because the mpi_type_parmesh size depend on the processor that define it???)
  MPI_Gather( parmesh,1,mpi_type_parmesh,allParmeshes,1,mpi_type_parmesh,0, comm);

  PMMG_free_MPI_ParMesh(&mpi_type_parmesh,
                        &mpi_type_int_node_comm,&mpi_type_ext_node_comm,
                        &mpi_type_grp,&mpi_type_mesh,
                        &mpi_type_point,&mpi_type_xpoint,
                        &mpi_type_tetra,&mpi_type_xtetra);

#warning to trash
  printf( " je suis le proc %d\n", parmesh->myrank);
  for ( k=0; k<nprocs; ++k ) {
    /** Tetra adjacency reconstruction */
    _MMG5_SAFE_FREE(allParmeshes[k]->listgrp[0].mesh->adja);
    mesh = allParmeshes[k]->listgrp[0].mesh;
    if ( !MMG3D_hashTetra(mesh,1) ) {
      fprintf(stderr,"  ## PMMG Hashing problem (1). Exit program.\n");
      return(0);
    }

    /* Pack the mesh */
    if ( !_MMG3D_packMesh(mesh,NULL,NULL ) ) {
      fprintf(stderr,"  ## PMMG Packing problem (1). Exit program.\n");
      return(0);
    }

    sprintf(filename,"PPproc%d.mesh",rank);
    MMG3D_saveMesh(mesh,filename);
    return(0);
  }

  /** Step 4: Procs travel through their external communicators to recover the
   * numbering of the points shared with a lower proc. The other points are
   * concatenated with the proc 0. */
#warning may be optimized if we know the number of non-renumbered points of each mesh

   // TO IMPLEMENT


  return(1);
}

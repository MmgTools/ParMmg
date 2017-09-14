/**
 * \file mergemesh.c
 * \brief Merge the mesh.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 */
#include "parmmg.h"
#include "mpitypes_pmmg.h"

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return
 *         PMMG_FAILURE
 *         PMMG_SUCCESS
 *
 * Merge all meshes (mesh elements + internal communicator) of a group into the
 * first mesh of the group
 *
 * \remark not yet tested
 *
 * \warning the groups meshes must be packed.
 */
int PMMG_mergeGrps( PMMG_pParMesh parmesh )
{
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh0,mesh;
  MMG5_pSol      met0,met;
  MMG5_pTetra    pt0,pt;
  MMG5_pxTetra   pxt0,pxt;
  MMG5_pPoint    ppt0,ppt;
  MMG5_pxPoint   pxp0,pxp;
  PMMG_pext_comm ext_node_comm;
  PMMG_pint_comm int_node_comm;
  int            nitem_int_node_comm,nitem_int_node_comm0,*intvalues;
  int            *node2int_node_comm_index1,*node2int_node_comm0_index1;
  int            *node2int_node_comm_index2,*node2int_node_comm0_index2;
  int            poi_id_int, poi_id_glo, idx, np,imsh,k,i,ie,ip;
  int            new_nitem_int_node_comm,new_nitem_ext_node_comm;
  int            ret_val = PMMG_SUCCESS;

  if ( parmesh->ngrp == 1 )
    return PMMG_SUCCESS;

  grp  = parmesh->listgrp;

  /* Use parmesh internal communicator buffer (intvalues) to remember the
   * location in the merged mesh of all the interface points */
  int_node_comm              = parmesh->int_node_comm;
  PMMG_CALLOC(parmesh,int_node_comm->intvalues,int_node_comm->nitem + 1,int,
              "node communicator",return PMMG_FAILURE);

  intvalues                  = int_node_comm->intvalues;

  mesh0                      = grp[0].mesh;
  met0                       = grp[0].met;
  nitem_int_node_comm0       = grp[0].nitem_int_node_comm;
  node2int_node_comm0_index1 = grp[0].node2int_node_comm_index1;
  node2int_node_comm0_index2 = grp[0].node2int_node_comm_index2;

//DEBUGGING: saveGrpsToMeshes( parmesh->listgrp, parmesh->ngrp, parmesh->myrank, "BeforeMergeGrp" );

  /** First step: Merge interface points from all meshes into mesh0->points */

  /* Use the tmp field of points in meshes to remember the id in the merged mesh
   * of points that have already been added to the merged mesh or 0 if they
   * haven't been merged yet */
  for ( imsh = 0; imsh < parmesh->ngrp; ++imsh )
    for ( k = 1; k <= grp[imsh].mesh->np; ++k )
      grp[imsh].mesh->point[k].tmp = 0;

  /* mesh0's points are already in mesh0. Add a reference to their location in
   * mesh0(ie n2i_n_c0_idx1) from the merging communicator buffer(ie intvalues) */
  for ( k = 0; k < grp[0].nitem_int_node_comm; ++k ) {
    poi_id_int = grp[0].node2int_node_comm_index1[k];
    poi_id_glo = grp[0].node2int_node_comm_index2[k];

    intvalues[ poi_id_glo ] = poi_id_int;
  }

  /* Add points referenced by the rest of the groups meshes' internal
   * communicators to the merging mesh (mesh0) and in intvalues */
  for ( imsh=1; imsh<parmesh->ngrp; ++imsh ) {
    mesh                      = grp[imsh].mesh;
    met                       = grp[imsh].met;
    nitem_int_node_comm       = grp[imsh].nitem_int_node_comm;
    node2int_node_comm_index1 = grp[imsh].node2int_node_comm_index1;
    node2int_node_comm_index2 = grp[imsh].node2int_node_comm_index2;

    for ( k=0; k<nitem_int_node_comm; ++k ) {
      poi_id_glo = node2int_node_comm_index2[ k ];
      assert ( ( 0 <= poi_id_glo ) && ( poi_id_glo < (int_node_comm->nitem+1) ) && "weird..." );

      // location in currently working mesh where point data actually are
      assert ( node2int_node_comm_index1[ k ] < mesh->npmax );
      ppt = &mesh->point[ node2int_node_comm_index1[ k ] ];

      /* point ppt is not found in the merged mesh. add it */
      if ( !intvalues[ poi_id_glo ] ) {
        ip = _MMG3D_newPt(mesh0,ppt->c,ppt->tag);
        if ( !ip ) {
          /* reallocation of point table */
#warning NIKOS: The correct law to pass to _MMG5_POINT_REALLOC is the commented code(the incompletely created array needs to be freed before exiting):
//                              PMMG_DEL_MEM(parmesh,
//                                           intvalues,int_node_comm->nitem+1,int,
//                                           "freeing temp buffer",
//                                           return PMMG_FAILURE),
          _MMG5_POINT_REALLOC(mesh0,met0,ip,mesh0->gap,
                              printf("  ## Error: unable to merge group points\n");
                              _MMG5_INCREASE_MEM_MESSAGE();
                              return PMMG_FAILURE;,
                              ppt->c,ppt->tag,PMMG_FAILURE);

        }
        assert( (ip < grp[ 0 ].mesh->npmax) && "run out of points" );
        //NIKOS TODO: Remember location in merged mesh in mesh, mesh0 and communicator tmp fields
        ppt->tmp = ip;
        intvalues[ poi_id_glo ] = ip;

        /* Add xpoint if needed */
        if ( ppt->xp ) {
          pxp  = &mesh->xpoint[ppt->xp];
          ppt0 = &mesh0->point[ip];
          pxp0 = &mesh0->xpoint[ppt0->xp];
          assert( (ppt0->xp < mesh0->xpmax) && "increase xpoints" );
          assert( (ppt0->xp > 0) && "negative xpoints?!?!?!?!?" );
          memcpy(pxp0,pxp,sizeof(MMG5_xPoint));
        }
        if ( met0->m ) {
          assert(met->m);
          assert( (met->size < 7) && "met size issues" );
          memcpy(&met0->m[met0->size*ip],&met->m[met->size*k],met->size*sizeof(double));
        }
      } else {
        /* point already exists in merged mesh. update his tmp field to point to
         * its mesh0 index */
        ppt->tmp = intvalues[ poi_id_glo ];
      }
    }
  }


  /* Second step: merge internal points and tetras of all meshes into mesh0 */
#warning Do we need a smart fit of the mesh size?
  np = mesh0->np;
  for ( imsh=1; imsh<parmesh->ngrp; ++imsh ) {
    mesh = grp[imsh].mesh;
    met  = grp[imsh].met;

    /* Loop over points and add the ones that are not already in the merdeg
     * mesh (mesh0) */
    for ( k=1; k<=mesh->np; k++ ) {
      ppt = &mesh->point[k];

      if ( !ppt->tmp ) {
        ip = _MMG3D_newPt(mesh0,ppt->c,ppt->tag);
        if ( !ip ) {
          /* reallocation of point table */
#warning NIKOS: The correct law to pass to _MMG5_POINT_REALLOC is the commented code:
//                              PMMG_DEL_MEM(parmesh,
//                                           intvalues,int_node_comm->nitem+1,int,
//                                           "freeing temp buffer",
//                                           return PMMG_FAILURE),
          _MMG5_POINT_REALLOC(mesh0,met0,ip,mesh0->gap,
                              printf("  ## Error: unable to merge group points\n");
                              _MMG5_INCREASE_MEM_MESSAGE();
                              return PMMG_FAILURE;,
                              ppt->c,ppt->tag,PMMG_FAILURE);
        }
        ppt->tmp = ip;
        /* Add xpoint if needed */
        if ( ppt->xp ) {
#warning add only the "true" xtetras (not those linked to the interfaces)
          pxp  = &mesh->xpoint[ppt->xp];
          ppt0 = &mesh0->point[ip];
          pxp0 = &mesh0->xpoint[ppt0->xp];
          memcpy(pxp0,pxp,sizeof(MMG5_xPoint));
        }
        if ( met0->m ) {
          assert(met->m);
          memcpy(&met0->m[met0->size*ip],&met->m[met->size*k],met->size*sizeof(double));
        }
      }
    }

    /* Add current meshs' tetras to the merged mesh (mesh0) */
    for ( k=1; k<=mesh->ne; k++ ) {
      pt  = &mesh->tetra[k];
      ie  = _MMG3D_newElt(mesh0);
      if ( !ie ) {
#warning NIKOS: The correct law to pass to _MMG5_TETRA_REALLOC is the commented code:
//                            PMMG_DEL_MEM(parmesh,
//                                         intvalues,int_node_comm->nitem,int,
//                                         "freeing temp buffer",
//                                         return PMMG_FAILURE),
        _MMG5_TETRA_REALLOC(mesh0,ie,mesh0->gap,
                            fprintf(stderr,"  ## Error: unable to merge group elts.\n");
                            _MMG5_INCREASE_MEM_MESSAGE();
                            return PMMG_SUCCESS;,
                            PMMG_FAILURE);
      }
      pt0 = &mesh0->tetra[ie];

      for ( i=0; i<4; ++i )
        pt0->v[i] = mesh->point[pt->v[i]].tmp;
      pt0->ref = pt->ref;

#warning add adjacency update

      /* Add xtetra if needed */
      if ( pt->xt ) {
#warning add only the "true" xtetras (not those linked to the interfaces)
        pxt = &mesh->xtetra[pt->xt];

        mesh0->xt++;
        if ( mesh0->xt > mesh0->xtmax ) {
          PMMG_RECALLOC(mesh0, mesh0->xtetra, 1.2 * mesh0->xtmax + 1,
                        mesh0->xtmax + 1, MMG5_xTetra,
                        "larger xtetra table", mesh0->xt--; goto fail_ncomm);
          mesh0->xtmax = 1.2 * mesh0->xtmax;
        }
        pt0->xt = mesh0->xt;
        pxt0 = &mesh0->xtetra[pt0->xt];
        memcpy(pxt0,pxt,sizeof(MMG5_xTetra));
      }
    }

    /* Free the mesh */
    MMG3D_Free_all(MMG5_ARG_start,
                   MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&met,
                   MMG5_ARG_end);
  }

  /** Update the communicators */
  /* Reset the tmp field of the point: it will be used to store the position of
   * a point in the internal communicator */
  for ( k=1; k<=mesh0->np; k++ )
    mesh0->point[k].tmp = 0;

  /* Travel through the external communicators and udpate all the communicators */
  poi_id_int = 0;
  poi_id_glo = 0;
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    // currently working external communicator
    ext_node_comm = &parmesh->ext_node_comm[k];

    // initialize new_nitem: used when increasing comm size
    new_nitem_ext_node_comm = ext_node_comm->nitem;

    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx = ext_node_comm->int_comm_index[i];
      ip  = intvalues[idx];
      ppt = &mesh0->point[ip];

      /* New point in the internal communicator */
      if ( !ppt->tmp ) {
        if ( poi_id_int == grp[0].nitem_int_node_comm ) {
          new_nitem_int_node_comm = (int)(1.2*grp[0].nitem_int_node_comm);
          PMMG_REALLOC(parmesh,grp[0].node2int_node_comm_index1,
                       new_nitem_int_node_comm,
                       grp[0].nitem_int_node_comm,int,
                       "(mergeGrps) node2int_node_comm_index1",
                       grp[0].nitem_int_node_comm = new_nitem_int_node_comm;
                       goto fail_ncomm);
          PMMG_REALLOC(parmesh,grp[0].node2int_node_comm_index2,
                       new_nitem_int_node_comm,
                       grp[0].nitem_int_node_comm,int,
                       "(mergeGrps) node2int_node_comm_index2",
                       grp[0].nitem_int_node_comm = new_nitem_int_node_comm;
                       goto fail_ncomm);
          grp[0].nitem_int_node_comm = new_nitem_int_node_comm;
          node2int_node_comm0_index1 = grp[0].node2int_node_comm_index1;
          node2int_node_comm0_index2 = grp[0].node2int_node_comm_index2;
        }
        node2int_node_comm0_index1[ poi_id_int ]     = ip;
        node2int_node_comm0_index2[ poi_id_int ]     = poi_id_int;

        if ( poi_id_glo >= new_nitem_ext_node_comm ) {
          new_nitem_ext_node_comm *= 1.2;
          PMMG_REALLOC(parmesh,ext_node_comm->int_comm_index,
                       new_nitem_ext_node_comm,
                       ext_node_comm->nitem,int,"(mergeGrps) ext_node_comm",
                       goto fail_ncomm);
        }
        ext_node_comm->int_comm_index[ poi_id_glo++ ]= poi_id_int;
        ppt->tmp                                     = poi_id_int + k*mesh0->np;
        poi_id_int++;
      /* The point already has a position in the internal comm */
      } else {
        if ( ppt->tmp < k * mesh0->np ) {
          /* The point has been stored in the internal comm by another external
           * comm: update its position in our external comm */
          if ( poi_id_glo >= new_nitem_ext_node_comm ) {

            new_nitem_ext_node_comm *= 1.2;
            PMMG_REALLOC(parmesh,ext_node_comm->int_comm_index,
                         new_nitem_ext_node_comm,
                         ext_node_comm->nitem,int,"(mergeGrps) ext_node_comm",
                         goto fail_ncomm);
          }
          ext_node_comm->int_comm_index[ poi_id_glo++] = ppt->tmp;
        }
      }
    }
    assert( (poi_id_glo != 0)  && "empty communicator?????" );
    PMMG_REALLOC(parmesh,ext_node_comm->int_comm_index, poi_id_glo,
                 new_nitem_ext_node_comm,int,"(mergeGrps) ext_node_comm",
                 ext_node_comm->nitem = poi_id_glo;goto fail_ncomm);
    ext_node_comm->nitem = poi_id_glo;
  }

  PMMG_REALLOC(parmesh,grp[0].node2int_node_comm_index1,poi_id_int,
               grp[0].nitem_int_node_comm,int,"(mergeGrps) node2int_node_comm_index1",
               goto fail_ncomm);
  PMMG_REALLOC(parmesh,grp[0].node2int_node_comm_index2,poi_id_int,
               grp[0].nitem_int_node_comm,int,"(mergeGrps) node2int_node_comm_index2",
               goto fail_ncomm);
  int_node_comm->nitem       = poi_id_int;
  grp[0].nitem_int_node_comm = poi_id_int;

  PMMG_DEL_MEM(parmesh,parmesh->int_node_comm->intvalues,parmesh->int_node_comm->nitem,
            int,"release int_n_comm intvalues");

  parmesh->int_node_comm->nitem = poi_id_int;

#warning NIKOS: REPLACE THIS WITH PMMG. THIS CAUSES A VALGRIND ERROR THAT I HAVNET YET UNDERSTOOD. I GUESS IT WILL BE FIXED BY REWRITTING THE ALLOCATIONS IN MERGEMESH/SPLITGRP
//  _MMG5_SAFE_REALLOC(grp,1,PMMG_Grp,"(mergeGrps) listgrp",0);

  parmesh->ngrp = 1;
  return PMMG_SUCCESS;

fail_ncomm:
  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,int_node_comm->nitem+1,int,"node communicator");
  return PMMG_FAILURE;
}


/**
 * \param parmesh pointer toward the parmesh structure.
 * \param merge 1 if the groups needs to be merged on the processor.
 * \return 0 if fail, 1 otherwise.
 *
 * Merge the mesh through the processors: The processors send their parmesh to
 * proc 0 that merge all the parmeshes.
 *
 * \warning the meshes must be packed before calling this procedure.
 */
int PMMG_mergeParMesh(PMMG_pParMesh parmesh, int merge) {
  PMMG_pGrp      grp;
  MMG5_pMesh     mesh;
  MMG5_pPoint    rcv_point,point_1, point_2,ppt;
  MMG5_pxPoint   rcv_xpoint,xpoint,pxp;
  MMG5_pTetra    rcv_tetra,tetra,pt;
  MMG5_pxTetra   rcv_xtetra,xtetra,pxt;
  MMG5_pSol      met;
  PMMG_pint_comm int_node_comm;
  PMMG_pext_comm ext_node_comm;
  MPI_Comm       comm;
  MPI_Datatype   mpi_point,mpi_xpoint,mpi_tetra,mpi_xtetra;
  double         *rcv_met,*met_1;
  int            *rcv_np,*rcv_ne,*rcv_xp,*rcv_xt,*rcv_nmet;
  int            *point_displs,*xpoint_displs,*tetra_displs,*xtetra_displs;
  int            *met_displs,nmet_tot,ne_tot;
  int            *intval_displs,*rcv_intvalues,nitem_int_node_comm_tot;
  int            *rcv_node2int_node_comm_index1,*rcv_node2int_node_comm_index2;
  int            *node2int_node_comm_index1,*node2int_node_comm_index2;
  int            *node2int_node_comm_index1_2,*node2int_node_comm_index2_2;
  int            *rcv_nitem_int_node_comm,*rcv_next_node_comm;
  int            *ext_comm_displs,*color_in_tab,*color_out_tab,*nitem_ext_tab;
  int            *rcv_color_in_tab,*rcv_color_out_tab,*rcv_nitem_ext_tab;
  int            *int_comm_index_displs,*int_comm_index,*nitems_ext_idx;
  int            *rcv_int_comm_index;
  int            *intvalues_1, *intvalues_2,*int_comm_index_2;
  int            ext_comm_displs_tot,nitem_ext_tot,nitems_1,nitems_2;
  int            nprocs,rank,k,i,j,idx,idx_2,cursor,color_in,color_out;
  int            np,ne,xt_tot,np_tot,xp_tot,nnpar,l,isMet;
  char           filename[32];

  rank   = parmesh->myrank;
  nprocs = parmesh->nprocs;
  comm   = parmesh->comm;


//DEBUGGING: saveGrpsToMeshes( PMMG_pGrp listgrp, 0, parmesh->myrank, "mesh" );

  /** Step 1: merge the groups over each procs and return 1 group per proc.
   * This group contains a packed mesh where the triangle and edges are not
   * reconstructed (the mesh contains tetra and xtetra). */
  if ( merge && !PMMG_mergeGrps(parmesh) ) return(0);

  /** Step 2: store the boundary entities into the xpoints and xtetra entities */
  grp  = &parmesh->listgrp[0];
  mesh  = grp->mesh;
  met   = grp->met;
  isMet = (met->m)? 1 : 0;

  /** Step 3: Allocate internal and external communicators and fill its: the
   *  intvalues array contains the indices of the matching nodes on the proc. */
  int_node_comm = parmesh->int_node_comm;

  if ( !int_node_comm->nitem ) return 1;

  _MMG5_SAFE_CALLOC(int_node_comm->intvalues,int_node_comm->nitem,int,0);

  assert(int_node_comm->nitem == grp->nitem_int_node_comm);

  for ( k=0; k<grp->nitem_int_node_comm; ++k ) {
    idx = grp->node2int_node_comm_index2[k];
    int_node_comm->intvalues[idx] = grp->node2int_node_comm_index1[k];
  }

  /** Step 4: Procs send their parmeshes to Proc 0 and Proc 0 recieve the data */
  _MMG5_SAFE_CALLOC(rcv_np        ,nprocs,int,0);
  _MMG5_SAFE_CALLOC(rcv_ne        ,nprocs,int,0);
  _MMG5_SAFE_CALLOC(rcv_xp        ,nprocs,int,0);
  _MMG5_SAFE_CALLOC(rcv_xt        ,nprocs,int,0);
  if(isMet)
    _MMG5_SAFE_CALLOC(rcv_nmet    ,nprocs,int,0);
  _MMG5_SAFE_CALLOC(point_displs  ,nprocs,int,0);
  _MMG5_SAFE_CALLOC(xpoint_displs ,nprocs,int,0);
  _MMG5_SAFE_CALLOC(tetra_displs  ,nprocs,int,0);
  _MMG5_SAFE_CALLOC(xtetra_displs ,nprocs,int,0);
  if(isMet)
    _MMG5_SAFE_CALLOC(met_displs    ,nprocs,int,0);

  _MMG5_SAFE_CALLOC(rcv_nitem_int_node_comm,nprocs,int,0);
  _MMG5_SAFE_CALLOC(rcv_next_node_comm,nprocs,int,0);

  /* if ( parmesh->ddebug ) { */
  /*   sprintf(filename,"Before_Gather_proc%d.mesh",rank); */
  /*   _MMG3D_bdryBuild(parmesh->listgrp[0].mesh); */
  /*   MMG3D_saveMesh(mesh,filename); */
  /*   MMG3D_saveSol(mesh,met,filename); */
  /* } */

  PMMG_create_MPI_Point (mesh->point,  &mpi_point );
  PMMG_create_MPI_xPoint(mesh->xpoint, &mpi_xpoint);
  PMMG_create_MPI_Tetra (mesh->tetra,  &mpi_tetra );
  PMMG_create_MPI_xTetra(mesh->xtetra, &mpi_xtetra);

  /* Gather parmesh size infos on proc 0 */
#warning try to compare with non-blocking comms (Igather(v))
  MPI_Gather(&mesh->np,1,MPI_INT,rcv_np,1,MPI_INT,0,comm);
  MPI_Gather(&mesh->ne,1,MPI_INT,rcv_ne,1,MPI_INT,0,comm);
  MPI_Gather(&mesh->xp,1,MPI_INT,rcv_xp,1,MPI_INT,0,comm);
  MPI_Gather(&mesh->xt,1,MPI_INT,rcv_xt,1,MPI_INT,0,comm);
  MPI_Gather(&grp->nitem_int_node_comm,1,MPI_INT,
              rcv_nitem_int_node_comm,1,MPI_INT,0,comm);
  MPI_Gather(&parmesh->next_node_comm,1,MPI_INT,
              rcv_next_node_comm,1,MPI_INT,0,comm);

  /* Gather meshes on proc 0 */
  /* Points */
  for ( k=1; k<=mesh->np; ++k ) mesh->point[k].tmp = 0;

  point_displs[0] = 0;
  for ( k=1; k<nprocs; ++k ) {
    point_displs[k] += point_displs[k-1] + rcv_np[k-1];
  }
  np_tot = point_displs[nprocs-1]+rcv_np[nprocs-1];
  rcv_point = (MMG5_pPoint)calloc(np_tot+1,sizeof(MMG5_Point));

  MPI_Gatherv(&mesh->point[1],mesh->np,mpi_point,&rcv_point[1],
               rcv_np,point_displs,mpi_point,0,comm);

  /* Solutions */
  if(isMet) {
    for ( k=0; k<nprocs; ++k ) {
      rcv_nmet[k] = met->size*rcv_np[k];
    }

    met_displs[0] = 0;
    for ( k=1; k<nprocs; ++k ) {
      met_displs[k] += met_displs[k-1] + rcv_nmet[k-1];
    }
    nmet_tot = met_displs[nprocs-1]+rcv_nmet[nprocs-1];
    rcv_met = (double*)calloc(nmet_tot+met->size,sizeof(double));

    MPI_Gatherv(&met->m[met->size],mesh->np*met->size,MPI_DOUBLE,&rcv_met[met->size],
                rcv_nmet,met_displs,MPI_DOUBLE,0,comm);
  }
  /* Tetra */
  tetra_displs[0] = 0;
  for ( k=1; k<nprocs; ++k ) {
    tetra_displs[k] += tetra_displs[k-1] + rcv_ne[k-1];
  }
  ne_tot    = tetra_displs[nprocs-1]+rcv_ne[nprocs-1];
  rcv_tetra = (MMG5_pTetra)calloc(ne_tot+1,sizeof(MMG5_Tetra));
  MPI_Gatherv(&mesh->tetra[1],mesh->ne,mpi_tetra,&rcv_tetra[1],
              rcv_ne,tetra_displs,mpi_tetra,0,comm);

  /* xPoints */
  xpoint_displs[0] = 0;
  for ( k=1; k<nprocs; ++k ) {
    xpoint_displs[k] += xpoint_displs[k-1] + rcv_xp[k-1];
  }
  xp_tot = xpoint_displs[nprocs-1]+rcv_xp[nprocs-1];
  rcv_xpoint = (MMG5_pxPoint)calloc(xp_tot+1,sizeof(MMG5_xPoint));

  MPI_Gatherv(&mesh->xpoint[1],mesh->xp,mpi_xpoint,&rcv_xpoint[1],
              rcv_xp,xpoint_displs,mpi_xpoint,0,comm);

  /* xTetra */
  xtetra_displs[0] = 0;
  for ( k=1; k<nprocs; ++k ) {
    xtetra_displs[k] += xtetra_displs[k-1] + rcv_xt[k-1];
  }
  xt_tot = xtetra_displs[nprocs-1]+rcv_xt[nprocs-1];
  rcv_xtetra = (MMG5_pxTetra)calloc(xt_tot+1,sizeof(MMG5_xTetra));

  MPI_Gatherv(&mesh->xtetra[1],mesh->xt,mpi_xtetra,&rcv_xtetra[1],
              rcv_xt,xtetra_displs,mpi_xtetra,0,comm);


  /* Internal communicator */
  _MMG5_SAFE_CALLOC(intval_displs,nprocs,int,0);
  intval_displs[0] = 0;
  for ( k=1; k<nprocs; ++k ) {
    intval_displs[k] += intval_displs[k-1] + rcv_nitem_int_node_comm[k-1];
  }
  nitem_int_node_comm_tot = intval_displs[nprocs-1]
    + rcv_nitem_int_node_comm[nprocs-1];

  rcv_intvalues                 = malloc(nitem_int_node_comm_tot*sizeof(int));
  rcv_node2int_node_comm_index1 = malloc(nitem_int_node_comm_tot*sizeof(int));
  rcv_node2int_node_comm_index2 = malloc(nitem_int_node_comm_tot*sizeof(int));

  MPI_Gatherv(int_node_comm->intvalues,int_node_comm->nitem,MPI_INT,
               rcv_intvalues,rcv_nitem_int_node_comm,intval_displs,MPI_INT,
              0,comm);
  MPI_Gatherv(grp->node2int_node_comm_index1,int_node_comm->nitem,MPI_INT,
              rcv_node2int_node_comm_index1,rcv_nitem_int_node_comm,
              intval_displs,MPI_INT,0,comm);
  MPI_Gatherv(grp->node2int_node_comm_index2,int_node_comm->nitem,MPI_INT,
              rcv_node2int_node_comm_index2,rcv_nitem_int_node_comm,
              intval_displs,MPI_INT,0,comm);


  /* External communicator */
  /* For each proc, we concatenate the external communicators fields into arrays */
  color_in_tab  =(int*)malloc(parmesh->next_node_comm*sizeof(int));
  color_out_tab =(int*)malloc(parmesh->next_node_comm*sizeof(int));
  nitem_ext_tab =(int*)malloc(parmesh->next_node_comm*sizeof(int));

  nitem_ext_tot = 0;
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm    = &parmesh->ext_node_comm[k];
    color_in_tab[k ] = ext_node_comm->color_in;
    color_out_tab[k] = ext_node_comm->color_out;
    nitem_ext_tab[k] = ext_node_comm->nitem;
    nitem_ext_tot += ext_node_comm->nitem;
  }

  int_comm_index = (int*)malloc(nitem_ext_tot*sizeof(int));
  idx = 0;
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    memcpy(&int_comm_index[idx],ext_node_comm->int_comm_index,
           ext_node_comm->nitem*sizeof(int));
    idx += ext_node_comm->nitem;
  }

  _MMG5_SAFE_CALLOC(ext_comm_displs,nprocs,int,0);
  ext_comm_displs[0] = 0;
  for ( k=1; k<nprocs; ++k ) {
    ext_comm_displs[k] += ext_comm_displs[k-1] + rcv_next_node_comm[k-1];
  }
  ext_comm_displs_tot = ext_comm_displs[nprocs-1] + rcv_next_node_comm[nprocs-1];

  rcv_color_in_tab  = (int*)malloc(nitem_ext_tot*sizeof(int));
  rcv_color_out_tab = (int*)malloc(nitem_ext_tot*sizeof(int));
  rcv_nitem_ext_tab = (int*)malloc(nitem_ext_tot*sizeof(int));

  MPI_Gatherv(color_in_tab,parmesh->next_node_comm,MPI_INT,
              rcv_color_in_tab,rcv_next_node_comm,ext_comm_displs,MPI_INT,
              0,comm);
  MPI_Gatherv(color_out_tab,parmesh->next_node_comm,MPI_INT,
              rcv_color_out_tab,rcv_next_node_comm,ext_comm_displs,MPI_INT,
              0,comm);
  MPI_Gatherv(nitem_ext_tab,parmesh->next_node_comm,MPI_INT,
              rcv_nitem_ext_tab,rcv_next_node_comm,ext_comm_displs,MPI_INT,
              0,comm);

  _MMG5_SAFE_CALLOC(int_comm_index_displs,nprocs,int,0);
  _MMG5_SAFE_CALLOC(nitems_ext_idx,nprocs,int,0);

  int_comm_index_displs[0] = 0;
  nitems_ext_idx[0] = 0;
  idx = 0;
  for ( i=0; i<rcv_next_node_comm[0]; ++i ) {
    nitems_ext_idx[0] += rcv_nitem_ext_tab[idx++];
  }
  for ( k=1; k<nprocs; ++k ) {
    nitems_ext_idx[k] = 0;
    for ( i=0; i<rcv_next_node_comm[k]; ++i ) {
      nitems_ext_idx[k] += rcv_nitem_ext_tab[idx++];
    }
    int_comm_index_displs[k] = int_comm_index_displs[k-1] + nitems_ext_idx[k-1];
  }

  rcv_int_comm_index = (int*)malloc((int_comm_index_displs[nprocs-1]+
                                     nitems_ext_idx[nprocs-1])*sizeof(int));

  MPI_Gatherv(int_comm_index,nitem_ext_tot,MPI_INT,
              rcv_int_comm_index,nitems_ext_idx,int_comm_index_displs,MPI_INT,
              0,comm);

  /* Free structures */
  _MMG5_SAFE_FREE(color_in_tab);
  _MMG5_SAFE_FREE(color_out_tab);
  _MMG5_SAFE_FREE(nitem_ext_tab);
  _MMG5_SAFE_FREE(int_comm_index);

  _MMG5_DEL_MEM(mesh,mesh->point,(mesh->npmax+1)*sizeof(MMG5_Point));
  _MMG5_DEL_MEM(mesh,mesh->tetra,(mesh->nemax+1)*sizeof(MMG5_Tetra));
  _MMG5_DEL_MEM(mesh,mesh->xpoint,(mesh->xpmax+1)*sizeof(MMG5_xPoint));
  _MMG5_DEL_MEM(mesh,mesh->xtetra,(mesh->xtmax+1)*sizeof(MMG5_Tetra));
  if(isMet)
    _MMG5_DEL_MEM(mesh,met->m,(met->npmax+1)*met->size*sizeof(double));

  _MMG5_SAFE_FREE(grp->node2int_edge_comm_index1);
  _MMG5_SAFE_FREE(grp->node2int_edge_comm_index2);
  _MMG5_SAFE_FREE(parmesh->int_node_comm->intvalues);

  for ( i=0; i<parmesh->next_node_comm; ++i ) {
    _MMG5_SAFE_FREE(parmesh->ext_node_comm->int_comm_index);
  }
  _MMG5_SAFE_FREE(parmesh->ext_node_comm);

  /** Step 5: Proc 0 merges the meshes: We travel through the external
   * communicators to recover the numbering of the points shared with a lower
   * proc. The other points are concatenated with the proc 0. */
  if ( !rank ) {
    np = 0;

    /* Mesh renumbering to have the same indices at points shared by multiple
     * processors; The new indices are stored in the tmp field of the MMG5_Point
     * structure. */
    for ( k=0; k<nprocs; ++k ) {
      /* Points and internal communicators for proc \a k */
      point_1     = &rcv_point[point_displs[k]];
      cursor      = intval_displs[k];
      intvalues_1 = &rcv_intvalues[cursor];
      node2int_node_comm_index1 = &rcv_node2int_node_comm_index1[cursor];
      node2int_node_comm_index2 = &rcv_node2int_node_comm_index2[cursor];

      /* Travel through the external communicators that lists the points at the
       * interface of the procs color_in and color_out: if color_in<color_out,
       * create new indices for the points over the color_in proc. If
       * color_out<color_in, get the indices from the matching points over the
       * proc color_out. */
      idx = int_comm_index_displs[k];

      for ( i=0; i<rcv_next_node_comm[k]; ++i ) {
        cursor   = ext_comm_displs[k]+i;
        nitems_1 = rcv_nitem_ext_tab[cursor];
        int_comm_index = &rcv_int_comm_index[idx];
        idx     += nitems_1;

        /* External communicator k->color_out */
        color_in    = rcv_color_in_tab[cursor];
        color_out   = rcv_color_out_tab[cursor];

        assert( color_in==k );
        assert( color_in!=color_out );

        if ( color_in < color_out ) {
          /* New point */
          for ( j=0; j<nitems_1; ++j ) {
            if ( !point_1[intvalues_1[int_comm_index[j]]].tmp ) {
              point_1[intvalues_1[int_comm_index[j]]].tmp = ++np;
            }
          }
        }
        else {
          /* Get the point index from the neighbouring proc */

          /* Find the matching external communicator over proc color_out (->k) */
          idx_2 =  int_comm_index_displs[color_out];
          for ( j=0; j<rcv_next_node_comm[color_out]; ++j ) {
            cursor   = ext_comm_displs[color_out]+j;
            nitems_2 = rcv_nitem_ext_tab[cursor];
            int_comm_index_2 = &rcv_int_comm_index[idx_2];
            idx_2   += nitems_2;

            if ( rcv_color_in_tab[cursor]  == color_out &&
                 rcv_color_out_tab[cursor] == color_in  )  break;
          }
          assert(j<rcv_next_node_comm[color_out]);
          assert(nitems_1 == nitems_2);


          /* Points and internal communicators for proc \a color_out */
          point_2     = &rcv_point[point_displs[color_out]];

          cursor      = intval_displs[color_out];
          intvalues_2 = &rcv_intvalues[cursor];
          node2int_node_comm_index1_2 = &rcv_node2int_node_comm_index1[cursor];
          node2int_node_comm_index2_2 = &rcv_node2int_node_comm_index2[cursor];

          /* Update point indices (stored in the tmp field) */
          for ( j=0; j<nitems_1; ++j ) {
            point_1[intvalues_1[int_comm_index[j]]].tmp =
              point_2[intvalues_2[int_comm_index_2[j]]].tmp;
          }
        }
      }

      /* Create new indices for the points that haven't been seen. */
      for ( i=1; i<=rcv_np[k]; ++i ) {
        if ( !point_1[i].tmp ) {
          point_1[i].tmp = ++np;
        }
      }
    }

    /* Tetra + xTetra */
    mesh->nemax = mesh->ne = ne_tot;
    mesh->nenil = 0;
    _MMG5_SAFE_CALLOC(mesh->xtetra,xt_tot+1,MMG5_xTetra,0);
    _MMG5_SAFE_CALLOC(mesh->tetra,mesh->nemax+1,MMG5_Tetra,0);
    mesh->xt = xt_tot;
    ne = idx = 0;
    for ( k=0; k<nprocs; ++k ) {
      xtetra      = &rcv_xtetra[xtetra_displs[k]];
      tetra       = &rcv_tetra[tetra_displs[k]];
      point_1     = &rcv_point[point_displs[k]];

      for ( i=1; i<=rcv_ne[k]; ++i ) {
        pt = &mesh->tetra[tetra_displs[k]+i];

        for ( j=0; j<4; ++j ) {
          tetra[i].v[j] = point_1[tetra[i].v[j]].tmp;
        }

        memcpy(pt,&tetra[i],sizeof(MMG5_Tetra));

        if ( tetra[i].xt ) {
          nnpar = 0;
          pxt = &xtetra[tetra[i].xt];
          for ( l=0; l<4; ++l ) {
            if ( pxt->ftag[l] && !(pxt->ftag[l] & MG_PARBDY) ) ++nnpar;
#warning doing this we loose the truely required entities
            if ( (pxt->ftag[l] & MG_PARBDY) && (pxt->ftag[l] & MG_REQ) )
              pxt->ftag[l] &= ~MG_REQ;
          }

          if ( !nnpar ) {
            pt->xt = 0;
          }
          else {
            memcpy(&mesh->xtetra[++ne],pxt,sizeof(MMG5_xTetra));
            pt->xt = ne;
          }
        }
      }
    }
    mesh->xt=ne;

    /* Points and solutions */
    mesh->np = met->np = np;
    mesh->npmax = met->npmax = mesh->np;
    mesh->npnil = 0;
    _MMG5_SAFE_CALLOC(mesh->point,mesh->npmax+1,MMG5_Point,0);
    if(isMet)
      _MMG5_SAFE_CALLOC(met->m,(met->npmax+1)*met->size,double,0);

    for ( i=1; i<=mesh->np; ++i ) mesh->point[i].tag = MG_NUL;

    np = 0;
    for ( k=0; k<nprocs; ++k ) {
      point_1     = &rcv_point[point_displs[k]];
      if(isMet)
        met_1       = &rcv_met[met_displs[k]];

      for ( i=1; i<=rcv_np[k]; ++i ) {
        idx = point_1[i].tmp;
        assert(idx);
        ppt = &mesh->point[idx];

        if ( MG_VOK(ppt) ) {
          point_1[i].tmp = 0;
          continue;
        }

        memcpy(ppt,&point_1[i],sizeof(MMG5_Point));
        ppt->tmp = 0;

        if(isMet)
          memcpy( &met->m[idx*met->size],
                  &met_1[i*met->size],
                  met->size*sizeof(double) );

        if ( point_1[i].xp ) ++np;
      }
    }

    /* xPoints */
    _MMG5_SAFE_CALLOC(mesh->xpoint,np+1,MMG5_xPoint,0);
    np = 0;
    for ( k=0; k<nprocs; ++k ) {
      point_1     = &rcv_point[point_displs[k]];
      xpoint      = &rcv_xpoint[xpoint_displs[k]];

      for ( i=1; i<=rcv_np[k]; ++i ) {
        idx = point_1[i].tmp;

        if ( !idx ) continue;
        ppt = &mesh->point[idx];

        if ( !point_1[i].xp ) continue;

        pxp = &mesh->xpoint[++np];
        memcpy(pxp, &xpoint[point_1[i].xp],sizeof(MMG5_xPoint));
        ppt->xp = np;
      }
    }
  }

  /* Free memory */
  MPI_Type_free(&mpi_point);
  MPI_Type_free(&mpi_xpoint);
  MPI_Type_free(&mpi_tetra);
  MPI_Type_free(&mpi_xtetra);

  _MMG5_SAFE_FREE(rcv_np);
  _MMG5_SAFE_FREE(rcv_ne);
  _MMG5_SAFE_FREE(rcv_xp);
  _MMG5_SAFE_FREE(rcv_xt);
  if(isMet)
    _MMG5_SAFE_FREE(rcv_nmet);
  _MMG5_SAFE_FREE(rcv_nitem_int_node_comm);
  _MMG5_SAFE_FREE(rcv_next_node_comm);
  _MMG5_SAFE_FREE(point_displs);
  _MMG5_SAFE_FREE(tetra_displs);
  _MMG5_SAFE_FREE(xpoint_displs);
  _MMG5_SAFE_FREE(xtetra_displs);
  if(isMet)
    _MMG5_SAFE_FREE(met_displs);
  _MMG5_SAFE_FREE(intval_displs);
  _MMG5_SAFE_FREE(rcv_point);
  _MMG5_SAFE_FREE(rcv_tetra);
  _MMG5_SAFE_FREE(rcv_xpoint);
  _MMG5_SAFE_FREE(rcv_xtetra);
  if(isMet)
    _MMG5_SAFE_FREE(rcv_met);
  _MMG5_SAFE_FREE(rcv_intvalues);
  _MMG5_SAFE_FREE(rcv_node2int_node_comm_index1);
  _MMG5_SAFE_FREE(rcv_node2int_node_comm_index2);
  _MMG5_SAFE_FREE(ext_comm_displs);
  _MMG5_SAFE_FREE(rcv_color_in_tab);
  _MMG5_SAFE_FREE(rcv_color_out_tab);
  _MMG5_SAFE_FREE(rcv_nitem_ext_tab);
  _MMG5_SAFE_FREE(int_comm_index_displs);
  _MMG5_SAFE_FREE(nitems_ext_idx);
  _MMG5_SAFE_FREE(rcv_int_comm_index);

  return(1);
}

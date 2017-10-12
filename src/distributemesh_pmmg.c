/**
 * \file distributemesh.c
 * \brief Distribute the mesh on the processors.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */
#include <mpi.h>
#include "parmmg.h"
#include "mpitypes_pmmg.h"
#include "metis_pmmg.h"

/**
 * \param xtetra pointer toward a table containing the xtetra structures.
 * \param *perm pointer toward the permutation table (to perform in place
 * permutations).
 * \param ind1 index of the first xtetra to swap.
 * \param ind2 index of the second xtetra to swap.
 *
 * Swap two xtetra in the table of xtetrahedras.
 */
static void swapxTetra( MMG5_pxTetra xtetra, int* perm, int ind1, int ind2 )
{
  MMG5_xTetra pxttmp;
  int         tmp;

  /** 1- swap the xtetra */
  memcpy(&pxttmp      ,&xtetra[ind2],sizeof(MMG5_xTetra));
  memcpy(&xtetra[ind2],&xtetra[ind1],sizeof(MMG5_xTetra));
  memcpy(&xtetra[ind1],&pxttmp      ,sizeof(MMG5_xTetra));

  /** 2- swap the permutation table */
  tmp        = perm[ind2];
  perm[ind2] = perm[ind1];
  perm[ind1] = tmp;
}

/**
 * \param xpoint pointer toward a table containing the xpoint structures.
 * \param *perm pointer toward the permutation table (to perform in place
 * permutations).
 * \param ind1 index of the first xpoint to swap.
 * \param ind2 index of the second xpoint to swap.
 *
 * Swap two xpoint in the table of xpoints.
 */
static void swapxPoint( MMG5_pxPoint xpoint, int* perm, int ind1, int ind2 )
{
  MMG5_xPoint pxptmp;
  int         tmp;

  /** 1- swap the xpoint */
  memcpy(&pxptmp      ,&xpoint[ind2],sizeof(MMG5_xPoint));
  memcpy(&xpoint[ind2],&xpoint[ind1],sizeof(MMG5_xPoint));
  memcpy(&xpoint[ind1],&pxptmp      ,sizeof(MMG5_xPoint));

  /** 2- swap the permutation table */
  tmp        = perm[ind2];
  perm[ind2] = perm[ind1];
  perm[ind1] = tmp;
}

/**
 * \param point pointer toward a table containing the point structures.
 * \param met pointer toward a table containing the metric structure.
 * \param *perm pointer toward the permutation table (to perform in place
 * permutations).
 * \param ind1 index of the first xpoint to swap.
 * \param ind2 index of the second xpoint to swap.
 * \param metsize size of the metric (1=iso,6=aniso).
 *
 * Swap two points in the table of points.
 */
static void swapPoint( MMG5_pPoint point, double* met,int* perm,
                       int ind1, int ind2, int metsiz )
{
  MMG5_Point ppttmp;
  MMG5_Sol   mettmp;
  int        tmp,addr2,addr1;

  /** 1- swap the xpoint */
  memcpy(&ppttmp      ,&point[ind2], sizeof(MMG5_Point));
  memcpy(&point[ind2] ,&point[ind1], sizeof(MMG5_Point));
  memcpy(&point[ind1] ,&ppttmp      ,sizeof(MMG5_Point));

  /** 2- swap the mets */
  if ( met ) {
    addr1 = ind1*metsiz;
    addr2 = ind2*metsiz;
    memcpy(&mettmp    ,&met[addr2],metsiz*sizeof(double));
    memcpy(&met[addr2],&met[addr1],metsiz*sizeof(double));
    memcpy(&met[addr1],&mettmp    ,metsiz*sizeof(double));
  }

  /** 3- swap the permutation table */
  tmp        = perm[ind2];
  perm[ind2] = perm[ind1];
  perm[ind1] = tmp;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \return PMMG_FAILURE
 *         PMMG_SUCCESS
 *
 * Send the initial mesh from proc 0 toward the other procs.
 */
#warning NIKOS TODO: could add a parameter in PMMG_bcastMesh to select which node broadcasts the mesh
#warning NIKOS TODO: perhaps change to: int PMMG_bcastMesh(MMG5_pMesh msh, int myid, int comm, int source);
int PMMG_bcastMesh( PMMG_pParMesh parmesh )
{
  PMMG_pGrp    grp;
  MMG5_pMesh   mesh;
  MMG5_pSol    met;
  MPI_Datatype mpi_light_point, mpi_light_tetra, mpi_tria,mpi_edge;
  int          rank;

  /** Proc 0 send the mesh to the other procs */
  grp    = &parmesh->listgrp[0];
  mesh   = grp->mesh;
  met    = grp->met;
  rank   = parmesh->myrank;

  /* Mesh */
  MPI_Bcast( &mesh->np,     1, MPI_INT,       0, parmesh->comm );
  MPI_Bcast( &mesh->ne,     1, MPI_INT,       0, parmesh->comm );
  MPI_Bcast( &mesh->nt,     1, MPI_INT,       0, parmesh->comm );
  MPI_Bcast( &mesh->na,     1, MPI_INT,       0, parmesh->comm );
  MPI_Bcast( &mesh->ntmax,  1, MPI_INT,       0, parmesh->comm );
  MPI_Bcast( &mesh->memMax, 1, MPI_LONG_LONG, 0, parmesh->comm );

  mesh->nemax = mesh->nei = mesh->ne;
  mesh->nenil = 0;
  mesh->npmax = mesh->npi = mesh->np;
  mesh->npnil = 0;
  mesh->nti   = mesh->nt;
  mesh->xtmax = mesh->ntmax;

  /* Metric */
  MPI_Bcast( &met->size,  1, MPI_INT, 0, parmesh->comm );
  MPI_Bcast( &met->type,  1, MPI_INT, 0, parmesh->comm );
  MPI_Bcast( &met->npmax, 1, MPI_INT, 0, parmesh->comm );
  MPI_Bcast( &met->np,    1, MPI_INT, 0, parmesh->comm );

  met->npi = met->np;
  met->ver = mesh->ver;
  met->dim = mesh->dim;

  if ( rank ) {
#warning NIKOS: DO WE NEED TO CLEANUP mesh member allocations or are they handled in mesh deallocation?
    PMMG_CALLOC(mesh,mesh->point,mesh->npmax+1,MMG5_Point,"initial vertices", return PMMG_FAILURE);

    PMMG_CALLOC(mesh,mesh->tetra,mesh->nemax+1,MMG5_Tetra,"initial tetrahedra",return PMMG_FAILURE);

    if ( mesh->nt )
      PMMG_CALLOC(mesh,mesh->tria,mesh->nt+1,MMG5_Tria,"initial triangles",return PMMG_SUCCESS);

    if ( mesh->na )
      PMMG_CALLOC(mesh,mesh->edge,mesh->na+1,MMG5_Edge,"initial edges",return PMMG_FAILURE);

    if ( met->npmax )
      PMMG_CALLOC(mesh,met->m,met->size*(met->npmax+1),double,"initial edges",return PMMG_FAILURE);
  }

  PMMG_create_MPI_lightPoint( mesh->point, &mpi_light_point );
  PMMG_create_MPI_lightTetra( mesh->tetra, &mpi_light_tetra );
  if ( mesh->nt )
    PMMG_create_MPI_Tria( mesh->tria, &mpi_tria );
  if ( mesh->na )
    PMMG_create_MPI_Edge( mesh->edge, &mpi_edge );

  MPI_Bcast( mesh->point, mesh->np + 1, mpi_light_point, 0, parmesh->comm );
  MPI_Bcast( mesh->tetra, mesh->ne + 1, mpi_light_tetra, 0, parmesh->comm );
  if ( mesh->nt )
    MPI_Bcast( mesh->tria, mesh->nt + 1, mpi_tria, 0, parmesh->comm );
  if ( mesh->na )
    MPI_Bcast( mesh->edge, mesh->na + 1, mpi_edge, 0, parmesh->comm );
  if ( met->m )
    MPI_Bcast( met->m, met->size * (met->npmax + 1), MPI_DOUBLE, 0, parmesh->comm );

  MPI_Type_free( &mpi_light_point );
  MPI_Type_free( &mpi_light_tetra );
  if ( mesh->nt )
    MPI_Type_free( &mpi_tria );
  if ( mesh->na )
    MPI_Type_free( &mpi_edge );

  return PMMG_SUCCESS;
}

/**
 * \param parmesh pointer toward the mesh structure.
 * \param part pointer toward an array of int containing the partitions.
 *
 * \return PMMG_FAILURE
 *         PMMG_SUCCESS
 *
 * Delete parts of the mesh not on the processor.
 */
int PMMG_distributeMesh( PMMG_pParMesh parmesh )
{
  PMMG_pGrp      grp = NULL;
  MMG5_pMesh     mesh = NULL;
  MMG5_pSol      met = NULL;
  MMG5_pTetra    pt = NULL, ptnew = NULL;
  MMG5_pxTetra   pxt = NULL;
  MMG5_pPoint    ppt = NULL;
  PMMG_pext_comm pext_comm = NULL;
  idx_t          *part = NULL;
  int            nprocs = 0 ,rank = 0, np = 0, ne = 0, nxt = 0, nxp = 0;
  int            ip = 0, iploc = 0, ifac = 0, i = 0, j = 0, k = 0, *idx = NULL;
  int            kvois = 0, rankVois =0;
  int            *node2int_node_comm_index1 = NULL;
  int            *node2int_node_comm_index2 = NULL;
  int            nitem_int_node_comm = 0, next_node_comm = 0;
  int            inIntComm = 0, nbl = 0;
  int            *shared_pt = NULL;
  int            *pointPerm = NULL, *xTetraPerm = NULL, *xPointPerm = NULL;
  int8_t         *pointRanks = NULL;
  int            ret_val = PMMG_SUCCESS;
  int            old_val = 0;
  MPI_Datatype   metis_dt;

  /** Proc 0 send the mesh to the other procs */
  nprocs = parmesh->nprocs;
  grp    = &parmesh->listgrp[0];
  mesh   = grp->mesh;
  met    = grp->met;
  rank   = parmesh->myrank;

  /** Call metis for partionning*/
  PMMG_CALLOC(parmesh,part,mesh->ne,idx_t,"allocate metis buffer",
              ret_val = PMMG_FAILURE;goto fail_alloc0);

#warning Perhaps I could change this 0 to be user configurable (via PMMG_distributeMesh function argument)
  if ( (!parmesh->myrank) && nprocs > 1 ) {
    if (    PMMG_partition_metis( parmesh, part, parmesh->nprocs )
         != PMMG_SUCCESS ) {
      ret_val = PMMG_FAILURE;
      goto fail_alloc1;
    }
  }
  if ( IDXTYPEWIDTH == 32 )
    metis_dt = MPI_INT32_T;
  else if ( IDXTYPEWIDTH == 64 )
    metis_dt = MPI_INT64_T;
  else {
    printf("  ## Error: %s: unable to detect the metis integer width (%d).\n",
           __func__,IDXTYPEWIDTH);
    goto fail_alloc1;
  }

  MPI_Bcast( &part[0], mesh->ne, metis_dt, 0, parmesh->comm );

  /** Remove the part of the mesh that are not on the proc rank */
  PMMG_CALLOC(parmesh,shared_pt,nprocs,int,"dist Mesh buffer0 ",
              ret_val = PMMG_FAILURE;goto fail_alloc1);
  PMMG_CALLOC(parmesh,pointRanks,nprocs*mesh->np,int8_t,"dist mesh buffer1",
              ret_val = PMMG_FAILURE; goto fail_alloc2);
  PMMG_CALLOC(parmesh,pointPerm,mesh->np+1,int,"dist mesh buffer2",
              ret_val = PMMG_FAILURE;goto fail_alloc3);
  PMMG_CALLOC(parmesh,xTetraPerm,mesh->xtmax+1,int,"dist mesh buffer3",
              ret_val = PMMG_FAILURE; goto fail_alloc4);
  PMMG_CALLOC(parmesh,xPointPerm,mesh->xp+1,int,"dist mesh buffer4",
              ret_val = PMMG_FAILURE;goto fail_alloc5);

  nxp = 0;
  nxt = 0;
  np  = 0;
  ne  = 0;
  nitem_int_node_comm = 0;

  /* Reset the tmp field of points */
  for ( k=1; k<=mesh->np; k++ )
    mesh->point[k].tmp = 0;

  /** Mark mesh entities that will stay on the proc and count the number of
   * point that must be communicated to the other procs  */
  for ( k=1; k<=mesh->ne; k++ ) {
    pt = &mesh->tetra[k];
    pt->mark = part[k-1];

    if ( pt->mark != rank )
      continue;

    for ( ifac=0; ifac<4; ifac++ ) {
      kvois = mesh->adja[4*k-3+ifac]/4;

      if ( kvois )
        rankVois = part[kvois-1];
      else
        rankVois = rank;

      /* Mark the interfaces between two procs */
      if ( rank != rankVois ) {
        if ( !pt->xt ) {
          if ( (mesh->xt + 1) > mesh->xtmax ) {
            /* realloc of xtetras table */
            PMMG_RECALLOC(mesh,mesh->xtetra,1.2*mesh->xtmax+1,mesh->xtmax+1,int,
                          "larger xtetra ",
                          ret_val = PMMG_FAILURE;goto fail_alloc6);
            PMMG_RECALLOC(parmesh,xTetraPerm,1.2*mesh->xtmax+1,mesh->xtmax+1,
                          int,"larger tetra permutation table ",
                          ret_val = PMMG_FAILURE; goto fail_alloc6);
            mesh->xtmax = 1.2 * mesh->xtmax;
          }
          ++mesh->xt;
          pt->xt = mesh->xt;
        }
        pxt = &mesh->xtetra[pt->xt];
        /* Parallel face */
        pxt->ftag[ifac] |= (MG_PARBDY + MG_BDY + MG_REQ + MG_NOSURF);

        /* Parallel edges */
        for ( j=0; j<3; ++j )
          pxt->tag[_MMG5_iarf[ifac][j]] |= (MG_PARBDY + MG_BDY + MG_REQ + MG_NOSURF);
      }

      for ( j=0; j<3; ++j ) {
        iploc = _MMG5_idir[ifac][j];
        ip    = pt->v[iploc];
        ppt   = &mesh->point[ip];

        /* Count (and mark) each point that will be shared between me and the
         * proc rankVois */
        if ( rankVois != rank && !pointRanks[nprocs*(ip-1)+rankVois] ) {
          pointRanks[nprocs*(ip-1)+rankVois] = 1;
          ++shared_pt[rankVois];

          /* Mark parallel vertex */
          ppt->tag |= (MG_PARBDY + MG_BDY + MG_REQ + MG_NOSURF);
        }

        if ( !ppt->tmp ) {
          /* Mark the new point index and update the table of permutation */
          ppt->tmp = ++np;
          pointPerm[ip] = np;

          /* update the table of permutation for xPoint if needed */
          if ( ppt->xp ) {
            xPointPerm[ppt->xp] = ++nxp;
            ppt->xp             = nxp;
          }
        }
      }
    }

    /* Update the table of permutation of the Tetra and the tetra vertices indices */
    for ( j=0; j<4; ++j )
      pt->v[j] = mesh->point[pt->v[j]].tmp;

    /* update the table of permutation for xTetra if needed */
    if ( !pt->xt )
      continue;

    xTetraPerm[pt->xt] = ++nxt;
    pt->xt = nxt;
  }
  PMMG_DEL_MEM(mesh,mesh->adja,4*mesh->nemax+5,int,"dealloc mesh adja");

  /** Count the number of external node communicators and initialize it */
  next_node_comm = 0;
  for ( k=0; k<nprocs; ++k )
    if ( shared_pt[k] )
      ++next_node_comm;

  old_val = parmesh->next_node_comm;
  parmesh->next_node_comm = next_node_comm;
  PMMG_CALLOC(parmesh,parmesh->ext_node_comm,next_node_comm,PMMG_ext_comm,
              "allocate node comm ",
              parmesh->next_node_comm = old_val; ret_val = PMMG_FAILURE; goto fail_alloc6);

  next_node_comm = 0;
  for ( k=0; k<nprocs; ++k ) {
    if ( shared_pt[k] ) {
      pext_comm = &parmesh->ext_node_comm[next_node_comm];
      pext_comm->color_in  = rank;
      pext_comm->color_out = k;
      old_val              = pext_comm->nitem;
      pext_comm->nitem     = shared_pt[k];
      PMMG_CALLOC(parmesh,pext_comm->int_comm_index,pext_comm->nitem,int,
                  "allocate comm idx",
                  pext_comm->nitem = old_val; ret_val = PMMG_FAILURE; goto fail_alloc6);
      /* Use shared_pt to store the idx of the external communicator me->k */
      shared_pt[k] = next_node_comm++;
    }
  }

  /** Initialize the internal node communicator */
  nitem_int_node_comm = 0;
  for ( k=1; k<=mesh->np; k++ ) {

    if ( !mesh->point[k].tmp )
      continue;

    for ( j=0; j<nprocs; ++j ) {
      if ( pointRanks[nprocs*(k-1)+j] ) {
        ++nitem_int_node_comm;
        break;
      }
    }
  }

  old_val = grp->nitem_int_node_comm;
  grp->nitem_int_node_comm = nitem_int_node_comm;
  PMMG_CALLOC(parmesh,grp->node2int_node_comm_index1,nitem_int_node_comm,int,
              "alloc n2i_n_c_idx1 ",
              grp->nitem_int_node_comm = old_val; ret_val = PMMG_FAILURE; goto fail_alloc6);
  // For the error handling to be complete at this point we have to:
  //   1) reverse the reallocation of n2i_n_c_idx1
  //   2) reset nitem_int_node_comm, set ret_val, etc
  PMMG_CALLOC(parmesh,grp->node2int_node_comm_index2, nitem_int_node_comm,int,
              "alloc n2i_n_c_idx2 ",
              PMMG_REALLOC(parmesh,grp->node2int_node_comm_index1,old_val,nitem_int_node_comm,
                           int, "realloc n2i_n_c_idx1 ", );
              grp->nitem_int_node_comm = old_val;
              ret_val = PMMG_FAILURE;
              goto fail_alloc6);

  /** Travel through the mesh and fill the communicators */
  i = 0;
  node2int_node_comm_index1 = grp->node2int_node_comm_index1;
  node2int_node_comm_index2 = grp->node2int_node_comm_index2;

  /* Idx is used to store the external communicator cursor */
  PMMG_CALLOC(parmesh,idx,parmesh->next_node_comm,int,"allocating idx",
              ret_val = PMMG_FAILURE; goto fail_alloc6);

  for ( k=1; k<=mesh->np; k++ ) {
    if ( !mesh->point[k].tmp )
      continue;

    inIntComm = 0;
    for ( j=0; j<nprocs; ++j ) {
      pext_comm = &parmesh->ext_node_comm[shared_pt[j]];

      if ( pointRanks[nprocs*(k-1)+j] == 1 ) {
        /* Add point in external communicator */
        pext_comm->int_comm_index[idx[shared_pt[j]]++] = i;

        if ( !inIntComm ) {
          /* Add point in internal communicator */
          inIntComm = 1;
          node2int_node_comm_index1[i] = mesh->point[k].tmp;
          node2int_node_comm_index2[i] = i;
        }
      }
    }
    /* Increment internal comm cursor */
    if ( inIntComm )
      ++i;
  }
#warning NIKOS: I need some help with managing error handling here: am I consistently managing the communicator deallocations? are they in a consistent state if an error happens and this returns?
  PMMG_CALLOC(parmesh,parmesh->int_node_comm,1,PMMG_int_comm,"allocating idx",
              ret_val = PMMG_FAILURE; goto fail_alloc7);
  /* We have 1 Grp per proc, thus : int_node_comm.nitem : nitem_int_node_comm */
  parmesh->int_node_comm->nitem = nitem_int_node_comm;

  /** Compact tetrahedra on the proc */
  ne  = 0;
  nbl = 1;
  for ( k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];

    if ( pt->mark != rank )
      continue;
    ++ne;

    if ( k!=nbl ) {
      ptnew = &mesh->tetra[nbl];
      memcpy(ptnew,pt,sizeof(MMG5_Tetra));
    }
    ++nbl;
  }
  mesh->ne = ne;

  /** Compact xtetra on the proc */
  for ( k=1; k<=mesh->xt; ++k )
    while ( xTetraPerm[k] != k && xTetraPerm[k] )
      swapxTetra(mesh->xtetra,xTetraPerm,k,xTetraPerm[k]);
  mesh->xt = nxt;

  /** Compact vertices on the proc: in place permutations */
  for ( k=1; k<=mesh->np; ++k )
    while ( pointPerm[k] != k && pointPerm[k] )
      swapPoint(mesh->point,met->m,pointPerm,k,pointPerm[k],met->size);
  mesh->np = np;
  met->np = np;

  /** Compact xpoint on the proc: in place permutations */
  for ( k=1; k<=mesh->xp; ++k )
    while ( xPointPerm[k] != k && xPointPerm[k] )
      swapxPoint(mesh->xpoint,xPointPerm,k,xPointPerm[k]);
  mesh->xp = nxp;

  /** Update xtetra edge tags */
  if ( PMMG_SUCCESS != PMMG_bdryUpdate( mesh ) ) {
    ret_val = PMMG_FAILURE;
    goto fail_alloc7;
  }


  /** Adjacency reconstruction */
  if ( 1 != MMG3D_hashTetra( parmesh->listgrp[0].mesh, 0 ) ) {
    ret_val = PMMG_FAILURE;
    goto fail_alloc7;
  }

//  if ( parmesh->ddebug ) {
//    grplst_meshes_to_saveMesh( parmesh->listgrp, 1, parmesh->myrank, "End_distributeMesh_proc");
//    if ( met )
//      PMMG_saveSol( parmesh, filename );
//  }

fail_alloc7:
  PMMG_DEL_MEM(parmesh,idx,parmesh->next_node_comm,int,"deallocating idx");
fail_alloc6:
  PMMG_DEL_MEM(parmesh,xPointPerm,mesh->xp+1,int,"deallocate metis buffer5");
fail_alloc5:
  PMMG_DEL_MEM(parmesh,xTetraPerm,mesh->xtmax+1,int,"deallocate metis buffer4");
fail_alloc4:
  PMMG_DEL_MEM(parmesh,pointPerm,mesh->np+1,int,"deallocate metis buffer3");
fail_alloc3:
  PMMG_DEL_MEM(parmesh,pointRanks,nprocs*mesh->np,int8_t,"deallocate metis buffer2");
fail_alloc2:
  PMMG_DEL_MEM(parmesh,shared_pt,nprocs,int,"deallocate metis buffer1");
fail_alloc1:
  PMMG_DEL_MEM(parmesh,part,mesh->ne,idx_t,"deallocate metis buffer0");
fail_alloc0:
  return ret_val;
}

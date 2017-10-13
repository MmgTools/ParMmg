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
 * \param mesh pointer toward a MMG5 mesh structure
 * \param rank MPI rank
 *
 * \return 1 if success, 0 if fail
 *
 * Pack the tetrahedra and remove those ones that are not on the processor (the
 * proc index is stored in pt->mark).
 *
 */
static inline
PMMG_packTetra(MMG5_pMesh mesh, int rank) {
  MMG5_pTetra pt,ptnew;
  int         ne,nbl,k;

  ne  = 0;
  nbl = 1;
  for ( k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];

    if ( (!MG_EOK(pt)) || (pt->mark != rank) )
      continue;

    pt->v[0] = mesh->point[pt->v[0]].tmp;
    pt->v[1] = mesh->point[pt->v[1]].tmp;
    pt->v[2] = mesh->point[pt->v[2]].tmp;
    pt->v[3] = mesh->point[pt->v[3]].tmp;

    ++ne;

    if ( k!=nbl ) {
      ptnew = &mesh->tetra[nbl];
      memcpy(ptnew,pt,sizeof(MMG5_Tetra));
    }
    ++nbl;
  }
  mesh->ne = ne;

  return 1;
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
  MMG5_pTetra    pt = NULL;
  MMG5_pxTetra   pxt = NULL;
  MMG5_pPoint    ppt = NULL;
  PMMG_pext_comm pext_node_comm,pext_face_comm;
  idx_t          *part = NULL;
  int            nprocs = 0 ,rank = 0, np = 0, nxt = 0, nxp = 0;
  int            ip = 0, iploc = 0, ifac = 0, i = 0, j = 0, k = 0, *idx = NULL;
  int            kvois = 0, rankCur,rankVois,ifacVois;
  int            *node2int_node_comm_index1,*node2int_node_comm_index2;
  int            *node2int_face_comm_index1,*node2int_face_comm_index2;
  int            nitem_int_node_comm,nitem_int_face_comm;
  int            next_node_comm,next_face_comm;
  int            inIntComm = 0;
  int            *shared_pt,*shared_face;
  int            *pointPerm = NULL, *xTetraPerm = NULL, *xPointPerm = NULL;
  int8_t         *seen_shared_pt = NULL;
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
  shared_pt = shared_face = NULL;
  PMMG_CALLOC(parmesh,shared_pt,nprocs,int,"shared_pt array",
              ret_val = PMMG_FAILURE;goto fail_alloc1);
  PMMG_CALLOC(parmesh,seen_shared_pt,nprocs*mesh->np,int8_t,"seen_shared_pt array",
              ret_val = PMMG_FAILURE; goto fail_alloc2);
  PMMG_CALLOC(parmesh,shared_face,nprocs,int,"shared_face array",
              ret_val = PMMG_FAILURE;goto fail_alloc3);
  PMMG_CALLOC(parmesh,pointPerm,mesh->np+1,int,"dist mesh buffer2",
              ret_val = PMMG_FAILURE;goto fail_alloc4);
  PMMG_CALLOC(parmesh,xTetraPerm,mesh->xtmax+1,int,"dist mesh buffer3",
              ret_val = PMMG_FAILURE; goto fail_alloc5);
  PMMG_CALLOC(parmesh,xPointPerm,mesh->xp+1,int,"dist mesh buffer4",
              ret_val = PMMG_FAILURE;goto fail_alloc6);

  nxp = 0;
  nxt = 0;
  np  = 0;
  nitem_int_node_comm = 0;

  /* Reset the tmp field of points (it will be used to store the local index of
   * the point on each proc) */
  for ( k=1; k<=mesh->np; k++ )
    mesh->point[k].tmp = 0;

  /* Reset the base field of tetras (it will be used to store the local index of
   * the tetra on each proc) */
  mesh->base = 0;
  for ( k=1; k<=mesh->ne; k++ )
    mesh->tetra[k].base = 0;

  /** Mark mesh entities that will stay on the proc and count the number of
   * point that must be communicated to the other procs  */
  for ( k=1; k<=mesh->ne; k++ ) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    pt->mark = part[k-1];

    if ( pt->mark != rank )
      continue;

    pt->base = ++mesh->base;

    for ( ifac=0; ifac<4; ifac++ ) {
      kvois = mesh->adja[4*k-3+ifac]/4;

      if ( kvois )
        rankVois = part[kvois-1];
      else
        rankVois = rank;

      /* Mark the interfaces between two procs */
      if ( rank != rankVois ) {
        ++shared_face[rankVois];

        if ( !pt->xt ) {
          if ( (mesh->xt + 1) > mesh->xtmax ) {
            /* realloc of xtetras table */
            PMMG_RECALLOC(mesh,mesh->xtetra,1.2*mesh->xtmax+1,mesh->xtmax+1,int,
                          "larger xtetra ",
                          ret_val = PMMG_FAILURE;goto fail_alloc7);
            PMMG_RECALLOC(parmesh,xTetraPerm,1.2*mesh->xtmax+1,mesh->xtmax+1,
                          int,"larger tetra permutation table ",
                          ret_val = PMMG_FAILURE; goto fail_alloc7);
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
        if ( rankVois != rank && !seen_shared_pt[nprocs*(ip-1)+rankVois] ) {
          seen_shared_pt[nprocs*(ip-1)+rankVois] = 1;
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

    /* update the table of permutation for xTetra if needed */
    if ( !pt->xt )
      continue;

    xTetraPerm[pt->xt] = ++nxt;
    pt->xt = nxt;
  }

  /** Count the number of external node/face communicators and initialize it */
  next_node_comm = next_face_comm = 0;
  for ( k=0; k<nprocs; ++k ) {
    if ( shared_pt[k] ) {
      ++next_node_comm;
    }
    if ( shared_face[k] ) {
      ++next_face_comm;
    }
  }

  PMMG_CALLOC(parmesh,parmesh->ext_node_comm,next_node_comm,PMMG_ext_comm,
              "allocate ext_node_comm ",ret_val=PMMG_FAILURE; goto fail_alloc7);
  parmesh->next_node_comm = next_node_comm;

  PMMG_CALLOC(parmesh,parmesh->ext_face_comm,next_face_comm,PMMG_ext_comm,
              "allocate ext_face_comm ",ret_val=PMMG_FAILURE; goto fail_alloc7);
  parmesh->next_face_comm = next_face_comm;

  /** Count internal communicators items before erasing the shared_face array */
  /* Internal node comm */
  nitem_int_node_comm = 0;
  for ( k=1; k<=mesh->np; k++ ) {

    if ( !mesh->point[k].tmp )
      continue;

    for ( j=0; j<nprocs; ++j ) {
      if ( seen_shared_pt[nprocs*(k-1)+j] ) {
        ++nitem_int_node_comm;
        break;
      }
    }
  }

  /* Internal face comm */
  nitem_int_face_comm = 0;
  for ( k=0; k<parmesh->nprocs; k++ ) {
    nitem_int_face_comm += shared_face[k];
  }


  /** External communicators allocation */
  next_node_comm = next_face_comm = 0;
  for ( k=0; k<nprocs; ++k ) {
    /* Node ext comm */
    if ( shared_pt[k] ) {
      pext_node_comm = &parmesh->ext_node_comm[next_node_comm];
      pext_node_comm->color_in  = rank;
      pext_node_comm->color_out = k;
      PMMG_CALLOC(parmesh,pext_node_comm->int_comm_index,shared_pt[k],int,
                  "allocate comm idx",ret_val = PMMG_FAILURE; goto fail_alloc7);
      pext_node_comm->nitem     = shared_pt[k];
      /* Use shared_pt to store the idx of the external communicator me->k */
      shared_pt[k] = next_node_comm++;
    }

    /* Face ext comm */
    if ( shared_face[k] ) {
      pext_face_comm            = &parmesh->ext_face_comm[next_face_comm];
      pext_face_comm->color_in  = rank;
      pext_face_comm->color_out = k;
      PMMG_CALLOC(parmesh,pext_face_comm->int_comm_index,shared_face[k],int,
                  "allocate comm idx",ret_val = PMMG_FAILURE; goto fail_alloc7);
      pext_face_comm->nitem     = shared_face[k];
      /* Use shared_pt to store the idx of the external communicator me->k */
      shared_face[k] = next_face_comm++;
    }
  }

  /** Internal communicators allocation */
  /* Internal node comm */
  assert ( !grp->nitem_int_node_comm );
  PMMG_CALLOC(parmesh,grp->node2int_node_comm_index1,nitem_int_node_comm,int,
              "node2int_node_comm_index1 ",ret_val = PMMG_FAILURE;
              goto fail_alloc7);
  grp->nitem_int_node_comm = nitem_int_node_comm;

  // For the error handling to be complete at this point we have to:
  //   1) free n2i_n_c_idx1
  //   2) reset nitem_int_node_comm, set ret_val, etc
  PMMG_CALLOC(parmesh,grp->node2int_node_comm_index2, nitem_int_node_comm,int,
              "alloc node2int_node_comm_index2 ",
              PMMG_DEL_MEM(parmesh,grp->node2int_node_comm_index1,
                           nitem_int_node_comm,int,
                           "free node2int_node_comm_index1 ");
              grp->nitem_int_node_comm=0;
              ret_val = PMMG_FAILURE;
              goto fail_alloc7);

  /* Internal face comm */
  assert ( !grp->nitem_int_face_comm );
  grp->nitem_int_face_comm = nitem_int_face_comm;
  PMMG_CALLOC(parmesh,grp->node2int_face_comm_index1,nitem_int_face_comm,int,
              "alloc node2int_face_comm_index1 ",ret_val = PMMG_FAILURE;
              goto fail_alloc7);
  // For the error handling to be complete at this point we have to:
  //   1) reverse the reallocation of n2i_n_c_idx1
  //   2) reset nitem_int_face_comm, set ret_val, etc
  PMMG_CALLOC(parmesh,grp->node2int_face_comm_index2,nitem_int_face_comm,int,
              "alloc node2int_face_comm_index2 ",
              PMMG_DEL_MEM(parmesh,grp->node2int_face_comm_index1,
                           nitem_int_face_comm,int,
                           "free node2int_face_comm_index1 ");
              grp->nitem_int_face_comm = 0;
              ret_val = PMMG_FAILURE;
              goto fail_alloc7);

  /** Travel through the mesh and fill the communicators */
  node2int_node_comm_index1 = grp->node2int_node_comm_index1;
  node2int_node_comm_index2 = grp->node2int_node_comm_index2;
  node2int_face_comm_index1 = grp->node2int_face_comm_index1;
  node2int_face_comm_index2 = grp->node2int_face_comm_index2;

  /* Node Communicators */
  /* Idx is used to store the external communicator cursor */
  PMMG_CALLOC(parmesh,idx,parmesh->next_node_comm,int,"allocating idx",
              ret_val = PMMG_FAILURE; goto fail_alloc7);

  i = 0;
  for ( k=1; k<=mesh->np; k++ ) {
    if ( !mesh->point[k].tmp )
      continue;

    inIntComm = 0;
    for ( j=0; j<nprocs; ++j ) {
      pext_node_comm = &parmesh->ext_node_comm[shared_pt[j]];

      if ( seen_shared_pt[nprocs*(k-1)+j] == 1 ) {
        /* Add point in external communicator */
        pext_node_comm->int_comm_index[idx[shared_pt[j]]++] = i;

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
  PMMG_DEL_MEM(parmesh,idx,parmesh->next_node_comm,int,"deallocating idx");

  /* Face Communicators */
  /* Idx is used to store the external communicator cursor */
  PMMG_CALLOC(parmesh,idx,parmesh->next_face_comm,int,"allocating idx",
              ret_val = PMMG_FAILURE; goto fail_alloc7);

  i = 0;
  for ( k=1; k<=mesh->ne; k++ ) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    rankCur = pt->mark;

    for ( ifac=0; ifac<4; ifac++ ) {
      kvois = mesh->adja[4*k-3+ifac]/4;

      if ( (!kvois) || k>kvois ) continue;

      rankVois = part[kvois-1];
      if ( rankCur == rankVois ) continue;
        ifacVois = mesh->adja[4*k-3+ifac]%4;

      if ( rankCur == rank ) {
        /* Add the elt k to communicators */
        pext_face_comm = &parmesh->ext_face_comm[shared_face[rankVois]];
        pext_face_comm->int_comm_index[idx[shared_face[rankVois]]++] = i;

        node2int_face_comm_index1[i] = 4*(pt->base-1) + ifac;
        node2int_face_comm_index2[i] = i;
        ++i;
      }
      else if ( rankVois == rank ) {
        /* Add the elt kvois to communicators */
        ifacVois = mesh->adja[4*k-3+ifac]%4;

        pext_face_comm = &parmesh->ext_face_comm[shared_face[rankCur]];
        pext_face_comm->int_comm_index[idx[shared_face[rankCur]]++] = i;

        node2int_face_comm_index1[i] = 4*(mesh->tetra[kvois].base-1)+ifacVois;
        node2int_face_comm_index2[i] = i;
        ++i;
      }
    }
  }
  PMMG_DEL_MEM(mesh,mesh->adja,4*mesh->nemax+5,int,"dealloc mesh adja");
  PMMG_DEL_MEM(parmesh,idx,parmesh->next_face_comm,int,"deallocating idx");

#warning NIKOS: I need some help with managing error handling here: am I consistently managing the communicator deallocations? are they in a consistent state if an error happens and this returns?
  PMMG_CALLOC(parmesh,parmesh->int_node_comm,1,PMMG_int_comm,
              "allocating int_node_comm",ret_val=PMMG_FAILURE;goto fail_alloc7);
  PMMG_CALLOC(parmesh,parmesh->int_face_comm,1,PMMG_int_comm,
              "allocating int_face_comm",ret_val=PMMG_FAILURE;goto fail_alloc7);

  /* We have 1 Grp per proc, thus : int_*_comm->nitem : nitem_int_*_comm */
  parmesh->int_node_comm->nitem = nitem_int_node_comm;
  parmesh->int_face_comm->nitem = nitem_int_face_comm;

  /** Compact tetrahedra on the proc */
  if ( !PMMG_packTetra(mesh,rank) ) goto fail_alloc7;

  /** Compact xtetra on the proc: in place permutations */
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
  PMMG_DEL_MEM(parmesh,xPointPerm,mesh->xp+1,int,"deallocate metis buffer5");
fail_alloc6:
  PMMG_DEL_MEM(parmesh,xTetraPerm,mesh->xtmax+1,int,"deallocate metis buffer4");
fail_alloc5:
  PMMG_DEL_MEM(parmesh,pointPerm,mesh->np+1,int,"deallocate metis buffer3");
fail_alloc4:
  PMMG_DEL_MEM(parmesh,shared_face,nprocs,int,"deallocate shared_face");
fail_alloc3:
  PMMG_DEL_MEM(parmesh,seen_shared_pt,nprocs*mesh->np,int8_t,
               "deallocate seen_shared_pt");
fail_alloc2:
  PMMG_DEL_MEM(parmesh,shared_pt,nprocs,int,"deallocate shared_pt");
fail_alloc1:
  PMMG_DEL_MEM(parmesh,part,mesh->ne,idx_t,"deallocate metis buffer0");
fail_alloc0:
  return ret_val;
}

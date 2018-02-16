/**
 * \file communicators_pmmg.c
 * \brief Functions related to the communicators construction
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria Soft)
 * \author Nikos Pattakos (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#include "linkedlist_pmmg.h"
#include "coorcell_pmmg.h"

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * \return 1 if success, 0 if fail.
 *
 * Build the node communicators (externals and internals) from the faces ones.
 *
 */
// BUG quand on a plusieurs groupes et plusieurs procs (voir yales2->grid_comm_m.f90).
int PMMG_build_nodeCommFromFaces( PMMG_pParMesh parmesh ) {
  int ier, ier_glob;

  assert ( PMMG_check_extFaceComm ( parmesh ) );
  assert ( PMMG_check_intFaceComm ( parmesh ) );

  /** Build the internal node communicator from the faces ones */
  ier = PMMG_build_intNodeComm(parmesh);
  if ( !ier ) {
    fprintf(stderr,"\n  ## Error: %s: unable to build the internal node"
            " communicators from the internal faces communicators.\n",__func__);
  }
  else {
    /* The internal comm is filled, try to begin to fill the external one */
    assert ( PMMG_check_intNodeComm ( parmesh ) );

    /** Build the external node communicators from the faces ones without looking
     * at the possibility to miss a processor that is not visible from the face
     * communicators: _0_\1/_2_ configuration ( from 0 it is not possible to say
     * that the node at edge intersction belongs to 2 */
    ier = PMMG_build_simpleExtNodeComm(parmesh);
    if ( !ier ) {
      fprintf(stderr,"\n  ## Error: %s: unable to build the simple externals node"
              " communicators from the external faces communicators.\n",__func__);
    }
  }

  /* Check that all steps have successed until here (because the next function
   * involves MPI comms) */
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);
  if ( !ier_glob ) return 0;

  /** Fill the external node communicator */
  ier = PMMG_build_completeExtNodeComm(parmesh);
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);
  if ( !ier ) {
    fprintf(stderr,"\n  ## Error: %s: unable to complete the external node"
            " communicators.\n",__func__);
  }
  if ( !ier_glob ) return 0;

  assert ( PMMG_check_extNodeComm ( parmesh ) );

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * \return 1 if success, 0 if fail.
 *
 * Build the external node communicators from the faces ones without looking
 * at the possibility to miss a processor that is not visible from the face
 * communicators: _0_\1/_2_ configuration ( from 0 it is not possible to say
 * that the node at edge intersction belongs to 2.
 *
 */
int PMMG_build_simpleExtNodeComm( PMMG_pParMesh parmesh ) {
  PMMG_pGrp       grp;
  PMMG_pext_comm  ext_node_comm,ext_face_comm;
  PMMG_pint_comm  int_node_comm,int_face_comm;
  MMG5_pMesh      mesh;
  MMG5_pTetra     pt;
  MMG5_pPoint     ppt;
  int             *face_vertices,*flag,idx1,idx2;
  int             *node2int_node_comm_index1,*node2int_node_comm_index2;
  int             *face2int_face_comm_index1,*face2int_face_comm_index2;
  int             next_face_comm,next_node_comm,nitem_ext_comm;
  int             color_out,ier,i,j,k,iel,ifac,ip,iploc,grpid;

  ier = 0;
  assert ( !parmesh->ext_node_comm && "external comms must be deleted" );

  next_face_comm = parmesh->next_face_comm;
  next_node_comm = next_face_comm;

  PMMG_CALLOC(parmesh,parmesh->ext_node_comm,next_node_comm,PMMG_ext_comm,
              "ext_node_comm ",return 0);
  parmesh->next_node_comm = next_node_comm;

  /** Step 1: For each face, store the position of its vertices in the internal
   * communicator */
  int_face_comm = parmesh->int_face_comm;
  PMMG_CALLOC(parmesh,face_vertices,3*int_face_comm->nitem,int,
              "position of face vertices in int_node_comm",goto end);

  for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
    grp  = &parmesh->listgrp[grpid];
    mesh = grp->mesh;

    for ( k=1; k<=mesh->np; ++k )
      mesh->point[k].tmp = PMMG_UNSET;

    node2int_node_comm_index1 = grp->node2int_node_comm_index1;
    node2int_node_comm_index2 = grp->node2int_node_comm_index2;
    for ( i=0; i<grp->nitem_int_node_comm; ++i ) {
      idx1 = node2int_node_comm_index1[i];
      idx2 = node2int_node_comm_index2[i];

      assert ( mesh->point[idx1].tmp < 0 && "node is not unique in the comm" );
      mesh->point[idx1].tmp = idx2;
    }

    face2int_face_comm_index1 = grp->face2int_face_comm_index1;
    face2int_face_comm_index2 = grp->face2int_face_comm_index2;
    for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
      iel   =  face2int_face_comm_index1[i]/12;
      ifac  = (face2int_face_comm_index1[i]%12)/3;
      iploc = (face2int_face_comm_index1[i]%12)%3;

      idx2 = face2int_face_comm_index2[i];

      pt = &mesh->tetra[iel];
      for ( j=0; j<3; ++j ) {
        ip  = pt->v[_MMG5_idir[ifac][(j+iploc)%3]];
        ppt = &mesh->point[ip];
        assert ( ppt->tmp>=0 && "missing node in the communicator" );

        face_vertices[3*idx2+j] = ppt->tmp;
      }
    }
  }

  /** Step 2: fill the external communicators */
  int_node_comm = parmesh->int_node_comm;
  PMMG_CALLOC(parmesh,flag,int_node_comm->nitem,int,"node flag",goto end);

  for ( k=0; k<next_face_comm; ++k ) {
    ext_face_comm = &parmesh->ext_face_comm[k];
    color_out = ext_face_comm->color_out;

    parmesh->ext_node_comm[k].color_in  = ext_face_comm->color_in;
    parmesh->ext_node_comm[k].color_out = color_out;

    if ( ext_face_comm->nitem )
      /* Reset the node flags */
      for ( j=0; j<int_node_comm->nitem; ++j ) flag[j] = 0;

    /* Count the number of nodes in the communicator */
    nitem_ext_comm = 0;
    for ( i=0; i<ext_face_comm->nitem; ++i ) {
      idx2  = ext_face_comm->int_comm_index[i];

      for ( j=0; j<3; ++j ) {
        ip = face_vertices[3*idx2+j];
        assert ( ip >=0 );

        if ( flag[ip] > 0 ) continue;

        flag[ip] = 1;
        ++nitem_ext_comm;
      }
    }

    /* External communicator allocation */
    PMMG_CALLOC(parmesh,parmesh->ext_node_comm[k].int_comm_index,nitem_ext_comm,
                int,"external node communicator",goto end);
    parmesh->ext_node_comm[k].nitem = nitem_ext_comm;
    ext_node_comm = &parmesh->ext_node_comm[k];

    /* Process the external face communicator and fill the node one */
    if ( ext_face_comm->nitem )
      /* Reset the node flags */
      for ( j=0; j<int_node_comm->nitem; ++j ) flag[j] = 0;

    nitem_ext_comm = 0;

    if ( color_out > parmesh->ext_node_comm[k].color_in ) {
      /* Travel the face vertices in the natural direction */
      for ( i=0; i<ext_face_comm->nitem; ++i ) {
        idx2  = ext_face_comm->int_comm_index[i];

        for  ( j=0; j<3; ++j ) {
          ip = face_vertices[3*idx2+j];
          assert ( ip>=0 );

          if ( flag[ip] > 0 ) continue;
          flag[ip] = 1;
          ext_node_comm->int_comm_index[nitem_ext_comm++] = ip;
        }
      }
    }
    else {
      /* Travel the face vertices in the opposite direction */
      for ( i=0; i<ext_face_comm->nitem; ++i ) {
        idx2  = ext_face_comm->int_comm_index[i];

        for  ( j=0; j<3; ++j ) {
          ip = face_vertices[3*idx2+(3-j)%3];
          assert ( ip>=0 );

          if ( flag[ip] > 0 ) continue;
          flag[ip] = 1;
          ext_node_comm->int_comm_index[nitem_ext_comm++] = ip;
        }
      }
    }
    assert ( nitem_ext_comm == ext_node_comm->nitem );
  }

  /* Success */
  ier = 1;

end:
  if ( !ier ) {
    /* Deallocation of the external comms because we have failed */
    if ( next_node_comm ) {
      for ( k=0; k<next_node_comm; ++k ) {
        PMMG_DEL_MEM(parmesh,parmesh->ext_node_comm[k].int_comm_index,
                     ext_node_comm->nitem,int,"external node communicator");
      }
    }
    PMMG_DEL_MEM(parmesh,parmesh->ext_node_comm,parmesh->next_node_comm
                 ,PMMG_ext_comm, "ext_node_comm ");
  }

  PMMG_DEL_MEM(parmesh,flag,int_node_comm->nitem,int,"node flag");

  PMMG_DEL_MEM(parmesh,face_vertices,3*parmesh->int_face_comm->nitem,int,
               "position of face vertices in int_node_comm");

  return ier;
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * \return 1 if success, 0 if fail.
 *
 * Build the internal node communicators from the faces ones.
 *
 */
int PMMG_build_intNodeComm( PMMG_pParMesh parmesh ) {
  PMMG_pGrp       grp;
  PMMG_coorCell   *coor_list;
  MMG5_pMesh      mesh;
  MMG5_pTetra     pt;
  MMG5_pPoint     ppt;
  double          bb_min[3],bb_max[3],delta,dist[3],dd;
  int             *face2int_face_comm_index1;
  int             *node2int_node_comm_index1,*node2int_node_comm_index2;
  int             *shared_fac,*new_pos,nitem_node,first_nitem_node,pos;
  int             *face_vertices,ier,i,j,iel,ifac,ip,iploc,grpid,idx,fac_idx;
  int8_t          update;

  ier = 0;

  /** Step 1: give a unique position in the internal communicator for each mesh
   * point but don't care about the unicity of the position for a point shared
   * by multiple groups */
  assert ( !parmesh->int_node_comm->nitem );
  nitem_node = 0;

  PMMG_MALLOC(parmesh,shared_fac,parmesh->int_face_comm->nitem,int,
              "Faces shared by 2 groups",goto end);
  for ( i=0; i<parmesh->int_face_comm->nitem; ++i )
    shared_fac[i] = PMMG_UNSET;

  for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
    /* Reset the flag field of the nodes of the group */
    grp  = &parmesh->listgrp[grpid];
    mesh = grp->mesh;

    for ( ip=1; ip<=mesh->np; ++ip )
      mesh->point[ip].tmp  = PMMG_UNSET;

    face2int_face_comm_index1 = grp->face2int_face_comm_index1;

    first_nitem_node = nitem_node;
    for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
      idx  = face2int_face_comm_index1[i];
      iel  =  idx/12;
      ifac = (idx%12)/3;

      assert ( iel && iel<=mesh->ne );
      assert ( 0<=ifac && ifac<4 );

      pt    = &mesh->tetra[iel];
      for ( j=0; j<3; ++j ) {
        ip  = pt->v[_MMG5_idir[ifac][j]];

        assert ( ip && ip<=mesh->np );

        ppt = &mesh->point[ip];
        if ( ppt->tmp < 0 )
          /** Give a position int the internal communicator to this point */
          ppt->tmp = nitem_node++;
      }

      /* Mark the face as seen */
      ++shared_fac[grp->face2int_face_comm_index2[i]];
    }

    /* Allocations of the node2int_node arrays */
    PMMG_CALLOC(parmesh,grp->node2int_node_comm_index1,nitem_node-first_nitem_node,
                int,"node2int_node_comm_index1",goto end);

    grp->nitem_int_node_comm = nitem_node-first_nitem_node;
    PMMG_CALLOC(parmesh,grp->node2int_node_comm_index2,nitem_node-first_nitem_node,
                int,"node2int_node_comm_index2",goto end);

    /* Fill this arrays */
    node2int_node_comm_index1 = grp->node2int_node_comm_index1;
    node2int_node_comm_index2 = grp->node2int_node_comm_index2;
    idx = 0;
    for ( ip=1; ip<=mesh->np; ++ip ) {
      ppt = &mesh->point[ip];
      if ( ppt->tmp < 0 ) continue;

      node2int_node_comm_index1[idx] = ip;
      node2int_node_comm_index2[idx] = first_nitem_node+idx;
      ++idx;
    }
    assert ( idx == nitem_node-first_nitem_node && "missing item");
  }


  /** Step 2: remove some of the multiple positions and pack communicators */
  PMMG_MALLOC(parmesh,new_pos,nitem_node,int,"new pos in int_node_comm",goto end);
  PMMG_MALLOC(parmesh,face_vertices,3*parmesh->int_face_comm->nitem,int,
              "pos of face vertices in int_node_comm",goto end);

  /* Initialisation of the new_pos array with the current position */
 for ( i=0; i<nitem_node; ++i )
    new_pos[i] = i;

 /* Initialisation of the face_vertices array */
 for ( i=0; i<3*parmesh->int_face_comm->nitem; ++i )
   face_vertices[i]   = PMMG_UNSET;


 PMMG_CALLOC(parmesh,coor_list,nitem_node,PMMG_coorCell,"node coordinates",
             goto end);

#ifndef NDEBUG
 /* Store the node coordinates of the nodes */
 for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
   grp  = &parmesh->listgrp[grpid];
   mesh = grp->mesh;
   assert ( mesh->info.delta && "unable to unscale the mesh");

   for ( i=0; i<grp->nitem_int_node_comm; ++i ) {
     ip      = grp->node2int_node_comm_index1[i];
     idx     = grp->node2int_node_comm_index2[i];
     for ( j=0; j<3; ++j )
       coor_list[idx].coor[j] = mesh->info.delta*mesh->point[ip].c[j]+mesh->info.min[j];
   }
 }
 /* Scale the coordinates depending to the bounding box ofthe internal comm */
 if ( !PMMG_scale_coorCellList(coor_list,nitem_node,bb_min,bb_max,&delta) )
   goto end;

#endif


  do {
    update = 0;
    /* Mark all the interface faces as not seen */
    for ( i=0; i<parmesh->int_face_comm->nitem; ++i )
      shared_fac[i] = -abs(shared_fac[i]);

    for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
      grp  = &parmesh->listgrp[grpid];
      mesh = grp->mesh;

      /* Reset tmp to the old position in the communicator */
      for ( i=0; i<grp->nitem_int_node_comm; ++i ) {
        ip = grp->node2int_node_comm_index1[i];
        mesh->point[ip].tmp  = grp->node2int_node_comm_index2[i];;
      }

      for ( i=0; i<grp->nitem_int_face_comm; ++i ) {
        idx   = grp->face2int_face_comm_index1[i];
        iel   =  idx/12;
        ifac  = (idx%12)/3;
        iploc = (idx%12)%3;

        fac_idx = grp->face2int_face_comm_index2[i];

        pt = &mesh->tetra[iel];
        if ( shared_fac[fac_idx]<=0 ) {
          /* First time that we see this face: store it in the natural direction
           * and mark it as seen */
          shared_fac[fac_idx] = abs(shared_fac[fac_idx]);
          for ( j=0; j<3; ++j ) {
            ip  = pt->v[_MMG5_idir[ifac][j]];
            ppt = &mesh->point[ip];
            assert ( ppt->tmp>=0 );

            pos = 3*fac_idx+j;
            if ( face_vertices[pos] < 0 ) {
              face_vertices[pos] = new_pos[mesh->point[ip].tmp];
              update = 1;
            }
            else if ( new_pos[face_vertices[pos]]>new_pos[mesh->point[ip].tmp] ) {
              new_pos[face_vertices[pos]] = new_pos[mesh->point[ip].tmp];
              update = 1;
            }
            else if ( new_pos[face_vertices[pos]]<new_pos[mesh->point[ip].tmp] ) {
              new_pos[mesh->point[ip].tmp] = new_pos[face_vertices[pos]];
              update = 1;
            }

#ifndef NDEBUG
            double scaled_coor[3];

            /* Unscale the corrdinates that has been scaled by the mesh */
            for ( j=0; j<3; ++j )
              scaled_coor[j] = mesh->info.delta*mesh->point[ip].c[j]
                + mesh->info.min[j];

            /* Scale it using the internal node comm scaling data */
            assert ( delta );
            dd = 1./delta;
            for ( j=0; j<3; ++j )
              scaled_coor[j] = dd*(scaled_coor[j]-bb_min[j]);

            /* Compute the distance between the points */
            dd = 0;
            for ( j=0; j<3; ++j ) {
              dist[j] = coor_list[new_pos[mesh->point[ip].tmp]].coor[j]
                - scaled_coor[j];
              dd += dist[j]*dist[j];
            }
            assert ( dd<_MMG5_EPSD );
#endif
          }
        }
        else {
          /* Second time we see this face: store it in the opposite direction */
          for ( j=0; j<3; ++j ) {
            ip  = pt->v[_MMG5_idir[ifac][(iploc+3-j)%3]];
            ppt = &mesh->point[ip];
            assert ( ppt->tmp>=0 );

            pos = 3*fac_idx+j;
            if ( face_vertices[pos] < 0 ) {
              face_vertices[pos] = new_pos[mesh->point[ip].tmp];
              update = 1;
            }
            else if ( new_pos[face_vertices[pos]]>new_pos[mesh->point[ip].tmp] ) {
              new_pos[face_vertices[pos]] = new_pos[mesh->point[ip].tmp];
              update = 1;
            }
            else if ( new_pos[face_vertices[pos]]<new_pos[mesh->point[ip].tmp] ) {
              new_pos[mesh->point[ip].tmp] = new_pos[face_vertices[pos]];
              update = 1;
            }
#ifndef NDEBUG
            double scaled_coor[3];

            /* Unscale the corrdinates that has been scaled by the mesh */
            for ( j=0; j<3; ++j )
              scaled_coor[j] = mesh->info.delta*mesh->point[ip].c[j]
                + mesh->info.min[j];

            /* Scale it using the internal node comm scaling data */
            assert ( delta );
            dd = 1./delta;
            for ( j=0; j<3; ++j )
              scaled_coor[j] = dd*(scaled_coor[j]-bb_min[j]);

            /* Compute the distance between the points */
            dd = 0;
            for ( j=0; j<3; ++j ) {
              dist[j] = coor_list[new_pos[mesh->point[ip].tmp]].coor[j]
                -scaled_coor[j];
              dd += dist[j]*dist[j];
            }
            assert ( dd<_MMG5_EPSD );
#endif
          }
        }
      }
    }
  } while ( update );

  /* Update node2int_node_comm arrays */
  for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
    grp  = &parmesh->listgrp[grpid];

    for ( i=0; i<grp->nitem_int_node_comm; ++i ) {
      idx = grp->node2int_node_comm_index2[i];
      grp->node2int_node_comm_index2[i] = new_pos[idx];
    }
  }

  /* Repove the empty positions in int_node_comm */
  for ( i=0; i<nitem_node; ++i )
    new_pos[i] = PMMG_UNSET;

  j = 0;
  for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
    grp  = &parmesh->listgrp[grpid];

    for ( i=0; i<grp->nitem_int_node_comm; ++i ) {

      idx = grp->node2int_node_comm_index2[i];
      if ( new_pos[idx]<0 ) new_pos[idx] = j++;

      grp->node2int_node_comm_index2[i] = new_pos[idx];
    }
  }
  nitem_node = j;

  /** Step 3: Remove the last doublon by comparing the listed coordinates */
  /* Store the node coordinates of the nodes */
  for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
    grp  = &parmesh->listgrp[grpid];
    mesh = grp->mesh;
    assert ( mesh->info.delta && "unable to unscale the mesh");

    for ( i=0; i<grp->nitem_int_node_comm; ++i ) {
      ip      = grp->node2int_node_comm_index1[i];
      idx     = grp->node2int_node_comm_index2[i];
      for ( j=0; j<3; ++j )
        coor_list[idx].coor[j] = mesh->info.delta*mesh->point[ip].c[j]+mesh->info.min[j];
    }
  }
  /* Scale the coordinates depending to the bounding box ofthe internal comm */
  if ( !PMMG_scale_coorCellList(coor_list,nitem_node,bb_min,bb_max,&delta) )
    goto end;

  /* Store the point position in the internal communicator */
  for ( i=0; i<nitem_node; ++i ) {
    new_pos[i]       = i;
    coor_list[i].pos = i;
  }

  /* Sort coor_list depending on its coordinates */
  qsort(coor_list,nitem_node,sizeof(PMMG_coorCell),PMMG_compare_coorCell);

  /* Travel the list and remove the identic nodes */
  idx = 0;
  new_pos[coor_list[0].pos] = 0;
  for ( i=1; i<nitem_node; ++i ) {
    if ( PMMG_compare_coorCell(&coor_list[i],&coor_list[idx]) ) {
      ++idx;
      if ( idx != i ) {
        coor_list[idx] = coor_list[i];
      }
      new_pos[coor_list[i].pos] = idx;
    }
    else
      new_pos[coor_list[i].pos] = new_pos[coor_list[idx].pos];
  }
  nitem_node = idx+1;

  /* Update node2int_node_comm arrays */
  for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
    grp  = &parmesh->listgrp[grpid];

    for ( i=0; i<grp->nitem_int_node_comm; ++i ) {

      idx = grp->node2int_node_comm_index2[i];
      grp->node2int_node_comm_index2[i] = new_pos[idx];
    }
  }

  /* Repove the empty positions in int_node_comm */
  for ( i=0; i<nitem_node; ++i )
    new_pos[i] = PMMG_UNSET;

  j = 0;
  for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
    grp  = &parmesh->listgrp[grpid];

    for ( i=0; i<grp->nitem_int_node_comm; ++i ) {

      idx = grp->node2int_node_comm_index2[i];
      if ( new_pos[idx]<0 ) new_pos[idx] = j++;

      grp->node2int_node_comm_index2[i] = new_pos[idx];
    }
  }
  nitem_node = j;

  /** Step 4: Update the number of items in the internal node communicator */
  parmesh->int_node_comm->nitem = nitem_node;

  /* Success */
  ier = 1;

end:
  if ( !ier ) {
    /* Dealloc of the communicators because we have failed */
    for ( grpid=0; grpid<parmesh->ngrp; ++grpid ) {
      grp  = &parmesh->listgrp[grpid];

      PMMG_DEL_MEM(parmesh,grp->node2int_node_comm_index1,
                   grp->nitem_int_node_comm,int,"node2int_node_comm_index1");
      PMMG_DEL_MEM(parmesh,grp->node2int_node_comm_index2,
                   grp->nitem_int_node_comm,int,"node2int_node_comm_index2");
    }
  }

  PMMG_DEL_MEM(parmesh,new_pos,nitem_node,int,"new pos in int_node_comm");
  PMMG_DEL_MEM(parmesh,face_vertices,3*parmesh->int_face_comm->nitem,int,
              "pos of face vertices in int_node_comm");
  PMMG_DEL_MEM(parmesh,shared_fac,parmesh->int_face_comm->nitem,int,
               "Faces shared by 2 groups");

  PMMG_DEL_MEM(parmesh,coor_list,nitem_node,PMMG_coorCell,"node coordinates");

  return ier;
}

/**
 * \param parmesh pointer toward a parmesh structure
 *
 * \return 1 if success, 0 if fail.
 *
 * Complete the external communicators by travelling through the processors and
 * faces to detect all the processors to which each node belongs.
 *
 */
int PMMG_build_completeExtNodeComm( PMMG_pParMesh parmesh ) {
  PMMG_pext_comm  ext_node_comm,*comm_ptr;
  PMMG_pint_comm  int_node_comm;
  PMMG_lnkdList   **proclists,list;
  int             *intvalues,nitem,nproclists,ier,ier2,k,i,j,idx,pos,rank,color;
  int             *itosend,*itorecv,*i2send_size,*i2recv_size,nitem2comm;
  int             *nitem_ext_comm,next_comm,val1_i,val2_i,val1_j,val2_j;;
  int8_t          glob_update,loc_update;
  MPI_Request     *request;
  MPI_Status      *status;

  ier = 0;

  rank          = parmesh->myrank;
  int_node_comm = parmesh->int_node_comm;
  nitem         = int_node_comm->nitem;

  PMMG_CALLOC(parmesh,int_node_comm->intvalues,nitem,int,"node communicator",
    return 0);
  intvalues     = int_node_comm->intvalues;

  /** Use intvalues to flag the already treated points of the communicator:
   * initialization to 0.  */
  for ( k=0; k<nitem; ++k ) intvalues[k] = 0;

  proclists = NULL;
  PMMG_CALLOC(parmesh,proclists,nitem,PMMG_lnkdList*,"array of linked lists",
              goto end);

  /* Reallocation of the list of external comms at maximal size (nprocs) to
   * avoid tricky treatment when filling it.*/
  PMMG_REALLOC(parmesh,parmesh->ext_node_comm,parmesh->nprocs,
               parmesh->next_node_comm,PMMG_ext_comm,
               "list of external communicators",goto end);
  next_comm = parmesh->next_node_comm;
  parmesh->next_node_comm = parmesh->nprocs;

  for ( k=next_comm; k<parmesh->nprocs; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    ext_node_comm->nitem          = 0;
    ext_node_comm->int_comm_index = NULL;
    ext_node_comm->color_in       = rank;
    ext_node_comm->color_out      = -1;
  }

  comm_ptr = NULL;
  PMMG_CALLOC(parmesh,comm_ptr,parmesh->nprocs,PMMG_pext_comm,
              "array of pointers toward the external communicators",
              goto end);

  /** Step 1: initialization of the list of the procs to which a point belongs
   * by the value of the current mpi rank */
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];

    if ( !ext_node_comm->nitem ) continue;

    comm_ptr[ext_node_comm->color_out] = ext_node_comm;

    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      idx = ext_node_comm->int_comm_index[i];

      if ( intvalues[idx] ) continue;

      PMMG_CALLOC(parmesh,proclists[idx],1,PMMG_lnkdList,"linked list pointer",
                  goto end);
      if ( !PMMG_lnkdListNew(parmesh,proclists[idx],PMMG_LISTSIZE) ) goto end;

      if ( !PMMG_add_cell2lnkdList( parmesh,proclists[idx],rank,idx ) )
        goto end;

      intvalues[idx] = 1;
    }
  }
  /* Fill the missing pointer toward the empty external communicators */
  for ( k=0; k<parmesh->nprocs; ++k ) {
    if ( !comm_ptr[k] ) {
      comm_ptr[k] = &parmesh->ext_node_comm[next_comm++];
      comm_ptr[k]->color_out = k;
    }
  }

  /** Step 2: While at least the proc list of 1 node is modified, send and
   * recieve the proc list of all the nodes to/from the other processors. At the
   * end of this loop, each node has the entire list of the proc to which it
   * belongs */
  PMMG_MALLOC(parmesh,request,parmesh->next_node_comm,MPI_Request,
              "mpi request array",goto end);

  PMMG_MALLOC(parmesh,status,parmesh->next_node_comm,MPI_Status,
              "mpi status array",goto end);

  PMMG_CALLOC(parmesh,i2send_size,parmesh->next_node_comm,int,
              "size of the i2send array",goto end);
  PMMG_CALLOC(parmesh,i2recv_size,parmesh->next_node_comm,int,
              "size of the i2recv array",goto end);

  if ( !PMMG_lnkdListNew(parmesh,&list,PMMG_LISTSIZE) ) goto end;

  do {
    glob_update = loc_update = 0;

    /** Send the list of procs to which belong each point of the communicator */
    for ( k=0; k<parmesh->next_node_comm; ++k ) {

      request[k] = MPI_REQUEST_NULL;
      ext_node_comm = &parmesh->ext_node_comm[k];

      /* Computation of the number of data to send to the other procs (we want
       * to send the the size of the linked list and the val1 and val2 fields of
       * each cell ) */
      nitem2comm = 0;
      for ( i=0; i<ext_node_comm->nitem; ++i ) {
        idx         = ext_node_comm->int_comm_index[i];
        nitem2comm += proclists[idx]->nitem*2+1;
      }

      if ( i2send_size[k] < nitem2comm ) {
        PMMG_REALLOC(parmesh,ext_node_comm->itosend,nitem2comm,i2send_size[k],int,
                     "itosend",goto end);
        i2send_size[k] = nitem2comm;
      }

      /* Filling of the array to send */
      pos     = 0;
      itosend = ext_node_comm->itosend;
      color   = ext_node_comm->color_out;
      for ( i=0; i<ext_node_comm->nitem; ++i ) {
        idx  = ext_node_comm->int_comm_index[i];
        pos += PMMG_packInArray_lnkdList(proclists[idx],&itosend[pos]);
      }
      assert ( pos==nitem2comm );
      MPI_CHECK( MPI_Isend(itosend,nitem2comm,MPI_INT,color,0,parmesh->comm,
                           &request[color]),goto end );
    }

    /** Recv the list of procs to which belong each point of the communicator */
    for ( k=0; k<parmesh->next_node_comm; ++k ) {
      ext_node_comm = &parmesh->ext_node_comm[k];

      if ( !ext_node_comm->nitem ) continue;

      color         = ext_node_comm->color_out;

      MPI_CHECK( MPI_Probe(color,0,parmesh->comm,&status[0] ),goto end);
      MPI_CHECK( MPI_Get_count(&status[0],MPI_INT,&nitem2comm),goto end);

      if ( i2recv_size[k] < nitem2comm ) {
        PMMG_REALLOC(parmesh,ext_node_comm->itorecv,nitem2comm,i2recv_size[k],
                     int,"itorecv",goto end);
        i2recv_size[k] = nitem2comm;
      }
      itorecv       = ext_node_comm->itorecv;
      MPI_CHECK( MPI_Recv(itorecv,nitem2comm,MPI_INT,color,0,parmesh->comm,
                          &status[0]), goto end );

      pos     = 0;
      for ( i=0; i<ext_node_comm->nitem; ++i ) {
        idx  = ext_node_comm->int_comm_index[i];
        assert ( idx>=0 );

        PMMG_reset_lnkdList( parmesh,&list );
        ier2 = PMMG_unpackArray_inLnkdList(parmesh,&list,&itorecv[pos]);

        if ( ier2 < 0 ) goto end;
        pos += ier2;

        ier2 = PMMG_merge_lnkdList(parmesh,proclists[idx],&list);
        if ( !ier2 ) goto end;

        loc_update |= (ier2%2);
      }
    }

    MPI_CHECK( MPI_Waitall(parmesh->next_node_comm,request,status), goto end );
    MPI_CHECK( MPI_Allreduce(&loc_update,&glob_update,1,MPI_INT8_T,MPI_LOR,
                             parmesh->comm),goto end);

  } while ( glob_update );

  /** Step 3: Cancel the old external communicator and build it again from the
   * list of proc of each node */
  PMMG_CALLOC(parmesh,nitem_ext_comm,parmesh->nprocs,int,
              "number of items in each external communicator",goto end);

  /* Remove the empty proc lists */
  nproclists = 0;
  for ( i=0; i<nitem; ++i ) {
    if ( intvalues[i] )
      proclists[nproclists++] = proclists[i];
  }

  /* Sort the list of procs to which each node belong to ensure that we will
   * build the external communicators in the same order on both processors
   * involved in an external comm */
  qsort(proclists,nproclists,sizeof(PMMG_lnkdList*),PMMG_compare_lnkdList);

  /* Remove the non unique paths */
  if ( nproclists ) {
    idx = 0;
    for ( i=1; i<nproclists; ++i ) {
      assert ( proclists[i]->nitem );
      if ( PMMG_compare_lnkdList(&proclists[i],&proclists[idx]) ) {
        ++idx;
        if ( idx != i ) {
          proclists[idx] = proclists[i];
        }
      }
    }
    nproclists = idx+1;
  }

  /* Double loop over the list of procs to which belong each point to add the
   * couples to the suitable external communicators */
  for ( k=0; k<nproclists; ++k ) {
    for ( i=0; i<proclists[k]->nitem; ++i ) {
      val1_i = proclists[k]->item[i].val1;
      val2_i = proclists[k]->item[i].val2;

      for ( j=i+1; j<proclists[k]->nitem; ++j ) {
        val1_j = proclists[k]->item[j].val1;
        val2_j = proclists[k]->item[j].val2;

        assert ( val1_i != val1_j );

        if ( val1_i == rank ) {
          ext_node_comm = comm_ptr[val1_j];
          assert ( ext_node_comm );

          if ( !ext_node_comm->nitem ) continue;

          if ( nitem_ext_comm[val1_j] == ext_node_comm->nitem ) {
            /* Reallocation */
            PMMG_REALLOC(parmesh,ext_node_comm->int_comm_index,
                         (int)((1.+PMMG_GAP)*ext_node_comm->nitem)+1,
                         ext_node_comm->nitem,int,
                         "external communicator",goto end);
            ext_node_comm->nitem = (int)((1.+PMMG_GAP)*ext_node_comm->nitem)+1;
            comm_ptr[val1_j] = ext_node_comm;
          }
          ext_node_comm->int_comm_index[nitem_ext_comm[val1_j]++] = val2_i;
        }
        else if ( val1_j == rank ) {
          ext_node_comm = comm_ptr[val1_i];
          assert ( ext_node_comm );
          if ( nitem_ext_comm[val1_j] == ext_node_comm->nitem ) {
            /* Reallocation */
            PMMG_REALLOC(parmesh,ext_node_comm->int_comm_index,
                          (int)((1.+PMMG_GAP)*ext_node_comm->nitem)+1,
                         ext_node_comm->nitem,int,
                         "external communicator",goto end);
            ext_node_comm->nitem = (int)((1.+PMMG_GAP)*ext_node_comm->nitem)+1;
            comm_ptr[val1_i] = ext_node_comm;
          }
          ext_node_comm->int_comm_index[nitem_ext_comm[val1_i]++] = val2_j;
        }
      }
    }
  }

  /* Pack and reallocate the external communicators */
  // The pack may be removed if we allow to have always nprocs external comms
  next_comm = 0;
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];
    if ( !ext_node_comm->nitem ) {
      /* Empty communicator */
      continue;
    }
    else {
      assert ( ext_node_comm->color_out>=0 && ext_node_comm->color_out!=rank );
      assert ( nitem_ext_comm[ext_node_comm->color_out] );

      PMMG_REALLOC(parmesh,ext_node_comm->int_comm_index,
                   nitem_ext_comm[ext_node_comm->color_out],
                   ext_node_comm->nitem,int,
                   "external communicator",goto end);
      ext_node_comm->nitem = nitem_ext_comm[ext_node_comm->color_out];

      if ( next_comm != k )
        parmesh->ext_node_comm[next_comm] = *ext_node_comm;
    }
    ++next_comm;
  }
  PMMG_REALLOC(parmesh,parmesh->ext_node_comm,
               next_comm,parmesh->next_node_comm,PMMG_ext_comm,
               "list of external communicator",goto end);
  parmesh->next_node_comm = next_comm;

  /* Success */
  ier = 1;

end:
  PMMG_DEL_MEM(parmesh,int_node_comm->intvalues,nitem,int,"node communicator");

  if ( proclists ) {
    for ( k=0; k<nitem; ++k ) {
      if ( proclists[k] ) {
        PMMG_DEL_MEM(parmesh,proclists[k]->item,proclists[k]->nitem_max,
                     PMMG_lnkdCell,"linked list array");

      }
    }
    PMMG_DEL_MEM(parmesh,proclists,nitem,PMMG_lnkdList*,"array of linked lists");
  }
  PMMG_DEL_MEM(parmesh,comm_ptr,parmesh->nprocs,PMMG_pext_comm,
              "array of pointers toward the external communicators");

  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];

    // Change this and add to the external comm the possibility to not
    // unalloc/realloc every time, thus, here, we will be able to reset the
    // communicators without unallocated it
    PMMG_DEL_MEM(parmesh,ext_node_comm->itosend,i2send_size[k],int,"i2send");
    PMMG_DEL_MEM(parmesh,ext_node_comm->itorecv,i2recv_size[k],int,"i2recv");
  }

  PMMG_DEL_MEM(parmesh,request,parmesh->next_node_comm,MPI_Request,
               "mpi request array");
  PMMG_DEL_MEM(parmesh,status,parmesh->next_node_comm,MPI_Status,
               "mpi status array");

  PMMG_DEL_MEM(parmesh,i2send_size,parmesh->next_node_comm,int,
               "size of the i2send array");

  PMMG_DEL_MEM(parmesh,i2recv_size,parmesh->next_node_comm,int,
               "size of the i2recv array");

  PMMG_DEL_MEM(parmesh,nitem_ext_comm,parmesh->nprocs,int,
               "number of items in each external communicator");

  PMMG_DEL_MEM(parmesh,list.item,list.nitem_max,PMMG_lnkdCell,
               "linked list array");

  return ier;
}

/**
 * \file libparmmg1.c
 * \brief Wrapper for the parallel remeshing library.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (InriaSoft)
 * \author Nikos Pattakos (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Internal function that perform the parallel remeshing.
 *
 */
#include "parmmg.h"

/**
 * \param grp pointer toward the group in which we want to update the list of
 * nodes that are in the internal communicator.
 *
 * \return 1 if success, 0 otherwise
 *
 * Update the list of vertices indices with the pack point index stored in the
 * tmp field of points.
 *
 */
int PMMG_update_node2intPackedVertices( PMMG_pGrp grp ) {
  MMG5_pPoint ppt;
  int         *node2int_node_comm_index1;
  int         k,iadr;

  node2int_node_comm_index1 = grp->node2int_node_comm_index1;

  for (k=0; k<grp->nitem_int_node_comm; ++k) {
    iadr = node2int_node_comm_index1[k];

    ppt = &grp->mesh->point[iadr];
    assert ( MG_VOK(ppt) );

    node2int_node_comm_index1[k] = ppt->tmp;
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure (unused).
 * \param ne pointer toward the number of packed tetra
 *
 * \return 1 if success, 0 if fail.
 *
 * Count the number of packed tetra and store the packed tetra index in flag.
 *
 */

int PMMG_mark_packedTetra(MMG5_pMesh mesh,int *ne) {
  MMG5_pTetra   pt;
  int           k;

  (*ne) = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;
    pt->flag = ++(*ne);
  }
  return 1;
}

/**
 * \param grp pointer toward the group in which we want to update the list of
 * faces that are in the internal communicator.
 *
 * \return 1 if success, 0 otherwise
 *
 * Update the list of tetra indices with the pack tetra index stored in the
 * flag field of tetras.
 *
 */
int PMMG_update_face2intPackedTetra( PMMG_pGrp grp ) {
  MMG5_pTetra pt;
  int         *face2int_face_comm_index1;
  int         k,iel,ifac,iploc;

  face2int_face_comm_index1 = grp->face2int_face_comm_index1;

  for (k=0; k<grp->nitem_int_face_comm; ++k) {
    iel   =  face2int_face_comm_index1[k]/12;
    ifac  = (face2int_face_comm_index1[k]%12)/3;
    iploc = (face2int_face_comm_index1[k]%12)%3;

    pt = &grp->mesh->tetra[iel];
    assert ( MG_EOK(pt) && pt->flag );

    face2int_face_comm_index1[k] = 12*pt->flag+3*ifac+iploc;
  }

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param igrp index of the group that we want to treat
 *
 * \return 1 if success, 0 if fail.
 *
 * Pack the tetra+adja arrays and update the face communicator with the new
 * tetra indices.
 *
 */
int PMMG_packTetra( PMMG_pParMesh parmesh, int igrp ) {
  PMMG_pGrp   grp;
  MMG5_pMesh  mesh;
  int         ne;

  assert ( igrp < parmesh->ngrp );

  grp                       = &parmesh->listgrp[igrp];
  mesh                      = grp->mesh;

  if ( !mesh ) return 1;

  /** Store in flag the pack index of each tetra */
  if ( !PMMG_mark_packedTetra(mesh,&ne) ) return 0;

  /** Update the tetra indices in the face communicator */
  if ( !PMMG_update_face2intPackedTetra(grp) ) return 0;

  /* compact tetrahedra */
  if ( mesh->adja ) {
    if ( !MMG3D_pack_tetraAndAdja(mesh) ) return 0;
  }
  else {
    if ( !MMG3D_pack_tetra(mesh) ) return 0;
  }
  assert ( ne==mesh->ne );

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 otherwise
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
  int         np,nc,igrp;

  for ( igrp=0; igrp<parmesh->ngrp; ++igrp ) {
    grp                       = &parmesh->listgrp[igrp];
    mesh                      = grp->mesh;
    met                       = grp->met;
    disp                      = grp->disp;

    if ( !mesh ) continue;

    /* Pack tetrahedra */
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

    /** Store in tmp the pack index of each point and count the corner*/
    if ( !MMG3D_mark_packedPoints(mesh,&np,&nc) ) return 0;

    /* node index update in internal communicator */
    if ( !PMMG_update_node2intPackedVertices( grp ) ) return 0;

    /** Update the element vertices indices */
    if ( !MMG3D_update_eltsVertices(mesh) ) return 0;

    if ( MMG3D_pack_pointArray(mesh) < 0 ) return 0;

    /* create prism adjacency */
    if ( !MMG3D_hashPrism(mesh) ) {
      fprintf(stderr,"\n  ## Error: %s: prism hashing problem. Exit program.\n",
              __func__);
      return 0;
    }

    /* Remove the MG_REQ tags added by the nosurf option */
    MMG3D_unset_reqBoundaries(mesh);

    if ( mesh->info.imprim > PMMG_VERB_VERSION ) {
      fprintf(stdout,"     NUMBER OF VERTICES   %8d   CORNERS %8d\n",mesh->np,nc);
      fprintf(stdout,"     NUMBER OF ELEMENTS   %8d\n",mesh->ne);
    }

    /* to could save the mesh, the adjacency have to be correct */
    if ( mesh->info.ddebug ) {
      if ( (!mesh->adja) && !MMG3D_hashTetra(mesh,1) ) {
        fprintf(stderr,"\n  ## Error: %s: tetra hashing problem. Exit program.\n",
                __func__);
        return 0;
      }
      if ( !MMG5_chkmsh(mesh,1,1) ) {
        fprintf(stderr,"  ##  Problem. Invalid mesh.\n");
        return 0;
      }
    }
  }

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param igrp index of the group that we want to treat
 * \param facesData pointer toward the allocatable list of the node indices of
 * the interface faces present in the list of interface triangles of the group.
 *
 * \return 1 if success, 0 if fail.
 *
 * Store in \a facesVertices the nodes indices of the interface faces.  The \f$
 * i^th \f$ face is stored at the \f$ [ 3*i;3$i+2 ] positions of the \a
 * facesData array. The \a facesData array is allocated in this function.
 *
 */
static inline
int PMMG_store_faceVerticesInIntComm( PMMG_pParMesh parmesh, int igrp,
                                      int **facesData) {
  PMMG_pGrp   grp;
  MMG5_pMesh  mesh;
  MMG5_pTetra pt;
  int         *face2int_face_comm_index1,nitem_int_face_comm;
  int         iel,ifac,iploc,ia,ib,ic;
  int         k;

  assert ( igrp < parmesh->ngrp );

  grp                 = &parmesh->listgrp[igrp];
  nitem_int_face_comm = grp->nitem_int_face_comm;

  PMMG_MALLOC(parmesh,*facesData,3*nitem_int_face_comm,int,"facesData",return 0);

  face2int_face_comm_index1 = grp->face2int_face_comm_index1;
  mesh                      = parmesh->listgrp[igrp].mesh;
  for ( k=0; k<nitem_int_face_comm; ++k ) {
    /** Get the vertices indices of the interface triangles */
    iel   =  face2int_face_comm_index1[k]/12;
    ifac  = (face2int_face_comm_index1[k]%12)/3;
    iploc = (face2int_face_comm_index1[k]%12)%3;


    pt = &mesh->tetra[iel];

    assert( MG_EOK(pt) );
    assert( pt->xt && ( mesh->xtetra[pt->xt].ftag[ifac] & MG_PARBDY ) );

    ia = pt->v[MMG5_idir[ifac][iploc]];
    ib = pt->v[MMG5_idir[ifac][(iploc+1)%3]];
    ic = pt->v[MMG5_idir[ifac][(iploc+2)%3]];

    /** Store the face vertices */
    (*facesData)[3*k]   = ia;
    (*facesData)[3*k+1] = ib;
    (*facesData)[3*k+2] = ic;
  }

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure.
 * \param igrp index of the group that we want to treat
 * \param facesData list the node vertices of the interface faces
 * present in the list of interface triangles of the group.
 *
 * \return 1 if success, 0 if fail.
 *
 * Find the index of the interface tetras from the data stored in the
 * \a facesData array (by the \ref store_faceVerticesInIntComm function) and
 * update the face2int_face_comm_index1 array. This function frees the \a
 * facesData array.
 *
 */
static inline
int  PMMG_update_face2intInterfaceTetra( PMMG_pParMesh parmesh, int igrp,
                                         int *facesData ) {
  PMMG_pGrp    grp;
  MMG5_pMesh   mesh;
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_Hash   hash;
  int          *face2int_face_comm_index1,nitem;
  int          hashVal,iel,ifac,iploc,ia,ib,ic;
  int          k,i,ier;

  assert ( igrp < parmesh->ngrp );

  grp   = &parmesh->listgrp[igrp];
  nitem = grp->nitem_int_face_comm;
  ier     = 1;

  /** Step 1: Hash the MG_PARBDY faces */
  mesh = parmesh->listgrp[igrp].mesh;
  if ( !MMG5_hashNew(mesh,&hash,0.51*nitem,1.51*nitem) ) {
    ier = 0;
    goto facesData;
  }

  for ( k=1; k<=mesh->ne; ++k ) {
    pt = &mesh->tetra[k];
    if ( (!MG_EOK(pt)) || !pt->xt ) continue;

    pxt = &mesh->xtetra[pt->xt];
    for ( i=0; i<4; ++i ) {
      if ( !(pxt->ftag[i] & MG_PARBDY) ) continue;

      ia = pt->v[MMG5_idir[i][0]];
      ib = pt->v[MMG5_idir[i][1]];
      ic = pt->v[MMG5_idir[i][2]];
      if ( !MMG5_hashFace(mesh,&hash,ia,ib,ic,12*k+3*i)  ) {
        ier = 0;
        goto hash;
      }
    }
  }

  /** Step 2: Travel through the \a facesData array, get int the hash table the
   * index of the element to which belong the face and update the face
   * communicator */
  face2int_face_comm_index1 = grp->face2int_face_comm_index1;
  for ( k=0; k<nitem; ++k ) {

    /* Get the interface triangle vertices */
    ia = facesData[3*k  ];
    ib = facesData[3*k+1];
    ic = facesData[3*k+2];

    hashVal = MMG5_hashGetFace(&hash,ia,ib,ic);
    assert( hashVal );

    iel  =  hashVal/12;
    ifac = (hashVal%12)/3;

    assert( MG_EOK(&mesh->tetra[iel]) );

    pt = &mesh->tetra[iel];

    for ( iploc=0; iploc<3; ++iploc )
      if ( pt->v[MMG5_idir[ifac][iploc]] == ia ) break;
    assert ( iploc < 3 );

    /* Update the face communicator */
    face2int_face_comm_index1[k] = hashVal+iploc;
  }

hash:
  MMG5_DEL_MEM(mesh,hash.item);

facesData:
  PMMG_DEL_MEM(parmesh,facesData,int,"facesData");

  return ier;
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
  MMG3D_pPROctree *q;
  size_t     oldMemMax,available;
  int        it,ier,ier_end,ieresult,i,k, *facesData;
  mytime     ctim[TIMEMAX];
  int8_t     tim;
  char       stim[32];


  tminit(ctim,TIMEMAX);

  ier_end = PMMG_SUCCESS;

  /** Groups creation */
  if ( parmesh->info.imprim > PMMG_VERB_QUAL ) {
    tim = 0;
    chrono(ON,&(ctim[tim]));
  }

  ier = PMMG_split_grps( parmesh,PMMG_GRPSPL_MMG_TARGET,0 );

  MPI_CHECK ( MPI_Allreduce( &ier,&ieresult,1,MPI_INT,MPI_MIN,parmesh->comm ),
              PMMG_CLEAN_AND_RETURN(parmesh,PMMG_LOWFAILURE) );

  if ( parmesh->info.imprim > PMMG_VERB_STEPS ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"       group splitting                   %s\n",stim);
  }

  if ( !ieresult ) {
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_LOWFAILURE);
  }
  else if ( ieresult<0 ) {
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_STRONGFAILURE);
  }

  //DEBUGGING: PMMG_grplst_meshes_to_saveMesh(parmesh->listgrp, 1, parmesh->myrank, "Begin_libparmmg1_proc");

  /** Reset the boundary fields between the old mesh size and the new one (Mmg
   * uses this fields assiming they are setted to 0)/ */
  for ( i=0; i<parmesh->ngrp; ++i ) {
    mesh         = parmesh->listgrp[i].mesh;

    if ( !mesh ) continue;

    memset(&mesh->xtetra[mesh->xt+1],0,(mesh->xtmax-mesh->xt)*sizeof(MMG5_xTetra));
    memset(&mesh->xpoint[mesh->xp+1],0,(mesh->xpmax-mesh->xp)*sizeof(MMG5_xPoint));

    /* Uncomment to debug tag errors */
    /* if(!mesh->ntmax) mesh->ntmax = mesh->xtmax;*/
    /* if ( !MMG3D_analys(mesh) ) { PMMG_CLEAN_AND_RETURN(parmesh,PMMG_STRONGFAILURE); } */
  }

  /** Mesh adaptation */
  for ( it = 0; it < parmesh->niter; ++it ) {

    if ( parmesh->info.imprim > PMMG_VERB_STEPS ) {
      tim = 1;
      if ( it > 0 ) {
        chrono(OFF,&(ctim[tim]));
      }
      if ( parmesh->info.imprim > PMMG_VERB_ITWAVES ) {
        fprintf(stdout,"\n" );
      }

      printim(ctim[tim].gdif,stim);
      chrono(ON,&(ctim[tim]));
      fprintf(stdout,"\r       adaptation: iter %d   cumul. timer %s",it+1,stim);fflush(stdout);
    }

    /** Update old groups for metrics interpolation */
    PMMG_update_oldGrps( parmesh );
    PMMG_CALLOC( parmesh, q, parmesh->ngrp, MMG3D_pPROctree, "octrees",
                 PMMG_CLEAN_AND_RETURN(parmesh,PMMG_STRONGFAILURE) );

    tim = 4;
    if ( parmesh->info.imprim > PMMG_VERB_ITWAVES ) {
      chrono(RESET,&(ctim[tim]));
      chrono(ON,&(ctim[tim]));
    }

    for ( i=0; i<parmesh->ngrp; ++i ) {
      mesh         = parmesh->listgrp[i].mesh;
      met          = parmesh->listgrp[i].met;
      /* Reset the value of the fem mode */
      mesh->info.fem = parmesh->info.fem;

      if ( (!mesh->np) && (!mesh->ne) ) {
        /* Empty mesh */
        continue;
      }

      PMMG_TRANSFER_AVMEM_TO_PARMESH(parmesh,available,oldMemMax);

      /** Store the vertices of interface faces in the internal communicator */
      if ( !(ier = PMMG_store_faceVerticesInIntComm(parmesh,i,&facesData) ) ) {
        /* We are not able to remesh */
        fprintf(stderr,"\n  ## Interface faces storage problem."
                " Exit program.\n");
        break;
      }
      else {
        /* We can remesh */

        PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,parmesh->listgrp[i].mesh,
                                               available,oldMemMax);

        /* Mark reinitialisation in order to be able to remesh all the mesh */
        mesh->mark = 0;
        mesh->base = 0;
        for ( k=1 ; k<=mesh->nemax ; k++ ) {
          mesh->tetra[k].mark = mesh->mark;
          mesh->tetra[k].flag = mesh->base;
        }

        /** Call the remesher */
        /* Here we need to scale the mesh */
        if ( !MMG5_scaleMesh(mesh,met) ) { goto strong_failed; }

        if ( !mesh->adja ) {
          if ( !MMG3D_hashTetra(mesh,0) ) {
            fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
            goto strong_failed;
          }
        }

        /** Build tetras octree */
        if( !PMMG_initPROctree( mesh, &q[i], 1000 ) ) goto strong_failed;

#ifdef PATTERN
        ier = MMG5_mmg3d1_pattern( mesh, met );
#else
        ier = MMG5_mmg3d1_delone( mesh, met );
#endif

        if ( !ier ) {
          fprintf(stderr,"\n  ## MMG remeshing problem. Exit program.\n");
        }

        /** Pack the tetra */
        if ( mesh->adja )
          PMMG_DEL_MEM(mesh,mesh->adja,int,"adja table");

        if ( !MMG5_paktet(mesh) ) {
          fprintf(stderr,"\n  ## Tetra packing problem. Exit program.\n");
          goto strong_failed;
        }

        /** Update interface tetra indices in the face communicator */
        if ( ! PMMG_update_face2intInterfaceTetra(parmesh,i,facesData) ) {
          fprintf(stderr,"\n  ## Interface tetra updating problem. Exit program.\n");
          goto strong_failed;
        }
        if ( !MMG5_unscaleMesh(mesh,met) ) { goto strong_failed; }

        PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,parmesh->listgrp[i].mesh,
                                               available,oldMemMax);

        if ( !ier ) { break; }
      }
#warning Luca: Check modifications in Mmg
      /* Reset the mesh->gap field in case Mmg have modified it */
      mesh->gap = 0.2;
    }

    MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
    if ( parmesh->info.imprim > PMMG_VERB_ITWAVES ) {
      chrono(OFF,&(ctim[tim]));
      printim(ctim[tim].gdif,stim);
      fprintf(stdout,"\n       mmg                               %s\n",stim);
    }

    if ( !ieresult )
      goto failed_handling;

    /** Interpolate metrics */
    if ( parmesh->info.imprim > PMMG_VERB_ITWAVES ) {
      tim = 2;
      chrono(RESET,&(ctim[tim]));
      chrono(ON,&(ctim[tim]));
    }

    ier = PMMG_interpMetrics( parmesh );

    MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
    if ( parmesh->info.imprim > PMMG_VERB_ITWAVES ) {
      chrono(OFF,&(ctim[tim]));
      printim(ctim[tim].gdif,stim);
      fprintf(stdout,"       metric interpolation              %s\n",stim);
    }

    if ( !ieresult ) {
      if ( !parmesh->myrank )
        fprintf(stderr,"\n  ## Metrics interpolation problem. Try to save the mesh and exit program.\n");
      PMMG_CLEAN_AND_RETURN(parmesh,PMMG_STRONGFAILURE);
    }

    /* Destroy octrees */
    PMMG_freePROctrees( parmesh, q );
    PMMG_DEL_MEM( parmesh, q, MMG3D_pPROctree, "octrees" );

    /* Compute quality in the interpolated metrics */
    ier = PMMG_tetraQual( parmesh,0 );

    /** load Balancing at group scale and communicators reconstruction */
    tim = 3;
    if ( parmesh->info.imprim > PMMG_VERB_ITWAVES ) {
      chrono(RESET,&(ctim[tim]));
      chrono(ON,&(ctim[tim]));
    }

    ier = PMMG_loadBalancing(parmesh);

    MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
   if ( parmesh->info.imprim > PMMG_VERB_ITWAVES ) {
      chrono(OFF,&(ctim[tim]));
      printim(ctim[tim].gdif,stim);
      fprintf(stdout,"       load balancing                    %s\n",stim);
    }

    if ( !ieresult ) {
      if ( !parmesh->myrank )
        fprintf(stderr,"\n  ## Load balancing problem. Try to save the mesh and exit program.\n");
      goto failed_handling;
    } else if ( ieresult < 0 ) {
      if ( !parmesh->myrank )
        fprintf(stderr,"\n  ## Load balancing problem. Exit program.\n");
      PMMG_CLEAN_AND_RETURN(parmesh,PMMG_STRONGFAILURE);
    }
  }

  if ( parmesh->info.imprim > PMMG_VERB_STEPS ) {
    printf("\n");
  }

  ier = PMMG_qualhisto( parmesh, PMMG_OUTQUA );

  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    ier_end = PMMG_LOWFAILURE;
  }

  if ( parmesh->info.imprim > PMMG_VERB_STEPS ) {
    tim = 4;
    chrono(ON,&(ctim[tim]));
  }

  ier = PMMG_packParMesh(parmesh);
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( parmesh->info.imprim > PMMG_VERB_STEPS ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"\n       mesh packing                      %s\n",stim);
  }

  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Parallel mesh packing problem. Exit program.\n");
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_STRONGFAILURE);
  }

  PMMG_listgrp_free( parmesh, &parmesh->old_listgrp, parmesh->nold_grp);

  if ( parmesh->info.imprim > PMMG_VERB_STEPS ) {
    tim = 5;
    chrono(ON,&(ctim[tim]));
  }

  ier = PMMG_merge_grps(parmesh);
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );

  if ( parmesh->info.imprim > PMMG_VERB_STEPS ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"       group merging                     %s\n",stim);
  }

  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Groups merging problem. Exit program.\n");
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_STRONGFAILURE);
  }

  /* Give memory to Mmg for the edge length computation */
  PMMG_TRANSFER_AVMEM_TO_MESHES(parmesh);

  if ( parmesh->info.imprim0 > PMMG_VERB_ITWAVES && !parmesh->info.iso ) {
    PMMG_prilen(parmesh,0);
  }

  PMMG_CLEAN_AND_RETURN(parmesh,ier_end);

  /** mmg3d1_delone failure */
strong_failed:
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  PMMG_CLEAN_AND_RETURN(parmesh,PMMG_STRONGFAILURE);

failed_handling:
  if ( parmesh->info.imprim > PMMG_VERB_STEPS ) {
    tim = 4;
    chrono(ON,&(ctim[tim]));
  }
  if ( !PMMG_packParMesh(parmesh) ) {
    fprintf(stderr,"\n  ## Parmesh packing problem. Exit program.\n");
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_STRONGFAILURE);
  }
  if ( parmesh->info.imprim > PMMG_VERB_STEPS ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"\n       mesh packing                      %s\n",stim);
  }

  if ( parmesh->info.imprim > PMMG_VERB_STEPS ) {
    chrono(ON,&(ctim[5]));
  }
  if ( !PMMG_merge_grps(parmesh) ) {
    fprintf(stderr,"\n  ## Groups merging problem. Exit program.\n");
    PMMG_CLEAN_AND_RETURN(parmesh,PMMG_STRONGFAILURE);
  }
  if ( parmesh->info.imprim > PMMG_VERB_STEPS ) {
    chrono(OFF,&(ctim[5]));
    printim(ctim[5].gdif,stim);
    fprintf(stdout,"       group merging                     %s\n",stim);
  }

  PMMG_CLEAN_AND_RETURN(parmesh,PMMG_LOWFAILURE);
}

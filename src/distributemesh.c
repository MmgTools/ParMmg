/**
 * \file distribute.c
 * \brief Distribute the mesh on the processors.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */


#include "libparmmg.h"
/**
 * \param parmesh pointer toward the mesh structure.
 * \param part pointer toward an array of int containing the partitions (proc master).
 * \return 0 if the file is not found, -1 if we detect mismatch parameters,
 * 1 otherwise.
 *
 * Proc 0 send at all the processors a part of the mesh.
 *
 */
int PMMG_distributeMesh(PMMG_pParMesh parmesh,int *part) {
  PMMG_pGrp    grp;
  MMG5_pMesh   mesh;
#warning add solution/metric communication
  MMG5_pSol    sol;
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  MMG5_pPoint  ppt;
  MMG5_pxPoint pxp;
  double       *pointcomm,*xPtcomm;
  int          i,j,k,iadr,*nbTet,*nbxTet,*nbPt,*nbxPt;
  int          *targetpoint,*tetracomm,*xTetcomm,*xpcomm,*xtcomm;
  int          iproc,icurc,itet,ixTet,nxp,nxt,ixPt;
  MPI_Status   status;

  _MMG5_SAFE_CALLOC(nbTet,parmesh->nprocs,int);
  _MMG5_SAFE_CALLOC(nbPt,parmesh->nprocs,int);
  _MMG5_SAFE_CALLOC(nbxTet,parmesh->nprocs,int);
  _MMG5_SAFE_CALLOC(nbxPt,parmesh->nprocs,int);

  if ( !parmesh->myrank ) {

    grp = parmesh->listgrp;
    mesh = grp[0].mesh;

    /*mark tetra with target proc*/
    /*count number of elt per proc*/
    for(k=1 ; k<=mesh->ne ; k++) {
      pt = &mesh->tetra[k];
      pt->mark = part[k-1];
      nbTet[part[k-1]]++;
      if ( pt->xt ) nbxTet[part[k-1]]++;
    }
  } else {
    /* In a first stage, we have only 1 group per processor */
    parmesh->ngrp = 1;
    _MMG5_SAFE_CALLOC(parmesh->listgrp,parmesh->ngrp,PMMG_Grp);

    grp = &parmesh->listgrp[0];

    grp->mesh = NULL;
    grp->sol  = NULL;

    MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&grp->mesh,
                    MMG5_ARG_ppMet,&grp->sol,
                    MMG5_ARG_end);
    mesh = grp->mesh;

  }

  /*send at each proc its number of elt*/
  MPI_Bcast(&nbTet[0],parmesh->nprocs,MPI_INT,0,parmesh->comm);
  MPI_Bcast(&nbxTet[0],parmesh->nprocs,MPI_INT,0,parmesh->comm);

  printf("je suis %d, recu ok %d\n",parmesh->myrank,nbTet[parmesh->myrank]);

  if ( !parmesh->myrank ) {
    for ( iproc=1 ; iproc<parmesh->nprocs ; iproc++ ) {
      _MMG5_SAFE_CALLOC(tetracomm,5*nbTet[iproc],int);
      _MMG5_SAFE_CALLOC(xtcomm,nbxTet[iproc],int);
      _MMG5_SAFE_CALLOC(xTetcomm,10*nbxTet[iproc],int);
      _MMG5_SAFE_CALLOC(pointcomm,12*nbTet[iproc],double);
      _MMG5_SAFE_CALLOC(xpcomm,4*nbxTet[iproc],int);
      _MMG5_SAFE_CALLOC(xPtcomm,24*nbxTet[iproc],double);
      _MMG5_SAFE_CALLOC(targetpoint,mesh->np+1,int);

      icurc = 0;
      itet  = 0;
      ixTet = 0;
      ixPt  = 0;
      nxp   = 0;
      nxt   = 0;

      for ( k=1; k<=mesh->ne; k++ ) {
        pt = &mesh->tetra[k];

        if ( pt->mark != iproc ) continue;

        for ( j=0; j<4; j++ ) {
          ppt = &mesh->point[pt->v[j]];

          if ( !targetpoint[pt->v[j]] ) {
            targetpoint[pt->v[j]] = icurc/3;
            pointcomm[icurc++] = ppt->c[0];
            pointcomm[icurc++] = ppt->c[1];
            pointcomm[icurc++] = ppt->c[2];

            if ( icurc == 12*nbTet[iproc] ) {
              printf("memory alloc problem in distribute mesh\n");
              exit(EXIT_FAILURE);
            }
            nbPt[iproc]++;

            /* Store xpoint if needed */
            if ( ppt->xp ) {
              nbxPt[iproc]++;
              pxp = &mesh->xpoint[ppt->xp];

              xpcomm[nxp++] = targetpoint[pt->v[j]];

              for ( i=0; i<3; ++i )
                xPtcomm[ixPt++] = pxp->n1[i];

              for ( i=0; i<3; ++i )
                xPtcomm[ixPt++] = pxp->n2[i];
            }
          }
          tetracomm[itet++] = targetpoint[pt->v[j]];
        }
        tetracomm[itet++] = pt->ref;

        /* Store xtetra if needed */
        if ( !pt->xt ) continue;

        pxt = &mesh->xtetra[pt->xt];
        xtcomm[nxt++] = (itet-1)/4;

        for ( j=0; j<4; ++j )
          xTetcomm[ixTet++] = pxt->ref[j];

        for ( j=0; j<6; ++j )
          xTetcomm[ixTet++] = pxt->edg[j];

      }

      //assert(itet==nbTet[iproc]);
      /* Send number of points and xpoints */
      MPI_Send(&nbPt[iproc],1,MPI_INT,iproc,0,parmesh->comm);
      MPI_Send(&nbxPt[iproc],1,MPI_INT,iproc,0,parmesh->comm);

      /* Send mesh data */
      MPI_Send(&tetracomm[0],5*nbTet[iproc],MPI_INT,iproc,0,parmesh->comm);
      MPI_Send(&xtcomm[0],nbxTet[iproc],MPI_INT,iproc,0,parmesh->comm);
      MPI_Send(&xTetcomm[0],10*nbxTet[iproc],MPI_INT,iproc,0,parmesh->comm);
      MPI_Send(&pointcomm[0],3*nbPt[iproc],MPI_DOUBLE,iproc,0,parmesh->comm);
      MPI_Send(&xpcomm[0],nbxPt[iproc],MPI_INT,iproc,0,parmesh->comm);
      MPI_Send(&xPtcomm[0],6*nbxPt[iproc],MPI_DOUBLE,iproc,0,parmesh->comm);

      _MMG5_SAFE_FREE(targetpoint);
      _MMG5_SAFE_FREE(tetracomm);
      _MMG5_SAFE_FREE(pointcomm);
      _MMG5_SAFE_FREE(xtcomm);
      _MMG5_SAFE_FREE(xTetcomm);
      _MMG5_SAFE_FREE(xpcomm);
      _MMG5_SAFE_FREE(xPtcomm);
    }

    for ( k=1; k<=mesh->ne; k++) {
      pt = &mesh->tetra[k];
      if ( pt->mark!=0 ) _MMG3D_delElt(mesh,k);
    }

    /* create tetra adjacency */
    if ( !MMG3D_hashTetra(mesh,1) ) {
      fprintf(stderr,"  ## PMMG Hashing problem (1). Exit program.\n");
      return(0);
    }

    /* Pack the mesh */
    if ( !_MMG3D_packMesh(mesh,NULL,NULL ) ) {
      fprintf(stderr,"  ## PMMG Packing problem (1). Exit program.\n");
      return(0);
    }
    MMG3D_saveMesh(mesh,"proc0.mesh");

  }
  else {

    /* Recieve number of points and xpoints */
    MPI_Recv(&mesh->np,1,MPI_INT,0,0,parmesh->comm,MPI_STATUS_IGNORE);
    MPI_Recv(&mesh->xp,1,MPI_INT,0,0,parmesh->comm,MPI_STATUS_IGNORE);

    /* Recieve mesh data */
    _MMG5_SAFE_CALLOC(tetracomm,5*nbTet[parmesh->myrank],int);
    _MMG5_SAFE_CALLOC(xtcomm,nbxTet[parmesh->myrank],int);
    _MMG5_SAFE_CALLOC(xTetcomm,10*nbxTet[parmesh->myrank],int);
    _MMG5_SAFE_CALLOC(pointcomm,3*mesh->np,double);
    _MMG5_SAFE_CALLOC(xpcomm,mesh->xp,int);
    _MMG5_SAFE_CALLOC(xPtcomm,6*mesh->xp,double);

    MPI_Recv(&tetracomm[0],5*nbTet[parmesh->myrank],MPI_INT,0,0,parmesh->comm,MPI_STATUS_IGNORE);
    MPI_Recv(&xtcomm[0],nbxTet[parmesh->myrank],MPI_INT,0,0,parmesh->comm,MPI_STATUS_IGNORE);
    MPI_Recv(&xTetcomm[0],10*nbxTet[parmesh->myrank],MPI_INT,0,0,parmesh->comm,MPI_STATUS_IGNORE);
    MPI_Recv(&pointcomm[0],3*mesh->np,MPI_DOUBLE,0,0,parmesh->comm,MPI_STATUS_IGNORE);
    MPI_Recv(&xpcomm[0],mesh->xp,MPI_INT,0,0,parmesh->comm,MPI_STATUS_IGNORE);
    MPI_Recv(&xPtcomm[0],6*mesh->xp,MPI_DOUBLE,0,0,parmesh->comm,MPI_STATUS_IGNORE);

    /*MMG5_mesh creation*/
    if ( MMG3D_Set_meshSize(mesh,mesh->np,nbTet[parmesh->myrank],0,0,0,0) != 1 )
      return(0);

    /* Points */
    icurc = 0;
    for ( k=1; k<=mesh->np; k++ ) {
      ppt = &mesh->point[k];
      for(i=0 ; i<3 ; i++) {
        ppt->c[i] = pointcomm[icurc++];
      }
    }

    /* xPoints */
    _MMG5_SAFE_CALLOC(mesh->xpoint,mesh->xp+1,MMG5_xPoint);
    icurc = 0;
    for ( k=1; k<=mesh->xp; ++k ) {
      pxp = &mesh->xpoint[k];

      for ( i=0; i<3; ++i )
        pxp->n1[i] = xPtcomm[icurc++];

      for ( i=0; i<3; ++i )
        pxp->n2[i] = xPtcomm[icurc++];

      mesh->point[xpcomm[k-1]+1].xp = k;
    }

    /* Tetra */
    mesh->xt = nbxTet[parmesh->myrank];
    _MMG5_SAFE_CALLOC(mesh->xtetra,mesh->xt+1,MMG5_xTetra);
    icurc = 0;
    for(k=1 ; k<=nbTet[parmesh->myrank] ; k++) {
      pt = &mesh->tetra[k];
      for(i=0 ; i<4 ; i++) {
        pt->v[i] = tetracomm[icurc++]+1;
      }
      pt->ref = tetracomm[icurc++];
    }

    /* xTetra */
    icurc = 0;
    for ( k=1; k<=mesh->xt; ++k ) {
      pxt = &mesh->xtetra[k];

      for ( i=0; i<4; ++i )
        pxt->ref[i] = xTetcomm[icurc++];

      for ( i=0; i<6; ++i )
        pxt->edg[i] = xTetcomm[icurc++];

      mesh->tetra[xtcomm[k-1]].xt = k;
    }

    /* create tetra adjacency */
    if ( !MMG3D_hashTetra(mesh,1) ) {
      fprintf(stderr,"  ## PMMG Hashing problem (1). Exit program.\n");
      return(0);
    }

    /* Pack the mesh */
    if ( !_MMG3D_packMesh(mesh,NULL,NULL ) ) {
      fprintf(stderr,"  ## PMMG Packing problem (1). Exit program.\n");
      return(0);
    }
    MMG3D_saveMesh(mesh,"proc1.mesh");

    _MMG5_SAFE_FREE(tetracomm);
    _MMG5_SAFE_FREE(pointcomm);
    _MMG5_SAFE_FREE(xtcomm);
    _MMG5_SAFE_FREE(xTetcomm);
    _MMG5_SAFE_FREE(xpcomm);
    _MMG5_SAFE_FREE(xPtcomm);
  }

  _MMG5_SAFE_FREE(nbTet);
  _MMG5_SAFE_FREE(nbPt);
  _MMG5_SAFE_FREE(nbxTet);
  _MMG5_SAFE_FREE(nbxPt);

  return(1);
}

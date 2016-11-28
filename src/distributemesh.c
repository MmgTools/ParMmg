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
  PMMG_pGrp   grp;
  MMG5_pMesh  mesh;
  MMG5_pTetra pt;
  MMG5_pPoint ppt;
  double      *pointcomm;
  int         i,j,k,iadr,*nbelt,*nbpt,*targetpoint,*tetracomm;
  int         iproc,icurc,itet;
  MPI_Status  status;

  _MMG5_SAFE_CALLOC(nbelt,parmesh->nprocs,int);
  _MMG5_SAFE_CALLOC(nbpt,parmesh->nprocs,int);

  if ( !parmesh->myrank ) {

    grp = parmesh->listgrp;
    mesh = grp[0].mesh;

    /*mark tetra with target proc*/
    /*count number of elt per proc*/
    for(k=1 ; k<=mesh->ne ; k++) {
      pt = &mesh->tetra[k];
      pt->mark = part[k-1];
      nbelt[part[k-1]]++;
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
  MPI_Bcast(&nbelt[0],parmesh->nprocs,MPI_INT,0,parmesh->comm);

  printf("je suis %d, recu ok %d\n",parmesh->myrank,nbelt[parmesh->myrank]);

  if(!parmesh->myrank) {
    for(iproc=1 ; iproc<parmesh->nprocs ; iproc++) {
      _MMG5_SAFE_CALLOC(tetracomm,5*nbelt[iproc],int);
      _MMG5_SAFE_CALLOC(pointcomm,12*nbelt[iproc],double);
      _MMG5_SAFE_CALLOC(targetpoint,mesh->np+1,int);
      icurc = 0;
      itet = 0;
      for(k=1 ; k<=mesh->ne ; k++) {
        pt = &mesh->tetra[k];
        if(pt->mark != iproc) continue;
        for(j=0 ; j<4 ; j++) {
          ppt = &mesh->point[pt->v[j]];
          if(!targetpoint[pt->v[j]]) {
            targetpoint[pt->v[j]] = icurc/3;
            pointcomm[icurc++] = ppt->c[0];
            pointcomm[icurc++] = ppt->c[1];
            pointcomm[icurc++] = ppt->c[2];
            if(icurc == 12*nbelt[iproc]) printf("memory alloc problem in distribute mesh\n");
            nbpt[iproc]++;
          }
          tetracomm[itet++] = targetpoint[pt->v[j]];
        }
        tetracomm[itet++] = pt->ref;
      }
      //assert(itet==nbelt[iproc]);
      MPI_Send(&tetracomm[0],5*nbelt[iproc],MPI_INT,iproc,0,parmesh->comm);
      MPI_Send(&pointcomm[0],12*nbelt[iproc],MPI_DOUBLE,iproc,0,parmesh->comm);
      _MMG5_SAFE_FREE(targetpoint);
      _MMG5_SAFE_FREE(tetracomm);
      _MMG5_SAFE_FREE(pointcomm);
    }
    for(k=1 ; k<=mesh->ne ; k++) {
      pt = &mesh->tetra[k];
      if(pt->mark!=0) _MMG3D_delElt(mesh,k);
    }
    MMG3D_saveMesh(mesh,"proc0.mesh");
  } else {
    _MMG5_SAFE_CALLOC(tetracomm,5*nbelt[parmesh->myrank],int);
    _MMG5_SAFE_CALLOC(pointcomm,12*nbelt[iproc],double);
    MPI_Recv(&tetracomm[0],5*nbelt[parmesh->myrank],MPI_INT,0,0,parmesh->comm,MPI_STATUS_IGNORE);
    MPI_Recv(&pointcomm[0],12*nbelt[parmesh->myrank],MPI_DOUBLE,0,0,parmesh->comm,MPI_STATUS_IGNORE);

    /*MMG5_mesh creation*/
    if ( MMG3D_Set_meshSize(mesh,4*nbelt[parmesh->myrank],nbelt[parmesh->myrank],0,0,0,0) != 1 )  exit(EXIT_FAILURE);

    icurc = 0;
    for(k=1 ; k<=4*nbelt[parmesh->myrank] ; k++) {
      ppt = &mesh->point[k];
      for(i=0 ; i<3 ; i++) {
        ppt->c[i] = pointcomm[icurc++];
      }
    }
    icurc = 0;
    for(k=1 ; k<=nbelt[parmesh->myrank] ; k++) {
      pt = &mesh->tetra[k];
      for(i=0 ; i<4 ; i++) {
        pt->v[i] = tetracomm[icurc++]+1;
      }
      pt->ref = tetracomm[icurc++];
    }
    MMG3D_saveMesh(mesh,"proc1.mesh");

    _MMG5_SAFE_FREE(tetracomm);
    _MMG5_SAFE_FREE(pointcomm);
  }

  _MMG5_SAFE_FREE(nbelt);
  _MMG5_SAFE_FREE(nbpt);


  return(1);
}

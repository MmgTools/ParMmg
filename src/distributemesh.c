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
  int         j,k,iadr,*nbelt;
  
 
  _MMG5_SAFE_CALLOC(nbelt,parmesh->nprocs,int);

  if(!parmesh->myrank) {
    grp = parmesh->listgrp;
    mesh = grp[0].mesh;
    /*mark tetra with target proc*/
    /*count number of elt per proc*/

    for(k=1 ; k<=mesh->ne ; k++) {
      pt = &mesh->tetra[k];
      pt->mark = part[k-1];
      nbelt[part[k-1]]++;
    }
    for(k=0 ; k<parmesh->nprocs ; k++) {
      printf("in proc %d : %d tet\n",k,nbelt[k]);
    }
  } else {
    parmesh->ngrp = 1;
    _MMG5_SAFE_CALLOC(parmesh->listgrp,parmesh->ngrp,PMMG_Grp);

    grp = &parmesh->listgrp[0];
    
    grp->mesh = NULL;
    grp->sol  = NULL;
    
    MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,&grp->mesh,
		    MMG5_ARG_ppMet,&grp->sol,
		    MMG5_ARG_end);
  }

  /*send at each proc its number of elt*/
  MPI_Bcast(&nbelt[0],parmesh->nprocs,MPI_INT,0,parmesh->comm);

  printf("je suis %d, recu ok %d\n",parmesh->myrank,nbelt[parmesh->myrank]);

  _MMG5_SAFE_FREE(nbelt);
  return(1);
}

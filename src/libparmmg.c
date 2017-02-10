#include "parmmg.h"


int PMMG_parmmglib(PMMG_pParMesh parmesh) {
  MMG5_pMesh       mesh;
  MMG5_pSol        sol;
  int              it,i,ier,niter;

#warning niter must be a param setted by the user
  niter = 1;


  /** Mesh adaptation */
  for ( it=0; it<niter; ++it ) {
    for ( i=0; i<parmesh->ngrp; ++i ) {
      mesh = parmesh->listgrp[i].mesh;
      sol  = parmesh->listgrp[i].sol;

      if ( !MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_nosurf,1 ) )
        return PMMG_STRONGFAILURE;

#warning for debugging purposes
      if ( !MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_noinsert,1 ) )
        return PMMG_STRONGFAILURE;

      if ( !MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_noswap,1 ) )
        return PMMG_STRONGFAILURE;

      if ( !MMG3D_Set_iparameter(mesh,sol,MMG3D_IPARAM_nomove,1 ) )
        return PMMG_STRONGFAILURE;

      ier = 0;//MMG3D_mmg3dlib(mesh,sol);

      if ( ier == MMG5_STRONGFAILURE ) return PMMG_STRONGFAILURE;

#warning Do we need to update the communicators? Does Mmg renum the boundary nodes with -nosurf option?

      /** load Balancing and communicators reconstruction */

    }
  }

  /** Merge all the meshes on the proc 0 */
  if ( !PMMG_mergeParMesh(parmesh) )  return(PMMG_STRONGFAILURE);

  _MMG3D_packMesh(parmesh->listgrp[0].mesh,parmesh->listgrp[0].sol,NULL);

  if ( parmesh->ddebug &&  !parmesh->myrank ) {
    MMG3D_saveMesh(parmesh->listgrp[0].mesh,"End_libparmmg.mesh");
    MMG3D_saveSol(parmesh->listgrp[0].mesh,parmesh->listgrp[0].sol,"End_libparmmg.sol");
  }
  return(PMMG_SUCCESS);
}

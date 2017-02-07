#include "libparmmg.h"


int PMMG_parmmglib(PMMG_pParMesh parmesh) {
  MMG5_pMesh       mesh;
  MMG5_pSol        sol;
  int              *part,it,i,ier,niter;

#warning niter must be a param setted by the user
  niter = 1;

  _MMG5_SAFE_CALLOC(part,(parmesh->listgrp[0].mesh)->ne,int);

  /** Call metis for partionning*/
  if(!PMMG_metispartitioning(parmesh,part)) return(PMMG_STRONGFAILURE);

  /** Mesh analysis: compute ridges, singularities, normals... and store the
      triangles into the xTetra structure */
  if ( !_MMG3D_analys(parmesh->listgrp[0].mesh) ) return(PMMG_STRONGFAILURE);

  /** Send mesh partionning to other proc*/
  if ( !PMMG_distributeMesh(parmesh,part) ) return(PMMG_STRONGFAILURE);
  _MMG5_SAFE_FREE(part);

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
  if ( !PMMG_mergeMesh(parmesh) )  return(PMMG_STRONGFAILURE);

  return(PMMG_SUCCESS);
}

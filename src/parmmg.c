#include "parmmg.h"

/**
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \return \ref PMMG_SUCCESS if success.
 * \return \ref PMMG_LOWFAILURE if failed but a conform mesh is saved.
 * \return \ref PMMG_STRONGFAILURE if failed and we can't save the mesh.
 *
 * Main program for PARMMG executable: perform parallel mesh adaptation.
 *
 */
int main(int argc,char *argv[]) {
  PMMG_pParMesh    parmesh=NULL;
  PMMG_pGrp        grp;
  MMG5_pMesh       mesh;
  MMG5_pSol        sol;
  int              i,ier;

  /** Init MPI */
  MPI_Init(&argc, &argv);

  /** Assign default values */
  if ( !PMMG_Init_parMesh(PMMG_ARG_start,
                          PMMG_ARG_ppParMesh,&parmesh,
                          PMMG_ARG_end) ) return PMMG_STRONGFAILURE;

  if ( !parmesh->myrank ) {
    fprintf(stdout,"  -- PARMMG3d, Release %s (%s) \n",PMMG_VER,PMMG_REL);
    fprintf(stdout,"     %s\n",PMMG_CPY);
    fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);

    /** Read sequential mesh */
#warning : for the moment, we only read a mesh named m.mesh

    if ( PMMG_loadMesh(parmesh,"m.mesh") < 1 ) return(PMMG_STRONGFAILURE);
    if ( PMMG_loadSol(parmesh,"m.sol") < 0 ) return(PMMG_STRONGFAILURE);
  }

  /** Send mesh partionning to other proc*/
  if ( !PMMG_distributeMesh(parmesh) ) return(PMMG_STRONGFAILURE);

  ier = _PMMG_parmmglib1(parmesh);

  if ( ier!= PMMG_STRONGFAILURE ) {

    /** Merge all the meshes on the proc 0 */
    if ( !PMMG_mergeParMesh(parmesh,0) )  return PMMG_STRONGFAILURE;

    if ( !parmesh->myrank ) {
      if ( _MMG3D_bdryBuild(parmesh->listgrp[0].mesh) < 0 )
        return(PMMG_STRONGFAILURE);

      /* Write mesh */
#warning : for the moment, we only write a mesh named out.mesh
      if ( !PMMG_saveMesh(parmesh,"out.mesh") ) return(PMMG_STRONGFAILURE);
      if ( !PMMG_saveSol(parmesh,"out.sol") ) return(PMMG_STRONGFAILURE);
    }
  }

  /* Free structures */
#warning todo: create an API function to free a whole parmesh?
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;
  sol  = grp->sol;

  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&sol,
                 MMG5_ARG_end);

  _MMG5_SAFE_FREE(parmesh);

  /** Finalize MPI */
  MPI_Finalize();

  return(ier);
}

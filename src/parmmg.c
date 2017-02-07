#include "libparmmg.h"

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
  }

  /** Read sequential mesh */
#warning : for the moment, we only read a mesh named m.mesh
#warning Algiane: with lot of procs mpi process may fail to read the same file at the same time so maybe we will need to read the mesh over one unique proc and to broadcast it over the others...
  if ( !parmesh->myrank ) {
    if ( !PMMG_loadMesh(parmesh,"m.mesh") ) return(PMMG_STRONGFAILURE);
    if ( PMMG_loadSol(parmesh,"m.sol") < 0 ) return(PMMG_STRONGFAILURE);
  }

  ier = PMMG_parmmglib(parmesh);

  if ( ier!= PMMG_STRONGFAILURE && !parmesh->myrank ) {
    /*write mesh*/
#warning : for the moment, we only write a mesh named out.mesh
    if ( !PMMG_saveMesh(parmesh,"out.mesh") ) return(PMMG_STRONGFAILURE);
    if ( !PMMG_saveSol(parmesh,"out.sol") ) return(PMMG_STRONGFAILURE);
  }

  /*free structures*/
#warning create an API function to free a whole parmesh
  grp  = &parmesh->listgrp[0];
  mesh = grp->mesh;
  sol  = grp->sol;

  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mesh,MMG5_ARG_ppMet,&sol,
                 MMG5_ARG_end);
  _MMG5_SAFE_FREE(grp->node2int_edge_comm_index1);
  _MMG5_SAFE_FREE(grp->node2int_edge_comm_index2);
  _MMG5_SAFE_FREE(parmesh->int_node_comm->intvalues);

  for ( i=0; i<parmesh->next_node_comm; ++i ) {
    _MMG5_SAFE_FREE(parmesh->ext_node_comm->int_comm_index);
  }
  _MMG5_SAFE_FREE(parmesh->ext_node_comm);

  /*Finalize MPI*/
  MPI_Finalize();

  return(ier);
}

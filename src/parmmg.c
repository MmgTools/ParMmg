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
  PMMG_pParMesh    parmesh;
  int              *part;

  if ( parmesh )  _MMG5_SAFE_FREE(parmesh);
  _MMG5_SAFE_CALLOC(parmesh,1,PMMG_ParMesh);

  /* Init MPI*/
  parmesh->comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(parmesh->comm, &parmesh->nprocs);
  MPI_Comm_rank(parmesh->comm, &parmesh->myrank);

#warning CECILE : do we need a variadic Init for Parmmg ?
  //???PMMG_Init(parmesh);

  if ( !parmesh->myrank ) {
    fprintf(stdout,"  -- PARMMG3d, Release %s (%s) \n",PMMG_VER,PMMG_REL);
    fprintf(stdout,"     %s\n",PMMG_CPY);
    fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);
  }

  /* Read sequential mesh*/
#warning : for the moment, we only read a mesh named m.mesh
  if(!PMMG_loadMesh(parmesh,"m.mesh")) return(PMMG_STRONGFAILURE);

  _MMG5_SAFE_CALLOC(part,(parmesh->listgrp[0].mesh)->ne,int);
  if ( !parmesh->myrank ) {
    /* Call metis for partionning*/
    if(!PMMG_metispartitioning(parmesh,part)) return(PMMG_STRONGFAILURE);
  }

  /* Send at each proc the graph */
  MPI_Bcast(&part[0],(parmesh->listgrp[0].mesh)->ne,MPI_INT,0,parmesh->comm);


  /* Mesh analysis: compute ridges, singularities, normals... and store the
     triangles into the xTetra structure */
  if ( !_MMG3D_analys(parmesh->listgrp[0].mesh) ) return(PMMG_STRONGFAILURE);

  /*send mesh partionning to other proc*/
  if ( !PMMG_distributeMesh(parmesh,part) ) return(PMMG_STRONGFAILURE);
  _MMG5_SAFE_FREE(part);


  /*perform mesh adaptation*/



  if(!parmesh->myrank) {
    /*receive mesh*/

    /*write mesh*/
#warning : for the moment, we only write a mesh named out.mesh
    PMMG_saveMesh(parmesh,"out.mesh");
  } else {
    /*send mesh*/
  }

  /*free structures*/

  /*Finalize MPI*/
  MPI_Finalize();

  return(PMMG_SUCCESS);
}

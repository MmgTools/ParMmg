/**
 * \file API_functions.c
 * \brief C API functions definitions for PARMMG library.
 * \copyright GNU Lesser General Public License.
 *
 * C API for PARMMG library.
 *
 */

#include "libparmmg.h"

/**
 * \param .
 *
 * Allocation and initialization of parmmg structures.
 *
 */
int PMMG_InitParMesh() {
  return(1);
}
/**
 * \param .
 *
 *
 */
int PMMG_Set_grpSize(PMMG_pGrp grp, int np, int ne, int nprism, int nt,
                     int nquad, int na, int typEntity, int typSol, int typMet){
  // Set mesh size
  if(!MMG3D_Set_meshSize(grp->mesh, np, ne, nprism, nt, nquad, na)){
    fprintf(stderr,"  ## Error in setting group mesh size.\n");
    return(0);
  }
  // Set physical solution size
  if(!MMG3D_Set_solSize(grp->mesh, grp->sol, typEntity, np, typSol)){
    fprintf(stderr,"  ## Error in setting group solution size.\n");
    return(0);
  }
  // Set metrics size
  if(!MMG3D_Set_solSize(grp->mesh, grp->met, typEntity, np, typMet)){
    fprintf(stderr,"  ## Error in setting group metrics size.\n");
    return(0);
  }
  return(1);
}

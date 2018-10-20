/**
 * \file scalem_pmmg.c
 * \brief Mesh scaling/unscaling.
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */
#include "parmmg.h"
#include "mpitypes_pmmg.h"
#include "metis_pmmg.h"


/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the metric or solution structure.
 * \return 1 if success, 0 if fail (computed bounding box too small
 * or one af the anisotropic input metric is not valid).
 *
 * Scale the mesh and the size informations between 0 and 1.
 * Compute a default value for the hmin/hmax parameters if needed.
 *
 */

int PMMG_scaleMesh(MMG5_pMesh mesh,MMG5_pSol met) {
  double     dd;
  
  if ( !MMG5_scaleMesh(mesh,met) ) return 0;

  if ( mesh->info.hsiz > 0 || mesh->info.optim ) {
    dd = 1.0 / mesh->info.delta;
    mesh->info.hmin *= dd;
    mesh->info.hmax *= dd;
  }
  
  return 1;
}

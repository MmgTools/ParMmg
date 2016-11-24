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


  return(1);
}

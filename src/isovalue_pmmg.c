/* =============================================================================
**  This file is part of the parmmg software package for parallel tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux, 2017-
**
**  parmmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  parmmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with parmmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the parmmg distribution only if you accept them.
** =============================================================================
*/

/**
 * \file isovalue_pmmg.c
 * \brief Create implicit surface in distribuited mesh.
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (InriaSoft)
 * \author Laetitia Mottet (UBordeaux)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Main library functions (parallel remeshing starting from centralized or
 * distributed data.
 *
 */

#include "parmmg.h"
#include "mmgexterns_private.h"

/**
 * \param parmesh pointer toward a parmesh structure
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set values.
 * \param met pointer toward a metric (non-mandatory).
 * \return 1 if success, 0 otherwise.
 *
 * Proceed to discretization of the implicit function carried by sol into mesh,
 * once values of sol have been snapped/checked
 *
 */
int PMMG_cuttet_ls(PMMG_pParMesh parmesh,MMG5_pMesh mesh, MMG5_pSol sol,MMG5_pSol met){
  fprintf(stdout,"\n  ## TODO:: PMMG_cuttet_ls.\n");
  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param mesh pointer toward the mesh.
 *
 * Reset mesh->info.isoref vertex and tetra references to 0.
 *
 * \warning to improve: for now, entities linked to the old ls (corners,required
 * points, normals/tangents, triangles and edges) are deleted in loadMesh. It
 * would be better to analyze wich entities must be kept and which one must be
 * deleted depending on the split/nosplit infos.
 */
int PMMG_resetRef_ls(PMMG_pParMesh parmesh,MMG5_pMesh mesh) {
  fprintf(stdout,"\n  ## TODO:: PMMG_resetRef_ls.\n");
  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set values.
 * \return 1.
 *
 * Set references to tets according to the sign of the level set function.
 *
 */
int PMMG_setref_ls(PMMG_pParMesh parmesh,MMG5_pMesh mesh, MMG5_pSol sol) {
  fprintf(stdout,"\n  ## TODO:: PMMG_setref_ls.\n");
  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set function.
 * \return 1 if success, 0 if fail.
 *
 * Snap values of the level set function very close to 0 to exactly 0,
 * and prevent nonmanifold patterns from being generated.
 *
 */
int PMMG_snpval_ls(PMMG_pParMesh parmesh,MMG5_pMesh mesh,MMG5_pSol sol) {
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n      ## TODO:: PMMG_snpval_ls.\n");
  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set.
 * \param met pointer toward  a metric (optionnal).
 * \return 0 if fail, 1 otherwise.
 *
 * Create implicit surface in mesh.
 *
 */
int PMMG_ls(PMMG_pParMesh parmesh,MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pSol met) {
  char str[16]="";

  /* Set function pointers */
  // if ( mesh->info.isosurf ) {
  //   fprintf(stdout,"  ** TODO - ls on surface\n");
  // }
  // else {
  //   PMMG_snpval   = PMMG_snpval_ls;
  //   PMMG_resetRef = PMMG_resetRef_ls;
  //   PMMG_cuttet   = PMMG_cuttet_ls;
  //   PMMG_setref   = PMMG_setref_ls;
  // }

  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"  ** ISOSURFACE EXTRACTION %s\n",str);

  if ( mesh->nprism || mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: Isosurface extraction not available with"
            " hybrid meshes. Exit program.\n");
    return 0;
  }

  /* Snap values of level set function if need be */
  if ( !PMMG_snpval_ls(parmesh,mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem with implicit function. Exit program.\n");
    return 0;
  }

//   if ( !MMG3D_hashTetra(mesh,1) ) {
//     fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
//     return 0;
//   }

//   /* Compatibility triangle orientation w/r tetras */
//   if ( !MMG5_bdryPerm(mesh) ) {
//     fprintf(stderr,"\n  ## Boundary orientation problem. Exit program.\n");
//     return 0;
//   }

//   if ( !MMG5_chkBdryTria(mesh) ) {
//     fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
//     return 0;
//   }

//   /* Build hash table for initial edges */
//   if ( !MMG5_hGeom(mesh) ) {
//     fprintf(stderr,"\n  ## Hashing problem (0). Exit program.\n");
//     return 0;
//   }

//   if ( !MMG5_bdrySet(mesh) ) {
//     fprintf(stderr,"\n  ## Problem in setting boundary. Exit program.\n");
//     return 0;
//   }

//   /* Reset the mesh->info.isoref field everywhere it appears */
//   if ( !MMG3D_resetRef(mesh) ) {
//     fprintf(stderr,"\n  ## Problem in resetting references. Exit program.\n");
//     return 0;
//   }

//   /* Removal of small parasitic components */
//   if ( mesh->info.iso ) {
//     if ( mesh->info.rmc > 0. && !MMG3D_rmc(mesh,sol) ) {
//       fprintf(stderr,"\n  ## Error in removing small parasitic components."
//               " Exit program.\n");
//       return 0;
//     }
//   }
//   else {
//     /* RMC : on verra */
//     if ( mesh->info.rmc > 0 ) {
//       fprintf(stdout,"\n  ## Warning: rmc option not implemented for boundary"
//               " isosurface extraction.\n");
//     }
//   }

// #ifdef USE_POINTMAP
//   /* Initialize source point with input index */
//   MMG5_int ip;
//   for( ip = 1; ip <= mesh->np; ip++ )
//     mesh->point[ip].src = ip;
// #endif

//   if ( !MMG3D_cuttet(mesh,sol,met) ) {
//     fprintf(stderr,"\n  ## Problem in discretizing implicit function. Exit program.\n");
//     return 0;
//   }

//   MMG5_DEL_MEM(mesh,mesh->adja);
//   MMG5_DEL_MEM(mesh,mesh->adjt);
//   MMG5_DEL_MEM(mesh,mesh->tria);

//   mesh->nt = 0;

//   if ( !MMG3D_setref(mesh,sol) ) {
//     fprintf(stderr,"\n  ## Problem in setting references. Exit program.\n");
//     return 0;
//   }

//   /* Clean old bdy analysis */
//   for ( MMG5_int k=1; k<=mesh->np; ++k ) {
//     if ( mesh->point[k].tag & MG_BDY ) {
//       mesh->point[k].tag &= ~MG_BDY;
//     }
//   }

//   /* Clean memory */
//   MMG5_DEL_MEM(mesh,sol->m);

  return 1;
}

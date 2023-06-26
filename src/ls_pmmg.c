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
 *
 * \return 1 if success, 0 otherwise.
 *
 * \todo Fill the funtion
 *
 * Proceed to discretization of the implicit function carried by sol into mesh,
 * once values of sol have been snapped/checked
 *
 */
int PMMG_cuttet_ls(PMMG_pParMesh parmesh,MMG5_pMesh mesh, MMG5_pSol sol,MMG5_pSol met){
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n      ## PMMG_cuttet_ls:: Under development.\n");

  MMG5_pTetra   pt;
  MMG5_pxTetra  pxt;
  MMG5_pPoint   p0,p1;
  MMG5_Hash     hash  ;
  double        c[3],v0,v1,s;
  int           ier;
  MMG5_int      vx[6],k,ip0,ip1,np,nb,ns,ne,src,refext,refint;
  int8_t        ia,j,npneg;
  static int8_t mmgWarn = 0;

  /*********************/
  /* Reset point flags */
  /*********************/
  for (k=1; k<=mesh->np; k++)
    mesh->point[k].flag = 0;

  /*********************************************************/
  /* Compute the number nb of intersection points on edges */
  /*********************************************************/
  nb = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    for (ia=0; ia<6; ia++) {
      ip0 = pt->v[MMG5_iare[ia][0]];
      ip1 = pt->v[MMG5_iare[ia][1]];
      p0  = &mesh->point[ip0];
      p1  = &mesh->point[ip1];
      if ( p0->flag && p1->flag )  continue;
      v0  = sol->m[ip0]-mesh->info.ls;
      v1  = sol->m[ip1]-mesh->info.ls;
      if ( fabs(v0) > MMG5_EPSD2 && fabs(v1) > MMG5_EPSD2 && v0*v1 < 0.0 ) {
        if ( !p0->flag ) {
          p0->flag = ++nb;
        }
        if ( !p1->flag ) {
          p1->flag = ++nb;
        }
      }
    }
  }
  if ( ! nb )  return 1;


  /* TODO :: test if the number of point proc by proc is correct */
#ifndef NDEBUG
  /* TODO */
#endif

  /*************************************************************************/
  /* Create the hash table for the edges created by the LS discretization  */
  /*************************************************************************/
  if ( !MMG5_hashNew(mesh,&hash,nb,7*nb) ) return 0;

  /********************************************/
  /* Create intersection points at isovalue 0 */
  /* and set flags to tetras                  */
  /********************************************/
  /* Loop over all boundaries and required edges,
     and put ip = -1 in hash structure */
  // Loop over tetra
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) )  continue;

    /* Avoid split of edges belonging to a required tet */
    if ( pt->tag & MG_REQ ) {
      for (ia=0; ia<6; ia++) {
        ip0 = pt->v[MMG5_iare[ia][0]];
        ip1 = pt->v[MMG5_iare[ia][1]];
        np  = -1;
        if ( !MMG5_hashEdge(mesh,&hash,ip0,ip1,np) )  return -1;
      }
      continue;
    }

    // If the xtetra associated to this tetra exists
    if ( !pt->xt ) continue;

    // Point towards the xtetra corresponding to the tetra...
    pxt = &mesh->xtetra[pt->xt];
    // ... then loop over the faces
    for (ia=0; ia<4; ia++) {
      if ( !(pxt->ftag[ia] & MG_BDY) ) continue;

      for (j=0; j<3; j++) {
        if ( !(pxt->tag[ MMG5_iarf[ia][j] ] & MG_REQ) ) continue;

        ip0 = pt->v[MMG5_idir[ia][MMG5_inxt2[j]]];
        ip1 = pt->v[MMG5_idir[ia][MMG5_iprv2[j]]];
        np  = -1;
        if ( !MMG5_hashEdge(mesh,&hash,ip0,ip1,np) )  return -1;
      }
    }
  }


//   for (k=1; k<=mesh->ne; k++) {
//     pt = &mesh->tetra[k];
//     if ( !MG_EOK(pt) )  continue;

//     for (ia=0; ia<6; ia++) {
//       ip0 = pt->v[MMG5_iare[ia][0]];
//       ip1 = pt->v[MMG5_iare[ia][1]];
//       np  = MMG5_hashGet(&hash,ip0,ip1);

//       if ( np>0 )  continue;

//       if ( !MMG5_isSplit(mesh,pt->ref,&refint,&refext) ) continue;

//       p0 = &mesh->point[ip0];
//       p1 = &mesh->point[ip1];
//       v0 = sol->m[ip0]-mesh->info.ls;
//       v1 = sol->m[ip1]-mesh->info.ls;
//       if ( fabs(v0) < MMG5_EPSD2 || fabs(v1) < MMG5_EPSD2 )  continue;
//       else if ( MG_SMSGN(v0,v1) )  continue;
//       else if ( !p0->flag || !p1->flag )  continue;

//       npneg = (np<0);

//       s = v0 / (v0-v1);

//       s = MG_MAX(MG_MIN(s,1.0-MMG5_EPS),MMG5_EPS);
//       c[0] = p0->c[0] + s*(p1->c[0]-p0->c[0]);
//       c[1] = p0->c[1] + s*(p1->c[1]-p0->c[1]);
//       c[2] = p0->c[2] + s*(p1->c[2]-p0->c[2]);

// #ifdef USE_POINTMAP
//       src = p0->src;
// #else
//       src = 1;
// #endif
//       np = MMG3D_newPt(mesh,c,0,src);
//       if ( !np ) {
//         MMG5_int oldnpmax = mesh->npmax;
//         MMG3D_POINT_REALLOC(mesh,sol,np,MMG5_GAP,
//                              fprintf(stderr,"\n  ## Error: %s: unable to"
//                                      " allocate a new point\n",__func__);
//                              MMG5_INCREASE_MEM_MESSAGE();
//                              return 0
//                              ,c,0,src);
//         if( met ) {
//           if( met->m ) {
//             MMG5_ADD_MEM(mesh,(met->size*(mesh->npmax-met->npmax))*sizeof(double),
//                          "larger solution",
//                          MMG5_SAFE_RECALLOC(mesh->point,mesh->npmax+1,oldnpmax+1,MMG5_Point,,);
//                          mesh->memCur -= (mesh->npmax - oldnpmax)*sizeof(MMG5_Point);
//                          mesh->npmax = oldnpmax;
//                          mesh->np = mesh->npmax-1;
//                          mesh->npnil = 0;
//                          return 0);
//             MMG5_SAFE_REALLOC(met->m,met->size*(met->npmax+1),
//                               met->size*(mesh->npmax+1),
//                               double,"larger solution",
//                               MMG5_SAFE_RECALLOC(mesh->point,mesh->npmax+1,oldnpmax+1,MMG5_Point,,);
//                               mesh->memCur -= (mesh->npmax - oldnpmax)*sizeof(MMG5_Point);
//                               mesh->npmax = oldnpmax;
//                               mesh->np = mesh->npmax-1;
//                               mesh->npnil = 0;
//                               return 0);
//           }
//           met->npmax = mesh->npmax;
//         }
//       }
//       sol->m[np] = mesh->info.ls;
//       /* If user provide a metric, interpolate it at the new point */
//       if ( met && met->m ) {
//         if ( met->size > 1 ) {
//           ier = MMG3D_intmet33_ani(mesh,met,k,ia,np,s);
//         }
//         else {
//           ier = MMG5_intmet_iso(mesh,met,k,ia,np,s);
//         }
//         if ( ier <= 0 ) {
//           // Unable to compute the metric
//           fprintf(stderr,"\n  ## Error: %s: unable to"
//                   " interpolate the metric during the level-set"
//                   " discretization\n",__func__);
//           return 0;
//         }
//       }

//       if ( npneg ) {
//         /* We split a required edge */
//         if ( !mmgWarn ) {
//           mmgWarn = 1;
//           fprintf(stderr,"  ## Warning: %s: the level-set intersect at least"
//                   " one required entity. Required entity ignored.\n\n",__func__);
//         }
//         MMG5_hashUpdate(&hash,ip0,ip1,np);
//       }
//       else
//         MMG5_hashEdge(mesh,&hash,ip0,ip1,np);
//     }
//   }

//   /* Proceed to splitting, according to flags to tets */
//   ne  = mesh->ne;
//   ns  = 0;
//   ier = 1;
  // for (k=1; k<=ne; k++) {
  //   pt = &mesh->tetra[k];
  //   if ( !MG_EOK(pt) )  continue;
  //   pt->flag = 0;
  //   memset(vx,0,6*sizeof(MMG5_int));
  //   for (ia=0; ia<6; ia++) {
  //     vx[ia] = MMG5_hashGet(&hash,pt->v[MMG5_iare[ia][0]],pt->v[MMG5_iare[ia][1]]);
  //     if ( vx[ia] > 0 )  MG_SET(pt->flag,ia);
  //   }
  //   switch (pt->flag) {
  //   case 1: case 2: case 4: case 8: case 16: case 32: /* 1 edge split */
  //     ier = MMG5_split1(mesh,met,k,vx,1);
  //     ns++;
  //     break;

  //   case 48: case 24: case 40: case 6: case 34: case 36:
  //   case 20: case 5: case 17: case 9: case 3: case 10: /* 2 edges (same face) split */
  //     ier = MMG5_split2sf(mesh,met,k,vx,1);
  //     ns++;
  //     break;

  //   case 7: case 25: case 42: case 52: /* 3 edges on conic configuration splitted */
  //     ier = MMG5_split3cone(mesh,met,k,vx,1);
  //     ns++;
  //     break;

  //   case 30: case 45: case 51:
  //     ier = MMG5_split4op(mesh,met,k,vx,1);
  //     ns++;
  //     break;

  //   default :
  //     assert(pt->flag == 0);
  //     break;
  //   }
  //   if ( !ier ) return 0;
  // }
  // if ( (mesh->info.ddebug || abs(mesh->info.imprim) > 5) && ns > 0 )
  //   fprintf(stdout,"     %7" MMG5_PRId " splitted\n",ns);

  // MMG5_DEL_MEM(mesh,hash.item);
  // return ns;
  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param mesh pointer toward the mesh
 * \param sol pointer toward the level-set
 *
 * \return 1 if success, 0 otherwise
 *
 * \todo Fill the funtion
 *
 * Removal of small parasitic components (bubbles of material, etc) with volume
 * less than mesh->info.rmc (default VOLFRAC) * volume of the mesh.
 *
 */
int PMMG_rmc(PMMG_pParMesh parmesh,MMG5_pMesh mesh,MMG5_pSol sol){
  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"\n      ## TODO:: PMMG_rmc.\n");
  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure
 * \param mesh pointer toward the mesh structure.
 * \param sol pointer toward the level-set function.
 *
 * \return 1 if success, 0 if fail.
 *
 * \todo Fill the funtion
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
 *
 * \return 0 if fail, 1 otherwise.
 *
 * Create implicit surface in mesh.
 *
 */
int PMMG_ls(PMMG_pParMesh parmesh, MMG5_pMesh mesh,MMG5_pSol sol,MMG5_pSol met) {
  char str[16]="";

  /* Set function pointers */
  /** \todo TODO :: Surface ls and alias functions */
  if ( mesh->info.isosurf ) {
    fprintf(stderr," ## Error: Splitting boundaries on isovalue not yet"
            " implemented. Exit program.\n");
    return 0;
  }

  if ( parmesh->info.imprim > PMMG_VERB_VERSION )
    fprintf(stdout,"  ** ISOSURFACE EXTRACTION %s\n",str);

  if ( mesh->nprism || mesh->nquad ) {
    fprintf(stderr,"\n  ## Error: Isosurface extraction not available with"
            " hybrid meshes. Exit program.\n");
    return 0;
  }

  /** \todo TODO :: Snap values of level set function if needed */
  if ( !PMMG_snpval_ls(parmesh,mesh,sol) ) {
    fprintf(stderr,"\n  ## Problem with implicit function. Exit program.\n");
    return 0;
  }

  /* OK :: Create table of adjacency for tetra */
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  /* OK :: Check the compatibility of triangle orientation with tetra faces */
  if ( !MMG5_bdryPerm(mesh) ) {
    fprintf(stderr,"\n  ## Boundary orientation problem. Exit program.\n");
    return 0;
  }

  /* TO BE CHECKED :: Check behaviour with PMMG_APIDISTRIB_nodes
    Identify surface mesh
    Clean triangle array - remove useless or double triangles
    and add the missing ones */
  if ( !MMG5_chkBdryTria(mesh) ) {
    fprintf(stderr,"\n  ## Boundary problem. Exit program.\n");
    return 0;
  }

  /* OK :: Build hash table for initial edges */
  if ( !MMG5_hGeom(mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem (0). Exit program.\n");
    return 0;
  }

  /* OK :: Set the triangles references to the tetrahedra faces and edges */
  if ( !MMG5_bdrySet(mesh) ) {
    fprintf(stderr,"\n  ## Problem in setting boundary. Exit program.\n");
    return 0;
  }

  /* OK :: Reset the mesh->info.isoref field everywhere */
  if ( !MMG3D_resetRef_ls(mesh) ) {
    fprintf(stderr,"\n  ## Problem in resetting references. Exit program.\n");
    return 0;
  }

  /** \todo TODO :: Removal of small parasitic components */
  if ( mesh->info.rmc > 0 ) {
    PMMG_rmc(parmesh,mesh,sol);
    fprintf(stdout,"\n  ## Warning: rmc option not implemented yet for ParMmg\n");
    return 0;
  }

#ifdef USE_POINTMAP
  /* OK :: Initialize source point with input index */
  MMG5_int ip;
  for( ip = 1; ip <= mesh->np; ip++ )
    mesh->point[ip].src = ip;
#endif

  /* OK :: Compute vertices and triangles global numerotation */
  if ( !PMMG_Compute_verticesGloNum( parmesh,parmesh->comm ) ) {
    if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
      fprintf(stdout,"\n\n\n  -- WARNING: IMPOSSIBLE TO COMPUTE NODE GLOBAL NUMBERING\n\n\n");
      PMMG_RETURN_AND_FREE( parmesh, PMMG_LOWFAILURE );
    }
  }

  if ( !PMMG_Compute_trianglesGloNum( parmesh,parmesh->comm ) ) {
    if ( parmesh->info.imprim > PMMG_VERB_VERSION ) {
      fprintf(stdout,"\n\n\n  -- WARNING: IMPOSSIBLE TO COMPUTE TRIANGLE GLOBAL NUMBERING\n\n\n");
      PMMG_RETURN_AND_FREE( parmesh, PMMG_LOWFAILURE );
    }
  }

  /** \todo TODO :: Discretization of the implicit funtion - Cut tetra */
  if ( !PMMG_cuttet_ls(parmesh,mesh,sol,met) ) {
    fprintf(stderr,"\n  ## Problem in discretizing implicit function. Exit program.\n");
    return 0;
  }

  /* Not sure which function to be used to deallocate memory */
  MMG5_DEL_MEM(mesh,mesh->adja);
  MMG5_DEL_MEM(mesh,mesh->adjt);
  MMG5_DEL_MEM(mesh,mesh->tria);

  mesh->nt = 0;

  /* OK :: Set ref to tetra according to the sign of the level-set */
  /* Comment for now as the level-set has not been performed yet it fails */
  // if ( !MMG3D_setref_ls(mesh,sol) ) {
  //   fprintf(stderr,"\n  ## Problem in setting references. Exit program.\n");
  //   return 0;
  // }

  /* Clean old bdy analysis */
  for ( MMG5_int k=1; k<=mesh->np; ++k ) {
    if ( mesh->point[k].tag & MG_BDY ) {
      mesh->point[k].tag &= ~MG_BDY;
    }
  }

  /* Clean memory */
  MMG5_DEL_MEM(mesh,sol->m);

  return 1;
}

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
 * \file interpmesh_pmmg.c
 * \brief Interpolate data from a background mesh to the current mesh.
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \author Nikos Pattakos (Inria)
 * \author Luca Cirrottola (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */
#include "parmmg.h"
#include "locate_pmmg.h"

/**
 * \param m pointer to the 2x2 symmetric matrix
 * \param im pointer to the inverse matrix
 *
 * \return 0 if fail, 1 if success
 *
 *  Invert 2x2 symmetric matrix.
 */
int PMMG_invmat22( double *m, double *im ) {
  double det;

  det = m[0]*m[2] - m[1]*m[1];
  if ( fabs(det) < MMG5_EPS*MMG5_EPS ) {
    fprintf(stderr,"\n  ## Error: %s: null metric det : %E \n",
            __func__,det);
    return 0;
  }
  det = 1.0 / det;

  im[0] =  det*m[3];
  im[1] = -det*m[1];
  im[2] =  det*m[0];

  return 1;
}

/**
 * \param mesh pointer to the current mesh
 * \param met pointer to the current metrics
 * \param oldMet pointer to the background metrics
 * \param ptr pointer to the target background triangle
 * \param ip index of the current point
 * \param l local index of the edge on the background triangle
 * \param barycoord barycentric coordinates of the point to be interpolated
 *
 * \return 0 if fail, 1 if success
 *
 *  Linearly interpolate point metrics on a target background edge.
 */
int PMMG_interp2bar_iso( MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,
                         MMG5_pTria ptr,int ip,int l,PMMG_barycoord *barycoord ) {
  double alpha;

  assert( met->size == 1 );

  alpha = barycoord[l].val;

  /** Linear interpolation of the squared size */
  met->m[ip] =      alpha *oldMet->m[ptr->v[MMG5_inxt2[l]]] +
               (1.0-alpha)*oldMet->m[ptr->v[MMG5_iprv2[l]]];

  return 1;
}

/**
 * \param mesh pointer to the current mesh
 * \param met pointer to the current metrics
 * \param oldMet pointer to the background metrics
 * \param ptr pointer to the target background triangle
 * \param ip index of the current point
 * \param l local index of the edge on the background triangle
 * \param barycoord barycentric coordinates of the point to be interpolated
 *
 * \return 0 if fail, 1 if success
 *
 *  Linearly interpolate point the metrics inverse on a target background
 *  edge.
 *
 */
int PMMG_interp2bar_ani( MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,
                         MMG5_pTria ptr,int ip,int l,PMMG_barycoord *barycoord ) {
  double alpha,dm[3][3],mi[3][3],m[3];
  int    iloc,i,isize,nsize,ier;

  assert( met->size == 3 );
  nsize  = met->size;

  alpha = barycoord[l].val;

  for( i=0; i<2; i++ ) {
    for(isize = 0; isize < nsize; isize++ )
      dm[i][isize] = oldMet->m[nsize*ptr->v[i]+isize];
    if( !PMMG_invmat22(dm[i],mi[i]) ) return 0;
  }

  /** Linear interpolation of the metrics */
  for( isize = 0; isize < nsize; isize++ ) {
    m[isize] =      alpha *mi[0][isize]+
               (1.0-alpha)*mi[1][isize];
  }

  if( !PMMG_invmat22(m,mi[0]) ) return 0;
  for( isize = 0; isize < nsize; isize++ )
    met->m[nsize*ip+isize] = mi[0][isize];

  return 1;
}

/**
 * \param mesh pointer to the current mesh
 * \param met pointer to the current metrics
 * \param oldMet pointer to the background metrics
 * \param ptr pointer to the target background triangle
 * \param ip index of the current point
 * \param barycoord barycentric coordinates of the point to be interpolated
 *
 * \return 0 if fail, 1 if success
 *
 *  Linearly interpolate point metrics on a target background triangle.
 *  This function is analogous to the MMG5_interp4bar_iso() function.
 */
int PMMG_interp3bar_iso( MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,
                         MMG5_pTria ptr,int ip,PMMG_barycoord *barycoord ) {
  double phi[3];
  int iloc,i,ier;

  assert( met->size == 1 );

  PMMG_barycoord_get( phi, barycoord, 3 );

  /** Linear interpolation of the squared size */
  met->m[ip] = phi[0]*oldMet->m[ptr->v[0]] +
               phi[1]*oldMet->m[ptr->v[1]] + 
               phi[2]*oldMet->m[ptr->v[2]];

  return 1;
}

/**
 * \param mesh pointer to the current mesh
 * \param met pointer to the current metrics
 * \param oldMet pointer to the background metrics
 * \param ptr pointer to the target background triangle
 * \param ip index of the current point
 * \param barycoord barycentric coordinates of the point to be interpolated
 *
 * \return 0 if fail, 1 if success
 *
 *  Linearly interpolate point the metrics inverse on a target background
 *  triangle.
 *  This function is analogous to the MMG5_interp4barintern() function.
 *
 */
int PMMG_interp3bar_ani( MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,
                         MMG5_pTria ptr,int ip,PMMG_barycoord *barycoord ) {
  double phi[3],dm[3][3],mi[3][3],m[3];
  int    iloc,i,isize,nsize,ier;

  assert( met->size == 3 );
  nsize  = met->size;

  PMMG_barycoord_get( phi, barycoord, 3 );

  for( i=0; i<3; i++ ) {
    for(isize = 0; isize < nsize; isize++ )
      dm[i][isize] = oldMet->m[nsize*ptr->v[i]+isize];
    if( !PMMG_invmat22(dm[i],mi[i]) ) return 0;
  }

  /** Linear interpolation of the metrics */
  for( isize = 0; isize < nsize; isize++ ) {
    m[isize] = phi[0]*mi[0][isize]+
               phi[1]*mi[1][isize]+
               phi[2]*mi[2][isize];
  }

  if( !PMMG_invmat22(m,mi[0]) ) return 0;
  for( isize = 0; isize < nsize; isize++ )
    met->m[nsize*ip+isize] = mi[0][isize];

  return 1;
}

/**
 * \param mesh pointer to the current mesh
 * \param met pointer to the current metrics
 * \param oldMet pointer to the background metrics
 * \param pt pointer to the target background tetrahedron
 * \param ip index of the current point
 * \param barycoord barycentric coordinates of the point to be interpolated
 *
 * \return 0 if fail, 1 if success
 *
 *  Linearly interpolate point metrics on a target background tetrahedron.
 *  This function is analogous to the MMG5_interp4bar_iso() function.
 */
int PMMG_interp4bar_iso( MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,
                         MMG5_pTetra pt,int ip,PMMG_barycoord *barycoord ) {
  double phi[4];
  int iloc,i,ier;

  assert( met->size == 1 );

  PMMG_barycoord_get( phi, barycoord, 4 );

  /** Linear interpolation of the squared size */
  met->m[ip] = phi[0]*oldMet->m[pt->v[0]]+
               phi[1]*oldMet->m[pt->v[1]]+
               phi[2]*oldMet->m[pt->v[2]]+
               phi[3]*oldMet->m[pt->v[3]];

  return 1;
}

/**
 * \param mesh pointer to the current mesh
 * \param met pointer to the current metrics
 * \param oldMet pointer to the background metrics
 * \param pt pointer to the target background tetrahedron
 * \param ip index of the current point
 * \param barycentric barycentric coordinates of the point to be interpolated
 *
 * \return 0 if fail, 1 if success
 *
 *  Linearly interpolate point the metrics inverse on a target background
 *  tetrahedron.
 *  This function is analogous to the MMG5_interp4barintern() function.
 *
 */
int PMMG_interp4bar_ani( MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,
                         MMG5_pTetra pt,int ip,PMMG_barycoord *barycoord ) {
  double phi[4],dm[4][6],mi[4][6],m[6];
  int    iloc,i,isize,nsize,ier;

  assert( met->size == 6 );
  nsize  = met->size;

  PMMG_barycoord_get( phi, barycoord, 4 );

  for( i=0; i<4; i++ ) {
    for(isize = 0; isize < nsize; isize++ )
      dm[i][isize] = oldMet->m[nsize*pt->v[i]+isize];
    if( !MMG5_invmat(dm[i],mi[i]) ) return 0;
  }

  /** Linear interpolation of the metrics */
  for( isize = 0; isize < nsize; isize++ ) {
    m[isize] = phi[0]*mi[0][isize] + phi[1]*mi[1][isize] +
               phi[2]*mi[2][isize] + phi[3]*mi[3][isize];
  }

  if( !MMG5_invmat(m,mi[0]) ) return 0;
  for( isize = 0; isize < nsize; isize++ )
    met->m[nsize*ip+isize] = mi[0][isize];

  return 1;
}

/**
 * \param mesh pointer to the current mesh.
 * \param met pointer to the current metrics.
 * \param oldMesh pointer to the background mesh.
 * \param oldMet pointer to the background metrics.
 * \param idest index of the target point.
 * \param isrc index of the source point.
 *
 * \return 0 if fail, 1 if success
 *
 * Copy the metric of a point.
 *
 */
int PMMG_copyMetrics( MMG5_pMesh mesh,MMG5_pSol met,MMG5_pMesh oldMesh,
                      MMG5_pSol oldMet,int idest,int isrc ) {
  int isize,nsize;

  nsize = met->size;

  for( isize = 0; isize<nsize; isize++ ) {
    met->m[nsize*idest+isize] = oldMet->m[nsize*isrc+isize];
  }

  return 1;
}

/**
 * \param mesh pointer to the current mesh.
 * \param oldMesh pointer to the background mesh.
 * \param met pointer to the current metrics.
 * \param oldMet pointer to the background metrics.
 * \param permNodGlob permutation array for nodes.
 * \param inputMet 1 if user provided metric.
 *
 * \return 0 if fail, 1 if success
 *
 * Copy the metric of a freezed interface point.
 *
 */
int PMMG_copyMetrics_point( MMG5_pMesh mesh,MMG5_pMesh oldMesh,
                            MMG5_pSol met,MMG5_pSol oldMet,int* permNodGlob,
                            unsigned char inputMet) {
  MMG5_pPoint    ppt;
  int            isize,nsize,ip;

  if ( (!inputMet) || mesh->info.hsiz > 0.0 ) return 1;

  nsize   = met->size;

#warning Luca: when surface adapt will be ready, distinguish BDY from PARBDY

  /** Freezed points: Copy the  metrics  */
  if ( (!oldMesh->info.renum) || !permNodGlob ) {
    /* No permutation array: simple copy */
    for( ip = 1; ip <= oldMesh->np; ++ip ) {
      ppt = &oldMesh->point[ip];
      if( !MG_VOK(ppt) ) continue;

      ppt->flag = mesh->base;
      if( ppt->tag & MG_REQ ) {
        for( isize = 0; isize<nsize; isize++ ) {
          met->m[nsize*ip+isize] = oldMet->m[nsize*ip+isize];
        }
      }
    }
  }
  else {
    /* Due to the scotch renumbering we must copy from the old mesh to the
     * new one */
    for( ip = 1; ip <= oldMesh->np; ++ip ) {
      ppt = &oldMesh->point[ip];
      if( !MG_VOK(ppt) ) continue;

      ppt->flag = mesh->base;
      if( ppt->tag & MG_REQ ) {
        for( isize = 0; isize<nsize; isize++ ) {
          met->m[nsize*permNodGlob[ip]+isize] = oldMet->m[nsize*ip+isize];
        }
      }
    }
  }

  return 1;
}

/**
 * \param mesh pointer to the current mesh structure.
 * \param triaNormals pointer to the array of non-normalized triangle normals.
 *
 * \return 1.
 *
 *  Precompute non-normalized triangle normals.
 *
 */
int PMMG_precompute_triaNormals( MMG5_pMesh mesh,double *triaNormals ) {
  MMG5_pTria  ptr;
  double      *normal,dd;
  int         k,ia,ib,ic;
  int         ier;

  for( k = 1; k <= mesh->nt; k++ ) {
    ptr = &mesh->tria[k];
    ia = ptr->v[0];
    ib = ptr->v[1];
    ic = ptr->v[2];
    normal = &triaNormals[3*k];
    /* Store triangle unit normal and volume */
    ier = MMG5_nonUnitNorPts( mesh,ia,ib,ic,normal );
    ptr->qual = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
    dd = 1.0/ptr->qual;
    normal[0] *= dd;
    normal[1] *= dd;
    normal[2] *= dd;
  }

  return 1;
}

/**
 * \param mesh pointer to the current mesh structure.
 * \param faceAreas pointer to the array of oriented face areas.
 *
 * \return 1.
 *
 *  Precompute oriented face areas on tetrahedra.
 *
 */
int PMMG_precompute_faceAreas( MMG5_pMesh mesh,double *faceAreas ) {
  MMG5_pTetra pt;
  double      *normal;
  int         ie,ifac,ia,ib,ic;
  int         ier;

  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    /* Store tetra volume in the qual field */
    pt->qual = MMG5_orvol( mesh->point, pt->v );
    /* Store oriented face normals */
    for( ifac = 0; ifac < 4; ifac++ ) {
      normal = &faceAreas[12*ie+3*ifac];
      ia = pt->v[MMG5_idir[ifac][0]];
      ib = pt->v[MMG5_idir[ifac][1]];
      ic = pt->v[MMG5_idir[ifac][2]];
      ier = MMG5_nonUnitNorPts( mesh,ia,ib,ic,normal );
    }
  }

  return 1;
}

/**
 * \param mesh pointer to the current mesh structure.
 * \param oldMesh pointer to the background mesh structure.
 * \param met pointer to the current metrics structure.
 * \param oldMet pointer to the background metrics structure.
 * \param faceAreas pointer to the array of oriented face areas.
 * \param triaNormals pointer to the array of non-normalized triangle normals.
 * \param permNodGlob permutation array of nodes.
 * \param inputMet 1 if user provided metric.
 * \param myrank process rank.
 * \param igrp current mesh group.
 * \param locStats pointer to the localization statistics structure.
 *
 * \return 0 if fail, 1 if success
 *
 *  Interpolate metrics for all groups from background to current meshes.
 *  Do nothing if no metrics is provided (inputMet == 0), otherwise:
 *  - if the metrics is constant, recompute it;
 *  - else, interpolate the non-constant metrics.
 *
 *  Oriented face areas are pre-computed in this function before proceeding
 *  with the localization.
 *
 */
int PMMG_interpMetrics_mesh( MMG5_pMesh mesh,MMG5_pMesh oldMesh,
                             MMG5_pSol met,MMG5_pSol oldMet,
                             double *faceAreas,double *triaNormals,
                             int *permNodGlob,unsigned char inputMet,
                             int myrank,int igrp,PMMG_locateStats *locStats ) {
  MMG5_pTetra pt;
  MMG5_pPoint ppt;
  PMMG_barycoord barycoord[4];
  double      *normal,dd;
  int         istartTetra,istartTria,ifoundTetra,ifoundTria;
  int         ifoundEdge,ifoundVertex;
  int         ip,ie,ifac,k,ia,ib,ic,iloc;
  int         ier;
  static int  mmgWarn=0;

  if( inputMet != 1 ) {

    /* Nothing to do */
    return 1;

  } else {

    if( mesh->info.hsiz > 0.0 ) {

      /* Compute constant metrics */
      if ( !MMG3D_Set_constantSize(mesh,met) ) return 0;

    } else {

      /** Pre-compute oriented face areas */
      ier = PMMG_precompute_faceAreas( oldMesh,faceAreas );

      /** Pre-compute surface unit normals */
      ier = PMMG_precompute_triaNormals( oldMesh,triaNormals );

      /** Interpolate metrics */
      oldMesh->base = 0;
      for ( ie = 1; ie < oldMesh->ne+1; ie++ ) {
        pt = &oldMesh->tetra[ie];
        if ( !MG_EOK(pt) ) continue;
        pt->flag = oldMesh->base;
      }

#ifdef NO_POINTMAP
      ifoundTetra = ifoundTria = 1;
#else
      PMMG_locate_setStart( mesh,oldMesh );
#endif
      /* Loop on new tetrahedra, and localize their vertices in the old mesh */
      mesh->base++;
      for( ie = 1; ie <= mesh->ne; ie++ ) {
        pt = &mesh->tetra[ie];
        if( !MG_EOK(pt) ) continue;
        for( iloc = 0; iloc < 4; iloc++ ) {
          ip = pt->v[iloc];
          ppt = &mesh->point[ip];
          if( !MG_VOK(ppt) ) continue;

          /* Skip already interpolated points */
          if( ppt->flag == mesh->base ) continue;

          if( ppt->tag & MG_REQ ) {
            /* Flag point as interpolated */
            ppt->flag = mesh->base;
            continue; // treated by copyMetric_points
          } else if ( ppt->tag & MG_BDY ) {

#ifdef NO_POINTMAP
            istartTria = ifoundTria;
#else
            istartTria = ppt->s;
#endif
            /** Locate point in the old mesh */
            ifoundTria = PMMG_locatePointBdy( oldMesh, ppt, istartTria,
                                              triaNormals, barycoord, ip,
                                              &ifoundEdge, &ifoundVertex );

            ifoundTria = PMMG_locatePoint_errorCheck( mesh,ip,ifoundTria,
                                                      myrank,igrp );

            /** Interpolate point metrics */
            if( ifoundVertex != PMMG_UNSET ) {
              ier = PMMG_copyMetrics( mesh,met,oldMesh,oldMet,ip,
                                      oldMesh->tria[ifoundTria].v[ifoundVertex] );
            } else if( ifoundEdge != PMMG_UNSET ) {
              ier = PMMG_interp2bar( mesh,met,oldMet,&oldMesh->tria[ifoundTria],
                                     ip,ifoundEdge,barycoord );
            } else {
            ier = PMMG_interp3bar(mesh,met,oldMet,&oldMesh->tria[ifoundTria],ip,
                                  barycoord);
            }

            /* Flag point as interpolated */
            ppt->flag = mesh->base;

          } else {

#ifdef NO_POINTMAP
            istartTetra = ifoundTetra;
#else
            istartTetra = ppt->s;
#endif
            /** Locate point in the old volume mesh */
            ifoundTetra = PMMG_locatePointVol( oldMesh, ppt, istartTetra,
                                               faceAreas, barycoord, ip );

            ifoundTetra = PMMG_locatePoint_errorCheck( mesh,ip,ifoundTetra,
                                                       myrank,igrp );

            /** Interpolate volume point metrics */
            ier = PMMG_interp4bar(mesh,met,oldMet,&oldMesh->tetra[ifoundTetra],ip,
                                  barycoord);

            /* Flag point as interpolated */
            ppt->flag = mesh->base;

          }
        }
      }
#ifndef NDEBUG
      PMMG_locate_postprocessing( mesh,oldMesh,locStats );
#endif
    }
  }

  return 1;
}

/**
 * \param parmesh pointer to the parmesh structure.
 * \param permNodGlob permutation array of nodes.
 *
 * \return 0 if fail, 1 if success
 *
 *  Interpolate metrics for all groups from background to current meshes.
 *  Do nothing if no metrics is provided (info.inputMet == 0), otherwise:
 *  - if the metrics is constant, recompute it;
 *  - else, interpolate the non-constant metrics.
 *
 */
int PMMG_interpMetrics( PMMG_pParMesh parmesh,int *permNodGlob ) {
  PMMG_pGrp        grp,oldGrp;
  MMG5_pMesh       mesh,oldMesh;
  MMG5_pSol        met,oldMet;
  MMG5_Hash        hash;
  PMMG_locateStats *locStats;
  size_t           memAv,oldMemMax;
  double           *faceAreas,*triaNormals;
  int              igrp,ier;

  PMMG_TRANSFER_AVMEM_TO_PARMESH(parmesh,memAv,oldMemMax);

  PMMG_MALLOC( parmesh,locStats,parmesh->ngrp,PMMG_locateStats,"locStats",return 0 );

  /** Loop on current groups */
  ier = 1;
  for( igrp = 0; igrp < parmesh->ngrp; igrp++ ) {

    grp  = &parmesh->listgrp[igrp];
    mesh = grp->mesh;
    met  = grp->met;

    oldGrp  = &parmesh->old_listgrp[igrp];
    oldMesh = oldGrp->mesh;
    oldMet  = oldGrp->met;

    /** Pre-allocate oriented face areas and surface unit normals */
    if( ( parmesh->info.inputMet == 1 ) && ( mesh->info.hsiz <= 0.0 ) ) {
      PMMG_MALLOC( parmesh,faceAreas,12*(oldMesh->ne+1),double,"faceAreas",return 0 );
      PMMG_MALLOC( parmesh,triaNormals,3*(oldMesh->nt+1),double,"triaNormals",return 0 );
    }

    if( !PMMG_interpMetrics_mesh( mesh,oldMesh,met,oldMet,faceAreas,triaNormals,
                                  permNodGlob,parmesh->info.inputMet,
                                  parmesh->myrank,igrp,&locStats[igrp] ) ) ier = 0;

    /** Deallocate oriented face areas and surface unit normals */
    if( ( parmesh->info.inputMet == 1 ) && ( mesh->info.hsiz <= 0.0 ) ) {
      PMMG_DEL_MEM(parmesh,faceAreas,double,"faceAreas");
      PMMG_DEL_MEM(parmesh,triaNormals,double,"triaNormals");
    }

  }

#ifndef NDEBUG
  if( ( parmesh->info.inputMet == 1 ) && ( mesh->info.hsiz <= 0.0 ) )
    PMMG_locate_print( locStats,parmesh->ngrp,parmesh->myrank );
#endif
  PMMG_DEL_MEM(parmesh,locStats,PMMG_locateStats,"locStats");
  return ier;
}

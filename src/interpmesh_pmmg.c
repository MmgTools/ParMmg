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
#include "interpmesh_pmmg.h"

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
  int    i0,i1;
  double phi[3];

  assert( met->size == 1 );

  PMMG_barycoord_get( phi, barycoord, 3 );

  i0 = MMG5_inxt2[l];
  i1 = MMG5_iprv2[l];

  /** Linear interpolation of the squared size */
  met->m[ip] = phi[i0]*oldMet->m[ptr->v[i0]] +
               phi[i1]*oldMet->m[ptr->v[i1]];

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
  int    i0,i1;
  double phi[3],mi[2][6],mint[6];
  int    iloc,i,isize,nsize,ier;

  assert( met->size == 6 );
  nsize  = met->size;

  PMMG_barycoord_get( phi, barycoord, 3 );

  i0 = MMG5_inxt2[l];
  i1 = MMG5_iprv2[l];

  if( !MMG5_invmat( &oldMet->m[nsize*ptr->v[i0]], mi[0] ) ) return 0;
  if( !MMG5_invmat( &oldMet->m[nsize*ptr->v[i1]], mi[1] ) ) return 0;

  /** Linear interpolation of the metrics */
  for( isize = 0; isize < nsize; isize++ ) {
    mint[isize] = phi[i0]*mi[0][isize]+
                  phi[i1]*mi[1][isize];
  }

  if( !MMG5_invmat( mint, &met->m[nsize*ip] ) ) return 0;

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
  int iadr,i,j;

  iadr = met->size *ip;

  assert (mesh->npmax==met->npmax );

  PMMG_barycoord_get( phi, barycoord, 3 );

  /** Linear interpolation of the squared size */
  for ( j=0; j<met->size; ++j ) {
    met->m [ iadr + j ] = 0.0;
  }

  for( i=0; i<3; i++ ) {
    for ( j=0; j<met->size; ++j ) {
      /* Barycentric coordinates could be permuted */
      met->m[iadr+j] += phi[i]*oldMet->m[ptr->v[i]*met->size+j];
    }
  }

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
  double phi[3],mi[3][6],mint[6];
  int    iloc,i,isize,nsize,ier;

  assert( met->size == 6 );
  nsize  = met->size;

  PMMG_barycoord_get( phi, barycoord, 3 );

  for( i=0; i<3; i++ ) {
    if( !MMG5_invmat( &oldMet->m[nsize*ptr->v[i]], mi[i]) ) return 0;
  }

  /** Linear interpolation of the metrics */
  for( isize = 0; isize < nsize; isize++ ) {
    mint[isize] = phi[0]*mi[0][isize]+
                  phi[1]*mi[1][isize]+
                  phi[2]*mi[2][isize];
  }

  if( !MMG5_invmat( mint, &met->m[nsize*ip] ) ) return 0;

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
 *  Linearly interpolate point data (metric or solution) on a target background
 *  tetrahedron. This function is analogous to the MMG5_interp4bar_iso()
 *  function.
 */
int PMMG_interp4bar_iso( MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,
                         MMG5_pTetra pt,int ip,PMMG_barycoord *barycoord ) {
  double phi[4];
  int i,j,iadr;

  PMMG_barycoord_get( phi, barycoord, 4 );

  /** Linear interpolation of the squared size */
  iadr = met->size *ip;

  assert (mesh->npmax==met->npmax );

  for ( j=0; j<met->size; ++j ) {
    met->m [ iadr + j ] = 0.0;
  }

  for( i=0; i<4; i++ ) {
    for ( j=0; j<met->size; ++j ) {
      /* Barycentric coordinates could be permuted */
      met->m[iadr+j] += phi[i]*oldMet->m[pt->v[i]*met->size+j];
    }
  }

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
  double phi[4],mi[4][6],mint[6];
  int    i,isize,nsize;

  assert( met->size == 6 );
  nsize  = met->size;

  PMMG_barycoord_get( phi, barycoord, 4 );

  for( i=0; i<4; i++ ) {
    if( !MMG5_invmat( &oldMet->m[nsize*pt->v[i]], mi[i]) ) return 0;
  }

  /** Linear interpolation of the metrics */
  for( isize = 0; isize < nsize; isize++ ) {
    mint[isize] = phi[0]*mi[0][isize] + phi[1]*mi[1][isize] +
                  phi[2]*mi[2][isize] + phi[3]*mi[3][isize];
  }

  if( !MMG5_invmat( mint, &met->m[nsize*ip] ) ) return 0;

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
 * Copy the data stored in a solution structure of a freezed interface point.
 *
 */
static inline
int PMMG_copySol_point( MMG5_pMesh mesh,MMG5_pMesh oldMesh,
                        MMG5_pSol sol,MMG5_pSol oldSol,int* permNodGlob) {
  MMG5_pPoint    ppt;
  int            isize,nsize,ip;

  nsize   = sol->size;

#warning Luca: when surface adapt will be ready, distinguish BDY from PARBDY

  /** Freezed points: Copy the data stored in solution structure  */
  if ( (!oldMesh->info.renum) || !permNodGlob ) {
    /* No permutation array: simple copy */
    for( ip = 1; ip <= oldMesh->np; ++ip ) {
      ppt = &oldMesh->point[ip];
      if( !MG_VOK(ppt) ) continue;

      /* Remark: ppt->flag is setted multiple times to mesh->base if copySol is
       * called multiple times (for example if we copy both a metric and a
       * solution fields) but it is not a problem since we use the same mesh  */
      ppt->flag = mesh->base;
      if( ppt->tag & MG_REQ ) {
        for( isize = 0; isize<nsize; isize++ ) {
          sol->m[nsize*ip+isize] = oldSol->m[nsize*ip+isize];
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

      /* Remark: ppt->flag is setted multiple times to mesh->base if copySol is
       * called multiple times (for example if we copy both a metric and a
       * solution fields) but it is not a problem since we use the same mesh  */
      ppt->flag = mesh->base;
      if( ppt->tag & MG_REQ ) {
        for( isize = 0; isize<nsize; isize++ ) {
          sol->m[nsize*permNodGlob[ip]+isize] = oldSol->m[nsize*ip+isize];
        }
      }
    }
  }

  return 1;
}

/**
 * \param mesh pointer to the current mesh.
 * \param oldMesh pointer to the background mesh.
 * \param met pointer to the current metrics.
 * \param oldMet pointer to the background metrics.
 * \param permNodGlob permutation array for nodes.
 *
 * \return 0 if fail, 1 if success
 *
 * Copy the metric of a freezed interface point.
 *
 */
static
int PMMG_copyMetrics_point( MMG5_pMesh mesh,MMG5_pMesh oldMesh,
                            MMG5_pSol met,MMG5_pSol oldMet,int* permNodGlob,
                            uint8_t inputMet ) {
  int            ier;

  if ( !inputMet || mesh->info.hsiz > 0.0 ) return 1;

  ier =  PMMG_copySol_point( mesh, oldMesh,met,oldMet,permNodGlob);

  return ier;
}

/**
 * \param mesh pointer to the current mesh.
 * \param oldMesh pointer to the background mesh.
 * \param field pointer to the current fields.
 * \param oldField pointer to the background fields.
 * \param permNodGlob permutation array for nodes.
 *
 * \return 0 if fail, 1 if success
 *
 * Copy the metric of a freezed interface point.
 *
 */
static
int PMMG_copyFields_point( MMG5_pMesh mesh,MMG5_pMesh oldMesh,
                           MMG5_pSol field,MMG5_pSol oldField,int* permNodGlob) {
  MMG5_pSol      psl,oldPsl;
  int            j,ier;

  if ( !mesh->nsols ) return 1;

  for ( j=0; j<mesh->nsols; ++j ) {
    psl    =    field + j;
    oldPsl = oldField + j;
    ier =  PMMG_copySol_point( mesh, oldMesh,psl,oldPsl,permNodGlob);
    if ( !ier ) {
      return 0;
    }
  }

  return 1;
}

/**
 * \param mesh pointer to the current mesh.
 * \param oldMesh pointer to the background mesh.
 * \param met pointer to the current metrics.
 * \param oldMet pointer to the background metrics.
 * \param field pointer to the current fields.
 * \param oldField pointer to the background fields.
 * \param permNodGlob permutation array for nodes.
 * \param inputMet 1 if user provided metric.
 *
 * \return 0 if fail, 1 if success
 *
 * Copy the metric and fields of a freezed interface point.
 *
 */
int PMMG_copyMetricsAndFields_point( MMG5_pMesh mesh ,MMG5_pMesh oldMesh,
                                     MMG5_pSol  met  ,MMG5_pSol  oldMet,
                                     MMG5_pSol  field,MMG5_pSol  oldField,
                                     int* permNodGlob,uint8_t inputMet) {
  int ier;

  ier = PMMG_copyMetrics_point(mesh,oldMesh,met,oldMet,permNodGlob,inputMet);
  if ( !ier ) {
    return 0;
  }

  ier = PMMG_copyFields_point(mesh,oldMesh,field,oldField,permNodGlob);

  return ier;
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
 * Interpolate metrics and solution fields for all groups from background
 * to current meshes.
 * For the metric:
 *  Do nothing if no metrics is provided (inputMet == 0), otherwise:
 *  - if the metrics is constant, recompute it;
 *  - else, interpolate the non-constant metrics.
 *
 * For the solution fields: Do nothing if no solution field is provided
 *   (mesh->nsols == 0), interpolate the non-constant field otherwise.
 *
 *  Oriented face areas are pre-computed in this function before proceeding
 *  with the localization.
 *
 */
static
int PMMG_interpMetricsAndFields_mesh( MMG5_pMesh mesh,MMG5_pMesh oldMesh,
                                      MMG5_pSol met,MMG5_pSol oldMet,
                                      MMG5_pSol field,MMG5_pSol oldField,
                                      double *faceAreas,double *triaNormals,int *nodeTrias,
                                      int *permNodGlob,uint8_t inputMet,
                                      int myrank,int igrp,PMMG_locateStats *locStats ) {
  MMG5_pTetra pt;
  MMG5_pPoint ppt;
  MMG5_pSol   psl,oldPsl;
  PMMG_barycoord barycoord[4];
  double      *normal,dd;
  int         istartTetra,istartTria,ifoundTetra,ifoundTria;
  int         ifoundEdge,ifoundVertex;
  int         ip,ie,ifac,k,ia,ib,ic,iloc,nsols;
  int         ismet,ier,j;
  static int  mmgWarn=0;

  nsols = mesh->nsols;

  ismet = 1;
  if( inputMet != 1 ) {
    ismet = 0;
  }
  else  if( mesh->info.hsiz > 0.0 ) {

    /* Compute constant metrics */
    if ( !MMG3D_Set_constantSize(mesh,met) ) return 0;
    ismet = 0;

  }
  if ( (!ismet) && (!nsols) ) {

    /* Nothing to do */
    return 1;
  }

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

#ifndef USE_POINTMAP
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

#ifdef USE_POINTMAP
        ifoundTria = ppt->s;
#endif
        /** Locate point in the old mesh */
        ier = PMMG_locatePointBdy( oldMesh, ppt,
                                   triaNormals, nodeTrias, barycoord,
                                   &ifoundTria,&ifoundEdge, &ifoundVertex );

        if( mesh->info.imprim > PMMG_VERB_ITWAVES )
          PMMG_locatePoint_errorCheck( mesh,ip,ier,myrank,igrp );

        /** Interpolate point metrics */
        if( ismet ) {
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
        }

#warning Luca: make this part consistent with metrics interpolation
        /** Field interpolation */
        if ( mesh->nsols ) {
          for ( j=0; j<mesh->nsols; ++j ) {
            psl    = field + j;
            oldPsl = oldField + j;
            if ( oldPsl->size == 6 ) {
              /* Tensor field */
              ier = PMMG_interp3bar_ani(mesh,psl,oldPsl,
                                        &oldMesh->tria[ifoundTria],
                                        ip,barycoord);
            }
            else {
              /* Scalar or vector field */
              ier = PMMG_interp3bar_iso(mesh,psl,oldPsl,
                                        &oldMesh->tria[ifoundTria],
                                        ip,barycoord);
            }
          }
        }

        /* Flag point as interpolated */
        ppt->flag = mesh->base;

      } else {

#ifdef USE_POINTMAP
        ifoundTetra = ppt->s;
#endif
        /** Locate point in the old volume mesh */
        ier = PMMG_locatePointVol( oldMesh, ppt,
                                   faceAreas, barycoord, &ifoundTetra );

        if( mesh->info.imprim > PMMG_VERB_ITWAVES )
          PMMG_locatePoint_errorCheck( mesh,ip,ier,myrank,igrp );

        /** Interpolate volume point metrics */
        if( ismet ) {
          ier = PMMG_interp4bar(mesh,met,oldMet,&oldMesh->tetra[ifoundTetra],ip,
                                barycoord);
        }

        /** Field interpolation */
        if ( mesh->nsols ) {
          for ( j=0; j<mesh->nsols; ++j ) {
            psl    = field + j;
            oldPsl = oldField + j;
            if ( oldPsl->size == 6 ) {
              /* Tensor field */
              ier = PMMG_interp4bar_ani(mesh,psl,oldPsl,
                                        &oldMesh->tetra[ifoundTetra],
                                        ip,barycoord);
            }
            else {
              /* Scalar or vector field */
              ier = PMMG_interp4bar_iso(mesh,psl,oldPsl,
                                        &oldMesh->tetra[ifoundTetra],
                                        ip,barycoord);
            }
          }
        }

        /* Flag point as interpolated */
        ppt->flag = mesh->base;

      }
    }
  }
#ifndef NDEBUG
  PMMG_locate_postprocessing( mesh,oldMesh,locStats );
#endif

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
int PMMG_interpMetricsAndFields( PMMG_pParMesh parmesh,int *permNodGlob ) {
  PMMG_pGrp        grp,oldGrp;
  MMG5_pMesh       mesh,oldMesh;
  MMG5_pSol        met,oldMet,field,oldField;
  MMG5_Hash        hash;
  PMMG_locateStats *locStats,*mylocStats;
  double           *faceAreas,*triaNormals;
  int              *nodeTrias;
  int              igrp,ier;
  int8_t           allocated;

  locStats = NULL;
#ifndef NDEBUG
  /* Check that alloc size does not overflow PTRDIFF_MAX (max size
   * authorized by malloc) */
  if ( parmesh->ngrp*sizeof(PMMG_locateStats) + sizeof(size_t) > PTRDIFF_MAX ) {
    printf(" ## Warning (dev): %s: alloc size (%zu) for locStats exceeds the"
           " maximum object size (%zu).\n",__func__,
           parmesh->ngrp*sizeof(PMMG_locateStats) + sizeof(size_t),PTRDIFF_MAX);
  }
  else {
    PMMG_MALLOC( parmesh,locStats,parmesh->ngrp,PMMG_locateStats,"locStats", );
  }
#endif

  /** Loop on current groups */
  ier = 1;
  for( igrp = 0; igrp < parmesh->ngrp; igrp++ ) {

    grp  = &parmesh->listgrp[igrp];
    mesh = grp->mesh;
    met  = grp->met;
    field = grp->field;

    oldGrp  = &parmesh->old_listgrp[igrp];
    oldMesh = oldGrp->mesh;
    oldMet  = oldGrp->met;
    oldField = oldGrp->field;

    /** Pre-allocate oriented face areas and surface unit normals */
    allocated = 0;
    if ( mesh->nsols || (( parmesh->info.inputMet == 1 ) && ( mesh->info.hsiz <= 0.0 )) ) {
      PMMG_MALLOC( parmesh,faceAreas,12*(oldMesh->ne+1),double,"faceAreas",return 0 );
      PMMG_MALLOC( parmesh,triaNormals,3*(oldMesh->nt+1),double,"triaNormals",return 0 );
      PMMG_precompute_nodeTrias( parmesh,oldMesh,&nodeTrias );
      allocated = 1;
    }

    mylocStats = NULL;
    if ( locStats ) {
      mylocStats = locStats + igrp;
    }
    if( !PMMG_interpMetricsAndFields_mesh( mesh,oldMesh,met,oldMet,
                                           field,oldField,
                                           faceAreas,triaNormals,nodeTrias,
                                           permNodGlob,parmesh->info.inputMet,
                                           parmesh->myrank,igrp,mylocStats ) ) {
      ier = 0;
    }

    /** Deallocate oriented face areas and surface unit normals */
    if( allocated ) {
      PMMG_DEL_MEM(parmesh,faceAreas,double,"faceAreas");
      PMMG_DEL_MEM(parmesh,triaNormals,double,"triaNormals");
      PMMG_DEL_MEM(parmesh,nodeTrias,int,"nodeTrias");
    }

  }

#ifndef NDEBUG
  if( ((( parmesh->info.inputMet == 1 ) && ( mesh->info.hsiz <= 0.0 )) || mesh->nsols)
      && (parmesh->info.imprim0 > PMMG_VERB_DETQUAL) ) {
    PMMG_locate_print( locStats,parmesh->ngrp,parmesh->myrank );
  }
  PMMG_DEL_MEM(parmesh,locStats,PMMG_locateStats,"locStats");
#endif

  return ier;
}

/**
 * \file grpsplit_pmmg.c
 * \brief Split groups into sub groups.
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \author Nikos Pattakos (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */
#include "parmmg.h"
#include "metis_pmmg.h"

/**
 * \param nelem number of elements in the initial group
 * \param target_mesh_size wanted number of elements per group
 *
 * \return the needed number of groups
 *
 *  Compute the needed number of groups to create groups of \a target_mesh_size
 *  elements from a group of nelem elements.
 *
 */
static int PMMG_howManyGroups ( const int nelem, const int target_mesh_size )
{
  int ngrp = nelem / target_mesh_size;

  if ( ngrp == 0 )
    return ( 1 );
  else if ( ngrp * target_mesh_size < nelem )
    return ( ngrp + 1 );

  return ( ngrp );
}


/**
 * \param to       mesh to copy xtetra item to
 * \param from     mesh to copy xtetra item from
 * \param location location of xtetra to copy in original mesh
 *
 * \return 1 if tetra is successfully appended
 *         0 if tetra is not successfully appended
 *
 *  Append a tetraedron to the mesh, increasing the array if the allocated
 *  array is not big enough.
 */
static int PMMG_xtetraAppend( MMG5_pMesh to, MMG5_pMesh from, int location )
{
  int newsize;

  /* Adjust the value of scale to reallocate memory more or less agressively */
  const float scale = 2.f;
  if ( (to->xt + 1) >= to->xtmax ) {
    newsize = MG_MAX(scale*to->xtmax,to->xtmax+1);
    PMMG_REALLOC(to, to->xtetra, newsize+1, to->xtmax+1,MMG5_xTetra,
                  "larger xtetra table",return 0);
    to->xtmax = newsize;
  }
  ++to->xt;
  memcpy( &to->xtetra[ to->xt ],
          &from->xtetra[ from->tetra[ location ].xt ], sizeof(MMG5_xTetra) );
  return 1;
}

/**
 * \param to       mesh to copy xpoint item to
 * \param from     mesh to copy xpoint item from
 * \param location location of xpoint to copy in original mesh
 * \param point    xpoint's number in tetrahedron
 *
 * \return 1 if tetra is successfully appended
 *         0 if tetra is not successfully appended
 *
 *  Append a point in the point list, increasing the array if the allocated
 *  array is not big enough.
 */
static int PMMG_xpointAppend( MMG5_pMesh to, MMG5_pMesh from, int tetrahedron, int point )
{
  int newsize;

  // Adjust the value of scale to reallocate memory more or less agressively
  const float scale = 2.f;
  if ( (to->xp + 1) >= to->xpmax + 1 ) {
    newsize = MG_MAX(scale*to->xpmax,to->xpmax+1);
    PMMG_RECALLOC(to, to->xpoint, newsize+1, to->xpmax+1,MMG5_xPoint,
                  "larger xpoint table",return 0);
    to->xpmax = newsize;
  }
  ++to->xp;
  memcpy( &to->xpoint[ to->xp ],
          &from->xpoint[ from->point[ from->tetra[tetrahedron].v[point] ].xp ],
          sizeof(MMG5_xPoint) );
  return 1;
}
/**
 * \param parmesh parmmg struct pointer
 * \param grp     working group pointer
 * \param max     pointer toward the size of the node2int_node_comm arrays
 * \param idx1    index1 value
 * \param idx2    index2 value
 *
 * \return 1 if tetra is successfully appended
 *         0 if tetra is not successfully appended
 *
 *  Append new values in  group's internal communitor, resizing the buffers if
 *  required
 */
static int PMMG_n2incAppend( PMMG_pParMesh parmesh, PMMG_pGrp grp, int *max,
                             int idx1, int idx2 )
{
  int newsize;

  assert( (max != 0) && "null pointer passed" );
  // Adjust the value of scale to reallocate memory more or less agressively
  const float scale = 2.f;
  if ( (grp->nitem_int_node_comm + 1) >= *max ) {
    newsize = MG_MAX(scale * (*max),(*max)+1);
    PMMG_RECALLOC(parmesh, grp->node2int_node_comm_index1, newsize, *max, int,
                  "increasing node2int_node_comm_index1",return 0);
    PMMG_RECALLOC(parmesh, grp->node2int_node_comm_index2, newsize, *max, int,
                  "increasing node2int_node_comm_index2",return 0);
    *max  = newsize;
  }
  grp->node2int_node_comm_index1[ grp->nitem_int_node_comm ] = idx1;
  grp->node2int_node_comm_index2[ grp->nitem_int_node_comm ] = idx2;
  ++grp->nitem_int_node_comm;
  return 1;
}

/**
 * \param parmesh pointer toward the parmmg structure
 * \param grp     pointer toward the working group
 * \param max     pointer toward the size of the node2int_face_comm arrays
 * \param idx1    \f$ 12\times iel + 4\times ifac + iploc \f$ with \a iel the
 * index of the element to which the face belong, \a ifac the local index of the
 * face in \a iel and \a iploc as starting point in \a ifac (to be able to build
 * the node communicators from the face ones)
 * \param idx2    position of the face in the internal face communicator.
 *
 * \return 1 if tetra is successfully appended
 *         0 if tetra is not successfully appended
 *
 *  Append new values in the face internal communicator, resizing the buffers if
 *  required.
 */
static int PMMG_f2ifcAppend( PMMG_pParMesh parmesh, PMMG_pGrp grp, int *max,
                        int idx1, int idx2 )
{
  int newsize;

  assert( (max != 0) && "null pointer passed" );

  // Adjust the value of scale to reallocate memory more or less agressively
  const float scale = 2.f;

  if ( (grp->nitem_int_face_comm + 1) >= *max ) {
    newsize = MG_MAX(scale * (*max),(*max)+1);
    PMMG_RECALLOC(parmesh, grp->face2int_face_comm_index1,newsize, *max, int,
                  "increasing face2int_face_comm_index1",return 0);
    PMMG_RECALLOC(parmesh, grp->face2int_face_comm_index2,newsize, *max, int,
                  "increasing face2int_face_comm_index2",return 0);
    *max  = newsize;
  }

  grp->face2int_face_comm_index1[ grp->nitem_int_face_comm ] = idx1;
  grp->face2int_face_comm_index2[ grp->nitem_int_face_comm ] = idx2;
  ++grp->nitem_int_face_comm;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param np number of vertices.
 * \param ne number of tetrahedra.
 * \param nprism number of prisms.
 * \param xp number of boundary point
 * \param xt number of boundary tetra
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Check the input mesh size and assign their values to the mesh.
 *
 */
int PMMG_grpSplit_setMeshSize_initData(MMG5_pMesh mesh, int np, int ne,
                                       int nprism, int xp, int xt ) {

  if ( ( (mesh->info.imprim > PMMG_VERB_DETQUAL) || mesh->info.ddebug ) &&
       ( mesh->point || mesh->xpoint || mesh->tetra || mesh->xtetra) )
    fprintf(stderr,"\n  ## Warning: %s: old mesh deletion.\n",__func__);

  if ( !np ) {
    fprintf(stderr,"  ** MISSING DATA:\n");
    fprintf(stderr,"     Your mesh must contains at least points.\n");
    return(0);
  }
  if ( !ne && (mesh->info.imprim > PMMG_VERB_DETQUAL || mesh->info.ddebug) ) {
    fprintf(stderr,"  ** WARNING:\n");
    fprintf(stderr,"     Your mesh don't contains tetrahedra.\n");
  }

  if ( mesh->point )
    MMG5_DEL_MEM(mesh,mesh->point);
  if ( mesh->tetra )
    MMG5_DEL_MEM(mesh,mesh->tetra);
  if ( mesh->prism )
    MMG5_DEL_MEM(mesh,mesh->prism);
  if ( mesh->tria )
    MMG5_DEL_MEM(mesh,mesh->tria);
  if ( mesh->quadra )
    MMG5_DEL_MEM(mesh,mesh->quadra);
  if ( mesh->edge )
    MMG5_DEL_MEM(mesh,mesh->edge);

  mesh->np  = np;
  mesh->ne  = ne;
  mesh->xp  = xp;
  mesh->xt  = xt;
  mesh->nprism = nprism;

  mesh->npi = mesh->np;
  mesh->nei = mesh->ne;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Allocation of the array fields of the mesh for the given xpmax, npmax,xtmax,
 * nemax.
 *
 */
int PMMG_grpSplit_setMeshSize_alloc( MMG5_pMesh mesh ) {

  PMMG_CALLOC(mesh,mesh->point,mesh->npmax+1,MMG5_Point,
              "vertices array", return 0);

  PMMG_CALLOC(mesh,mesh->xpoint,mesh->xpmax+1,MMG5_xPoint,
              "boundary vertices array", return 0);

  PMMG_CALLOC(mesh,mesh->tetra,mesh->nemax+1,MMG5_Tetra,
              "tetra array", return 0);

  PMMG_CALLOC(mesh,mesh->xtetra,mesh->xtmax+1,MMG5_xTetra,
              "boundary tetra array", return 0);

  if ( mesh->nprism ) {
    PMMG_CALLOC(mesh,mesh->prism,mesh->nprism+1,MMG5_Prism,
                "prisms array", return 0);
  }
  if ( mesh->xpr ) {
    PMMG_CALLOC(mesh,mesh->xprism,mesh->xpr+1,MMG5_xPrism,
                "boundary prisms array", return 0);
  }

  return ( PMMG_link_mesh( mesh ) );
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param np number of vertices.
 * \param ne number of tetrahedra.
 * \param nprism number of prisms.
 * \param xp number of boundary points.
 * \param xt number of boundary tetra.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Check the input mesh size and assign their values to the mesh.
 *
 */
int PMMG_grpSplit_setMeshSize(MMG5_pMesh mesh,int np,int ne,
                              int nprism,int xp,int xt ) {

  /* Check input data and set mesh->ne/na/np/nt to the suitable values */
  if ( !PMMG_grpSplit_setMeshSize_initData(mesh,np,ne,nprism,xp,xt) )
    return 0;

  mesh->npmax  = mesh->np;
  mesh->xpmax  = mesh->xp;
  mesh->nemax  = mesh->ne;
  mesh->xtmax  = mesh->xt;

  /* Mesh allocation and linkage */
  if ( !PMMG_grpSplit_setMeshSize_alloc( mesh ) ) return 0;

  return(1);

}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param igrp index of the group to create
 *
 * \return 0 if fail, 1 if success
 *
 * Creation of a new group for the background mesh (without communication
 * structures (info.inputMet == 1 if a metrics is provided by the user).
 *
 */
int PMMG_oldGrps_newGroup( PMMG_pParMesh parmesh,int igrp ) {
  MMG5_pMesh const meshOld= parmesh->listgrp[igrp].mesh;
  MMG5_pSol  const metOld = parmesh->listgrp[igrp].met;
  PMMG_pGrp        grp;
  MMG5_pMesh       mesh;
  MMG5_pSol        met;
  size_t           oldMemMax,memAv;

  grp = &parmesh->old_listgrp[igrp];
  grp->mesh = NULL;
  grp->met  = NULL;

  MMG3D_Init_mesh( MMG5_ARG_start,
                   MMG5_ARG_ppMesh, &grp->mesh,
                   MMG5_ARG_ppMet, &grp->met,
                   MMG5_ARG_end );

  mesh = grp->mesh;
  met  = grp->met;

  /* Give all the available memory to the mesh */
  oldMemMax = parmesh->memCur;
  memAv     = parmesh->memMax-oldMemMax;
  PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh,memAv,oldMemMax);

  /* Copy the mesh filenames */
  if ( !MMG5_Set_inputMeshName(  mesh,meshOld->namein) )      return 0;
  if ( !MMG5_Set_inputSolName(   mesh,met,metOld->namein ) )  return 0;
  if ( !MMG5_Set_outputMeshName( mesh,meshOld->nameout ) )    return 0;
  if ( !MMG5_Set_outputSolName(  mesh,met,metOld->nameout ) ) return 0;

  /* Set sizes and allocate new mesh */
  if ( !PMMG_grpSplit_setMeshSize( mesh,meshOld->np,meshOld->ne,0,0,0) )
    return 0;

  PMMG_CALLOC(mesh,mesh->adja,4*mesh->nemax+5,int,"adjacency table",return 0);

  /* Set metrics size */
  if ( meshOld->info.inputMet == 1 )
    if ( !MMG3D_Set_solSize(mesh,met,MMG5_Vertex,meshOld->np,metOld->type) )
      return 0;

  /* Copy the info structure of the initial mesh: it contains the remeshing
   * options */
  memcpy(&(mesh->info),&(meshOld->info),sizeof(MMG5_Info) );


  /* Give the available memory to the parmesh */
  PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,mesh,memAv,oldMemMax);

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param group pointer toward the new group to create
 * \param igrp index of the old group which is splitted
 * \param memAv available mem for the mesh allocation
 * \param ne number of elements in the new group mesh
 * \param f2ifc_max maximum number of elements in the face2int_face_comm arrays
 * \param n2inc_max maximum number of elements in the node2int_node_comm arrays
 *
 * \return 0 if fail, 1 if success
 *
 * Creation of the new group \a grp: allocation and initialization of the mesh
 * and communicator structures. Set the f2ifc_max variable at the size at which
 * we allocate the face2int_face_comm arrays and the n2inc_max variable at the
 * size at which the node2int_node_comm arrays are allocated.
 *
 */
static int
PMMG_splitGrps_newGroup( PMMG_pParMesh parmesh,PMMG_pGrp grp,int igrp,
                         size_t *memAv,int ne,int *f2ifc_max,int *n2inc_max ) {
  PMMG_pGrp  const grpOld = &parmesh->listgrp[igrp];
  MMG5_pMesh const meshOld= parmesh->listgrp[igrp].mesh;
  MMG5_pMesh       mesh;
  size_t           oldMemMax;

  grp->mesh = NULL;
  grp->met  = NULL;
  grp->disp = NULL;

  MMG3D_Init_mesh( MMG5_ARG_start, MMG5_ARG_ppMesh, &grp->mesh,
                   MMG5_ARG_ppMet, &grp->met, MMG5_ARG_end );

  mesh      = grp->mesh;

  /* Give all the available memory to the mesh */
  mesh->memMax = *memAv;

  /* Copy the mesh filenames */
  if ( !MMG5_Set_inputMeshName( mesh,meshOld->namein) )                return 0;
  if ( !MMG5_Set_inputSolName(  mesh,grp->met,grpOld->met->namein ) )  return 0;
  if ( !MMG5_Set_outputMeshName(mesh, meshOld->nameout ) )             return 0;
  if ( !MMG5_Set_outputSolName( mesh,grp->met,grpOld->met->nameout ) ) return 0;

  /* Uses the Euler-poincare formulae to estimate the number of entities (np =
   * ne/6, nt=ne/3 */
  if ( !PMMG_grpSplit_setMeshSize( mesh,MG_MAX(ne/6,4),ne,0,MG_MAX(ne/6,3),MG_MAX(ne/3,1)) ) return 0;

  PMMG_CALLOC(mesh,mesh->adja,4*mesh->nemax+5,int,"adjacency table",return 0);

  grp->mesh->np = 0;
  grp->mesh->npi = 0;

  if ( grpOld->met->m ) {
    if ( grpOld->met->size == 1 )
      grp->met->type = MMG5_Scalar;
    else if ( grpOld->met->size == 6 )
      grp->met->type = MMG5_Tensor;

    /** If we have an initial metric, force the metric allocation (even if for
     * now, we don't know the number of point that will be stored in it) */
    int allocMetric = 1;
    if ( !MMG3D_Set_solSize(grp->mesh,grp->met,MMG5_Vertex,allocMetric,grp->met->type) )
      return 0;
  }

  /* Copy the info structure of the initial mesh: it contains the remeshing
   * options */
  if ( !PMMG_copy_mmgInfo ( &meshOld->info,&grp->mesh->info ) ) return 0;


  /* Give the available memory to the parmesh */
  grp->mesh->memMax = grp->mesh->memCur;
  if ( *memAv < grp->mesh->memMax ) {
      fprintf(stderr,"\n  ## Error: %s: not enough memory to allocate a new"
              " mesh.\n",__func__);
      return 0;
  }

  *memAv           -= grp->mesh->memMax;

  oldMemMax         = parmesh->memMax;
  parmesh->memMax  += *memAv;

  *n2inc_max = ne/3;
  assert( (grp->nitem_int_node_comm == 0 ) && "non empty comm" );
  PMMG_CALLOC(parmesh,grp->node2int_node_comm_index1,*n2inc_max,int,
              "subgroup internal1 communicator ",return 0);
  PMMG_CALLOC(parmesh,grp->node2int_node_comm_index2,*n2inc_max,int,
              "subgroup internal2 communicator ",return 0);

  *f2ifc_max = mesh->xtmax;
  PMMG_CALLOC(parmesh,grp->face2int_face_comm_index1,*f2ifc_max,int,
              "face2int_face_comm_index1 communicator",return 0);
  PMMG_CALLOC(parmesh,grp->face2int_face_comm_index2,*f2ifc_max,int,
              "face2int_face_comm_index2 communicator",return 0);

  /* Update the available memory */
  parmesh->memMax = parmesh->memCur;
  assert ( parmesh->memMax >= oldMemMax );
  if ( *memAv < (parmesh->memMax - oldMemMax) ) {
      fprintf(stderr,"\n  ## Error: %s: not enough memory.\n",__func__);
      return 0;
  }

  *memAv -= (parmesh->memMax - oldMemMax);

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param igrp index of the group to fill
 *
 * Fill the background mesh with the current mesh on group igrp
 * (info.inputMet == 1 if a  metrics is provided by the user).
 *
 */
int PMMG_oldGrps_fillGroup( PMMG_pParMesh parmesh,int igrp ) {

  MMG5_pMesh const meshOld= parmesh->listgrp[igrp].mesh;
  MMG5_pSol  const metOld = parmesh->listgrp[igrp].met;
  MMG5_pMesh       mesh;
  MMG5_pSol        met;
  MMG5_pTetra      pt,ptCur;
  MMG5_pPoint      ppt,pptCur;
  int              *adja,*oldAdja,ie,ip;

  mesh = parmesh->old_listgrp[igrp].mesh;
  met  = parmesh->old_listgrp[igrp].met;

  assert( mesh->ne == meshOld->ne );
  assert( mesh->np == meshOld->np );

  /* Loop on tetras */
  for ( ie = 1; ie < meshOld->ne+1; ++ie ) {
    pt = &meshOld->tetra[ie];
    ptCur = &mesh->tetra[ie];
 
    if ( !MG_EOK(pt) ) continue;

    /* Copy tetra */
    memcpy( ptCur, pt, sizeof(MMG5_Tetra) );

    /* Copy element's adjacency */
    if( meshOld->adja ) {
      adja    =    &mesh->adja[ 4*( ie-1 )+1 ];
      oldAdja = &meshOld->adja[ 4*( ie-1 )+1 ];
      memcpy( adja, oldAdja, 4*sizeof(int) );
    }

    /* Skip xtetra */
    ptCur->xt = 0;

  }

  /* Loop on points */
  for ( ip = 1; ip < meshOld->np+1; ++ip ) {
    ppt = &meshOld->point[ip];
    pptCur = &mesh->point[ip];
 
    if ( !MG_VOK(ppt) ) continue;

    /* Copy point */
    memcpy( pptCur, ppt, sizeof(MMG5_Point) );

    /* Copy metrics */
    if ( mesh->info.inputMet == 1 )
      memcpy( &met->m[ ip*met->size ], &metOld->m[ip*met->size], met->size*sizeof(double) );

    /* Skip xpoint */
    pptCur->xp = 0;

  }
  
  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param group pointer toward the new group to fill
 * \param mesh pointer toward the new mesh to fill
 * \param meshOld pointer toward the old mesh
 * \param tetraCur pointer to the current tetra
 * \param pt pointer to the old tetra
 * \param tet index of the old tetra
 * \param grpId intex of the current group
 * \param n2ifc_max maximum number of nodes in the node2int_node_comm arrays
 * \param part partition array for the old tetrahedra
 * \param memAv pointer to the available memory
 * \param oldMemMax pointer to the old max memory
 *
 * \return 0 if fail, 1 if success.
 *
 * Fill the node communicator with new interface nodes.
 *
 */
static int PMMG_splitGrps_updateNodeCommNew( PMMG_pParMesh parmesh,PMMG_pGrp grp,
    MMG5_pMesh mesh, MMG5_pMesh meshOld,MMG5_pTetra tetraCur,MMG5_pTetra pt,
    int tet,int grpId,int *n2inc_max,int *part,size_t *memAv,size_t *oldMemMax ) {
  MMG5_pPoint ppt;
  int *adja,fac,poi,adjidx;

  adja = &meshOld->adja[ 4 * ( tet - 1 ) + 1 ];
  for ( fac = 0; fac < 4; ++fac ) {
    adjidx = adja[ fac ] / 4;

    if ( adjidx && grpId != part[ adjidx - 1 ] ) {
      for ( poi = 0; poi < 3; ++poi ) {
        ppt = & mesh->point[ tetraCur->v[ MMG5_idir[fac][poi] ] ];
        if ( ppt->tmp == -1 ) {
          if (   !PMMG_n2incAppend( parmesh, grp, n2inc_max,
                                   tetraCur->v[ MMG5_idir[fac][poi] ],
                                   parmesh->int_node_comm->nitem + 1 ) ) {
            return 0;
          }

          meshOld->point[ pt->v[ MMG5_idir[fac][poi]  ] ].tmp =
            parmesh->int_node_comm->nitem + 1;
          ppt->tmp = parmesh->int_node_comm->nitem + 1;

          ++parmesh->int_node_comm->nitem;
        }
      }
    }
  }
  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param group pointer toward the new group to fill
 * \param mesh pointer toward the new mesh to fill
 * \param ppt pointer to the current point
 * \param it index of the current point
 * \param fac intex of the current face
 * \param n2ifc_max maximum number of nodes in the node2int_node_comm arrays
 * \param memAv pointer to the available memory
 * \param oldMemMax pointer to the old max memory
 * \param pos pointer to the tetra+face index
 *
 * \return 0 if fail, 1 if success.
 *
 * Fill the node communicator if the node was already parallel.
 *
 */
static int PMMG_splitGrps_updateNodeCommOld( PMMG_pParMesh parmesh,PMMG_pGrp grp,
    MMG5_pMesh mesh,MMG5_pPoint ppt,int ip,int *n2inc_max,size_t *memAv,size_t *oldMemMax ) {

  /* Give the available memory to the parmesh */
  PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,mesh,*memAv,*oldMemMax);

  /* Add point in subgroup's communicator if it already was in group's
     communicator */
  if ( ppt->tmp != PMMG_UNSET ) {
    if (  !PMMG_n2incAppend( parmesh, grp, n2inc_max,
                             ip, ppt->tmp ) ) {
      return 0;
    }
    ++parmesh->int_node_comm->nitem;
  }

  /* Give back the memory to the mesh */
  PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh,*memAv,*oldMemMax);

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param group pointer toward the new group to fill
 * \param mesh pointer toward the new mesh to fill
 * \param meshOld pointer toward the old mesh
 * \param fac intex of the current face
 * \param adjidx index of the old tetra
 * \param vidx face of the old tetra
,* \param posInIntFaceComm position of each tetra face in the internal face
 * \param iplocFaceComm starting index to list the vertices of the faces in the
 * \param f2ifc_max maximum number of elements in the face2int_face_comm arrays
 * face2int_face arrays (to be able to build the node communicators from the
 * face ones).
 * communicator (-1 if not in the internal face comm)
 * \param tetPerGrp number of tetra in the group
 * \param memAv pointer to the available memory
 * \param oldMemMax pointer to the old max memory
 * \param pos tetra+face index
 *
 * \return 0 if fail, 1 if success.
 *
 * Fill the face communicator with new parallel faces.
 *
 */
static int PMMG_splitGrps_updateFaceCommNew( PMMG_pParMesh parmesh,
    PMMG_pGrp grp, MMG5_pMesh mesh,MMG5_pMesh meshOld,MMG5_pTetra pt,int fac,int adjidx,int vidx,
    int *posInIntFaceComm, int *iplocFaceComm,int *f2ifc_max,int tetPerGrp,
    size_t *memAv,size_t *oldMemMax,int pos ) {
  MMG5_pTetra ptadj;
  int ip,iploc,iplocadj;

  /* Give the available memory to the parmesh */
  PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,mesh,*memAv,*oldMemMax);

  /** Build the internal communicator for the boundary faces */
  /* 1) Check if this face has already a position in the internal face
   * communicator, and if not, increment the internal face comm and
   * store the face position in posInIntFaceComm */
  if ( posInIntFaceComm[4*(adjidx-1)+1+vidx]<0 ) {

    posInIntFaceComm[pos]                 = parmesh->int_face_comm->nitem;
    posInIntFaceComm[4*(adjidx-1)+1+vidx] = parmesh->int_face_comm->nitem;

    /* Find a common starting point inside the face for both tetra */
    iploc = 0; // We impose the starting point in tet
    ip    = pt->v[MMG5_idir[fac][iploc]];

    ptadj = &meshOld->tetra[adjidx];
    for ( iplocadj=0; iplocadj < 3; ++iplocadj )
      if ( ptadj->v[MMG5_idir[vidx][iplocadj]] == ip ) break;
    assert ( iplocadj < 3 );

    iplocFaceComm[pos]                 = iploc;
    iplocFaceComm[4*(adjidx-1)+1+vidx] = iplocadj;

    /* 2) Add the face in the list of interface faces of the group */
    if ( !PMMG_f2ifcAppend( parmesh, grp, f2ifc_max,12*tetPerGrp+3*fac+iploc,
                           posInIntFaceComm[pos] ) ) {

      return 0;
    }
    ++parmesh->int_face_comm->nitem;
  }
  else {
    assert ( posInIntFaceComm[pos] >=0 );
    assert ( iplocFaceComm   [pos] >=0 );

    /* 2) Add the face in the list of interface faces of the group */
    iploc = iplocFaceComm[pos];
    if ( !PMMG_f2ifcAppend( parmesh, grp, f2ifc_max,12*tetPerGrp+3*fac+iploc,
                           posInIntFaceComm[pos] ) ) {

      return 0;
    }
  }

  PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh,*memAv,*oldMemMax);

  return 1;

}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param group pointer toward the new group to fill
 * \param mesh pointer toward the new mesh to fill
 * \param adja pointer to the current tetra adjacency
 * \param tet index of the current tetra
 * \param fac intex of the current face
,* \param posInIntFaceComm position of each tetra face in the internal face
 * \param iplocFaceComm starting index to list the vertices of the faces in the
 * \param f2ifc_max maximum number of elements in the face2int_face_comm arrays
 * face2int_face arrays (to be able to build the node communicators from the
 * face ones).
 * communicator (-1 if not in the internal face comm)
 * \param tetPerGrp number of tetra in the group
 * \param memAv pointer to the available memory
 * \param oldMemMax pointer to the old max memory
 * \param pos pointer to the tetra+face index
 *
 * \return 0 if fail, 1 if success, -1 if you can continue to the next tetra
 *
 * Fill the face communicator if the face was already parallel.
 *
 */
static int PMMG_splitGrps_updateFaceCommOld( PMMG_pParMesh parmesh,
    PMMG_pGrp grp, MMG5_pMesh mesh,int *adja,int tet,int fac,
    int *posInIntFaceComm, int *iplocFaceComm,int *f2ifc_max,int tetPerGrp,
    size_t *memAv,size_t *oldMemMax,int *pos ) {
  int iploc;
 
  *pos   = 4*(tet-1)+1+fac;

  /* Give the available memory to the parmesh */
  PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,mesh,*memAv,*oldMemMax);

  if ( adja[ fac ] == 0 ) {
    /** Build the internal communicator for the parallel faces */
    /* 1) Check if this face has a position in the internal face
     * communicator. If not, the face is not parallel and we have nothing
     * to do */
    if ( posInIntFaceComm[*pos]<0 ) {
      /* Give the available memory to the mesh */
      PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh,*memAv,*oldMemMax);
      return -1;
    }

    /* 2) Add the point in the list of interface faces of the group */
    iploc = iplocFaceComm[*pos];
    assert ( iploc >=0 );

    if ( !PMMG_f2ifcAppend( parmesh, grp, f2ifc_max,12*tetPerGrp+3*fac+iploc,
                           posInIntFaceComm[*pos] ) ) {
      return 0;
    }
    /* Give the available memory to the mesh */
    PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh,*memAv,*oldMemMax);
    return -1;
  }

  /* Give the available memory to the mesh */
  PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh,*memAv,*oldMemMax);

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param group pointer toward the new group to fill
 * \param grpIdOld index of the group that is splitted in the old list of groups
 * \param grpId index of the group that we create in the list of groups
 * \param ne number of elements in the new group mesh
 * \param np pointer toward number of points in the new group mesh
 * \param f2ifc_max maximum number of elements in the face2int_face_comm arrays
 * \param n2inc_max maximum number of elements in the node2int_node_comm arrays
 * \param part metis partition
 * \param posInIntFaceComm position of each tetra face in the internal face
 * \param iplocFaceComm starting index to list the vertices of the faces in the
 * face2int_face arrays (to be able to build the node communicators from the
 * face ones).
 * communicator (-1 if not in the internal face comm)
 *
 * \return 0 if fail, 1 if success
 *
 * Fill the mesh and communicators of the new group \a grp.
 *
 */
static int
PMMG_splitGrps_fillGroup( PMMG_pParMesh parmesh,PMMG_pGrp grp,int grpIdOld,int grpId,int ne,
                          int *np,int *f2ifc_max,int *n2inc_max,idx_t *part,
                          int* posInIntFaceComm,int* iplocFaceComm,size_t *memAv ) {

  PMMG_pGrp  const grpOld = &parmesh->listgrp[grpIdOld];
  MMG5_pMesh const meshOld= parmesh->listgrp[grpIdOld].mesh;
  MMG5_pMesh       mesh;
  MMG5_pSol        met;
  MMG5_pTetra      pt,ptadj,tetraCur;
  MMG5_pxTetra     pxt;
  MMG5_pPoint      ppt;
  size_t           oldMemMax;
  int              *adja,adjidx,vidx,fac,pos,ip,iploc,iplocadj;
  int              ie,tetPerGrp,tet,poi,j,newsize;
  int              ier;

  mesh = grp->mesh;
  met  = grp->met;

  /** Give the available memory to the mesh */
  oldMemMax = mesh->memMax;
  mesh->memMax += *memAv;

  /** Reinitialize to the value that n2i_n_c arrays are initially allocated
      Otherwise grp #1,2,etc will incorrectly use the values that the previous
      grps assigned to n2inc_max */
  (*n2inc_max) = ne/3;

  /* Loop over tetras and choose the ones to add in the submesh being constructed */
  tetPerGrp = 0;
  *np       = 0;
  for ( tet = 1; tet < meshOld->ne + 1; ++tet ) {
    pt = &meshOld->tetra[tet];

    if ( !MG_EOK(pt) ) continue;

    /* Skip elements that do not belong in the group processed in this iteration */
    if ( grpId != part[ tet - 1 ] )
      continue;

    /** MMG3D_Tetra.flag is used to update adjacency vector:
       if the tetra belongs to the group we store the local tetrahedron id
       in the tet.flag */
    ++tetPerGrp;
    assert( ( tetPerGrp <= mesh->nemax ) && "overflowing tetra array?" );
    tetraCur = mesh->tetra + tetPerGrp;
    pt->flag = tetPerGrp;

    /* add tetrahedron to subgroup (copy from original group) */
    memcpy( tetraCur, pt, sizeof(MMG5_Tetra) );
    tetraCur->base = 0;
    tetraCur->flag = tet;

    /* xTetra: this element was already an xtetra (in meshOld) */
    if ( tetraCur->xt != 0 ) {
      if ( !PMMG_xtetraAppend( mesh, meshOld, tet ) )
        return 0;
      tetraCur->xt = mesh->xt;
    }

    /* Add tetrahedron vertices in points struct and
       adjust tetrahedron vertices indices */
    for ( poi = 0; poi < 4 ; ++poi ) {
      if ( meshOld->point[ pt->v[poi] ].s != grpId ) {
        /* 1st time that this point is seen in this subgroup
           Add point in subgroup point array */
        ++(*np);

        if ( *np > mesh->npmax ) {
          newsize = MG_MAX((int)((1+mesh->gap)*mesh->npmax),mesh->npmax+1);
          PMMG_RECALLOC(mesh,mesh->point,newsize+1,mesh->npmax+1,MMG5_Point,
                        "point array",return 0);
          mesh->npmax = newsize;
          mesh->npnil = (*np)+1;;
          for (j=mesh->npnil; j<mesh->npmax-1; j++)
            mesh->point[j].tmp  = j+1;

          if ( met->m ) {
            PMMG_REALLOC(mesh,met->m,met->size*(mesh->npmax+1),
                         met->size*(met->npmax+1),double,
                         "metric array",return 0);
          }
          met->npmax = mesh->npmax;
          assert ( *np<=mesh->npmax );
        }
        memcpy( mesh->point+(*np),&meshOld->point[pt->v[poi]],
                sizeof(MMG5_Point) );
        if ( met->m ) {
          memcpy( &met->m[ (*np) * met->size ],
                  &grpOld->met->m[pt->v[poi] * met->size],
                  met->size * sizeof( double ) );
        }

        /* Update tetra vertex index */
        tetraCur->v[poi] = (*np);

        /* Store the point group */
        meshOld->point[ pt->v[poi] ].s = grpId;

        /* Store the point id */
        meshOld->point[ pt->v[poi] ].flag = (*np);

        ppt = &mesh->point[*np];
        /* xPoints: this was already a boundary point */
        if ( mesh->point[*np].xp != 0 ) {
          if ( !PMMG_xpointAppend(mesh,meshOld,tet,poi) ) {
            return 0;
          }
          mesh->point[*np].xp = mesh->xp;
        }

        if( !PMMG_splitGrps_updateNodeCommOld( parmesh,grp,mesh,ppt,*np,
              n2inc_max,memAv,&oldMemMax ) ) return 0;

      } else {
        // point is already included in this subgroup, update current tetra
        // vertex reference
        tetraCur->v[poi] = meshOld->point[ pt->v[poi] ].flag;
      }
    }


    /* Copy element's vertices adjacency from old mesh and update them to the
     * new mesh values */
    assert( ((4*(tetPerGrp-1)+4)<(4*(mesh->ne-1)+5)) && "adja overflow" );
    adja = &mesh->adja[ 4 * ( tetPerGrp - 1 ) + 1 ];

    assert( (4 *(tet-1)+1+3) < (4*(meshOld->ne-1)+1+5) && "mesh->adja overflow" );
    memcpy( adja, &meshOld->adja[ 4 * ( tet - 1 ) + 1 ], 4 * sizeof(int) );

    /* Update element's adjaceny to elements in the new mesh */
    for ( fac = 0; fac < 4; ++fac ) {

      ier = PMMG_splitGrps_updateFaceCommOld( parmesh,grp,mesh,adja,tet,fac,
          posInIntFaceComm,iplocFaceComm,f2ifc_max,tetPerGrp,memAv,&oldMemMax,&pos );
      if( ier == 0 )  return 0;
      if( ier == -1 ) continue;

      adjidx = adja[ fac ] / 4;
      vidx   = adja[ fac ] % 4;

      /* new boundary face: set to 0, add xtetra and set tags */
      assert( ((adjidx - 1) < meshOld->ne ) && "part[adjaidx] would overflow" );
      if ( part[ adjidx - 1 ] != grpId ) {
        adja[ fac ] = 0;

        /* creation of the interface faces : ref 0 and tag MG_PARBDY */
        if( !mesh->tetra[tetPerGrp].xt ) {
          if ( !PMMG_xtetraAppend( mesh, meshOld, tet ) ) {
            return 0;
          }
          tetraCur->xt = mesh->xt;
        }
        pxt = &mesh->xtetra[tetraCur->xt];
        pxt->ref[fac] = 0;
        /* If already boundary, make it recognizable as a "true" boundary */
        if( pxt->ftag[fac] & MG_BDY ) pxt->ftag[fac] |= MG_PARBDYBDY; 
        pxt->ftag[fac] |= (MG_PARBDY + MG_BDY + MG_REQ + MG_NOSURF);

        if( !PMMG_splitGrps_updateFaceCommNew( parmesh,grp,mesh,meshOld,pt,fac,adjidx,vidx,posInIntFaceComm,iplocFaceComm,f2ifc_max,tetPerGrp,memAv,&oldMemMax,pos ) ) return 0;

        for ( j=0; j<3; ++j ) {
          /* Update the face and face vertices tags */
          pxt->tag[MMG5_iarf[fac][j]] |= (MG_PARBDY + MG_BDY + MG_REQ + MG_NOSURF);
          ppt = &mesh->point[tetraCur->v[MMG5_idir[fac][j]]];
          ppt->tag |= (MG_PARBDY + MG_BDY + MG_REQ + MG_NOSURF);

          /** Add an xPoint if needed */
// TO REMOVE WHEN MMG WILL BE READY
          if ( !ppt->xp ) {
            if ( (mesh->xp+1) > mesh->xpmax ) {
              /* realloc of xtetras table */
              newsize = MG_MAX((int)((1+mesh->gap)*mesh->xpmax),mesh->xpmax+1);
              PMMG_RECALLOC(mesh,mesh->xpoint,newsize+1,mesh->xpmax+1,MMG5_xPoint,
                            "larger xpoint ",return 0);
              mesh->xpmax = newsize;
            }
            ++mesh->xp;
            ppt->xp = mesh->xp;
          }
// TO REMOVE WHEN MMG WILL BE READY
        }

        // if the adjacent number is already processed
      } else if ( adjidx < tet ) {
        adja[ fac ] =  4 * meshOld->tetra[ adjidx ].flag  + vidx;
        // Now that we know the local gr_idx for both tetra we should also
        // update the adjacency entry for the adjacent face
        assert(     (4 * (meshOld->tetra[adjidx].flag-1) + 1 + vidx )
                    < (4 * (mesh->ne -1 ) + 1 + 4)
                    && "adja overflow" );
        mesh->adja[4*(meshOld->tetra[adjidx].flag-1)+1+vidx] =
          4*tetPerGrp+fac;
      }
    }

  }
  assert( (mesh->ne == tetPerGrp) && "Error in the tetra count" );

  /* Give the available memory to the parmesh */
  PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,mesh,*memAv,oldMemMax);

  for ( ie = 1; ie <= mesh->ne; ie++ ) {
    /* Get tetra in the new mesh */
    tetraCur = &mesh->tetra[ie];

    /* Get tetra in the old mesh */
    tet = tetraCur->flag;
    pt = &meshOld->tetra[tet];

    if( !PMMG_splitGrps_updateNodeCommNew( parmesh,grp,mesh,meshOld,tetraCur,pt,
          tet,grpId,n2inc_max,part,memAv,&oldMemMax ) ) return 0;
  }

  PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh,*memAv,oldMemMax);

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param meshOld pointer toward the old mesh
 * \param part partition array for the tetra in the old mesh
 * \param tetPerGrp array of number of tetra for each new group
 *
 * Count the number of new tetra in each group, and store the new tetra ID in
 * the old tetra flag.
 *
 */
static void
PMMG_splitGrps_countTetPerGrp( PMMG_pParMesh parmesh,MMG5_pMesh meshOld,
                               idx_t *part,int *tetPerGrp) {
  MMG5_pTetra pt;
  int ie;

  for( ie = 1; ie <= meshOld->ne; ie++ ) {
    pt = &meshOld->tetra[ie];
    pt->flag = ++tetPerGrp[ part[ie-1] ];
  }
}

/**
 * \param mesh pointer toward an MMG5 mesh structure
 * \param met pointer toward an MMG5 metric structure
 * \param np number of points in the mesh
 *
 * \return 0 if fail, 1 if success
 *
 * Clean the mesh filled by the \a split_grps function to make it valid:
 *   - reallocate the mesh at it exact size
 *   - set the np/ne/npi/nei/npnil/nenil fields to suitables value and keep
 *   track of empty link
 *   - update the edge tags in all the xtetra of the edge shell
 *
 */
static inline
int PMMG_splitGrps_cleanMesh( MMG5_pMesh mesh,MMG5_pSol met,int np )
{
  MMG5_pTetra  pt;
  MMG5_pxTetra pxt;
  int          k,i;

  /* Mesh reallocation at the smallest possible size */
  PMMG_REALLOC(mesh,mesh->point,np+1,mesh->npmax+1,
               MMG5_Point,"fitted point table",return 0);
  mesh->npmax = np;
  mesh->npnil = 0;
  mesh->nenil = 0;

  PMMG_REALLOC(mesh,mesh->xpoint,mesh->xp+1,mesh->xpmax+1,
               MMG5_xPoint,"fitted xpoint table",return 0);
  mesh->xpmax = mesh->xp;
  PMMG_REALLOC(mesh,mesh->tetra,mesh->ne+1,mesh->nemax+1,
               MMG5_Tetra,"fitted tetra table",return 0);
  mesh->nemax = mesh->ne;
  PMMG_REALLOC(mesh,mesh->xtetra,mesh->xt+1,mesh->xtmax+1,
               MMG5_xTetra,"fitted xtetra table",return 0);
  mesh->xtmax = mesh->xt;
  if ( met->m )
    PMMG_REALLOC(mesh,met->m,met->size*(np+1),met->size*(met->npmax+1),
                 double,"fitted metric table",return 0);
  met->npmax = mesh->npmax;

  /* Set memMax to the smallest possible value */
  mesh->memMax = mesh->memCur;

  // Update the empty points' values as per the convention used in MMG3D
  mesh->np  = np;
  mesh->npi = np;

  if ( met->m ) {
    met->np  = np;
    met->npi = np;
  }

  /* Udate tags and refs of tetra edges (if we have 2 boundary tetra in the
   * shell of an edge, it is possible that one of the xtetra has set the edge
   * as MG_PARBDY. In this case, this tag must be reported in the second
   * xtetra) */
  for ( k = 1; k < mesh->ne +1; ++k ) {
    pt = &mesh->tetra[k];

    if ( !MG_EOK(pt) ) continue;
    if ( !pt->xt )     continue;

    pxt = &mesh->xtetra[pt->xt];

    for ( i=0; i<6; ++i ) {
      if ( !(pxt->tag[i] & MG_PARBDY ) )  continue;
      // Algiane: if this step is too long, try to hash the updated edges to not
      // update twice the same shell (PMMG_bdryUpdate function).
      MMG5_settag(mesh,k,i,pxt->tag[i],pxt->edg[i]);
    }
  }

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 if success
 *
 * Copy all groups from the current to the background list.
 *
 */
int PMMG_update_oldGrps( PMMG_pParMesh parmesh ) {
  int grpId;

  PMMG_listgrp_free(parmesh, &parmesh->old_listgrp, parmesh->nold_grp);

  /* Allocate list of subgroups struct and allocate memory */
  parmesh->nold_grp = parmesh->ngrp;
  PMMG_CALLOC(parmesh,parmesh->old_listgrp,parmesh->nold_grp,PMMG_Grp,
              "old group list ",return 0);

  /** Copy every group */
  for ( grpId = 0; grpId < parmesh->ngrp; ++grpId ) {

    /* New group initialisation */
    if ( !PMMG_oldGrps_newGroup( parmesh, grpId ) ) {
      fprintf(stderr,"\n  ## Error: %s: unable to initialize new background"
              " group (%d).\n",__func__,grpId);
      return 0;
    }

    /* Fill group */
    if ( !PMMG_oldGrps_fillGroup( parmesh, grpId ) ) {
      fprintf(stderr,"\n  ## Error: %s: unable to fill new background"
              "group (%d).\n",__func__,grpId);
      return 0;
    }
  }

  return 1;
}

int PMMG_split_eachGrp( PMMG_pParMesh parmesh,int grpIdOld,PMMG_pGrp grpsNew,idx_t ngrp,int *countPerGrp,idx_t *part ) {
  PMMG_pGrp grpOld,grpCur;
  MMG5_pMesh meshOld,meshCur;
  /** size of allocated node2int_node_comm_idx. when comm is ready trim to
   *  actual node2int_node_comm */
  int n2inc_max,f2ifc_max;
  n2inc_max = f2ifc_max = 0;
  size_t memAv,oldMemMax;
  int poiPerGrp = 0;
  int *posInIntFaceComm,*iplocFaceComm;
  int i, grpId, poi, fac, ie;
  int ret_val = 1;
 
  grpOld  = &parmesh->listgrp[grpIdOld];
  meshOld = grpOld->mesh;

  /* Use the posInIntFaceComm array to remember the position of the tetra faces
   * in the internal face communicator */
  posInIntFaceComm = NULL;
  iplocFaceComm    = NULL;
  PMMG_MALLOC(parmesh,posInIntFaceComm,4*meshOld->ne+1,int,
              "array of faces position in the internal face commmunicator ",
              ret_val = 0;goto fail_facePos);
  for ( i=0; i<=4*meshOld->ne; ++i )
    posInIntFaceComm[i] = PMMG_UNSET;

  PMMG_MALLOC(parmesh,iplocFaceComm,4*meshOld->ne+1,int,
              "starting vertices of the faces of face2int_face_comm_index1",
              ret_val = 0;goto fail_facePos);
  for ( i=0; i<=4*meshOld->ne; ++i )
    iplocFaceComm[i] = PMMG_UNSET;

  for ( i=0; i<grpOld->nitem_int_face_comm; ++i ) {
    ie  =  grpOld->face2int_face_comm_index1[i]/12;
    fac = (grpOld->face2int_face_comm_index1[i]%12)/3;
    posInIntFaceComm[4*(ie-1)+1+fac] = grpOld->face2int_face_comm_index2[i];
    iplocFaceComm[4*(ie-1)+1+fac] = (grpOld->face2int_face_comm_index1[i]%12)%3;
  }

  /* 
   * Use point[].tmp field to store index in internal communicator of
     vertices. specifically: place a copy of vertices' node2index2 position at
     point[].tmp field or -1 if they are not in the comm.
     Use point[].s field to store assigned point subgroup.
     Use point[].flag field to "remember" assigned local(in subgroup) numbering
     (no need to initialize it).
   */
  for ( poi = 1; poi < meshOld->np + 1; ++poi ) {
    meshOld->point[poi].tmp  = PMMG_UNSET;
    meshOld->point[poi].s    = PMMG_UNSET;
  }

  for ( i = 0; i < grpOld->nitem_int_node_comm; i++ )
    meshOld->point[ grpOld->node2int_node_comm_index1[ i ] ].tmp =
      grpOld->node2int_node_comm_index2[ i ];

  /* Available memory to create the groups */
  parmesh->memMax = parmesh->memCur;
  parmesh->listgrp[grpIdOld].mesh->memMax = parmesh->listgrp[grpIdOld].mesh->memCur;
  memAv = parmesh->memGloMax-parmesh->memMax-parmesh->listgrp[grpIdOld].mesh->memCur;

  for ( grpId = 0; grpId < ngrp; ++grpId ) {
    /** New group */
    grpCur  = &grpsNew[grpId];

    /** New group initialisation */
    if ( !PMMG_splitGrps_newGroup(parmesh,&grpsNew[grpId],grpIdOld,&memAv,
                                  countPerGrp[grpId],&f2ifc_max,&n2inc_max) ) {
      fprintf(stderr,"\n  ## Error: %s: unable to initialize new"
              " group (%d).\n",__func__,grpId);
      ret_val = -1;
      goto fail_sgrp;
    }
  }

  /* Store index of the old tetra in the new tetra field, for all grps */
  for ( ie = 1; ie <= meshOld->ne; ie++ ) {
    grpId = part[ ie-1 ];
    grpCur  = &grpsNew[grpId];
    meshCur = grpCur->mesh;
    meshCur->tetra[ meshOld->tetra[ie].flag ].flag = ie;
  }

  for ( grpId = 0; grpId < ngrp; ++grpId ) {
    /** New group filling */
    grpCur  = &grpsNew[grpId];
    meshCur = grpCur->mesh;

    if ( !PMMG_splitGrps_fillGroup(parmesh,&grpsNew[grpId],grpIdOld,grpId,
                                   countPerGrp[grpId],&poiPerGrp,&f2ifc_max,
                                   &n2inc_max,part,posInIntFaceComm,
                                   iplocFaceComm,&memAv) ) {
      fprintf(stderr,"\n  ## Error: %s: unable to fill new group (%d).\n",
              __func__,grpId);
      ret_val = -1;
      goto fail_sgrp;
    }

    /* Mesh cleaning in the new group */
    if ( !PMMG_splitGrps_cleanMesh(meshCur,grpCur->met,poiPerGrp) ) {
      fprintf(stderr,"\n  ## Error: %s: unable to clean the mesh of"
              " new group (%d).\n",__func__,grpId);
      ret_val = -1;
      goto fail_sgrp;
    }
    /* Remove the memory used by this mesh from the available memory */
    if ( memAv < meshCur->memMax ) {
      fprintf(stderr,"\n  ## Error: %s: not enough memory to allocate a new"
              " mesh.\n",__func__);
      goto fail_sgrp;
    }
    memAv -= meshCur->memMax;

    /* Give the memory to the parmesh */
    oldMemMax        = parmesh->memMax;
    parmesh->memMax += memAv;

    /* Fitting of the communicator sizes */
    PMMG_RECALLOC(parmesh, grpCur->node2int_node_comm_index1,
                  grpCur->nitem_int_node_comm, n2inc_max, int,
                  "subgroup internal1 communicator ",
                  ret_val = -1;goto fail_sgrp );
    PMMG_RECALLOC(parmesh, grpCur->node2int_node_comm_index2,
                  grpCur->nitem_int_node_comm, n2inc_max, int,
                  "subgroup internal2 communicator ",
                  ret_val = -1;goto fail_sgrp );
    n2inc_max = grpCur->nitem_int_node_comm;
    PMMG_RECALLOC(parmesh, grpCur->face2int_face_comm_index1,
                  grpCur->nitem_int_face_comm, f2ifc_max, int,
                  "subgroup interface faces communicator ",
                  ret_val = -1;goto fail_sgrp );
    PMMG_RECALLOC(parmesh, grpCur->face2int_face_comm_index2,
                  grpCur->nitem_int_face_comm, f2ifc_max, int,
                  "subgroup interface faces communicator ",
                  ret_val = -1;goto fail_sgrp );

    parmesh->memMax = parmesh->memCur;
    assert ( parmesh->memMax <= oldMemMax );
    if ( memAv < oldMemMax - parmesh->memMax ) {
      fprintf(stderr,"\n  ## Error: %s: not enough memory.\n",__func__);
      goto fail_sgrp;
    }
    memAv          -= (parmesh->memMax - oldMemMax);
  }

  /* No error so far, skip deallocation of lstgrps */
  goto fail_facePos;

  /* fail_sgrp deallocates any mesh that has been allocated in listgroup.
     Should be executed only if an error has occured */
fail_sgrp:
  for ( grpId = 0; grpId < ngrp; ++grpId ) {
    grpCur  = &grpsNew[grpId];
    meshCur = grpCur->mesh;

    /* internal comm for nodes */
    if ( grpCur->node2int_node_comm_index2 != NULL )
      PMMG_DEL_MEM(parmesh,grpCur->node2int_node_comm_index2,int,
                   "subgroup internal2 communicator ");
    if ( grpCur->node2int_node_comm_index1 != NULL )
      PMMG_DEL_MEM(parmesh,grpCur->node2int_node_comm_index1,int,
                   "subgroup internal1 communicator ");

    /* internal communicator for faces */
    if ( grpCur->face2int_face_comm_index1 )
      PMMG_DEL_MEM(parmesh,grpCur->face2int_face_comm_index1,int,
                   "face2int_face_comm_index1 communicator ");
    if ( grpCur->face2int_face_comm_index2 )
      PMMG_DEL_MEM(parmesh,grpCur->face2int_face_comm_index2,int,
                   "face2int_face_comm_index1 communicator ");

    /* mesh */
    if ( meshCur ) {
      if ( meshCur->adja != NULL )
        PMMG_DEL_MEM(meshCur,meshCur->adja,int,"adjacency table");
      if ( meshCur->xpoint != NULL )
        PMMG_DEL_MEM(meshCur,meshCur->xpoint,MMG5_xPoint,"boundary points");
      if ( meshCur->xtetra != NULL )
        PMMG_DEL_MEM(meshCur,meshCur->xtetra,MMG5_xTetra,"msh boundary tetra");
    }
  }

  /* these labels should be executed as part of normal code execution before
     returning as well as error handling */
fail_facePos:
  PMMG_DEL_MEM(parmesh,iplocFaceComm,int,
               "starting vertices of the faces of face2int_face_comm_index1");

  PMMG_DEL_MEM(parmesh,posInIntFaceComm,int,
               "array to store faces positions in internal face communicator");

  return ret_val;
}


/**
 * \param parmesh pointer toward the parmesh structure.
 * \param target software for which we split the groups
 * (\a PMMG_GRPSPL_METIS_TARGET or \a PMMG_GRPSPL_MMG_TARGET)
 * \param fitMesh alloc the meshes at their exact sizes
 *
 * \return -1 : no possibility to save the mesh
 *         0  : failed but the mesh is correct
 *         1  : success
 *
 * if the existing group of only one mesh is too big, split it into into several
 * meshes.
 *
 * \warning tetra must be packed.
 *
 */
int PMMG_split_grps( PMMG_pParMesh parmesh,int target,int fitMesh)
{
  PMMG_pGrp grpOld;
  PMMG_pGrp grpsNew = NULL;
  MMG5_pMesh meshOld;
  int *countPerGrp = NULL;
  int ret_val = 1;
  /** remember meshOld->ne to correctly free the metis buffer */
  int meshOld_ne = 0;
  idx_t ngrp = 1;
  idx_t *part = NULL;
  int grpId,grpIdOld, tet;
  int ne_all[parmesh->nprocs],ngrps_all[parmesh->nprocs];
  int npmax,nemax,xpmax,xtmax;

  if ( !parmesh->ngrp ) goto end;

  /* We are splitting group 0 */
  grpIdOld = 0;
  grpOld = &parmesh->listgrp[grpIdOld];
  meshOld = parmesh->listgrp[grpIdOld].mesh;

  assert ( (parmesh->ngrp == 1) && " split_grps can not split m groups to n");

  if ( !meshOld ) goto end;

  ngrp = PMMG_howManyGroups( meshOld->ne,abs(parmesh->info.target_mesh_size) );
  if ( parmesh->info.target_mesh_size < 0 ) {
    /* default value : do not authorize large number of groups */
    ngrp = MG_MIN ( PMMG_REMESHER_NGRPS_MAX, ngrp );
  }

  if ( target == PMMG_GRPSPL_METIS_TARGET ) {
    /* Compute the number of metis nodes from the number of groups */
    ngrp = MG_MIN( ngrp*abs(parmesh->info.metis_ratio), meshOld->ne/PMMG_METIS_NELEM_MIN+1 );
    if ( parmesh->info.metis_ratio < 0 ) {
      /* default value : do not authorize large number of groups */
      if ( ngrp > PMMG_METIS_NGRPS_MAX ) {
        printf("  ## Warning: %s: too much metis nodes needed...\n"
               "     Partitions may remains freezed. Try to use more processors.\n",
               __func__);
        ngrp = PMMG_METIS_NGRPS_MAX;
      }
    }
    if ( ngrp > meshOld->ne ) {
      /* Correction if it leads to more groups than elements */
      printf("  ## Warning: %s: too much metis nodes needed...\n"
             "     Partitions may remains freezed. Try to reduce the number of processors.\n",
             __func__);
      ngrp = MG_MIN ( meshOld->ne, ngrp );
    }
  }

  /* Print split info */
  MPI_CHECK( MPI_Gather(&ngrp,1,MPI_INT,
                        ngrps_all,1,MPI_INT,0,parmesh->comm),
             PMMG_CLEAN_AND_RETURN(parmesh,PMMG_LOWFAILURE) );
  MPI_CHECK( MPI_Gather(&meshOld->ne,1,MPI_INT,
                        ne_all,1,MPI_INT,0,parmesh->comm),
             PMMG_CLEAN_AND_RETURN(parmesh,PMMG_LOWFAILURE) );
  if ( parmesh->info.imprim > PMMG_VERB_STEPS ) {
    int i;
    for( i=0; i<parmesh->nprocs; i++ )
      fprintf(stdout,"         rank %d splitting %d elts into %d grps\n",
              i,ne_all[i],ngrps_all[i]);
  }

  /* Does the group need to be further subdivided to subgroups or not? */
  if ( ngrp == 1 )  {
    if ( parmesh->ddebug )
      fprintf( stdout,
               "[%d-%d]: %d group is enough, no need to create sub groups.\n",
               parmesh->myrank+1, parmesh->nprocs, ngrp );
    goto end;

  } else {
    if ( parmesh->ddebug )
      fprintf( stdout,
               "[%d-%d]: %d groups required, splitting into sub groups...\n",
               parmesh->myrank+1, parmesh->nprocs, ngrp );
  }

  /* Crude check whether there is enough free memory to allocate the new group */
  if ( parmesh->memCur+2*meshOld->memCur>parmesh->memGloMax ) {
    npmax = meshOld->npmax;
    nemax = meshOld->nemax;
    xpmax = meshOld->xpmax;
    xtmax = meshOld->xtmax;
    meshOld->npmax = meshOld->np;
    meshOld->nemax = meshOld->ne;
    meshOld->xpmax = meshOld->xp;
    meshOld->xtmax = meshOld->xt;
    if ( (!PMMG_setMemMax_realloc( meshOld, npmax, xpmax, nemax, xtmax )) ||
         parmesh->memCur+2*meshOld->memCur>parmesh->memGloMax ) {
      fprintf( stderr, "Not enough memory to create listgrp struct\n" );
      return 0;
    }
  }

  /* use metis to partition the mesh into the computed number of groups needed
     part array contains the groupID computed by metis for each tetra */
  PMMG_CALLOC(parmesh,part,meshOld->ne,idx_t,"metis buffer ", return 0);
  meshOld_ne = meshOld->ne;

  if ( !PMMG_part_meshElts2metis(parmesh, part, ngrp) ) {
    ret_val = 0;
    goto fail_part;
  }

  /* count_per_grp: count new tetra per group, and store new ID in the old
   * tetra flag */
  PMMG_CALLOC(parmesh,countPerGrp,ngrp,int,"counter buffer ",
              ret_val = 0;goto fail_part);
  PMMG_splitGrps_countTetPerGrp( parmesh,meshOld,part,countPerGrp );
  
  /* Allocate list of subgroups struct and allocate memory */
  PMMG_CALLOC(parmesh,grpsNew,ngrp,PMMG_Grp,"subgourp list ",
              ret_val = 0; goto fail_counters);


  /** Perform group splitting */
  ret_val = PMMG_split_eachGrp( parmesh,grpIdOld,grpsNew,ngrp,countPerGrp,part );
  if( ret_val != 1) goto fail_counters;
 
//DEBUGGING:  saveGrpsToMeshes( grpsNew, ngrp, parmesh->myrank, "AfterSplitGrp" );

  PMMG_listgrp_free(parmesh, &parmesh->listgrp, parmesh->ngrp);
  parmesh->listgrp = grpsNew;
  parmesh->ngrp = ngrp;

  /** Check grps contiguity */
  ret_val = PMMG_checkAndReset_grps_contiguity( parmesh );

  if ( PMMG_parmesh_updateMemMax(parmesh, 5, fitMesh) ) {
    /* No error so far, skip deallocation of lstgrps */
    goto fail_counters;
  }
  else
    ret_val = -1;

fail_counters:
  PMMG_DEL_MEM(parmesh,countPerGrp,int,"counter buffer ");
fail_part:
  PMMG_DEL_MEM(parmesh,part,idx_t,"free metis buffer ");

end:
  /* Check the communicators */
  assert ( PMMG_check_intNodeComm(parmesh) && "Wrong internal node comm" );
  assert ( PMMG_check_intFaceComm(parmesh) && "Wrong internal face comm" );
  assert ( PMMG_check_extNodeComm(parmesh) && "Wrong external node comm" );
  assert ( PMMG_check_extFaceComm(parmesh) && "Wrong external face comm" );

  return ret_val;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param target software for which we split the groups
 * (\a PMMG_GRPSPL_METIS_TARGET or \a PMMG_GRPSPL_MMG_TARGET)
 * \param fitMesh alloc the meshes at their exact sizes
 *
 * \return 0 if fail, 1 if success, -1 if the mesh is not correct
 *
 * Redistribute the n groups of listgrps into \a target_mesh_size groups.
 *
 */
int PMMG_split_n2mGrps(PMMG_pParMesh parmesh,int target,int fitMesh) {
  int     ier,ier1;
#ifndef NDEBUG
  int     ier_glob;
#endif
  int     tim;
  mytime  ctim[3];
  char    stim[32];

  assert ( PMMG_check_intFaceComm ( parmesh ) );
  assert ( PMMG_check_extFaceComm ( parmesh ) );
  assert ( PMMG_check_intNodeComm ( parmesh ) );
  assert ( PMMG_check_extNodeComm ( parmesh ) );

  if ( parmesh->info.imprim > PMMG_VERB_DETQUAL ) {
      tminit(ctim,3);
      tim = 0;
      chrono(ON,&(ctim[tim]));
  }

  /** Merge the parmesh groups into 1 group */
  ier = PMMG_merge_grps(parmesh);
  if ( !ier ) {
    fprintf(stderr,"\n  ## Merge groups problem.\n");
  }

  if ( parmesh->info.imprim > PMMG_VERB_DETQUAL ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"\n                   merge group           %s\n",stim);

    tim = 1;
    chrono(ON,&(ctim[tim]));
  }

  /** Pack the tetra and update the face communicator */
  ier1 = 1;
  if ( parmesh->ngrp ) {
    ier1 = PMMG_packTetra(parmesh,0);
    if ( !ier1 ) {
      fprintf(stderr,"\n  ## Pack tetrahedra and face communicators problem.\n");
    }
  }
  ier = MG_MIN( ier, ier1 );

  if ( parmesh->info.imprim > PMMG_VERB_DETQUAL ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"                   pack tetra            %s\n",stim);
  }

#ifndef NDEBUG
  /* In debug mode we have mpi comm in split_grps, thus, if 1 proc fails
   * and the other not we will deadlock */
  MPI_Allreduce( &ier, &ier_glob, 1, MPI_INT, MPI_MIN, parmesh->comm);
  if ( ier_glob <=0 ) return ier_glob;
#endif

  if ( parmesh->ddebug ) {
    size_t dummy1,dummy2;

    PMMG_TRANSFER_AVMEM_TO_MESHES(parmesh);

    PMMG_qualhisto( parmesh, PMMG_INQUA );
    PMMG_prilen( parmesh, 0 );

    PMMG_TRANSFER_AVMEM_TO_PARMESH(parmesh,dummy1,dummy2);

  }

  if ( parmesh->info.imprim > PMMG_VERB_DETQUAL ) {
    tim = 2;
    chrono(ON,&(ctim[tim]));
  }

  /** Split the group into the suitable number of groups */
  if ( ier )
    ier = PMMG_split_grps(parmesh,target,fitMesh);

  if ( parmesh->info.imprim > PMMG_VERB_DETQUAL ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"                   group split           %s\n",stim);
  }

  if ( ier<=0 )
    fprintf(stderr,"\n  ## Split group problem.\n");

  assert ( PMMG_check_intFaceComm ( parmesh ) );
  assert ( PMMG_check_extFaceComm ( parmesh ) );
  assert ( PMMG_check_intNodeComm ( parmesh ) );
  assert ( PMMG_check_extNodeComm ( parmesh ) );

  return ier;
}

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

/**
 * \param .
 *
 *
 */
int PMMG_Set_vertex(PMMG_pGrp grp, double c0, double c1, double c2,
                    int ref, int pos){
  if(!MMG3D_Set_vertex(grp->mesh, c0, c1, c2, ref, pos)){
    fprintf(stderr,"  ## Error in setting mesh vertex.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_vertices(PMMG_pGrp grp, double *vertices,int *refs){
  if(!MMG3D_Set_vertices(grp->mesh, vertices, refs)){
    fprintf(stderr,"  ## Error in setting mesh vertices.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_tetrahedron(PMMG_pGrp grp, int v0, int v1, int v2, int v3,
                         int ref, int pos){
  if(!MMG3D_Set_tetrahedron(grp->mesh, v0, v1, v2, v3, ref, pos)){
    fprintf(stderr,"  ## Error in setting mesh tetrahedron.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_tetrahedra(PMMG_pGrp grp, int *tetra, int *refs){
  if(!MMG3D_Set_tetrahedra(grp->mesh, tetra, refs)){
    fprintf(stderr,"  ## Error in setting mesh tetrahedra.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_prism(PMMG_pGrp grp, int v0, int v1, int v2,
                   int v3, int v4, int v5, int ref, int pos){
  if(!MMG3D_Set_prism(grp->mesh, v0, v1, v2, v3, v4, v5, ref, pos)){
    fprintf(stderr,"  ## Error in setting mesh prism.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_prisms(PMMG_pGrp grp, int *prisms, int *refs){
  if(!MMG3D_Set_prisms(grp->mesh, prisms, refs)){
    fprintf(stderr,"  ## Error in setting mesh prisms.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_triangle(PMMG_pGrp grp, int v0, int v1, int v2,
                      int ref,int pos){
  if(!MMG3D_Set_triangle(grp->mesh, v0, v1, v2, ref, pos)){
    fprintf(stderr,"  ## Error in setting mesh triangle.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_triangles(PMMG_pGrp grp, int *tria, int *refs){
  if(!MMG3D_Set_triangles(grp->mesh, tria, refs)){
    fprintf(stderr,"  ## Error in setting mesh triangles.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_quadrilateral(PMMG_pGrp grp, int v0, int v1,
                           int v2, int v3, int ref,int pos){
  if(!MMG3D_Set_quadrilateral(grp->mesh, v0, v1, v2, v3, ref, pos)){
    fprintf(stderr,"  ## Error in setting mesh quadrilateral.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_quadrilaterals(PMMG_pGrp grp, int *quads, int *refs){
  if(!MMG3D_Set_quadrilaterals(grp->mesh, quads, refs)){
    fprintf(stderr,"  ## Error in setting mesh quadrilaterals.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_edge(PMMG_pGrp grp, int v0, int v1, int ref, int pos){
  if(!MMG3D_Set_edge(grp->mesh, v0, v1, ref, pos)){
    fprintf(stderr,"  ## Error in setting mesh edge.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_corner(PMMG_pGrp grp, int k){
  if(!MMG3D_Set_corner(grp->mesh, k)){
    fprintf(stderr,"  ## Error in setting corner.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_requiredVertex(PMMG_pGrp grp, int k){
  if(!MMG3D_Set_requiredVertex(grp->mesh, k)){
    fprintf(stderr,"  ## Error in setting required vertex.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_requiredTetrahedron(PMMG_pGrp grp, int k){
  if(!MMG3D_Set_requiredTetrahedron(grp->mesh, k)){
    fprintf(stderr,"  ## Error in setting required tetrahedron.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_requiredTetrahedra(PMMG_pGrp grp, int *reqIdx, int nreq){
  if(!MMG3D_Set_requiredTetrahedra(grp->mesh, reqIdx, nreq)){
    fprintf(stderr,"  ## Error in setting required tetrahedra.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_requiredTriangle(PMMG_pGrp grp, int k){
  if(!MMG3D_Set_requiredTriangle(grp->mesh, k)){
    fprintf(stderr,"  ## Error in setting required triangle.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_requiredTriangles(PMMG_pGrp grp, int *reqIdx, int nreq){
  if(!MMG3D_Set_requiredTriangles(grp->mesh, reqIdx, nreq)){
    fprintf(stderr,"  ## Error in setting required triangles.\n");
    return(0);
  }
  return(1);
}


/**
 * \param .
 *
 *
 */
int PMMG_Set_ridge(PMMG_pGrp grp, int k){
  if(!MMG3D_Set_ridge(grp->mesh, k)){
    fprintf(stderr,"  ## Error in setting ridge.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_requiredEdge(PMMG_pGrp grp, int k){
  if(!MMG3D_Set_requiredEdge(grp->mesh, k)){
    fprintf(stderr,"  ## Error in setting required edge.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_normalAtVertex(PMMG_pGrp grp, int k, double n0, double n1,
                            double n2){
  if(!MMG3D_Set_normalAtVertex(grp->mesh, k, n0, n1, n2)){
    fprintf(stderr,"  ## Error in setting normal at vertex.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_scalarSol(PMMG_pGrp grp, double s,int pos){
  if(!MMG3D_Set_scalarSol(grp->sol, s, pos)){
    fprintf(stderr,"  ## Error in setting scalar solution at vertex.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_scalarSols(PMMG_pGrp grp, double *sol){
  if(!MMG3D_Set_scalarSols(grp->sol, sol)){
    fprintf(stderr,"  ## Error in setting scalar solution at vertices.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_vectorSol(PMMG_pGrp grp, double vx,double vy, double vz,
                       int pos){
  if(!MMG3D_Set_vectorSol(grp->sol, vx, vy, vz, pos)){
    fprintf(stderr,"  ## Error in setting vectorial solution at vertex.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_vectorSols(PMMG_pGrp grp, double *sols){
  if(!MMG3D_Set_vectorSols(grp->sol, sols)){
    fprintf(stderr,"  ## Error in setting vectorial solution at vertices.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_tensorSol(PMMG_pGrp grp, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos){
  if(!MMG3D_Set_tensorSol(grp->sol, m11, m12, m13, m22, m23, m33, pos)){
    fprintf(stderr,"  ## Error in setting tensorial solution at vertex.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_tensorSols(PMMG_pGrp grp, double *sols){
  if(!MMG3D_Set_tensorSols(grp->sol, sols)){
    fprintf(stderr,"  ## Error in setting tensorial solution at vertices.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_scalarMet(PMMG_pGrp grp, double m,int pos){
  if(!MMG3D_Set_scalarSol(grp->met, m, pos)){
    fprintf(stderr,"  ## Error in setting scalar metrics at vertex.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_scalarMets(PMMG_pGrp grp, double *met){
  if(!MMG3D_Set_scalarSols(grp->met, met)){
    fprintf(stderr,"  ## Error in setting scalar metrics at vertices.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_vectorMet(PMMG_pGrp grp, double vx,double vy, double vz,
                       int pos){
  if(!MMG3D_Set_vectorSol(grp->met, vx, vy, vz, pos)){
    fprintf(stderr,"  ## Error in setting vectorial metrics at vertex.\n");
    return(0);
  }

  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_vectorMets(PMMG_pGrp grp, double *mets){
  if(!MMG3D_Set_vectorSols(grp->met, mets)){
    fprintf(stderr,"  ## Error in setting vectorial metrics at vertices.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_tensorMet(PMMG_pGrp grp, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos){
  if(!MMG3D_Set_tensorSol(grp->met, m11, m12, m13, m22, m23, m33, pos)){
    fprintf(stderr,"  ## Error in setting tensorial metrics at vertex.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Set_tensorMets(PMMG_pGrp grp, double *mets){
  if(!MMG3D_Set_tensorSols(grp->met, mets)){
    fprintf(stderr,"  ## Error in setting tensorial metrics at vertices.\n");
    return(0);
  }
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_grpSize(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_vertex(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_vertices(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_tetrahedron(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_tetrahedra(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_prism(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_prisms(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_triangle(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_triangles(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_quadrilateral(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_quadrilaterals(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_edge(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_normalAtVertex(PMMG_pGrp grp){
  return(1);
}


/**
 * \param .
 *
 *
 */
int PMMG_Get_scalarSol(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_scalarSols(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_vectorSol(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_vectorSols(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_tensorSol(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_tensorSols(PMMG_pGrp grp){
  return(1);
}


/**
 * \param .
 *
 *
 */
int PMMG_Get_scalarMet(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_scalarMets(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_vectorMet(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_vectorMets(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_tensorMet(PMMG_pGrp grp){
  return(1);
}

/**
 * \param .
 *
 *
 */
int PMMG_Get_tensorMets(PMMG_pGrp grp){
  return(1);
}
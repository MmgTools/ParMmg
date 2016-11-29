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
 * Set group size.
 */
int PMMG_Set_grpSize(PMMG_pGrp grp, int np, int ne, int nprism, int nt,
                     int nquad, int na, int typEntity, int typSol, int typMet,
                     int *loc2glob){
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
 * Set group mesh vertex.
 */
int PMMG_Set_vertex(PMMG_pGrp grp, double c0, double c1, double c2,
                    int ref, int pos){
  return(MMG3D_Set_vertex(grp->mesh, c0, c1, c2, ref, pos));
}

/**
 * Set group mesh vertices.
 */
int PMMG_Set_vertices(PMMG_pGrp grp, double *vertices,int *refs){
  return(MMG3D_Set_vertices(grp->mesh, vertices, refs));
}

/**
 * Set group mesh tetrahedron.
 */
int PMMG_Set_tetrahedron(PMMG_pGrp grp, int v0, int v1, int v2, int v3,
                         int ref, int pos){
  return(MMG3D_Set_tetrahedron(grp->mesh, v0, v1, v2, v3, ref, pos));
}

/**
 * Set group mesh tetrahedra.
 */
int PMMG_Set_tetrahedra(PMMG_pGrp grp, int *tetra, int *refs){
  return(MMG3D_Set_tetrahedra(grp->mesh, tetra, refs));
}

/**
 * Set group mesh prism.
 */
int PMMG_Set_prism(PMMG_pGrp grp, int v0, int v1, int v2,
                   int v3, int v4, int v5, int ref, int pos){
  return(MMG3D_Set_prism(grp->mesh, v0, v1, v2, v3, v4, v5, ref, pos));
}

/**
 * Set group mesh prisms.
 */
int PMMG_Set_prisms(PMMG_pGrp grp, int *prisms, int *refs){
  return(MMG3D_Set_prisms(grp->mesh, prisms, refs));
}

/**
 * Set group mesh triangle.
 */
int PMMG_Set_triangle(PMMG_pGrp grp, int v0, int v1, int v2,
                      int ref,int pos){
  return(MMG3D_Set_triangle(grp->mesh, v0, v1, v2, ref, pos));
}

/**
 * Set group mesh triangles.
 */
int PMMG_Set_triangles(PMMG_pGrp grp, int *tria, int *refs){
  return(MMG3D_Set_triangles(grp->mesh, tria, refs));
}

/**
 * Set group mesh quadrilateral.
 */
int PMMG_Set_quadrilateral(PMMG_pGrp grp, int v0, int v1,
                           int v2, int v3, int ref,int pos){
  return(MMG3D_Set_quadrilateral(grp->mesh, v0, v1, v2, v3, ref, pos));
}

/**
 * Set group mesh quadrilaterals.
 */
int PMMG_Set_quadrilaterals(PMMG_pGrp grp, int *quads, int *refs){
  return(MMG3D_Set_quadrilaterals(grp->mesh, quads, refs));
}

/**
 * Set group mesh edge.
 */
int PMMG_Set_edge(PMMG_pGrp grp, int v0, int v1, int ref, int pos){
  return(MMG3D_Set_edge(grp->mesh, v0, v1, ref, pos));
}

/**
 * Set group mesh corner.
 */
int PMMG_Set_corner(PMMG_pGrp grp, int k){
  return(MMG3D_Set_corner(grp->mesh, k));
}

/**
 * Set group mesh required vertex.
 */
int PMMG_Set_requiredVertex(PMMG_pGrp grp, int k){
  return(MMG3D_Set_requiredVertex(grp->mesh, k));
}

/**
 * Set group mesh required tetrahedron.
 */
int PMMG_Set_requiredTetrahedron(PMMG_pGrp grp, int k){
  return(MMG3D_Set_requiredTetrahedron(grp->mesh, k));
}

/**
 * Set group mesh required tetrahedra.
 */
int PMMG_Set_requiredTetrahedra(PMMG_pGrp grp, int *reqIdx, int nreq){
  return(MMG3D_Set_requiredTetrahedra(grp->mesh, reqIdx, nreq));
}

/**
 * Set group mesh required triangle.
 */
int PMMG_Set_requiredTriangle(PMMG_pGrp grp, int k){
  return(MMG3D_Set_requiredTriangle(grp->mesh, k));
}

/**
 * Set group mesh required triangles.
 */
int PMMG_Set_requiredTriangles(PMMG_pGrp grp, int *reqIdx, int nreq){
  return(MMG3D_Set_requiredTriangles(grp->mesh, reqIdx, nreq));
}


/**
 * Set group mesh ridge.
 */
int PMMG_Set_ridge(PMMG_pGrp grp, int k){
  return(MMG3D_Set_ridge(grp->mesh, k));
}

/**
 * Set group mesh required edge.
 */
int PMMG_Set_requiredEdge(PMMG_pGrp grp, int k){
  return(MMG3D_Set_requiredEdge(grp->mesh, k));
}

/**
 * Set group mesh normal at vertex.
 */
int PMMG_Set_normalAtVertex(PMMG_pGrp grp, int k, double n0, double n1,
                            double n2){
  return(MMG3D_Set_normalAtVertex(grp->mesh, k, n0, n1, n2));
}

/**
 * Set group scalar solution at vertex.
 */
int PMMG_Set_scalarSol(PMMG_pGrp grp, double s,int pos){
  return(MMG3D_Set_scalarSol(grp->sol, s, pos));
}

/**
 * Set group scalar solution at vertices.
 */
int PMMG_Set_scalarSols(PMMG_pGrp grp, double *sol){
  return(MMG3D_Set_scalarSols(grp->sol, sol));
}

/**
 * Set group vectorial solution at vertex.
 */
int PMMG_Set_vectorSol(PMMG_pGrp grp, double vx,double vy, double vz,
                       int pos){
  return(MMG3D_Set_vectorSol(grp->sol, vx, vy, vz, pos));
}

/**
 * Set group vectorial solution at vertices.
 */
int PMMG_Set_vectorSols(PMMG_pGrp grp, double *sols){
  return(MMG3D_Set_vectorSols(grp->sol, sols));
}

/**
 * Set group tensorial solution at vertex.
 */
int PMMG_Set_tensorSol(PMMG_pGrp grp, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos){
  return(MMG3D_Set_tensorSol(grp->sol, m11, m12, m13, m22, m23, m33, pos));
}

/**
 * Set group tensorial solution at vertices.
 */
int PMMG_Set_tensorSols(PMMG_pGrp grp, double *sols){
  return(MMG3D_Set_tensorSols(grp->sol, sols));
}

/**
 * Set group scalar metrics at vertex.
 */
int PMMG_Set_scalarMet(PMMG_pGrp grp, double m,int pos){
  return(MMG3D_Set_scalarSol(grp->met, m, pos));
}

/**
 * Set group scalar metrics at vertices.
 */
int PMMG_Set_scalarMets(PMMG_pGrp grp, double *met){
  return(MMG3D_Set_scalarSols(grp->met, met));
}

/**
 * Set group vectorial metrics at vertex.
 */
int PMMG_Set_vectorMet(PMMG_pGrp grp, double vx,double vy, double vz,
                       int pos){
  return(MMG3D_Set_vectorSol(grp->met, vx, vy, vz, pos));
}

/**
 * Set group vectorial metrics at vertices.
 */
int PMMG_Set_vectorMets(PMMG_pGrp grp, double *mets){
  return(MMG3D_Set_vectorSols(grp->met, mets));
}

/**
 * Set group tensorial metrics at vertex.
 */
int PMMG_Set_tensorMet(PMMG_pGrp grp, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos){
  return(MMG3D_Set_tensorSol(grp->met, m11, m12, m13, m22, m23, m33, pos));
}

/**
 * Set group tensorial metrics at vertices.
 */
int PMMG_Set_tensorMets(PMMG_pGrp grp, double *mets){
  return(MMG3D_Set_tensorSols(grp->met, mets));
}

/**
 * Get group sizes.
 */
int PMMG_Get_grpSize(PMMG_pGrp grp, int* np, int* ne,int *nprism, int* nt,
                     int* nquad, int* na, int* typEntity,
                     int* typSol, int* typMet){
  if(!MMG3D_Get_meshSize(grp->mesh, np, ne, nprism, nt, nquad, na)){
    fprintf(stderr,"  ## Error in getting group mesh size.\n");
    return(0);
  }
  if(!MMG3D_Get_solSize(grp->mesh, grp->sol, typEntity, np, typSol)){
    // FIXME np is read again.
    fprintf(stderr,"  ## Error in getting group solution size.\n");
    return(0);
  }
  if(!MMG3D_Get_solSize(grp->mesh, grp->met, typEntity, np, typMet)){
    // FIXME np is read again, as typEntity.
    fprintf(stderr,"  ## Error in getting group metrics size.\n");
    return(0);
  }
  return(1);
}

/**
 * Get group mesh vertex.
 */
int PMMG_Get_vertex(PMMG_pGrp grp, double* c0, double* c1, double* c2,
                    int* ref,int* isCorner, int* isRequired){
  return(MMG3D_Get_vertex(grp->mesh, c0, c1, c2, ref, isCorner, isRequired));
}

/**
 * Get group mesh vertices.
 */
int PMMG_Get_vertices(PMMG_pGrp grp, double* vertices, int* refs,
                      int* areCorners, int* areRequired){
  return(MMG3D_Get_vertices(grp->mesh, vertices, refs, areCorners, areRequired));
}

/**
 * Get group mesh tetrahedron.
 */
int PMMG_Get_tetrahedron(PMMG_pGrp grp, int* v0, int* v1, int* v2,
                         int* v3,int* ref, int* isRequired){
  return(MMG3D_Get_tetrahedron(grp->mesh, v0, v1, v2, v3, ref, isRequired));
}

/**
 * Get group mesh tetrahedra.
 */
int PMMG_Get_tetrahedra(PMMG_pGrp grp, int* tetra,int* refs,
                        int* areRequired){
  return(MMG3D_Get_tetrahedra(grp->mesh, tetra, refs, areRequired));
}

/**
 * Get group mesh prism.
 */
int PMMG_Get_prism(PMMG_pGrp grp, int* v0, int* v1, int* v2,
                   int* v3,int* v4,int* v5,int* ref, int* isRequired){
  return(MMG3D_Get_prism(grp->mesh, v0, v1, v2, v3, v4, v5, ref, isRequired));
}

/**
 * Get group mesh prisms.
 */
int PMMG_Get_prisms(PMMG_pGrp grp, int* prisms,int* refs,
                    int* areRequired){
  return(MMG3D_Get_prisms(grp->mesh, prisms, refs, areRequired));
}

/**
 * Get group mesh triangle.
 */
int PMMG_Get_triangle(PMMG_pGrp grp, int* v0, int* v1, int* v2, int* ref,
                      int* isRequired){
  return(MMG3D_Get_triangle(grp->mesh, v0, v1, v2, ref, isRequired));
}

/**
 * Get group mesh triangles.
 */
int PMMG_Get_triangles(PMMG_pGrp grp, int* tria, int* refs,
                       int* areRequired){
  return(MMG3D_Get_triangles(grp->mesh, tria, refs, areRequired));
}

/**
 * Get group mesh quadrilateral.
 */
int PMMG_Get_quadrilateral(PMMG_pGrp grp, int* v0, int* v1, int* v2,int* v3,
                            int* ref, int* isRequired){
  return(MMG3D_Get_quadrilateral(grp->mesh, v0, v1, v2, v3, ref, isRequired));
}

/**
 * Get group mesh quadrilaterals.
 */
int PMMG_Get_quadrilaterals(PMMG_pGrp grp, int* quads, int* refs,
                            int* areRequired){
  return(MMG3D_Get_quadrilaterals(grp->mesh, quads, refs, areRequired));
}

/**
 * Get group mesh edge.
 */
int PMMG_Get_edge(PMMG_pGrp grp, int* e0, int* e1, int* ref,
                   int* isRidge, int* isRequired){
  return(MMG3D_Get_edge(grp->mesh, e0, e1, ref, isRidge, isRequired));
}

/**
 * Get group mesh normal at vertex.
 */
int PMMG_Get_normalAtVertex(PMMG_pGrp grp, int k, double *n0, double *n1,
                              double *n2){
  return(MMG3D_Get_normalAtVertex(grp->mesh, k, n0, n1, n2));
}


/**
 * Get group scalar solution at vertex.
 */
int PMMG_Get_scalarSol(PMMG_pGrp grp, double* s){
  return(MMG3D_Get_scalarSol(grp->sol, s));
}

/**
 * Get group scalar solution at vertices.
 */
int PMMG_Get_scalarSols(PMMG_pGrp grp, double* s){
  return(MMG3D_Get_scalarSols(grp->sol, s));
}

/**
 * Get group vectorial solution at vertex.
 */
int PMMG_Get_vectorSol(PMMG_pGrp grp, double* vx, double* vy, double* vz){
  return(MMG3D_Get_vectorSol(grp->sol, vx, vy, vz));
}

/**
 * Get group vectorial solution at vertices.
 */
int PMMG_Get_vectorSols(PMMG_pGrp grp, double* sols){
  return(MMG3D_Get_vectorSols(grp->sol, sols));
}

/**
 * Get group tensorial solution at vertex.
 */
int PMMG_Get_tensorSol(PMMG_pGrp grp, double *m11,double *m12, double *m13,
                       double *m22,double *m23, double *m33){
  return(MMG3D_Get_tensorSol(grp->sol, m11, m12, m13, m22, m23, m33));
}

/**
 * Get group tensorial solution at vertices.
 */
int PMMG_Get_tensorSols(PMMG_pGrp grp, double *sols){
  return(MMG3D_Get_tensorSols(grp->sol, sols));
}


/**
 * Get group scalar metrics at vertex.
 */
int PMMG_Get_scalarMet(PMMG_pGrp grp, double* m){
  return(MMG3D_Get_scalarSol(grp->met, m));
}

/**
 * Get group scalar metrics at vertices.
 */
int PMMG_Get_scalarMets(PMMG_pGrp grp, double* m){
  return(MMG3D_Get_scalarSols(grp->met, m));
}

/**
 * Get group vectorial metrics at vertex.
 */
int PMMG_Get_vectorMet(PMMG_pGrp grp, double* vx, double* vy, double* vz){
  return(MMG3D_Get_vectorSol(grp->met, vx, vy, vz));
}

/**
 * Get group vectorial metrics at vertices.
 */
int PMMG_Get_vectorMets(PMMG_pGrp grp, double* mets){
  return(MMG3D_Get_vectorSols(grp->met, mets));
}

/**
 * Get group tensorial metrics at vertex.
 */
int PMMG_Get_tensorMet(PMMG_pGrp grp, double *m11,double *m12, double *m13,
                       double *m22,double *m23, double *m33){
  return(MMG3D_Get_tensorSol(grp->met, m11, m12, m13, m22, m23, m33));
}

/**
 * Get group tensorial metrics at vertices.
 */
int PMMG_Get_tensorMets(PMMG_pGrp grp, double *mets){
  return(MMG3D_Get_tensorSols(grp->met, mets));
}

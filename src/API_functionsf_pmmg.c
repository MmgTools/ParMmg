/**
 * \file API_functionsf_pmmg.c
 * \brief Fortran API functions for PARMMG library.
 *
 * Define the Fortran API functions for PARMMG library: adds function
 * definitions with upcase, underscore and double underscore to match
 * any fortran compiler.
 *
 */
#include "libparmmg.h"
#include "mmgcommon.h"

/**
 * See \ref PMMG_Init_parMesh function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_INIT_PARMESH,pmmg_init_parmesh,
             (PMMG_pParMesh *parmesh,MPI_Comm *comm,int *retval),
             (parmesh,comm,retval)) {

  *retval = PMMG_Init_parMesh(parmesh,*comm);

  return;
}

/**
 * See \ref PMMG_Init_parameters function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_INIT_PARAMETERS,pmmg_init_parameters,
             (PMMG_pParMesh *parmesh,MPI_Comm *comm),
             (parmesh,comm)) {

  PMMG_Init_parameters(*parmesh,*comm);

  return;
}

/**
 * See \ref PMMG_Set_inputMeshName function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_INPUTMESHNAME, pmmg_set_inputmeshname,
             (PMMG_pParMesh *parmesh, char* meshin, int *strlen, int* retval),
             (parmesh,meshin,strlen,retval)) {
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,meshin,*strlen);
  tmp[*strlen] = '\0';
  *retval = PMMG_Set_inputMeshName(*parmesh,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref PMMG_Set_inputSolName function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_INPUTSOLNAME, pmmg_set_inputsolname,
             (PMMG_pParMesh *parmesh,char* solin, int* strlen, int* retval),
             (parmesh,solin,strlen,retval)) {

  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,solin,*strlen);
  tmp[*strlen] = '\0';
  *retval = PMMG_Set_inputSolName(*parmesh,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref PMMG_Set_inputMetName function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_INPUTMETNAME, pmmg_set_inputmetname,
             (PMMG_pParMesh *parmesh,char* metin, int* strlen, int* retval),
             (parmesh,metin,strlen,retval)) {

  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,metin,*strlen);
  tmp[*strlen] = '\0';
  *retval = PMMG_Set_inputMetName(*parmesh,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref PMMG_Set_outputMeshName function in libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_OUTPUTMESHNAME,pmmg_set_outputmeshname,
             (PMMG_pParMesh *parmesh, char* meshout, int* strlen,int* retval),
             (parmesh,meshout,strlen,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,meshout,*strlen);
  tmp[*strlen] = '\0';
  *retval = PMMG_Set_outputMeshName(*parmesh, tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref PMMG_Set_outputSolName function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_OUTPUTSOLNAME,pmmg_set_outputsolname,
             (PMMG_pParMesh *parmesh, char* solout,int* strlen, int* retval),
             (parmesh,solout,strlen,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,solout,*strlen);
  tmp[*strlen] = '\0';
  *retval = PMMG_Set_outputSolName(*parmesh,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref PMMG_Set_outputMetName function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_OUTPUTMETNAME,pmmg_set_outputmetname,
             (PMMG_pParMesh *parmesh, char* metout,int* strlen, int* retval),
             (parmesh,metout,strlen,retval)){
  char *tmp = NULL;

  tmp = (char*)malloc((*strlen+1)*sizeof(char));
  strncpy(tmp,metout,*strlen);
  tmp[*strlen] = '\0';
  *retval = PMMG_Set_outputMetName(*parmesh,tmp);
  _MMG5_SAFE_FREE(tmp);

  return;
}

/**
 * See \ref PMMG_Set_grpSize function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_GRPSIZE, pmmg_set_grpsize,
    (PMMG_pGrp *grp, int* np, int* ne, int* nprism, int* nt,
     int* nquad, int* na, int* typEntity, int* typSol, int* typMet,
     int* loc2glob, int* retval),
    (grp, np, ne, nprism, nt, nquad, na, typEntity, typSol, typMet,
     loc2glob, retval)) {
  *retval = PMMG_Set_grpSize(*grp, *np, *ne, *nprism, *nt, *nquad,
                             *na, *typEntity, *typSol, *typMet, loc2glob);
  return;
}

/**
 * See \ref PMMG_Set_vertex function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_VERTEX,pmmg_set_vertex,
             (PMMG_pGrp *grp, double* c0, double* c1, double* c2, int* ref,
              int* pos, int* retval),
             (grp,c0,c1,c2,ref,pos,retval)) {

  *retval = PMMG_Set_vertex(*grp,*c0,*c1,*c2,*ref,*pos);
  return;
}

/**
 * See \ref PMMG_Set_vertices function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_VERTICES,pmmg_set_vertices,
             (PMMG_pGrp *grp, double* vertices, int* refs, int* retval),
             (grp,vertices,refs,retval)) {

  *retval = PMMG_Set_vertices(*grp,vertices,refs);
  return;
}

/**
 * See \ref PMMG_Set_tetrahedron function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_TETRAHEDRON,pmmg_set_tetrahedron,
             (PMMG_pGrp *grp, int *v0, int *v1, int *v2, int *v3, int *ref,
              int *pos, int* retval),
             (grp,v0,v1,v2,v3,ref,pos,retval)){
  *retval = PMMG_Set_tetrahedron(*grp,*v0,*v1,*v2,*v3,*ref,*pos);
  return;
}

/**
 * See \ref PMMG_Set_tetrahedra function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_TETRAHEDRA,pmmg_set_tetrahedra,
             (PMMG_pGrp *grp, int *tetra, int *refs, int* retval),
             (grp,tetra,refs,retval)){
  *retval = PMMG_Set_tetrahedra(*grp,tetra,refs);
  return;
}

/**
 * See \ref PMMG_Set_prism function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_PRISM,pmmg_set_prism,
             (PMMG_pGrp *grp, int *v0, int *v1, int *v2, int *v3,
              int *v4, int *v5,int *ref,int *pos, int* retval),
             (grp,v0,v1,v2,v3,v4,v5,ref,pos,retval)){
  *retval = PMMG_Set_prism(*grp,*v0,*v1,*v2,*v3,*v4,*v5,*ref,*pos);
  return;
}

/**
 * See \ref PMMG_Set_prisms function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_PRISMS,pmmg_set_prisms,
             (PMMG_pGrp *grp, int *prisms, int *refs, int* retval),
             (grp,prisms,refs,retval)){
  *retval = PMMG_Set_prisms(*grp,prisms,refs);
  return;
}

/**
 * See \ref PMMG_Set_triangle function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_TRIANGLE,pmmg_set_triangle,
             (PMMG_pGrp *grp, int* v0, int* v1, int* v2, int* ref,int* pos,
              int* retval),
             (grp,v0,v1,v2,ref,pos,retval)) {
  *retval = PMMG_Set_triangle(*grp, *v0, *v1, *v2, *ref, *pos);
  return;
}

/**
 * See \ref PMMG_Set_triangles function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_TRIANGLES,pmmg_set_triangles,
             (PMMG_pGrp *grp, int* tria, int* refs,
              int* retval),
             (grp,tria,refs,retval)) {
  *retval = PMMG_Set_triangles(*grp, tria, refs);
  return;
}

/**
 * See \ref PMMG_Set_quadrilateral function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_QUADRILATERAL,pmmg_set_quadrilateral,
             (PMMG_pGrp *grp, int* v0, int* v1, int* v2,int *v3,
              int* ref,int* pos,int* retval),
             (grp,v0,v1,v2,v3,ref,pos,retval)) {
  *retval = PMMG_Set_quadrilateral(*grp, *v0, *v1, *v2, *v3,*ref, *pos);
  return;
}

/**
 * See \ref PMMG_Set_quadrilaterals function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_QUADRILATERALS,pmmg_set_quadrilaterals,
             (PMMG_pGrp *grp, int* quads, int* refs,
              int* retval),
             (grp,quads,refs,retval)) {
  *retval = PMMG_Set_quadrilaterals(*grp, quads, refs);
  return;
}

/**
 * See \ref PMMG_Set_corner function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_CORNER,pmmg_set_corner,(PMMG_pGrp *grp, int *k, int* retval),
             (grp,k,retval)) {
  *retval =  PMMG_Set_corner(*grp,*k);
  return;
}

/**
 * See \ref PMMG_Set_requiredVertex function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_REQUIREDVERTEX,pmmg_set_requiredvertex,
             (PMMG_pGrp *grp, int *k, int* retval),
             (grp,k,retval)) {
  *retval =  PMMG_Set_requiredVertex(*grp,*k);
  return;
}

/**
 * See \ref PMMG_Set_requiredTetrahedron function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_REQUIREDTETRAHEDRON,pmmg_set_requiredtetrahedron,
             (PMMG_pGrp *grp, int *k, int* retval),
             (grp,k,retval)) {
  *retval = PMMG_Set_requiredTetrahedron(*grp,*k);
  return;
}

/**
 * See \ref PMMG_Set_requiredTetrahedra function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_REQUIREDTETRAHEDRA,pmmg_set_requiredtetrahedra,
             (PMMG_pGrp *grp, int *reqIdx, int *nreq, int* retval),
             (grp,reqIdx,nreq,retval)) {
  *retval = PMMG_Set_requiredTetrahedra(*grp,reqIdx, *nreq);
  return;
}

/**
 * See \ref PMMG_Set_requiredTriangle function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_REQUIREDTRIANGLE,pmmg_set_requiredtriangle,
             (PMMG_pGrp *grp, int *k, int* retval),
             (grp,k,retval)) {
  *retval = PMMG_Set_requiredTriangle(*grp, *k);
  return;
}

/**
 * See \ref PMMG_Set_requiredTriangles function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_REQUIREDTRIANGLES,pmmg_set_requiredtriangles,
             (PMMG_pGrp *grp, int *reqIdx, int *nreq, int* retval),
             (grp,reqIdx,nreq,retval)) {
  *retval = PMMG_Set_requiredTriangles(*grp, reqIdx, *nreq);
  return;
}

/**
 * See \ref PMMG_Set_ridge function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_RIDGE,pmmg_set_ridge,
             (PMMG_pGrp *grp, int *k, int* retval),
             (grp,k,retval)) {
  *retval = PMMG_Set_ridge(*grp,*k);
  return;
}

/**
 * See \ref PMMG_Set_requiredEdge function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_REQUIREDEDGE,pmmg_set_requirededge,
             (PMMG_pGrp *grp, int *k, int* retval),
             (grp,k,retval)) {
  *retval = PMMG_Set_requiredEdge(*grp,*k);
  return;
}

/**
 * See \ref PMMG_Set_normalAtVertex function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_NORMALATVERTEX,pmmg_set_normalatvertex,
             (PMMG_pGrp *grp, int *k, double* n0, double* n1, double* n2,int* retval),
             (grp,k,n0,n1,n2,retval)) {
  *retval = PMMG_Set_normalAtVertex(*grp,*k, *n0, *n1, *n2);
  return;
}

/**
 * See \ref PMMG_Set_scalarSol function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_SCALARSOL,pmmg_set_scalarsol,
             (PMMG_pGrp *grp, double *s, int *pos, int* retval),
             (grp,s,pos,retval)) {
  *retval = PMMG_Set_scalarSol(*grp,*s,*pos);
  return;
}

/**
 * See \ref PMMG_Set_scalarSols function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_SCALARSOLS,pmmg_set_scalarsols,
             (PMMG_pGrp *grp, double *s, int* retval),
             (grp,s,retval)) {
  *retval = PMMG_Set_scalarSols(*grp,s);
  return;
}

/**
 * See \ref PMMG_Set_vectorSol function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_VECTORSOL,pmmg_set_vectorsol,
             (PMMG_pGrp *grp, double *vx, double *vy, double *vz,
              int *pos, int* retval),
             (grp,vx,vy,vz,pos,retval)) {
  *retval = PMMG_Set_vectorSol(*grp,*vx,*vy,*vz,*pos);
  return;
}

/**
 * See \ref PMMG_Set_vectorSols function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_VECTORSOLS,pmmg_set_vectorsols,
             (PMMG_pGrp *grp, double *sols, int* retval),
             (grp,sols,retval)) {
  *retval = PMMG_Set_vectorSols(*grp,sols);
  return;
}

/**
 * See \ref PMMG_Set_tensorSol function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_TENSORSOL,pmmg_set_tensorsol,
             (PMMG_pGrp *grp, double* m11,double *m12, double *m13,
              double* m22,double *m23, double *m33, int *pos, int* retval),
             (grp,m11,m12,m13,m22,m23,m33,pos,retval)) {
  *retval = PMMG_Set_tensorSol(*grp,*m11,*m12,*m13,*m22,*m23,*m33,*pos);
  return;
}

/**
 * See \ref PMMG_Set_tensorSols function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_TENSORSOLS,pmmg_set_tensorsols,
             (PMMG_pGrp *grp, double* sols,int* retval),
             (grp,sols,retval)) {
  *retval = PMMG_Set_tensorSols(*grp,sols);
  return;
}

/**
 * See \ref PMMG_Set_scalarMet function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_SCALARMET,pmmg_set_scalarmet,
             (PMMG_pGrp *grp, double *m, int *pos, int* retval),
             (grp,m,pos,retval)) {
  *retval = PMMG_Set_scalarMet(*grp,*m,*pos);
  return;
}

/**
 * See \ref PMMG_Set_scalaMets function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_SCALARMETS,pmmg_set_scalarmets,
             (PMMG_pGrp *grp, double *m, int* retval),
             (grp,m,retval)) {
  *retval = PMMG_Set_scalarMets(*grp,m);
  return;
}

/**
 * See \ref PMMG_Set_vectorMet function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_VECTORMET,pmmg_set_vectormet,
             (PMMG_pGrp *grp, double *vx, double *vy, double *vz,
              int *pos, int* retval),
             (grp,vx,vy,vz,pos,retval)) {
  *retval = PMMG_Set_vectorMet(*grp,*vx,*vy,*vz,*pos);
  return;
}

/**
 * See \ref PMMG_Set_vectorMets function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_VECTORMETS,pmmg_set_vectormets,
             (PMMG_pGrp *grp, double *mets, int* retval),
             (grp,mets,retval)) {
  *retval = PMMG_Set_vectorMets(*grp,mets);
  return;
}

/**
 * See \ref PMMG_Set_tensorMet function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_TENSORMET,pmmg_set_tensormet,
             (PMMG_pGrp *grp, double* m11,double *m12, double *m13,
              double* m22,double *m23, double *m33, int *pos, int* retval),
             (grp,m11,m12,m13,m22,m23,m33,pos,retval)) {
  *retval = PMMG_Set_tensorMet(*grp,*m11,*m12,*m13,*m22,*m23,*m33,*pos);
  return;
}

/**
 * See \ref PMMG_Set_tensorMets function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_TENSORMETS,pmmg_set_tensormets,
             (PMMG_pGrp *grp, double* mets,int* retval),
             (grp,mets,retval)) {
  *retval = PMMG_Set_tensorMets(*grp,mets);
  return;
}

/**
 * See \ref PMMG_Get_grpSize function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_GRPSIZE, pmmg_get_grpsize,
    (PMMG_pGrp *grp, int* np, int* ne, int* nprism, int* nt,
     int* nquad, int* na, int* typEntity, int* typSol, int* typMet,
     int* retval),
    (grp, np, ne, nprism, nt, nquad, na, typEntity, typSol, typMet,
     retval)) {
  *retval = PMMG_Get_grpSize(*grp, np, ne, nprism, nt, nquad, na,
                             typEntity, typSol, typMet);
  return;
}

/**
 * See \ref PMMG_Get_vertex function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_VERTEX,pmmg_get_vertex,
             (PMMG_pGrp *grp, double* c0, double* c1, double* c2, int* ref,
              int* isCorner, int* isRequired, int* retval),
             (grp,c0,c1,c2,ref,isCorner,isRequired, retval)) {
  *retval = PMMG_Get_vertex(*grp,c0,c1,c2,ref,isCorner,isRequired);
  return;
}

/**
 * See \ref PMMG_Get_vertices function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_VERTICES,pmmg_get_vertices,
             (PMMG_pGrp *grp, double* vertices, int* refs,
              int* areCorners, int* areRequired, int* retval),
             (grp,vertices,refs,areCorners,areRequired, retval)) {
  *retval = PMMG_Get_vertices(*grp,vertices,refs,areCorners,areRequired);
  return;
}

/**
 * See \ref PMMG_Get_tetrahedron function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_TETRAHEDRON,pmmg_get_tetrahedron,
             (PMMG_pGrp *grp, int* v0, int* v1, int* v2, int *v3,
              int* ref, int* isRequired, int* retval),
             (grp,v0,v1,v2,v3,ref,isRequired,retval)) {
  *retval = PMMG_Get_tetrahedron(*grp,v0,v1,v2,v3,ref,isRequired);
  return;
}

/**
 * See \ref PMMG_Get_tetrahedra function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_TETRAHEDRA,pmmg_get_tetrahedra,
             (PMMG_pGrp *grp, int* tetra, int* refs, int* areRequired,
              int* retval),
             (grp,tetra,refs,areRequired,retval)) {
  *retval = PMMG_Get_tetrahedra(*grp,tetra,refs,areRequired);
  return;
}

/**
 * See \ref PMMG_Get_prism function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_PRISM,pmmg_get_prism,
             (PMMG_pGrp *grp, int* v0, int* v1, int* v2, int *v3,
              int *v4, int* v5,int* ref, int* isRequired, int* retval),
             (grp,v0,v1,v2,v3,v4,v5,ref,isRequired,retval)) {
  *retval = PMMG_Get_prism(*grp,v0,v1,v2,v3,v4,v5,ref,isRequired);
  return;
}

/**
 * See \ref PMMG_Get_prisms function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_PRISMS,pmmg_get_prisms,
             (PMMG_pGrp *grp, int* prisms, int* refs, int* areRequired,
              int* retval),
             (grp,prisms,refs,areRequired,retval)) {
  *retval = PMMG_Get_prisms(*grp,prisms,refs,areRequired);
  return;
}

/**
 * See \ref PMMG_Get_triangle function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_TRIANGLE,pmmg_get_triangle,
             (PMMG_pGrp *grp, int* v0, int* v1, int* v2, int* ref
              ,int* isRequired, int* retval),
             (grp,v0,v1,v2,ref,isRequired,retval)) {
  *retval = PMMG_Get_triangle(*grp,v0,v1,v2,ref,isRequired);
  return;
}

/**
 * See \ref PMMG_Get_triangles function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_TRIANGLES,pmmg_get_triangles,
             (PMMG_pGrp *grp, int* tria, int* refs,int* areRequired,
              int* retval),
             (grp,tria,refs,areRequired,retval)) {
  *retval = PMMG_Get_triangles(*grp,tria,refs,areRequired);
  return;
}

/**
 * See \ref PMMG_Get_quadrilateral function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_QUADRILATERAL,pmmg_get_quadrilateral,
             (PMMG_pGrp *grp, int* v0, int* v1, int* v2,int *v3,
               int* ref,int* isRequired, int* retval),
             (grp,v0,v1,v2,v3,ref,isRequired,retval)) {
  *retval = PMMG_Get_quadrilateral(*grp,v0,v1,v2,v3,ref,isRequired);
  return;
}

/**
 * See \ref PMMG_Get_quadrilaterals function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_QUADRILATERALS,pmmg_get_quadrilaterals,
             (PMMG_pGrp *grp, int* quads, int* refs,int* areRequired,
              int* retval),
             (grp,quads,refs,areRequired,retval)) {
  *retval = PMMG_Get_quadrilaterals(*grp,quads,refs,areRequired);
  return;
}

/**
 * See \ref PMMG_Get_edge function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_EDGE,pmmg_get_edge,(PMMG_pGrp *grp, int* e0, int* e1, int* ref
                                          ,int* isRidge, int* isRequired, int* retval),
             (grp,e0,e1,ref,isRidge,isRequired,retval)) {
  *retval = PMMG_Get_edge(*grp,e0,e1,ref,isRidge,isRequired);
  return;
}

/**
 * See \ref PMMG_Get_normalAtVertex function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_NORMALATVERTEX,pmmg_get_normalatvertex,
             (PMMG_pGrp *grp, int *k, double* n0, double* n1, double* n2,int* retval),
             (grp,k,n0,n1,n2,retval)) {
  *retval = PMMG_Get_normalAtVertex(*grp,*k, n0, n1, n2);
  return;
}

/**
 * See \ref PMMG_Get_scalarSol function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_SCALARSOL,pmmg_get_scalarsol,
             (PMMG_pGrp *grp, double* s, int* retval),
             (grp,s,retval)) {
  *retval = PMMG_Get_scalarSol(*grp,s);
  return;
}

/**
 * See \ref PMMG_Get_scalarSols function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_SCALARSOLS,pmmg_get_scalarsols,
             (PMMG_pGrp *grp, double* s, int* retval),
             (grp,s,retval)) {
  *retval = PMMG_Get_scalarSols(*grp,s);
  return;
}

/**
 * See \ref PMMG_Get_vectorSol function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_VECTORSOL,pmmg_get_vectorsol,
             (PMMG_pGrp *grp, double* vx,double *vy, double *vz, int* retval),
             (grp,vx,vy,vz,retval)) {
  *retval = PMMG_Get_vectorSol(*grp,vx,vy,vz);
  return;
}

/**
 * See \ref PMMG_Get_vectorSols function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_VECTORSOLS,pmmg_get_vectorsols,
             (PMMG_pGrp *grp, double* sols, int* retval),
             (grp,sols,retval)) {
  *retval = PMMG_Get_vectorSols(*grp,sols);
  return;
}

/**
 * See \ref PMMG_Get_tensorSol function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_TENSORSOL,pmmg_get_tensorsol,
             (PMMG_pGrp *grp, double* m11,double *m12, double *m13,
              double* m22,double *m23, double *m33, int* retval),
             (grp,m11,m12,m13,m22,m23,m33,retval)) {
  *retval = PMMG_Get_tensorSol(*grp,m11,m12,m13,m22,m23,m33);
  return;
}

/**
 * See \ref PMMG_Get_tensorSols function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_TENSORSOLS,pmmg_get_tensorsols,
             (PMMG_pGrp *grp, double* sols, int* retval),
             (grp,sols,retval)) {
  *retval = PMMG_Get_tensorSols(*grp,sols);
  return;
}

/**
 * See \ref PMMG_Get_scalarMet function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_SCALARMET,pmmg_get_scalarmet,
             (PMMG_pGrp *grp, double* m, int* retval),
             (grp,m,retval)) {
  *retval = PMMG_Get_scalarMet(*grp,m);
  return;
}

/**
 * See \ref PMMG_Get_scalarMets function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_SCALARMETS,pmmg_get_scalarmets,
             (PMMG_pGrp *grp, double* m, int* retval),
             (grp,m,retval)) {
  *retval = PMMG_Get_scalarMets(*grp,m);
  return;
}

/**
 * See \ref PMMG_Get_vectorMet function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_VECTORMET,pmmg_get_vectormet,
             (PMMG_pGrp *grp, double* vx,double *vy, double *vz, int* retval),
             (grp,vx,vy,vz,retval)) {
  *retval = PMMG_Get_vectorMet(*grp,vx,vy,vz);
  return;
}

/**
 * See \ref PMMG_Get_vectorMets function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_VECTORMETS,pmmg_get_vectormets,
             (PMMG_pGrp *grp, double* mets, int* retval),
             (grp,mets,retval)) {
  *retval = PMMG_Get_vectorMets(*grp,mets);
  return;
}

/**
 * See \ref PMMG_Get_tensorMet function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_TENSORMET,pmmg_get_tensormet,
             (PMMG_pGrp *grp, double* m11,double *m12, double *m13,
              double* m22,double *m23, double *m33, int* retval),
             (grp,m11,m12,m13,m22,m23,m33,retval)) {
  *retval = PMMG_Get_tensorMet(*grp,m11,m12,m13,m22,m23,m33);
  return;
}

/**
 * See \ref PMMG_Get_tensorMets function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_TENSORMETS,pmmg_get_tensormets,
             (PMMG_pGrp *grp, double* mets, int* retval),
             (grp,mets,retval)) {
  *retval = PMMG_Get_tensorMets(*grp,mets);
  return;
}

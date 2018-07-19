/**
 * \file API_functionsf_pmmg.c
 * \brief Fortran API functions for PARMMG library.
 *
 * Define the Fortran API functions for PARMMG library: adds function
 * definitions with upcase, underscore and double underscore to match
 * any fortran compiler.
 *
 */
#include "parmmg.h"

/**
 * See \ref PMMG_Init_parmesh function in \ref libparmmg.h file.
 */
FORTRAN_VARIADIC ( PMMG_INIT_PARMESH, pmmg_init_parmesh,
                   (const int starter, ... ),
                   va_list argptr;
                   int     ier;

                   va_start(argptr, starter);

                   ier = PMMG_Init_parMesh_var(argptr);

                   va_end(argptr);

                   if ( !ier ) exit(EXIT_FAILURE);

                   return;
  )

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
 * See \ref PMMG_Init_parameters function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_INIT_PARAMETERS,pmmg_init_parameters,
             (PMMG_pParMesh *parmesh,MPI_Comm *comm),
             (parmesh,comm)) {

  PMMG_Init_parameters(*parmesh,*comm);

  return;
}

/**
 * See \ref PMMG_Set_meshSize function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_MESHSIZE, pmmg_set_meshsize,
    (PMMG_pParMesh *parmesh, int* np, int* ne, int* nprism, int* nt,
     int* nquad, int* na, int* retval),
    (parmesh, np, ne, nprism, nt, nquad, na, retval)) {

  *retval = PMMG_Set_meshSize(*parmesh, *np, *ne, *nprism, *nt, *nquad,*na);

  return;
}

/**
 * See \ref PMMG_Set_allSolsSizes function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_ALLSOLSSIZES, pmmg_set_allsolssizes,
    (PMMG_pParMesh *parmesh,int *nsol,int* typEntity,int *np,int* typSol,
     int* retval),
    (parmesh, nsol, typEntity, np, typSol, retval)) {
  *retval = PMMG_Set_allSolsSizes(*parmesh,*nsol,typEntity,*np,typSol);
  return;
}

/**
 * See \ref PMMG_Set_metSize function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_METSIZE, pmmg_set_metsize,
    (PMMG_pParMesh *parmesh,int* typEntity,int *np,int* typMet,
     int* retval),
    (parmesh, typEntity, np, typMet, retval)) {
  *retval = PMMG_Set_metSize(*parmesh,*typEntity,*np,*typMet);
  return;
}

/**
 * See \ref PMMG_Set_vertex function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_VERTEX,pmmg_set_vertex,
             (PMMG_pParMesh *parmesh, double* c0, double* c1, double* c2, int* ref,
              int* pos, int* retval),
             (parmesh,c0,c1,c2,ref,pos,retval)) {

  *retval = PMMG_Set_vertex(*parmesh,*c0,*c1,*c2,*ref,*pos);
  return;
}

/**
 * See \ref PMMG_Set_vertices function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_VERTICES,pmmg_set_vertices,
             (PMMG_pParMesh *parmesh, double* vertices, int* refs, int* retval),
             (parmesh,vertices,refs,retval)) {

  *retval = PMMG_Set_vertices(*parmesh,vertices,refs);
  return;
}

/**
 * See \ref PMMG_Set_tetrahedron function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_TETRAHEDRON,pmmg_set_tetrahedron,
             (PMMG_pParMesh *parmesh, int *v0, int *v1, int *v2, int *v3, int *ref,
              int *pos, int* retval),
             (parmesh,v0,v1,v2,v3,ref,pos,retval)){
  *retval = PMMG_Set_tetrahedron(*parmesh,*v0,*v1,*v2,*v3,*ref,*pos);
  return;
}

/**
 * See \ref PMMG_Set_tetrahedra function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_TETRAHEDRA,pmmg_set_tetrahedra,
             (PMMG_pParMesh *parmesh, int *tetra, int *refs, int* retval),
             (parmesh,tetra,refs,retval)){
  *retval = PMMG_Set_tetrahedra(*parmesh,tetra,refs);
  return;
}

/**
 * See \ref PMMG_Set_prism function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_PRISM,pmmg_set_prism,
             (PMMG_pParMesh *parmesh, int *v0, int *v1, int *v2, int *v3,
              int *v4, int *v5,int *ref,int *pos, int* retval),
             (parmesh,v0,v1,v2,v3,v4,v5,ref,pos,retval)){
  *retval = PMMG_Set_prism(*parmesh,*v0,*v1,*v2,*v3,*v4,*v5,*ref,*pos);
  return;
}

/**
 * See \ref PMMG_Set_prisms function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_PRISMS,pmmg_set_prisms,
             (PMMG_pParMesh *parmesh, int *prisms, int *refs, int* retval),
             (parmesh,prisms,refs,retval)){
  *retval = PMMG_Set_prisms(*parmesh,prisms,refs);
  return;
}

/**
 * See \ref PMMG_Set_triangle function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_TRIANGLE,pmmg_set_triangle,
             (PMMG_pParMesh *parmesh, int* v0, int* v1, int* v2, int* ref,int* pos,
              int* retval),
             (parmesh,v0,v1,v2,ref,pos,retval)) {
  *retval = PMMG_Set_triangle(*parmesh, *v0, *v1, *v2, *ref, *pos);
  return;
}

/**
 * See \ref PMMG_Set_triangles function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_TRIANGLES,pmmg_set_triangles,
             (PMMG_pParMesh *parmesh, int* tria, int* refs,
              int* retval),
             (parmesh,tria,refs,retval)) {
  *retval = PMMG_Set_triangles(*parmesh, tria, refs);
  return;
}

/**
 * See \ref PMMG_Set_quadrilateral function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_QUADRILATERAL,pmmg_set_quadrilateral,
             (PMMG_pParMesh *parmesh, int* v0, int* v1, int* v2,int *v3,
              int* ref,int* pos,int* retval),
             (parmesh,v0,v1,v2,v3,ref,pos,retval)) {
  *retval = PMMG_Set_quadrilateral(*parmesh, *v0, *v1, *v2, *v3,*ref, *pos);
  return;
}

/**
 * See \ref PMMG_Set_quadrilaterals function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_QUADRILATERALS,pmmg_set_quadrilaterals,
             (PMMG_pParMesh *parmesh, int* quads, int* refs,
              int* retval),
             (parmesh,quads,refs,retval)) {
  *retval = PMMG_Set_quadrilaterals(*parmesh, quads, refs);
  return;
}

/**
 * See \ref PMMG_Set_corner function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_CORNER,pmmg_set_corner,(PMMG_pParMesh *parmesh, int *k, int* retval),
             (parmesh,k,retval)) {
  *retval =  PMMG_Set_corner(*parmesh,*k);
  return;
}

/**
 * See \ref PMMG_Set_requiredVertex function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_REQUIREDVERTEX,pmmg_set_requiredvertex,
             (PMMG_pParMesh *parmesh, int *k, int* retval),
             (parmesh,k,retval)) {
  *retval =  PMMG_Set_requiredVertex(*parmesh,*k);
  return;
}

/**
 * See \ref PMMG_Set_requiredTetrahedron function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_REQUIREDTETRAHEDRON,pmmg_set_requiredtetrahedron,
             (PMMG_pParMesh *parmesh, int *k, int* retval),
             (parmesh,k,retval)) {
  *retval = PMMG_Set_requiredTetrahedron(*parmesh,*k);
  return;
}

/**
 * See \ref PMMG_Set_requiredTetrahedra function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_REQUIREDTETRAHEDRA,pmmg_set_requiredtetrahedra,
             (PMMG_pParMesh *parmesh, int *reqIdx, int *nreq, int* retval),
             (parmesh,reqIdx,nreq,retval)) {
  *retval = PMMG_Set_requiredTetrahedra(*parmesh,reqIdx, *nreq);
  return;
}

/**
 * See \ref PMMG_Set_requiredTriangle function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_REQUIREDTRIANGLE,pmmg_set_requiredtriangle,
             (PMMG_pParMesh *parmesh, int *k, int* retval),
             (parmesh,k,retval)) {
  *retval = PMMG_Set_requiredTriangle(*parmesh, *k);
  return;
}

/**
 * See \ref PMMG_Set_requiredTriangles function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_REQUIREDTRIANGLES,pmmg_set_requiredtriangles,
             (PMMG_pParMesh *parmesh, int *reqIdx, int *nreq, int* retval),
             (parmesh,reqIdx,nreq,retval)) {
  *retval = PMMG_Set_requiredTriangles(*parmesh, reqIdx, *nreq);
  return;
}

/**
 * See \ref PMMG_Set_ridge function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_RIDGE,pmmg_set_ridge,
             (PMMG_pParMesh *parmesh, int *k, int* retval),
             (parmesh,k,retval)) {
  *retval = PMMG_Set_ridge(*parmesh,*k);
  return;
}

/**
 * See \ref PMMG_Set_requiredEdge function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_REQUIREDEDGE,pmmg_set_requirededge,
             (PMMG_pParMesh *parmesh, int *k, int* retval),
             (parmesh,k,retval)) {
  *retval = PMMG_Set_requiredEdge(*parmesh,*k);
  return;
}

/**
 * See \ref PMMG_Set_normalAtVertex function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_NORMALATVERTEX,pmmg_set_normalatvertex,
             (PMMG_pParMesh *parmesh, int *k, double* n0, double* n1, double* n2,int* retval),
             (parmesh,k,n0,n1,n2,retval)) {
  *retval = PMMG_Set_normalAtVertex(*parmesh,*k, *n0, *n1, *n2);
  return;
}

/**
 * See \ref PMMG_Set_scalarSol function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_SCALARSOL,pmmg_set_scalarsol,
             (PMMG_pParMesh *parmesh, double *s, int *pos, int* retval),
             (parmesh,s,pos,retval)) {
  *retval = PMMG_Set_scalarSol(*parmesh,*s,*pos);
  return;
}

/**
 * See \ref PMMG_Set_scalarSols function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_SCALARSOLS,pmmg_set_scalarsols,
             (PMMG_pParMesh *parmesh, double *s, int* retval),
             (parmesh,s,retval)) {
  *retval = PMMG_Set_scalarSols(*parmesh,s);
  return;
}

/**
 * See \ref PMMG_Set_vectorSol function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_VECTORSOL,pmmg_set_vectorsol,
             (PMMG_pParMesh *parmesh, double *vx, double *vy, double *vz,
              int *pos, int* retval),
             (parmesh,vx,vy,vz,pos,retval)) {
  *retval = PMMG_Set_vectorSol(*parmesh,*vx,*vy,*vz,*pos);
  return;
}

/**
 * See \ref PMMG_Set_vectorSols function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_VECTORSOLS,pmmg_set_vectorsols,
             (PMMG_pParMesh *parmesh, double *sols, int* retval),
             (parmesh,sols,retval)) {
  *retval = PMMG_Set_vectorSols(*parmesh,sols);
  return;
}

/**
 * See \ref PMMG_Set_tensorSol function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_TENSORSOL,pmmg_set_tensorsol,
             (PMMG_pParMesh *parmesh, double* m11,double *m12, double *m13,
              double* m22,double *m23, double *m33, int *pos, int* retval),
             (parmesh,m11,m12,m13,m22,m23,m33,pos,retval)) {
  *retval = PMMG_Set_tensorSol(*parmesh,*m11,*m12,*m13,*m22,*m23,*m33,*pos);
  return;
}

/**
 * See \ref PMMG_Set_tensorSols function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_TENSORSOLS,pmmg_set_tensorsols,
             (PMMG_pParMesh *parmesh, double* sols,int* retval),
             (parmesh,sols,retval)) {
  *retval = PMMG_Set_tensorSols(*parmesh,sols);
  return;
}

/**
 * See \ref PMMG_Set_scalarMet function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_SCALARMET,pmmg_set_scalarmet,
             (PMMG_pParMesh *parmesh, double *m, int *pos, int* retval),
             (parmesh,m,pos,retval)) {
  *retval = PMMG_Set_scalarMet(*parmesh,*m,*pos);
  return;
}

/**
 * See \ref PMMG_Set_scalaMets function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_SCALARMETS,pmmg_set_scalarmets,
             (PMMG_pParMesh *parmesh, double *m, int* retval),
             (parmesh,m,retval)) {
  *retval = PMMG_Set_scalarMets(*parmesh,m);
  return;
}

/**
 * See \ref PMMG_Set_vectorMet function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_VECTORMET,pmmg_set_vectormet,
             (PMMG_pParMesh *parmesh, double *vx, double *vy, double *vz,
              int *pos, int* retval),
             (parmesh,vx,vy,vz,pos,retval)) {
  *retval = PMMG_Set_vectorMet(*parmesh,*vx,*vy,*vz,*pos);
  return;
}

/**
 * See \ref PMMG_Set_vectorMets function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_VECTORMETS,pmmg_set_vectormets,
             (PMMG_pParMesh *parmesh, double *mets, int* retval),
             (parmesh,mets,retval)) {
  *retval = PMMG_Set_vectorMets(*parmesh,mets);
  return;
}

/**
 * See \ref PMMG_Set_tensorMet function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_TENSORMET,pmmg_set_tensormet,
             (PMMG_pParMesh *parmesh, double* m11,double *m12, double *m13,
              double* m22,double *m23, double *m33, int *pos, int* retval),
             (parmesh,m11,m12,m13,m22,m23,m33,pos,retval)) {
  *retval = PMMG_Set_tensorMet(*parmesh,*m11,*m12,*m13,*m22,*m23,*m33,*pos);
  return;
}

/**
 * See \ref PMMG_Set_tensorMets function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_SET_TENSORMETS,pmmg_set_tensormets,
             (PMMG_pParMesh *parmesh, double* mets,int* retval),
             (parmesh,mets,retval)) {
  *retval = PMMG_Set_tensorMets(*parmesh,mets);
  return;
}

/**
 * See \ref PMMG_Get_meshSize function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_MESHSIZE, pmmg_get_meshsize,
    (PMMG_pParMesh *parmesh, int* np, int* ne, int* nprism, int* nt,
     int* nquad, int* na, int* retval),
    (parmesh, np, ne, nprism, nt, nquad, na, retval)) {

  *retval = PMMG_Get_meshSize(*parmesh, np, ne, nprism, nt, nquad, na);

  return;
}

/**
 * See \ref PMMG_Get_allSolsSizes function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_ALLSOLSSIZES, pmmg_get_allsolssizes,
    (PMMG_pParMesh *parmesh,int *nsol,int* typEntity,int *np,int* typSol,
     int* retval),
    (parmesh, nsol, typEntity, np, typSol, retval)) {
  *retval = PMMG_Get_allSolsSizes(*parmesh,nsol,typEntity,np,typSol);
  return;
}

/**
 * See \ref PMMG_Get_metSize function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_METSIZE, pmmg_get_metsize,
    (PMMG_pParMesh *parmesh,int* typEntity,int *np,int* typMet,
     int* retval),
    (parmesh, typEntity, np, typMet, retval)) {
  *retval = PMMG_Get_metSize(*parmesh,typEntity,np,typMet);
  return;
}

/**
 * See \ref PMMG_Get_vertex function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_VERTEX,pmmg_get_vertex,
             (PMMG_pParMesh *parmesh, double* c0, double* c1, double* c2, int* ref,
              int* isCorner, int* isRequired, int* retval),
             (parmesh,c0,c1,c2,ref,isCorner,isRequired, retval)) {
  *retval = PMMG_Get_vertex(*parmesh,c0,c1,c2,ref,isCorner,isRequired);
  return;
}

/**
 * See \ref PMMG_Get_vertices function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_VERTICES,pmmg_get_vertices,
             (PMMG_pParMesh *parmesh, double* vertices, int* refs,
              int* areCorners, int* areRequired, int* retval),
             (parmesh,vertices,refs,areCorners,areRequired, retval)) {
  *retval = PMMG_Get_vertices(*parmesh,vertices,refs,areCorners,areRequired);
  return;
}

/**
 * See \ref PMMG_Get_tetrahedron function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_TETRAHEDRON,pmmg_get_tetrahedron,
             (PMMG_pParMesh *parmesh, int* v0, int* v1, int* v2, int *v3,
              int* ref, int* isRequired, int* retval),
             (parmesh,v0,v1,v2,v3,ref,isRequired,retval)) {
  *retval = PMMG_Get_tetrahedron(*parmesh,v0,v1,v2,v3,ref,isRequired);
  return;
}

/**
 * See \ref PMMG_Get_tetrahedra function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_TETRAHEDRA,pmmg_get_tetrahedra,
             (PMMG_pParMesh *parmesh, int* tetra, int* refs, int* areRequired,
              int* retval),
             (parmesh,tetra,refs,areRequired,retval)) {
  *retval = PMMG_Get_tetrahedra(*parmesh,tetra,refs,areRequired);
  return;
}

/**
 * See \ref PMMG_Get_prism function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_PRISM,pmmg_get_prism,
             (PMMG_pParMesh *parmesh, int* v0, int* v1, int* v2, int *v3,
              int *v4, int* v5,int* ref, int* isRequired, int* retval),
             (parmesh,v0,v1,v2,v3,v4,v5,ref,isRequired,retval)) {
  *retval = PMMG_Get_prism(*parmesh,v0,v1,v2,v3,v4,v5,ref,isRequired);
  return;
}

/**
 * See \ref PMMG_Get_prisms function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_PRISMS,pmmg_get_prisms,
             (PMMG_pParMesh *parmesh, int* prisms, int* refs, int* areRequired,
              int* retval),
             (parmesh,prisms,refs,areRequired,retval)) {
  *retval = PMMG_Get_prisms(*parmesh,prisms,refs,areRequired);
  return;
}

/**
 * See \ref PMMG_Get_triangle function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_TRIANGLE,pmmg_get_triangle,
             (PMMG_pParMesh *parmesh, int* v0, int* v1, int* v2, int* ref
              ,int* isRequired, int* retval),
             (parmesh,v0,v1,v2,ref,isRequired,retval)) {
  *retval = PMMG_Get_triangle(*parmesh,v0,v1,v2,ref,isRequired);
  return;
}

/**
 * See \ref PMMG_Get_triangles function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_TRIANGLES,pmmg_get_triangles,
             (PMMG_pParMesh *parmesh, int* tria, int* refs,int* areRequired,
              int* retval),
             (parmesh,tria,refs,areRequired,retval)) {
  *retval = PMMG_Get_triangles(*parmesh,tria,refs,areRequired);
  return;
}

/**
 * See \ref PMMG_Get_quadrilateral function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_QUADRILATERAL,pmmg_get_quadrilateral,
             (PMMG_pParMesh *parmesh, int* v0, int* v1, int* v2,int *v3,
               int* ref,int* isRequired, int* retval),
             (parmesh,v0,v1,v2,v3,ref,isRequired,retval)) {
  *retval = PMMG_Get_quadrilateral(*parmesh,v0,v1,v2,v3,ref,isRequired);
  return;
}

/**
 * See \ref PMMG_Get_quadrilaterals function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_QUADRILATERALS,pmmg_get_quadrilaterals,
             (PMMG_pParMesh *parmesh, int* quads, int* refs,int* areRequired,
              int* retval),
             (parmesh,quads,refs,areRequired,retval)) {
  *retval = PMMG_Get_quadrilaterals(*parmesh,quads,refs,areRequired);
  return;
}

/**
 * See \ref PMMG_Get_edge function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_EDGE,pmmg_get_edge,(PMMG_pParMesh *parmesh, int* e0, int* e1, int* ref
                                          ,int* isRidge, int* isRequired, int* retval),
             (parmesh,e0,e1,ref,isRidge,isRequired,retval)) {
  *retval = PMMG_Get_edge(*parmesh,e0,e1,ref,isRidge,isRequired);
  return;
}

/**
 * See \ref PMMG_Get_normalAtVertex function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_NORMALATVERTEX,pmmg_get_normalatvertex,
             (PMMG_pParMesh *parmesh, int *k, double* n0, double* n1, double* n2,int* retval),
             (parmesh,k,n0,n1,n2,retval)) {
  *retval = PMMG_Get_normalAtVertex(*parmesh,*k, n0, n1, n2);
  return;
}

/**
 * See \ref PMMG_Get_scalarSol function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_SCALARSOL,pmmg_get_scalarsol,
             (PMMG_pParMesh *parmesh, double* s, int* retval),
             (parmesh,s,retval)) {
  *retval = PMMG_Get_scalarSol(*parmesh,s);
  return;
}

/**
 * See \ref PMMG_Get_scalarSols function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_SCALARSOLS,pmmg_get_scalarsols,
             (PMMG_pParMesh *parmesh, double* s, int* retval),
             (parmesh,s,retval)) {
  *retval = PMMG_Get_scalarSols(*parmesh,s);
  return;
}

/**
 * See \ref PMMG_Get_vectorSol function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_VECTORSOL,pmmg_get_vectorsol,
             (PMMG_pParMesh *parmesh, double* vx,double *vy, double *vz, int* retval),
             (parmesh,vx,vy,vz,retval)) {
  *retval = PMMG_Get_vectorSol(*parmesh,vx,vy,vz);
  return;
}

/**
 * See \ref PMMG_Get_vectorSols function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_VECTORSOLS,pmmg_get_vectorsols,
             (PMMG_pParMesh *parmesh, double* sols, int* retval),
             (parmesh,sols,retval)) {
  *retval = PMMG_Get_vectorSols(*parmesh,sols);
  return;
}

/**
 * See \ref PMMG_Get_tensorSol function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_TENSORSOL,pmmg_get_tensorsol,
             (PMMG_pParMesh *parmesh, double* m11,double *m12, double *m13,
              double* m22,double *m23, double *m33, int* retval),
             (parmesh,m11,m12,m13,m22,m23,m33,retval)) {
  *retval = PMMG_Get_tensorSol(*parmesh,m11,m12,m13,m22,m23,m33);
  return;
}

/**
 * See \ref PMMG_Get_tensorSols function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_TENSORSOLS,pmmg_get_tensorsols,
             (PMMG_pParMesh *parmesh, double* sols, int* retval),
             (parmesh,sols,retval)) {
  *retval = PMMG_Get_tensorSols(*parmesh,sols);
  return;
}

/**
 * See \ref PMMG_Get_scalarMet function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_SCALARMET,pmmg_get_scalarmet,
             (PMMG_pParMesh *parmesh, double* m, int* retval),
             (parmesh,m,retval)) {
  *retval = PMMG_Get_scalarMet(*parmesh,m);
  return;
}

/**
 * See \ref PMMG_Get_scalarMets function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_SCALARMETS,pmmg_get_scalarmets,
             (PMMG_pParMesh *parmesh, double* m, int* retval),
             (parmesh,m,retval)) {
  *retval = PMMG_Get_scalarMets(*parmesh,m);
  return;
}

/**
 * See \ref PMMG_Get_vectorMet function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_VECTORMET,pmmg_get_vectormet,
             (PMMG_pParMesh *parmesh, double* vx,double *vy, double *vz, int* retval),
             (parmesh,vx,vy,vz,retval)) {
  *retval = PMMG_Get_vectorMet(*parmesh,vx,vy,vz);
  return;
}

/**
 * See \ref PMMG_Get_vectorMets function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_VECTORMETS,pmmg_get_vectormets,
             (PMMG_pParMesh *parmesh, double* mets, int* retval),
             (parmesh,mets,retval)) {
  *retval = PMMG_Get_vectorMets(*parmesh,mets);
  return;
}

/**
 * See \ref PMMG_Get_tensorMet function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_TENSORMET,pmmg_get_tensormet,
             (PMMG_pParMesh *parmesh, double* m11,double *m12, double *m13,
              double* m22,double *m23, double *m33, int* retval),
             (parmesh,m11,m12,m13,m22,m23,m33,retval)) {
  *retval = PMMG_Get_tensorMet(*parmesh,m11,m12,m13,m22,m23,m33);
  return;
}

/**
 * See \ref PMMG_Get_tensorMets function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_GET_TENSORMETS,pmmg_get_tensormets,
             (PMMG_pParMesh *parmesh, double* mets, int* retval),
             (parmesh,mets,retval)) {
  *retval = PMMG_Get_tensorMets(*parmesh,mets);
  return;
}

/**
 * See \ref PMMG_Free_all function in \ref mmg3d/libmmg3d.h file.
 */
FORTRAN_VARIADIC(PMMG_FREE_ALL,pmmg_free_all,
                 (const int starter,...),
                 va_list argptr;
                 int     ier;

                 va_start(argptr, starter);

                 ier = PMMG_Free_all_var(argptr);

                 va_end(argptr);

                 if ( !ier ) exit(EXIT_FAILURE);
                 return;
  )


/**
 * See \ref PMMG_parmmglib_distributed function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_PARMMGLIB_DISTRIBUTED,pmmg_parmmglib_distributed,
             (PMMG_pParMesh *parmesh,int* retval),
             (parmesh,retval)) {
  *retval = PMMG_parmmglib_distributed(*parmesh);
  return;
}

/**
 * See \ref PMMG_parmmglib_centralized function in \ref libparmmg.h file.
 */
FORTRAN_NAME(PMMG_PARMMGLIB_CENTRALIZED,pmmg_parmmglib_centralized,
             (PMMG_pParMesh *parmesh,int* retval),
             (parmesh,retval)) {
  *retval = PMMG_parmmglib_centralized(*parmesh);
  return;
}

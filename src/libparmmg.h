#include "libparmmgtypes.h"
#include "mmg3d.h"
#include "libmmg3d.h"

#include "mpi.h"



#define PMMG_VER   "1.0.0"
#define PMMG_REL   "2016"
#define PMMG_CPY   "Copyright (c) Bx INP/INRIA, 2016-"
#define PMMG_STR   "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"


/*inout_3d.c*/
int PMMG_saveMesh(PMMG_pParMesh ,const char *);
int PMMG_loadMesh(PMMG_pParMesh ,const char *);

/*metisfunctions.c*/
int PMMG_metispartitioning(PMMG_pParMesh ,int *);

/*distributemesh*/
int PMMG_distributeMesh(PMMG_pParMesh ,int *);

/**
 * \param grp       Pointer towards the group structure.
 * \param np        Number of vertices.
 * \param ne        Number of tetrahedra.
 * \param nprism    Number of prisms.
 * \param nt        Number of triangles.
 * \param nquad     Number of quads.
 * \param na        Number of edges.
 * \param typEntity Type of solution/metrics entities (vertices, triangles...).
 * \param typSol    Type of solution (scalar, vectorial...).
 * \param typMet    Type of metrics (scalar, vectorial...).
 * \return          0 if failed, 1 otherwise.
 *
 * Set the group mesh, solution, metrics size by calling MMG3D functions.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_GRPSIZE(grp, np, ne, nprism,n t, nquad, na, &\n
 * >                               typEntity, typSol, typMet, retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     INTEGER, INTENT(IN)           :: np, ne, nprism, nt, nquad, na, &\n
 * >                                      typEntity, typSol, typMet\n
 * >     INTEGER, INTENT(INOUT)        :: retval\n
 *
 */
int PMMG_Set_grpSize(PMMG_pGrp grp, int np, int ne, int nprism, int nt,
                     int nquad, int na, int typEntity, int typSol, int typMet);

/**
 * \param grp  pointer toward the group structure.
 * \param c0   coordinate of the point along the first dimension.
 * \param c1   coordinate of the point along the second dimension.
 * \param c2   coordinate of the point along the third dimension.
 * \param ref  point reference.
 * \param pos  position of the point in the mesh.
 * \return 1.
 *
 * Set vertex of coordinates \a c0, \a c1,\a c2 and reference \a ref
 * at position \a pos in mesh structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_VERTEX(grp,c0,c1,c2,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     REAL(KIND=8), INTENT(IN)      :: c0,c1,c2\n
 * >     INTEGER, INTENT(IN)           :: ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_vertex(PMMG_pGrp grp, double c0, double c1, double c2,
                    int ref, int pos);

/**
 * \param grp      pointer toward the group structure.
 * \param vertices table of the points coor.
 * The coordinates of the \f$i^{th}\f$ point are stored in
 * vertices[(i-1)*3]\@3.
 * \param refs     table of points references.
 * The ref of the \f$i^th\f$ point is stored in refs[i-1].
 * \return 1.
 *
 * Set vertices coordinates and references in mesh structure
 * (wrapper for MMG3D function).
 *
 * \remark Fortran interface: (commentated in order to allow to pass
 * \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE PMMG_SET_VERTICES(grp,vertices,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * > !    REAL(KIND=8), INTENT(IN)      :: vertices(*)\n
 * > !    INTEGER,INTENT(IN)            :: refs(*)\n
 * > !    INTEGER, INTENT(OUT)          :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int PMMG_Set_vertices(PMMG_pGrp grp, double *vertices, int *refs);

/**
 * \param grp pointer toward the group structure.
 * \param v0  first vertex of tetrahedron.
 * \param v1  second vertex of tetrahedron.
 * \param v2  third vertex of tetrahedron.
 * \param v3  fourth vertex of tetrahedron.
 * \param ref tetrahedron reference.
 * \param pos tetrahedron position in the mesh.
 * \return 0  if failed, 1 otherwise.
 *
 * Set tetrahedra of vertices \a v0, \a v1,\a v2,\a v3 and reference
 * \a ref at position \a pos in mesh structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_TETRAHEDRON(grp,v0,v1,v2,v3,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     INTEGER, INTENT(IN)           :: v0,v1,v2,v3,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_tetrahedron(PMMG_pGrp grp, int v0, int v1, int v2, int v3,
                         int ref, int pos);

/**
 * \param grp   pointer toward the group structure.
 * \param tetra vertices of the tetras of the mesh
 * Vertices of the \f$i^{th}\f$ tetra are stored in tetra[(i-1)*4]\@4.
 * \param refs  table of the tetrahedra references.
 * References of the \f$i^{th}\f$ tetra is stored in refs[i-1].
 * \return 0    if failed, 1 otherwise.
 *
 * Set vertices and references of the mesh tetrahedra.
 *
 * \remark Fortran interface: (commentated in
 * order to allow to pass \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE PMMG_SET_TETRAHEDRA(grp,tetra,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)     :: grp\n
 * > !    INTEGER, DIMENSION(*), INTENT(IN) :: tetra,refs\n
 * > !    INTEGER, INTENT(OUT)              :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int PMMG_Set_tetrahedra(PMMG_pGrp grp, int *tetra, int *refs);

/**
 * \param grp pointer toward the group structure.
 * \param v0  first vertex of prism.
 * \param v1  second vertex of prism.
 * \param v2  third vertex of prism.
 * \param v3  fourth vertex of prism.
 * \param v4  fifth vertex of prism.
 * \param v5  sixth vertex of prism.
 * \param ref prism reference.
 * \param pos prism position in the mesh.
 * \return 0  if failed, 1 otherwise.
 *
 * Set prisms of vertices \a v0, \a v1,\a v2,\a v3,\a v4,\a v5 and reference
 * \a ref at position \a pos in mesh structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_PRISM(grp,v0,v1,v2,v3,v4,v5,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     INTEGER, INTENT(IN)           :: v0,v1,v2,v3,v4,v5,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_prism(PMMG_pGrp grp, int v0, int v1, int v2,
                   int v3, int v4, int v5, int ref, int pos);

/**
 * \param grp    pointer toward the mesh structure.
 * \param prisms vertices of the prisms of the mesh
 * Vertices of the \f$i^{th}\f$ prism are stored in prism[(i-1)*6]\@6.
 * \param refs   table of the prisms references.
 * References of the \f$i^{th}\f$ prisms is stored in refs[i-1].
 * \return 0     if failed, 1 otherwise.
 *
 * Set vertices and references of the mesh prisms (wrapper for MMG3D
 * function).
 *
 * \remark Fortran interface: (commentated in
 * order to allow to pass \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE PMMG_SET_PRISMS(grp,prisms,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)     :: grp\n
 * > !    INTEGER, DIMENSION(*), INTENT(IN) :: prisms,refs\n
 * > !    INTEGER, INTENT(OUT)              :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int PMMG_Set_prisms(PMMG_pGrp grp, int *prisms, int *refs);

/**
 * \param grp pointer toward the group structure.
 * \param v0  first vertex of triangle.
 * \param v1  second vertex of triangle.
 * \param v2  third vertex of triangle.
 * \param ref triangle reference.
 * \param pos triangle position in the mesh.
 * \return 0  if failed, 1 otherwise.
 *
 * Set triangle of vertices \a v0, \a v1, \a v2 and reference \a ref
 * at position \a pos in mesh structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_TRIANGLE(grp,v0,v1,v2,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     INTEGER, INTENT(IN)           :: v0,v1,v2,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_triangle(PMMG_pGrp grp, int v0, int v1, int v2,
                      int ref,int pos);

/**
 * \param grp  pointer toward the group structure.
 * \param tria pointer toward the table of the tria vertices
 * Vertices of the \f$i^{th}\f$ tria are stored in tria[(i-1)*3]\@3.
 * \param refs pointer toward the table of the triangle references.
 * refs[i-1] is the ref of the \f$i^{th}\f$ tria.
 * \return 0   if failed, 1 otherwise.
 *
 * Set vertices and references of the mesh triangles (wrapper for MMG3D
 * function).
 *
 * \remark Fortran interface: (commentated in
 * order to allow to pass \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE PMMG_SET_TRIANGLES(grp,tria,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)    :: grp\n
 * > !    INTEGER,DIMENSION(*), INTENT(IN) :: tria,refs\n
 * > !    INTEGER, INTENT(OUT)             :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int PMMG_Set_triangles(PMMG_pGrp grp, int *tria, int *refs);

/**
 * \param grp pointer toward the group structure.
 * \param v0  first vertex of quadrilateral.
 * \param v1  second vertex of quadrilateral.
 * \param v2  third vertex of quadrilateral.
 * \param v3  fourth vertex of quadrilateral.
 * \param ref quadrilateral reference.
 * \param pos quadrilateral position in the mesh.
 * \return 0  if failed, 1 otherwise.
 *
 * Set quadrilateral of vertices \a v0, \a v1, \a v2, \a v3 and reference \a ref
 * at position \a pos in mesh structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_QUADRILATERAL(grp,v0,v1,v2,v3,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     INTEGER, INTENT(IN)           :: v0,v1,v2,v3,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_quadrilateral(PMMG_pGrp grp, int v0, int v1, int v2, int v3,
                           int ref,int pos);

/**
 * \param grp   pointer toward the group structure.
 * \param quads pointer toward the table of the quads vertices
 * Vertices of the \f$i^{th}\f$ quad are stored in quads[(i-1)*3]\@3.
 * \param refs  pointer toward the table of the quadrilateral references.
 * refs[i-1] is the ref of the \f$i^{th}\f$ quad.
 * \return 0    if failed, 1 otherwise.
 *
 * Set vertices and references of the mesh quadrilaterals
 * (wrapper for MMG3D function).
 *
 * \remark Fortran interface: (commentated in
 * order to allow to pass \%val(0) instead of the refs array)
 *
 * > !  SUBROUTINE PMMG_SET_QUADRILATERALS(grp,quads,refs,retval)\n
 * > !    MMG5_DATA_PTR_T,INTENT(INOUT)    :: grp\n
 * > !    INTEGER,DIMENSION(*), INTENT(IN) :: quads,refs\n
 * > !    INTEGER, INTENT(OUT)             :: retval\n
 * > !  END SUBROUTINE\n
 *
 */
int PMMG_Set_quadrilaterals(PMMG_pGrp grp, int *quads, int *refs);

/**
 * \param grp pointer toward the grp structure.
 * \param v0  first extremity of the edge.
 * \param v1  second extremity of the edge.
 * \param ref edge reference.
 * \param pos edge position in the mesh.
 * \return 0  if failed, 1 otherwise.
 *
 * Set edges of extremities \a v0, \a v1 and reference \a ref at
 * position \a pos in mesh structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_EDGE(grp,v0,v1,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     INTEGER, INTENT(IN)           :: v0,v1,ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_edge(PMMG_pGrp grp, int v0, int v1, int ref, int pos);

/**
 * \param grp pointer toward the group structure.
 * \param k   vertex index.
 * \return 1.
 *
 * Set corner at point \a pos.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_CORNER(grp,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_corner(PMMG_pGrp grp, int k);

/**
 * \param grp pointer toward the group structure.
 * \param k   vertex index.
 * \return 1.
 *
 * Set point \a k as required (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_REQUIREDVERTEX(grp,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_requiredVertex(PMMG_pGrp grp, int k);

/**
 * \param grp pointer toward the group structure.
 * \param k   element index.
 * \return 1.
 *
 * Set element \a k as required (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_REQUIREDTETRAHEDRON(grp,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_requiredTetrahedron(PMMG_pGrp grp, int k);

/**
 * \param grp    pointer toward the group structure.
 * \param reqIdx table of the indices of the required elements.
 * \param nreq   number of required elements
 * \return 1.
 *
 * Set the required Tetra (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_REQUIREDTETRAHEDRA(grp,reqIdx,nreq,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     INTEGER, DIMENSION(*),INTENT(IN) :: reqIdx\n
 * >     INTEGER, INTENT(IN)           :: nreq\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_requiredTetrahedra(PMMG_pGrp grp, int *reqIdx, int nreq);

/**
 * \param grp pointer toward the group structure.
 * \param k   triangle index.
 * \return 1.
 *
 * Set triangle \a k as required (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_REQUIREDTRIANGLE(grp,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_requiredTriangle(PMMG_pGrp grp, int k);

/**
 * \param grp    pointer toward the group structure.
 * \param reqIdx table of the indices of the required trias.
 * \param nreq   number of required trias
 * \return 1.
 *
 * Set the required triangles (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_REQUIREDTRIANGLES(grp,reqIdx,nreq,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)    :: grp\n
 * >     INTEGER, DIMENSION(*),INTENT(IN) :: reqIdx\n
 * >     INTEGER, INTENT(IN)              :: nreq\n
 * >     INTEGER, INTENT(OUT)             :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_requiredTriangles(PMMG_pGrp grp, int *reqIdx, int nreq);

/**
 * \param grp pointer toward the group structure.
 * \param k   edge index.
 * \return 1.
 *
 * Set ridge at edge \a k (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_RIDGE(grp,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_ridge(PMMG_pGrp grp, int k);

/**
 * \param grp pointer toward the group structure.
 * \param k   edge index.
 * \return 1.
 *
 * Set edge \a k as required (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_REQUIREDEDGE(grp,k,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_requiredEdge(PMMG_pGrp grp, int k);

/**
 * \param grp pointer toward the group structure.
 * \param k   point index
 * \param n0  x componant of the normal at point \a k.
 * \param n1  y componant of the normal at point \a k.
 * \param n2  z componant of the normal at point \a k.
 *
 * \return 1 if success.
 *
 * Set normals (n0,n1,n2) at point \a k (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_NORMALATVERTEX(grp,k,n0,n1,n2,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     INTEGER, INTENT(IN)           :: k\n
 * >     REAL(KIND=8), INTENT(IN)      :: n0,n1,n2\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_normalAtVertex(PMMG_pGrp grp, int k, double n0, double n1,
                              double n2);

/**
 * \param grp pointer toward the group structure.
 * \param s   solution scalar value.
 * \param pos position of the solution in the mesh.
 * \return 0  if failed, 1 otherwise.
 *
 * Set scalar value \a s at position \a pos in solution structure (wrapper
 * for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_SCALARSOL(grp,s,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     REAL(KIND=8), INTENT(IN)      :: s\n
 * >     INTEGER, INTENT(IN)           :: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_scalarSol(PMMG_pGrp grp, double s,int pos);

/**
 * \param grp pointer toward the group structure.
 * \param s   table of the scalar solutions values.
 * s[i-1] is the solution at vertex i.
 * \return 0  if failed, 1 otherwise.
 *
 * Set scalar solutions at mesh vertices (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_SCALARSOLS(grp,s,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: s\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_scalarSols(PMMG_pGrp grp, double *sol);

/**
 * \param grp pointer toward the group structure.
 * \param vx  x value of the vectorial solution.
 * \param vy  y value of the vectorial solution.
 * \param vz  z value of the vectorial solution.
 * \param pos position of the solution in the mesh (begin to 1).
 * \return 0  if failed, 1 otherwise.
 *
 * Set vectorial value \f$(v_x,v_y,v_z)\f$ at position \a pos in solution
 * structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_VECTORSOL(grp,vx,vy,vz,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     REAL(KIND=8), INTENT(IN)      :: vx,vy,vz\n
 * >     INTEGER, INTENT(IN)           :: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_vectorSol(PMMG_pGrp grp, double vx,double vy, double vz,
                        int pos);

/**
 * \param grp  pointer toward the group structure.
 * \param sols table of the vectorial solutions
 * sols[3*(i-1)]\@3 is the solution at vertex i
 * \return 0   if failed, 1 otherwise.
 *
 * Set vectorial solutions at mesh vertices (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_VECTORSOLS(grp,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN)      :: sols\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_vectorSols(PMMG_pGrp grp, double *sols);

/**
 * \param grp pointer toward the group structure.
 * \param m11 value of the tensorial solution at position (1,1) in the tensor
 * \param m12 value of the tensorial solution at position (1,2) in the tensor
 * \param m13 value of the tensorial solution at position (1,3) in the tensor
 * \param m22 value of the tensorial solution at position (2,2) in the tensor
 * \param m23 value of the tensorial solution at position (2,3) in the tensor
 * \param m33 value of the tensorial solution at position (3,3) in the tensor
 * \param pos position of the solution in the mesh (begin to 1).
 * \return 0 if failed, 1 otherwise.
 *
 * Set tensorial values at position \a pos in solution
 * structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_TENSORSOL(grp,m11,m12,m13,m22,m23,m33,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     REAL(KIND=8), INTENT(IN)      :: m11,m12,m13,m22,m23,m33\n
 * >     INTEGER, INTENT(IN)           :: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_tensorSol(PMMG_pGrp grp, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos);

/**
 * \param grp  pointer toward the group structure.
 * \param sols table of the tensorial solutions.
 * sols[6*(i-1)]\@6 is the solution at vertex i
 * \return 0   if failed, 1 otherwise.
 *
 * Set tensorial values at position \a pos in solution
 * structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_TENSORSOLS(grp,sols,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)         :: grp\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: sols\n
 * >     INTEGER, INTENT(OUT)                  :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_tensorSols(PMMG_pGrp grp, double *sols);

/**
 * \param grp pointer toward the group structure.
 * \param m   metrics scalar value.
 * \param pos position of the solution in the mesh.
 * \return 0  if failed, 1 otherwise.
 *
 * Set scalar value \a m at position \a pos in metrics structure (wrapper
 * for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_SCALARMET(grp,m,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     REAL(KIND=8), INTENT(IN)      :: m\n
 * >     INTEGER, INTENT(IN)           :: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_scalarMet(PMMG_pGrp grp, double m, int pos);

/**
 * \param grp pointer toward the group structure.
 * \param m   table of the metrics solutions values.
 * m[i-1] is the solution at vertex i.
 * \return 0  if failed, 1 otherwise.
 *
 * Set scalar metrics at mesh vertices (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_SCALARMETS(grp,m,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: m\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_scalarMets(PMMG_pGrp grp, double *m);

/**
 * \param grp pointer toward the group structure.
 * \param vx  x value of the vectorial solution.
 * \param vy  y value of the vectorial solution.
 * \param vz  z value of the vectorial solution.
 * \param pos position of the solution in the mesh (begin to 1).
 * \return 0  if failed, 1 otherwise.
 *
 * Set vectorial value \f$(v_x,v_y,v_z)\f$ at position \a pos in metrics
 * structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_VECTORMET(grp,vx,vy,vz,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     REAL(KIND=8), INTENT(IN)      :: vx,vy,vz\n
 * >     INTEGER, INTENT(IN)           :: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_vectorMet(PMMG_pGrp grp, double vx,double vy, double vz,
                        int pos);

/**
 * \param grp  pointer toward the group structure.
 * \param mets table of the vectorial metrics
 * mets[3*(i-1)]\@3 is the solution at vertex i
 * \return 0   if failed, 1 otherwise.
 *
 * Set vectorial solutions at mesh vertices (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_VECTORMETS(grp,mets,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)         :: grp\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: mets\n
 * >     INTEGER, INTENT(OUT)                  :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_vectorMets(PMMG_pGrp grp, double *mets);

/**
 * \param grp pointer toward the group structure.
 * \param m11 value of the tensorial metrics at position (1,1) in the tensor
 * \param m12 value of the tensorial metrics at position (1,2) in the tensor
 * \param m13 value of the tensorial metrics at position (1,3) in the tensor
 * \param m22 value of the tensorial metrics at position (2,2) in the tensor
 * \param m23 value of the tensorial metrics at position (2,3) in the tensor
 * \param m33 value of the tensorial metrics at position (3,3) in the tensor
 * \param pos position of the metrics in the mesh (begin to 1).
 * \return 0 if failed, 1 otherwise.
 *
 * Set tensorial values at position \a pos in metrics
 * structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_TENSORMET(grp,m11,m12,m13,m22,m23,m33,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * >     REAL(KIND=8), INTENT(IN)      :: m11,m12,m13,m22,m23,m33\n
 * >     INTEGER, INTENT(IN)           :: pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_tensorMet(PMMG_pGrp grp, double m11,double m12, double m13,
                       double m22,double m23, double m33, int pos);

/**
 * \param grp  pointer toward the group structure.
 * \param mets table of the tensorial solutions.
 * mets[6*(i-1)]\@6 is the solution at vertex i
 * \return 0   if failed, 1 otherwise.
 *
 * Set tensorial values at position \a pos in solution
 * structure (wrapper for MMG3D function).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_SET_TENSORMETS(grp,mets,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)         :: grp\n
 * >     REAL(KIND=8),DIMENSION(*), INTENT(IN) :: mets\n
 * >     INTEGER, INTENT(OUT)                  :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG_Set_tensorMets(PMMG_pGrp grp, double *mets);

/**
 * \param grp Pointer towards the group structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_GRPSIZE(grp)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 *
 */
int PMMG_Get_grpSize(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_VERTEX(grp)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * 
 */
int PMMG_Get_vertex(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG_GET_VERTICES(grp)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: grp\n
 * 
 */
int PMMG_Get_vertices(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_tetrahedron(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_tetrahedra(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_prism(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_prisms(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_triangle(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_triangles(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_quadrilateral(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_quadrilaterals(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_edge(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_normalAtVertex(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_scalarSol(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_scalarSols(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_vectorSol(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_vectorSols(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_tensorSol(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_tensorSols(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_scalarMet(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_scalarMets(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_vectorMet(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_vectorMets(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_tensorMet(PMMG_pGrp grp);

/**
 * \param grp Pointer towards the group structure.
 *
 *
 */
int PMMG_Get_tensorMets(PMMG_pGrp grp);
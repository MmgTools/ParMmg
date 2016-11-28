#include "libparmmgtypes.h"
#include "mpi.h"



#define PMMG_VER   "1.0.0"
#define PMMG_REL   "2016"
#define PMMG_CPY   "Copyright (c) Bx INP/INRIA, 2016-"
#define PMMG_STR   "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"

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

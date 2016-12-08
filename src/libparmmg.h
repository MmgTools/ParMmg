#include "libparmmgtypes.h"
#include "mmg3d.h"
#include "libmmg3d.h"
#include "mpitypes.h"

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

/*mergeMesh*/
int PMMG_mergeMesh(PMMG_pParMesh);

#include "libparmmgtypes.h"
#include "mmg3d.h"
#include "mpi.h"



#define PMMG_VER   "1.0.0"
#define PMMG_REL   "2016"
#define PMMG_CPY   "Copyright (c) Bx INP/INRIA, 2016-"
#define PMMG_STR   "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"


int PMMG_saveMesh(PMMG_pParMesh ,const char *);
int PMMG_loadMesh(PMMG_pParMesh ,const char *);

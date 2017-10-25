#ifndef METIS_PMMG_H

#define METIS_PMMG_H

#include "parmmg.h"
#include <metis.h>

int PMMG_graph_meshElts2metis(PMMG_pParMesh,MMG5_pMesh,idx_t**,idx_t**,idx_t*,idx_t);
int PMMG_part_meshElts2metis( PMMG_pParMesh parmesh, idx_t* part, idx_t nproc );

#endif

#ifndef METIS_PMMG_H
#define METIS_PMMG_H
#include "parmmg.h"
#include <metis.h>
int PMMG_partition_metis( PMMG_pParMesh parmesh, idx_t* part, idx_t nproc );
#endif

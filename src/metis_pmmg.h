#ifndef METIS_PMMG_H
#define METIS_PMMG_H
#include "parmmg.h"
#include <metis.h>
int PMMG_mesh2metis( PMMG_pParMesh parmesh, idx_t** xadj, idx_t** adjncy,
                     int *adjsize );
int PMMG_metis_wrapper( idx_t* nelt, idx_t* ncon, idx_t* xadj, idx_t* adjncy,
                        idx_t* nproc, idx_t *objval, idx_t* part );
int PMMG_metispartitioning( PMMG_pParMesh parmesh, idx_t* part );
#endif

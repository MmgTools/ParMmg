#ifndef DEBUG_PMMG_H
#define DEBUG_PMMG_H
#include "parmmg.h"
int adja_idx_of_face( MMG5_pMesh mesh, int element, int face );
int adja_tetra_to_face( MMG5_pMesh mesh, int element, int face );
int adja_face_to_face( MMG5_pMesh mesh, int element, int face );

void grplst_meshes_to_txt( char *name, PMMG_pGrp grp, int ngrp );
void tetras_of_mesh_to_txt( char *name, MMG5_pMesh mesh, int num );
void find_tetras_referencing_null_points_to_txt( char *name, PMMG_pGrp grp, int nmsh );
void listgrp_meshes_adja_of_tetras_to_txt( char *name, PMMG_pGrp grp, int ngrp );
void grplst_meshes_to_saveMesh( PMMG_pGrp listgrp, int ngrp, int rank, char *basename );
void dump_malloc_allocator_info( char *msg, int id );
void check_mem_max_and_mem_cur( PMMG_pParMesh parmesh );
#endif

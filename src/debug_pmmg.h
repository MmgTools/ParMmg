#ifndef DEBUG_PMMG_H
#define DEBUG_PMMG_H
#include "parmmg.h"
int PMMG_adja_idx_of_face( MMG5_pMesh mesh, int element, int face );
int PMMG_adja_tetra_to_face( MMG5_pMesh mesh, int element, int face );
int PMMG_adja_face_to_face( MMG5_pMesh mesh, int element, int face );

void PMMG_grplst_meshes_to_txt( char *name, PMMG_pGrp grp, int ngrp );
void PMMG_tetras_of_mesh_to_txt( char *name, MMG5_pMesh mesh, int num );
void PMMG_find_tetras_referencing_null_points_to_txt( char *name, PMMG_pGrp grp, int nmsh );
void PMMG_listgrp_meshes_adja_of_tetras_to_txt( char *name, PMMG_pGrp grp, int ngrp );
int PMMG_grp_to_saveMesh( PMMG_pParMesh parmesh, int grpId, char *basename );
int PMMG_listgrp_to_saveMesh( PMMG_pParMesh parmesh, char *basename );
void PMMG_dump_malloc_allocator_info( char *msg, int id );
void PMMG_check_mem_max_and_mem_cur( PMMG_pParMesh parmesh, const char *msg );
#endif

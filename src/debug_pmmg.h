/* =============================================================================
**  This file is part of the parmmg software package for parallel tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux, 2017-
**
**  parmmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  parmmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with parmmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the parmmg distribution only if you accept them.
** =============================================================================
*/

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
int PMMG_grp_to_saveEdges( PMMG_pParMesh parmesh,int grpId,int16_t tag,char *basename );
int PMMG_grp_to_saveMesh( PMMG_pParMesh parmesh, int grpId, char *basename );
int PMMG_listgrp_to_saveMesh( PMMG_pParMesh parmesh, char *basename );
int PMMG_listgrp_quality_to_saveMesh( PMMG_pParMesh parmesh, char *basename );
int PMMG_listgrp_mark_to_saveMesh( PMMG_pParMesh parmesh, char *basename );
void PMMG_print_ext_comm( PMMG_pParMesh parmesh, PMMG_pInt_comm int_comm,
                          PMMG_pExt_comm ext_comm, int next_comm );
void PMMG_dump_malloc_allocator_info( char *msg, int id );
void PMMG_check_mem_max_and_mem_cur( PMMG_pParMesh parmesh, const char *msg );
#endif

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

/**
 * \brief functions to pack a group mesh into a char buffer (for mpi call)
 * \author Luca Cirrottola (Inria)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 */

#include "parmmg.h"

/**
 * \param grp pointer toward a PMMG group
 *
 * \return the computed size
 *
 * Compute the size of compressed mesh and solution data related to array sizes.
 *
 */
static
int PMMG_mpisizeof_meshSizes ( PMMG_pGrp grp ) {
  const MMG5_pMesh mesh = grp->mesh;
  const MMG5_pSol  met  = grp->met;
  const MMG5_pSol  ls   = grp->ls;
  const MMG5_pSol  disp = grp->disp;
  int              idx = 0;

  /** Mesh size */
  idx += sizeof(int); // mesh->np
  idx += sizeof(int); // mesh->xp
  idx += sizeof(int); // mesh->ne
  idx += sizeof(int); // mesh->xt
  idx += sizeof(int); // mesh->nsols

  /** Met size */
  idx += sizeof(int); // met && met->m
  if ( met && met->m ) {
    idx += 2*sizeof(int); // met->size/type;
  }

  /** ls size */
  idx += sizeof(int);
  if ( ls && ls->m ) {
    idx += 2*sizeof(int);
  }

  /** disp size */
  idx += sizeof(int);
  if ( disp && disp->m ) {
    idx += 2*sizeof(int);
  }

  /** Fields info and sizes */
  if ( mesh->nsols ) {
    idx += mesh->nsols*sizeof(int); // psl->size;
    idx += mesh->nsols*sizeof(int); // psl->type;
  }

  return idx;
}

/**
 * \param grp pointer toward a PMMG group
 *
 * \return the computed size
 *
 * Compute the size of compressed mesh and solution names.
 *
 */
static
int PMMG_mpisizeof_filenames ( PMMG_pGrp grp ) {
  const MMG5_pMesh mesh = grp->mesh;
  const MMG5_pSol  met  = grp->met;
  const MMG5_pSol  ls   = grp->ls;
  const MMG5_pSol  disp = grp->disp;
  MMG5_pSol        psl;
  int              is;
  int              idx = 0;

  /** Mesh names */
  idx += sizeof(int); // meshin
  idx += sizeof(int); // meshout
  idx += (strlen(mesh->namein) + 1) * sizeof(char);
  idx += (strlen(mesh->nameout) + 1) * sizeof(char);


  /* Met */
  idx += sizeof(int); // metin
  idx += sizeof(int); // metout
  if ( met && met->namein && met->nameout ) {
    idx += (strlen(met->namein) + 1) * sizeof(char);
    idx += (strlen(met->nameout) + 1) * sizeof(char);
  }

  /* ls */
  idx += sizeof(int);
  idx += sizeof(int);
  if ( ls && ls->namein && ls->nameout ) {
    idx += (strlen(ls->namein) + 1) * sizeof(char);
    idx += (strlen(ls->nameout) + 1) * sizeof(char);
  }

  /* disp */
  idx += sizeof(int); // metin
  idx += sizeof(int); // metout
  if ( disp && disp->namein && disp->nameout ) {
    idx += (strlen(disp->namein) + 1) * sizeof(char);
    idx += (strlen(disp->nameout) + 1) * sizeof(char);
  }

  /* Fields */
  if ( mesh->nsols ) {
    idx += mesh->nsols*sizeof(int);
    idx += mesh->nsols*sizeof(int);
    for ( is=0; is<mesh->nsols; ++is ) {
      psl = &grp->field[is];
      if ( psl && psl->namein && psl->nameout ) {
        idx += (strlen(psl->namein) + 1) * sizeof(char);
        idx += (strlen(psl->nameout) + 1) * sizeof(char);
      }
    }
  }

  return idx;
}

/**
 * \param info pointer toward a MMG5 info structure
 *
 * \return the computed size
 *
 * Compute the size of compressed info structure.
 *
 */
static
int PMMG_mpisizeof_infos ( MMG5_Info *info ) {
  int idx = 0;

  /** Mesh infos: warning, some "useless" info are not sended */
  idx += sizeof(double); // mesh->info.dhd
  idx += sizeof(double); // mesh->info.hmin
  idx += sizeof(double); // mesh->info.hmax
  idx += sizeof(double); // mesh->info.hsiz
  idx += sizeof(double); // mesh->info.hgrad
  idx += sizeof(double); // mesh->info.hgradreq
  idx += sizeof(double); // mesh->info.hausd

  idx += sizeof(double); // mesh->info.delta
  idx += sizeof(double); // mesh->info.min[0]
  idx += sizeof(double); // mesh->info.min[1]
  idx += sizeof(double); // mesh->info.min[2]

  idx += sizeof(double); // mesh->info.ls

  idx += sizeof(int); // npar
  idx += sizeof(int); // openbdy
  idx += sizeof(int); // renum
  idx += sizeof(int); // PROctree
  idx += sizeof(int); // nmat

  idx += sizeof(char); // nreg
  idx += sizeof(char); // imprim
  idx += sizeof(char); // ddebug
  idx += sizeof(char); // iso
  idx += sizeof(char); // lag
  idx += sizeof(char); // parTyp
  idx += sizeof(char); // optim
  idx += sizeof(char); // optimLES
  idx += sizeof(char); // noinsert
  idx += sizeof(char); // noswap
  idx += sizeof(char); // nomove
  idx += sizeof(char); // nosurf
  idx += sizeof(char); // sethmin
  idx += sizeof(char); // sethmax

  /* affectation of old refs in ls-mode */
  if ( info->nmat ) {
    assert( info->mat );
    idx += info->nmat*sizeof(char); // mat->dospl
    idx += info->nmat*sizeof(int); //  mat->ref
    idx += info->nmat*sizeof(int); //  mat->rin
    idx += info->nmat*sizeof(int); //  mat->rex
  }

  /* local parameters */
  if ( info->npar ) {
    assert( info->par );
    idx += info->npar*sizeof(double); // par->hmin
    idx += info->npar*sizeof(double); // par->hmax
    idx += info->npar*sizeof(double); // par->hausd
    idx += info->npar*sizeof(int); //  par->ref
    idx += info->npar*sizeof(char); // par->elt
  }
  return idx;
}

/**
 * \param grp pointer toward a PMMG group
 *
 * \return the computed size
 *
 * Compute the size of compressed mesh and solutions arrays.
 *
 */
static
int PMMG_mpisizeof_meshArrays ( PMMG_pGrp grp ) {
  const MMG5_pMesh mesh = grp->mesh;
  const MMG5_pSol  met  = grp->met;
  const MMG5_pSol  ls   = grp->ls;
  const MMG5_pSol  disp = grp->disp;
  MMG5_pSol        psl;
  int              idx = 0;
  int              is;

  /** Pack mesh points */
  /* Coordinates */
  idx += mesh->np*sizeof(double); // mesh->point[k].c[0];
  idx += mesh->np*sizeof(double); // mesh->point[k].c[1];
  idx += mesh->np*sizeof(double); // mesh->point[k].c[2];
  /* Tangent */
  idx += mesh->np*sizeof(double); // mesh->point[k].n[0];
  idx += mesh->np*sizeof(double); // mesh->point[k].n[1];
  idx += mesh->np*sizeof(double); // mesh->point[k].n[2];
  /* Pointer toward the boundary entity */
  idx += mesh->np*sizeof(int); // mesh->point[k].xp;
  /* Ref */
  idx += mesh->np*sizeof(int); // mesh->point[k].ref;
  /* Tag */
  idx += mesh->np*sizeof(int16_t); // mesh->point[k].tag;
  /* Src */
  idx += mesh->np*sizeof(int); // mesh->point[k].src;

  /** Pack mesh boundary points */
  /* First normal */
  idx += mesh->xp*sizeof(double); // mesh->xpoint[k].n1[0];
  idx += mesh->xp*sizeof(double); // mesh->xpoint[k].n1[1];
  idx += mesh->xp*sizeof(double); // mesh->xpoint[k].n1[2];
  /* Second normal */
  idx += mesh->xp*sizeof(double); // mesh->xpoint[k].n2[0];
  idx += mesh->xp*sizeof(double); // mesh->xpoint[k].n2[1];
  idx += mesh->xp*sizeof(double); // mesh->xpoint[k].n2[2];

  /** Pack mesh elements */
  /* Tetra vertices */
  idx += mesh->ne*sizeof(int); // mesh->tetra[k].v[0];
  idx += mesh->ne*sizeof(int); // mesh->tetra[k].v[1];
  idx += mesh->ne*sizeof(int); // mesh->tetra[k].v[2];
  idx += mesh->ne*sizeof(int); // mesh->tetra[k].v[3];
  /* Pointer toward the boundary entity */
  idx += mesh->ne*sizeof(int); // mesh->tetra[k].xt;
  /* Ref */
  idx += mesh->ne*sizeof(int); // mesh->tetra[k].ref;
  /* Mark */
  idx += mesh->ne*sizeof(int); // mesh->point[k].mark;
  /* Tag */
  idx += mesh->ne*sizeof(int16_t); // mesh->tetra[k].tag;
  /* Quality */
  idx += mesh->ne*sizeof(double); // mesh->tetra[k].qual;

  /** Pack mesh boundary tetra */
  /* Faces references  */
  idx += mesh->xt*sizeof(int); // mesh->xtetra[k].ref[0];
  idx += mesh->xt*sizeof(int); // mesh->xtetra[k].ref[1];
  idx += mesh->xt*sizeof(int); // mesh->xtetra[k].ref[2];
  idx += mesh->xt*sizeof(int); // mesh->xtetra[k].ref[3];
  /* Faces tags */
  idx += mesh->xt*sizeof(int16_t); // mesh->xtetra[k].ftag[0];
  idx += mesh->xt*sizeof(int16_t); // mesh->xtetra[k].ftag[1];
  idx += mesh->xt*sizeof(int16_t); // mesh->xtetra[k].ftag[2];
  idx += mesh->xt*sizeof(int16_t); // mesh->xtetra[k].ftag[3];
  /* Edges references  */
  idx += mesh->xt*sizeof(int); // mesh->xtetra[k].edg[0];
  idx += mesh->xt*sizeof(int); // mesh->xtetra[k].edg[1];
  idx += mesh->xt*sizeof(int); // mesh->xtetra[k].edg[2];
  idx += mesh->xt*sizeof(int); // mesh->xtetra[k].edg[3];
  idx += mesh->xt*sizeof(int); // mesh->xtetra[k].edg[4];
  idx += mesh->xt*sizeof(int); // mesh->xtetra[k].edg[5];
  /* Edges tags */
  idx += mesh->xt*sizeof(int16_t); // mesh->xtetra[k].tag[0];
  idx += mesh->xt*sizeof(int16_t); // mesh->xtetra[k].tag[1];
  idx += mesh->xt*sizeof(int16_t); // mesh->xtetra[k].tag[2];
  idx += mesh->xt*sizeof(int16_t); // mesh->xtetra[k].tag[3];
  idx += mesh->xt*sizeof(int16_t); // mesh->xtetra[k].tag[4];
  idx += mesh->xt*sizeof(int16_t); // mesh->xtetra[k].tag[5];
  /* Orientation of the triangles */
  idx += mesh->xt*sizeof(char); // mesh->xtetra[k].ori;

  /** Pack metric */
  if ( met && met->m ) {
    idx += met->size*met->np*sizeof(double); // met->m;
  }

  /** Pack ls */
  if ( ls && ls->m ) {
    idx += ls->size*ls->np*sizeof(double); // ls->m;
  }
  /** Pack disp */
  if ( disp && disp->m ) {
    idx += disp->size*disp->np*sizeof(double); // disp->m;
  }

  /** Pack Fields  */
  if ( mesh->nsols ) {
    assert ( grp->field );
    for ( is=0; is<mesh->nsols; ++is ) {
      psl = &grp->field[is];
      idx += psl->size*psl->np*sizeof(double); // psl->m;
    }
  }

  return idx;
}

/**
 * \param grp pointer toward a PMMG group
 *
 * \return the computed size
 *
 * Compute the size of compressed internal group communicators size and arrays.
 *
 */
static
int PMMG_mpisizeof_grpintcomm ( PMMG_pGrp grp ) {
  int              idx = 0;

  /** Pack communicators */
  /* Communicator sizes */
  idx += sizeof(int); // grp->nitem_int_node_comm;
  idx += sizeof(int); // grp->nitem_int_face_comm;

  /* Node communicators */
  idx += 2*grp->nitem_int_node_comm*sizeof(int); // grp->node2int_node_comm_index1-2

  /* Face communicators */
  idx += 2*grp->nitem_int_face_comm*sizeof(int); // grp->face2int_face_comm_index1-2

  return idx;
}

/**
 * \param parmesh pointer toward a PMMG parmesh
 *
 * \return the computed size
 *
 * Compute the size of compressed intvalue array of the nodal internal
 * communicator.
 *
 */
static
int PMMG_mpisizeof_nodeintvalues ( PMMG_pParMesh parmesh ) {
  int              idx = 0;

  /** Pack intvalues array of nodal communicator */
  /* Array size */
  idx += sizeof(int);

  /* Intvalues array */
  idx += parmesh->int_node_comm->nitem*sizeof(int);

  return idx;
}

/**
 * \param parmesh pointer toward a PMMG parmesh
 *
 * \return the computed size
 *
 * Compute the size of compressed nodal external communicator.
 *
 */
static
int PMMG_mpisizeof_extnodecomm ( PMMG_pParMesh parmesh ) {
  PMMG_pExt_comm ext_node_comm;
  int            k,idx = 0;

  /** Pack nodal external communicators */
  /* Number of external communicators */
  idx += sizeof(int);

  /* List of external communicators */
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];

    idx += 3*sizeof(int); // Colors + nitem
    idx += parmesh->ext_node_comm[k].nitem * sizeof(int);
  }

  return idx;
}


/**
 * \param grp pointer toward a PMMG_Grp structure.
 * \return the size (in char) of the packed group.
 *
 * \warning the mesh prisms are not treated.
 *
 * Compute the size of the compressed group.
 *
 */
int PMMG_mpisizeof_grp ( PMMG_pGrp grp ) {
  const MMG5_pMesh mesh = grp->mesh;

  int idx;

  /** Used or unused group */
  idx = sizeof(int);

  if ( !grp->mesh ) {
    /* unused group */
    return idx;
  }

  /** Size of mesh / metric / fields... arrays (np, met->siz...) */
  idx += PMMG_mpisizeof_meshSizes ( grp );

  /** Size of Info */
  idx += PMMG_mpisizeof_infos ( &mesh->info );

  /** Size of compressed points / tetra / metric / fields... */
  idx += PMMG_mpisizeof_meshArrays ( grp );

  /** Size of compressed internal group communicators */
  idx += PMMG_mpisizeof_grpintcomm ( grp );

  return idx;
}

/**
 * \param grp pointer toward a PMMG_Grp structure.
 * \return the size (in char) of the packed group.
 *
 * \warning the mesh prisms are not treated.
 *
 * Compute the size of the compressed parmesh (groups must have been merged
 * before entering this function).
 *
 */
int PMMG_mpisizeof_parmesh ( PMMG_pParMesh parmesh ) {
  PMMG_pGrp grp;
  int       idx;

  assert ( parmesh->ngrp < 2 ); // Check that groups are merged

  /** Used or unused group */
  idx = sizeof(int);

  if ( (!parmesh->ngrp) || (!parmesh->listgrp[0].mesh) ) {
    /* unused group */
    return idx;
  }

  grp = &parmesh->listgrp[0];

  /** Size of mesh / metric / fields... arrays (np, met->siz...) */
  idx += PMMG_mpisizeof_meshSizes ( grp );

  /** Size of Info */
  idx += PMMG_mpisizeof_infos ( &parmesh->listgrp[0].mesh->info );

  /** Size of compressed points / tetra / metric / fields... */
  idx += PMMG_mpisizeof_meshArrays ( grp );

  /** Size of intvalues array */
  idx += PMMG_mpisizeof_nodeintvalues ( parmesh );

  /** Size of external communicators */
  idx += PMMG_mpisizeof_extnodecomm ( parmesh );

  return idx;
}

/**
 * \param grp pointer toward a PMMG_Grp structure.
 * \param buffer pointer toward the buffer in which we pack the group
 *
 * \warning the mesh prisms are not treated.
 *
 * Pack the size of the mesh/met/fields... arrays to send and shift the buffer
 * pointer at the end of the written area.
 *
 */
static
void PMMG_mpipack_meshSizes ( PMMG_pGrp grp,char **buffer ) {
  const MMG5_pMesh mesh = grp->mesh;
  const MMG5_pSol  met  = grp->met;
  const MMG5_pSol  ls   = grp->ls;
  const MMG5_pSol  disp = grp->disp;
  MMG5_pSol        psl;

  int   is;
  char  *tmp;

  tmp = *buffer;

  /** Mesh size */
  *( (int *) tmp) = mesh->np; tmp += sizeof(int);
  *( (int *) tmp) = mesh->xp; tmp += sizeof(int);
  *( (int *) tmp) = mesh->ne; tmp += sizeof(int);
  *( (int *) tmp) = mesh->xt; tmp += sizeof(int);
  *( (int *) tmp) = mesh->nsols; tmp += sizeof(int);

  /** Metric info and sizes */
  *( (int *) tmp) = ( (met && met->m) ? 1 : 0 );  tmp += sizeof(int);
  if ( met && met->m ) {
    assert ( mesh->npmax == met->npmax );
    assert ( mesh->np    == met->np    );

    *( (int *) tmp) = met->size;          tmp += sizeof(int);
    *( (int *) tmp) = met->type;          tmp += sizeof(int);
  }

  /** Ls info and sizes */
  *( (int *) tmp) = ( (ls && ls->m) ? 1 : 0 );  tmp += sizeof(int);
  if ( ls && ls->m ) {
    assert ( mesh->npmax == ls->npmax );
    assert ( mesh->np    == ls->np    );

    *( (int *) tmp) = ls->size;          tmp += sizeof(int);
    *( (int *) tmp) = ls->type;          tmp += sizeof(int);
  }

  /** Disp info and sizes */
  *( (int *) tmp) = ( (disp && disp->m) ? 1 : 0 );  tmp += sizeof(int);
  if ( disp && disp->m ) {
    assert ( mesh->npmax == disp->npmax );
    assert ( mesh->np    == disp->np    );

    *( (int *) tmp) = disp->size;          tmp += sizeof(int);
    *( (int *) tmp) = disp->type;          tmp += sizeof(int);
  }

  /** Fields info and sizes */
  if ( mesh->nsols ) {
    for ( is=0; is<mesh->nsols; ++is ) {
      psl = &grp->field[is];
      assert ( mesh->npmax == psl->npmax );
      assert ( mesh->np    == psl->np    );

      *( (int *) tmp) = psl->size;          tmp += sizeof(int);
      *( (int *) tmp) = psl->type;          tmp += sizeof(int);
    }
  }

  *buffer = tmp;
}

/**
 * \param grp pointer toward a PMMG_Grp structure.
 * \param buffer pointer toward the buffer in which we pack the group
 *
 * Pack the filenames and shift the buffer
 * pointer at the end of the written area.
 *
 */
static
int PMMG_mpipack_filenames ( PMMG_pGrp grp,char **buffer ) {
  const MMG5_pMesh mesh = grp->mesh;
  const MMG5_pSol  met  = grp->met;
  const MMG5_pSol  ls   = grp->ls;
  const MMG5_pSol  disp = grp->disp;
  MMG5_pSol        psl;

  int   is,meshin_s,meshout_s,metin_s,metout_s;
  int   lsin_s,lsout_s,dispin_s,dispout_s,pslin_s,pslout_s;
  char  *tmp;

  tmp = *buffer;

  /** Mesh names */
  meshin_s  = (strlen(mesh->namein) + 1);
  meshout_s = (strlen(mesh->nameout) + 1);
  if ( meshin_s  > MMG5_FILENAME_LEN_MAX || meshout_s > MMG5_FILENAME_LEN_MAX ) {
    printf("  ## Error: input meshnames too long.");
    return 0;
  }
  *( (int *) tmp) = meshin_s; tmp += sizeof(int);
  *( (int *) tmp) = meshout_s; tmp += sizeof(int);
  strncpy( ( (char *) tmp), mesh->namein,meshin_s);
  tmp += meshin_s * sizeof(char);
  strncpy( ( (char *) tmp), mesh->nameout,meshout_s);
  tmp += meshout_s * sizeof(char);

  if ( met && met->namein && met->nameout ) {
    metin_s   = (strlen(met->namein) + 1);
    metout_s  = (strlen(met->nameout) + 1);
    if ( metin_s  > MMG5_FILENAME_LEN_MAX || metout_s > MMG5_FILENAME_LEN_MAX ) {
      printf("  ## Error: input metnames too long.");
      return 0;
    }
  }
  else {
    metin_s  = 0;
    metout_s = 0;
  }

  *( (int *) tmp) = metin_s; tmp += sizeof(int);
  *( (int *) tmp) = metout_s; tmp += sizeof(int);

  if ( metin_s ) {
    assert ( metout_s );
    strncpy( ( (char *) tmp), met->namein,metin_s);
    tmp += metin_s * sizeof(char);
    strncpy( ( (char *) tmp), met->nameout, metout_s);
    tmp += metout_s * sizeof(char);
  }

  if ( ls && ls->namein && ls->nameout ) {
    lsin_s   = (strlen(ls->namein) + 1);
    lsout_s  = (strlen(ls->nameout) + 1);
    if ( lsin_s  > MMG5_FILENAME_LEN_MAX || lsout_s > MMG5_FILENAME_LEN_MAX ) {
      printf("  ## Error: input level-set names too long.");
      return 0;
    }
  }
  else {
    lsin_s  = 0;
    lsout_s = 0;
  }

  *( (int *) tmp) = lsin_s; tmp += sizeof(int);
  *( (int *) tmp) = lsout_s; tmp += sizeof(int);

  if ( lsin_s ) {
    assert ( lsout_s );
    strncpy( ( (char *) tmp), ls->namein,lsin_s);
    tmp += lsin_s * sizeof(char);
    strncpy( ( (char *) tmp), ls->nameout, lsout_s);
    tmp += lsout_s * sizeof(char);
  }

  if ( disp && disp->namein && disp->nameout ) {
    dispin_s   = (strlen(disp->namein) + 1);
    dispout_s  = (strlen(disp->nameout) + 1);
    if ( dispin_s  > MMG5_FILENAME_LEN_MAX || dispout_s > MMG5_FILENAME_LEN_MAX ) {
      printf("  ## Error: input displacement names too long.");
      return 0;
    }
  }
  else {
    dispin_s  = 0;
    dispout_s = 0;
  }

  *( (int *) tmp) = dispin_s; tmp += sizeof(int);
  *( (int *) tmp) = dispout_s; tmp += sizeof(int);

  if ( dispin_s ) {
    assert ( dispout_s );
    strncpy( ( (char *) tmp), disp->namein,dispin_s);
    tmp += dispin_s * sizeof(char);
    strncpy( ( (char *) tmp), disp->nameout, dispout_s);
    tmp += dispout_s * sizeof(char);
  }


  for ( is=0; is<mesh->nsols; ++is ) {
    psl = &grp->field[is];
    if ( psl && psl->namein && psl->nameout ) {
      pslin_s   = (strlen(psl->namein) + 1);
      pslout_s  = (strlen(psl->nameout) + 1);
      if ( pslin_s  > MMG5_FILENAME_LEN_MAX || pslout_s > MMG5_FILENAME_LEN_MAX ) {
        printf("  ## Error: input displacement names too long.");
        return 0;
      }
    }
    else {
      pslin_s  = 0;
      pslout_s = 0;
    }

    *( (int *) tmp) = pslin_s; tmp += sizeof(int);
    *( (int *) tmp) = pslout_s; tmp += sizeof(int);

    if ( pslin_s ) {
      assert ( pslout_s );
      strncpy( ( (char *) tmp), psl->namein,pslin_s);
      tmp += pslin_s * sizeof(char);
      strncpy( ( (char *) tmp), psl->nameout, pslout_s);
      tmp += pslout_s * sizeof(char);
    }
  }

  *buffer = tmp;
  return 1;
}

/**
 * \param info pointer toward a MMG5 info structure
 * \param buffer pointer toward the buffer in which we pack the group
 *
 * \return 1 if success, 0 if fail
 *
 * Pack the info structure and shift the buffer
 * pointer at the end of the written area.
 *
 */
static
void PMMG_mpipack_infos ( MMG5_Info *info,char **buffer ) {
  int   k;
  char  *tmp;

  tmp = *buffer;

  /** Mesh infos */
  *( (double *) tmp) = info->dhd;      tmp += sizeof(double);
  *( (double *) tmp) = info->hmin;     tmp += sizeof(double);
  *( (double *) tmp) = info->hmax;     tmp += sizeof(double);
  *( (double *) tmp) = info->hsiz;     tmp += sizeof(double);
  *( (double *) tmp) = info->hgrad;    tmp += sizeof(double);
  *( (double *) tmp) = info->hgradreq; tmp += sizeof(double);
  *( (double *) tmp) = info->hausd;    tmp += sizeof(double);

  *( (double *) tmp) = info->delta;    tmp += sizeof(double);
  *( (double *) tmp) = info->min[0];   tmp += sizeof(double);
  *( (double *) tmp) = info->min[1];   tmp += sizeof(double);
  *( (double *) tmp) = info->min[2];   tmp += sizeof(double);

  *( (double *) tmp) = info->ls;       tmp += sizeof(double);

  *( (int *) tmp) = info->opnbdy;    tmp += sizeof(int);
  *( (int *) tmp) = info->renum;     tmp += sizeof(int);
  *( (int *) tmp) = info->PROctree;  tmp += sizeof(int);
  *( (int *) tmp) = info->npar;      tmp += sizeof(int);
  *( (int *) tmp) = info->nmat;      tmp += sizeof(int);

  *( (char *) tmp) = info->nreg;     tmp += sizeof(char);
  *( (char *) tmp) = info->imprim;   tmp += sizeof(char);
  *( (char *) tmp) = info->ddebug;   tmp += sizeof(char);
  *( (char *) tmp) = info->iso;      tmp += sizeof(char);
  *( (char *) tmp) = info->lag;      tmp += sizeof(char);
  *( (char *) tmp) = info->parTyp;   tmp += sizeof(char);
  *( (char *) tmp) = info->optim;    tmp += sizeof(char);
  *( (char *) tmp) = info->optimLES; tmp += sizeof(char);
  *( (char *) tmp) = info->noinsert; tmp += sizeof(char);
  *( (char *) tmp) = info->noswap;   tmp += sizeof(char);
  *( (char *) tmp) = info->nomove;   tmp += sizeof(char);
  *( (char *) tmp) = info->nosurf;   tmp += sizeof(char);
  *( (char *) tmp) = info->sethmin;  tmp += sizeof(char);
  *( (char *) tmp) = info->sethmax;  tmp += sizeof(char);

  /* affectation of old refs in ls-mode */
  if ( info->nmat ) {
    assert( info->mat );
    for ( k=0; k<info->nmat; ++k ) {
      *( (char *) tmp) = info->mat[k].dospl; tmp += sizeof(char);
      *( (int *) tmp)  = info->mat[k].ref; tmp += sizeof(int);
      *( (int *) tmp)  = info->mat[k].rin; tmp += sizeof(int);
      *( (int *) tmp)  = info->mat[k].rex; tmp += sizeof(int);
    }
  }

  /* local parameters */
  if ( info->npar ) {
    assert( info->par );
    for ( k=0; k<info->npar; ++k ) {
      *( (double *) tmp) =  info->par[k].hmin; tmp += sizeof(double);
      *( (double *) tmp) =  info->par[k].hmax; tmp += sizeof(double);
      *( (double *) tmp) =  info->par[k].hausd; tmp += sizeof(double);
      *( (int *)    tmp) =  info->par[k].ref; tmp += sizeof(int);
      *( (char *)   tmp) =  info->par[k].elt; tmp += sizeof(char);
    }
  }

  *buffer = tmp;
}

/**
 * \param grp pointer toward a PMMG_Grp structure.
 * \param buffer pointer toward the buffer in which we pack the group
 *
 * Pack the mesh and solutions (metric, ls, disp and fields) arrays.
 *
 */
static
void PMMG_mpipack_meshArrays ( PMMG_pGrp grp,char **buffer ) {
  const MMG5_pMesh mesh = grp->mesh;
  const MMG5_pSol  met  = grp->met;
  const MMG5_pSol  ls   = grp->ls;
  const MMG5_pSol  disp = grp->disp;
  MMG5_pSol        psl;

  int   k,i,is;
  char  *tmp;

  tmp = *buffer;

  /** Pack mesh points */
  for ( k=1; k<=mesh->np; ++k ) {
    /* Coordinates */
    *( (double *) tmp) = mesh->point[k].c[0]; tmp += sizeof(double);
    *( (double *) tmp) = mesh->point[k].c[1]; tmp += sizeof(double);
    *( (double *) tmp) = mesh->point[k].c[2]; tmp += sizeof(double);
    /* Tangent */
    *( (double *) tmp) = mesh->point[k].n[0]; tmp += sizeof(double);
    *( (double *) tmp) = mesh->point[k].n[1]; tmp += sizeof(double);
    *( (double *) tmp) = mesh->point[k].n[2]; tmp += sizeof(double);
    /* Pointer toward the boundary entity */
    *( (int *) tmp) = mesh->point[k].xp; tmp += sizeof(int);
    /* Ref */
    *( (int *) tmp) = mesh->point[k].ref; tmp += sizeof(int);
    /* Tag */
    *( (int16_t *) tmp) = mesh->point[k].tag; tmp += sizeof(int16_t);
    /* Src */
    *( (int *) tmp) = mesh->point[k].src; tmp += sizeof(int);
  }

  /** Pack mesh boundary points */
  for ( k=1; k<=mesh->xp; ++k ) {
    /* First normal */
    *( (double *) tmp) = mesh->xpoint[k].n1[0]; tmp += sizeof(double);
    *( (double *) tmp) = mesh->xpoint[k].n1[1]; tmp += sizeof(double);
    *( (double *) tmp) = mesh->xpoint[k].n1[2]; tmp += sizeof(double);
    /* Second normal */
    *( (double *) tmp) = mesh->xpoint[k].n2[0]; tmp += sizeof(double);
    *( (double *) tmp) = mesh->xpoint[k].n2[1]; tmp += sizeof(double);
    *( (double *) tmp) = mesh->xpoint[k].n2[2]; tmp += sizeof(double);
  }

  /** Pack mesh elements */
  for ( k=1; k<=mesh->ne; ++k ) {
    /* Tetra vertices */
    *( (int *) tmp) = mesh->tetra[k].v[0]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->tetra[k].v[1]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->tetra[k].v[2]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->tetra[k].v[3]; tmp += sizeof(int);
    /* Pointer toward the boundary entity */
    *( (int *) tmp) = mesh->tetra[k].xt; tmp += sizeof(int);
    /* Ref */
    *( (int *) tmp) = mesh->tetra[k].ref; tmp += sizeof(int);
    /* Mark */
    *( (int *) tmp) = mesh->tetra[k].mark; tmp += sizeof(int);
    /* Tag */
    *( (int16_t *) tmp) = mesh->tetra[k].tag; tmp += sizeof(int16_t);
    /* Quality */
    *( (double *) tmp) = mesh->tetra[k].qual; tmp += sizeof(double);
  }

  /** Pack mesh boundary tetra */
  for ( k=1; k<=mesh->xt; ++k ) {
    /* Faces references  */
    *( (int *) tmp) = mesh->xtetra[k].ref[0]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->xtetra[k].ref[1]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->xtetra[k].ref[2]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->xtetra[k].ref[3]; tmp += sizeof(int);
    /* Faces tags */
    *( (int16_t *) tmp) = mesh->xtetra[k].ftag[0]; tmp += sizeof(int16_t);
    *( (int16_t *) tmp) = mesh->xtetra[k].ftag[1]; tmp += sizeof(int16_t);
    *( (int16_t *) tmp) = mesh->xtetra[k].ftag[2]; tmp += sizeof(int16_t);
    *( (int16_t *) tmp) = mesh->xtetra[k].ftag[3]; tmp += sizeof(int16_t);
    /* Edges references  */
    *( (int *) tmp) = mesh->xtetra[k].edg[0]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->xtetra[k].edg[1]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->xtetra[k].edg[2]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->xtetra[k].edg[3]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->xtetra[k].edg[4]; tmp += sizeof(int);
    *( (int *) tmp) = mesh->xtetra[k].edg[5]; tmp += sizeof(int);
    /* Edges tags */
    *( (int16_t *) tmp) = mesh->xtetra[k].tag[0]; tmp += sizeof(int16_t);
    *( (int16_t *) tmp) = mesh->xtetra[k].tag[1]; tmp += sizeof(int16_t);
    *( (int16_t *) tmp) = mesh->xtetra[k].tag[2]; tmp += sizeof(int16_t);
    *( (int16_t *) tmp) = mesh->xtetra[k].tag[3]; tmp += sizeof(int16_t);
    *( (int16_t *) tmp) = mesh->xtetra[k].tag[4]; tmp += sizeof(int16_t);
    *( (int16_t *) tmp) = mesh->xtetra[k].tag[5]; tmp += sizeof(int16_t);
    /* Orientation of the triangles */
    *( (char *) tmp) = mesh->xtetra[k].ori; tmp += sizeof(char);
  }

  /** Pack metric */
  if ( met && met->m ) {
    for ( k=1; k<=met->np; ++k ) {
      for ( i=0; i<met->size; ++i ) {
        *( (double *) tmp) = met->m[met->size*k + i]; tmp += sizeof(double);
      }
    }
  }

  /** Pack ls */
  if ( ls && ls->m ) {
    for ( k=1; k<=ls->np; ++k ) {
      for ( i=0; i<ls->size; ++i ) {
        *( (double *) tmp) = ls->m[ls->size*k + i]; tmp += sizeof(double);
      }
    }
  }

  /** Pack disp */
  if ( disp && disp->m ) {
    for ( k=1; k<=disp->np; ++k ) {
      for ( i=0; i<disp->size; ++i ) {
        *( (double *) tmp) = disp->m[disp->size*k + i]; tmp += sizeof(double);
      }
    }
  }

  /** Pack Fields */
  if ( mesh->nsols ) {
    for ( is=0; is<mesh->nsols; ++is ) {
      psl = &grp->field[is];
      for ( k=1; k<=psl->np; ++k ) {
        for ( i=0; i<psl->size; ++i ) {
          *( (double *) tmp) = psl->m[psl->size*k + i]; tmp += sizeof(double);
        }
      }
    }
  }

  *buffer = tmp;
}

/**
 * \param grp pointer toward a PMMG_Grp structure.
 * \param buffer pointer toward the buffer in which we pack the group
 *
 * Pack the internal communicators of the group
 *
 */
static
void PMMG_mpipack_grpintcomm ( PMMG_pGrp grp,char **buffer ) {
  int   k;
  char  *tmp;

  tmp = *buffer;

  /** Pack communicators */
  /* Communicator sizes */
  *( (int *) tmp) = grp->nitem_int_node_comm; tmp += sizeof(int);
  *( (int *) tmp) = grp->nitem_int_face_comm; tmp += sizeof(int);

  /* Node communicator */
  for ( k=0; k<grp->nitem_int_node_comm; ++k ) {
    *( (int *) tmp) = grp->node2int_node_comm_index1[k]; tmp += sizeof(int);
  }
  for ( k=0; k<grp->nitem_int_node_comm; ++k ) {
    *( (int *) tmp) = grp->node2int_node_comm_index2[k]; tmp += sizeof(int);
  }
  /* Face communicator */
  for ( k=0; k<grp->nitem_int_face_comm; ++k ) {
    *( (int *) tmp) = grp->face2int_face_comm_index1[k]; tmp += sizeof(int);
  }
  for ( k=0; k<grp->nitem_int_face_comm; ++k ) {
    *( (int *) tmp) = grp->face2int_face_comm_index2[k]; tmp += sizeof(int);
  }

  *buffer = tmp;
}

/**
 * \param parmesh pointer toward a PMMG_pParMesh structure.
 * \param buffer pointer toward the buffer in which we pack the group
 *
 * Pack the nodal intvalues array.
 *
 */
static
void PMMG_mpipack_nodeintvalues ( PMMG_pParMesh parmesh,char **buffer ) {
  PMMG_pInt_comm int_node_comm = parmesh->int_node_comm;
  int            k;
  char           *tmp;

  tmp = *buffer;

  /** Pack intvalues array of nodal communicator */
  /* Array size */
  *( (int *) tmp) = int_node_comm->nitem; tmp += sizeof(int);

  /* Intvalues array */
  for ( k=0; k<int_node_comm->nitem; ++k ) {
    *( (int *) tmp) = int_node_comm->intvalues[k]; tmp += sizeof(int);
  }

  *buffer = tmp;
}

/**
 * \param parmesh pointer toward a PMMG_pParMesh structure.
 * \param buffer pointer toward the buffer in which we pack the group
 *
 * Pack the nodal external communicators
 *
 */
static
void PMMG_mpipack_extnodecomm ( PMMG_pParMesh parmesh,char **buffer ) {
  PMMG_pExt_comm ext_node_comm;
  int   k,i;
  char  *tmp;

  tmp = *buffer;

  /** Pack nodal external communicators */
  /* Number of external communicators */
  *( (int *) tmp) = parmesh->next_node_comm; tmp += sizeof(int);

  /* List of external communicators */
  for ( k=0; k<parmesh->next_node_comm; ++k ) {
    ext_node_comm = &parmesh->ext_node_comm[k];

    /* Color in */
    *( (int *) tmp) = ext_node_comm->color_in; tmp += sizeof(int);
    /* Color out */
    *( (int *) tmp) = ext_node_comm->color_out; tmp += sizeof(int);
    /* nitem */
    *( (int *) tmp) = ext_node_comm->nitem; tmp += sizeof(int);
    /* int_comm_index */
    for ( i=0; i<ext_node_comm->nitem; ++i ) {
      *( (int *) tmp) = ext_node_comm->int_comm_index[i]; tmp += sizeof(int);
    }
  }

  *buffer = tmp;
}

/**
 * \param grp pointer toward a PMMG_Grp structure.
 * \param buffer pointer toward the buffer in which we pack the group
 *
 * \return 1 if success, 0 if fail
 *
 * \warning the mesh prisms are not treated.
 *
 * \remark To store the correct value of a variable of type "x" into buffer, we
 * need to cast the void pointer toward the buffer into a pointer toward a
 * buffer of type "x". Then we can store the variable value by dereferencing the
 * adress of the buffer.
 *
 * Pack a group into a buffer (to allow mpi communication) and shift the buffer
 * pointer at the end of the written area.
 *
 */
int PMMG_mpipack_grp ( PMMG_pGrp grp,char **buffer ) {
  int   ier;
  char  *tmp;

  ier = 1;

  tmp = *buffer;

  if ( !grp->mesh ) {
    /* unused group */
    *( (int *) tmp ) = 0; tmp += sizeof(int);

    return ier;
  }

  /* used group */
  *( (int *) tmp) = 1; tmp += sizeof(int);

  *buffer = tmp;

  PMMG_mpipack_meshSizes(grp,buffer);

  PMMG_mpipack_infos(&(grp->mesh->info),buffer);

  PMMG_mpipack_meshArrays(grp,buffer);

  PMMG_mpipack_grpintcomm(grp,buffer);

  return ier;
}

/**
 * \param parmesh pointer toward a PMMG_ParMesh structure.
 * \param buffer pointer toward the buffer in which we pack the group
 *
 * \return 1 if success, 0 if fail
 *
 * \warning the mesh prisms are not treated.
 *
 * \remark To store the correct value of a variable of type "x" into buffer, we
 * need to cast the void pointer toward the buffer into a pointer toward a
 * buffer of type "x". Then we can store the variable value by dereferencing the
 * adress of the buffer.
 *
 * Pack a parmesh into a buffer (to allow mpi communication) and shift the
 * buffer pointer at the end of the written area. The parmesh groups must have
 * been merged before entering this function.
 *
 */
int PMMG_mpipack_parmesh ( PMMG_pParMesh parmesh ,char **buffer ) {
  PMMG_pGrp grp;
  int       ier;
  char      *tmp;

  ier = 1;

  /* Check that the groups are merged */
  assert ( parmesh->ngrp < 2 );

  tmp = *buffer;

  if ( (!parmesh->ngrp) || (!parmesh->listgrp[0].mesh) ) {
    /* unused group */
    *( (int *) tmp ) = 0; tmp += sizeof(int);

    return ier;
  }

  /* used group */
  grp = &parmesh->listgrp[0];

  *( (int *) tmp) = 1; tmp += sizeof(int);

  *buffer = tmp;

  PMMG_mpipack_meshSizes(grp,buffer);

  PMMG_mpipack_infos(&(grp->mesh->info),buffer);

  PMMG_mpipack_meshArrays(grp,buffer);

  PMMG_mpipack_nodeintvalues(parmesh,buffer);

  PMMG_mpipack_extnodecomm(parmesh,buffer);

  return ier;
}

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
 * \brief functions to unpack a group mesh from a char buffer (for mpi call)
 * \author Luca Cirrottola (Inria)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 *
 */

#include "parmmg.h"


/**
 * \param listgrp pointer toward a PMMG_Grp structure array.
 * \param igrp index of the group to handle.
 * \param buffer pointer toward the buffer in which we unpack the group
 * \param np number of point    in the mesh
 * \param ne number of elements in the mesh
 * \param xp number of boundary points in the mesh
 * \param xt number of boundary elements in the mesh
 * \param ier_mesh  1 if the mesh         is allocated, 0 otherwise
 * \param npmet number of points in the metric
 * \param ier_met   1 if the metric       is allocated, 0 otherwise
 * \param metsize size of the metric
 * \param npls number of points in the level-set
 * \param ier_ls    1 if the level-set    is allocated, 0 otherwise
 * \param lssize size of the level-set
 * \param npdisp number of points in the displacement
 * \param ier_disp  1 if the displacement is allocated, 0 otherwise
 * \param dispsize size of the displacement
 * \param nsols number of solution fields
 * \param ier_field 1 if the sol fields  are allocated, 0 otherwise
 * \param fieldsize size of the solution fields.
 * \return 0 if fail, 1 if success.
 *
 * \warning the mesh prisms are not treated.
 *
 * Upack the mesh and solutions sizes from a double buffer and shift the pointer
 * toward the buffer. Fill the ier_<x> values depending on the success of the
 * allocation of the <x> structure.
 * This function gives all the available memory to the mesh.
 *
 */
static
int PMMG_mpiunpack_meshSizes ( PMMG_pParMesh parmesh,PMMG_pGrp listgrp,int igrp,
                                char **buffer,
                                int *np,int *ne,int *xp,int *xt,
                                int *ier_mesh,int *npmet,int *ier_met,int *metsize,
                                int *npls,int *ier_ls,int *lssize,
                                int *npdisp,int *ier_disp,int *dispsize,
                                int *nsols,int* ier_field,
                                int *fieldsize ) {
  PMMG_pGrp  const grp = &listgrp[igrp];
  MMG5_pMesh mesh;
  MMG5_pSol  met,ls,disp;
  int const  nprocs = parmesh->nprocs;
  int        is,ier_grp;
  int        type[MMG5_NSOLS_MAX];
  int        ismet,isls,isdisp;
  int        ier = 1;

  /** Mesh size */
  grp->mesh  = NULL;
  grp->met   = NULL;
  grp->ls    = NULL;
  grp->disp  = NULL;
  grp->field = NULL;

  ier_grp = MMG3D_Init_mesh(MMG5_ARG_start,
                            MMG5_ARG_ppMesh,&(grp->mesh),
                            MMG5_ARG_ppMet ,&(grp->met),
                            MMG5_ARG_end);

  mesh = grp->mesh;
  met  = grp->met;

  /** Get the mesh maximal authorized memory */
  (*np) = *( (int *) *buffer); *buffer += sizeof(int);
  (*xp) = *( (int *) *buffer); *buffer += sizeof(int);
  (*ne) = *( (int *) *buffer); *buffer += sizeof(int);
  (*xt) = *( (int *) *buffer); *buffer += sizeof(int);

  if ( ier_grp ) {
    /** Set the mesh size */
    (*ier_mesh) = PMMG_grpSplit_setMeshSize( mesh,*np,*ne,0,*xp,*xt );
  }
  else ier = (*ier_mesh) = 0;

  /** Number of fields */
  (*nsols) = *( (int *) *buffer); *buffer += sizeof(int);

  /** Metric info and sizes */
  ismet     = *( (int *) *buffer); *buffer += sizeof(int);

  (*npmet)   = 0;
  (*ier_met) = 1;
  if ( ismet ) {
    *metsize = *( (int *) *buffer); *buffer += sizeof(int);
    *npmet   = *np;

    if ( ier_grp ) {
      met->size = *metsize;
      /** Metric type */
      met->type = *( (int *) *buffer); *buffer += sizeof(int);

      /** Set the metric size */
      *ier_met = MMG3D_Set_solSize(mesh,met,MMG5_Vertex,*npmet,met->type);

      ier = MG_MIN ( ier, *ier_met );
    }
    else {
      ier = *ier_met = 0;
      /** Metric type */
      *buffer += sizeof(int);
    }
  }

  /** Ls info and sizes */
  isls      = *( (int *) *buffer); *buffer += sizeof(int);

  *npls   = 0;
  *ier_ls = 1;
  if ( isls ) {
    *lssize    = *( (int *) *buffer); *buffer += sizeof(int);

    *npls = *np;
    if ( ier_grp && (!grp->ls) ) {
      PMMG_CALLOC(grp->mesh,grp->ls,1,MMG5_Sol,"ls",ier_grp=ier=0);
    }

    if ( ier_grp ) {
      ls = grp->ls;

      ls->size = *lssize;
      /** Ls type */
      ls->type = *( (int *) *buffer); *buffer += sizeof(int);

      /** Set the ls size */
      *ier_ls = MMG3D_Set_solSize(mesh,ls,MMG5_Vertex,*npls,ls->type);

      ier = MG_MIN ( ier, *ier_ls );
    }
    else if ( *npls ) {
      ier = *ier_ls = 0;
      /** ls type */
      *buffer += sizeof(int);
    }
  }

  /** Disp info and sizes */
  isdisp    = *( (int *) *buffer); *buffer += sizeof(int);

  *npdisp   = 0;
  *ier_disp = 1;

  if ( isdisp ) {
    *dispsize  = *( (int *) *buffer); *buffer += sizeof(int);

    *npdisp = *np;

    if ( ier_grp && (!grp->disp) ) {
      PMMG_CALLOC(grp->mesh,grp->disp,1,MMG5_Sol,"disp",ier_grp=ier= 0);
    }

    if ( ier_grp ) {
      disp = grp->disp;

      disp->size = *dispsize;
      /** Disp type */
      disp->type = *( (int *) *buffer); *buffer += sizeof(int);

      /** Set the disp size */
      *ier_disp = MMG3D_Set_solSize(mesh,disp,MMG5_Vertex,*npdisp,disp->type);

      ier = MG_MIN ( ier, *ier_disp );
    }
    else if ( *npdisp ) {
      ier = *ier_disp = 0;
      /** disp type */
      *buffer += sizeof(int);
    }
  }

  *ier_field = 1;
  if ( *nsols ) {
    for ( is=0; is<*nsols; ++is ) {
      /** Store fields size */
      fieldsize[is] = *( (int *) *buffer); *buffer += sizeof(int);

      /** Store fields types */
      type[is] = *( (int *) *buffer); *buffer += sizeof(int);
    }

    if ( ier_grp ) {
      /** Set the fields sizes */
      *ier_field = MMG3D_Set_solsAtVerticesSize( mesh,&grp->field,*nsols,*np,type);
    }
    else {
      ier = *ier_field = 0;
    }
  }

  return ier;
}


/**
 * \param grp pointer toward a PMMG_Grp structure.
 * \param buffer pointer toward the buffer in which we pack the group
 * \param ier pointer toward the error variable: setted to 0 if this function fail.
 * \param ier_mesh  1 if the mesh         is allocated, 0 otherwise
 * \param ier_met   1 if the metric       is allocated, 0 otherwise
 * \param ier_ls    1 if the level-set    is allocated, 0 otherwise
 * \param ier_disp  1 if the displacement is allocated, 0 otherwise
 * \param nsols number of solution fields
 * \param ier_field 1 if the sol fields  are allocated, 0 otherwise
 *
 * \return 0 if the filename overflow the \a MMG5_FILENAME_LEN_MAX value, 1 otherwise.
 *
 * Unpack the filenames and shift the buffer pointer at the end of the readed
 * area.
 *
 */
static
int PMMG_mpiunpack_filenames ( PMMG_pGrp grp,char **buffer,int *ier,int ier_mesh,
                               int ier_ls,int ier_disp,int nsols,int ier_field ) {
  const MMG5_pMesh mesh = grp->mesh;
  const MMG5_pSol  met  = grp->met;
  const MMG5_pSol  ls   = grp->ls;
  const MMG5_pSol  disp = grp->disp;
  const MMG5_pSol  field= grp->field;
  int              is;
  int              meshin_s,metin_s,meshout_s,metout_s;
  int              lsin_s,lsout_s,dispin_s,dispout_s;
  int              pslin_s,pslout_s;
  char             chaine[256];

  /* File names */
  meshin_s  = *( (int *) *buffer); *buffer += sizeof(int);
  meshout_s = *( (int *) *buffer); *buffer += sizeof(int);
  if ( meshin_s > MMG5_FILENAME_LEN_MAX || meshout_s > MMG5_FILENAME_LEN_MAX ) {
    printf("  ## Error: input or output mesh name too long. Exit");
    return 0;
  }
  if ( meshin_s ) {
    strncpy(chaine,( (char *) *buffer),meshin_s); *buffer += meshin_s * sizeof(char);
    if ( !MMG5_Set_inputMeshName( mesh, chaine ) ) { *ier=0; }
  }
  if ( meshout_s ) {
    strncpy(chaine,( (char *) *buffer),meshout_s); *buffer += meshout_s * sizeof(char);
    if ( !MMG5_Set_outputMeshName( mesh, chaine ) ) { *ier = 0; }
  }

  metin_s  = *( (int *) *buffer); *buffer += sizeof(int);
  metout_s = *( (int *) *buffer); *buffer += sizeof(int);

  if ( metin_s > MMG5_FILENAME_LEN_MAX || metout_s > MMG5_FILENAME_LEN_MAX ) {
    printf("  ## Error: input or output metric name too long. Exit");
    return 0;
  }
  if ( metin_s ) {
    strncpy(chaine,( (char *) *buffer),metin_s); *buffer += metin_s * sizeof(char);
    if ( ier_mesh ) {
      if ( !MMG5_Set_inputSolName( mesh, met, chaine ) ) { *ier = 0; }
    }
  }
  if ( metout_s ) {
    strncpy(chaine,( (char *) *buffer),metout_s); *buffer += metout_s * sizeof(char);
    if ( ier_mesh ) {
      if ( !MMG5_Set_outputSolName( mesh,met, chaine ) ) { *ier = 0; }
    }
  }

  lsin_s  = *( (int *) *buffer); *buffer += sizeof(int);
  lsout_s = *( (int *) *buffer); *buffer += sizeof(int);

  if ( lsin_s > MMG5_FILENAME_LEN_MAX || lsout_s > MMG5_FILENAME_LEN_MAX ) {
    printf("  ## Error: input or output levelset name too long. Exit");
    return 0;
  }
  if ( lsin_s ) {
    strncpy(chaine,( (char *) *buffer),lsin_s); *buffer += lsin_s * sizeof(char);
    if ( ier_ls ) {
      if ( !MMG5_Set_inputSolName( mesh, ls, chaine ) ) { *ier = 0; }
    }
  }
  if ( lsout_s ) {
    strncpy(chaine,( (char *) *buffer),lsout_s); *buffer += lsout_s * sizeof(char);
    if ( ier_ls ) {
      if ( !MMG5_Set_outputSolName( mesh,ls, chaine ) ) { *ier = 0; }
    }
  }

  dispin_s  = *( (int *) *buffer); *buffer += sizeof(int);
  dispout_s = *( (int *) *buffer); *buffer += sizeof(int);

  if ( dispin_s > MMG5_FILENAME_LEN_MAX || dispout_s > MMG5_FILENAME_LEN_MAX ) {
    printf("  ## Error: input or output displacement name too long. Exit");
    return 0;
  }
  if ( dispin_s ) {
    strncpy(chaine,( (char *) *buffer),dispin_s); *buffer += dispin_s * sizeof(char);
    if ( ier_disp ) {
      if ( !MMG5_Set_inputSolName( mesh, disp, chaine ) ) { *ier = 0; }
    }
  }
  if ( dispout_s ) {
    strncpy(chaine,( (char *) *buffer),dispout_s); *buffer += dispout_s * sizeof(char);
    if ( ier_disp ) {
      if ( !MMG5_Set_outputSolName( mesh,disp, chaine ) ) { *ier = 0; }
    }
  }

  for ( is=0; is<nsols; ++is ) {
    pslin_s  = *( (int *) *buffer); *buffer += sizeof(int);
    pslout_s = *( (int *) *buffer); *buffer += sizeof(int);

    if ( pslin_s > MMG5_FILENAME_LEN_MAX || pslout_s > MMG5_FILENAME_LEN_MAX ) {
      printf("  ## Error: input or output field name too long. Exit");
      return 0;
    }
    if ( pslin_s ) {
      strncpy(chaine,( (char *) *buffer),pslin_s); *buffer += pslin_s * sizeof(char);
      if ( ier_field && field ) {
        if ( !MMG5_Set_inputSolName( mesh, &field[is], chaine ) ) { *ier = 0; }
      }
    }
    if ( pslout_s ) {
      strncpy(chaine,( (char *) *buffer),pslout_s); *buffer += pslout_s * sizeof(char);
      if ( ier_field && field ) {
        if ( !MMG5_Set_outputSolName( mesh,&field[is], chaine ) ) { *ier = 0; }
      }
    }
  }

  return 1;
}

/**
 * \param parmesh pointer toward a PMMG_pParMesh structure.
 * \param grp pointer toward a PMMG_Grp structure.
 * \param ier pointer toward the error variable: setted to 0 if this function fail.
 * \param ier_mesh  1 if the mesh         is allocated, 0 otherwise
 * \param ier_met   1 if the metric       is allocated, 0 otherwise
 * \param ier_ls    1 if the level-set    is allocated, 0 otherwise
 * \param ier_disp  1 if the displacement is allocated, 0 otherwise
 * \param nsols number of solution fields
 * \param ier_field 1 if the sol fields  are allocated, 0 otherwise
 *
 * \return 0 if the filename overflow the \a MMG5_FILENAME_LEN_MAX value, 1 otherwise.
 *
 * Unpack the filenames and shift the buffer pointer at the end of the readed
 * area.
 *
 */
static
int PMMG_copy_filenames ( PMMG_pParMesh parmesh,PMMG_pGrp grp,int *ier,int ier_mesh,
                          int ier_ls,int ier_disp,int nsols,int ier_field ) {
  const MMG5_pMesh mesh = grp->mesh;
  const MMG5_pSol  met  = grp->met;
  const MMG5_pSol  ls   = grp->ls;
  const MMG5_pSol  disp = grp->disp;
  const MMG5_pSol  field= grp->field;
  int              is;

  if ( !ier_mesh ) {
    /* mesh and metrics are not allocated */
    *ier=0;
    return 1;
  }

  /* File names */
  if ( parmesh->meshin && *parmesh->meshin ) {
    if ( !MMG5_Set_inputMeshName ( mesh, parmesh->meshin ) )  { *ier=0; }
  }
  if ( parmesh->meshout && *parmesh->meshout ) {
    if ( !MMG5_Set_outputMeshName( mesh, parmesh->meshout ) ) { *ier = 0; }
  }

  if ( parmesh->metin && *parmesh->metin ) {
    if ( !MMG5_Set_inputSolName( mesh, met, parmesh->metin ) ) { *ier = 0; }
  }
  if ( parmesh->metout && *parmesh->metout ) {
    if ( !MMG5_Set_outputSolName( mesh,met, parmesh->metout ) ) { *ier = 0; }
  }

  if ( ier_ls ) {
    /* ls structure is allocated */
    if ( parmesh->lsin && *parmesh->lsin ) {
      if ( !MMG5_Set_inputSolName( mesh, ls, parmesh->lsin ) ) { *ier = 0; }
    }
  }

  if ( ier_disp ) {
    /* disp structure is allocated */
    if ( parmesh->dispin && *parmesh->dispin ) {
      if ( !MMG5_Set_inputSolName( mesh, disp, parmesh->dispin ) ) { *ier = 0; }
    }
  }

  if ( ier_field ) {
    /* field structure is allocated */
    for ( is=0; is<nsols; ++is ) {
      if ( parmesh->fieldin ) {
        if ( !MMG5_Set_inputSolName( mesh, &field[is], parmesh->fieldin ) ) {
          *ier = 0;
        }
      }
      if ( parmesh->fieldout ) {
        if ( !MMG5_Set_outputSolName( mesh, &field[is], parmesh->fieldout ) ) {
          *ier = 0;
        }
      }
    }
  }

  return 1;
}

/**
 * \param info pointer toward a MMG5 info structure
 * \param buffer pointer toward the buffer in which we pack the group
 * \param ier pointer toward the error variable: setted to 0 if this function fail.
 * \param ier_mesh 1 if the mesh is allocated, 0 otherwise
 *
 * \return 1 if success, 0 if fail
 *
 * Unpack the info structure and shift the buffer
 * pointer at the end of the readed area.
 *
 */
static
void PMMG_mpiunpack_infos ( MMG5_Info *info,char **buffer,int *ier,int ier_mesh ) {
  int   k,nmat,npar;

  if ( ier_mesh ) {
    /** Mesh infos */
    info->dhd       = *( (double *) *buffer); *buffer += sizeof(double);
    info->hmin      = *( (double *) *buffer); *buffer += sizeof(double);
    info->hmax      = *( (double *) *buffer); *buffer += sizeof(double);
    info->hsiz      = *( (double *) *buffer); *buffer += sizeof(double);
    info->hgrad     = *( (double *) *buffer); *buffer += sizeof(double);
    info->hgradreq  = *( (double *) *buffer); *buffer += sizeof(double);
    info->hausd     = *( (double *) *buffer); *buffer += sizeof(double);

    info->delta  = *( (double *) *buffer); *buffer += sizeof(double);
    info->min[0] = *( (double *) *buffer); *buffer += sizeof(double);
    info->min[1] = *( (double *) *buffer); *buffer += sizeof(double);
    info->min[2] = *( (double *) *buffer); *buffer += sizeof(double);

    info->ls     = *( (double *) *buffer); *buffer += sizeof(double);

    info->opnbdy    = *( (int *) *buffer); *buffer += sizeof(int);
    info->renum     = *( (int *) *buffer); *buffer += sizeof(int);
    info->PROctree  = *( (int *) *buffer); *buffer += sizeof(int);
    info->npar      = *( (int *) *buffer); *buffer += sizeof(int);
    info->nmat      = *( (int *) *buffer); *buffer += sizeof(int);
    info->nsd       = *( (int *) *buffer); *buffer += sizeof(int);

    info->nreg      = *( (int8_t *) *buffer); *buffer += sizeof(int8_t);
    info->imprim    = *( (int8_t *) *buffer); *buffer += sizeof(int8_t);
    info->ddebug    = *( (int8_t *) *buffer); *buffer += sizeof(int8_t);
    info->iso       = *( (int8_t *) *buffer); *buffer += sizeof(int8_t);
    info->lag       = *( (int8_t *) *buffer); *buffer += sizeof(int8_t);
    info->parTyp    = *( (int8_t *) *buffer); *buffer += sizeof(int8_t);
    info->sethmin   = *( (int8_t *) *buffer); *buffer += sizeof(int8_t);
    info->sethmax   = *( (int8_t *) *buffer); *buffer += sizeof(int8_t);

    info->optim     = *( (uint8_t *) *buffer); *buffer += sizeof(uint8_t);
    info->optimLES  = *( (uint8_t *) *buffer); *buffer += sizeof(uint8_t);
    info->noinsert  = *( (uint8_t *) *buffer); *buffer += sizeof(uint8_t);
    info->noswap    = *( (uint8_t *) *buffer); *buffer += sizeof(uint8_t);
    info->nomove    = *( (uint8_t *) *buffer); *buffer += sizeof(uint8_t);
    info->nosurf    = *( (uint8_t *) *buffer); *buffer += sizeof(uint8_t);
    info->nosizreq  = *( (uint8_t *) *buffer); *buffer += sizeof(uint8_t);

    /* affectation of old refs in ls-mode */
    if ( info->nmat ) {

      MMG5_SAFE_CALLOC(info->mat,info->nmat,MMG5_Mat, *ier = 0);

      if ( *ier ) {
        for ( k=0; k<info->nmat; ++k ) {
          info->mat[k].dospl = *( (int8_t *) *buffer); *buffer += sizeof(int8_t);
          info->mat[k].ref   = *( (int *) *buffer); *buffer += sizeof(int);
          info->mat[k].rin   = *( (int *) *buffer); *buffer += sizeof(int);
          info->mat[k].rex   = *( (int *) *buffer); *buffer += sizeof(int);
        }
      }
      else {
        *buffer += sizeof(int8_t);
        *buffer += 3*sizeof(int);
      }
    }

    /* local parameters */
    if ( info->npar ) {

      MMG5_SAFE_CALLOC(info->par,info->npar,MMG5_Par, *ier = 0);

      if ( *ier ) {
        for ( k=0; k<info->npar; ++k ) {
          info->par[k].hmin = *( (double *) *buffer); *buffer += sizeof(double);
          info->par[k].hmax = *( (double *) *buffer); *buffer += sizeof(double);
          info->par[k].hausd = *( (double *) *buffer); *buffer += sizeof(double);
          info->par[k].ref = *( (int *) *buffer); *buffer += sizeof(int);
          info->par[k].elt = *( (int8_t *) *buffer); *buffer += sizeof(int8_t);
        }
      }
      else {
        *buffer += 3*sizeof(double) + sizeof(int) + sizeof(int8_t);
      }
    }
  }
  else {
    ier = 0;

    /** Mesh infos */
    *buffer += 12*sizeof(double);
    *buffer += 3 *sizeof(int);

    npar = *( (int *) *buffer); *buffer += sizeof(int);
    nmat = *( (int *) *buffer); *buffer += sizeof(int);

    *buffer += 8*sizeof(int8_t);
    *buffer += 7*sizeof(uint8_t);

    if ( nmat ) {
      *buffer += nmat*sizeof(int8_t);
      *buffer += nmat*sizeof(int);
      *buffer += nmat*sizeof(int);
      *buffer += nmat*sizeof(int);
    }

    /* local parameters */
    if ( npar ) {
      *buffer += npar*sizeof(double);
      *buffer += npar*sizeof(double);
      *buffer += npar*sizeof(double);
      *buffer += npar*sizeof(int);
      *buffer += npar*sizeof(int8_t);
    }
  }
}

/**
 * \param listgrp pointer toward a PMMG_Grp structure array.
 * \param igrp index of the group to handle.
 * \param buffer pointer toward the buffer in which we unpack the group
 * \param np number of point    in the mesh
 * \param ne number of elements in the mesh
 * \param xp number of boundary points in the mesh
 * \param xt number of boundary elements in the mesh
 * \param ier_mesh  1 if the mesh         is allocated, 0 otherwise
 * \param npmet number of points in the metric
 * \param ier_met   1 if the metric       is allocated, 0 otherwise
 * \param metsize size of the metric
 * \param npls number of points in the level-set
 * \param ier_ls    1 if the level-set    is allocated, 0 otherwise
 * \param lssize size of the level-set
 * \param npdisp number of points in the displacement
 * \param ier_disp  1 if the displacement is allocated, 0 otherwise
 * \param dispsize size of the displacement
 * \param nsols number of solution fields
 * \param ier_field 1 if the sol fields  are allocated, 0 otherwise
 * \param fieldsize size of the solution fields (array allocated inside this function)
 * \return 0 if fail, 1 if success.
 *
 * \warning the mesh prisms are not treated.
 *
 * Unpack the mesh and solutions (metric, ls, disp and fields) arraysand shift
 * the buffer pointer toward the end of the readed area.
 *
 */
static
int PMMG_mpiunpack_meshArrays ( PMMG_pParMesh parmesh,PMMG_pGrp listgrp,int igrp,
                                char **buffer,
                                int np,int ne,int xp,int xt,
                                int ier_mesh,int npmet,int ier_met,int metsize,
                                int npls,int ier_ls,int lssize,
                                int npdisp,int ier_disp,int dispsize,
                                int nsols,int ier_field,
                                int *fieldsize ) {
  const PMMG_pGrp  grp   = &listgrp[igrp];
  const MMG5_pMesh mesh  = grp->mesh;
  const MMG5_pSol  met   = grp->met;
  const MMG5_pSol  ls    = grp->ls;
  const MMG5_pSol  disp  = grp->disp;
  MMG5_pSol        psl;
  int const        nprocs = parmesh->nprocs;
  int              ier = 1;

  int   k,i,is;

  if ( !mesh ) ier = 0;


  if ( ier_mesh ) {
    /** Get mesh points */
    for ( k=1; k<=mesh->np; ++k ) {
      /* Coordinates */
      mesh->point[k].c[0] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->point[k].c[1] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->point[k].c[2] = *( (double *) *buffer); *buffer += sizeof(double);
      /* Tangent */
      mesh->point[k].n[0] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->point[k].n[1] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->point[k].n[2] = *( (double *) *buffer); *buffer += sizeof(double);
      /* Pointer toward the boundary entity */
      mesh->point[k].xp = *( (int *) *buffer); *buffer += sizeof(int);
      /* Ref */
      mesh->point[k].ref = *( (int *) *buffer); *buffer += sizeof(int);
      /* Tag */
      mesh->point[k].tag = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
#ifdef USE_POINTMAP
      /* Src */
      mesh->point[k].src = *( (int *) *buffer); *buffer += sizeof(int);
#endif
    }

    /** Unpack mesh boundary points */
    for ( k=1; k<=mesh->xp; ++k ) {
      /* First normal */
      mesh->xpoint[k].n1[0] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->xpoint[k].n1[1] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->xpoint[k].n1[2] = *( (double *) *buffer); *buffer += sizeof(double);
      /* Second normal */
      mesh->xpoint[k].n2[0] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->xpoint[k].n2[1] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->xpoint[k].n2[2] = *( (double *) *buffer); *buffer += sizeof(double);
      /* nnor */
      mesh->xpoint[k].nnor  = *( (int8_t *) *buffer); *buffer += sizeof(int8_t);
    }

    /** Unpack mesh elements */
    for ( k=1; k<=mesh->ne; ++k ) {
      /* Tetra vertices */
      mesh->tetra[k].v[0] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->tetra[k].v[1] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->tetra[k].v[2] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->tetra[k].v[3] = *( (int *) *buffer); *buffer += sizeof(int);
      /* Pointer toward the boundary entity */
      mesh->tetra[k].xt = *( (int *) *buffer); *buffer += sizeof(int);
      /* Ref */
      mesh->tetra[k].ref = *( (int *) *buffer); *buffer += sizeof(int);
      /* Mark */
      mesh->tetra[k].mark = *( (int *) *buffer); *buffer += sizeof(int);
      /* Tag */
      mesh->tetra[k].tag = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      /* Quality */
      mesh->tetra[k].qual = *( (double *) *buffer); *buffer += sizeof(double);
    }

    /** Unpack mesh boundary tetra */
    for ( k=1; k<=mesh->xt; ++k ) {
      /* Faces references  */
      mesh->xtetra[k].ref[0] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->xtetra[k].ref[1] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->xtetra[k].ref[2] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->xtetra[k].ref[3] = *( (int *) *buffer); *buffer += sizeof(int);
      /* Faces tags */
      mesh->xtetra[k].ftag[0] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      mesh->xtetra[k].ftag[1] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      mesh->xtetra[k].ftag[2] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      mesh->xtetra[k].ftag[3] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      /* Edges references  */
      mesh->xtetra[k].edg[0] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->xtetra[k].edg[1] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->xtetra[k].edg[2] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->xtetra[k].edg[3] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->xtetra[k].edg[4] = *( (int *) *buffer); *buffer += sizeof(int);
      mesh->xtetra[k].edg[5] = *( (int *) *buffer); *buffer += sizeof(int);
      /* Edges tags */
      mesh->xtetra[k].tag[0] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      mesh->xtetra[k].tag[1] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      mesh->xtetra[k].tag[2] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      mesh->xtetra[k].tag[3] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      mesh->xtetra[k].tag[4] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      mesh->xtetra[k].tag[5] = *( (int16_t *) *buffer); *buffer += sizeof(int16_t);
      /* Orientation of the triangles */
      mesh->xtetra[k].ori = *( (int8_t *) *buffer); *buffer += sizeof(int8_t);
    }
  }
  else {
    /* The mesh can't be allocated */
   /** Get mesh points */
    /* Coordinates */
    *buffer += np*sizeof(double);
    *buffer += np*sizeof(double);
    *buffer += np*sizeof(double);
    /* Tangent */
    *buffer += np*sizeof(double);
    *buffer += np*sizeof(double);
    *buffer += np*sizeof(double);
    /* Pointer toward the boundary entity */
    *buffer += np*sizeof(int);
    /* Ref */
    *buffer += np*sizeof(int);
    /* Tag */
    *buffer += np*sizeof(int16_t);

    /** Unpack mesh boundary points */
    /* First normal */
    *buffer += xp*sizeof(double);
    *buffer += xp*sizeof(double);
    *buffer += xp*sizeof(double);
    /* Second normal */
    *buffer += xp*sizeof(double);
    *buffer += xp*sizeof(double);
    *buffer += xp*sizeof(double);
    /* nnor */
    *buffer += xp*sizeof(int8_t);

    /** Unpack mesh elements */
    /* Tetra vertices */
    *buffer += ne*sizeof(int);
    *buffer += ne*sizeof(int);
    *buffer += ne*sizeof(int);
    *buffer += ne*sizeof(int);
    /* Pointer toward the boundary entity */
    *buffer += ne*sizeof(int);
    /* Ref */
    *buffer += ne*sizeof(int);
    /* Tag */
    *buffer += ne*sizeof(int16_t);

    /** Unpack mesh boundary tetra */
    /* Faces references  */
    *buffer += xt*sizeof(int);
    *buffer += xt*sizeof(int);
    *buffer += xt*sizeof(int);
    *buffer += xt*sizeof(int);
    /* Faces tags */
    *buffer += xt*sizeof(int16_t);
    *buffer += xt*sizeof(int16_t);
    *buffer += xt*sizeof(int16_t);
    *buffer += xt*sizeof(int16_t);
    /* Edges references  */
    *buffer += xt*sizeof(int);
    *buffer += xt*sizeof(int);
    *buffer += xt*sizeof(int);
    *buffer += xt*sizeof(int);
    *buffer += xt*sizeof(int);
    *buffer += xt*sizeof(int);
    /* Edges tags */
    *buffer += xt*sizeof(int16_t);
    *buffer += xt*sizeof(int16_t);
    *buffer += xt*sizeof(int16_t);
    *buffer += xt*sizeof(int16_t);
    *buffer += xt*sizeof(int16_t);
    *buffer += xt*sizeof(int16_t);
    /* Triangle orientation */
    *buffer += xt*sizeof(int8_t);
  }

  /** Unpack metric */
  if( npmet ) {/* only if the metrics size is non null, i.e. not default metrics */
    if ( ier_met ) {
      for ( k=1; k<=np; ++k ) {
        for ( i=0; i<metsize; ++i ) {
          met->m[metsize*k + i] = *( (double *) *buffer);
          *buffer += sizeof(double);
        }
      }
    }
    else {
      /* The metric array can't be allocated */
      *buffer += np*metsize*sizeof(double);
    }
  }

  /** Unpack ls */
  if( npls ) {/* only if the ls size is non null */
    if ( ier_ls ) {
      for ( k=1; k<=np; ++k ) {
        for ( i=0; i<lssize; ++i ) {
          ls->m[lssize*k + i] = *( (double *) *buffer);
          *buffer += sizeof(double);
        }
      }
    }
    else {
      /* The ls array can't be allocated */
      *buffer += np*lssize*sizeof(double);
    }
  }

  /** Unpack metric */
  if( npdisp ) {/* only if the displacement size is non null */
    if ( ier_disp ) {
      for ( k=1; k<=np; ++k ) {
        for ( i=0; i<dispsize; ++i ) {
          disp->m[dispsize*k + i] = *( (double *) *buffer);
          *buffer += sizeof(double);
        }
      }
    }
    else {
      /* The metric array can't be allocated */
      *buffer += np*dispsize*sizeof(double);
    }
  }

  if ( nsols ) {
    if ( ier_field ) {
      for ( is=0; is<nsols; ++is ) {
        psl = &grp->field[is];
        for ( k=1; k<=np; ++k ) {
          for ( i=0; i<psl->size; ++i ) {
            psl->m[psl->size*k + i] = *( (double *) *buffer);
            *buffer += sizeof(double);
          }
        }
      }
    }
    else {
      for ( is=0; is<nsols; ++is ) {
        *buffer += np*fieldsize[is]*sizeof(double);
      }
    }
  }

  return ier;
}

/**
 * \param parmesh pointer toward a parmesh structure.
 * \param grp pointer toward a PMMG_Grp structure.
 * \param buffer pointer toward the buffer in which we pack the group
 * \param ier pointer toward the error value (setted to 0 if we fail)
 *
 * Unpack the internal communicators of the group and move the buffer pointer at
 * the end of the readed area. Set the parmesh memory to the minimal amount needed.
 *
 */
static
void PMMG_mpiunpack_grpintcomm ( PMMG_pParMesh parmesh,PMMG_pGrp grp,
                                 char **buffer,int *ier ) {

  int   k,ier_comm;

  /** Unpack communicators */
  ier_comm = 1;

  /** Communicator sizes */
  grp->nitem_int_node_comm = *( (int *) *buffer); *buffer += sizeof(int);
  grp->nitem_int_face_comm = *( (int *) *buffer); *buffer += sizeof(int);

  /* Node communicator */
  PMMG_MALLOC(parmesh,grp->node2int_node_comm_index1,grp->nitem_int_node_comm,
              int,"node2int_node_comm_index1",*ier = ier_comm = 0);
  PMMG_MALLOC(parmesh,grp->node2int_node_comm_index2,grp->nitem_int_node_comm,
              int,"node2int_node_comm_index1",*ier = ier_comm = 0);

  if ( ier_comm ) {
    for ( k=0; k<grp->nitem_int_node_comm; ++k ) {
      grp->node2int_node_comm_index1[k] = *( (int *) *buffer);
      *buffer += sizeof(int);
    }
    for ( k=0; k<grp->nitem_int_node_comm; ++k ) {
      grp->node2int_node_comm_index2[k] = *( (int *) *buffer);
      *buffer += sizeof(int);
    }
  }
  else {
    /* The communicators arrays can't be allocated */
    *buffer += 2*grp->nitem_int_node_comm*sizeof(int);
  }

  ier_comm=1;
  /* Face communicator */
  PMMG_MALLOC(parmesh,grp->face2int_face_comm_index1,grp->nitem_int_face_comm,
              int,"face2int_face_comm_index1",*ier = ier_comm = 0);
  PMMG_MALLOC(parmesh,grp->face2int_face_comm_index2,grp->nitem_int_face_comm,
              int,"face2int_face_comm_index1",*ier = ier_comm = 0);


  if ( ier_comm ) {
    for ( k=0; k<grp->nitem_int_face_comm; ++k ) {
      grp->face2int_face_comm_index1[k] = *( (int *) *buffer);
      *buffer += sizeof(int);
    }
    for ( k=0; k<grp->nitem_int_face_comm; ++k ) {
      grp->face2int_face_comm_index2[k] = *( (int *) *buffer);
      *buffer += sizeof(int);
    }
  }
  else {
    /* The communicators arrays can't be allocated */
    *buffer += 2*grp->nitem_int_face_comm*sizeof(int);
  }
}

/**
 * \param parmesh pointer toward a parmesh structure.
 * \param int_node_comm pointer toward an internal communicator.
 * \param buffer pointer toward the buffer in which we pack the group
 * \param ier pointer toward the error value (setted to 0 if we fail)
 *
 * Unpack the nodal intvalues array
 *
 */
static
void PMMG_mpiunpack_nodeintvalues ( PMMG_pParMesh parmesh,
                                    PMMG_pInt_comm int_node_comm,
                                    char **buffer,int *ier ) {
  int            k,ier_comm;

  /** Unpack intvalues array */
  ier_comm = 1;

  /** Size of intvalue array */
  int_node_comm->nitem = *( (int *) *buffer); *buffer += sizeof(int);

  /* Allocation of intvalues */
  PMMG_MALLOC(parmesh,int_node_comm->intvalues,int_node_comm->nitem,
              int,"int_node_comm->intvalues",*ier = ier_comm = 0);

  if ( ier_comm ) {
    for ( k=0; k<int_node_comm->nitem; ++k ) {
      int_node_comm->intvalues[k] = *( (int *) *buffer);
      *buffer += sizeof(int);
    }
  }
  else {
    /* The intvalue array can't be allocated */
    *buffer += int_node_comm->nitem*sizeof(int);
  }
}

/**
 * \param parmesh pointer toward a parmesh structure.
 * \param next_node_comm pointer toward the number of external communicators.
 * \param ext_node_comm pointer toward an array of external communicators.
 * \param buffer pointer toward the buffer in which we pack the group
 * \param ier pointer toward the error value (setted to 0 if we fail)
 *
 * Unpack the nodal external communicators.
 *
 */
static
void PMMG_mpiunpack_extnodecomm ( PMMG_pParMesh parmesh,
                                  int *next_node_comm,
                                  PMMG_pExt_comm *ext_node_comm,
                                  char **buffer,int *ier ) {

  int k,i,ier_comm,nitem_unread;

  /** Unpack communicators */
  ier_comm = 1;

  /** Communicator sizes */
  *next_node_comm = *( (int *) *buffer); *buffer += sizeof(int);

  /* External communicators */
  PMMG_MALLOC(parmesh,*ext_node_comm,*next_node_comm,
              PMMG_Ext_comm,"ext_node_comm",*ier = ier_comm = 0);

  if ( ier_comm ) {
    for ( k=0; k<*next_node_comm; ++k ) {
      (*ext_node_comm)[k].color_in  = *( (int *) *buffer); *buffer += sizeof(int);
      (*ext_node_comm)[k].color_out = *( (int *) *buffer); *buffer += sizeof(int);
      (*ext_node_comm)[k].nitem     = *( (int *) *buffer); *buffer += sizeof(int);

      if ( ier_comm ) {
        PMMG_MALLOC(parmesh,(*ext_node_comm)[k].int_comm_index,
                    (*ext_node_comm)[k].nitem,
                    int,"ext_node_comm->int_comm_index",*ier = ier_comm = 0);
      }
      if ( ier_comm ) {
        for ( i=0; i<(*ext_node_comm)[k].nitem; ++i ) {
          (*ext_node_comm)[k].int_comm_index[i] = *( (int *) *buffer);
          *buffer += sizeof(int);
        }
      }
      else {
        /* Skip unreadable items */
        *buffer += (*ext_node_comm)[k].nitem*sizeof(int);
      }
    }
  }
  else {
    for ( k=0 ; k<*next_node_comm; ++k ) {
      /* Travel the remaining communicator to skip unreadable items */
      *buffer += 2*sizeof(int); //colors

      nitem_unread = *( (int *) *buffer); *buffer += sizeof(int);
      *buffer += nitem_unread*sizeof(int); // ext_node_comm->int_comm_index
    }
  }

}


/**
 * \param parmesh pointer toward a parmesh structure.
 * \param listgrp pointer toward a PMMG_Grp structure array.
 * \param igrp index of the group to handle.
 * \param buffer pointer toward the buffer in which we unpack the group
 *
 * \return 0 if fail, 1 otherwise
 *
 * \warning the mesh prisms are not treated.
 *
 * Upack a group from a double buffer and shift the pointer toward the buffer to
 * point to the next group stored in the buffer.
 *
 * \remark  To get the correct value of a variable
 * of type "x", we need to cast the void pointer toward the buffer into a
 * pointer toward a buffer of type "x". Then we can get the variable value by
 * dereferencing the adress of the buffer.
 *
 */
int PMMG_mpiunpack_grp ( PMMG_pParMesh parmesh,PMMG_pGrp listgrp,int igrp,char **buffer ) {
  PMMG_pGrp const grp = &listgrp[igrp];
  int        ier,ier_mesh,ier_met,ier_ls,ier_disp,ier_field;
  int        np,npmet,npdisp,npls,xp,ne,xt;
  int        metsize,lssize,dispsize,fieldsize[MMG5_NSOLS_MAX];
  int        nsols,used;

  ier = 1;

  used = *( (int *) *buffer ); *buffer += sizeof(int);
  if ( !used ) {
    /* unused group */
    return ier;
  }

  ier = PMMG_mpiunpack_meshSizes ( parmesh,listgrp,igrp,buffer,&np,&ne,&xp,&xt,
                             &ier_mesh,&npmet,&ier_met,&metsize,
                             &npls,&ier_ls,&lssize,&npdisp,&ier_disp,&dispsize,
                             &nsols,&ier_field,fieldsize );

  PMMG_copy_filenames ( parmesh,grp,&ier,ier_mesh,ier_ls,ier_disp,nsols,ier_field );

  PMMG_mpiunpack_infos(&(grp->mesh->info),buffer,&ier,ier_mesh);

  ier = PMMG_mpiunpack_meshArrays( parmesh,listgrp,igrp,buffer,np,ne,xp,xt,ier_mesh,
                             npmet,ier_met,metsize,npls,ier_ls,lssize,npdisp,
                             ier_disp,dispsize,nsols,ier_field,fieldsize );


  PMMG_mpiunpack_grpintcomm ( parmesh,grp,buffer,&ier);

  return ier;
}

/**
 * \param parmesh pointer toward a parmesh structure.
 * \param listgrp pointer toward a PMMG_Grp structure array.
 * \param igrp index of the group to handle.
 * \param int_node_comm pointer toward an internal communicator.
 * \param next_node_comm pointer toward the number of external communicators.
 * \param ext_node_comm pointer toward an array of external communicators.
 * \param buffer pointer toward the buffer in which we unpack the group
 *
 * \return 0 if fail, 1 otherwise
 *
 * \warning the mesh prisms are not treated.
 *
 * Upack a  from a double buffer and shift the pointer toward the buffer to
 * point to the next group stored in the buffer.
 *
 * \remark  To get the correct value of a variable
 * of type "x", we need to cast the void pointer toward the buffer into a
 * pointer toward a buffer of type "x". Then we can get the variable value by
 * dereferencing the adress of the buffer.
 *
 */
int PMMG_mpiunpack_parmesh ( PMMG_pParMesh parmesh,PMMG_pGrp listgrp,int igrp,
                             PMMG_pInt_comm int_node_comm,
                             int *next_node_comm,
                             PMMG_pExt_comm *ext_node_comm,
                             char **buffer ) {
  PMMG_pGrp const grp = &listgrp[igrp];
  int        ier,ier_mesh,ier_met,ier_ls,ier_disp,ier_field;
  int        np,npmet,npdisp,npls,xp,ne,xt;
  int        metsize,lssize,dispsize,fieldsize[MMG5_NSOLS_MAX];
  int        nsols,used;

  ier = 1;

  used = *( (int *) *buffer ); *buffer += sizeof(int);
  if ( !used ) {
    /* unused group */
    return ier;
  }

  ier = PMMG_mpiunpack_meshSizes ( parmesh,listgrp,igrp,buffer,&np,&ne,&xp,&xt,
                             &ier_mesh,&npmet,&ier_met,&metsize,
                             &npls,&ier_ls,&lssize,&npdisp,&ier_disp,&dispsize,
                             &nsols,&ier_field,fieldsize );

  PMMG_mpiunpack_infos(&(grp->mesh->info),buffer,&ier,ier_mesh);

  ier = PMMG_mpiunpack_meshArrays( parmesh,listgrp,igrp,buffer,np,ne,xp,xt,ier_mesh,
                             npmet,ier_met,metsize,npls,ier_ls,lssize,npdisp,
                             ier_disp,dispsize,nsols,ier_field,fieldsize );


  PMMG_mpiunpack_nodeintvalues ( parmesh,int_node_comm,buffer,&ier);

  PMMG_mpiunpack_extnodecomm ( parmesh,next_node_comm,ext_node_comm,buffer,
                               &ier);

  return ier;
}

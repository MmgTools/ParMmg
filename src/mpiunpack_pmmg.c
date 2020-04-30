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
 * \param grp pointer toward a PMMG_Grp structure.
 * \param buffer pointer toward the buffer in which we unpack the group
 * \param memAv pointer toward the available memory whose value is updated.
 * \param ier pointer toward the error variable: setted to 0 if this function fail.
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
void PMMG_mpiunpack_meshSizes ( PMMG_pGrp grp,char **buffer,
                                size_t *memAv,int *ier,
                                int *np,int *ne,int *xp,int *xt,
                                int *ier_mesh,int *npmet,int *ier_met,int *metsize,
                                int *npls,int *ier_ls,int *lssize,
                                int *npdisp,int *ier_disp,int *dispsize,
                                int *nsols,int* ier_field,
                                int **fieldsize ) {
  MMG5_pMesh mesh;
  MMG5_pSol  met,ls,disp;
  int        is,ier_grp,ier_field0;
  int        *type;
  int        ismet,isls,isdisp;

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
    /* Give all the available memory to the mesh */
    mesh->memMax = *memAv;

    /** Set the mesh size */
    (*ier_mesh) = PMMG_grpSplit_setMeshSize( mesh,*np,*ne,0,*xp,*xt );
  }
  else (*ier) = (*ier_mesh) = 0;

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

      *ier = MG_MIN ( *ier, *ier_met );
    }
    else {
      *ier = *ier_met = 0;
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
      PMMG_CALLOC(grp->mesh,grp->ls,1,MMG5_Sol,"ls",ier_grp=(*ier)=0);
    }

    if ( ier_grp ) {
      ls = grp->ls;

      ls->size = *lssize;
      /** Ls type */
      ls->type = *( (int *) *buffer); *buffer += sizeof(int);

      /** Set the ls size */
      *ier_ls = MMG3D_Set_solSize(mesh,ls,MMG5_Vertex,*npls,ls->type);

      *ier = MG_MIN ( *ier, *ier_ls );
    }
    else if ( *npls ) {
      *ier = *ier_ls = 0;
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
      PMMG_CALLOC(grp->mesh,grp->disp,1,MMG5_Sol,"disp",ier_grp=(*ier)= 0);
    }

    if ( ier_grp ) {
      disp = grp->disp;

      disp->size = *dispsize;
      /** Disp type */
      disp->type = *( (int *) *buffer); *buffer += sizeof(int);

      /** Set the disp size */
      *ier_disp = MMG3D_Set_solSize(mesh,disp,MMG5_Vertex,*npdisp,disp->type);

      *ier = MG_MIN ( *ier, *ier_disp );
    }
    else if ( *npdisp ) {
      *ier = *ier_disp = 0;
      /** disp type */
      *buffer += sizeof(int);
    }
  }

  ier_field0 = *ier_field = 1;
  *fieldsize = type = NULL;
  if ( *nsols ) {
    MMG5_SAFE_MALLOC (       type,*nsols,int,ier_field0 = *ier_field = 0 );
    MMG5_SAFE_MALLOC ( *fieldsize,*nsols,int,ier_field0 = *ier_field = 0 );

    if ( ier_grp && (!grp->field) ) {
      PMMG_CALLOC(grp->mesh,grp->field,*nsols,MMG5_Sol,"field",ier_grp=(*ier)=0);
    }

    if ( ier_field0 ) {
      for ( is=0; is<*nsols; ++is ) {
        /** Store fields size */
        (*fieldsize)[is] = *( (int *) *buffer); *buffer += sizeof(int);

        /** Store fields types */
        type[is] = *( (int *) *buffer); *buffer += sizeof(int);
      }
    }
    else {
      /** Skip field type and sizes */
      *buffer += 2*(*nsols)*sizeof(int);
    }

    if ( ier_grp && *ier_field ) {
      /** Set the fields sizes */
      *ier_field = MMG3D_Set_solsAtVerticesSize( mesh,&grp->field,*nsols,*np,type);
    }
    else {
      *ier = *ier_field = 0;
    }
    if ( type ) {
      MMG5_SAFE_FREE(type);
    }
  }
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
  else {
    strncpy(chaine,( (char *) *buffer),meshin_s); *buffer += meshin_s * sizeof(char);
    if ( !MMG5_Set_inputMeshName( mesh, chaine ) ) { *ier=0; }
    strncpy(chaine,( (char *) *buffer),meshout_s); *buffer += meshout_s * sizeof(char);
    if ( !MMG5_Set_outputMeshName( mesh, chaine ) ) { *ier = 0; }
  }

  metin_s  = *( (int *) *buffer); *buffer += sizeof(int);
  metout_s = *( (int *) *buffer); *buffer += sizeof(int);

  if ( metin_s > MMG5_FILENAME_LEN_MAX || metout_s > MMG5_FILENAME_LEN_MAX ) {
    printf("  ## Error: input or output metric name too long. Exit");
    return 0;
  }
  else if ( metin_s || metout_s ) {
    /* Either both names or none */
    assert ( metin_s && metout_s );

    strncpy(chaine,( (char *) *buffer),metin_s); *buffer += metin_s * sizeof(char);
    if ( ier_mesh ) {
      if ( !MMG5_Set_inputSolName( mesh, met, chaine ) ) { *ier = 0; }
    }
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
  else if ( lsin_s || lsout_s ) {
    /* Either both names or none */
    assert ( lsin_s && lsout_s );

    strncpy(chaine,( (char *) *buffer),lsin_s); *buffer += lsin_s * sizeof(char);
    if ( ier_ls ) {
      if ( !MMG5_Set_inputSolName( mesh, ls, chaine ) ) { *ier = 0; }
    }
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
  else if ( dispin_s || dispout_s ) {
    /* Either both names or none */
    assert ( dispin_s && dispout_s );

    strncpy(chaine,( (char *) *buffer),dispin_s); *buffer += dispin_s * sizeof(char);
    if ( ier_disp ) {
      if ( !MMG5_Set_inputSolName( mesh, disp, chaine ) ) { *ier = 0; }
    }
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
    else if ( pslin_s || pslout_s ) {
      /* Either both names or none */
      assert ( pslin_s && pslout_s );

      strncpy(chaine,( (char *) *buffer),pslin_s); *buffer += pslin_s * sizeof(char);
      if ( ier_field && field ) {
        if ( !MMG5_Set_inputSolName( mesh, &field[is], chaine ) ) { *ier = 0; }
      }
      strncpy(chaine,( (char *) *buffer),pslout_s); *buffer += pslout_s * sizeof(char);
      if ( ier_field && field ) {
        if ( !MMG5_Set_outputSolName( mesh,&field[is], chaine ) ) { *ier = 0; }
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

    info->npar      = *( (int *) *buffer); *buffer += sizeof(int);
    info->opnbdy    = *( (int *) *buffer); *buffer += sizeof(int);
    info->renum     = *( (int *) *buffer); *buffer += sizeof(int);
    info->PROctree  = *( (int *) *buffer); *buffer += sizeof(int);
    info->nmat      = *( (int *) *buffer); *buffer += sizeof(int);

    info->nreg      = *( (char *) *buffer); *buffer += sizeof(char);
    info->imprim    = *( (char *) *buffer); *buffer += sizeof(char);
    info->ddebug    = *( (char *) *buffer); *buffer += sizeof(char);
    info->iso       = *( (char *) *buffer); *buffer += sizeof(char);
    info->lag       = *( (char *) *buffer); *buffer += sizeof(char);
    info->parTyp    = *( (char *) *buffer); *buffer += sizeof(char);
    info->optim     = *( (char *) *buffer); *buffer += sizeof(char);
    info->optimLES  = *( (char *) *buffer); *buffer += sizeof(char);
    info->noinsert  = *( (char *) *buffer); *buffer += sizeof(char);
    info->noswap    = *( (char *) *buffer); *buffer += sizeof(char);
    info->nomove    = *( (char *) *buffer); *buffer += sizeof(char);
    info->nosurf    = *( (char *) *buffer); *buffer += sizeof(char);

    /* affectation of old refs in ls-mode */
    if ( info->nmat ) {

      MMG5_SAFE_CALLOC(info->mat,info->nmat,MMG5_Mat, *ier = 0);

      if ( *ier ) {
        for ( k=0; k<info->nmat; ++k ) {
          info->mat[k].dospl = *( (char *) *buffer); *buffer += sizeof(char);
          info->mat[k].ref   = *( (int *) *buffer); *buffer += sizeof(int);
          info->mat[k].rin   = *( (int *) *buffer); *buffer += sizeof(int);
          info->mat[k].rex   = *( (int *) *buffer); *buffer += sizeof(int);
        }
      }
      else {
        *buffer += sizeof(char);
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
          info->par[k].elt = *( (char *) *buffer); *buffer += sizeof(char);
        }
      }
      else {
        *buffer += 3*sizeof(double) + sizeof(int) + sizeof(char);
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

    *buffer += 13*sizeof(char);

    if ( nmat ) {
      *buffer += nmat*sizeof(char);
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
      *buffer += npar*sizeof(char);
    }
  }
}

/**
 * \param grp pointer toward a PMMG_Grp structure.
 * \param buffer pointer toward the buffer in which we unpack the group
 * \param memAv pointer toward the available memory whose value is updated.
 * \param ier pointer toward the error variable: setted to 0 if this function fail.
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
 *
 * \warning the mesh prisms are not treated.
 *
 * Unpack the mesh and solutions (metric, ls, disp and fields) arraysand shift
 * the buffer pointer toward the end of the readed area.
 *
 */
static
void PMMG_mpiunpack_meshArrays ( PMMG_pGrp grp,char **buffer,
                                 size_t *memAv,int *ier,
                                 int np,int ne,int xp,int xt,
                                 int ier_mesh,int npmet,int ier_met,int metsize,
                                 int npls,int ier_ls,int lssize,
                                 int npdisp,int ier_disp,int dispsize,
                                 int nsols,int ier_field,
                                 int **fieldsize ) {
  const MMG5_pMesh mesh  = grp->mesh;
  const MMG5_pSol  met   = grp->met;
  const MMG5_pSol  ls    = grp->ls;
  const MMG5_pSol  disp  = grp->disp;
  MMG5_pSol        psl;

  int   k,i,is;

  if ( mesh ) {
    /* Use exactly the amount of needed memory for this mesh and metric */
    mesh->memMax = mesh->memCur;

    /* Update the available memory count */
    *memAv -= mesh->memMax;
  }
  else {
    *ier = 0;
  }

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
    }

    /** Pack mesh boundary points */
    for ( k=1; k<=mesh->xp; ++k ) {
      /* First normal */
      mesh->xpoint[k].n1[0] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->xpoint[k].n1[1] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->xpoint[k].n1[2] = *( (double *) *buffer); *buffer += sizeof(double);
      /* Second normal */
      mesh->xpoint[k].n2[0] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->xpoint[k].n2[1] = *( (double *) *buffer); *buffer += sizeof(double);
      mesh->xpoint[k].n2[2] = *( (double *) *buffer); *buffer += sizeof(double);
    }

    /** Pack mesh elements */
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

    /** Pack mesh boundary tetra */
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

    /** Pack mesh boundary points */
    /* First normal */
    *buffer += xp*sizeof(double);
    *buffer += xp*sizeof(double);
    *buffer += xp*sizeof(double);
    /* Second normal */
    *buffer += xp*sizeof(double);
    *buffer += xp*sizeof(double);
    *buffer += xp*sizeof(double);

    /** Pack mesh elements */
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

    /** Pack mesh boundary tetra */
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

  if ( nsols && *fieldsize ) {
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
        *buffer += np*(*fieldsize)[is]*sizeof(double);
      }
    }
    MMG5_SAFE_FREE(*fieldsize);
  }

}

/**
 * \param parmesh pointer toward a parmesh structure.
 * \param grp pointer toward a PMMG_Grp structure.
 * \param buffer pointer toward the buffer in which we pack the group
 * \param memAv pointer toward the available memory whose value is updated.
 * \param ier pointer toward the error value (setted to 0 if we fail)
 * \param ptrok 1 if the buffer pointer is at the good position (ie fieldsize has been allocated)
 *
 * Unpack the internal communicators of the group and move the buffer pointer at
 * the end of the readed area. Set the parmesh memory to the minimal amount needed.
 *
 */
static
void PMMG_mpiunpack_grpintcomm ( PMMG_pParMesh parmesh,PMMG_pGrp grp,
                                 char **buffer,size_t *memAv,
                                 int *ier,int ptrok ) {

  int   k,ier_comm;

  /* Give all the available mem to the communicators */
  parmesh->memMax = *memAv;

  if ( !ptrok ) {
    /* We were not able to read the field sizes so we can't know the span to
     * be skipped and we can't continue to read the mesh */
    ier_comm = 0;
  }
  else {
    /** Pack communicators */
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

    /* Use the minimal memory needed */
    parmesh->memMax = parmesh->memCur;
    *memAv         -= parmesh->memMax;

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

}

/**
 * \param parmesh pointer toward a parmesh structure.
 * \param grp pointer toward a PMMG_Grp structure.
 * \param buffer pointer toward the buffer in which we unpack the group
 * \param memAv pointer toward the available memory whose value is updated.
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
int PMMG_mpiunpack_grp ( PMMG_pParMesh parmesh,PMMG_pGrp grp,char **buffer,
                         size_t *memAv) {
  int        ier,ier_mesh,ier_met,ier_ls,ier_disp,ier_field;
  int        np,npmet,npdisp,npls,xp,ne,xt;
  int        metsize,lssize,dispsize,*fieldsize;
  int        nsols,used,ptrok;

  ier = 1;

  used = *( (int *) *buffer ); *buffer += sizeof(int);
  if ( !used ) {
    /* unused group */
    return ier;
  }

  PMMG_mpiunpack_meshSizes ( grp,buffer,memAv,&ier,&np,&ne,&xp,&xt,
                             &ier_mesh,&npmet,&ier_met,&metsize,
                             &npls,&ier_ls,&lssize,&npdisp,&ier_disp,&dispsize,
                             &nsols,&ier_field,&fieldsize );

  ptrok =  nsols ? (fieldsize!=NULL) : 1;

  if ( !PMMG_mpiunpack_filenames ( grp,buffer,&ier,ier_mesh,ier_ls,ier_disp,
                                   nsols,ier_field ) ) {
    return 0;
  }

  PMMG_mpiunpack_infos(&(grp->mesh->info),buffer,&ier,ier_mesh);

  PMMG_mpiunpack_meshArrays( grp,buffer,memAv,&ier,np,ne,xp,xt,ier_mesh,
                             npmet,ier_met,metsize,npls,ier_ls,lssize,npdisp,
                             ier_disp,dispsize,nsols,ier_field,&fieldsize );


  PMMG_mpiunpack_grpintcomm ( parmesh,grp,buffer,memAv,&ier,ptrok);

  return ier;
}

/**
 * \param parmesh pointer toward a parmesh structure.
 * \param grp pointer toward a PMMG_Grp structure.
 * \param buffer pointer toward the buffer in which we unpack the group
 * \param memAv pointer toward the available memory whose value is updated.
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
int PMMG_mpiunpack_grp4finalmerge ( PMMG_pParMesh parmesh,PMMG_pGrp grp,char **buffer,
                                    size_t *memAv) {
  int        ier,ier_mesh,ier_met,ier_ls,ier_disp,ier_field;
  int        np,npmet,npdisp,npls,xp,ne,xt;
  int        metsize,lssize,dispsize,*fieldsize;
  int        nsols,used,ptrok;

  ier = 1;

  used = *( (int *) *buffer ); *buffer += sizeof(int);
  if ( !used ) {
    /* unused group */
    return ier;
  }

  PMMG_mpiunpack_meshSizes ( grp,buffer,memAv,&ier,&np,&ne,&xp,&xt,
                             &ier_mesh,&npmet,&ier_met,&metsize,
                             &npls,&ier_ls,&lssize,&npdisp,&ier_disp,&dispsize,
                             &nsols,&ier_field,&fieldsize );

  ptrok =  nsols ? (fieldsize!=NULL) : 1;

#ifndef NDEBUG
  if ( !PMMG_mpiunpack_filenames ( grp,buffer,&ier,ier_mesh,ier_ls,ier_disp,
                                   nsols,ier_field ) ) {
    return 0;
  }
#endif

  PMMG_mpiunpack_meshArrays( grp,buffer,memAv,&ier,np,ne,xp,xt,ier_mesh,
                             npmet,ier_met,metsize,npls,ier_ls,lssize,npdisp,
                             ier_disp,dispsize,nsols,ier_field,&fieldsize );


  PMMG_mpiunpack_grpintcomm ( parmesh,grp,buffer,memAv,&ier,ptrok);

  return ier;
}

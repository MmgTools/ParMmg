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
 * \file inout_pmmg.c
 * \brief io for the parmmg software
 * \author Algiane Froehly (InriaSoft)
 * \version 1
 * \date 07 2018
 * \copyright GNU Lesser General Public License.
 *
 * input/outputs for parmmg.
 *
 */

#include "parmmg.h"

/**
 * \param n integer for which we want to know the number of digits
 *
 * \return the number of digits of n.
 *
 */
static inline
int PMMG_count_digits(int n) {

  int count = 0;
  while (n != 0) {
    n /= 10;
    ++count;
  }

  return count;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param inm pointer to the mesh file.
 * \param bin binary (1) or ascii (0) file.
 * \param iswp perform byte swapping (1) or not (0).
 * \param pos position of the communicators in the binary file.
 * \param ncomm number of communicators to read.
 * \param nitem_comm pointer to the nb of items in each communicator.
 * \param color pointer to the color of each communicator.
 * \param idx_loc pointer to the local indices of entities in each communicator.
 * \param idx_glo pointer to the global indices of entities in each communicator.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Load parallel communicator in Medit format (only one group per process is
 * allowed).
 *
 */
int PMMG_loadCommunicator( PMMG_pParMesh parmesh,FILE *inm,int bin,int iswp,
                           int pos,int ncomm,int *nitem_comm,int *color,
                           int **idx_loc,int **idx_glo ) {
  char chaine[MMG5_FILESTR_LGTH],strskip[MMG5_FILESTR_LGTH];
  int binch,bpos;
  int *inxt;
  int ntot,k,idxl,idxg,icomm,i;

  PMMG_CALLOC(parmesh,inxt,ncomm,int,"inxt",return 0);

  rewind(inm);
  fseek(inm,pos,SEEK_SET);
  /* Read color and nb of items */
  ntot = 0;
  if(!bin) {
    for( icomm = 0; icomm < ncomm; icomm++ ) {
      MMG_FSCANF(inm,"%d %d",&color[icomm],&nitem_comm[icomm]);
      ntot += nitem_comm[icomm];
    }
  }
  else {
    for( icomm = 0; icomm < ncomm; icomm++ ) {
      MMG_FREAD(&k,MMG5_SW,1,inm);
      if(iswp) k=MMG5_swapbin(k);
      color[icomm] = k;
      MMG_FREAD(&k,MMG5_SW,1,inm);
      if(iswp) k=MMG5_swapbin(k);
      nitem_comm[icomm] = k;
      ntot += nitem_comm[icomm];
    }
  }
  /* Allocate indices arrays */
  for( icomm = 0; icomm < ncomm; icomm++ ) {
    PMMG_CALLOC(parmesh,idx_loc[icomm],nitem_comm[icomm],int,
                "idx_loc",return 0);
    PMMG_CALLOC(parmesh,idx_glo[icomm],nitem_comm[icomm],int,
                "idx_glo",return 0);
  }

  rewind(inm);
  if (!bin) {
    strcpy(chaine,"D");
    while(fscanf(inm,"%127s",&chaine[0])!=EOF && strncmp(chaine,"End",strlen("End")) ) {
      if ( chaine[0] == '#' ) {
        fgets(strskip,MMG5_FILESTR_LGTH,inm);
        continue;
      }
      if( parmesh->info.API_mode == PMMG_APIDISTRIB_faces ) {
        if(!strncmp(chaine,"ParallelCommunicatorTriangles",strlen("ParallelCommunicatorTriangles"))) {
          pos = ftell(inm);
          break;
        }
      } else if( parmesh->info.API_mode == PMMG_APIDISTRIB_nodes ) {
         if(!strncmp(chaine,"ParallelCommunicatorVertices",strlen("ParallelCommunicatorVertices"))) {
          pos = ftell(inm);
          break;
        }
      }
    }
  } else { //binary file
    while(fread(&binch,MMG5_SW,1,inm)!=0 && binch!=54 ) {
      if(iswp) binch=MMG5_swapbin(binch);
      if(binch==54) break;
      if(!ncomm && binch==72) { // ParallelCommunicatorTriangles
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        pos = ftell(inm);
        break; // if parallel triangles are found, ignore parallel nodes
      } else if(!ncomm && binch==73) { // ParallelCommunicatorVertices
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        pos = ftell(inm);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else {
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
      }
    }
  }


  /* Read indices */
  if(!bin) {
    for( i = 0; i < ntot; i++ ) {
      MMG_FSCANF(inm,"%d %d %d",&idxl,&idxg,&icomm);
      idx_loc[icomm][inxt[icomm]] = idxl;
      idx_glo[icomm][inxt[icomm]] = idxg;
      inxt[icomm]++;
    }
  } else {
    for( i = 0; i < ntot; i++ ) {
      MMG_FREAD(&k,MMG5_SW,1,inm);
      if(iswp) k=MMG5_swapbin(k);
      idxl = k;
      MMG_FREAD(&k,MMG5_SW,1,inm);
      if(iswp) k=MMG5_swapbin(k);
      idxg = k;
      MMG_FREAD(&k,MMG5_SW,1,inm);
      if(iswp) k=MMG5_swapbin(k);
      icomm = k;
      idx_loc[icomm][inxt[icomm]] = idxl;
      idx_glo[icomm][inxt[icomm]] = idxg;
      inxt[icomm]++;
    }
  }

  PMMG_DEL_MEM(parmesh,inxt,int,"inxt");
  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of the file to load the mesh from.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Load parallel communicators in Medit format (only one group per process is
 * allowed).
 *
 */
int PMMG_loadCommunicators( PMMG_pParMesh parmesh,const char *filename ) {
  MMG5_pMesh  mesh;
  int         meshver;
  int         API_mode,icomm,ier;
  int         ncomm,*nitem_comm,*color;
  int         **idx_loc,**idx_glo;
  FILE        *inm;
  int         bin;
  long        pos;
  int         iswp;
  int         binch,bpos;
  char        chaine[MMG5_FILESTR_LGTH],strskip[MMG5_FILESTR_LGTH];

  assert( parmesh->ngrp == 1 );
  mesh = parmesh->listgrp[0].mesh;

  /** Open mesh file */
  ier = MMG3D_openMesh(mesh->info.imprim,filename,&inm,&bin,"rb","rb");

  /** Read communicators */
  pos = 0;
  ncomm = 0;
  iswp = 0;
  API_mode = PMMG_UNSET;

  rewind(inm);
  if (!bin) {
    strcpy(chaine,"D");
    while(fscanf(inm,"%127s",&chaine[0])!=EOF && strncmp(chaine,"End",strlen("End")) ) {
      if ( chaine[0] == '#' ) {
        fgets(strskip,MMG5_FILESTR_LGTH,inm);
        continue;
      }

      if(!strncmp(chaine,"ParallelTriangleCommunicators",strlen("ParallelTriangleCommunicators"))) {
        MMG_FSCANF(inm,"%d",&ncomm);
        pos = ftell(inm);
        API_mode = PMMG_APIDISTRIB_faces;
        break;
      } else if(!strncmp(chaine,"ParallelVertexCommunicators",strlen("ParallelVertexCommunicators"))) {
        MMG_FSCANF(inm,"%d",&ncomm);
        pos = ftell(inm);
        API_mode = PMMG_APIDISTRIB_nodes;
        break;
      }
    }
  } else { //binary file
    MMG_FREAD(&meshver,MMG5_SW,1,inm);
    iswp=0;
    if(meshver==16777216)
      iswp=1;
    else if(meshver!=1) {
      fprintf(stderr,"BAD FILE ENCODING\n");
    }

    int endcount = 0;
    while(fread(&binch,MMG5_SW,1,inm)!=0 && endcount != 2 ) {
      if(iswp) binch=MMG5_swapbin(binch);
      if(binch==54) break;
      if(!ncomm && binch==70) { // ParallelTriangleCommunicators
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&ncomm,MMG5_SW,1,inm);
        if(iswp) ncomm=MMG5_swapbin(ncomm);
        pos = ftell(inm);
        API_mode = PMMG_APIDISTRIB_faces;
        break; // if parallel triangles are found, ignore parallel nodes
      } else if(!ncomm && binch==71) { // ParallelVertexCommunicators
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        MMG_FREAD(&ncomm,MMG5_SW,1,inm);
        if(iswp) ncomm=MMG5_swapbin(ncomm);
        pos = ftell(inm);
        API_mode = PMMG_APIDISTRIB_nodes;
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
        continue;
      } else if ( binch==54 ) {
        /* The end keyword will be present twice */
        ++endcount;
      } else {
        MMG_FREAD(&bpos,MMG5_SW,1,inm); //NulPos
        if(iswp) bpos=MMG5_swapbin(bpos);
        rewind(inm);
        fseek(inm,bpos,SEEK_SET);
      }
    }
  }

  /* Set API mode */
  if( API_mode == PMMG_UNSET ) {
    fprintf(stderr,"### Error: No parallel communicators provided on rank %d!\n",parmesh->myrank);
    return 0;
  } else if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_APImode, API_mode ) ) {
    return 0;
  }

  /* memory allocation */
  PMMG_CALLOC(parmesh,nitem_comm,ncomm,int,"nitem_comm",return 0);
  PMMG_CALLOC(parmesh,color,ncomm,int,"color",return 0);
  PMMG_CALLOC(parmesh,idx_loc,ncomm,int*,"idx_loc pointer",return 0);
  PMMG_CALLOC(parmesh,idx_glo,ncomm,int*,"idx_glo pointer",return 0);

  /* Load the communicator */
  if( !PMMG_loadCommunicator( parmesh,inm,bin,iswp,pos,ncomm,nitem_comm,color,
                              idx_loc,idx_glo ) ) return 0;

  /* Set triangles or nodes interfaces depending on API mode */
  switch( API_mode ) {

    case PMMG_APIDISTRIB_faces :

      /* Set the number of interfaces */
      ier = PMMG_Set_numberOfFaceCommunicators(parmesh, ncomm);

      /* Loop on each interface (proc pair) seen by the current rank) */
      for( icomm = 0; icomm < ncomm; icomm++ ) {

        /* Set nb. of entities on interface and rank of the outward proc */
        ier = PMMG_Set_ithFaceCommunicatorSize(parmesh, icomm,
                                               color[icomm],
                                               nitem_comm[icomm]);

        /* Set local and global index for each entity on the interface */
        ier = PMMG_Set_ithFaceCommunicator_faces(parmesh, icomm,
                                                 idx_loc[icomm],
                                                 idx_glo[icomm], 1 );
      }
      break;

    case PMMG_APIDISTRIB_nodes :

      /* Set the number of interfaces */
      ier = PMMG_Set_numberOfNodeCommunicators(parmesh, ncomm);

      /* Loop on each interface (proc pair) seen by the current rank) */
      for( icomm = 0; icomm < ncomm; icomm++ ) {

        /* Set nb. of entities on interface and rank of the outward proc */
        ier = PMMG_Set_ithNodeCommunicatorSize(parmesh, icomm,
                                               color[icomm],
                                               nitem_comm[icomm]);

        /* Set local and global index for each entity on the interface */
        ier = PMMG_Set_ithNodeCommunicator_nodes(parmesh, icomm,
                                                 idx_loc[icomm],
                                                 idx_glo[icomm], 1 );
      }
      break;
  }

  /* Release memory and return */
  PMMG_DEL_MEM(parmesh,nitem_comm,int,"nitem_comm");
  PMMG_DEL_MEM(parmesh,color,int,"color");
  for( icomm = 0; icomm < ncomm; icomm++ ) {
    PMMG_DEL_MEM(parmesh,idx_loc[icomm],int,"idx_loc");
    PMMG_DEL_MEM(parmesh,idx_glo[icomm],int,"idx_glo");
  }
  PMMG_DEL_MEM(parmesh,idx_loc,int*,"idx_loc pointer");
  PMMG_DEL_MEM(parmesh,idx_glo,int*,"idx_glo pointer");

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param endame string to allocate to store the final filename
 * \param initname initial file name in which we want to insert rank index
 * \param ASCIIext extension to search for ASCII format
 * \param binext extension to search for binary format
 *
 * Allocate the endname string and copy the initname string with the mpir rank
 * index before the file extension.
 *
 */
static inline
void PMMG_insert_rankIndex(PMMG_pParMesh parmesh,char **endname,const char *initname,
                           char *ASCIIext, char *binext) {
  int    lenmax;
  int8_t fmt;
  char   *ptr;

  lenmax = PMMG_count_digits ( parmesh->nprocs );

  /* Check for pointer validity */
  if ( (!endname) || (!initname) ) {
    return;
  }

  MMG5_SAFE_CALLOC(*endname,strlen(initname)+lenmax+7,char,return);

  strcpy(*endname,initname);

  ptr = strstr(*endname,binext);

  fmt = 0; /* noext */
  if ( ptr ) {
    *ptr = '\0';
    fmt = 1; /* binary */
  }
  else {
    ptr = strstr(*endname,ASCIIext);
    if ( ptr ) {
      *ptr = '\0';
      fmt = 2; /* ASCII */
    }
  }
  sprintf(*endname, "%s.%d", *endname, parmesh->myrank );
  if ( fmt==1 ) {
    strcat ( *endname, binext );
  }
  else if ( fmt==2 ) {
    strcat ( *endname, ASCIIext );
  }

  return;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of the file to load the mesh from.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Load a distributed mesh with parallel communicators in Medit format (only one
 * group per process is allowed). The rank index is inserted in the input file
 * name.
 *
 */
int PMMG_loadMesh_distributed(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh  mesh;
  int         ier;
  char*       data = NULL;

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }

  mesh = parmesh->listgrp[0].mesh;

  /* Add rank index to mesh name */
  if ( filename ) {
    PMMG_insert_rankIndex(parmesh,&data,filename,".mesh", ".meshb");
  }
  else if ( parmesh->meshin ) {
    PMMG_insert_rankIndex(parmesh,&data,parmesh->meshin,".mesh", ".meshb");
  }
  else if ( mesh->namein ) {
    PMMG_insert_rankIndex(parmesh,&data,mesh->namein,".mesh", ".meshb");
  }

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  ier = MMG3D_loadMesh(mesh,data);

  /* Restore the mmg verbosity to its initial value */
  mesh->info.imprim = parmesh->info.mmg_imprim;

  if ( ier < 1 ) {
    MMG5_SAFE_FREE(data);
    return ier;
  }

  /* Load parallel communicators */
  ier = PMMG_loadCommunicators( parmesh,data );

  MMG5_SAFE_FREE(data);

  if ( 1 != ier ) return 0;

  return 1;
}

int PMMG_loadMesh_centralized(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;
  int        ier;
  const char *data;

  if ( parmesh->myrank!=parmesh->info.root ) {
    return 1;
  }

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }
  mesh = parmesh->listgrp[0].mesh;

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  if ( filename ) {
    data = filename;
  }
  else if ( parmesh->meshin ) {
    data = parmesh->meshin;
  }
  else if ( mesh->namein ) {
    data = mesh->namein;
  }
  else {
    data = NULL;
  }
  ier = MMG3D_loadMesh(mesh,data);

  /* Restore the mmg verbosity to its initial value */
  mesh->info.imprim = parmesh->info.mmg_imprim;

  return ier;
}

int PMMG_loadMet_centralized(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  int        ier;
  const char *data;

  if ( parmesh->myrank!=parmesh->info.root ) {
    return 1;
  }

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }
  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  if ( filename ) {
    data = filename;
  }
  else if ( parmesh->metin ) {
    data = parmesh->metin;
  }
  else if ( met->namein ) {
    data = met->namein;
  }
  else {
    data = NULL;
  }
  ier = MMG3D_loadSol(mesh,met,data);

  /* Restore the mmg verbosity to its initial value */
  mesh->info.imprim = parmesh->info.mmg_imprim;

  return ier;
}

int PMMG_loadMet_distributed(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  int        ier;
  char       *data = NULL;

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }

  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  /* Add rank index to mesh name */
  if ( filename ) {
    PMMG_insert_rankIndex(parmesh,&data,filename,".sol", ".sol");
  }
  else if ( parmesh->metin ) {
    PMMG_insert_rankIndex(parmesh,&data,parmesh->metin,".sol", ".sol");
  }
  else if ( met->namein ) {
    PMMG_insert_rankIndex(parmesh,&data,met->namein,".sol", ".sol");
  }
  else if ( parmesh->meshin ) {
    PMMG_insert_rankIndex(parmesh,&data,parmesh->meshin,".mesh", ".meshb");
  }
  else if ( mesh->namein ) {
    PMMG_insert_rankIndex(parmesh,&data,mesh->namein,".mesh", ".meshb");
  }

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  ier = MMG3D_loadSol(mesh,met,data);

  /* Restore the mmg verbosity to its initial value */
  mesh->info.imprim = parmesh->info.mmg_imprim;

  MMG5_SAFE_FREE(data);

  return ier;
}

int PMMG_loadLs_centralized(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;
  MMG5_pSol  ls;
  int        ier;
  const char *data;

  if ( parmesh->myrank!=parmesh->info.root ) {
    return 1;
  }

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }
  mesh = parmesh->listgrp[0].mesh;
  ls   = parmesh->listgrp[0].ls;

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  if ( filename ) {
    data = filename;
  }
  else if ( parmesh->lsin ) {
    data = parmesh->lsin;
  }
  else if ( ls->namein ) {
    data = ls->namein;
  }
  else {
    data = NULL;
  }
  ier = MMG3D_loadSol(mesh,ls,data);

  /* Restore the mmg verbosity to its initial value */
  mesh->info.imprim = parmesh->info.mmg_imprim;

  return ier;
}

int PMMG_loadDisp_centralized(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;
  MMG5_pSol  disp;
  int        ier;
  const char *data;

  if ( parmesh->myrank!=parmesh->info.root ) {
    return 1;
  }

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }
  mesh = parmesh->listgrp[0].mesh;
  disp = parmesh->listgrp[0].disp;

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  if ( filename ) {
    data = filename;
  }
  else if ( parmesh->dispin ) {
    data = parmesh->dispin;
  }
  else if ( disp->namein ) {
    data = disp->namein;
  }
  else {
    data = NULL;
  }
  ier = MMG3D_loadSol(mesh,disp,data);

  /* Restore the mmg verbosity to its initial value */
  mesh->info.imprim = parmesh->info.mmg_imprim;

  return ier;
}

int PMMG_loadSol_centralized(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;
  MMG5_pSol  sol;
  int        ier;
  const char *namein;

  if ( parmesh->myrank!=parmesh->info.root ) {
    return 1;
  }

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }
  mesh = parmesh->listgrp[0].mesh;

  /* For each mode: pointer over the solution structure to load */
  if ( mesh->info.lag >= 0 ) {
    sol = parmesh->listgrp[0].disp;
  }
  else if ( mesh->info.iso ) {
    sol = parmesh->listgrp[0].ls;
  }
  else {
    sol = parmesh->listgrp[0].met;
  }

  if ( !filename ) {
    namein = sol->namein;
  }
  else {
    namein = filename;
  }
  assert ( namein );

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  ier = MMG3D_loadSol(mesh,sol,namein);

  /* Restore the mmg verbosity to its initial value */
  mesh->info.imprim = parmesh->info.mmg_imprim;

  return ier;
}

int PMMG_loadAllSols_centralized(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;
  MMG5_pSol  *sol;
  int        ier;
  const char *data;

  if ( parmesh->myrank!=parmesh->info.root ) {
    return 1;
  }

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }
  mesh = parmesh->listgrp[0].mesh;
  sol  = &parmesh->listgrp[0].field;

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  if ( filename ) {
    data = filename;
  }
  else if ( parmesh->fieldin ) {
    data = parmesh->fieldin;
  }
  else {
    data = NULL;
  }
  ier = MMG3D_loadAllSols(mesh,sol,data);

  /* Restore the mmg verbosity to its initial value */
  mesh->info.imprim = parmesh->info.mmg_imprim;

  return ier;

}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of the file to load the mesh from.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Save a distributed mesh with parallel communicators in Medit format (only one
 * group per process is allowed).
 *
 */
int PMMG_saveMesh_distributed(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh  mesh;
  int         ier;
  char       *data = NULL;

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }

  mesh = parmesh->listgrp[0].mesh;

  /* Add rank index to mesh name */
  if ( filename ) {
    PMMG_insert_rankIndex(parmesh,&data,filename,".mesh", ".meshb");
  }
  else if ( parmesh->meshout ) {
    PMMG_insert_rankIndex(parmesh,&data,parmesh->meshout,".mesh", ".meshb");
  }
  else if ( mesh->nameout ) {
    PMMG_insert_rankIndex(parmesh,&data,mesh->nameout,".mesh", ".meshb");
  }

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  ier = MMG3D_saveMesh(mesh,data);

  /* Restore the mmg verbosity to its initial value */
  mesh->info.imprim = parmesh->info.mmg_imprim;

  if ( ier < 1 ) {
    MMG5_SAFE_FREE(data);
    return ier;
  }

  /* Load parallel communicators */
  ier = PMMG_printCommunicator ( parmesh,data );

  MMG5_SAFE_FREE ( data );

  if ( 1 != ier ) return 0;

  return 1;
}


int PMMG_saveMesh_centralized(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;
  int        ier;

  if ( parmesh->myrank!=parmesh->info.root ) {
    return 1;
  }

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }
  mesh = parmesh->listgrp[0].mesh;

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  if ( filename && *filename ) {
    ier = MMG3D_saveMesh(mesh,filename);
  }
  else {
    ier = MMG3D_saveMesh(mesh,parmesh->meshout);
  }

  /* Restore the mmg verbosity to its initial value */
  mesh->info.imprim = parmesh->info.mmg_imprim;

  return ier;
}

int PMMG_saveMet_centralized(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  int        ier;

  if ( parmesh->myrank!=parmesh->info.root ) {
    return 1;
  }

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }
  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  if ( filename && *filename ) {
    ier =  MMG3D_saveSol(mesh,met,filename);
  }
  else {
    ier =  MMG3D_saveSol(mesh,met,parmesh->metout);
  }

  /* Restore the mmg verbosity to its initial value */
  mesh->info.imprim = parmesh->info.mmg_imprim;

  return ier;
}

int PMMG_saveMet_distributed(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  int        ier;
  char       *data = NULL;

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }

  mesh = parmesh->listgrp[0].mesh;
  met  = parmesh->listgrp[0].met;

  /* Add rank index to mesh name */
  if ( filename ) {
    PMMG_insert_rankIndex(parmesh,&data,filename,".sol", ".sol");
  }
  else if ( parmesh->metout ) {
    PMMG_insert_rankIndex(parmesh,&data,parmesh->metout,".sol", ".sol");
  }
  else if ( met->nameout ) {
    PMMG_insert_rankIndex(parmesh,&data,met->nameout,".sol", ".sol");
  }
  else if ( parmesh->meshout ) {
    PMMG_insert_rankIndex(parmesh,&data,parmesh->meshout,".mesh", ".meshb");
  }
  else if ( mesh->nameout ) {
    PMMG_insert_rankIndex(parmesh,&data,mesh->nameout,".mesh", ".meshb");
  }

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  ier =  MMG3D_saveSol(mesh,met,data);

  /* Restore the mmg verbosity to its initial value */
  mesh->info.imprim = parmesh->info.mmg_imprim;

  MMG5_SAFE_FREE ( data );

  return ier;
}

int PMMG_saveAllSols_centralized(PMMG_pParMesh parmesh,const char *filename) {
  MMG5_pMesh mesh;
  MMG5_pSol  sol;
  int        ier;

  if ( parmesh->myrank!=parmesh->info.root ) {
    return 1;
  }

  if ( parmesh->ngrp != 1 ) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }
  mesh = parmesh->listgrp[0].mesh;
  sol  = parmesh->listgrp[0].field;

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  if ( filename && *filename ) {
    ier = MMG3D_saveAllSols(mesh,&sol,filename);
  }
  else {
    ier = MMG3D_saveAllSols(mesh,&sol,parmesh->fieldout);
  }

  /* Restore the mmg verbosity to its initial value */
  mesh->info.imprim = parmesh->info.mmg_imprim;

  return ier;
}

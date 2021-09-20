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
  int         iswp,k;
  int         binch,bpos;
  char        chaine[MMG5_FILESTR_LGTH],strskip[MMG5_FILESTR_LGTH];

  assert( parmesh->ngrp == 1 );
  mesh = parmesh->listgrp[0].mesh;

  /* A non-// tria may be marked as // in Medit serial I/O (if its 3 edges are
   * //): as we can infer // triangles from communicators, reset useless (and
   * maybe erroneous) tags */
  for ( k=1; k<=mesh->nt; ++k ) {
    if ( (mesh->tria[k].tag[0] & MG_PARBDY) &&
         (mesh->tria[k].tag[1] & MG_PARBDY) &&
         (mesh->tria[k].tag[2] & MG_PARBDY) ) {
      mesh->tria[k].tag[0] &= ~MG_PARBDY;
      mesh->tria[k].tag[1] &= ~MG_PARBDY;
      mesh->tria[k].tag[2] &= ~MG_PARBDY;
    }
  }

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

  /* Close the file */
  if( filename ) {
    fclose(inm);
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

int PMMG_Set_defaultIOEntities_hdf5(int *io_entities) {
  /* Default: save/load everything */
  for (int i = 0 ; i < PMMG_NTYPENTITIES ; i++) io_entities[i] = 1;

  return 1;
}

int PMMG_Set_requiredEntitiesIO_hdf5(int *io_entities, int val) {
  io_entities[PMMG_IO_Req] = val;
  io_entities[PMMG_IO_EdReq] = val;
  io_entities[PMMG_IO_TriaReq] = val;
  io_entities[PMMG_IO_QuadReq] = val;
  io_entities[PMMG_IO_TetReq] = val;

  return 1;
}

int PMMG_Set_parallelEntitiesIO_hdf5(int *io_entities, int val) {
  io_entities[PMMG_IO_Par] = val;
  io_entities[PMMG_IO_EdPar] = val;
  io_entities[PMMG_IO_TriaPar] = val;
  io_entities[PMMG_IO_QuadPar] = val;
  io_entities[PMMG_IO_TetPar] = val;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param nentities array of size PMMG_NTYPENTITIES * nprocs that will contain the number of entities of every proc.
 * \param nentitiesl array of size PMMG_NTYP_ENTITIES that will contain the local number of entities.
 * \param nentitiesg array of size PMMG_NTYP_ENTITIES that will contain the global number of entities.
 *
 * \return 1, there is no reason for this function to fail.
 *
 * Count the local and global number of entities in the parallel mesh (only one
 * group per process is allowed).
 *
 */
static int PMMG_countEntities(PMMG_pParMesh parmesh, hsize_t *nentities, hsize_t *nentitiesl, hsize_t* nentitiesg, int *io_entities) {
  /* MMG variables */
  PMMG_pGrp   grp;
  MMG5_pMesh  mesh;
  MMG5_pPoint ppt;
  MMG5_pEdge  pa;
  MMG5_pTria  pt;
  MMG5_pQuad  pq;
  MMG5_pTetra pe;
  MMG5_pPrism pp;
  /* Local number of entities */
  hsize_t ne, np, nt, na, nquad, nprism;       /* Tetra, points, triangles, edges, quads, prisms */
  hsize_t nc, npreq, nppar;                      /* Corners, required and parallel vertices */
  hsize_t nr, nedreq, nedpar;                  /* Ridges, required and parallel edges */
  hsize_t ntreq, ntpar;                        /* Required and parallel triangles */
  hsize_t nqreq, nqpar;                        /* Required and parallel quads */
  hsize_t nereq, nepar;                        /* Required and parallel tetra */
  hsize_t nnor, ntan;                          /* Normals and Tangents */
  /* MPI variables */
  int rank, root, nprocs;

  /* Only one group is allowed */
  assert( parmesh->ngrp == 1 );

  /* Set MPI variables */
  nprocs = parmesh->nprocs;
  rank = parmesh->myrank;
  root = parmesh->info.root;

  /* Set mesh size to 0 */
  np = na = nt = nquad = ne = nprism = 0;
  nc = npreq = nppar = 0;
  nr = nedreq = nedpar = 0;
  ntreq = ntpar = 0;
  nqreq = nqpar = 0;
  nereq = nepar = 0;
  nnor = ntan = 0;

  /* Set ParMmg variables */
  grp = &parmesh->listgrp[0];
  mesh = grp->mesh;
  ppt = NULL;
  pa = NULL;
  pt = NULL;
  pq = NULL;
  pe = NULL;
  pp = NULL;

  /* Check arguments */
  assert( nentities && "\n  ## Error: %s: nentities array not allocated.\n" );
  assert( nentitiesl && "\n  ## Error: %s: nentitiesl array not allocated.\n" );
  assert( nentitiesg && "\n  ## Error: %s: nentitiesg array not allocated.\n" );

  /* Count local entities */

  /* Vertices, normals and tangents */
  for (int k = 1 ; k <= mesh->np ; k++) {
    ppt = &mesh->point[k];
    if (MG_VOK(ppt)) {
      ppt->tmp = ++np;
      ppt->flag = 0;
      if (ppt->tag & MG_CRN) nc++;
      if (ppt->tag & MG_REQ) npreq++;
      if (ppt->tag & MG_PARBDY) nppar++;
      if (MG_SIN(ppt->tag)) continue;
      if ((ppt->tag & MG_BDY) && (!(ppt->tag & MG_GEO) || (ppt->tag & MG_NOM))) nnor++;
      if (MG_EDG(ppt->tag) || (ppt->tag & MG_NOM)) ntan++;
    }
  }

  /* Edges */
  if (mesh->na) {
    for (int k = 1 ; k <= mesh->na ; k++) {
      pa = &mesh->edge[k];
      na++;
      if (pa->tag & MG_GEO) nr++;
      if (pa->tag & MG_REQ) nedreq++;
      if (pa->tag & MG_PARBDY) nedpar++;
    }
  }

  /* Triangles */
  if (mesh->nt) {
    for (int k = 1 ; k <= mesh->nt ; k++) {
      pt = &mesh->tria[k];
      nt++;
      if (pt->tag[0] & MG_REQ && pt->tag[1] & MG_REQ && pt->tag[2] & MG_REQ) ntreq++;
      if (pt->tag[0] & MG_PARBDY && pt->tag[1] & MG_PARBDY && pt->tag[2] & MG_PARBDY) ntpar++;
    }
  }

  /* Quadrilaterals */
  if (mesh->nquad) {
    for (int k = 1 ; k <= mesh->nquad ; k++) {
      pq = &mesh->quadra[k];
      nquad++;
      if (pq->tag[0] & MG_REQ && pq->tag[1] & MG_REQ &&
          pq->tag[2] & MG_REQ && pq->tag[3] & MG_REQ) {
        nqreq++;
      }
      if (pq->tag[0] & MG_PARBDY && pq->tag[1] & MG_PARBDY &&
          pq->tag[2] & MG_PARBDY && pq->tag[3] & MG_PARBDY) {
        nqpar++;
      }
    }
  }

  /* Tetrahedra */
  if (mesh->ne) {
    for (int k = 1 ; k <= mesh->ne ; k++) {
      pe = &mesh->tetra[k];
      if (!MG_EOK(pe)) {
        continue;
      }
      ne++;
      if (pe->tag & MG_REQ) nereq++;
      if (pe->tag & MG_PARBDY) nepar++;
    }
  }

  /* Prisms */
  if (mesh->nprism) {
    for (int k = 1 ; k <= mesh->nprism ; k++) {
      pp = &mesh->prism[k];
      if (!MG_EOK(pp)){
        continue;
      }
      nprism++;
    }
  }

  /* Gather the number of entities in the nentities array */
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_Vertex]  = np;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_Edge]    = na;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_Tria]    = nt;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_Quad]    = nquad;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_Tetra]   = ne;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_Prism]   = nprism;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_Corner]  = nc;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_Req]     = npreq;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_Par]     = nppar;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_Ridge]   = nr;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_EdReq]   = nedreq;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_EdPar]   = nedpar;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_TriaReq] = ntreq;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_TriaPar] = ntpar;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_QuadReq] = nqreq;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_QuadPar] = nqpar;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_TetReq]  = nereq;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_TetPar]  = nepar;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_Normal]  = nnor;
  nentities[PMMG_NTYPENTITIES * rank + PMMG_IO_Tangent] = ntan;

  /* Set the number of all entities that are not saved to zero */
  for (int typent = 0 ; typent < PMMG_NTYPENTITIES ; typent++) {
    if (!io_entities[typent]) {
      nentities[PMMG_NTYPENTITIES * rank + typent] = 0;
    }
  }

  /* Count local entities */
  for (int typent = 0 ; typent < PMMG_NTYPENTITIES ; typent++) {
    nentitiesl[typent] = nentities[PMMG_NTYPENTITIES * rank + typent];
  }

  MPI_Allgather(&nentities[PMMG_NTYPENTITIES * rank], PMMG_NTYPENTITIES, MPI_UNSIGNED_LONG_LONG,
                nentities                           , PMMG_NTYPENTITIES, MPI_UNSIGNED_LONG_LONG, parmesh->comm);

  /* Count global entities */
  for (int k = 0 ; k < nprocs ; k++) {
    for (int typent = 0 ; typent < PMMG_NTYPENTITIES ; typent++) {
      nentitiesg[typent] += nentities[PMMG_NTYPENTITIES * k + typent];
    }
  }

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param nentities array of size PMMG_NTYPENTITIES * nprocs that contains the number of entities of every proc.
 * \param offset array of size 2 * PMMG_NTYPENTITIES that will contain the offsets for each type of entity.
 *
 * \return 1
 *
 * Compute the offset for parallel writing/reading in an HDF5 file.
 *
 */
static inline int PMMG_computeHDFoffset(PMMG_pParMesh parmesh, hsize_t *nentities, hsize_t *offset) {
  for (int k = 0 ; k < parmesh->myrank ; k++) {
    for (int typent = 0 ; typent < PMMG_NTYPENTITIES ; typent++) {
      offset[2 * typent] += nentities[PMMG_NTYPENTITIES * k + typent];
    }
  }
  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param file_id identifier of the HDF5 file in which to write the mesh
 *
 * \return 1, there is no reason for this function to fail.
 *
 * Save the version and the dimension of \a parmesh aswell as the number of
 * partitions and the API mode into the opened HDF5 file \a file_id.
 *
 */
static int PMMG_saveHeader_hdf5(PMMG_pParMesh parmesh, hid_t file_id) {
  MMG5_pMesh mesh;
  hid_t      dspace_id;
  hid_t      attr_id;
  int        rank, root;

  mesh = parmesh->listgrp[0].mesh;
  rank = parmesh->myrank;
  root = parmesh->info.root;

  dspace_id = H5Screate(H5S_SCALAR);

  attr_id = H5Acreate(file_id, "MeshVersionFormatted", H5T_NATIVE_INT, dspace_id, H5P_DEFAULT, H5P_DEFAULT);
  if (rank == root)
    H5Awrite(attr_id, H5T_NATIVE_INT, &mesh->ver);
  H5Aclose(attr_id);

  attr_id = H5Acreate(file_id, "Dimension", H5T_NATIVE_INT, dspace_id, H5P_DEFAULT, H5P_DEFAULT);
  if (rank == root)
    H5Awrite(attr_id, H5T_NATIVE_INT, &mesh->dim);
  H5Aclose(attr_id);

  attr_id = H5Acreate(file_id, "NumberOfPartitions", H5T_NATIVE_INT, dspace_id, H5P_DEFAULT, H5P_DEFAULT);
  if (rank == root)
    H5Awrite(attr_id, H5T_NATIVE_INT, &parmesh->nprocs);
  H5Aclose(attr_id);

  attr_id = H5Acreate(file_id, "API_mode", H5T_NATIVE_INT, dspace_id, H5P_DEFAULT, H5P_DEFAULT);
  if (rank == root)
    H5Awrite(attr_id, H5T_NATIVE_INT, &parmesh->info.API_mode);
  H5Aclose(attr_id);

  H5Sclose(dspace_id);

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param grp_entities_id identifier of the HDF5 group in which to write the mesh entities.
 * \param dcpl_id identifier of the dataset creation property list (no fill value).
 * \param dxpl_id identifier of the dataset transfer property list (MPI-IO).
 * \param nentitiesl array of size PMMG_NTYP_ENTITIES containing the local number of entities.
 * \param nentitiesg array of size PMMG_NTYP_ENTITIES containing the global number of entities.
 * \param offset array of size PMMG_NTYP_ENTITIES containing the offset for parallel writing.
 * \param save_entities array of size PMMG_NTYP_ENTITIES to tell which entities are to be saved.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Save the mesh entities in the \a grp_entities_id group of an HDF5 file (only
 * one group per process is allowed).
 *
 */
static int PMMG_saveMeshEntities_hdf5(PMMG_pParMesh parmesh, hid_t grp_entities_id, hid_t dcpl_id, hid_t dxpl_id,
                                      hsize_t *nentitiesl, hsize_t *nentitiesg, hsize_t *offset, int *save_entities) {
  /* MMG variables */
  PMMG_pGrp    grp;
  MMG5_pMesh   mesh;
  MMG5_pPoint  ppt;
  MMG5_pEdge   pa;
  MMG5_pTria   pt;
  MMG5_pQuad   pq;
  MMG5_pTetra  pe;
  MMG5_pPrism  pp;
  MMG5_pxPoint pxp;

  /* Local mesh size */
  hsize_t ne, np, nt, na, nquad, nprism;       /* Tetra, points, triangles, edges, quads, prisms */
  hsize_t nc, npreq, nppar;                    /* Corners, required and parallel vertices */
  hsize_t nr, nedreq, nedpar;                  /* Ridges, required and parallel edges */
  hsize_t ntreq, ntpar;                        /* Required and parallel triangles */
  hsize_t nqreq, nqpar;                        /* Required and parallel quads */
  hsize_t nereq, nepar;                        /* Required and parallel tetra */
  hsize_t nnor, ntan;                          /* Normals and Tangents */
  /* Global mesh size */
  hsize_t neg, npg, ntg, nag, nquadg, nprismg; /* Tetra, points, triangles, edges, quads, prisms */
  hsize_t ncg, npreqg, npparg;                 /* Corners, required and parallel vertices */
  hsize_t nrg, nedreqg, nedparg;               /* Ridges, required and parallel edges */
  hsize_t ntreqg, ntparg;                      /* Required and parallel triangles */
  hsize_t nqreqg, nqparg;                      /* Required and parallel quads */
  hsize_t nereqg, neparg;                      /* Required and parallel tetra */
  hsize_t nnorg, ntang;                        /* Normals and Tangents */

  /* Mesh buffer arrays */
  /* 6 buffers is the minimum amount for what we have to do */
  double *ppoint;      /* Point coordinates */
  int    *pent;        /* Other entities : edges, trias, quads, tetra, prisms. */
  int    *pcr;         /* Corners and ridges */
  int    *preq, *ppar; /* Required and parallel entities */
  int    *pref;        /* References */

  /* Normals and tangents */
  /* We could reuse the previous buffers, but the names would be confusing */
  int    *pnorat, *ptanat; /* Normals and Tangents at vertices */
  double *pnor, *ptan;     /* Normals and Tangents */

  /* Counters for the corners/ridges, the required and parallel entities, the normals and the tangents */
  int crcount, reqcount, parcount, ncount, tcount;

  /* Flag to remember if the save_entities was NULL or not */
  int nullf = 0;

  /* MPI variables */
  hsize_t rank, root;
  hsize_t nprocs;

  /* HDF5 variables */
  hid_t dspace_mem_id, dspace_file_id;
  hid_t dset_id;

  /*------------------------- INIT -------------------------*/

  assert ( parmesh->ngrp == 1 );

  /* Set all buffers to NULL */
  ppoint = NULL;
  pent = NULL;
  pcr = NULL;
  preq = NULL; ppar = NULL;
  pref = NULL;
  pnor = NULL; ptan = NULL;
  pnorat = NULL; ptanat = NULL;

  /* Set MPI variables */
  nprocs = parmesh->nprocs;
  rank = parmesh->myrank;
  root = parmesh->info.root;

  /* Set ParMmg variables */
  grp = &parmesh->listgrp[0];
  mesh = grp->mesh;
  ppt = NULL;
  pa = NULL;
  pt = NULL;
  pq = NULL;
  pe = NULL;
  pp = NULL;

  /* Check the save_entities argument */
  if (save_entities == NULL) {
    nullf = 1;
    PMMG_CALLOC(parmesh, save_entities, PMMG_NTYPENTITIES, int, "save_entities", return 0);
    PMMG_Set_defaultIOEntities_hdf5(save_entities);
  }

  /* Get the number of entities */
  np     = nentitiesl[PMMG_IO_Vertex];
  na     = nentitiesl[PMMG_IO_Edge];
  nt     = nentitiesl[PMMG_IO_Tria];
  nquad  = nentitiesl[PMMG_IO_Quad];
  ne     = nentitiesl[PMMG_IO_Tetra];
  nprism = nentitiesl[PMMG_IO_Prism];
  nc     = nentitiesl[PMMG_IO_Corner];
  npreq  = nentitiesl[PMMG_IO_Req];
  nppar  = nentitiesl[PMMG_IO_Par];
  nr     = nentitiesl[PMMG_IO_Ridge];
  nedreq = nentitiesl[PMMG_IO_EdReq];
  nedpar = nentitiesl[PMMG_IO_EdPar];
  ntreq  = nentitiesl[PMMG_IO_TriaReq];
  ntpar  = nentitiesl[PMMG_IO_TriaPar];
  nqreq  = nentitiesl[PMMG_IO_QuadReq];
  nqpar  = nentitiesl[PMMG_IO_QuadPar];
  nereq  = nentitiesl[PMMG_IO_TetReq];
  nepar  = nentitiesl[PMMG_IO_TetPar];
  nnor   = nentitiesl[PMMG_IO_Normal];
  ntan   = nentitiesl[PMMG_IO_Tangent];

  npg     = nentitiesg[PMMG_IO_Vertex];
  nag     = nentitiesg[PMMG_IO_Edge];
  ntg     = nentitiesg[PMMG_IO_Tria];
  nquadg  = nentitiesg[PMMG_IO_Quad];
  neg     = nentitiesg[PMMG_IO_Tetra];
  nprismg = nentitiesg[PMMG_IO_Prism];
  ncg     = nentitiesg[PMMG_IO_Corner];
  npreqg  = nentitiesg[PMMG_IO_Req];
  npparg  = nentitiesg[PMMG_IO_Par];
  nrg     = nentitiesg[PMMG_IO_Ridge];
  nedreqg = nentitiesg[PMMG_IO_EdReq];
  nedparg = nentitiesg[PMMG_IO_EdPar];
  ntreqg  = nentitiesg[PMMG_IO_TriaReq];
  ntparg  = nentitiesg[PMMG_IO_TriaPar];
  nqreqg  = nentitiesg[PMMG_IO_QuadReq];
  nqparg  = nentitiesg[PMMG_IO_QuadPar];
  nereqg  = nentitiesg[PMMG_IO_TetReq];
  neparg  = nentitiesg[PMMG_IO_TetPar];
  nnorg   = nentitiesg[PMMG_IO_Normal];
  ntang   = nentitiesg[PMMG_IO_Tangent];

  /* Arrays for bidimensional dataspaces */
  hsize_t hnp[2]      = {np, 3};
  hsize_t hna[2]      = {na, 2};
  hsize_t hnt[2]      = {nt, 3};
  hsize_t hnquad[2]   = {nquad, 4};
  hsize_t hne[2]      = {ne, 4};
  hsize_t hnprism[2]  = {nprism, 2};
  hsize_t hnnor[2]    = {nnor, 3};
  hsize_t hntan[2]    = {ntan, 3};
  hsize_t hnpg[2]     = {npg, 3};
  hsize_t hnag[2]     = {nag, 2};
  hsize_t hntg[2]     = {ntg, 3};
  hsize_t hnquadg[2]  = {nquadg, 4};
  hsize_t hneg[2]     = {neg, 4};
  hsize_t hnprismg[2] = {nprismg, 2};
  hsize_t hnnorg[2]   = {nnorg, 3};
  hsize_t hntang[2]   = {ntang, 3};

  /* Vertices, Normals and Tangents */
  if (save_entities[PMMG_IO_Vertex] && npg) {

    PMMG_MALLOC(parmesh, ppoint, 3 * np, double, "ppoint", goto free_and_return );
    PMMG_MALLOC(parmesh, pref, np, int, "pref", goto free_and_return );

    if (save_entities[PMMG_IO_Corner]) PMMG_MALLOC(parmesh, pcr, nc, int, "pcr",
                                                   goto free_and_return );
    if (save_entities[PMMG_IO_Req])    PMMG_MALLOC(parmesh, preq, npreq, int, "preq",
                                                   goto free_and_return );
    if (save_entities[PMMG_IO_Par])    PMMG_MALLOC(parmesh, ppar, nppar, int, "ppar",
                                                   goto free_and_return );

    if (save_entities[PMMG_IO_Normal]) {
      PMMG_MALLOC(parmesh, pnor, 3 * nnor, double, "pnor", goto free_and_return );
      PMMG_MALLOC(parmesh, pnorat, nnor, int, "pnorat", goto free_and_return );
    }

    if (save_entities[PMMG_IO_Tangent]) {
      PMMG_MALLOC(parmesh, ptan, 3 * ntan, double, "ptan", goto free_and_return );
      PMMG_MALLOC(parmesh, ptanat, ntan, int, "ptanat", goto free_and_return );
    }

    crcount = reqcount = parcount = ncount = tcount = 0;

    for (int i = 0 ; i < mesh->np ; i++) {
      ppt = &mesh->point[i + 1];
      if (!MG_VOK(ppt))  continue;
      for (int j = 0 ; j < 3 ; j++) {
        ppoint[3 * (ppt->tmp - 1) + j] = ppt->c[j];
      }
      if (save_entities[PMMG_IO_Corner] && (ppt->tag & MG_CRN)) {
        pcr[crcount++] = ppt->tmp + offset[2 * PMMG_IO_Vertex] - 1;
      }
      if (save_entities[PMMG_IO_Req] && (ppt->tag & MG_REQ)) {
        preq[reqcount++] = ppt->tmp + offset[2 * PMMG_IO_Vertex] - 1;
      }
      if (save_entities[PMMG_IO_Par] && (ppt->tag & MG_PARBDY)) {
        ppar[parcount++] = ppt->tmp + offset[2 * PMMG_IO_Vertex] - 1;
      }
      if (!MG_SIN(ppt->tag)) {
        /* Normals */
        if (save_entities[PMMG_IO_Normal]) {
          if ((ppt->tag & MG_BDY) && (!(ppt->tag & MG_GEO) || (ppt->tag & MG_NOM))) {
            pxp = &mesh->xpoint[ppt->xp];
            for (int j = 0 ; j < 3 ; j++) {
              pnor[3 * ncount + j] = pxp->n1[j];
            }
            pnorat[ncount++] = ppt->tmp + offset[2 * PMMG_IO_Vertex] - 1;
          }
        }
        /* Tangents */
        if (save_entities[PMMG_IO_Tangent]) {
          if (MG_EDG(ppt->tag) || (ppt->tag & MG_NOM)) {
            for (int j = 0 ; j < 3 ; j++) {
              ptan[3 * tcount + j] = ppt->n[j];
            }
            ptanat[tcount++] = ppt->tmp + offset[2 * PMMG_IO_Vertex] - 1;
          }
        }
      }
      pref[ppt->tmp - 1] = abs(ppt->ref);
    }

    dspace_mem_id  = H5Screate_simple(2, hnp, NULL);
    dspace_file_id = H5Screate_simple(2, hnpg, NULL);
    H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Vertex], NULL, hnp, NULL);
    dset_id = H5Dcreate(grp_entities_id, "Vertices", H5T_NATIVE_DOUBLE, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dspace_mem_id, dspace_file_id, dxpl_id, ppoint);
    H5Dclose(dset_id);
    H5Sclose(dspace_mem_id);
    H5Sclose(dspace_file_id);
    PMMG_DEL_MEM(parmesh, ppoint, double, "ppoint");

    dspace_mem_id  = H5Screate_simple(1, hnp, NULL);
    dspace_file_id = H5Screate_simple(1, hnpg, NULL);
    H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Vertex], NULL, hnp, NULL);
    dset_id = H5Dcreate(grp_entities_id, "VerticesRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
    H5Dclose(dset_id);
    H5Sclose(dspace_mem_id);
    H5Sclose(dspace_file_id);
    PMMG_DEL_MEM(parmesh, pref, int, "pref");

    if (save_entities[PMMG_IO_Corner] && ncg) {
      dspace_mem_id  = H5Screate_simple(1, &nc, NULL);
      dspace_file_id = H5Screate_simple(1, &ncg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Corner], NULL, &nc, NULL);
      dset_id = H5Dcreate(grp_entities_id, "Corners", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
      H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pcr);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);
      PMMG_DEL_MEM(parmesh, pcr, int, "pcr");
    }

    if (save_entities[PMMG_IO_Req] && npreqg) {
      dspace_mem_id  = H5Screate_simple(1, &npreq, NULL);
      dspace_file_id = H5Screate_simple(1, &npreqg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Req], NULL, &npreq, NULL);
      dset_id = H5Dcreate(grp_entities_id, "RequiredVertices", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
      H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);
      PMMG_DEL_MEM(parmesh, preq, int, "preq");
    }

    if (save_entities[PMMG_IO_Par] && npparg) {
      dspace_mem_id  = H5Screate_simple(1, &nppar, NULL);
      dspace_file_id = H5Screate_simple(1, &npparg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Par], NULL, &nppar, NULL);
      dset_id = H5Dcreate(grp_entities_id, "ParallelVertices", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
      H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);
      PMMG_DEL_MEM(parmesh, ppar, int, "ppar");
    }

    if (save_entities[PMMG_IO_Normal] && nnorg) {
      dspace_mem_id  = H5Screate_simple(2, hnnor, NULL);
      dspace_file_id = H5Screate_simple(2, hnnorg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Normal], NULL, hnnor, NULL);
      dset_id = H5Dcreate(grp_entities_id, "Normals", H5T_NATIVE_DOUBLE, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dspace_mem_id, dspace_file_id, dxpl_id, pnor);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);
      PMMG_DEL_MEM(parmesh, pnor, double, "pnor");

      dspace_mem_id  = H5Screate_simple(1, hnnor, NULL);
      dspace_file_id = H5Screate_simple(1, hnnorg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Normal], NULL, hnnor, NULL);
      dset_id = H5Dcreate(grp_entities_id, "NormalsAtVertices", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
      H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pnorat);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);
      PMMG_DEL_MEM(parmesh, pnorat, int, "pnorat");
    }

    if (save_entities[PMMG_IO_Tangent] && ntang) {
      dspace_mem_id  = H5Screate_simple(2, hntan, NULL);
      dspace_file_id = H5Screate_simple(2, hntang, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Tangent], NULL, hntan, NULL);
      dset_id = H5Dcreate(grp_entities_id, "Tangents", H5T_NATIVE_DOUBLE, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
      H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dspace_mem_id, dspace_file_id, dxpl_id, ptan);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);
      PMMG_DEL_MEM(parmesh, ptan, double, "ptan");

      dspace_mem_id  = H5Screate_simple(1, hntan, NULL);
      dspace_file_id = H5Screate_simple(1, hntang, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Tangent], NULL, hntan, NULL);
      dset_id = H5Dcreate(grp_entities_id, "TangentsAtVertices", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
      H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ptanat);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);
      PMMG_DEL_MEM(parmesh, ptanat, int, "ptanat");
    }

  }

  /* Edges */
  if (save_entities[PMMG_IO_Edge] && nag) {

    PMMG_MALLOC(parmesh, pent, 2 * na, int, "pent", goto free_and_return);
    PMMG_MALLOC(parmesh, pref, na, int, "pref", goto free_and_return);
    if (save_entities[PMMG_IO_Ridge]) PMMG_MALLOC(parmesh, pcr , nr    , int, "pcr" , goto free_and_return);
    if (save_entities[PMMG_IO_EdReq]) PMMG_MALLOC(parmesh, preq, nedreq, int, "preq", goto free_and_return);
    if (save_entities[PMMG_IO_EdPar]) PMMG_MALLOC(parmesh, ppar, nedpar, int, "ppar", goto free_and_return);

    crcount = reqcount = parcount = 0;

    if (na) {
      na = 0;
      for (int i = 0 ; i < mesh->na ; i++) {
        pa = &mesh->edge[i + 1];
        pent[2 * i]     = mesh->point[pa->a].tmp + offset[2 * PMMG_IO_Vertex] - 1;
        pent[2 * i + 1] = mesh->point[pa->b].tmp + offset[2 * PMMG_IO_Vertex] - 1;
        pref[i] = pa->ref;
        if (save_entities[PMMG_IO_Ridge] && (pa->tag & MG_GEO))    pcr[crcount++]   = na + offset[2 * PMMG_IO_Edge];
        if (save_entities[PMMG_IO_EdReq] && (pa->tag & MG_REQ))    preq[reqcount++] = na + offset[2 * PMMG_IO_Edge];
        if (save_entities[PMMG_IO_EdPar] && (pa->tag & MG_PARBDY)) ppar[parcount++] = na + offset[2 * PMMG_IO_Edge];
        na++;
      }
    }

    dspace_mem_id  = H5Screate_simple(2, hna, NULL);
    dspace_file_id = H5Screate_simple(2, hnag, NULL);
    H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Edge], NULL, hna, NULL);
    dset_id = H5Dcreate(grp_entities_id, "Edges", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
    H5Dclose(dset_id);
    H5Sclose(dspace_mem_id);
    H5Sclose(dspace_file_id);
    PMMG_DEL_MEM(parmesh, pent, int, "pent");

    dspace_mem_id  = H5Screate_simple(1, hna, NULL);
    dspace_file_id = H5Screate_simple(1, hnag, NULL);
    H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Edge], NULL, hna, NULL);
    dset_id = H5Dcreate(grp_entities_id, "EdgesRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
    H5Dclose(dset_id);
    H5Sclose(dspace_mem_id);
    H5Sclose(dspace_file_id);
    PMMG_DEL_MEM(parmesh, pref, int, "pref");

    if (save_entities[PMMG_IO_Ridge] && nrg) {
      dspace_mem_id  = H5Screate_simple(1, &nr, NULL);
      dspace_file_id = H5Screate_simple(1, &nrg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Ridge], NULL, &nr, NULL);
      dset_id = H5Dcreate(grp_entities_id, "Ridges", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
      H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pcr);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);
      PMMG_DEL_MEM(parmesh, pcr, int, "pcr");
    }

    if (save_entities[PMMG_IO_EdReq] && nedreqg) {
      dspace_mem_id  = H5Screate_simple(1, &nedreq, NULL);
      dspace_file_id = H5Screate_simple(1, &nedreqg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_EdReq], NULL, &nedreq, NULL);
      dset_id = H5Dcreate(grp_entities_id, "RequiredEdges", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
      H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);
      PMMG_DEL_MEM(parmesh, preq, int, "preq");
    }

    if (save_entities[PMMG_IO_EdPar] && nedparg) {
      dspace_mem_id  = H5Screate_simple(1, &nedpar, NULL);
      dspace_file_id = H5Screate_simple(1, &nedparg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_EdPar], NULL, &nedpar, NULL);
      dset_id = H5Dcreate(grp_entities_id, "ParallelEdges", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
      H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);
      PMMG_DEL_MEM(parmesh, ppar, int, "ppar");
    }
  }

  /* Triangles */
  if (save_entities[PMMG_IO_Tria] && nt) {

    PMMG_MALLOC(parmesh, pent, 3 * nt, int, "pent", goto free_and_return);
    PMMG_MALLOC(parmesh, pref, nt, int, "pref", goto free_and_return);
    if (save_entities[PMMG_IO_TriaReq]) PMMG_MALLOC(parmesh, preq, ntreq, int, "preq", goto free_and_return);
    if (save_entities[PMMG_IO_TriaPar]) PMMG_MALLOC(parmesh, ppar, ntpar, int, "ppar", goto free_and_return);

    reqcount = parcount = 0;

    if (nt) {
      nt = 0;
      for (int i = 0 ; i < mesh->nt ; i++) {
        pt = &mesh->tria[i + 1];
        for (int j = 0 ; j < 3 ; j++) {
          pent[3 * i + j] = mesh->point[pt->v[j]].tmp + offset[2 * PMMG_IO_Vertex] - 1;
        }
        pref[i] = pt->ref;
        if (save_entities[PMMG_IO_TriaReq]) {
          if (pt->tag[0] & MG_REQ && pt->tag[1] & MG_REQ && pt->tag[2] & MG_REQ) {
            preq[reqcount++] = nt + offset[2 * PMMG_IO_Tria];
          }
        }
        if (save_entities[PMMG_IO_TriaPar]) {
          if (pt->tag[0] & MG_PARBDY && pt->tag[1] & MG_PARBDY && pt->tag[2] & MG_PARBDY) {
            ppar[parcount++] = nt + offset[2 * PMMG_IO_Tria];
          }
        }
        nt++;
      }
    }

    dspace_mem_id  = H5Screate_simple(2, hnt, NULL);
    dspace_file_id = H5Screate_simple(2, hntg, NULL);
    H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Tria], NULL, hnt, NULL);
    dset_id = H5Dcreate(grp_entities_id, "Triangles", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
    H5Dclose(dset_id);
    H5Sclose(dspace_mem_id);
    H5Sclose(dspace_file_id);
    PMMG_DEL_MEM(parmesh, pent, int, "pent");

    dspace_mem_id  = H5Screate_simple(1, hnt, NULL);
    dspace_file_id = H5Screate_simple(1, hntg, NULL);
    H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Tria], NULL, hnt, NULL);
    dset_id = H5Dcreate(grp_entities_id, "TrianglesRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
    H5Dclose(dset_id);
    H5Sclose(dspace_mem_id);
    H5Sclose(dspace_file_id);
    PMMG_DEL_MEM(parmesh, pref, int, "pref");

    if (save_entities[PMMG_IO_TriaReq] && ntreqg) {
      dspace_mem_id  = H5Screate_simple(1, &ntreq, NULL);
      dspace_file_id = H5Screate_simple(1, &ntreqg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_TriaReq], NULL, &ntreq, NULL);
      dset_id = H5Dcreate(grp_entities_id, "RequiredTriangles", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
      H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);
      PMMG_DEL_MEM(parmesh, preq, int, "preq");
    }

    if (save_entities[PMMG_IO_TriaPar] && ntparg) {
      dspace_mem_id  = H5Screate_simple(1, &ntpar, NULL);
      dspace_file_id = H5Screate_simple(1, &ntparg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_TriaPar], NULL, &ntpar, NULL);
      dset_id = H5Dcreate(grp_entities_id, "ParallelTriangles", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
      H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);
      PMMG_DEL_MEM(parmesh, ppar, int, "ppar");
    }
  }


  /* Quadrilaterals */
  if (save_entities[PMMG_IO_Quad] && nquadg) {

    PMMG_MALLOC(parmesh, pent, 4 * nquad, int, "pent", goto free_and_return);
    PMMG_MALLOC(parmesh, pref, nquad, int, "pref", goto free_and_return);
    if (save_entities[PMMG_IO_QuadReq]) PMMG_MALLOC(parmesh, preq, nqreq, int, "preq", goto free_and_return);
    if (save_entities[PMMG_IO_QuadPar]) PMMG_MALLOC(parmesh, ppar, nqpar, int, "ppar", goto free_and_return);

    reqcount = parcount = 0;

    if (nquad){
      nquad = 0;
      for (int i = 0 ; i < mesh->nquad ; i++) {
        pq = &mesh->quadra[i + 1];
        for (int j = 0 ; j < 4 ; j++) {
          pent[4 * i + j] = mesh->point[pq->v[j]].tmp + offset[2 * PMMG_IO_Vertex] - 1;
        }
        pref[i] = pq->ref;
        if (save_entities[PMMG_IO_QuadReq]) {
          if (pq->tag[0] & MG_REQ && pq->tag[1] & MG_REQ &&
              pq->tag[2] & MG_REQ && pq->tag[3] & MG_REQ) {
            preq[reqcount++] = nquad + offset[2 * PMMG_IO_Quad];
          }
        }
        if (save_entities[PMMG_IO_QuadReq]) {
          if (pq->tag[0] & MG_PARBDY && pq->tag[1] & MG_PARBDY &&
              pq->tag[2] & MG_PARBDY && pq->tag[3] & MG_PARBDY) {
            ppar[parcount++] = nquad + offset[2 * PMMG_IO_Quad];
          }
        }
        nquad++;
      }
    }

    dspace_mem_id  = H5Screate_simple(2, hnquad, NULL);
    dspace_file_id = H5Screate_simple(2, hnquadg, NULL);
    H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Quad], NULL, hnquad, NULL);
    dset_id = H5Dcreate(grp_entities_id, "Quadrilaterals", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
    H5Dclose(dset_id);
    H5Sclose(dspace_mem_id);
    H5Sclose(dspace_file_id);
    PMMG_DEL_MEM(parmesh, pent, int, "pent");

    dspace_mem_id  = H5Screate_simple(1, hnquad, NULL);
    dspace_file_id = H5Screate_simple(1, hnquadg, NULL);
    H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Quad], NULL, hnquad, NULL);
    dset_id = H5Dcreate(grp_entities_id, "QuadrilateralsRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
    H5Dclose(dset_id);
    H5Sclose(dspace_mem_id);
    H5Sclose(dspace_file_id);
    PMMG_DEL_MEM(parmesh, pref, int, "pref");

    if (save_entities[PMMG_IO_QuadReq] && nqreqg) {
      dspace_mem_id  = H5Screate_simple(1, &nqreq, NULL);
      dspace_file_id = H5Screate_simple(1, &nqreqg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_QuadReq], NULL, &nqreq, NULL);
      dset_id = H5Dcreate(grp_entities_id, "RequiredQuadrilaterals", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
      H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);
      PMMG_DEL_MEM(parmesh, preq, int, "preq");
    }

    if (save_entities[PMMG_IO_QuadPar] && nqparg) {
      dspace_mem_id  = H5Screate_simple(1, &nqpar, NULL);
      dspace_file_id = H5Screate_simple(1, &nqparg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_QuadPar], NULL, &nqpar, NULL);
      dset_id = H5Dcreate(grp_entities_id, "ParallelQuadrilaterals", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
      H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);
      PMMG_DEL_MEM(parmesh, ppar, int, "ppar");
    }
  }

  /* Tetrahedra */
  if (save_entities[PMMG_IO_Tetra] && neg) {

    PMMG_MALLOC(parmesh, pent, 4 * ne, int, "pent", goto free_and_return);
    PMMG_MALLOC(parmesh, pref, ne, int, "pref", goto free_and_return);
    if (save_entities[PMMG_IO_TetReq]) PMMG_MALLOC(parmesh, preq, nereq, int, "preq", goto free_and_return);
    if (save_entities[PMMG_IO_TetPar]) PMMG_MALLOC(parmesh, ppar, nepar, int, "ppar", goto free_and_return);

    reqcount = parcount = 0;

    if (ne) {
      ne = 0;
      for (int i = 0 ; i < mesh->ne ; i++) {
        pe = &mesh->tetra[i + 1];
        if (!MG_EOK(pe)) continue;
        for (int j = 0 ; j < 4 ; j++) {
          pent[4 * ne + j] = mesh->point[pe->v[j]].tmp + offset[2 * PMMG_IO_Vertex] - 1;
        }
        pref[i] = pe->ref;
        if (save_entities[PMMG_IO_TetReq] && (pe->tag & MG_REQ))    preq[reqcount++] = ne + offset[2 * PMMG_IO_Tetra];
        if (save_entities[PMMG_IO_TetPar] && (pe->tag & MG_PARBDY)) ppar[parcount++] = ne + offset[2 * PMMG_IO_Tetra];
        ne++;
      }
    }

    dspace_mem_id  = H5Screate_simple(2, hne, NULL);
    dspace_file_id = H5Screate_simple(2, hneg, NULL);
    H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Tetra], NULL, hne, NULL);
    dset_id = H5Dcreate(grp_entities_id, "Tetrahedra", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
    H5Dclose(dset_id);
    H5Sclose(dspace_mem_id);
    H5Sclose(dspace_file_id);
    PMMG_DEL_MEM(parmesh, pent, int, "pent");

    dspace_mem_id  = H5Screate_simple(1, hne, NULL);
    dspace_file_id = H5Screate_simple(1, hneg, NULL);
    H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Tetra], NULL, hne, NULL);
    dset_id = H5Dcreate(grp_entities_id, "TetrahedraRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
    H5Dclose(dset_id);
    H5Sclose(dspace_mem_id);
    H5Sclose(dspace_file_id);
    PMMG_DEL_MEM(parmesh, pref, int, "pref");

    if (save_entities[PMMG_IO_TetReq] && nereqg) {
      dspace_mem_id  = H5Screate_simple(1, &nereq, NULL);
      dspace_file_id = H5Screate_simple(1, &nereqg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_TetReq], NULL, &nereq, NULL);
      dset_id = H5Dcreate(grp_entities_id, "RequiredTetrahedra", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
      H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);
      PMMG_DEL_MEM(parmesh, preq, int, "preq");
    }

    if (save_entities[PMMG_IO_TetPar] && neparg) {
      dspace_mem_id  = H5Screate_simple(1, &nepar, NULL);
      dspace_file_id = H5Screate_simple(1, &neparg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_TetPar], NULL, &nepar, NULL);
      dset_id = H5Dcreate(grp_entities_id, "ParallelTetrahedra", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
      H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);
      PMMG_DEL_MEM(parmesh, ppar, int, "ppar");
    }
  }

  /* Prisms */
  if (save_entities[PMMG_IO_Prism] && nprismg) {
    PMMG_MALLOC(parmesh, pent, 6 * nprism, int, "pent", goto free_and_return);
    PMMG_MALLOC(parmesh, pref, nprism, int, "pref", goto free_and_return);

    if (nprism){
      for (int i = 0 ; i < mesh->nprism ; i++) {
        pp = &mesh->prism[i + 1];
        for (int j = 0 ; j < 6 ; j++) {
          pent[6 * i + j] = mesh->point[pp->v[j]].tmp + offset[2 * PMMG_IO_Vertex] - 1;
        }
        pref[i] = pp->ref;
      }
    }

    dspace_mem_id  = H5Screate_simple(2, hnprism, NULL);
    dspace_file_id = H5Screate_simple(2, hnprismg, NULL);
    H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Prism], NULL, hnprism, NULL);
    dset_id = H5Dcreate(grp_entities_id, "Prisms", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
    H5Dclose(dset_id);
    H5Sclose(dspace_mem_id);
    H5Sclose(dspace_file_id);
    PMMG_DEL_MEM(parmesh, pent, int, "pent");

    dspace_mem_id  = H5Screate_simple(1, hnprism, NULL);
    dspace_file_id = H5Screate_simple(1, hnprismg, NULL);
    H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Prism], NULL, hnprism, NULL);
    dset_id = H5Dcreate(grp_entities_id, "PrismsRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
    H5Dclose(dset_id);
    H5Sclose(dspace_mem_id);
    H5Sclose(dspace_file_id);
    PMMG_DEL_MEM(parmesh, pref, int, "pref");
  }

  if (nullf) PMMG_DEL_MEM(parmesh, save_entities, int, "save_entities");

  /* Print the number of mesh entities */
  if ( parmesh->info.imprim > PMMG_VERB_STEPS ) {
    fprintf(stdout,"     NUMBER OF VERTICES       %lld   CORNERS %lld"
            "   REQUIRED %lld\n",npg,ncg,npreqg);
    fprintf(stdout,"     NUMBER OF TETRAHEDRA     %lld   REQUIRED  %lld\n",
            neg,nereqg);
    if ( nprismg )
      fprintf(stdout,"     NUMBER OF PRISMS         %lld\n",nprismg);

    if ( nag )
      fprintf(stdout,"     NUMBER OF EDGES          %lld   RIDGES  %lld"
              "   REQUIRED  %lld\n",nag,nrg,nedreqg);
    if ( ntg )
      fprintf(stdout,"     NUMBER OF TRIANGLES      %lld   REQUIRED  %lld\n",
              ntg, ntreqg);
    if ( nquadg )
      fprintf(stdout,"     NUMBER OF QUADRILATERALS %lld   REQUIRED  %lld\n",
              nquadg,nqreqg);

    if ( npparg || nedparg || ntparg || neparg || nqparg ) {
      fprintf(stdout,"     NUMBER OF PARALLEL ENTITIES: \n");
      if ( npparg )
        fprintf(stdout,"                  VERTICES       %lld \n",npparg);
      if ( nedparg )
        fprintf(stdout,"                  EDGES          %lld \n",nedparg);
      if ( ntparg )
        fprintf(stdout,"                  TRIANGLES      %lld \n",ntparg);
      if ( nqparg )
        fprintf(stdout,"                  QUADRILATERALS %lld \n",nqparg);
      if ( neparg )
        fprintf(stdout,"                  TETRAHEDRA    %lld \n",neparg);
    }
  }

  return 1;

 free_and_return:
  PMMG_DEL_MEM(parmesh, ppoint, double, "ppoint");
  PMMG_DEL_MEM(parmesh, pent, int, "pent");
  PMMG_DEL_MEM(parmesh, pref, int, "pref");
  PMMG_DEL_MEM(parmesh, pcr, int, "pcr");
  PMMG_DEL_MEM(parmesh, preq, int, "preq");
  PMMG_DEL_MEM(parmesh, ppar, int, "ppar");
  PMMG_DEL_MEM(parmesh, pnor, double, "pnor");
  PMMG_DEL_MEM(parmesh, pnorat, int, "pnorat");
  PMMG_DEL_MEM(parmesh, ptan, double, "ptan");
  PMMG_DEL_MEM(parmesh, ptanat, int, "ptanat");
  if (nullf) PMMG_DEL_MEM(parmesh, save_entities, int, "save_entities");
  return 0;

}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param grp_part_id identifier of the HDF5 group in which to write the mesh partitioning.
 * \param dcpl_id identifier of the dataset creation property list (no fill value).
 * \param dxpl_id identifier of the dataset transfer property list (MPI-IO).
 * \param nentities array of size nprocs * PMMG_NTYP_ENTITIES containing the number of entities on every proc.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Save the mesh partitioning in the \a grp_part_id group of an HDF5 file (only
 * one group per process is allowed).
 *
 */
static int PMMG_savePartitioning_hdf5(PMMG_pParMesh parmesh, hid_t grp_part_id, hid_t dcpl_id, hid_t dxpl_id, hsize_t *nentities) {
  PMMG_pExt_comm comms;
  hsize_t        *ncomms, *nitem, *nitem_proc;
  hsize_t        ncommg, comm_offset, nitemg, item_offset, rank_offset;
  hsize_t        icomm;
  int            *colors;
  int            **idx_loc, **idx_glob, *loc_buf, *glob_buf;
  int            ier;
  int            rank, nprocs;
  hid_t          dspace_mem_id, dspace_file_id;
  hid_t          dset_id;

  /* Set pointers to NULL */
  ncomms = nitem = nitem_proc = NULL;
  colors = NULL;
  idx_loc = idx_glob = NULL;
  loc_buf = glob_buf = NULL;

  /* Init variables */
  rank = parmesh->myrank;
  nprocs = parmesh->nprocs;

  if (parmesh->info.API_mode == PMMG_APIDISTRIB_faces)
    comms = parmesh->ext_face_comm;
  else
    comms = parmesh->ext_node_comm;

  ncommg = nitemg = comm_offset = item_offset = rank_offset = 0;

  /* Write the number of entities per proc */
  hsize_t hn[2] = {nprocs, PMMG_NTYPENTITIES};
  dspace_file_id = H5Screate_simple(2, hn, NULL);
  dset_id = H5Dcreate(grp_part_id, "NumberOfEntities", H5T_NATIVE_HSIZE, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  if (rank == parmesh->info.root)
    H5Dwrite(dset_id, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nentities);
  H5Dclose(dset_id);
  H5Sclose(dspace_file_id);

  /* Don't save the communicators outside the adaptation loop */
  if (parmesh->iter == PMMG_UNSET) return 1;

  /* Dont try to save communicators if there isn't any */
  if (nprocs == 1) return 1;

  /* Count the number of communicators */
  PMMG_MALLOC(parmesh, ncomms, nprocs, hsize_t, "ncomms", goto free_and_return);

  if (parmesh->info.API_mode == PMMG_APIDISTRIB_faces)
    ncomms[rank] = parmesh->next_face_comm;
  else
    ncomms[rank] = parmesh->next_node_comm;

  MPI_CHECK( MPI_Allgather(&ncomms[rank], 1, MPI_LONG_LONG, ncomms, 1, MPI_LONG_LONG, parmesh->comm),
             goto free_and_return );

  for (int i = 0 ; i < nprocs ; i++) {
    ncommg += ncomms[i];
  }

  /* Create the buffers */
  PMMG_MALLOC(parmesh, colors, ncomms[rank], int, "colors", goto free_and_return);
  PMMG_MALLOC(parmesh, nitem, ncomms[rank], hsize_t, "nitem", goto free_and_return);
  PMMG_CALLOC(parmesh, nitem_proc, nprocs, hsize_t, "nitem_proc", goto free_and_return);

  for (icomm = 0 ; icomm < ncomms[rank] ; icomm++) {
    colors[icomm] = comms[icomm].color_out;
    nitem[icomm] = comms[icomm].nitem;
    nitem_proc[rank] += nitem[icomm];
  }

  MPI_CHECK( MPI_Allgather(&nitem_proc[rank], 1, MPI_LONG_LONG, nitem_proc, 1, MPI_LONG_LONG, parmesh->comm),
             goto free_and_return );

  /* Count the total number of items */
  for (int i = 0 ; i < nprocs ; i++) {
    nitemg += nitem_proc[i];
  }

  /* Count the offsets for parallel writing */
  for (int i = 0 ; i < rank ; i++) {
    comm_offset += ncomms[i];
    rank_offset += nitem_proc[i];
  }

  /* Write the things */
  /* Number of communicators */
  hsize_t hnprocs = nprocs;
  dspace_file_id = H5Screate_simple(1, &hnprocs, NULL);
  if (parmesh->info.API_mode == PMMG_APIDISTRIB_faces)
    dset_id = H5Dcreate(grp_part_id, "NumberOfFaceCommunicators", H5T_NATIVE_HSIZE, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  else
    dset_id = H5Dcreate(grp_part_id, "NumberOfNodeCommunicators", H5T_NATIVE_HSIZE, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  if (rank == parmesh->info.root)
    H5Dwrite(dset_id, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ncomms);
  H5Dclose(dset_id);
  H5Sclose(dspace_file_id);

  /* For each communicator, write the outward proc color and the number of faces */
  dspace_mem_id  = H5Screate_simple(1, &ncomms[rank], NULL);
  dspace_file_id = H5Screate_simple(1, &ncommg, NULL);
  H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &comm_offset, NULL, &ncomms[rank], NULL);

  dset_id = H5Dcreate(grp_part_id, "ColorsOut", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, colors);
  H5Dclose(dset_id);

  if (parmesh->info.API_mode == PMMG_APIDISTRIB_faces)
    dset_id = H5Dcreate(grp_part_id, "NumberOfCommunicatorFaces", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  else
    dset_id = H5Dcreate(grp_part_id, "NumberOfCommunicatorNodes", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);

  H5Dwrite(dset_id, H5T_NATIVE_HSIZE, dspace_mem_id, dspace_file_id, dxpl_id, nitem);
  H5Dclose(dset_id);

  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);

  /* Get the communicator items */
  PMMG_CALLOC(parmesh, idx_loc, ncomms[rank], int*, "idx_loc", goto free_and_return);
  PMMG_CALLOC(parmesh, idx_glob, ncomms[rank], int*, "idx_glob", goto free_and_return);
  for (icomm = 0 ; icomm < ncomms[rank] ; icomm++) {
    PMMG_CALLOC(parmesh, idx_loc[icomm], nitem[icomm], int, "idx_loc[icomm]", goto free_and_return);
    PMMG_CALLOC(parmesh, idx_glob[icomm], nitem[icomm], int, "idx_glob[icomm]", goto free_and_return);
  }

  if (parmesh->info.API_mode == PMMG_APIDISTRIB_faces) {
    ier = PMMG_Get_FaceCommunicator_faces(parmesh, idx_loc);

    if ( !ier ) {
      fprintf(stderr,"\n  ## Error: %s: unable to compute face communicators on rank %d.\n",
              __func__, rank);
      return ier;
    }
    ier = PMMG_Get_FaceCommunicator_owners(parmesh, NULL, idx_glob, NULL, NULL);
  }
  else {
    ier = PMMG_Get_NodeCommunicator_nodes(parmesh, idx_loc);

    if ( !ier ) {
      fprintf(stderr,"\n  ## Error: %s: unable to compute node communicators on rank %d.\n",
              __func__, rank);
      return ier;
    }
    ier = PMMG_Get_NodeCommunicator_owners(parmesh, NULL, idx_glob, NULL, NULL);
  }

  /* Make a unique buffer for each proc */
  PMMG_CALLOC(parmesh, loc_buf, nitem_proc[rank], int, "loc_buf", goto free_and_return);
  PMMG_CALLOC(parmesh, glob_buf, nitem_proc[rank], int, "glob_buf", goto free_and_return);

  for (icomm = 0 ; icomm < ncomms[rank] ; icomm++) {
    for (int k = 0 ; k < nitem[icomm] ; k++) {
      loc_buf[item_offset + k] = idx_loc[icomm][k];
      glob_buf[item_offset + k] = idx_glob[icomm][k];
    }
    item_offset += nitem[icomm];
  }

  /* Free the memory of idx_loc and idx_glob */
  for (icomm = 0 ; icomm < ncomms[rank] ; icomm++) {
    PMMG_DEL_MEM(parmesh, idx_loc[icomm], int, "idx_loc[icomms]");
    PMMG_DEL_MEM(parmesh, idx_glob[icomm], int, "idx_glob[icomm]");
  }
  PMMG_DEL_MEM(parmesh, idx_loc, int*, "idx_loc");
  PMMG_DEL_MEM(parmesh, idx_glob, int*, "idx_glob");

  /* Write the local indices */
  dspace_file_id = H5Screate_simple(1, &nitemg, NULL);

  if (parmesh->info.API_mode == PMMG_APIDISTRIB_faces)
    dset_id = H5Dcreate(grp_part_id, "LocalFaceIndices", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  else
    dset_id = H5Dcreate(grp_part_id, "LocalNodeIndices", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  dspace_mem_id = H5Screate_simple(1, &nitem_proc[rank], NULL);
  H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &rank_offset, NULL, &nitem_proc[rank], NULL);
  H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, loc_buf);
  H5Sclose(dspace_mem_id);
  H5Dclose(dset_id);

  /* Write the global indices */
  if (parmesh->info.API_mode == PMMG_APIDISTRIB_faces)
    dset_id = H5Dcreate(grp_part_id, "GlobalFaceIndices", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  else
    dset_id = H5Dcreate(grp_part_id, "GlobalNodeIndices", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  dspace_mem_id = H5Screate_simple(1, &nitem_proc[rank], NULL);
  H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &rank_offset, NULL, &nitem_proc[rank], NULL);
  H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, glob_buf);
  H5Sclose(dspace_mem_id);
  H5Dclose(dset_id);

  H5Sclose(dspace_file_id);

  /* Free the memory */
  PMMG_DEL_MEM(parmesh, loc_buf, int, "loc_buf");
  PMMG_DEL_MEM(parmesh, glob_buf, int, "glob_buf");
  PMMG_DEL_MEM(parmesh, ncomms, hsize_t, "ncomms");
  PMMG_DEL_MEM(parmesh, colors, int, "colors");
  PMMG_DEL_MEM(parmesh, nitem, int, "nitem");
  PMMG_DEL_MEM(parmesh, nitem_proc, hsize_t, "nitem_proc");

  return 1;

 free_and_return:
  PMMG_DEL_MEM(parmesh, loc_buf, int, "loc_buf");
  PMMG_DEL_MEM(parmesh, glob_buf, int, "glob_buf");
  PMMG_DEL_MEM(parmesh, ncomms, hsize_t, "ncomms");
  PMMG_DEL_MEM(parmesh, colors, int, "colors");
  PMMG_DEL_MEM(parmesh, nitem, int, "nitem");
  PMMG_DEL_MEM(parmesh, nitem_proc, hsize_t, "nitem_proc");

  if ( idx_loc && idx_glob ) {
    for (icomm = 0 ; icomm < ncomms[rank] ; icomm++) {
      PMMG_DEL_MEM(parmesh, idx_loc[icomm], int, "idx_loc[icomms]");
      PMMG_DEL_MEM(parmesh, idx_glob[icomm], int, "idx_glob[icomm]");
    }
    PMMG_DEL_MEM(parmesh, idx_loc, int*, "idx_loc");
    PMMG_DEL_MEM(parmesh, idx_glob, int*, "idx_glob");
  }

  return 0;

}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param grp_sols_id identifier of the HDF5 group in which to write the metric.
 * \param dcpl_id identifier of the dataset creation property list (no fill value).
 * \param dxpl_id identifier of the dataset transfer property list (MPI-IO).
 * \param nentitiesl array of size PMMG_NTYP_ENTITIES containing the local number of entities.
 * \param nentitiesg array of size PMMG_NTYP_ENTITIES containing the global number of entities.
 * \param offset array of size PMMG_NTYP_ENTITIES containing the offset for parallel writing.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Save the metric in the \a grp_sols_id group of an HDF5 file (only one group per process is allowed).
 *
 */
static int PMMG_saveMetric_hdf5(PMMG_pParMesh parmesh, hid_t grp_sols_id, hid_t dcpl_id, hid_t dxpl_id,
                                hsize_t *nentitiesl, hsize_t *nentitiesg, hsize_t *offset) {
  int         np, npg, mcount, isMet;
  MMG5_pMesh  mesh;
  MMG5_pSol   met;
  MMG5_pPoint ppt;
  double      *met_buf;
  hsize_t     met_offset[2] = {0, 0};
  hid_t       dspace_mem_id, dspace_file_id;
  hid_t       dset_id;

  assert ( parmesh->ngrp == 1 );

  mesh = parmesh->listgrp[0].mesh;
  met = parmesh->listgrp[0].met;
  np  = nentitiesl[PMMG_IO_Vertex];
  npg = nentitiesg[PMMG_IO_Vertex];

  /* Check the metric */
  isMet = (met && met->m);

  MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, &isMet, 1, MPI_INT, MPI_MAX, parmesh->comm), return 0);

  if (!isMet) {
    if (parmesh->myrank == parmesh->info.root) {
      fprintf(stderr, "\n  ## Warning: %s: No metric data to save \n", __func__);
    }
    return 1;
  }
  if (met->size != 1 && met->size != 6) {
    fprintf(stderr, "\n  ## Error: %s: Wrong metric size\n", __func__);
    return 0;
  }
  if (mesh->np != met->np) {
    fprintf(stderr, "\n  ## Error: %s: The metric vertices do not match with the mesh vertices \n", __func__);
    return 0;
  }

  /* Arrays for bidimensional dataspaces */
  hsize_t hns[2]  = {np, met->size};
  hsize_t hnsg[2] = {npg, met->size};

  /* Offset for parallel writing */
  met_offset[0] = offset[2 * PMMG_IO_Vertex];

  /* Fill the metric buffer */
  PMMG_MALLOC(parmesh, met_buf, np * met->size, double, "met_buf", return 0);
  mcount = 0;
  for (int k = 0 ; k < mesh->np ; k++) {
    ppt = &mesh->point[k + 1];
    if (!MG_VOK(ppt)) continue;
    for (int j = 0 ; j < met->size ; j++) {
      met_buf[mcount++] = met->m[1 + k * met->size + j];
    }
  }

  /* Write the buffer */
  dspace_mem_id = H5Screate_simple(2, hns, NULL);
  dspace_file_id = H5Screate_simple(2, hnsg, NULL);
  H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, met_offset, NULL, hns, NULL);
  dset_id = H5Dcreate(grp_sols_id, "MetricAtVertices", H5T_NATIVE_DOUBLE, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dspace_mem_id, dspace_file_id, dxpl_id, met_buf);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);

  /* Free the memory */
  PMMG_DEL_MEM(parmesh, met_buf, double, "met_buf");

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param grp_sols_id identifier of the HDF5 group in which to write the solutions.
 * \param dcpl_id identifier of the dataset creation property list (no fill value).
 * \param dxpl_id identifier of the dataset transfer property list (MPI-IO).
 * \param nentitiesl array of size PMMG_NTYP_ENTITIES containing the local number of entities.
 * \param nentitiesg array of size PMMG_NTYP_ENTITIES containing the global number of entities.
 * \param offset array of size PMMG_NTYP_ENTITIES containing the offset for parallel writing.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Save the solutions in the \a grp_sols_id group of an HDF5 file (only one group per process is allowed).
 *
 */
static int PMMG_saveAllSols_hdf5(PMMG_pParMesh parmesh, hid_t grp_sols_id, hid_t dcpl_id, hid_t dxpl_id,
                                 hsize_t *nentitiesl, hsize_t *nentitiesg, hsize_t *offset) {
  MMG5_pMesh  mesh;
  MMG5_pSol   *sols;
  MMG5_pPoint ppt;
  MMG5_pTetra pt;
  int         nsols, ndigits, np, npg, ne, neg, size, count, vcount, tcount, iwar = 1;
  char        *solname, *tmp;
  hsize_t     *sol_offset;
  double      *sol_buf;
  hid_t       dspace_mem_id, dspace_file_id;
  hid_t       dset_id, attr_id;

  /* Set ParMmg variables */
  mesh = parmesh->listgrp[0].mesh;
  sols = &(parmesh->listgrp[0].field);
  nsols = mesh->nsols;

  /* Get the local and global number of vertices/tetra */
  np  = nentitiesl[PMMG_IO_Vertex];
  ne  = nentitiesl[PMMG_IO_Tetra];
  npg = nentitiesg[PMMG_IO_Vertex];
  neg = nentitiesg[PMMG_IO_Tetra];

  /* Arrays for bidimensional dataspaces */
  hsize_t hns[2]  = {0, 0};
  hsize_t hnsg[2] = {0, 0};

  /* Count digits for the name of the datasets */
  ndigits = PMMG_count_digits(nsols);

  /* Initialize the counts */
  vcount = tcount = 0;

  for (int i = 0 ; i < nsols ; i++) {

    if ( !sols[i] || !sols[i]->m ) {
      iwar = 0;
    }
    MPI_Allreduce(MPI_IN_PLACE, &iwar, 1, MPI_INT, MPI_MIN, parmesh->comm);
    if ( !iwar ) {
      fprintf(stderr, "\n  ## Warning: %s: Skipping empty solution %d.\n", __func__, i);
      continue;
    }

    size = sols[i]->size;
    count = 0;

    if (sols[i]->entities == MMG5_Noentity || sols[i]->entities == MMG5_Vertex) {
      PMMG_MALLOC(parmesh, sol_buf, size * np, double, "sol_buf", goto free_buf);
      for (int k = 0 ; k < mesh->np ; k++) {
        ppt = &mesh->point[k + 1];
        if ( !MG_VOK(ppt) ) continue;
        for (int j = 0 ; j < size ; j++) {
          sol_buf[count++] = sols[i]->m[1 + k * size + j];
        }
      }
      hns[0] = np; hns[1] = size;
      hnsg[0] = npg; hnsg[1] = size;
      PMMG_CALLOC(parmesh, sol_offset, np * size, hsize_t, "sol_offset", goto free_buf);
      sol_offset[0] = offset[2 * PMMG_IO_Vertex];

      PMMG_CALLOC(parmesh, solname, strlen("SolAtVertices") + ndigits + 1, char, "solname", goto free_buf);
      PMMG_CALLOC(parmesh, tmp, ndigits + 1, char, "tmp", goto free_buf);
      strcpy(solname, "SolAtVertices");
      sprintf(tmp, "%d", vcount);
      strcat(solname, tmp);
      vcount++;
    }

    else if (sols[i]->entities == MMG5_Tetrahedron) {
      PMMG_MALLOC(parmesh, sol_buf, size * ne, double, "sol_buf", goto free_buf);
      for (int k = 0 ; k < mesh->ne ; k++) {
        pt = &mesh->tetra[k + 1];
        if ( !MG_EOK(pt) ) continue;
        for (int j = 0 ; j < size ; j++) {
          sol_buf[count++] = sols[i]->m[1 + k * size + j];
        }
      }
      hns[0] = ne; hns[1] = size;
      hnsg[0] = neg; hnsg[1] = size;
      PMMG_CALLOC(parmesh, sol_offset, ne * size, hsize_t, "sol_offset", goto free_buf);
      sol_offset[0] = offset[2 * PMMG_IO_Tetra];

      PMMG_CALLOC(parmesh, solname, strlen("SolAtTetrahedra") + ndigits + 1, char, "solname", goto free_buf);
      PMMG_CALLOC(parmesh, tmp, ndigits + 1, char, "tmp", goto free_buf);
      strcpy(solname, "SolAtTetrahedra");
      sprintf(tmp, "%d", vcount);
      strcat(solname, tmp);
      tcount++;
    }

    else {
      printf("\n  ## Warning: %s: unexpected entity type for solution %d: %s."
             "\n Ignored.\n",
             __func__, i, MMG5_Get_entitiesName(sols[i]->entities));
      continue;
    }

    dspace_mem_id = H5Screate_simple(2, hns, NULL);
    dspace_file_id = H5Screate_simple(2, hnsg, NULL);
    H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, sol_offset, NULL, hns, NULL);
    dset_id = H5Dcreate(grp_sols_id, solname, H5T_NATIVE_DOUBLE, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dspace_mem_id, dspace_file_id, dxpl_id, sol_buf);
    H5Dclose(dset_id);
    H5Sclose(dspace_mem_id);
    H5Sclose(dspace_file_id);

  free_buf:
    PMMG_DEL_MEM(parmesh, sol_offset, hsize_t, "sol_offset");
    PMMG_DEL_MEM(parmesh, sol_buf, double, "sol_buf");
    PMMG_DEL_MEM(parmesh, solname, char, "solname");
    PMMG_DEL_MEM(parmesh, tmp, char, "tmp");
  }

  /* Save the actual number of solutions as group attributes */
  dspace_file_id = H5Screate(H5S_SCALAR);

  attr_id = H5Acreate(grp_sols_id, "NSolsAtVertices", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr_id, H5T_NATIVE_INT, &vcount);
  H5Aclose(attr_id);

  attr_id = H5Acreate(grp_sols_id, "NSolsAtTetrahedra", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr_id, H5T_NATIVE_INT, &tcount);
  H5Aclose(attr_id);

  H5Sclose(dspace_file_id);

  return 1;
}


/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of the HDF5 file containing the heavy data.
 * \param xdmfname name of the XDMF file that will contain the light data.
 * \param nentitiesg array of size PMMG_NTYP_ENTITIES containing the global number of entities.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Create the XDMF file \a xdmfname and write light data describing the mesh that was saved in
 * the HDF5 file \a filename (only one group per process is allowed).
 *
 */
static int PMMG_writeXDMF(PMMG_pParMesh parmesh, const char *filename, const char *xdmfname, hsize_t *nentitiesg) {
  hsize_t neg, npg;
  PMMG_pGrp grp;
  MMG5_pSol met, *sols;
  int nsols, entities;
  FILE *xdmf_file = NULL;

  assert ( parmesh->ngrp == 1 );

  npg  = nentitiesg[PMMG_IO_Vertex];
  neg  = nentitiesg[PMMG_IO_Tetra];
  grp  = &parmesh->listgrp[0];
  met  = grp->met;
  sols = &grp->field;
  nsols = grp->mesh->nsols;

  if (parmesh->myrank == parmesh->info.root) {

    xdmf_file = fopen(xdmfname, "w");

    if ( !xdmf_file ) return 0;

    /* XDMF header */
    fprintf(xdmf_file, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
    fprintf(xdmf_file, "<Xdmf Version=\"3.0\">\n");
    fprintf(xdmf_file, "<Domain>\n");
    fprintf(xdmf_file, "    <Grid Name=\"3D Unstructured Mesh\" GridType=\"Uniform\">\n");

    /* Tetrahedra */
    fprintf(xdmf_file, "      <Topology TopologyType=\"Tetrahedron\" NumberOfElements=\"%llu\">\n", neg);
    fprintf(xdmf_file, "        <DataItem DataType=\"Int\"\n");
    fprintf(xdmf_file, "                  Format=\"HDF\"\n");
    fprintf(xdmf_file, "                  Dimensions=\"%llu 4\">\n", neg);
    fprintf(xdmf_file, "          %s:/Mesh/MeshEntities/Tetrahedra\n", filename);
    fprintf(xdmf_file, "        </DataItem>\n");
    fprintf(xdmf_file, "      </Topology>\n");

    /* Vertices */
    fprintf(xdmf_file, "      <Geometry GeometryType=\"XYZ\">\n");
    fprintf(xdmf_file, "        <DataItem DataType=\"Float\"\n");
    fprintf(xdmf_file, "                  Precision=\"8\"\n");
    fprintf(xdmf_file, "                  Format=\"HDF\"\n");
    fprintf(xdmf_file, "                  Dimensions=\"%llu 3\">\n", npg);
    fprintf(xdmf_file, "          %s:/Mesh/MeshEntities/Vertices\n", filename);
    fprintf(xdmf_file, "        </DataItem>\n");
    fprintf(xdmf_file, "      </Geometry>\n");

    /* Metric */
    if (met && met->m) {
      if (met->size == 6)
        fprintf(xdmf_file, "      <Attribute Center=\"Node\" Name=\"Metric\" AttributeType=\"Tensor6\">\n");
      else if (met->size == 1)
        fprintf(xdmf_file, "      <Attribute Center=\"Node\" Name=\"Metric\" AttributeType=\"Scalar\">\n");
      fprintf(xdmf_file, "        <DataItem DataType=\"Float\"\n");
      fprintf(xdmf_file, "                  Precision=\"8\"\n");
      fprintf(xdmf_file, "                  Format=\"HDF\"\n");
      fprintf(xdmf_file, "                  Dimensions=\"%lld %d\">\n", npg, met->size);
      fprintf(xdmf_file, "          %s:/Solutions/MetricAtVertices\n", filename);
      fprintf(xdmf_file, "        </DataItem>\n");
      fprintf(xdmf_file, "      </Attribute>\n");
    }

    /* Solutions */
    for (int i = 0 ; i < nsols ; i++) {

      /* Ignore invalid solutions */
      if ( !sols[i] || !sols[i]->m ) continue;

      entities = sols[i]->entities;

      if (entities != MMG5_Noentity && entities != MMG5_Vertex && entities != MMG5_Tetrahedron) continue;
      if (sols[i]->type == MMG5_Scalar) {
        fprintf(xdmf_file, "      <Attribute Center=\"Node\" Name=\"Sol%d\" AttributeType=\"Scalar\">\n", i);
      }
      else if (sols[i]->type == MMG5_Vector) {
        fprintf(xdmf_file, "      <Attribute Center=\"Node\" Name=\"Sol%d\" AttributeType=\"Vector\">\n", i);
      }
      else if (sols[i]->type == MMG5_Tensor) {
        fprintf(xdmf_file, "      <Attribute Center=\"Node\" Name=\"Sol%d\" AttributeType=\"Tensor\">\n", i);
      }
      fprintf(xdmf_file, "        <DataItem DataType=\"Float\"\n");
      fprintf(xdmf_file, "                  Precision=\"8\"\n");
      fprintf(xdmf_file, "                  Format=\"HDF\"\n");
      if (sols[i]->entities == MMG5_Noentity || sols[i]->entities == MMG5_Vertex) {
        fprintf(xdmf_file, "                  Dimensions=\"%lld %d\">\n", npg, sols[i]->size);
        fprintf(xdmf_file, "          %s:/Solutions/SolAtVertices%d\n", filename, i);
      }
      else if (sols[i]->entities == MMG5_Tetrahedron) {
        fprintf(xdmf_file, "                  Dimensions=\"%lld %d\">\n", neg, sols[i]->size);
        fprintf(xdmf_file, "          %s:/Solutions/SolAtTetrahedra%d\n", filename, i);
      }
      fprintf(xdmf_file, "        </DataItem>\n");
      fprintf(xdmf_file, "      </Attribute>\n");
    }

    /* End */
    fprintf(xdmf_file, "    </Grid>\n");
    fprintf(xdmf_file, "  </Domain>\n");
    fprintf(xdmf_file, "</Xdmf>\n");

    fclose(xdmf_file);
  }

  return 1;
}

int PMMG_saveParmesh_hdf5(PMMG_pParMesh parmesh, int *save_entities, const char *filename, const char *xdmfname) {

#ifndef USE_HDF5

  fprintf(stderr,"  ** HDF5 library not found. Unavailable file format.\n");
  return -1;

#else

  int        ier = 1;
  hsize_t    *nentities, *nentitiesl, *nentitiesg; /* Number of entities (on each proc/on the current proc/global) */
  hsize_t    *offset;                              /* Offset for the parallel writing with HDF5 */
  hid_t      file_id, grp_mesh_id, grp_part_id, grp_entities_id, grp_sols_id; /* HDF5 objects */
  hid_t      fapl_id, dxpl_id, dcpl_id;                                       /* HDF5 property lists */
  MPI_Info   info = MPI_INFO_NULL;
  mytime     ctim[TIMEMAX];
  int8_t     tim;
  char       stim[32];

  /* Check arguments */
  if (parmesh->ngrp != 1) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in your parmesh.\n",
            __func__);
    ier = 0;
  }
  if (!filename || !*filename) {
    fprintf(stderr,"  ## Error: %s: no HDF5 file name provided.\n",
            __func__);
    ier = 0;
  }
  if (!save_entities[PMMG_IO_Vertex] || !save_entities[PMMG_IO_Tetra]) {
    fprintf(stderr, "\n  ## Error: %s: save_entities: you must at least save the vertices and the tetra.\n",
            __func__);
    ier = 0;
  }

  MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, &ier, 1, MPI_INT, MPI_MIN, parmesh->comm ), return 0);
  if ( !ier ) {
    return 0;
  }

  tminit(ctim, TIMEMAX);
  chrono(ON, &ctim[0]);

  /* Set all pointers to NULL */
  nentities = nentitiesl = nentitiesg = offset = NULL;

  /** Count the number of entities on each proc and globally */
  tim = 1;
  chrono(ON, &ctim[tim]);

  PMMG_CALLOC(parmesh, nentities, PMMG_NTYPENTITIES * parmesh->nprocs, hsize_t, "nentities",
              goto free_and_return);
  PMMG_CALLOC(parmesh, nentitiesg, PMMG_NTYPENTITIES, hsize_t, "nentitiesg",
              goto free_and_return );
  PMMG_CALLOC(parmesh, nentitiesl, PMMG_NTYPENTITIES, hsize_t, "nentitiesl",
              goto free_and_return );

  ier = PMMG_countEntities(parmesh, nentities, nentitiesl, nentitiesg, save_entities);

  /* Check that the global mesh is not empty */
  if ( !nentitiesg[PMMG_IO_Vertex] ) {
    if (parmesh->myrank == parmesh->info.root)
      fprintf(stderr, "\n  ## Error: %s: can't save an empty mesh.\n", __func__);
    return 0;
  }
  if ( !nentitiesg[PMMG_IO_Tetra] ) {
    if (parmesh->myrank == parmesh->info.root)
      fprintf(stderr, "\n  ## Warning: %s: there is no tetra in your mesh.\n", __func__);
  }

  chrono(OFF, &ctim[tim]);
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    fprintf(stdout,"  -- Count entities.     %s\n",stim);
  }

  /** Compute the offset for parallel writing */
  tim = 2;
  chrono(ON, &ctim[tim]);

  PMMG_CALLOC(parmesh, offset, 2 * PMMG_NTYPENTITIES, hsize_t, "offset",
              goto free_and_return );

  ier = PMMG_computeHDFoffset(parmesh, nentities, offset);

  chrono(OFF, &ctim[tim]);
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    fprintf(stdout,"  -- Compute HDF5 write offsets.     %s\n",stim);
  }

  /*------------------------- HDF5 IOs START HERE -------------------------*/

  /* Shut HDF5 error stack */
  HDF_CHECK( H5Eset_auto(H5E_DEFAULT, NULL, NULL),
             goto free_and_return );

  /* TODO ? Pass MPI hints via the info struct */
  MPI_Info_create(&info);

  /* Create the property lists */
  fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl_id, parmesh->comm, info);
  H5Pset_coll_metadata_write(fapl_id, 1);
  dxpl_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
  dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER);

  /* Create the file */
  HDF_CHECK( file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id),
             fprintf(stderr,"\n  ## Error: %s: Rank %d could not create the hdf5 file %s.\n",
                     __func__, parmesh->myrank, filename);
             goto free_and_return );

  if (parmesh->info.imprim > PMMG_VERB_VERSION) {
    fprintf(stdout, "\n  %%%% %s OPENED \n", filename);
  }

  /* Save the attributes (Version and Dimension, and number of entities per proc) */
  ier = PMMG_saveHeader_hdf5(parmesh, file_id);

  MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, &ier, 1, MPI_INT, MPI_MIN, parmesh->comm),
             goto free_and_return );

  if (!ier) {
    if (parmesh->myrank == parmesh->info.root) {
      fprintf(stderr,"\n  ## Error: %s: Could not write the mesh attributes.\n",__func__);
    }
    goto free_and_return;
  }

  /* Open the mesh group */
  HDF_CHECK( grp_mesh_id = H5Gcreate(file_id, "Mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
             fprintf(stderr,"\n  ## Error: %s: Could not create the /Mesh group.\n",__func__);
             goto free_and_return );

  /** Write the partitionning information */
  tim = 3;
  chrono(ON, &ctim[tim]);

  HDF_CHECK( grp_part_id = H5Gcreate(grp_mesh_id, "Partitioning", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
             fprintf(stderr,"\n  ## Error: %s: Could not create the /Mesh/Partitioning group.\n",__func__);
             H5Gclose(grp_mesh_id);
             goto free_and_return );

  ier = PMMG_savePartitioning_hdf5(parmesh, grp_part_id, dcpl_id, dxpl_id, nentities);

  MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, &ier, 1, MPI_INT, MPI_MIN, parmesh->comm),
             H5Gclose(grp_part_id);
             H5Gclose(grp_mesh_id);
             goto free_and_return );

  if (!ier) {
    if (parmesh->myrank == parmesh->info.root) {
      fprintf(stderr,"\n  ## Error: %s: Could not write the partitioning information.\n",__func__);
    }
    H5Gclose(grp_part_id);
    H5Gclose(grp_mesh_id);
    goto free_and_return;
  }

  H5Gclose(grp_part_id);

  chrono(OFF, &ctim[tim]);
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    fprintf(stdout,"  -- Save partitioning info.     %s\n",stim);
  }

  /* Now that we have computed the offsets and written the number of entities,
     each proc only needs to know its local number of entities and the
     global number of entities */
  PMMG_DEL_MEM(parmesh, nentities, hsize_t, "nentities");

  /** Write the mesh entities */
  tim = 4;
  chrono(ON, &ctim[tim]);

  HDF_CHECK( grp_entities_id = H5Gcreate(grp_mesh_id, "MeshEntities", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
             fprintf(stderr,"\n  ## Error: %s: Could not create the /Mesh/MeshEntities group.\n",__func__);
             H5Gclose(grp_mesh_id);
             goto free_and_return );

  ier = PMMG_saveMeshEntities_hdf5(parmesh, grp_entities_id, dcpl_id, dxpl_id, nentitiesl, nentitiesg, offset, save_entities);

  MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, &ier, 1, MPI_INT, MPI_MIN, parmesh->comm),
             H5Gclose(grp_entities_id);
             H5Gclose(grp_mesh_id);
             goto free_and_return );

  if (!ier) {
    if (parmesh->myrank == parmesh->info.root) {
      fprintf(stderr,"\n  ## Error: %s: Could not write the mesh entities.\n",__func__);
    }
    H5Gclose(grp_entities_id);
    H5Gclose(grp_mesh_id);
    goto free_and_return;
  }

  H5Gclose(grp_entities_id);

  chrono(OFF, &ctim[tim]);
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    fprintf(stdout,"  -- Save mesh entities.     %s\n",stim);
  }

  /* Close the mesh group */
  H5Gclose(grp_mesh_id);

  /** Write the metric and the solutions */
  tim = 5;
  chrono(ON, &ctim[tim]);

  HDF_CHECK( grp_sols_id = H5Gcreate(file_id, "Solutions", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
             fprintf(stderr,"\n  ## Error: %s: Could not create the /Solutions group.\n",__func__);
             goto free_and_return );

  ier = PMMG_saveMetric_hdf5(parmesh, grp_sols_id, dcpl_id, dxpl_id, nentitiesl, nentitiesg, offset);

  MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, &ier, 1, MPI_INT, MPI_MIN, parmesh->comm),
             H5Gclose(grp_sols_id);
             goto free_and_return );

  if (!ier) {
    if (parmesh->myrank == parmesh->info.root) {
      fprintf(stderr,"\n  ## Error: %s: Could not write the metric.\n",__func__);
    }
    H5Gclose(grp_sols_id);
    goto free_and_return;
  }

  ier = PMMG_saveAllSols_hdf5(parmesh, grp_sols_id, dcpl_id, dxpl_id, nentitiesl, nentitiesg, offset);

  MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, &ier, 1, MPI_INT, MPI_MIN, parmesh->comm),
             H5Gclose(grp_sols_id);
             goto free_and_return );

  if (!ier) {
    if (parmesh->myrank == parmesh->info.root) {
      fprintf(stderr,"\n  ## Error: %s: Could not write the solutions.\n",__func__);
    }
    H5Gclose(grp_sols_id);
    goto free_and_return;
  }

  H5Gclose(grp_sols_id);

  chrono(OFF, &ctim[tim]);
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    fprintf(stdout,"  -- Save metric and solutions.     %s\n",stim);
  }

  /* Release all HDF5 IDs */
  H5Fclose(file_id);
  H5Pclose(fapl_id);
  H5Pclose(dxpl_id);
  H5Pclose(dcpl_id);

  if (parmesh->info.imprim > PMMG_VERB_VERSION) {
    fprintf(stdout, "  %%%% %s CLOSED \n", filename);
  }

  /* We no longer need the offset nor the local nuumber of entities */
  PMMG_DEL_MEM(parmesh, offset, hsize_t, "offset");
  PMMG_DEL_MEM(parmesh, nentitiesl, hsize_t, "nentitiesl");

  /** Write light data in XDMF file */
  if (!xdmfname || !*xdmfname)
    fprintf(stderr,"  ## Warning: %s: no XDMF file name provided.", __func__);
  else {
    ier = PMMG_writeXDMF(parmesh, filename, xdmfname, nentitiesg);

    MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, &ier, 1, MPI_INT, MPI_MIN, parmesh->comm),
              goto free_and_return );
    if (!ier) {
      if (parmesh->myrank == parmesh->info.root) {
        fprintf(stderr,"\n  ## Error: %s: Could not write the xdmf file %s.\n",
                __func__, xdmfname);
      }
      goto free_and_return;
    }
  }

  /* We no longer need the global number of entities */
  PMMG_DEL_MEM(parmesh, nentitiesg, hsize_t, "nentitiesg");

  chrono(OFF, &ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    fprintf(stdout,"\n   SAVE_PARMESH_HDF5: ELAPSED TIME  %s\n",stim);
  }
  return 1;

 free_and_return:
  PMMG_DEL_MEM(parmesh, nentities, hsize_t, "nentities");
  PMMG_DEL_MEM(parmesh, nentitiesg, hsize_t, "nentitiesg");
  PMMG_DEL_MEM(parmesh, nentitiesl, hsize_t, "nentitiesl");
  PMMG_DEL_MEM(parmesh, offset, hsize_t, "offset");
  H5Fclose(file_id);
  H5Pclose(fapl_id);
  H5Pclose(dxpl_id);
  H5Pclose(dcpl_id);
  return 0;

#endif
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param file_id identifier of the HDF5 file.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Load the version and the dimension of \a parmesh aswell as the number of
 * partitions and the API mode from the opened HDF5 file \a file_id.
 *
 */
static int PMMG_loadHeader_hdf5(PMMG_pParMesh parmesh, hid_t file_id) {
  MMG5_pMesh mesh;
  hid_t attr_id;

  assert(parmesh->ngrp == 1);

  mesh = parmesh->listgrp[0].mesh;

  attr_id = H5Aopen(file_id, "MeshVersionFormatted", H5P_DEFAULT);
  H5Aread(attr_id, H5T_NATIVE_INT, &mesh->ver);
  H5Aclose(attr_id);

  attr_id = H5Aopen(file_id, "Dimension", H5P_DEFAULT);
  H5Aread(attr_id, H5T_NATIVE_INT, &mesh->dim);
  H5Aclose(attr_id);

  attr_id = H5Aopen(file_id, "NumberOfPartitions", H5P_DEFAULT);
  H5Aread(attr_id, H5T_NATIVE_INT, &parmesh->info.npartin);
  H5Aclose(attr_id);

  attr_id = H5Aopen(file_id, "API_mode", H5P_DEFAULT);
  H5Aread(attr_id, H5T_NATIVE_INT, &parmesh->info.API_mode);
  H5Aclose(attr_id);

  if (mesh->dim != 3) {
    if (parmesh->myrank == parmesh->info.root)
      fprintf(stderr,"\n  ## Error: %s: Wrong mesh dimension: %d (expected 3)!\n", __func__, mesh->dim);
    return 0;
  }

  if (parmesh->info.API_mode == PMMG_UNSET) {
    if (parmesh->myrank == parmesh->info.root)
      fprintf(stderr,"\n  ## Error: %s: No APIDISTRIB mode provided!\n", __func__);
    return 0;
  }

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param grp_part_id identifier of the HDF5 group from which to read the mesh partitioning.
 * \param dxpl_id identifier of the dataset transfer property list (MPI-IO).
 * \param nentities array of size nprocs * PMMG_NTYPENTITIES that will contain the number of entities on each proc.
 * \param nentitiesl array of size PMMG_NTYPENTITIES that will contain the local number of entities.
 * \param nentitiesg array of size PMMG_NTYPENTITIES that will contain the global number of entities.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Load the mesh partitioning and the communicators from the \a grp_part_id group
 * of the HDF5 file (only one group per process is allowed).
 * Three situations can occur:
 *   1/ nprocs = npartitions - Each proc just loads its corresponding partition
 *                             from the HDF5 file.
 *   2/ nprocs > npartitions - The excess procs do not load anything, and load balancing
 *                             is performed before the remeshing loop.
 *   3/ nprocs < npartitions - Some partitions are merged into the same proc, and load
 *                             balancing is performed before the remeshing loop.
 *
 * \warning Only the cases 1/ and 2/ are actually implemented yet.
 *
 */
static int PMMG_loadPartitioning_hdf5(PMMG_pParMesh parmesh, hid_t grp_part_id, hid_t dxpl_id,
                                      hsize_t *nentities, hsize_t *nentitiesl, hsize_t *nentitiesg) {
  hsize_t        *nentities_read;
  hsize_t        *ncomms, *nitem, *nitem_part;
  hsize_t        ncommg, comm_offset, nitemg, item_offset, rank_offset;
  hsize_t        icomm;
  int            *colors;
  int            **idx_loc, **idx_glob, *loc_buf, *glob_buf;
  int            npartitions, nprocs, rank;
  hid_t          dspace_file_id, dspace_mem_id;
  hid_t          dset_id;

  assert ( parmesh->ngrp == 1 );

  /* Set pointers to NULL */
  ncomms = nitem = nitem_part = NULL;
  colors = NULL;
  idx_loc = idx_glob = NULL;
  loc_buf = glob_buf = NULL;

  /* Init */
  nprocs = parmesh->nprocs;
  rank = parmesh->myrank;
  npartitions = parmesh->info.npartin;

  ncommg = nitemg = comm_offset = item_offset = rank_offset = 0;

  if (nprocs < npartitions) {
    fprintf(stderr, "\n ## Error : cannot read %d partitions with %d procs. \n",
            npartitions, nprocs);
    return 0;
  }

  /* Read the number of entities per partition */
  PMMG_CALLOC(parmesh, nentities_read, npartitions * PMMG_NTYPENTITIES, hsize_t, "nentities_read", goto free_and_return);
  dset_id = H5Dopen(grp_part_id, "NumberOfEntities", H5P_DEFAULT);
  H5Dread(dset_id, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL, dxpl_id, nentities_read);
  H5Dclose(dset_id);
  for (int i = 0 ; i < npartitions ; i++) {
    for (int j = 0 ; j < PMMG_NTYPENTITIES ; j++) {
      nentities[PMMG_NTYPENTITIES * i + j] = nentities_read[PMMG_NTYPENTITIES * i + j];
    }
  }
  PMMG_DEL_MEM(parmesh, nentities_read, hsize_t, "nentities_read");

  /* Set at least 1 communicator for each proc (even the ones that wont read the mesh) */
  if (parmesh->info.API_mode == PMMG_APIDISTRIB_faces)
    PMMG_Set_numberOfFaceCommunicators(parmesh, 1);
  else
    PMMG_Set_numberOfNodeCommunicators(parmesh, 1);

  if (rank < npartitions) {
    /* Get the local number of entities */
    for (int i = 0 ; i < PMMG_NTYPENTITIES ; i++)
      nentitiesl[i] = nentities[PMMG_NTYPENTITIES * rank + i];
  }

  /* Get the global number of entities */
  for (int i = 0 ; i < npartitions ; i++)
    for (int j = 0 ; j < PMMG_NTYPENTITIES ; j++)
      nentitiesg[j] += nentities[PMMG_NTYPENTITIES * i + j];

  /* Read the number of comms */
  PMMG_CALLOC(parmesh, ncomms, nprocs, hsize_t, "ncomms", goto free_and_return);
  if (parmesh->info.API_mode == PMMG_APIDISTRIB_faces)
    dset_id = H5Dopen(grp_part_id, "NumberOfFaceCommunicators", H5P_DEFAULT);
  else
    dset_id = H5Dopen(grp_part_id, "NumberOfNodeCommunicators", H5P_DEFAULT);
  H5Dread(dset_id, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ncomms);
  H5Dclose(dset_id);

  /* Compute the total number of comms */
  for (int i = 0 ; i < npartitions ; i++) {
    ncommg += ncomms[i];
  }
  for (int i = 0 ; i < rank ; i++) {
    comm_offset += ncomms[i];
  }

  /* Read the colors and the number of items */
  PMMG_MALLOC(parmesh, colors, ncomms[rank], int, "colors", goto free_and_return);
  PMMG_MALLOC(parmesh, nitem, ncomms[rank], hsize_t, "nitem", goto free_and_return);

  dset_id = H5Dopen(grp_part_id, "ColorsOut", H5P_DEFAULT);
  dspace_file_id = H5Dget_space(dset_id);
  dspace_mem_id = H5Screate_simple(1, &ncomms[rank], NULL);
  H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &comm_offset, NULL, &ncomms[rank], NULL);
  H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, colors);
  H5Dclose(dset_id);
  H5Sclose(dspace_file_id);
  H5Sclose(dspace_mem_id);

  if (parmesh->info.API_mode == PMMG_APIDISTRIB_faces)
    dset_id = H5Dopen(grp_part_id, "NumberOfCommunicatorFaces", H5P_DEFAULT);
  else
    dset_id = H5Dopen(grp_part_id, "NumberOfCommunicatorNodes", H5P_DEFAULT);
  dspace_file_id = H5Dget_space(dset_id);
  dspace_mem_id = H5Screate_simple(1, &ncomms[rank], NULL);
  H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &comm_offset, NULL, &ncomms[rank], NULL);
  H5Dread(dset_id, H5T_NATIVE_HSIZE, dspace_mem_id, dspace_file_id, dxpl_id, nitem);
  H5Dclose(dset_id);
  H5Sclose(dspace_file_id);
  H5Sclose(dspace_mem_id);

  /* Compute the total number of items and the item offset for the parallel reading */

  PMMG_CALLOC(parmesh, nitem_part, nprocs, hsize_t, "nitem_part", goto free_and_return);

  for (icomm = 0 ; icomm < ncomms[rank] ; icomm++) {
    nitem_part[rank] += nitem[icomm];
  }

  MPI_Allgather(&nitem_part[rank], 1, MPI_LONG_LONG, nitem_part, 1, MPI_LONG_LONG, parmesh->comm);

  for (int i = 0 ; i < npartitions ; i++) {
    nitemg += nitem_part[i];
  }

  for (int i = 0 ; i < rank ; i++) {
    rank_offset += nitem_part[i];
  }

  /* Read the communicator items in one buffer */

  PMMG_MALLOC(parmesh, loc_buf, nitem_part[rank], int, "loc_buf", goto free_and_return);
  PMMG_MALLOC(parmesh, glob_buf, nitem_part[rank], int, "glob_buf", goto free_and_return);

  dspace_file_id = H5Screate_simple(1, &nitemg, NULL);

  /* Local indices */
  if (parmesh->info.API_mode == PMMG_APIDISTRIB_faces)
    dset_id = H5Dopen(grp_part_id, "LocalFaceIndices", H5P_DEFAULT);
  else
    dset_id = H5Dopen(grp_part_id, "LocalNodeIndices", H5P_DEFAULT);
  dspace_mem_id = H5Screate_simple(1, &nitem_part[rank], NULL);
  H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &rank_offset, NULL, &nitem_part[rank], NULL);
  H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, loc_buf);
  H5Sclose(dspace_mem_id);
  H5Dclose(dset_id);

  /* Global indices */
  if (parmesh->info.API_mode == PMMG_APIDISTRIB_faces)
    dset_id = H5Dopen(grp_part_id, "GlobalFaceIndices", H5P_DEFAULT);
  else
    dset_id = H5Dopen(grp_part_id, "GlobalNodeIndices", H5P_DEFAULT);
  dspace_mem_id = H5Screate_simple(1, &nitem_part[rank], NULL);
  H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &rank_offset, NULL, &nitem_part[rank], NULL);
  H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, glob_buf);
  H5Sclose(dspace_mem_id);
  H5Dclose(dset_id);
  H5Sclose(dspace_file_id);

  H5Sclose(dspace_file_id);

  PMMG_DEL_MEM(parmesh, nitem_part, hsize_t, "nitem_part");

  /* Set the communicator items */
  PMMG_CALLOC(parmesh, idx_loc, ncomms[rank], int*, "idx_loc", goto free_and_return);
  PMMG_CALLOC(parmesh, idx_glob, ncomms[rank], int*, "idx_glob", goto free_and_return);
  for (icomm = 0 ; icomm < ncomms[rank] ; icomm++) {
    PMMG_CALLOC(parmesh, idx_loc[icomm], nitem[icomm], int, "idx_loc[icomm]", goto free_and_return);
    PMMG_CALLOC(parmesh, idx_glob[icomm], nitem[icomm], int, "idx_glob[icomm]", goto free_and_return);
  }

  for (icomm = 0 ; icomm < ncomms[rank] ; icomm++) {
    for (int k = 0 ; k < nitem[icomm] ; k++) {
      idx_loc[icomm][k] = loc_buf[item_offset + k];
      idx_glob[icomm][k] = glob_buf[item_offset + k];
    }
    item_offset += nitem[icomm];
  }

  /* Free the buffers */
  PMMG_DEL_MEM(parmesh, loc_buf, int, "loc_buf");
  PMMG_DEL_MEM(parmesh, glob_buf, int, "glob_buf");

  /* Set the communicators */
  if (rank < parmesh->info.npartin) {
    if (parmesh->info.API_mode == PMMG_APIDISTRIB_faces) {

      if ( !PMMG_Set_numberOfFaceCommunicators(parmesh, ncomms[rank]) ) {
        fprintf(stderr,"\n  ## Error: %s: unable to set number of face communicators on rank %d.\n",
                __func__, rank);
        goto free_and_return;
      }

      for (icomm = 0 ; icomm < ncomms[rank] ; icomm++) {

        if ( !PMMG_Set_ithFaceCommunicatorSize(parmesh, icomm, colors[icomm], nitem[icomm]) ) {
          fprintf(stderr,"\n  ## Error: %s: unable to set %lld th face communicator size on rank %d.\n",
                  __func__, icomm, rank);
          goto free_and_return;
        }

        if ( !PMMG_Set_ithFaceCommunicator_faces(parmesh, icomm, idx_loc[icomm], idx_glob[icomm], 1) ) {
          fprintf(stderr,"\n  ## Error: %s: unable to set %lld th face communicator faces on rank %d.\n",
                  __func__, icomm, rank);
          goto free_and_return;
        }

      }
    }
    else {

      if ( !PMMG_Set_numberOfNodeCommunicators(parmesh, ncomms[rank]) ) {
        fprintf(stderr,"\n  ## Error: %s: unable to set number of node communicators on rank %d.\n",
                __func__, rank);
        goto free_and_return;
      }

      for (icomm = 0 ; icomm < ncomms[rank] ; icomm++) {

        if ( !PMMG_Set_ithNodeCommunicatorSize(parmesh, icomm, colors[icomm], nitem[icomm]) ) {
          fprintf(stderr,"\n  ## Error: %s: unable to set %lld th node communicator size on rank %d.\n",
                  __func__, icomm, rank);
          goto free_and_return;
        }

        if ( !PMMG_Set_ithNodeCommunicator_nodes(parmesh, icomm, idx_loc[icomm], idx_glob[icomm], 1) ) {
          fprintf(stderr,"\n  ## Error: %s: unable to set %lld th node communicator faces on rank %d.\n",
                  __func__, icomm, rank);
          goto free_and_return;
        }

      }
    }
  }

  /* Free all memory */
  for (icomm = 0 ; icomm < ncomms[rank] ; icomm++) {
    PMMG_DEL_MEM(parmesh, idx_loc[icomm], int, "idx_loc[icomm]");
    PMMG_DEL_MEM(parmesh, idx_glob[icomm], int, "idx_glob[icomm]");
  }
  PMMG_DEL_MEM(parmesh, idx_loc, int*, "idx_loc");
  PMMG_DEL_MEM(parmesh, idx_glob, int*, "idx_glob");

  PMMG_DEL_MEM(parmesh, ncomms, hsize_t, "ncomms");
  PMMG_DEL_MEM(parmesh, colors, int, "colors");
  PMMG_DEL_MEM(parmesh, nitem, int, "nitem");

  return 1;

 free_and_return:
  PMMG_DEL_MEM(parmesh, loc_buf, int, "loc_buf");
  PMMG_DEL_MEM(parmesh, glob_buf, int, "glob_buf");
  PMMG_DEL_MEM(parmesh, ncomms, hsize_t, "ncomms");
  PMMG_DEL_MEM(parmesh, colors, int, "colors");
  PMMG_DEL_MEM(parmesh, nitem, int, "nitem");

  if ( idx_loc && idx_glob ) {
    for (icomm = 0 ; icomm < ncomms[rank] ; icomm++) {
      PMMG_DEL_MEM(parmesh, idx_loc[icomm], int, "idx_loc[icomms]");
      PMMG_DEL_MEM(parmesh, idx_glob[icomm], int, "idx_glob[icomm]");
    }
    PMMG_DEL_MEM(parmesh, idx_loc, int*, "idx_loc");
    PMMG_DEL_MEM(parmesh, idx_glob, int*, "idx_glob");
  }

  return 0;

}

static int PMMG_loadMeshEntities_hdf5(PMMG_pParMesh parmesh, hid_t grp_entities_id, hid_t dxpl_id, hsize_t *nentitiesl, hsize_t *nentitiesg, hsize_t *offset, int *load_entities) {
  /* MMG variables */
  PMMG_pGrp    grp;
  MMG5_pMesh   mesh;
  MMG5_pPoint  ppt;
  MMG5_pEdge   pa;
  MMG5_pTria   pt;
  MMG5_pQuad   pq;
  MMG5_pTetra  pe;
  MMG5_pPrism  pp;

  /* Local mesh size */
  hsize_t ne, np, nt, na, nquad, nprism;       /* Tetra, points, triangles, edges, quads, prisms */
  hsize_t nc, npreq, nppar;                    /* Corners, required and parallel vertices */
  hsize_t nr, nedreq, nedpar;                  /* Ridges, required and parallel edges */
  hsize_t ntreq, ntpar;                        /* Required and parallel triangles */
  hsize_t nqreq, nqpar;                        /* Required and parallel quads */
  hsize_t nereq, nepar;                        /* Required and parallel tetra */
  hsize_t nnor, ntan;                          /* Normals and Tangents */
  /* Global mesh size */
  hsize_t neg, npg, ntg, nag, nquadg, nprismg; /* Tetra, points, triangles, edges, quads, prisms */
  hsize_t ncg, npreqg, npparg;                 /* Corners, required and parallel vertices */
  hsize_t nrg, nedreqg, nedparg;               /* Ridges, required and parallel edges */
  hsize_t ntreqg, ntparg;                      /* Required and parallel triangles */
  hsize_t nqreqg, nqparg;                      /* Required and parallel quads */
  hsize_t nereqg, neparg;                      /* Required and parallel tetra */
  hsize_t nnorg, ntang;                        /* Normals and Tangents */

  /* Mesh buffer arrays */
  /* 6 buffers is the minimum amount for what we have to do */
  double *ppoint;      /* Point coordinates */
  int    *pent;        /* Other entities : edges, trias, quads, tetra, prisms. */
  int    *pcr;         /* Corners and ridges */
  int    *preq, *ppar; /* Required and parallel entities */
  int    *pref;        /* References */

  /* Normals and tangents */
  /* We could reuse the previous buffers, but the names would be confusing */
  int    *pnorat, *ptanat; /* Normals and Tangents at vertices */
  double *pnor, *ptan;     /* Normals and Tangents */

  /* Flag to remember if the save_entities array was NULL or not */
  int nullf = 0;

  /* HDF5 variables */
  hid_t dspace_mem_id, dspace_file_id;
  hid_t dset_id;

  /*------------------------- INIT -------------------------*/

  /* Set all buffers to NULL */
  ppoint = NULL;
  pent = NULL;
  pcr = NULL;
  preq = NULL; ppar = NULL;
  pref = NULL;
  pnor = NULL; ptan = NULL;
  pnorat = NULL; ptanat = NULL;

  /* Set ParMmg variables */
  grp = &parmesh->listgrp[0];
  mesh = grp->mesh;
  ppt = NULL;
  pa = NULL;
  pt = NULL;
  pq = NULL;
  pe = NULL;
  pp = NULL;

  /* Check the load_entities argument */
  if (load_entities == NULL) {
    nullf = 1;
    PMMG_CALLOC(parmesh, load_entities, PMMG_NTYPENTITIES, int, "load_entities", return 0);
    PMMG_Set_defaultIOEntities_hdf5(load_entities);
  }

  if (parmesh->myrank < parmesh->info.npartin) {

    /* Get the number of entities */
    np     = nentitiesl[PMMG_IO_Vertex];
    na     = nentitiesl[PMMG_IO_Edge];
    nt     = nentitiesl[PMMG_IO_Tria];
    nquad  = nentitiesl[PMMG_IO_Quad];
    ne     = nentitiesl[PMMG_IO_Tetra];
    nprism = nentitiesl[PMMG_IO_Prism];
    nc     = nentitiesl[PMMG_IO_Corner];
    npreq  = nentitiesl[PMMG_IO_Req];
    nppar  = nentitiesl[PMMG_IO_Par];
    nr     = nentitiesl[PMMG_IO_Ridge];
    nedreq = nentitiesl[PMMG_IO_EdReq];
    nedpar = nentitiesl[PMMG_IO_EdPar];
    ntreq  = nentitiesl[PMMG_IO_TriaReq];
    ntpar  = nentitiesl[PMMG_IO_TriaPar];
    nqreq  = nentitiesl[PMMG_IO_QuadReq];
    nqpar  = nentitiesl[PMMG_IO_QuadPar];
    nereq  = nentitiesl[PMMG_IO_TetReq];
    nepar  = nentitiesl[PMMG_IO_TetPar];
    nnor   = nentitiesl[PMMG_IO_Normal];
    ntan   = nentitiesl[PMMG_IO_Tangent];

    npg     = nentitiesg[PMMG_IO_Vertex];
    nag     = nentitiesg[PMMG_IO_Edge];
    ntg     = nentitiesg[PMMG_IO_Tria];
    nquadg  = nentitiesg[PMMG_IO_Quad];
    neg     = nentitiesg[PMMG_IO_Tetra];
    nprismg = nentitiesg[PMMG_IO_Prism];
    ncg     = nentitiesg[PMMG_IO_Corner];
    npreqg  = nentitiesg[PMMG_IO_Req];
    npparg  = nentitiesg[PMMG_IO_Par];
    nrg     = nentitiesg[PMMG_IO_Ridge];
    nedreqg = nentitiesg[PMMG_IO_EdReq];
    nedparg = nentitiesg[PMMG_IO_EdPar];
    ntreqg  = nentitiesg[PMMG_IO_TriaReq];
    ntparg  = nentitiesg[PMMG_IO_TriaPar];
    nqreqg  = nentitiesg[PMMG_IO_QuadReq];
    nqparg  = nentitiesg[PMMG_IO_QuadPar];
    nereqg  = nentitiesg[PMMG_IO_TetReq];
    neparg  = nentitiesg[PMMG_IO_TetPar];
    nnorg   = nentitiesg[PMMG_IO_Normal];
    ntang   = nentitiesg[PMMG_IO_Tangent];

    /* Arrays for bidimensional dataspaces */
    hsize_t hnp[2]      = {np, 3};
    hsize_t hna[2]      = {na, 2};
    hsize_t hnt[2]      = {nt, 3};
    hsize_t hnquad[2]   = {nquad, 4};
    hsize_t hne[2]      = {ne, 4};
    hsize_t hnprism[2]  = {nprism, 2};
    hsize_t hnnor[2]    = {nnor, 3};
    hsize_t hntan[2]    = {ntan, 3};
    hsize_t hnpg[2]     = {npg, 3};
    hsize_t hnag[2]     = {nag, 2};
    hsize_t hntg[2]     = {ntg, 3};
    hsize_t hnquadg[2]  = {nquadg, 4};
    hsize_t hneg[2]     = {neg, 4};
    hsize_t hnprismg[2] = {nprismg, 2};
    hsize_t hnnorg[2]   = {nnorg, 3};
    hsize_t hntang[2]   = {ntang, 3};

    PMMG_Set_meshSize(parmesh, np, ne, nprism, nt, nquad, na);

    /* Vertices, Normals and Tangents */
    if (load_entities[PMMG_IO_Vertex] && npg) {

      PMMG_MALLOC(parmesh, ppoint, 3 * np, double, "ppoint", goto free_and_return);
      PMMG_MALLOC(parmesh, pref, np, int, "pref", goto free_and_return);

      dspace_mem_id  = H5Screate_simple(2, hnp, NULL);
      dspace_file_id = H5Screate_simple(2, hnpg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Vertex], NULL, hnp, NULL);
      dset_id = H5Dopen(grp_entities_id, "Vertices", H5P_DEFAULT);
      H5Dread(dset_id, H5T_NATIVE_DOUBLE, dspace_mem_id, dspace_file_id, dxpl_id, ppoint);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);

      dspace_mem_id  = H5Screate_simple(1, hnp, NULL);
      dspace_file_id = H5Screate_simple(1, hnpg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Vertex], NULL, hnp, NULL);
      dset_id = H5Dopen(grp_entities_id, "VerticesRef", H5P_DEFAULT);
      H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);

      PMMG_Set_vertices(parmesh, ppoint, pref);

      PMMG_DEL_MEM(parmesh, ppoint, double, "ppoint");
      PMMG_DEL_MEM(parmesh, pref, int, "pref");

      if (load_entities[PMMG_IO_Corner] && ncg) {

        PMMG_MALLOC(parmesh, pcr, nc, int, "pcr", goto free_and_return);

        dspace_mem_id  = H5Screate_simple(1, &nc, NULL);
        dspace_file_id = H5Screate_simple(1, &ncg, NULL);
        H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Corner], NULL, &nc, NULL);
        dset_id = H5Dopen(grp_entities_id, "Corners", H5P_DEFAULT);
        H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pcr);
        H5Dclose(dset_id);
        H5Sclose(dspace_mem_id);
        H5Sclose(dspace_file_id);

        for (int i = 0 ; i < nc ; i++) {
          PMMG_Set_corner(parmesh, pcr[i] - offset[2 * PMMG_IO_Vertex] + 1);
        }

        PMMG_DEL_MEM(parmesh, pcr, int, "pcr");

      }

      if (load_entities[PMMG_IO_Req] && npreqg) {

        PMMG_MALLOC(parmesh, preq, npreq, int, "preq", goto free_and_return);

        dspace_mem_id  = H5Screate_simple(1, &npreq, NULL);
        dspace_file_id = H5Screate_simple(1, &npreqg, NULL);
        H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Req], NULL, &npreq, NULL);
        dset_id = H5Dopen(grp_entities_id, "RequiredVertices", H5P_DEFAULT);
        H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
        H5Dclose(dset_id);
        H5Sclose(dspace_mem_id);
        H5Sclose(dspace_file_id);

        for (int i = 0 ; i < npreq ; i++) {
          PMMG_Set_requiredVertex(parmesh, preq[i] - offset[2 * PMMG_IO_Vertex] + 1);
        }

        PMMG_DEL_MEM(parmesh, preq, int, "preq");

      }

      if (load_entities[PMMG_IO_Par] && npparg) {

        PMMG_MALLOC(parmesh, ppar, nppar, int, "ppar", goto free_and_return);

        dspace_mem_id  = H5Screate_simple(1, &nppar, NULL);
        dspace_file_id = H5Screate_simple(1, &npparg, NULL);
        H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Par], NULL, &nppar, NULL);
        dset_id = H5Dopen(grp_entities_id, "ParallelVertices", H5P_DEFAULT);
        H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
        H5Dclose(dset_id);
        H5Sclose(dspace_mem_id);
        H5Sclose(dspace_file_id);

        for (int i = 0 ; i < nppar ; i++) {
          ppt = &mesh->point[ppar[i] - offset[2 * PMMG_IO_Vertex] + 1];
          ppt->tag |= MG_PARBDY;
        }

        PMMG_DEL_MEM(parmesh, ppar, int, "ppar");

      }

      if (load_entities[PMMG_IO_Normal] && nnorg) {

        PMMG_MALLOC(parmesh, pnor, 3 * nnor, double, "pnor", goto free_and_return);
        PMMG_MALLOC(parmesh, pnorat, nnor, int, "pnorat", goto free_and_return);

        dspace_mem_id  = H5Screate_simple(2, hnnor, NULL);
        dspace_file_id = H5Screate_simple(2, hnnorg, NULL);
        H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Normal], NULL, hnnor, NULL);
        dset_id = H5Dopen(grp_entities_id, "Normals", H5P_DEFAULT);
        H5Dread(dset_id, H5T_NATIVE_DOUBLE, dspace_mem_id, dspace_file_id, dxpl_id, pnor);
        H5Dclose(dset_id);
        H5Sclose(dspace_mem_id);
        H5Sclose(dspace_file_id);

        dspace_mem_id  = H5Screate_simple(1, hnnor, NULL);
        dspace_file_id = H5Screate_simple(1, hnnorg, NULL);
        H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Normal], NULL, hnnor, NULL);
        dset_id = H5Dopen(grp_entities_id, "NormalsAtVertices", H5P_DEFAULT);
        H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pnorat);
        H5Dclose(dset_id);
        H5Sclose(dspace_mem_id);
        H5Sclose(dspace_file_id);

        for (int i = 0 ; i < nnor ; i++) {
          PMMG_Set_normalAtVertex(parmesh, pnorat[i] - offset[2 * PMMG_IO_Vertex] + 1, pnor[3 * i], pnor[3 * i + 1], pnor[3 * i + 2]);
        }

        PMMG_DEL_MEM(parmesh, pnor, double, "pnor");
        PMMG_DEL_MEM(parmesh, pnorat, int, "pnorat");

      }

      if (load_entities[PMMG_IO_Tangent] && ntang) {

        PMMG_MALLOC(parmesh, ptan, 3 * ntan, double, "ptan", goto free_and_return);
        PMMG_MALLOC(parmesh, ptanat, ntan, int, "ptanat", goto free_and_return);

        dspace_mem_id  = H5Screate_simple(2, hntan, NULL);
        dspace_file_id = H5Screate_simple(2, hntang, NULL);
        H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Tangent], NULL, hntan, NULL);
        dset_id = H5Dopen(grp_entities_id, "Tangents", H5P_DEFAULT);
        H5Dread(dset_id, H5T_NATIVE_DOUBLE, dspace_mem_id, dspace_file_id, dxpl_id, ptan);
        H5Dclose(dset_id);
        H5Sclose(dspace_mem_id);
        H5Sclose(dspace_file_id);

        dspace_mem_id  = H5Screate_simple(1, hntan, NULL);
        dspace_file_id = H5Screate_simple(1, hntang, NULL);
        H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Tangent], NULL, hntan, NULL);
        dset_id = H5Dopen(grp_entities_id, "TangentsAtVertices", H5P_DEFAULT);
        H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ptanat);
        H5Dclose(dset_id);
        H5Sclose(dspace_mem_id);
        H5Sclose(dspace_file_id);

        PMMG_DEL_MEM(parmesh, ptan, double, "ptan");
        PMMG_DEL_MEM(parmesh, ptanat, int, "ptanat");

      }

    }

    /* Edges */
    if (load_entities[PMMG_IO_Edge] && nag) {

      PMMG_MALLOC(parmesh, pent, 2 * na, int, "pent", goto free_and_return);
      PMMG_MALLOC(parmesh, pref, na, int, "pref", goto free_and_return);

      dspace_mem_id  = H5Screate_simple(2, hna, NULL);
      dspace_file_id = H5Screate_simple(2, hnag, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Edge], NULL, hna, NULL);
      dset_id = H5Dopen(grp_entities_id, "Edges", H5P_DEFAULT);
      H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);

      dspace_mem_id  = H5Screate_simple(1, hna, NULL);
      dspace_file_id = H5Screate_simple(1, hnag, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Edge], NULL, hna, NULL);
      dset_id = H5Dopen(grp_entities_id, "EdgesRef", H5P_DEFAULT);
      H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);

      for (int i = 0 ; i < na ; i++) {
        PMMG_Set_edge(parmesh, pent[2 * i] - offset[PMMG_IO_Vertex] + 1, pent[2 * i + 1] - offset[PMMG_IO_Vertex] + 1, pref[i], i + 1);
      }

      PMMG_DEL_MEM(parmesh, pent, int, "pent");
      PMMG_DEL_MEM(parmesh, pref, int, "pref");

      if (load_entities[PMMG_IO_Ridge] && nrg) {

        PMMG_MALLOC(parmesh, pcr, nr, int, "pcr", goto free_and_return);

        dspace_mem_id  = H5Screate_simple(1, &nr, NULL);
        dspace_file_id = H5Screate_simple(1, &nrg, NULL);
        H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Ridge], NULL, &nr, NULL);
        dset_id = H5Dopen(grp_entities_id, "Ridges", H5P_DEFAULT);
        H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pcr);
        H5Dclose(dset_id);
        H5Sclose(dspace_mem_id);
        H5Sclose(dspace_file_id);

        for (int i = 0 ; i < nr ; i++) {
          PMMG_Set_ridge(parmesh, pcr[i] - offset[2 * PMMG_IO_Edge] + 1);
        }

        PMMG_DEL_MEM(parmesh, pcr, int, "pcr");

      }

      if (load_entities[PMMG_IO_EdReq] && nedreqg) {

        PMMG_MALLOC(parmesh, preq, nedreq, int, "preq", goto free_and_return);

        dspace_mem_id  = H5Screate_simple(1, &nedreq, NULL);
        dspace_file_id = H5Screate_simple(1, &nedreqg, NULL);
        H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_EdReq], NULL, &nedreq, NULL);
        dset_id = H5Dopen(grp_entities_id, "RequiredEdges", H5P_DEFAULT);
        H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
        H5Dclose(dset_id);
        H5Sclose(dspace_mem_id);
        H5Sclose(dspace_file_id);

        for (int i = 0 ; i < nedreq ; i++) {
          PMMG_Set_requiredEdge(parmesh, preq[i] - offset[2 * PMMG_IO_Edge] + 1);
        }

        PMMG_DEL_MEM(parmesh, preq, int, "preq");

      }

      if (load_entities[PMMG_IO_EdPar] && nedparg) {

        PMMG_MALLOC(parmesh, ppar, nedpar, int, "ppar", goto free_and_return);

        dspace_mem_id  = H5Screate_simple(1, &nedpar, NULL);
        dspace_file_id = H5Screate_simple(1, &nedparg, NULL);
        H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_EdPar], NULL, &nedpar, NULL);
        dset_id = H5Dopen(grp_entities_id, "ParallelEdges", H5P_DEFAULT);
        H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
        H5Dclose(dset_id);
        H5Sclose(dspace_mem_id);
        H5Sclose(dspace_file_id);

        for (int i = 0 ; i < nedpar ; i++) {
          pa = &mesh->edge[ppar[i] - offset[2 * PMMG_IO_Edge] + 1];
          pa->tag |= MG_PARBDY;
        }

        PMMG_DEL_MEM(parmesh, ppar, int, "ppar");

      }
    }

    /* Triangles */
    if (load_entities[PMMG_IO_Tria] && ntg) {

      PMMG_MALLOC(parmesh, pent, 3 * nt, int, "pent", goto free_and_return);
      PMMG_MALLOC(parmesh, pref, nt, int, "pref", goto free_and_return);

      dspace_mem_id  = H5Screate_simple(2, hnt, NULL);
      dspace_file_id = H5Screate_simple(2, hntg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Tria], NULL, hnt, NULL);
      dset_id = H5Dopen(grp_entities_id, "Triangles", H5P_DEFAULT);
      H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);

      dspace_mem_id  = H5Screate_simple(1, hnt, NULL);
      dspace_file_id = H5Screate_simple(1, hntg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Tria], NULL, hnt, NULL);
      dset_id = H5Dopen(grp_entities_id, "TrianglesRef", H5P_DEFAULT);
      H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);

      for (int i = 0 ; i < nt ; i++) {
        PMMG_Set_triangle(parmesh, pent[3 * i] - offset[2 * PMMG_IO_Vertex] + 1,
                          pent[3 * i + 1] - offset[2 * PMMG_IO_Vertex] + 1,
                          pent[3 * i + 2] - offset[2 * PMMG_IO_Vertex] + 1,
                          pref[i], i + 1);
      }

      PMMG_DEL_MEM(parmesh, pent, int, "pent");
      PMMG_DEL_MEM(parmesh, pref, int, "pref");

      if (load_entities[PMMG_IO_TriaReq] && ntreqg) {

        PMMG_MALLOC(parmesh, preq, ntreq, int, "preq", goto free_and_return);

        dspace_mem_id  = H5Screate_simple(1, &ntreq, NULL);
        dspace_file_id = H5Screate_simple(1, &ntreqg, NULL);
        H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_TriaReq], NULL, &ntreq, NULL);
        dset_id = H5Dopen(grp_entities_id, "RequiredTriangles", H5P_DEFAULT);
        H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
        H5Dclose(dset_id);
        H5Sclose(dspace_mem_id);
        H5Sclose(dspace_file_id);

        for (int i = 0 ; i < ntreq ; i++) {
          PMMG_Set_requiredTriangle(parmesh, preq[i] + offset[2 * PMMG_IO_Tria] + 1);
        }

        PMMG_DEL_MEM(parmesh, preq, int, "preq");

      }

      if (load_entities[PMMG_IO_TriaPar] && ntparg) {

        PMMG_MALLOC(parmesh, ppar, ntpar, int, "ppar", goto free_and_return);

        dspace_mem_id  = H5Screate_simple(1, &ntpar, NULL);
        dspace_file_id = H5Screate_simple(1, &ntparg, NULL);
        H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_TriaPar], NULL, &ntpar, NULL);
        dset_id = H5Dopen(grp_entities_id, "ParallelTriangles", H5P_DEFAULT);
        H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
        H5Dclose(dset_id);
        H5Sclose(dspace_mem_id);
        H5Sclose(dspace_file_id);

        for (int i = 0 ; i < ntpar ; i++) {
          pt = &mesh->tria[ppar[i] - offset[2 * PMMG_IO_Tria] + 1];
          for (int j = 0 ; j < 3 ; j++) {
            pt->tag[j] |= MG_PARBDY;
          }
        }

        PMMG_DEL_MEM(parmesh, ppar, int, "ppar");

      }
    }


    /* Quadrilaterals */
    if (load_entities[PMMG_IO_Quad] && nquadg) {

      PMMG_MALLOC(parmesh, pent, 4 * nquad, int, "pent", goto free_and_return);
      PMMG_MALLOC(parmesh, pref, nquad, int, "pref", goto free_and_return);

      dspace_mem_id  = H5Screate_simple(2, hnquad, NULL);
      dspace_file_id = H5Screate_simple(2, hnquadg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Quad], NULL, hnquad, NULL);
      dset_id = H5Dopen(grp_entities_id, "Quadrilaterals", H5P_DEFAULT);
      H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);

      dspace_mem_id  = H5Screate_simple(1, hnquad, NULL);
      dspace_file_id = H5Screate_simple(1, hnquadg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Quad], NULL, hnquad, NULL);
      dset_id = H5Dopen(grp_entities_id, "QuadrilateralsRef", H5P_DEFAULT);
      H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);

      for (int i = 0 ; i < nquad ; i++) {
        PMMG_Set_quadrilateral(parmesh, pent[4 * i] - offset[2 * PMMG_IO_Vertex] + 1,
                               pent[4 * i + 1] - offset[2 * PMMG_IO_Vertex] + 1,
                               pent[4 * i + 2] - offset[2 * PMMG_IO_Vertex] + 1,
                               pent[4 * i + 3] - offset[2 * PMMG_IO_Vertex] + 1,
                               pref[i], i + 1);
      }

      PMMG_DEL_MEM(parmesh, pent, int, "pent");
      PMMG_DEL_MEM(parmesh, pref, int, "pref");

      if (load_entities[PMMG_IO_QuadReq] && nqreqg) {

        PMMG_MALLOC(parmesh, preq, nqreq, int, "preq", goto free_and_return);

        dspace_mem_id  = H5Screate_simple(1, &nqreq, NULL);
        dspace_file_id = H5Screate_simple(1, &nqreqg, NULL);
        H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_QuadReq], NULL, &nqreq, NULL);
        dset_id = H5Dopen(grp_entities_id, "RequiredQuadrilaterals", H5P_DEFAULT);
        H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
        H5Dclose(dset_id);
        H5Sclose(dspace_mem_id);
        H5Sclose(dspace_file_id);

        for (int i = 0 ; i < nqreq ; i++) {
          pq = &mesh->quadra[preq[i] - offset[2 * PMMG_IO_Quad] + 1];
          for (int j = 0 ; j < 4 ; j++) {
            pq->tag[j] |= MG_REQ;
          }
        }

        PMMG_DEL_MEM(parmesh, preq, int, "preq");

      }

      if (load_entities[PMMG_IO_QuadPar] && nqparg) {

        PMMG_MALLOC(parmesh, ppar, nqpar, int, "ppar", goto free_and_return);

        dspace_mem_id  = H5Screate_simple(1, &nqpar, NULL);
        dspace_file_id = H5Screate_simple(1, &nqparg, NULL);
        H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_QuadPar], NULL, &nqpar, NULL);
        dset_id = H5Dopen(grp_entities_id, "ParallelQuadrilaterals", H5P_DEFAULT);
        H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
        H5Dclose(dset_id);
        H5Sclose(dspace_mem_id);
        H5Sclose(dspace_file_id);

        for (int i = 0 ; i < nqpar ; i++) {
          pq = &mesh->quadra[ppar[i] - offset[2 * PMMG_IO_Quad] + 1];
          for (int j = 0 ; j < 4 ; j++) {
            pq->tag[j] |= MG_PARBDY;
          }
        }

        PMMG_DEL_MEM(parmesh, ppar, int, "ppar");

      }
    }

    /* Tetrahedra */
    if (load_entities[PMMG_IO_Tetra] && neg) {

      PMMG_MALLOC(parmesh, pent, 4 * ne, int, "pent", goto free_and_return);
      PMMG_MALLOC(parmesh, pref, ne, int, "pref", goto free_and_return);

      dspace_mem_id  = H5Screate_simple(2, hne, NULL);
      dspace_file_id = H5Screate_simple(2, hneg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Tetra], NULL, hne, NULL);
      dset_id = H5Dopen(grp_entities_id, "Tetrahedra", H5P_DEFAULT);
      H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);

      dspace_mem_id  = H5Screate_simple(1, hne, NULL);
      dspace_file_id = H5Screate_simple(1, hneg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Tetra], NULL, hne, NULL);
      dset_id = H5Dopen(grp_entities_id, "TetrahedraRef", H5P_DEFAULT);
      H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);

      for (int i = 0 ; i < ne ; i++) {
        PMMG_Set_tetrahedron(parmesh, pent[4 * i] - offset[2 * PMMG_IO_Vertex] + 1,
                             pent[4 * i + 1] - offset[2 * PMMG_IO_Vertex] + 1,
                             pent[4 * i + 2] - offset[2 * PMMG_IO_Vertex] + 1,
                             pent[4 * i + 3] - offset[2 * PMMG_IO_Vertex] + 1,
                             pref[i], i + 1);
      }

      PMMG_DEL_MEM(parmesh, pent, int, "pent");
      PMMG_DEL_MEM(parmesh, pref, int, "pref");

      if (load_entities[PMMG_IO_TetReq] && nereqg) {

        PMMG_MALLOC(parmesh, preq, nereq, int, "preq", goto free_and_return);

        dspace_mem_id  = H5Screate_simple(1, &nereq, NULL);
        dspace_file_id = H5Screate_simple(1, &nereqg, NULL);
        H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_TetReq], NULL, &nereq, NULL);
        dset_id = H5Dopen(grp_entities_id, "RequiredTetrahedra", H5P_DEFAULT);
        H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
        H5Dclose(dset_id);
        H5Sclose(dspace_mem_id);
        H5Sclose(dspace_file_id);

        for (int i = 0 ; i < nereq ; i++) {
          PMMG_Set_requiredTetrahedron(parmesh, preq[i] - offset[2 * PMMG_IO_Tetra] + 1);
        }

        PMMG_DEL_MEM(parmesh, preq, int, "preq");

      }

      if (load_entities[PMMG_IO_TetPar] && neparg) {

        PMMG_MALLOC(parmesh, ppar, nepar, int, "ppar", goto free_and_return);

        dspace_mem_id  = H5Screate_simple(1, &nepar, NULL);
        dspace_file_id = H5Screate_simple(1, &neparg, NULL);
        H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_TetPar], NULL, &nepar, NULL);
        dset_id = H5Dopen(grp_entities_id, "ParallelTetrahedra", H5P_DEFAULT);
        H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
        H5Dclose(dset_id);
        H5Sclose(dspace_mem_id);
        H5Sclose(dspace_file_id);

        for (int i = 0 ; i < nepar ; i++) {
          pe = &mesh->tetra[ppar[i] - offset[2 * PMMG_IO_Tetra] + 1];
          pe->tag |= MG_PARBDY;
        }

        PMMG_DEL_MEM(parmesh, ppar, int, "ppar");

      }
    }

    /* Prisms */
    if (load_entities[PMMG_IO_Prism] && nprismg) {
      PMMG_MALLOC(parmesh, pent, 6 * nprism, int, "pent", goto free_and_return);
      PMMG_MALLOC(parmesh, pref, nprism, int, "pref", goto free_and_return);

      dspace_mem_id  = H5Screate_simple(2, hnprism, NULL);
      dspace_file_id = H5Screate_simple(2, hnprismg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Prism], NULL, hnprism, NULL);
      dset_id = H5Dopen(grp_entities_id, "Prisms", H5P_DEFAULT);
      H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);

      dspace_mem_id  = H5Screate_simple(1, hnprism, NULL);
      dspace_file_id = H5Screate_simple(1, hnprismg, NULL);
      H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_IO_Prism], NULL, hnprism, NULL);
      dset_id = H5Dopen(grp_entities_id, "PrismsRef", H5P_DEFAULT);
      H5Dread(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
      H5Dclose(dset_id);
      H5Sclose(dspace_mem_id);
      H5Sclose(dspace_file_id);

      for (int i = 0 ; i < nprism ; i++) {
        PMMG_Set_prism(parmesh, pent[6 * i] - offset[2 * PMMG_IO_Vertex] + 1,
                       pent[6 * i + 1] - offset[2 * PMMG_IO_Vertex] + 1,
                       pent[6 * i + 2] - offset[2 * PMMG_IO_Vertex] + 1,
                       pent[6 * i + 3] - offset[2 * PMMG_IO_Vertex] + 1,
                       pent[6 * i + 4] - offset[2 * PMMG_IO_Vertex] + 1,
                       pent[6 * i + 5] - offset[2 * PMMG_IO_Vertex] + 1,
                       pref[i], i + 1);
      }

      PMMG_DEL_MEM(parmesh, pent, int, "pent");
      PMMG_DEL_MEM(parmesh, pref, int, "pref");

    }

    /* Print the number of entities */
    if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
      fprintf(stdout,"     NUMBER OF VERTICES       %lld\n", npg);
      fprintf(stdout,"     NUMBER OF TETRAHEDRA     %lld\n", neg);
      if ( nprismg )
        fprintf(stdout,"     NUMBER OF PRISMS         %lld\n",nprismg);

      if ( nag ) {
        fprintf(stdout,"     NUMBER OF EDGES          %lld\n",nag);
        if ( nrg )
          fprintf(stdout,"     NUMBER OF RIDGES         %lld\n",nrg);
      }
      if ( ntg )
        fprintf(stdout,"     NUMBER OF TRIANGLES      %lld\n",ntg);
      if ( nquadg )
        fprintf(stdout,"     NUMBER OF QUADRILATERALS %lld\n",nquadg);


      if ( npreqg || nedreqg || ntreqg || nereqg || nqreqg ) {
        fprintf(stdout,"     NUMBER OF REQUIRED ENTITIES: \n");
        if ( npreqg )
          fprintf(stdout,"                  VERTICES       %lld \n",npreqg);
        if ( nedreqg )
          fprintf(stdout,"                  EDGES          %lld \n",nedreqg);
        if ( ntreqg )
          fprintf(stdout,"                  TRIANGLES      %lld \n",ntreqg);
        if ( nqreqg )
          fprintf(stdout,"                  QUADRILATERALS %lld \n",nqreqg);
        if ( nereqg )
          fprintf(stdout,"                  TETRAHEDRA    %lld \n",nereqg);
      }
      if( ncg )
        fprintf(stdout,"     NUMBER OF CORNERS        %lld \n",ncg);

      if ( npparg || nedparg || ntparg || neparg || nqparg ) {
        fprintf(stdout,"     NUMBER OF PARALLEL ENTITIES: \n");
        if ( npparg )
          fprintf(stdout,"                  VERTICES       %lld \n",npparg);
        if ( nedparg )
          fprintf(stdout,"                  EDGES          %lld \n",nedparg);
        if ( neparg )
          fprintf(stdout,"                  TRIANGLES      %lld \n",ntparg);
        if ( nqparg )
          fprintf(stdout,"                  QUADRILATERALS %lld \n",nqparg);
        if ( neparg )
          fprintf(stdout,"                  TETRAHEDRA     %lld \n",neparg);
      }
    }

  }

  if (nullf) PMMG_DEL_MEM(parmesh, load_entities, int, "load_entities");

  return 1;

 free_and_return:
  PMMG_DEL_MEM(parmesh, ppoint, double, "ppoint");
  PMMG_DEL_MEM(parmesh, pent, int, "pent");
  PMMG_DEL_MEM(parmesh, pref, int, "pref");
  PMMG_DEL_MEM(parmesh, pcr, int, "pcr");
  PMMG_DEL_MEM(parmesh, preq, int, "preq");
  PMMG_DEL_MEM(parmesh, ppar, int, "ppar");
  PMMG_DEL_MEM(parmesh, pnor, double, "pnor");
  PMMG_DEL_MEM(parmesh, pnorat, int, "pnorat");
  PMMG_DEL_MEM(parmesh, ptan, double, "ptan");
  PMMG_DEL_MEM(parmesh, ptanat, int, "ptanat");
  if (nullf) PMMG_DEL_MEM(parmesh, load_entities, int, "save_entities");

  return 0;

}

static int PMMG_loadMetric_hdf5(PMMG_pParMesh parmesh, hid_t grp_sols_id, hid_t dxpl_id,
                                hsize_t *nentitiesl, hsize_t *offset) {
  int np;
  MMG5_pMesh mesh;
  MMG5_pSol met;
  MMG5_pPoint ppt;
  double *met_buf;
  hsize_t met_offset[2] = {0, 0};
  hsize_t hnsg[2] = {0, 0};
  hid_t dspace_mem_id, dspace_file_id;
  hid_t dset_id;

  assert ( parmesh->ngrp == 1 );

  /* Init mmg variables*/
  mesh = parmesh->listgrp[0].mesh;
  met = parmesh->listgrp[0].met;
  np = nentitiesl[PMMG_IO_Vertex];

  /* Get the metric size */
  dset_id = H5Dopen(grp_sols_id, "MetricAtVertices", H5P_DEFAULT);
  if (dset_id < 0) return -1;
  dspace_file_id = H5Dget_space(dset_id);
  H5Sget_simple_extent_dims(dspace_file_id, hnsg, NULL);
  H5Sclose(dspace_file_id);
  H5Dclose(dset_id);

  /* Set the metric size */
  if (hnsg[1] == 1)
    PMMG_Set_metSize(parmesh, MMG5_Vertex, np, MMG5_Scalar);
  else if (hnsg[1] == 3 && mesh->info.lag != -1)
    PMMG_Set_metSize(parmesh, MMG5_Vertex, np, MMG5_Vector);
  else if (hnsg[1] == 6)
    PMMG_Set_metSize(parmesh, MMG5_Vertex, np, MMG5_Tensor);
  else {
    if (parmesh->myrank == parmesh->info.root) {
      fprintf(stderr, "\n  ## Error: %s: Wrong metric size/type \n", __func__);
    }
    return 0;
  }

  /* Compute the offset for parallel reading */
  hsize_t hns[2] = {np, met->size};
  met_offset[0] = offset[2 * PMMG_IO_Vertex];

  /* Read the metric buffer */
  PMMG_MALLOC(parmesh, met_buf, np * met->size, double, "met_buf", return 0);
  dset_id = H5Dopen(grp_sols_id, "MetricAtVertices", H5P_DEFAULT);
  dspace_file_id = H5Dget_space(dset_id);
  dspace_mem_id = H5Screate_simple(2, hns, NULL);
  H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, met_offset, NULL, hns, NULL);
  H5Dread(dset_id, H5T_NATIVE_DOUBLE, dspace_mem_id, dspace_file_id, dxpl_id, met_buf);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  H5Dclose(dset_id);

  /* Set the metric */
  for (int k = 0 ; k < mesh->np ; k++) {
    ppt = &mesh->point[k + 1];
    if (!MG_VOK(ppt)) continue;
    for (int j = 0 ; j < met->size ; j++) {
      met->m[1 + k * met->size + j] = met_buf[k * met->size + j];
    }
  }

  return 1;
}

static int  PMMG_loadAllSols_hdf5(PMMG_pParMesh parmesh, hid_t grp_sols_id, hid_t dcpl_id, hsize_t *nentitiesl, hsize_t *nentitiesg, hsize_t *offset) {
  return 1;
}

int PMMG_loadParmesh_hdf5(PMMG_pParMesh parmesh, int *load_entities, const char *filename) {

#ifndef USE_HDF5

  fprintf(stderr,"  ** HDF5 library not found. Unavailable file format.\n");
  return -1;

#else

  int ier = 1;
  hsize_t *nentities, *nentitiesl, *nentitiesg;
  hsize_t *offset;
  int npartitions;
  hid_t file_id, grp_mesh_id, grp_part_id, grp_entities_id, grp_sols_id; /* Objects */
  hid_t fapl_id, dxpl_id;                                                /* Property lists */
  MPI_Info info = MPI_INFO_NULL;
  MPI_Comm read_comm;
  int rank, nprocs, mpi_color;
  mytime ctim[TIMEMAX];
  int8_t tim;
  char   stim[32];

  /* Check arguments */
  if (parmesh->ngrp != 1) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in your parmesh.",
            __func__);
    return 0;
  }
  if (!filename || !*filename) {
    fprintf(stderr,"  ## Error: %s: no HDF5 file name provided.",
            __func__);
    return 0;
  }
  if (!load_entities[PMMG_IO_Vertex] || !load_entities[PMMG_IO_Tetra]) {
    fprintf(stderr, "\n  ## Error: %s: load_entities: you must at least load the vertices and the tetra.\n",
            __func__);
    return 0;
  }

  tminit(ctim, TIMEMAX);
  chrono(ON, &ctim[0]);

  /* Set all pointers to NULL */
  nentities = nentitiesl = nentitiesg = offset = NULL;

  /* Set MPI variables */
  nprocs = parmesh->nprocs;
  rank = parmesh->myrank;

  /* Store the input format in the parmesh->info.fmtout field */
  parmesh->info.fmtout = PMMG_FMT_HDF5;

  /* Shut HDF5 error stack */
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);

  /** Open the file a first time to read npartin */

  /* Create the property lists */
  fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl_id, parmesh->comm, info);  /* Parallel access to the file */
  H5Pset_all_coll_metadata_ops(fapl_id, 1);        /* Collective metadata read */
  dxpl_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE); /* Collective dataset xfer operations */

  /* Open the HDF5 file */
  HDF_CHECK( file_id = H5Fopen(filename, H5F_ACC_RDONLY, fapl_id),
             fprintf(stderr,"\n  ## Error: %s: Rank %d could not open the hdf5 file %s.\n",
                     __func__, rank, filename);
             goto free_and_return );

  /* Load the header (version, dimension, number of partitions and API mode) */
  ier = PMMG_loadHeader_hdf5(parmesh, file_id);

  MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, &ier, 1, MPI_INT, MPI_MIN, parmesh->comm),
             goto free_and_return );

  if (ier == 0) {
    if (rank == parmesh->info.root) {
      fprintf(stderr,"\n  ## Error: %s: Wrong mesh attributes in hdf5 file %s.\n",
              __func__, filename);
    }
    goto free_and_return;
  }

  /* Close the file and create a new communicator if there are less partitions
     than MPI processes */
  H5Fclose(file_id);
  H5Pclose(fapl_id);

  npartitions = parmesh->info.npartin;

  /* Reading more partitions than there are procs available is not supported yet */
  if (npartitions > nprocs) {
    if (rank == parmesh->info.root) {
      fprintf(stderr,"\n  ## Error: %s: Can't read %d partitions with %d procs yet.\n",
              __func__, npartitions, nprocs);
    }
    return 0;
  }

  /* Set the new communicator containing the procs reading the mesh */
  mpi_color = (rank < npartitions) ? 1 : 0;

  MPI_CHECK( MPI_Comm_split(parmesh->comm, mpi_color, rank, &read_comm),
             goto free_and_return );

  /* Set MPI error handling */
  MPI_CHECK( MPI_Comm_set_errhandler(read_comm, MPI_ERRORS_RETURN),
             goto free_and_return );

  parmesh->info.read_comm = read_comm;

  /** Open the file a second (and final) time to actually read the mesh */

  /* Set the file access property list with the new communicator */
  fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl_id, read_comm, info);    /* Parallel access to the file */
  H5Pset_all_coll_metadata_ops(fapl_id, 1);      /* Collective metadata read */

  /* Reopen the file with the new communicator */
  HDF_CHECK( file_id = H5Fopen(filename, H5F_ACC_RDONLY, fapl_id),
             fprintf(stderr,"\n  ## Error: %s: Rank %d could not open the hdf5 file %s.\n",
                     __func__, rank, filename);
             goto free_and_return );

  if (parmesh->info.imprim > PMMG_VERB_VERSION) {
    fprintf(stdout, "\n  %%%% %s OPENED \n", filename);
  }

  /* Open the mesh group */
  HDF_CHECK( grp_mesh_id = H5Gopen(file_id, "Mesh", H5P_DEFAULT),
             fprintf(stderr,"\n  ## Error: %s: Rank %d could not open the /Mesh group in file %s.\n",
                     __func__, rank, filename);
             goto free_and_return );

  /** Load the partitioning information*/
  tim = 1;
  chrono(ON, &ctim[tim]);

  /* Open the partitionning group */
  HDF_CHECK( grp_part_id = H5Gopen(grp_mesh_id, "Partitioning", H5P_DEFAULT),
             fprintf(stderr,"\n  ## Error: %s: Rank %d could not open the /Mesh/Partitioning group in file %s.\n",
                     __func__, rank, filename);
             H5Gclose(grp_mesh_id);
             goto free_and_return );

  /* Load the old partitioning of the mesh */
  PMMG_CALLOC(parmesh, nentities, PMMG_NTYPENTITIES * nprocs, hsize_t, "nentities",
              goto free_and_return );
  PMMG_CALLOC(parmesh, nentitiesl, PMMG_NTYPENTITIES, hsize_t, "nentitiesl",
              goto free_and_return );
  PMMG_CALLOC(parmesh, nentitiesg, PMMG_NTYPENTITIES, hsize_t, "nentitiesg",
              goto free_and_return );

  ier = PMMG_loadPartitioning_hdf5(parmesh, grp_part_id, dxpl_id, nentities, nentitiesl, nentitiesg);

  MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, &ier, 1, MPI_INT, MPI_MIN, read_comm),
             H5Gclose(grp_part_id);
             H5Gclose(grp_mesh_id);
             goto free_and_return );

  if (ier == 0) {
    if (rank == parmesh->info.root) {
      fprintf(stderr,"\n  ## Error: %s: Could not read mesh partitioning.\n", __func__);
    }
    H5Gclose(grp_part_id);
    H5Gclose(grp_mesh_id);
    goto free_and_return;
  }

  H5Gclose(grp_part_id);

  chrono(OFF, &ctim[tim]);
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    fprintf(stdout,"  -- Read mesh partitioning.     %s\n",stim);
  }

  /** Compute the offset for parallel reading */
  tim = 2;
  chrono(ON, &ctim[tim]);

  PMMG_CALLOC(parmesh, offset, 2 * PMMG_NTYPENTITIES, hsize_t, "offset",
              goto free_and_return );

  ier = PMMG_computeHDFoffset(parmesh, nentities, offset);

  /* We do not need the number of entities anymore */
  PMMG_DEL_MEM(parmesh, nentities, hsize_t, "nentities");

  chrono(OFF, &ctim[tim]);
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    fprintf(stdout,"  -- Compute HDF5 read offsets.     %s\n",stim);
  }

  /** Load the mesh entities */
  tim = 3;
  chrono(ON, &ctim[tim]);

  /* Each proc reads the part of the mesh that is assigned to him */
  HDF_CHECK( grp_entities_id = H5Gopen(grp_mesh_id, "MeshEntities", H5P_DEFAULT),
             fprintf(stderr,"\n  ## Error: %s: Rank %d could not open the /Mesh/MeshEntities group in file %s.\n",
                     __func__, parmesh->myrank, filename);
             H5Gclose(grp_mesh_id);
             goto free_and_return );

  ier = PMMG_loadMeshEntities_hdf5(parmesh, grp_entities_id, dxpl_id, nentitiesl, nentitiesg, offset, load_entities);

  MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, &ier, 1, MPI_INT, MPI_MIN, read_comm),
             H5Gclose(grp_entities_id);
             H5Gclose(grp_mesh_id);
             goto free_and_return );

  if (!ier) {
    if (parmesh->myrank == parmesh->info.root) {
      fprintf(stderr,"\n  ## Error: %s: Could not read the mesh entities.\n", __func__);
    }
    H5Gclose(grp_entities_id);
    H5Gclose(grp_mesh_id);
    goto free_and_return;
  }

  H5Gclose(grp_entities_id);

  chrono(OFF, &ctim[tim]);
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    fprintf(stdout,"  -- Read mesh entites.     %s\n",stim);
  }

  /* Close the mesh group */
  H5Gclose(grp_mesh_id);

  /** Load the metric and the solutions */
  tim = 4;
  chrono(ON, &ctim[tim]);

  /* Each proc reads the part of the solutions/metric that is assigned to him */
  HDF_CHECK( grp_sols_id = H5Gopen(file_id, "Solutions", H5P_DEFAULT),
             fprintf(stderr,"\n  ## Error: %s: Could not open the /Solutions group in file %s.\n",
                     __func__, filename);
             goto free_and_return );

  ier = PMMG_loadMetric_hdf5(parmesh, grp_sols_id, dxpl_id, nentitiesl, offset);

  MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, &ier, 1, MPI_INT, MPI_MIN, read_comm),
             H5Gclose(grp_sols_id);
             goto free_and_return );

  if ( ier == -1 ) {
    if (parmesh->myrank == parmesh->info.root) {
      fprintf(stderr,"\n  ## Warning: %s: Metric not found, using default metric.\n",__func__);
    }
  }
  if ( ier == 0 ) {
    if (parmesh->myrank == parmesh->info.root) {
      fprintf(stderr,"\n  ## Error: %s: Could not load the metric.\n",__func__);
    }
    H5Gclose(grp_sols_id);
    goto free_and_return;
  }

  ier = PMMG_loadAllSols_hdf5(parmesh, grp_sols_id, dxpl_id, nentitiesl, nentitiesg, offset);

  MPI_CHECK( MPI_Allreduce(MPI_IN_PLACE, &ier, 1, MPI_INT, MPI_MIN, read_comm),
             H5Gclose(grp_sols_id);
             goto free_and_return );

  if (!ier) {
    if (parmesh->myrank == parmesh->info.root) {
      fprintf(stderr,"\n  ## Error: %s: Could not load the solutions.\n",__func__);
    }
    H5Gclose(grp_sols_id);
    goto free_and_return;
  }

  H5Gclose(grp_sols_id);

  chrono(OFF, &ctim[tim]);
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    fprintf(stdout,"  -- Read metric and solutions.     %s\n",stim);
  }

  /*------------------------- RELEASE ALL HDF5 IDs AND MEMORY -------------------------*/

  H5Fclose(file_id);
  H5Pclose(fapl_id);
  H5Pclose(dxpl_id);
  PMMG_DEL_MEM(parmesh, nentitiesl, hsize_t, "nentitiesl");
  PMMG_DEL_MEM(parmesh, nentitiesg, hsize_t, "nentitiesg");
  PMMG_DEL_MEM(parmesh, offset, hsize_t, "offset");

  if (parmesh->info.imprim > PMMG_VERB_VERSION) {
    fprintf(stdout, "  %%%% %s CLOSED \n\n", filename);
  }

  /* Very ugly : if the rank is above the number of partitions of the input mesh,
     allocate an internal communicator of opposite type of API_mode. This is necessary
     because all those processes wont enter the PMMG_preprocessMesh_distributed function.
     Also set ngrp to 0 for loadBalancing to work properly. */
  if ( rank >= parmesh->info.npartin ) {
    parmesh->ngrp = 0;
    if (parmesh->info.API_mode == PMMG_APIDISTRIB_faces)
      PMMG_CALLOC(parmesh, parmesh->int_node_comm, 1, PMMG_Int_comm, "int_face_comm", goto free_and_return);
    else
      PMMG_CALLOC(parmesh, parmesh->int_face_comm, 1, PMMG_Int_comm, "int_face_comm", goto free_and_return);
  }

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( parmesh->info.imprim >= PMMG_VERB_STEPS ) {
    fprintf(stdout,"\n   LOAD_PARMESH_HDF5: ELAPSED TIME  %s\n",stim);
  }

  return 1;

 free_and_return:
  H5Fclose(file_id);
  H5Pclose(fapl_id);
  H5Pclose(dxpl_id);
  PMMG_DEL_MEM(parmesh, nentities, hsize_t, "nentities");
  PMMG_DEL_MEM(parmesh, nentitiesg, hsize_t, "nentitiesg");
  PMMG_DEL_MEM(parmesh, nentitiesl, hsize_t, "nentitiesl");
  PMMG_DEL_MEM(parmesh, offset, hsize_t, "offset");
  return 0;

#endif

}

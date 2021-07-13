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
#include "hdf5.h"

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

int PMMG_Get_parmeshSize(PMMG_pParMesh parmesh, int ntyp_entities,
                         hsize_t *np, hsize_t *na, hsize_t *nt, hsize_t *nquad, hsize_t *ne, hsize_t *nprism,
                         hsize_t *nc, hsize_t *nreq, hsize_t *npar,
                         hsize_t *nr, hsize_t *nedreq, hsize_t *nedpar,
                         hsize_t *ntreq, hsize_t *ntpar,
                         hsize_t *nqreq, hsize_t *nqpar,
                         hsize_t *nereq, hsize_t *nepar,
                         hsize_t *nnor, hsize_t *ntan,
                         hsize_t *npg, hsize_t *nag, hsize_t *ntg, hsize_t *nquadg, hsize_t *neg, hsize_t *nprismg,
                         hsize_t *ncg, hsize_t *nreqg, hsize_t *nparg,
                         hsize_t *nrg, hsize_t *nedreqg, hsize_t *nedparg,
                         hsize_t *ntreqg, hsize_t *ntparg,
                         hsize_t *nqreqg, hsize_t *nqparg,
                         hsize_t *nereqg, hsize_t *neparg,
                         hsize_t *nnorg, hsize_t *ntang,
                         hsize_t *point_offset, hsize_t *edge_offset, hsize_t *tria_offset,
                         hsize_t *quad_offset, hsize_t *tetra_offset, hsize_t *prism_offset,
                         hsize_t *corner_offset, hsize_t *ridge_offset,
                         hsize_t *required_offset, hsize_t *parallel_offset) {

  /* MMG variables */
  PMMG_pGrp grp   = &parmesh->listgrp[0];
  MMG5_pMesh mesh = grp->mesh;
  MMG5_pPoint ppt = NULL;
  MMG5_pEdge pa   = NULL;
  MMG5_pTria pt   = NULL;
  MMG5_pQuad pq   = NULL;
  MMG5_pTetra pe  = NULL;
  MMG5_pPrism pp  = NULL;
  hsize_t *nentities;

  /* MPI variables */
  int nprocs = parmesh->nprocs;
  int rank = parmesh->myrank;
  MPI_Comm comm = parmesh->comm;

  nentities = NULL;

  /*------------------------- COUNT LOCAL MESH ENTITIES -------------------------*/

  if ( !mesh->point ) {
    fprintf(stderr, "\n  ## Error: %s: points array not allocated.\n",
            __func__);
    return 0;
  }

  /* Vertices, normals and tangents */
  for (int k = 1 ; k <= mesh->np ; k++) {
    ppt = &mesh->point[k];
    if (MG_VOK(ppt)) {
      ppt->tmp = ++(*np);
      ppt->flag = 0;
      if (ppt->tag & MG_CRN) (*nc)++;
      if (ppt->tag & MG_REQ) (*nreq)++;
      if (ppt->tag & MG_PARBDY) (*npar)++;
      if ((!ppt->flag) || MG_SIN(ppt->tag)) continue;
      if (ppt->tag & MG_BDY && (!(ppt->tag & MG_GEO) || ppt->tag & MG_NOM)) (*nnor)++;
      if (MG_EDG(ppt->tag) || (ppt->tag & MG_NOM)) (*ntan)++;
    }
  }

  /* Edges */
  if (mesh->na) {
    for (int k = 1 ; k <= mesh->na ; k++) {
      pa = &mesh->edge[k];
      (*na)++;
      if (pa->tag & MG_GEO) (*nr)++;
      if (pa->tag & MG_REQ) (*nedreq)++;
      if (pa->tag & MG_PARBDY) (*nedpar)++;
    }
  }

  /* Triangles */
  if (mesh->nt) {
    for (int k = 1 ; k <= mesh->nt ; k++) {
      pt = &mesh->tria[k];
      (*nt)++;
      if (pt->tag[0] & MG_REQ && pt->tag[1] & MG_REQ && pt->tag[2] & MG_REQ) (*ntreq)++;
      if (pt->tag[0] & MG_PARBDY && pt->tag[1] & MG_PARBDY && pt->tag[2] & MG_PARBDY) (*ntpar)++;
    }
  }

  /* Quadrilaterals */
  if (mesh->nquad) {
    for (int k = 1 ; k <= mesh->nquad ; k++) {
      pq = &mesh->quadra[k];
      (*nquad)++;
      if (pq->tag[0] & MG_REQ && pq->tag[1] & MG_REQ &&
          pq->tag[2] & MG_REQ && pq->tag[3] & MG_REQ) {
        (*nqreq)++;
      }
      if (pq->tag[0] & MG_PARBDY && pq->tag[1] & MG_PARBDY &&
          pq->tag[2] & MG_PARBDY && pq->tag[3] & MG_PARBDY) {
        (*nqpar)++;
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
      (*ne)++;
      if (pe->tag & MG_REQ) (*nereq)++;
      if (pe->tag & MG_PARBDY) (*nepar)++;
    }
  } else {
    fprintf(stderr, "\n  ## Warning: %s: tetra array not allocated.\n",
            __func__);
  }

  /* Prisms */
  if (mesh->nprism) {
    for (int k = 1 ; k <= mesh->nprism ; k++) {
      pp = &mesh->prism[k];
      if (!MG_EOK(pp)){
        continue;
      }
      (*nprism)++;
    }
  }

  /*------------------------- COUNT GLOBAL MESH ENTITIES -------------------------*/

  nentities = (hsize_t*) calloc(ntyp_entities * nprocs, sizeof(hsize_t));

  nentities[ntyp_entities * rank]      = *np;
  nentities[ntyp_entities * rank + 1]  = *na;
  nentities[ntyp_entities * rank + 2]  = *nt;
  nentities[ntyp_entities * rank + 3]  = *nquad;
  nentities[ntyp_entities * rank + 4]  = *ne;
  nentities[ntyp_entities * rank + 5]  = *nprism;
  nentities[ntyp_entities * rank + 6]  = *nc;
  nentities[ntyp_entities * rank + 7]  = *nreq;
  nentities[ntyp_entities * rank + 8]  = *npar;
  nentities[ntyp_entities * rank + 9]  = *nr;
  nentities[ntyp_entities * rank + 10] = *nedreq;
  nentities[ntyp_entities * rank + 11] = *nedpar;
  nentities[ntyp_entities * rank + 12] = *ntreq;
  nentities[ntyp_entities * rank + 13] = *ntpar;
  nentities[ntyp_entities * rank + 14] = *nqreq;
  nentities[ntyp_entities * rank + 15] = *nqpar;
  nentities[ntyp_entities * rank + 16] = *nereq;
  nentities[ntyp_entities * rank + 17] = *nepar;
  nentities[ntyp_entities * rank + 18] = *nnor;
  nentities[ntyp_entities * rank + 19] = *ntan;

  MPI_Allgather(&nentities[ntyp_entities * rank], ntyp_entities, MPI_UNSIGNED_LONG_LONG,
                nentities                       , ntyp_entities, MPI_UNSIGNED_LONG_LONG, comm);

  for (int k = 0 ; k < nprocs ; k++) {
    *npg     += nentities[ntyp_entities * k];
    *nag     += nentities[ntyp_entities * k + 1];
    *ntg     += nentities[ntyp_entities * k + 2];
    *nquadg  += nentities[ntyp_entities * k + 3];
    *neg     += nentities[ntyp_entities * k + 4];
    *nprismg += nentities[ntyp_entities * k + 5];
    *ncg     += nentities[ntyp_entities * k + 6];
    *nreqg   += nentities[ntyp_entities * k + 7];
    *nparg   += nentities[ntyp_entities * k + 8];
    *nrg     += nentities[ntyp_entities * k + 9];
    *nedreqg += nentities[ntyp_entities * k + 10];
    *nedparg += nentities[ntyp_entities * k + 11];
    *ntreqg  += nentities[ntyp_entities * k + 12];
    *ntparg  += nentities[ntyp_entities * k + 13];
    *nqreqg  += nentities[ntyp_entities * k + 14];
    *nqparg  += nentities[ntyp_entities * k + 15];
    *nereqg  += nentities[ntyp_entities * k + 16];
    *neparg  += nentities[ntyp_entities * k + 17];
    *nnorg   += nentities[ntyp_entities * k + 18];
    *ntang   += nentities[ntyp_entities * k + 19];
  }

  /*------------------------- COMPUTE OFFSET ARRAYS -------------------------*/

  for (int k = 0 ; k < rank ; k++) {
    point_offset[0]    += nentities[ntyp_entities * k];
    edge_offset[0]     += nentities[ntyp_entities * k + 1];
    tria_offset[0]     += nentities[ntyp_entities * k + 2];
    quad_offset[0]     += nentities[ntyp_entities * k + 3];
    tetra_offset[0]    += nentities[ntyp_entities * k + 4];
    prism_offset[0]    += nentities[ntyp_entities * k + 5];
    corner_offset      += nentities[ntyp_entities * k + 6];
    required_offset[0] += nentities[ntyp_entities * k + 7];
    parallel_offset[0] += nentities[ntyp_entities * k + 8];
    ridge_offset       += nentities[ntyp_entities * k + 9];
    required_offset[1] += nentities[ntyp_entities * k + 10];
    parallel_offset[1] += nentities[ntyp_entities * k + 11];
    required_offset[2] += nentities[ntyp_entities * k + 12];
    parallel_offset[2] += nentities[ntyp_entities * k + 13];
    required_offset[3] += nentities[ntyp_entities * k + 14];
    parallel_offset[3] += nentities[ntyp_entities * k + 15];
    required_offset[4] += nentities[ntyp_entities * k + 16];
    parallel_offset[4] += nentities[ntyp_entities * k + 17];
  }

  /* We no longer need the number of entities array */
  free(nentities); nentities = NULL;

}


int PMMG_saveParmesh_hdf5(PMMG_pParMesh parmesh, const char *filename, const char *xdmfname) {
  /* MMG variables */
  int ier = 1;
  int nsols;
  PMMG_pGrp grp;
  PMMG_pExt_comm comms;
  MMG5_pSol met, sols;
  MMG5_pMesh mesh;
  MMG5_pPoint ppt;
  MMG5_pEdge pa;
  MMG5_pTria pt;
  MMG5_pQuad pq;
  MMG5_pTetra pe;
  MMG5_pPrism pp;

  /* Local mesh size */
  hsize_t ne, np, nt, na, nquad, nprism;       /* Tetra, points, triangles, edges, quads, prisms */
  hsize_t nc, nreq, npar;                      /* Corners, required and parallel vertices */
  hsize_t nr, nedreq, nedpar;                  /* Ridges, required and parallel edges */
  hsize_t ntreq, ntpar;                        /* Required and parallel triangles */
  hsize_t nqreq, nqpar;                        /* Required and parallel quads */
  hsize_t nereq, nepar;                        /* Required and parallel tetra */
  hsize_t nnor, ntan;                          /* Normals and Tangents */

  /* Global mesh size */
  hsize_t neg, npg, ntg, nag, nquadg, nprismg; /* Tetra, points, triangles, edges, quads, prisms */
  hsize_t ncg, nreqg, nparg;                   /* Corners, required and parallel vertices */
  hsize_t nrg, nedreqg, nedparg;               /* Ridges, required and parallel edges */
  hsize_t ntreqg, ntparg;                      /* Required and parallel triangles */
  hsize_t nqreqg, nqparg;                      /* Required and parallel quads */
  hsize_t nereqg, neparg;                      /* Required and parallel tetra */
  hsize_t nnorg, ntang;                        /* Normals and Tangents */

  /* Buffer arrays */
  /* 6 buffers is the minimum amount for what we have to do */
  double *ppoint;   /* Point coordinates */
  int *pent;        /* Other entities : edges, trias, quads, tetra, prisms. */
  int *pcr;         /* Corners and ridges */
  int *preq, *ppar; /* Required and parallel entities */
  int *pref;        /* References */

  /* Store the number of entities */
  /* This only serves as a buffer for the MPI communication */
  int ntyp_entities = 20;
  hsize_t *nentities;

  /* Counters for the corners/ridges, the required entities and the parallel entities */
  int crcount, reqcount, parcount;

  /* Offsets for parallel writing */
  hsize_t point_offset[3] = {0, 0, 0};
  hsize_t edge_offset[2]  = {0, 0};
  hsize_t tria_offset[3]  = {0, 0, 0};
  hsize_t quad_offset[4]  = {0, 0, 0, 0};
  hsize_t tetra_offset[4] = {0, 0, 0, 0};
  hsize_t prism_offset[6] = {0, 0, 0, 0, 0, 0};
  hsize_t required_offset[5] = {0, 0, 0, 0, 0};     /* Used for the required entities */
  hsize_t parallel_offset[5] = {0, 0, 0, 0, 0};     /* Used for the parallel entities */
  hsize_t corner_offset = 0;                        /* Used for the corners */
  hsize_t ridge_offset = 0;                         /* Used for the ridges */

  /* HDF5 variables */
  hid_t file_id, grp_mesh_id, grp_comm_id, grp_entities_id, grp_sols_id; /* Objects */
  hid_t fapl_id, dxpl_id, dcpl_id;                                       /* Property lists */
  hid_t attr_dim_id, attr_ver_id;                                        /* Attributes */
  hid_t dspace_mem_id, dspace_file_id;                                   /* Dataspaces */
  hid_t dset_id;                                                         /* Dataset */
  herr_t status;

  /* MPI variables */
  MPI_Info info = MPI_INFO_NULL;
  MPI_Comm comm = parmesh->comm;
  int rank, root, nprocs;

  /*------------------------- INIT -------------------------*/

  /* Set all buffers to NULL */
  ppoint = NULL;
  pent = NULL;
  pcr = NULL;
  preq = NULL; ppar = NULL;
  pref = NULL;
  nentities = NULL;

  /* Set MPI variables */
  nprocs = parmesh->nprocs;
  rank = parmesh->myrank;
  root = parmesh->info.root;

  /* Set mesh size to 0 */
  np = na = nt = nquad = ne = nprism = 0;
  nc = nreq = npar = 0;
  nr = nedreq = nedpar = 0;
  ntreq = ntpar = 0;
  nqreq = nqpar = 0;
  nereq = nepar = 0;
  nnor = ntan = 0;

  npg = nag = ntg = nquadg = neg = nprismg = 0;
  ncg = nreqg = nparg = 0;
  nrg = nedreqg = nedparg = 0;
  ntreqg = ntparg = 0;
  nqreqg = nqparg = 0;
  nereqg = neparg = 0;
  nnorg = ntang = 0;

  /* Check arguments */
  if (parmesh->ngrp != 1) {
    fprintf(stderr,"  ## Error: %s: you must have exactly 1 group in you parmesh.",
            __func__);
    return 0;
  }
  if (!filename || !*filename) {
    fprintf(stderr,"  ## Error: %s: no HDF5 file name provided.",
            __func__);
    return 0;
  }

  /* Set ParMmg variables */
  grp = &parmesh->listgrp[0];
  comms = parmesh->ext_face_comm;
  mesh = grp->mesh;
  met = grp->met;
  sols = grp->field;
  nsols = mesh->nsols;
  ppt = NULL;
  pa = NULL;
  pt = NULL;
  pq = NULL;
  pe = NULL;
  pp = NULL;

  /*------------------------- COUNT LOCAL MESH ENTITIES -------------------------*/

  if ( !mesh->point ) {
    fprintf(stderr, "\n  ## Error: %s: points array not allocated.\n",
            __func__);
    return 0;
  }

  /* /\* Vertices, normals and tangents *\/ */
  /* for (int k = 1 ; k <= mesh->np ; k++) { */
  /*   ppt = &mesh->point[k]; */
  /*   if (MG_VOK(ppt)) { */
  /*     ppt->tmp = ++np; */
  /*     ppt->flag = 0; */
  /*     if (ppt->tag & MG_CRN) nc++; */
  /*     if (ppt->tag & MG_REQ) nreq++; */
  /*     if (ppt->tag & MG_PARBDY) npar++; */
  /*     if ((!ppt->flag) || MG_SIN(ppt->tag)) continue; */
  /*     if (ppt->tag & MG_BDY && (!(ppt->tag & MG_GEO) || ppt->tag & MG_NOM)) nnor++; */
  /*     if (MG_EDG(ppt->tag) || (ppt->tag & MG_NOM)) ntan++; */
  /*   } */
  /* } */

  /* /\* Edges *\/ */
  /* if (mesh->na) { */
  /*   for (int k = 1 ; k <= mesh->na ; k++) { */
  /*     pa = &mesh->edge[k]; */
  /*     na++; */
  /*     if (pa->tag & MG_GEO) nr++; */
  /*     if (pa->tag & MG_REQ) nedreq++; */
  /*     if (pa->tag & MG_PARBDY) nedpar++; */
  /*   } */
  /* } */

  /* /\* Triangles *\/ */
  /* if (mesh->nt) { */
  /*   for (int k = 1 ; k <= mesh->nt ; k++) { */
  /*     pt = &mesh->tria[k]; */
  /*     nt++; */
  /*     if (pt->tag[0] & MG_REQ && pt->tag[1] & MG_REQ && pt->tag[2] & MG_REQ) ntreq++; */
  /*     if (pt->tag[0] & MG_PARBDY && pt->tag[1] & MG_PARBDY && pt->tag[2] & MG_PARBDY) ntpar++; */
  /*   } */
  /* } */

  /* /\* Quadrilaterals *\/ */
  /* if (mesh->nquad) { */
  /*   for (int k = 1 ; k <= mesh->nquad ; k++) { */
  /*     pq = &mesh->quadra[k]; */
  /*     nquad++; */
  /*     if (pq->tag[0] & MG_REQ && pq->tag[1] & MG_REQ && */
  /*         pq->tag[2] & MG_REQ && pq->tag[3] & MG_REQ) { */
  /*       nqreq++; */
  /*     } */
  /*     if (pq->tag[0] & MG_PARBDY && pq->tag[1] & MG_PARBDY && */
  /*         pq->tag[2] & MG_PARBDY && pq->tag[3] & MG_PARBDY) { */
  /*       nqpar++; */
  /*     } */
  /*   } */
  /* } */

  /* /\* Tetrahedra *\/ */
  /* if (mesh->ne) { */
  /*   for (int k = 1 ; k <= mesh->ne ; k++) { */
  /*     pe = &mesh->tetra[k]; */
  /*     if (!MG_EOK(pe)) { */
  /*       continue; */
  /*     } */
  /*     ne++; */
  /*     if (pe->tag & MG_REQ) nereq++; */
  /*     if (pe->tag & MG_PARBDY) nepar++; */
  /*   } */
  /* } else { */
  /*   fprintf(stderr, "\n  ## Warning: %s: tetra array not allocated.\n", */
  /*           __func__); */
  /* } */

  /* /\* Prisms *\/ */
  /* if (mesh->nprism) { */
  /*   for (int k = 1 ; k <= mesh->nprism ; k++) { */
  /*     pp = &mesh->prism[k]; */
  /*     if (!MG_EOK(pp)){ */
  /*       continue; */
  /*     } */
  /*     nprism++; */
  /*   } */
  /* } */

  /*------------------------- COUNT GLOBAL MESH ENTITIES -------------------------*/

  /* nentities = (hsize_t*) calloc(ntyp_entities * nprocs, sizeof(hsize_t)); */

  /* nentities[ntyp_entities * rank]      = np; */
  /* nentities[ntyp_entities * rank + 1]  = na; */
  /* nentities[ntyp_entities * rank + 2]  = nt; */
  /* nentities[ntyp_entities * rank + 3]  = nquad; */
  /* nentities[ntyp_entities * rank + 4]  = ne; */
  /* nentities[ntyp_entities * rank + 5]  = nprism; */
  /* nentities[ntyp_entities * rank + 6]  = nc; */
  /* nentities[ntyp_entities * rank + 7]  = nreq; */
  /* nentities[ntyp_entities * rank + 8]  = npar; */
  /* nentities[ntyp_entities * rank + 9]  = nr; */
  /* nentities[ntyp_entities * rank + 10] = nedreq; */
  /* nentities[ntyp_entities * rank + 11] = nedpar; */
  /* nentities[ntyp_entities * rank + 12] = ntreq; */
  /* nentities[ntyp_entities * rank + 13] = ntpar; */
  /* nentities[ntyp_entities * rank + 14] = nqreq; */
  /* nentities[ntyp_entities * rank + 15] = nqpar; */
  /* nentities[ntyp_entities * rank + 16] = nereq; */
  /* nentities[ntyp_entities * rank + 17] = nepar; */
  /* nentities[ntyp_entities * rank + 18] = nnor; */
  /* nentities[ntyp_entities * rank + 19] = ntan; */

  /* MPI_Allgather(&nentities[ntyp_entities * rank], ntyp_entities, MPI_UNSIGNED_LONG_LONG, */
  /*               nentities                       , ntyp_entities, MPI_UNSIGNED_LONG_LONG, comm); */

  /* for (int k = 0 ; k < nprocs ; k++) { */
  /*   npg     += nentities[ntyp_entities * k]; */
  /*   nag     += nentities[ntyp_entities * k + 1]; */
  /*   ntg     += nentities[ntyp_entities * k + 2]; */
  /*   nquadg  += nentities[ntyp_entities * k + 3]; */
  /*   neg     += nentities[ntyp_entities * k + 4]; */
  /*   nprismg += nentities[ntyp_entities * k + 5]; */
  /*   ncg     += nentities[ntyp_entities * k + 6]; */
  /*   nreqg   += nentities[ntyp_entities * k + 7]; */
  /*   nparg   += nentities[ntyp_entities * k + 8]; */
  /*   nrg     += nentities[ntyp_entities * k + 9]; */
  /*   nedreqg += nentities[ntyp_entities * k + 10]; */
  /*   nedparg += nentities[ntyp_entities * k + 11]; */
  /*   ntreqg  += nentities[ntyp_entities * k + 12]; */
  /*   ntparg  += nentities[ntyp_entities * k + 13]; */
  /*   nqreqg  += nentities[ntyp_entities * k + 14]; */
  /*   nqparg  += nentities[ntyp_entities * k + 15]; */
  /*   nereqg  += nentities[ntyp_entities * k + 16]; */
  /*   neparg  += nentities[ntyp_entities * k + 17]; */
  /*   nnorg   += nentities[ntyp_entities * k + 18]; */
  /*   ntang   += nentities[ntyp_entities * k + 19]; */
  /* } */

  /*------------------------- COMPUTE OFFSET ARRAYS -------------------------*/

  /* for (int k = 0 ; k < rank ; k++) { */
  /*   point_offset[0]    += nentities[ntyp_entities * k]; */
  /*   edge_offset[0]     += nentities[ntyp_entities * k + 1]; */
  /*   tria_offset[0]     += nentities[ntyp_entities * k + 2]; */
  /*   quad_offset[0]     += nentities[ntyp_entities * k + 3]; */
  /*   tetra_offset[0]    += nentities[ntyp_entities * k + 4]; */
  /*   prism_offset[0]    += nentities[ntyp_entities * k + 5]; */
  /*   corner_offset      += nentities[ntyp_entities * k + 6]; */
  /*   required_offset[0] += nentities[ntyp_entities * k + 7]; */
  /*   parallel_offset[0] += nentities[ntyp_entities * k + 8]; */
  /*   ridge_offset       += nentities[ntyp_entities * k + 9]; */
  /*   required_offset[1] += nentities[ntyp_entities * k + 10]; */
  /*   parallel_offset[1] += nentities[ntyp_entities * k + 11]; */
  /*   required_offset[2] += nentities[ntyp_entities * k + 12]; */
  /*   parallel_offset[2] += nentities[ntyp_entities * k + 13]; */
  /*   required_offset[3] += nentities[ntyp_entities * k + 14]; */
  /*   parallel_offset[3] += nentities[ntyp_entities * k + 15]; */
  /*   required_offset[4] += nentities[ntyp_entities * k + 16]; */
  /*   parallel_offset[4] += nentities[ntyp_entities * k + 17]; */
  /* } */

  /* /\* We no longer need the number of entities array *\/ */
  /* free(nentities); nentities = NULL; */

  PMMG_Get_parmeshSize(parmesh, ntyp_entities, &np, &na, &nt, &nquad, &ne, &nprism, &nc, &nreq, &npar,
                       &nr, &nedreq, &nedpar, &ntreq, &ntpar, &nqreq, &nqpar, &nereq, &nepar, &nnor, &ntan,
                       &npg, &nag, &ntg, &nquadg, &neg, &nprismg, &ncg, &nreqg, &nparg,
                       &nrg, &nedreqg, &nedparg, &ntreqg, &ntparg, &nqreqg, &nqparg, &nereqg, &neparg, &nnorg, &ntang,
                       point_offset, edge_offset, tria_offset, quad_offset, tetra_offset, prism_offset,
                       &corner_offset, &ridge_offset, required_offset, parallel_offset);

  /* Arrays for bidimensional dataspaces */
  hsize_t hnp[2]      = {np,3};
  hsize_t hna[2]      = {na,2};
  hsize_t hnt[2]      = {nt,3};
  hsize_t hnquad[2]   = {nquad,4};
  hsize_t hne[2]      = {ne,4};
  hsize_t hnprism[2]  = {nprism,2};
  hsize_t hnpg[2]     = {npg,3};
  hsize_t hnag[2]     = {nag,2};
  hsize_t hntg[2]     = {ntg,3};
  hsize_t hnquadg[2]  = {nquadg,4};
  hsize_t hneg[2]     = {neg,4};
  hsize_t hnprismg[2] = {nprismg,2};

  /*------------------------- COMMUNICATORS -------------------------*/

  hsize_t *ncomms, ncommg, comm_offset;
  int *colors, *nface;

  ncommg = comm_offset = 0;

  /* Count the number of communicators */
  ncomms = (hsize_t*) calloc(nprocs, sizeof(hsize_t));
  ncomms[rank] = parmesh->next_face_comm;
  MPI_Allgather(&ncomms[rank], 1, MPI_LONG_LONG, ncomms, 1, MPI_LONG_LONG, comm);

  for (int i = 0 ; i < nprocs ; i++) {
    ncommg += ncomms[i];
  }
  for (int i = 0 ; i < rank ; i++) {
    comm_offset += ncomms[i];
  }

  /* Create the buffers */
  colors = (int*) calloc(ncomms[rank], sizeof(int));
  nface = (int*) calloc(ncomms[rank], sizeof(int));

  for (int icomm = 0 ; icomm < ncomms[rank] ; icomm++) {
    colors[icomm] = comms[icomm].color_out;
    nface[icomm] = comms[icomm].nitem;
  }

  /*------------------------- HDF5 IOs START HERE -------------------------*/

  /* Shut HDF5 error stack */
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);

  /* Create the parallel file acces and the parallel dataset transfer property lists */
  fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(fapl_id, comm, info);
  dxpl_id = H5Pcreate(H5P_DATASET_XFER);
  status = H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);

  /* Create the dataset creation property list to tell we don't want to write any fill value */
  dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  status = H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER);

  /*------------------------- OPEN FILE AND WRITE DATA -------------------------*/

  /* Create the file */
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);

  /* Save the attributes (Version and Dimension) */
  dspace_file_id = H5Screate(H5S_SCALAR);
  attr_ver_id = H5Acreate(file_id, "MeshVersionFormatted", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, H5P_DEFAULT);
  if (rank == root)
    status = H5Awrite(attr_ver_id, H5T_NATIVE_INT, &mesh->ver);
  H5Aclose(attr_ver_id);
  attr_dim_id = H5Acreate(file_id, "Dimension", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, H5P_DEFAULT);
  if (rank == root)
    status = H5Awrite(attr_dim_id, H5T_NATIVE_INT, &mesh->dim);
  H5Aclose(attr_dim_id);
  H5Sclose(dspace_file_id);

  /*------------------------- WRITE MESH -------------------------*/

  grp_mesh_id = H5Gcreate(file_id, "Mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /*------------------------- WRITE COMMUNICATORS -------------------------*/

  grp_comm_id = H5Gcreate(grp_mesh_id, "FaceCommunicators", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* Number of communicators */
  hsize_t hnprocs = nprocs;
  dspace_file_id = H5Screate_simple(1, &hnprocs, NULL);
  dset_id = H5Dcreate(grp_comm_id, "NumberOfFaceCommunicators", H5T_NATIVE_LLONG, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  if (rank == root)
    status = H5Dwrite(dset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ncomms[rank]);
  H5Dclose(dset_id);
  H5Sclose(dspace_file_id);

  /* For each communicator, write the number of faces and the outward proc color */
  dspace_mem_id  = H5Screate_simple(1, &ncomms[rank], NULL);
  dspace_file_id = H5Screate_simple(1, &ncommg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &comm_offset, NULL, &ncomms[rank], NULL);

  dset_id = H5Dcreate(grp_comm_id, "ColorsOut", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, colors);
  H5Dclose(dset_id);

  dset_id = H5Dcreate(grp_comm_id, "NumberOfCommunicatorFaces", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, nface);
  H5Dclose(dset_id);

  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  H5Gclose(grp_comm_id);

  free(ncomms); ncomms = NULL;
  free(colors); colors = NULL;
  free(nface); nface = NULL;

  /*------------------------- WRITE MESH ENTITIES -------------------------*/

  grp_entities_id = H5Gcreate(grp_mesh_id, "MeshEntities", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* Vertices */
  ppoint = (double*) calloc(np, 3 * sizeof(double));
  pref = (int*) calloc(np, sizeof(int));
  pcr = (int*) calloc(nc, sizeof(int));
  preq = (int*) calloc(nreq, sizeof(int));
  ppar = (int*) calloc(npar, sizeof(int));

  crcount = reqcount = parcount = 0;

  for (int i = 0 ; i < mesh->np ; i++) {
    ppt = &mesh->point[i + 1];
    if (MG_VOK(ppt)){
      for (int j = 0 ; j < 3 ; j++) {
        ppoint[3 * (ppt->tmp - 1) + j] = ppt->c[j];
      }
      if (ppt->tag & MG_CRN)    pcr[crcount++] = ppt->tmp + point_offset[0] - 1;
      if (ppt->tag & MG_REQ)    preq[reqcount++]   = ppt->tmp + point_offset[0] - 1;
      if (ppt->tag & MG_PARBDY) ppar[parcount++]   = ppt->tmp - 1; /* Local index for parallel entities */
      pref[ppt->tmp - 1] = abs(ppt->ref);
    }
  }

  dspace_mem_id  = H5Screate_simple(2, hnp, NULL);
  dspace_file_id = H5Screate_simple(2, hnpg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, point_offset, NULL, hnp, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Vertices", H5T_NATIVE_DOUBLE, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dspace_mem_id, dspace_file_id, dxpl_id, ppoint);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(ppoint); ppoint = NULL;

  dspace_mem_id  = H5Screate_simple(1, hnp, NULL);
  dspace_file_id = H5Screate_simple(1, hnpg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, point_offset, NULL, hnp, NULL);
  dset_id = H5Dcreate(grp_entities_id, "VerticesRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(pref); pref = NULL;

  dspace_mem_id  = H5Screate_simple(1, &nc, NULL);
  dspace_file_id = H5Screate_simple(1, &ncg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &corner_offset, NULL, &nc, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Corners", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pcr);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(pcr) ; pcr = NULL;

  dspace_mem_id  = H5Screate_simple(1, &nreq, NULL);
  dspace_file_id = H5Screate_simple(1, &nreqg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &required_offset[0], NULL, &nreq, NULL);
  dset_id = H5Dcreate(grp_entities_id, "RequiredVertices", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(preq) ; preq = NULL;

  dspace_mem_id  = H5Screate_simple(1, &npar, NULL);
  dspace_file_id = H5Screate_simple(1, &nparg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &parallel_offset[0], NULL, &npar, NULL);
  dset_id = H5Dcreate(grp_entities_id, "ParallelVertices", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(ppar) ; ppar = NULL;

  /* Edges */
  pent = (int*) calloc(na, 2 * sizeof(int));
  pref = (int*) calloc(na, sizeof(int));
  pcr = (int*) calloc(nr, sizeof(int));
  preq = (int*) calloc(nedreq, sizeof(int));
  ppar = (int*) calloc(nedpar, sizeof(int));

  crcount = reqcount = parcount = 0;

  if (na) {
    na = 0;
    for (int i = 0 ; i < mesh->na ; i++) {
      pa = &mesh->edge[i + 1];
      pent[2 * i]     = mesh->point[pa->a].tmp + point_offset[0] - 1;
      pent[2 * i + 1] = mesh->point[pa->b].tmp + point_offset[0] - 1;
      pref[i] = pa->ref;
      if (pa->tag & MG_GEO)    pcr[crcount++] = na + edge_offset[0];
      if (pa->tag & MG_REQ)    preq[reqcount++] = na + edge_offset[0];
      if (pa->tag & MG_PARBDY) ppar[parcount++] = na; /* Local index for parallel entities */
      na++;
    }
  }

  dspace_mem_id  = H5Screate_simple(2, hna, NULL);
  dspace_file_id = H5Screate_simple(2, hnag, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, edge_offset, NULL, hna, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Edges", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(pent); pent = NULL;

  dspace_mem_id  = H5Screate_simple(1, hna, NULL);
  dspace_file_id = H5Screate_simple(1, hnag, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, edge_offset, NULL, hna, NULL);
  dset_id = H5Dcreate(grp_entities_id, "EdgesRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(pref); pref = NULL;

  dspace_mem_id  = H5Screate_simple(1, &nr, NULL);
  dspace_file_id = H5Screate_simple(1, &nrg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &ridge_offset, NULL, &nr, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Ridges", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pcr);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(pcr); pcr = NULL;

  dspace_mem_id  = H5Screate_simple(1, &nedreq, NULL);
  dspace_file_id = H5Screate_simple(1, &nedreqg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &required_offset[1], NULL, &nedreq, NULL);
  dset_id = H5Dcreate(grp_entities_id, "RequiredEdges", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(preq); preq = NULL;

  dspace_mem_id  = H5Screate_simple(1, &nedpar, NULL);
  dspace_file_id = H5Screate_simple(1, &nedparg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &parallel_offset[1], NULL, &nedpar, NULL);
  dset_id = H5Dcreate(grp_entities_id, "ParallelEdges", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(ppar); ppar = NULL;

  /* Triangles */
  pent = (int*) calloc(nt, 3 * sizeof(int));
  pref = (int*) calloc(nt, sizeof(int));
  preq = (int*) calloc(ntreq, sizeof(int));
  ppar = (int*) calloc(ntpar, sizeof(int));

  reqcount = parcount = 0;

  if (nt) {
    nt = 0;
    for (int i = 0 ; i < mesh->nt ; i++) {
      pt = &mesh->tria[i + 1];
      for (int j = 0 ; j < 3 ; j++) {
        pent[3 * i + j] = mesh->point[pt->v[j]].tmp + point_offset[0] - 1;
      }
      pref[i] = pt->ref;
      if (pt->tag[0] & MG_REQ && pt->tag[1] & MG_REQ && pt->tag[2] & MG_REQ) {
        preq[reqcount++] = nt + tria_offset[0];
      }
      if (pt->tag[0] & MG_PARBDY && pt->tag[1] & MG_PARBDY && pt->tag[2] & MG_PARBDY) {
        ppar[parcount++] = nt; /* Local index for parallel entities */
      }
      nt++;
    }
  }

  dspace_mem_id  = H5Screate_simple(2, hnt, NULL);
  dspace_file_id = H5Screate_simple(2, hntg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, tria_offset, NULL, hnt, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Triangles", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(pent); pent = NULL;

  dspace_mem_id  = H5Screate_simple(1, hnt, NULL);
  dspace_file_id = H5Screate_simple(1, hntg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, tria_offset, NULL, hnt, NULL);
  dset_id = H5Dcreate(grp_entities_id, "TrianglesRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(pref); pref = NULL;

  dspace_mem_id  = H5Screate_simple(1, &ntreq, NULL);
  dspace_file_id = H5Screate_simple(1, &ntreqg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &required_offset[2], NULL, &ntreq, NULL);
  dset_id = H5Dcreate(grp_entities_id, "RequiredTriangles", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(preq); preq = NULL;

  dspace_mem_id  = H5Screate_simple(1, &ntpar, NULL);
  dspace_file_id = H5Screate_simple(1, &ntparg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &parallel_offset[2], NULL, &ntpar, NULL);
  dset_id = H5Dcreate(grp_entities_id, "ParallelTriangles", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(ppar); ppar = NULL;

  /* Quadrilaterals */
  pent = (int*) calloc(nquad, 4 * sizeof(int));
  pref = (int*) calloc(nquad, sizeof(int));
  preq = (int*) calloc(nqreq, sizeof(int));
  ppar = (int*) calloc(nqpar, sizeof(int));

  reqcount = parcount = 0;

  if (nquad){
    nquad = 0;
    for (int i = 0 ; i < mesh->nquad ; i++) {
      pq = &mesh->quadra[i + 1];
      for (int j = 0 ; j < 4 ; j++) {
        pent[4 * i + j] = mesh->point[pq->v[j]].tmp + point_offset[0] - 1;
      }
      pref[i] = pq->ref;
      if (pq->tag[0] & MG_REQ && pq->tag[1] & MG_REQ &&
          pq->tag[2] & MG_REQ && pq->tag[3] & MG_REQ) {
        preq[reqcount++] = nquad + quad_offset[0];
      }
      if (pq->tag[0] & MG_PARBDY && pq->tag[1] & MG_PARBDY &&
          pq->tag[2] & MG_PARBDY && pq->tag[3] & MG_PARBDY) {
        ppar[parcount++] = nquad; /* Local index for parallel entities */
      }
      nquad++;
    }
  }

  dspace_mem_id  = H5Screate_simple(2, hnquad, NULL);
  dspace_file_id = H5Screate_simple(2, hnquadg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, quad_offset, NULL, hnquad, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Quadrilaterals", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(pent); pent = NULL;

  dspace_mem_id  = H5Screate_simple(1, hnquad, NULL);
  dspace_file_id = H5Screate_simple(1, hnquadg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, quad_offset, NULL, hnquad, NULL);
  dset_id = H5Dcreate(grp_entities_id, "QuadrilateralsRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(pref); pref = NULL;

  dspace_mem_id  = H5Screate_simple(1, &nqreq, NULL);
  dspace_file_id = H5Screate_simple(1, &nqreqg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &required_offset[3], NULL, &nqreq, NULL);
  dset_id = H5Dcreate(grp_entities_id, "RequiredQuadrilaterals", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(preq); preq = NULL;

  dspace_mem_id  = H5Screate_simple(1, &nqpar, NULL);
  dspace_file_id = H5Screate_simple(1, &nqparg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &parallel_offset[3], NULL, &nqpar, NULL);
  dset_id = H5Dcreate(grp_entities_id, "ParallelQuadrilaterals", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(ppar); ppar = NULL;

  /* Tetrahedra */
  pent = (int*) calloc(ne, 4 * sizeof(int));
  pref = (int*) calloc(ne, sizeof(int));
  preq = (int*) calloc(nereq, sizeof(int));
  ppar = (int*) calloc(nepar, sizeof(int));

  reqcount = parcount = 0;

  if (ne) {
    ne = 0;
    for (int i = 0 ; i < mesh->ne ; i++) {
      pe = &mesh->tetra[i + 1];
      if (MG_EOK(pe)) {
        for (int j = 0 ; j < 4 ; j++) {
          pent[4 * ne + j] = mesh->point[pe->v[j]].tmp + point_offset[0] - 1;
        }
      }
      pref[i] = pe->ref;
      if (pe->tag & MG_REQ)    preq[reqcount++] = ne + tetra_offset[0];
      if (pe->tag & MG_PARBDY) ppar[parcount++] = ne; /* Local index for parallel entities */
      ne++;
    }
  }

  dspace_mem_id  = H5Screate_simple(2, hne, NULL);
  dspace_file_id = H5Screate_simple(2, hneg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, tetra_offset, NULL, hne, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Tetrahedra", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(pent); pent = NULL;

  dspace_mem_id  = H5Screate_simple(1, hne, NULL);
  dspace_file_id = H5Screate_simple(1, hneg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, tetra_offset, NULL, hne, NULL);
  dset_id = H5Dcreate(grp_entities_id, "TetrahedraRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(pref); pref = NULL;

  dspace_mem_id  = H5Screate_simple(1, &nereq, NULL);
  dspace_file_id = H5Screate_simple(1, &nereqg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &required_offset[4], NULL, &nereq, NULL);
  dset_id = H5Dcreate(grp_entities_id, "RequiredTetrahedra", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(preq); preq = NULL;

  dspace_mem_id  = H5Screate_simple(1, &nepar, NULL);
  dspace_file_id = H5Screate_simple(1, &neparg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &parallel_offset[4], NULL, &nepar, NULL);
  dset_id = H5Dcreate(grp_entities_id, "ParallelTetrahedra", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(ppar); ppar = NULL;

  /* Prisms */
  pent = (int*) calloc(nprism, 6 * sizeof(int));
  pref = (int*) calloc(nprism, sizeof(int));
  if (nprism){
    for (int i = 0 ; i < mesh->nprism ; i++) {
      pp = &mesh->prism[i + 1];
      for (int j = 0 ; j < 6 ; j++) {
        pent[6 * i + j] = mesh->point[pp->v[j]].tmp + point_offset[0] - 1;
      }
      pref[i] = pp->ref;
    }
  }

  dspace_mem_id  = H5Screate_simple(2, hnprism, NULL);
  dspace_file_id = H5Screate_simple(2, hnprismg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, prism_offset, NULL, hnprism, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Prisms", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(pent); pent = NULL;

  dspace_mem_id  = H5Screate_simple(1, hnprism, NULL);
  dspace_file_id = H5Screate_simple(1, hnprismg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, prism_offset, NULL, hnprism, NULL);
  dset_id = H5Dcreate(grp_entities_id, "PrismsRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  free(pref); pref = NULL;

  /* Close the groups */
  status = H5Gclose(grp_entities_id);
  status = H5Gclose(grp_mesh_id);

  /*------------------------- WRITE METRIC AND SOLUTIONS -------------------------*/

  if (met->size != 1 && met->size != 6) {
    fprintf(stderr, "\n  ## Error: %s: Wrong metric size\n", __func__);
    return 0;
  }
  if (np != met->np) {
    fprintf(stderr, "\n  ## Error: %s: The metric vertices do not match with the mesh vertices \n", __func__);
    return 0;
  }

  /* Open the group */
  grp_sols_id = H5Gcreate(file_id, "Solutions", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /* Arrays for bidimensional dataspaces */
  hsize_t hns[2]  = {np, met->size};
  hsize_t hnsg[2] = {npg, met->size};

  dspace_mem_id = H5Screate_simple(2, hns, NULL);
  dspace_file_id = H5Screate_simple(2, hnsg, NULL);
  H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, point_offset, NULL, hns, NULL);
  dset_id = H5Dcreate(grp_sols_id, "MetricAtVertices", H5T_NATIVE_DOUBLE, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dspace_mem_id, dspace_file_id, dxpl_id, &(met->m[1]));
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);

  /* Close the group */
  status = H5Gclose(grp_sols_id);

  /* Release the remaining HDF5 IDs (property lists and the file) */
  status = H5Fclose(file_id);
  status = H5Pclose(fapl_id);
  status = H5Pclose(dxpl_id);
  status = H5Pclose(dcpl_id);

  /*------------------------- HDF5 IOs END HERE -------------------------*/

  /*------------------------- WRITE LIGHT DATA IN XDMF FILE -------------------------*/

  if (rank == root) {
    FILE *xdmf_file = NULL;
    xdmf_file = fopen(xdmfname, "w");
    fprintf(xdmf_file, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
    fprintf(xdmf_file, "<Xdmf Version=\"3.0\">\n");
    fprintf(xdmf_file, "<Domain>\n");
    fprintf(xdmf_file, "    <Grid Name=\"3D Unstructured Mesh\" GridType=\"Uniform\">\n");
    fprintf(xdmf_file, "      <Topology TopologyType=\"Tetrahedron\" NumberOfElements=\"%llu\">\n", neg);
    fprintf(xdmf_file, "        <DataItem DataType=\"Int\"\n");
    fprintf(xdmf_file, "                  Format=\"HDF\"\n");
    fprintf(xdmf_file, "                  Dimensions=\"%llu 4\">\n", neg);
    fprintf(xdmf_file, "          %s:/Mesh/MeshEntities/Tetrahedra\n", filename);
    fprintf(xdmf_file, "        </DataItem>\n");
    fprintf(xdmf_file, "      </Topology>\n");
    fprintf(xdmf_file, "      <Geometry GeometryType=\"XYZ\">\n");
    fprintf(xdmf_file, "        <DataItem DataType=\"Float\"\n");
    fprintf(xdmf_file, "                  Precision=\"8\"\n");
    fprintf(xdmf_file, "                  Format=\"HDF\"\n");
    fprintf(xdmf_file, "                  Dimensions=\"%llu 3\">\n", npg);
    fprintf(xdmf_file, "          %s:/Mesh/MeshEntities/Vertices\n", filename);
    fprintf(xdmf_file, "        </DataItem>\n");
    fprintf(xdmf_file, "      </Geometry>\n");
    if (met) {
      if (met->size == 6)
        fprintf(xdmf_file, "      <Attribute Center=\"Node\" Name=\"Metric\" AttributeType=\"Tensor6\">\n");
      else if (met->size == 1)
        fprintf(xdmf_file, "      <Attribute Center=\"Node\" Name=\"Metric\" AttributeType=\"Scalar\">\n");
      fprintf(xdmf_file, "        <DataItem DataType=\"Float\"\n");
      fprintf(xdmf_file, "                  Precision=\"8\"\n");
      fprintf(xdmf_file, "                  Format=\"HDF\"\n");
      fprintf(xdmf_file, "                  Dimensions=\"%lld %d\">\n", npg, grp->met->size);
      fprintf(xdmf_file, "          %s:/Solutions/MetricAtVertices\n", filename);
      fprintf(xdmf_file, "        </DataItem>\n");
      fprintf(xdmf_file, "      </Attribute>\n");
    }
    /* for (int i = 0 ; i < nsols ; i++) { */
    /*   if (sols->type == MMG5_Scalar) { */
    /*     fprintf(xdmf_file, "      <Attribute Center=\"Node\" Name=\"Sol%d\" AttributeType=\"Scalar\">\n", i); */
    /*   } */
    /*   else if (sols->type == MMG5_Vector) { */
    /*     fprintf(xdmf_file, "      <Attribute Center=\"Node\" Name=\"Sol%d\" AttributeType=\"Vector\">\n", i); */
    /*   } */
    /*   else if (sols->type == MMG5_Tensor) { */
    /*     fprintf(xdmf_file, "      <Attribute Center=\"Node\" Name=\"Sol%d\" AttributeType=\"Tensor\">\n", i); */
    /*   } */
    /*   fprintf(xdmf_file, "        <DataItem DataType=\"Float\"\n"); */
    /*   fprintf(xdmf_file, "                  Precision=\"8\"\n"); */
    /*   fprintf(xdmf_file, "                  Format=\"HDF\"\n"); */
    /*   fprintf(xdmf_file, "                  Dimensions=\"%d %d\">\n", npg, sols[i].size); */
    /*   fprintf(xdmf_file, "          %s:/sols_grp/SolAtVertices%d\n", fileout, i); */
    /*   fprintf(xdmf_file, "        </DataItem>\n"); */
    /*   fprintf(xdmf_file, "      </Attribute>\n"); */
    /* } */
    fprintf(xdmf_file, "    </Grid>\n");
    fprintf(xdmf_file, "  </Domain>\n");
    fprintf(xdmf_file, "</Xdmf>\n");
    fclose(xdmf_file);
  }

  /*------------------------- END -------------------------*/

  return ier;
}


int PMMG_loadParmesh_hdf5(PMMG_pParMesh parmesh, const char *filename) {

  return 1;
}

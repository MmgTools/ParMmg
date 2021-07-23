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

static int PMMG_countEntities(PMMG_pParMesh parmesh, int ntyp_entities, hsize_t *nentities, hsize_t *nentitiesl, hsize_t* nentitiesg) {
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
  MPI_Comm comm = parmesh->comm;
  int rank, root, nprocs;

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
  if (!nentities) {
    fprintf(stderr, "\n  ## Error: %s: nentities array not allocated.\n",
            __func__);
    return 0;
  }
  if (!nentitiesl) {
    fprintf(stderr, "\n  ## Error: %s: nentitiesl array not allocated.\n",
            __func__);
    return 0;
  }
  if (!nentitiesg) {
    fprintf(stderr, "\n  ## Error: %s: nentitiesg array not allocated.\n",
            __func__);
    return 0;
  }

  /* Check mesh */
  if ( !mesh->point ) {
    fprintf(stderr, "\n  ## Error: %s: points array not allocated.\n",
            __func__);
    return 0;
  }

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
      /* if ((!ppt->flag) || MG_SIN(ppt->tag)) continue; */
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
      nprism++;
    }
  }

  /* Count global entities */
  nentities[ntyp_entities * rank + PMMG_saveVertex]  = np;
  nentities[ntyp_entities * rank + PMMG_saveEdge]    = na;
  nentities[ntyp_entities * rank + PMMG_saveTria]    = nt;
  nentities[ntyp_entities * rank + PMMG_saveQuad]    = nquad;
  nentities[ntyp_entities * rank + PMMG_saveTetra]   = ne;
  nentities[ntyp_entities * rank + PMMG_savePrism]   = nprism;
  nentities[ntyp_entities * rank + PMMG_saveCorner]  = nc;
  nentities[ntyp_entities * rank + PMMG_saveReq]     = npreq;
  nentities[ntyp_entities * rank + PMMG_savePar]     = nppar;
  nentities[ntyp_entities * rank + PMMG_saveRidge]   = nr;
  nentities[ntyp_entities * rank + PMMG_saveEdReq]   = nedreq;
  nentities[ntyp_entities * rank + PMMG_saveEdPar]   = nedpar;
  nentities[ntyp_entities * rank + PMMG_saveTriaReq] = ntreq;
  nentities[ntyp_entities * rank + PMMG_saveTriaPar] = ntpar;
  nentities[ntyp_entities * rank + PMMG_saveQuadReq] = nqreq;
  nentities[ntyp_entities * rank + PMMG_saveQuadPar] = nqpar;
  nentities[ntyp_entities * rank + PMMG_saveTetReq]  = nereq;
  nentities[ntyp_entities * rank + PMMG_saveTetPar]  = nepar;
  nentities[ntyp_entities * rank + PMMG_saveNormal]  = nnor;
  nentities[ntyp_entities * rank + PMMG_saveTangent] = ntan;

  for (int typent = 0 ; typent < ntyp_entities ; typent++) {
    nentitiesl[typent] = nentities[ntyp_entities * rank + typent];
  }

  MPI_Allgather(&nentities[ntyp_entities * rank], ntyp_entities, MPI_UNSIGNED_LONG_LONG,
                nentities                       , ntyp_entities, MPI_UNSIGNED_LONG_LONG, comm);

  for (int k = 0 ; k < nprocs ; k++) {
    for (int typent = 0 ; typent < ntyp_entities ; typent++) {
      nentitiesg[typent] += nentities[ntyp_entities * k + typent];
    }
  }

  return 1;
}

static inline int PMMG_computeHDFoffset(PMMG_pParMesh parmesh, int ntyp_entities, hsize_t *nentities, hsize_t *offset) {
  for (int k = 0 ; k < parmesh->myrank ; k++) {
    for (int typent = 0 ; typent < ntyp_entities ; typent++) {
      offset[2 * typent] += nentities[ntyp_entities * k + typent];
    }
  }
  return 1;
}

static int PMMG_saveHeader_hdf5(PMMG_pParMesh parmesh, hid_t file_id) {
  MMG5_pMesh mesh;
  hid_t      dspace_id;
  hid_t      attr_id;
  int        rank, root;
  herr_t     status;

  mesh = parmesh->listgrp[0].mesh;
  rank = parmesh->myrank;
  root = parmesh->info.root;

  dspace_id = H5Screate(H5S_SCALAR);
  attr_id = H5Acreate(file_id, "MeshVersionFormatted", H5T_NATIVE_INT, dspace_id, H5P_DEFAULT, H5P_DEFAULT);
  if (rank == root)
    status = H5Awrite(attr_id, H5T_NATIVE_INT, &mesh->ver);
  H5Aclose(attr_id);
  attr_id = H5Acreate(file_id, "Dimension", H5T_NATIVE_INT, dspace_id, H5P_DEFAULT, H5P_DEFAULT);
  if (rank == root)
    status = H5Awrite(attr_id, H5T_NATIVE_INT, &mesh->dim);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  return 1;
}

static int PMMG_saveMeshEntities_hdf5(PMMG_pParMesh parmesh, hid_t grp_entities_id, hid_t dcpl_id, hid_t dxpl_id,
                                      hsize_t *nentitiesl, hsize_t *nentitiesg, hsize_t *offset) {
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
  hsize_t nc, npreq, nppar;                      /* Corners, required and parallel vertices */
  hsize_t nr, nedreq, nedpar;                  /* Ridges, required and parallel edges */
  hsize_t ntreq, ntpar;                        /* Required and parallel triangles */
  hsize_t nqreq, nqpar;                        /* Required and parallel quads */
  hsize_t nereq, nepar;                        /* Required and parallel tetra */
  hsize_t nnor, ntan;                          /* Normals and Tangents */
  /* Global mesh size */
  hsize_t neg, npg, ntg, nag, nquadg, nprismg; /* Tetra, points, triangles, edges, quads, prisms */
  hsize_t ncg, npreqg, npparg;                   /* Corners, required and parallel vertices */
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

  /* MPI variables */
  int rank, root, nprocs;

  /* HDF5 variables */
  hid_t dspace_mem_id, dspace_file_id;
  hid_t dset_id;
  herr_t status;

  /*------------------------- INIT -------------------------*/

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

  /* Get the number of entities */
  np     = nentitiesl[PMMG_saveVertex];
  na     = nentitiesl[PMMG_saveEdge];
  nt     = nentitiesl[PMMG_saveTria];
  nquad  = nentitiesl[PMMG_saveQuad];
  ne     = nentitiesl[PMMG_saveTetra];
  nprism = nentitiesl[PMMG_savePrism];
  nc     = nentitiesl[PMMG_saveCorner];
  npreq   = nentitiesl[PMMG_saveReq];
  nppar   = nentitiesl[PMMG_savePar];
  nr     = nentitiesl[PMMG_saveRidge];
  nedreq = nentitiesl[PMMG_saveEdReq];
  nedpar = nentitiesl[PMMG_saveEdPar];
  ntreq  = nentitiesl[PMMG_saveTriaReq];
  ntpar  = nentitiesl[PMMG_saveTriaPar];
  nqreq  = nentitiesl[PMMG_saveQuadReq];
  nqpar  = nentitiesl[PMMG_saveQuadPar];
  nereq  = nentitiesl[PMMG_saveTetReq];
  nepar  = nentitiesl[PMMG_saveTetPar];
  nnor   = nentitiesl[PMMG_saveNormal];
  ntan   = nentitiesl[PMMG_saveTangent];

  npg     = nentitiesg[PMMG_saveVertex];
  nag     = nentitiesg[PMMG_saveEdge];
  ntg     = nentitiesg[PMMG_saveTria];
  nquadg  = nentitiesg[PMMG_saveQuad];
  neg     = nentitiesg[PMMG_saveTetra];
  nprismg = nentitiesg[PMMG_savePrism];
  ncg     = nentitiesg[PMMG_saveCorner];
  npreqg   = nentitiesg[PMMG_saveReq];
  npparg   = nentitiesg[PMMG_savePar];
  nrg     = nentitiesg[PMMG_saveRidge];
  nedreqg = nentitiesg[PMMG_saveEdReq];
  nedparg = nentitiesg[PMMG_saveEdPar];
  ntreqg  = nentitiesg[PMMG_saveTetReq];
  ntparg  = nentitiesg[PMMG_saveTetPar];
  nqreqg  = nentitiesg[PMMG_saveQuadReq];
  nqparg  = nentitiesg[PMMG_saveQuadPar];
  nereqg  = nentitiesg[PMMG_saveTetReq];
  neparg  = nentitiesg[PMMG_saveTetPar];
  nnorg   = nentitiesg[PMMG_saveNormal];
  ntang   = nentitiesg[PMMG_saveTangent];

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
  PMMG_MALLOC(parmesh, ppoint, 3 * np, double, "ppoint", return 0);
  PMMG_MALLOC(parmesh, pref, np, int, "pref", return 0);
  PMMG_MALLOC(parmesh, pcr, nc, int, "pcr", return 0);
  PMMG_MALLOC(parmesh, preq, npreq, int, "preq", return 0);
  PMMG_MALLOC(parmesh, ppar, nppar, int, "ppar", return 0);
  PMMG_MALLOC(parmesh, pnor, 3 * nnor, double, "pnor", return 0);
  PMMG_MALLOC(parmesh, ptan, 3 * ntan, double, "ptan", return 0);
  PMMG_MALLOC(parmesh, pnorat, nnor, int, "pnorat", return 0);
  PMMG_MALLOC(parmesh, ptanat, ntan, int, "ptanat", return 0);

  crcount = reqcount = parcount = ncount = tcount = 0;

  for (int i = 0 ; i < mesh->np ; i++) {
    ppt = &mesh->point[i + 1];
    if (MG_VOK(ppt)){
      for (int j = 0 ; j < 3 ; j++) {
        ppoint[3 * (ppt->tmp - 1) + j] = ppt->c[j];
      }
      if (ppt->tag & MG_CRN) {
        pcr[crcount++] = ppt->tmp + offset[2 * PMMG_saveVertex] - 1;
      }
      if (ppt->tag & MG_REQ) {
        preq[reqcount++] = ppt->tmp + offset[2 * PMMG_saveVertex] - 1;
      }
      if (ppt->tag & MG_PARBDY) {
        ppar[parcount++] = ppt->tmp - 1; /* Local index for parallel entities */
      }
      if (!MG_SIN(ppt->tag)) {
        /* Normals */
        if ((ppt->tag & MG_BDY) && (!(ppt->tag & MG_GEO) || (ppt->tag & MG_NOM))) {
          pxp = &mesh->xpoint[ppt->xp];
          for (int j = 0 ; j < 3 ; j++) {
            pnor[3 * ncount + j] = pxp->n1[j];
          }
          pnorat[ncount++] = ppt->tmp + offset[2 * PMMG_saveVertex] - 1;
        }
        /* Tangents */
        if (MG_EDG(ppt->tag) || (ppt->tag & MG_NOM)) {
          for (int j = 0 ; j < 3 ; j++) {
            ptan[3 * tcount + j] = ppt->n[j];
          }
          ptanat[tcount++] = ppt->tmp + offset[2 * PMMG_saveVertex] - 1;
        }
      }
      pref[ppt->tmp - 1] = abs(ppt->ref);
    }
  }

  dspace_mem_id  = H5Screate_simple(2, hnp, NULL);
  dspace_file_id = H5Screate_simple(2, hnpg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveVertex], NULL, hnp, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Vertices", H5T_NATIVE_DOUBLE, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dspace_mem_id, dspace_file_id, dxpl_id, ppoint);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, ppoint, double, "ppoint");

  dspace_mem_id  = H5Screate_simple(1, hnp, NULL);
  dspace_file_id = H5Screate_simple(1, hnpg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveVertex], NULL, hnp, NULL);
  dset_id = H5Dcreate(grp_entities_id, "VerticesRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, pref, int, "pref");

  dspace_mem_id  = H5Screate_simple(1, &nc, NULL);
  dspace_file_id = H5Screate_simple(1, &ncg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveCorner], NULL, &nc, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Corners", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pcr);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, pcr, int, "pcr");

  dspace_mem_id  = H5Screate_simple(1, &npreq, NULL);
  dspace_file_id = H5Screate_simple(1, &npreqg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveReq], NULL, &npreq, NULL);
  dset_id = H5Dcreate(grp_entities_id, "RequiredVertices", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, preq, int, "preq");

  dspace_mem_id  = H5Screate_simple(1, &nppar, NULL);
  dspace_file_id = H5Screate_simple(1, &npparg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_savePar], NULL, &nppar, NULL);
  dset_id = H5Dcreate(grp_entities_id, "ParallelVertices", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, ppar, int, "ppar");

  dspace_mem_id  = H5Screate_simple(2, hnnor, NULL);
  dspace_file_id = H5Screate_simple(2, hnnorg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveNormal], NULL, hnnor, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Normals", H5T_NATIVE_DOUBLE, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dspace_mem_id, dspace_file_id, dxpl_id, pnor);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, pnor, double, "pnor");

  dspace_mem_id  = H5Screate_simple(1, hnnor, NULL);
  dspace_file_id = H5Screate_simple(1, hnnorg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveNormal], NULL, hnnor, NULL);
  dset_id = H5Dcreate(grp_entities_id, "NormalsAtVertices", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pnorat);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, pnorat, int, "pnorat");

  dspace_mem_id  = H5Screate_simple(2, hntan, NULL);
  dspace_file_id = H5Screate_simple(2, hntang, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveTangent], NULL, hntan, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Tangents", H5T_NATIVE_DOUBLE, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dspace_mem_id, dspace_file_id, dxpl_id, ptan);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, ptan, double, "ptan");

  dspace_mem_id  = H5Screate_simple(1, hntan, NULL);
  dspace_file_id = H5Screate_simple(1, hntang, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveTangent], NULL, hntan, NULL);
  dset_id = H5Dcreate(grp_entities_id, "TangentsAtVertices", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ptanat);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, ptanat, int, "ptanat");

  /* Edges */
  PMMG_MALLOC(parmesh, pent, 2 * na, int, "pent", return 0);
  PMMG_MALLOC(parmesh, pref, na, int, "pref", return 0);
  PMMG_MALLOC(parmesh, pcr, nr, int, "pcr", return 0);
  PMMG_MALLOC(parmesh, preq, nedreq, int, "preq", return 0);
  PMMG_MALLOC(parmesh, ppar, nedpar, int, "ppar", return 0);

  crcount = reqcount = parcount = 0;

  if (na) {
    na = 0;
    for (int i = 0 ; i < mesh->na ; i++) {
      pa = &mesh->edge[i + 1];
      pent[2 * i]     = mesh->point[pa->a].tmp + offset[2 * PMMG_saveVertex] - 1;
      pent[2 * i + 1] = mesh->point[pa->b].tmp + offset[2 * PMMG_saveVertex] - 1;
      pref[i] = pa->ref;
      if (pa->tag & MG_GEO)    pcr[crcount++] = na + offset[2 * PMMG_saveEdge];
      if (pa->tag & MG_REQ)    preq[reqcount++] = na + offset[2 * PMMG_saveEdge];
      if (pa->tag & MG_PARBDY) ppar[parcount++] = na; /* Local index for parallel entities */
      na++;
    }
  }

  dspace_mem_id  = H5Screate_simple(2, hna, NULL);
  dspace_file_id = H5Screate_simple(2, hnag, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveEdge], NULL, hna, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Edges", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, pent, int, "pent");

  dspace_mem_id  = H5Screate_simple(1, hna, NULL);
  dspace_file_id = H5Screate_simple(1, hnag, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveEdge], NULL, hna, NULL);
  dset_id = H5Dcreate(grp_entities_id, "EdgesRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, pref, int, "pref");

  dspace_mem_id  = H5Screate_simple(1, &nr, NULL);
  dspace_file_id = H5Screate_simple(1, &nrg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveRidge], NULL, &nr, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Ridges", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pcr);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, pcr, int, "pcr");

  dspace_mem_id  = H5Screate_simple(1, &nedreq, NULL);
  dspace_file_id = H5Screate_simple(1, &nedreqg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveEdReq], NULL, &nedreq, NULL);
  dset_id = H5Dcreate(grp_entities_id, "RequiredEdges", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, preq, int, "preq");

  dspace_mem_id  = H5Screate_simple(1, &nedpar, NULL);
  dspace_file_id = H5Screate_simple(1, &nedparg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveEdPar], NULL, &nedpar, NULL);
  dset_id = H5Dcreate(grp_entities_id, "ParallelEdges", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, ppar, int, "ppar");

  /* Triangles */
  PMMG_MALLOC(parmesh, pent, 3 * nt, int, "pent", return 0);
  PMMG_MALLOC(parmesh, pref, nt, int, "pref", return 0);
  PMMG_MALLOC(parmesh, preq, ntreq, int, "preq", return 0);
  PMMG_MALLOC(parmesh, ppar, ntpar, int, "ppar", return 0);

  reqcount = parcount = 0;

  if (nt) {
    nt = 0;
    for (int i = 0 ; i < mesh->nt ; i++) {
      pt = &mesh->tria[i + 1];
      for (int j = 0 ; j < 3 ; j++) {
        pent[3 * i + j] = mesh->point[pt->v[j]].tmp + offset[2 * PMMG_saveVertex] - 1;
      }
      pref[i] = pt->ref;
      if (pt->tag[0] & MG_REQ && pt->tag[1] & MG_REQ && pt->tag[2] & MG_REQ) {
        preq[reqcount++] = nt + offset[2 * PMMG_saveTria];
      }
      if (pt->tag[0] & MG_PARBDY && pt->tag[1] & MG_PARBDY && pt->tag[2] & MG_PARBDY) {
        ppar[parcount++] = nt; /* Local index for parallel entities */
      }
      nt++;
    }
  }

  dspace_mem_id  = H5Screate_simple(2, hnt, NULL);
  dspace_file_id = H5Screate_simple(2, hntg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveTria], NULL, hnt, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Triangles", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, pent, int, "pent");

  dspace_mem_id  = H5Screate_simple(1, hnt, NULL);
  dspace_file_id = H5Screate_simple(1, hntg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveTria], NULL, hnt, NULL);
  dset_id = H5Dcreate(grp_entities_id, "TrianglesRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, pref, int, "pref");

  dspace_mem_id  = H5Screate_simple(1, &ntreq, NULL);
  dspace_file_id = H5Screate_simple(1, &ntreqg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveTriaReq], NULL, &ntreq, NULL);
  dset_id = H5Dcreate(grp_entities_id, "RequiredTriangles", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, preq, int, "preq");

  dspace_mem_id  = H5Screate_simple(1, &ntpar, NULL);
  dspace_file_id = H5Screate_simple(1, &ntparg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveTriaPar], NULL, &ntpar, NULL);
  dset_id = H5Dcreate(grp_entities_id, "ParallelTriangles", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, ppar, int, "ppar");

  /* Quadrilaterals */
  PMMG_MALLOC(parmesh, pent, 4 * nquad, int, "pent", return 0);
  PMMG_MALLOC(parmesh, pref, nquad, int, "pref", return 0);
  PMMG_MALLOC(parmesh, preq, nqreq, int, "preq", return 0);
  PMMG_MALLOC(parmesh, ppar, nqpar, int, "ppar", return 0);

  reqcount = parcount = 0;

  if (nquad){
    nquad = 0;
    for (int i = 0 ; i < mesh->nquad ; i++) {
      pq = &mesh->quadra[i + 1];
      for (int j = 0 ; j < 4 ; j++) {
        pent[4 * i + j] = mesh->point[pq->v[j]].tmp + offset[2 * PMMG_saveVertex] - 1;
      }
      pref[i] = pq->ref;
      if (pq->tag[0] & MG_REQ && pq->tag[1] & MG_REQ &&
          pq->tag[2] & MG_REQ && pq->tag[3] & MG_REQ) {
        preq[reqcount++] = nquad + offset[2 * PMMG_saveQuad];
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
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveQuad], NULL, hnquad, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Quadrilaterals", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, pent, int, "pent");

  dspace_mem_id  = H5Screate_simple(1, hnquad, NULL);
  dspace_file_id = H5Screate_simple(1, hnquadg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveQuad], NULL, hnquad, NULL);
  dset_id = H5Dcreate(grp_entities_id, "QuadrilateralsRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, pref, int, "pref");

  dspace_mem_id  = H5Screate_simple(1, &nqreq, NULL);
  dspace_file_id = H5Screate_simple(1, &nqreqg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveQuadReq], NULL, &nqreq, NULL);
  dset_id = H5Dcreate(grp_entities_id, "RequiredQuadrilaterals", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, preq, int, "preq");

  dspace_mem_id  = H5Screate_simple(1, &nqpar, NULL);
  dspace_file_id = H5Screate_simple(1, &nqparg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveQuadPar], NULL, &nqpar, NULL);
  dset_id = H5Dcreate(grp_entities_id, "ParallelQuadrilaterals", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, ppar, int, "ppar");

  /* Tetrahedra */
  PMMG_MALLOC(parmesh, pent, 4 * ne, int, "pent", return 0);
  PMMG_MALLOC(parmesh, pref, ne, int, "pref", return 0);
  PMMG_MALLOC(parmesh, preq, nereq, int, "preq", return 0);
  PMMG_MALLOC(parmesh, ppar, nepar, int, "ppar", return 0);

  reqcount = parcount = 0;

  if (ne) {
    ne = 0;
    for (int i = 0 ; i < mesh->ne ; i++) {
      pe = &mesh->tetra[i + 1];
      if (MG_EOK(pe)) {
        for (int j = 0 ; j < 4 ; j++) {
          pent[4 * ne + j] = mesh->point[pe->v[j]].tmp + offset[2 * PMMG_saveVertex] - 1;
        }
      }
      pref[i] = pe->ref;
      if (pe->tag & MG_REQ)    preq[reqcount++] = ne + offset[2 * PMMG_saveTetra];
      if (pe->tag & MG_PARBDY) ppar[parcount++] = ne; /* Local index for parallel entities */
      ne++;
    }
  }

  dspace_mem_id  = H5Screate_simple(2, hne, NULL);
  dspace_file_id = H5Screate_simple(2, hneg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveTetra], NULL, hne, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Tetrahedra", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, pent, int, "pent");

  dspace_mem_id  = H5Screate_simple(1, hne, NULL);
  dspace_file_id = H5Screate_simple(1, hneg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveTetra], NULL, hne, NULL);
  dset_id = H5Dcreate(grp_entities_id, "TetrahedraRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, pref, int, "pref");

  dspace_mem_id  = H5Screate_simple(1, &nereq, NULL);
  dspace_file_id = H5Screate_simple(1, &nereqg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveTetReq], NULL, &nereq, NULL);
  dset_id = H5Dcreate(grp_entities_id, "RequiredTetrahedra", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, preq);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, preq, int, "preq");

  dspace_mem_id  = H5Screate_simple(1, &nepar, NULL);
  dspace_file_id = H5Screate_simple(1, &neparg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_saveTetPar], NULL, &nepar, NULL);
  dset_id = H5Dcreate(grp_entities_id, "ParallelTetrahedra", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, ppar);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, ppar, int, "ppar");

  /* Prisms */
  PMMG_MALLOC(parmesh, pent, 6 * nprism, int, "pent", return 0);
  PMMG_MALLOC(parmesh, pref, nprism, int, "pref", return 0);

  if (nprism){
    for (int i = 0 ; i < mesh->nprism ; i++) {
      pp = &mesh->prism[i + 1];
      for (int j = 0 ; j < 6 ; j++) {
        pent[6 * i + j] = mesh->point[pp->v[j]].tmp + offset[2 * PMMG_saveVertex] - 1;
      }
      pref[i] = pp->ref;
    }
  }

  dspace_mem_id  = H5Screate_simple(2, hnprism, NULL);
  dspace_file_id = H5Screate_simple(2, hnprismg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_savePrism], NULL, hnprism, NULL);
  dset_id = H5Dcreate(grp_entities_id, "Prisms", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pent);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, pent, int, "pent");

  dspace_mem_id  = H5Screate_simple(1, hnprism, NULL);
  dspace_file_id = H5Screate_simple(1, hnprismg, NULL);
  status = H5Sselect_hyperslab(dspace_file_id, H5S_SELECT_SET, &offset[2 * PMMG_savePrism], NULL, hnprism, NULL);
  dset_id = H5Dcreate(grp_entities_id, "PrismsRef", H5T_NATIVE_INT, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, dspace_mem_id, dspace_file_id, dxpl_id, pref);
  H5Dclose(dset_id);
  H5Sclose(dspace_mem_id);
  H5Sclose(dspace_file_id);
  PMMG_DEL_MEM(parmesh, pref, int, "pref");

  return 1;
}

static int PMMG_saveCommunicators_hdf5(PMMG_pParMesh parmesh, hid_t grp_comm_id, hid_t dcpl_id, hid_t dxpl_id) {
  PMMG_pExt_comm comms;
  hsize_t        *ncomms, ncommg, comm_offset;
  int            *colors, *nface;
  MPI_Comm       comm;
  int            rank, nprocs, root;
  hid_t          dspace_mem_id, dspace_file_id;
  hid_t          dset_id;
  herr_t         status;

  /* Init variables */
  rank = parmesh->myrank;
  nprocs = parmesh->nprocs;
  root = parmesh->info.root;
  comm = parmesh->comm;
  comms = parmesh->ext_face_comm;

  ncommg = comm_offset = 0;

  /* Count the number of communicators */
  PMMG_MALLOC(parmesh, ncomms, nprocs, hsize_t, "ncomms", return 0);
  ncomms[rank] = parmesh->next_face_comm;
  MPI_Allgather(&ncomms[rank], 1, MPI_LONG_LONG, ncomms, 1, MPI_LONG_LONG, comm);

  for (int i = 0 ; i < nprocs ; i++) {
    ncommg += ncomms[i];
  }
  for (int i = 0 ; i < rank ; i++) {
    comm_offset += ncomms[i];
  }

  /* Create the buffers */
  PMMG_MALLOC(parmesh, colors, ncomms[rank], int, "colors", return 0);
  PMMG_MALLOC(parmesh, nface, ncomms[rank], int, "nface", return 0);
  for (int icomm = 0 ; icomm < ncomms[rank] ; icomm++) {
    colors[icomm] = comms[icomm].color_out;
    nface[icomm] = comms[icomm].nitem;
  }

  /* Write the things */
  /* Number of communicators */
  hsize_t hnprocs = nprocs;
  dspace_file_id = H5Screate_simple(1, &hnprocs, NULL);
  dset_id = H5Dcreate(grp_comm_id, "NumberOfFaceCommunicators", H5T_NATIVE_HSIZE, dspace_file_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  if (rank == root)
    status = H5Dwrite(dset_id, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ncomms[rank]);
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

  status = H5Sclose(dspace_mem_id);
  status = H5Sclose(dspace_file_id);

  PMMG_DEL_MEM(parmesh, ncomms, hsize_t, "ncomms");
  PMMG_DEL_MEM(parmesh, colors, int, "colors");
  PMMG_DEL_MEM(parmesh, nface, int, "nface");

  return 1;
}

static int PMMG_saveMetric_hdf5(PMMG_pParMesh parmesh, hid_t grp_sols_id, hid_t dcpl_id, hid_t dxpl_id,
                                hsize_t *nentitiesl, hsize_t *nentitiesg, hsize_t *offset) {
  int         np, npg, mcount;
  MMG5_pMesh  mesh;
  MMG5_pSol   met;
  MMG5_pPoint ppt;
  double      *met_buf;
  hsize_t     *met_offset;
  hid_t       dspace_mem_id, dspace_file_id;
  hid_t       dset_id;

  mesh = parmesh->listgrp[0].mesh;
  met = parmesh->listgrp[0].met;
  np  = nentitiesl[PMMG_saveVertex];
  npg = nentitiesg[PMMG_saveVertex];

  /* Check the metric */
  if (met->size != 1 && met->size != 6) {
    fprintf(stderr, "\n  ## Error: %s: Wrong metric size\n", __func__);
    return 0;
  }
  if (mesh->np != met->np) {
    fprintf(stderr, "\n  ## Error: %s: The metric vertices do not match with the mesh vertices \n", __func__);
    return 0;
  }
  if (!met->m) {
    fprintf(stderr, "\n  ## Warning: %s: No metric data to save \n", __func__);
    return 0;
  }

  /* Arrays for bidimensional dataspaces */
  hsize_t hns[2]  = {np, met->size};
  hsize_t hnsg[2] = {npg, met->size};

  PMMG_CALLOC(parmesh, met_offset, met->size, hsize_t, "met_offset", return 0);
  met_offset[0] = offset[2 * PMMG_saveVertex];

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
  PMMG_DEL_MEM(parmesh, met_offset, hsize_t, "met_offset");

  return 1;
}

static int PMMG_saveAllSols_hdf5(PMMG_pParMesh parmesh, hid_t grp_sols_id, hid_t dcpl_id, hid_t dxpl_id,
                                 hsize_t *nentitiesl, hsize_t *nentitiesg, hsize_t *offset) {
  int nsols, ndigits, np, npg, ne, neg, size, count, vcount, tcount;
  char        *solname;
  PMMG_pGrp   grp;
  MMG5_pMesh  mesh;
  MMG5_pSol   *sols;
  MMG5_pPoint ppt;
  MMG5_pTetra pt;
  hsize_t     *sol_offset;
  double      *sol_buf;
  hid_t       dspace_mem_id, dspace_file_id;
  hid_t       dset_id;

  /* Set ParMmg variables */
  grp = &parmesh->listgrp[0];
  mesh = grp->mesh;
  sols = &grp->field;
  nsols = mesh->nsols;

  /* Get the local and global number of vertices/tetra */
  np  = nentitiesl[PMMG_saveVertex];
  ne  = nentitiesl[PMMG_saveTetra];
  npg = nentitiesg[PMMG_saveVertex];
  neg = nentitiesg[PMMG_saveTetra];

  /* Arrays for bidimensional dataspaces */
  hsize_t hns[2]  = {0, 0};
  hsize_t hnsg[2] = {0, 0};

  /* Count digits for the name of the datasets */
  ndigits = PMMG_count_digits(nsols);

  /* Initialize the counts */
  vcount = tcount = 0;

  for (int i = 0 ; i < nsols ; i++) {
    size = sols[i]->size;
    count = 0;

    if (sols[i]->entities == MMG5_Noentity || sols[i]->entities == MMG5_Vertex) {
      PMMG_MALLOC(parmesh, sol_buf, size * np, double, "sol_buf", return 0);
      for (int k = 0 ; k < mesh->np ; k++) {
        ppt = &mesh->point[k + 1];
        if ( !MG_VOK(ppt)) continue;
        for (int j = 0 ; j < size ; j++) {
          sol_buf[count++] = sols[i]->m[1 + k * size + j];
        }
      }
      hns[0] = np; hns[1] = size;
      hnsg[0] = npg; hnsg[1] = size;
      PMMG_CALLOC(parmesh, sol_offset, np * size, hsize_t, "sol_offset", return 0);
      sol_offset[0] = offset[2 * PMMG_saveVertex];

      solname = (char*) malloc((strlen("SolAtVertices") + ndigits) * sizeof(char));
      strcpy(solname, "SolAtVertices");
      strcat(solname, (char*)&vcount);
      vcount++;
    }

    else if (sols[i]->entities == MMG5_Tetrahedron) {
      PMMG_MALLOC(parmesh, sol_buf, size * ne, double, "sol_buf", return 0);
      for (int k = 0 ; k < mesh->ne ; k++) {
        pt = &mesh->tetra[k + 1];
        if ( !MG_EOK(pt)) continue;
        for (int j = 0 ; j < size ; j++) {
          sol_buf[count++] = sols[i]->m[1 + k * size + j];
        }
      }
      hns[0] = ne; hns[1] = size;
      hnsg[0] = neg; hnsg[1] = size;
      PMMG_CALLOC(parmesh, sol_offset, ne * size, hsize_t, "sol_offset", return 0);
      sol_offset[0] = offset[2 * PMMG_saveTetra];

      solname = (char*) malloc((strlen("SolAtTetrahedra") + ndigits) * sizeof(char));
      strcpy(solname, "SolAtTetrahedra");
      strcat(solname, (char*)&tcount);
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

    PMMG_DEL_MEM(parmesh, sol_offset, hsize_t, "sol_offset");
    PMMG_DEL_MEM(parmesh, sol_buf, double, "sol_buf");
  }

  return 1;
}

static int PMMG_writeXDMF(PMMG_pParMesh parmesh, const char *filename, const char *xdmfname, hsize_t *nentitiesg) {
  hsize_t neg, npg;
  PMMG_pGrp grp;
  MMG5_pSol met, *sols;
  int nsols, rank, root, entities;

  npg  = nentitiesg[PMMG_saveVertex];
  neg  = nentitiesg[PMMG_saveTetra];
  grp  = &parmesh->listgrp[0];
  met  = grp->met;
  sols = &grp->field;
  nsols = grp->mesh->nsols;
  rank = parmesh->myrank;
  root = parmesh->info.root;

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
      fprintf(xdmf_file, "                  Dimensions=\"%lld %d\">\n", npg, met->size);
      fprintf(xdmf_file, "          %s:/Solutions/MetricAtVertices\n", filename);
      fprintf(xdmf_file, "        </DataItem>\n");
      fprintf(xdmf_file, "      </Attribute>\n");
    }
    for (int i = 0 ; i < nsols ; i++) {
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
    fprintf(xdmf_file, "    </Grid>\n");
    fprintf(xdmf_file, "  </Domain>\n");
    fprintf(xdmf_file, "</Xdmf>\n");
    fclose(xdmf_file);
  }

  return 1;
}

int PMMG_saveParmesh_hdf5(PMMG_pParMesh parmesh, const char *filename, const char *xdmfname) {
  int ier = 1;
  int ntyp_entities = 20;
  hsize_t *nentities, *nentitiesl, *nentitiesg; /* Number of entities (on each proc/on the current proc/global) */
  hsize_t *offset;                              /* Offset for the parallel writing with HDF5 */
  hid_t file_id, grp_mesh_id, grp_comm_id, grp_entities_id, grp_sols_id; /* HDF5 objects */
  hid_t fapl_id, dxpl_id, dcpl_id;                                       /* HDF5 property lists */
  herr_t status;
  MPI_Info info = MPI_INFO_NULL;
  MPI_Comm comm = parmesh->comm;
  int nprocs;

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

  /* Set all buffers to NULL */
  nentities = NULL;
  nentitiesl = NULL;
  nentitiesg = NULL;

  /* Set MPI variables */
  nprocs = parmesh->nprocs;

  PMMG_CALLOC(parmesh, nentities, ntyp_entities * nprocs, hsize_t, "nentities", return 0);
  PMMG_CALLOC(parmesh, nentitiesg, ntyp_entities, hsize_t, "nentitiesg", return 0);
  PMMG_CALLOC(parmesh, nentitiesl, ntyp_entities, hsize_t, "nentitiesl", return 0);

  /* Count the number of entities on each proc and globally */
  PMMG_countEntities(parmesh, ntyp_entities, nentities, nentitiesl, nentitiesg);

  /* Compute the offset for parallel writing */
  PMMG_CALLOC(parmesh, offset, 2 * ntyp_entities, hsize_t, "offset", return 0);
  PMMG_computeHDFoffset(parmesh, ntyp_entities, nentities, offset);

  /* Now that we have computed the offsets, each proc only needs to know its
     local number of entities and the global number of entities */
  PMMG_DEL_MEM(parmesh, nentities, hsize_t, "nentities");

  /*------------------------- HDF5 IOs START HERE -------------------------*/

  /* Shut HDF5 error stack */
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);

  /* Create the property lists */
  fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(fapl_id, comm, info);
  status = H5Pset_coll_metadata_write(fapl_id, 1);
  dxpl_id = H5Pcreate(H5P_DATASET_XFER);
  status = H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);
  dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
  status = H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER);

  /* Create the file */
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
  if (file_id < 0) {
    fprintf(stderr,"\n  ## Error: %s: Could not create the hdf5 file.\n",
            __func__);
    return 0;
  }

  /* Save the attributes (Version and Dimension) */
  PMMG_saveHeader_hdf5(parmesh, file_id);

  /* Open the mesh group */
  grp_mesh_id = H5Gcreate(file_id, "Mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (grp_mesh_id < 0) {
    fprintf(stderr,"\n  ## Error: %s: Could not create the mesh group.\n",
            __func__);
    return 0;
  }

  /* Write the communicators */
  grp_comm_id = H5Gcreate(grp_mesh_id, "FaceCommunicators", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (grp_comm_id < 0) {
    fprintf(stderr,"\n  ## Error: %s: Could not create the communicators group.\n",
            __func__);
    return 0;
  }
  PMMG_saveCommunicators_hdf5(parmesh, grp_comm_id, dcpl_id, dxpl_id);
  H5Gclose(grp_comm_id);

  /* Write the mesh entities */
  grp_entities_id = H5Gcreate(grp_mesh_id, "MeshEntities", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (grp_entities_id < 0) {
    fprintf(stderr,"\n  ## Error: %s: Could not create the mesh entities group.\n",
            __func__);
    return 0;
  }
  PMMG_saveMeshEntities_hdf5(parmesh, grp_entities_id, dcpl_id, dxpl_id, nentitiesl, nentitiesg, offset);
  H5Gclose(grp_entities_id);

  /* Close the mesh group */
  H5Gclose(grp_mesh_id);

  /* Write the metric and the solutions */
  grp_sols_id = H5Gcreate(file_id, "Solutions", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (grp_sols_id < 0) {
    fprintf(stderr,"\n  ## Error: %s: Could not create the solutions group.\n",
            __func__);
    return 0;
  }
  PMMG_saveMetric_hdf5(parmesh, grp_sols_id, dcpl_id, dxpl_id, nentitiesl, nentitiesg, offset);
  PMMG_saveAllSols_hdf5(parmesh, grp_sols_id, dcpl_id, dxpl_id, nentitiesl, nentitiesg, offset);
  H5Gclose(grp_sols_id);

  /* We no longer need the offset */
  PMMG_DEL_MEM(parmesh, offset, hsize_t, "offset");

  /* Release all HDF5 IDs */
  status = H5Fclose(file_id);
  status = H5Pclose(fapl_id);
  status = H5Pclose(dxpl_id);
  status = H5Pclose(dcpl_id);

  /* Write light data in XDMF file */
  if (!xdmfname || !*xdmfname)
    fprintf(stderr,"  ## Warning: %s: no XDMF file name provided.", __func__);
  else
    PMMG_writeXDMF(parmesh, filename, xdmfname, nentitiesg);

  /* We no longer need the number of entities */
  PMMG_DEL_MEM(parmesh, nentitiesl, hsize_t, "nentitiesl");
  PMMG_DEL_MEM(parmesh, nentitiesg, hsize_t, "nentitiesg");

  return ier;
}

static int PMMG_loadHeader_hdf5(PMMG_pParMesh parmesh, hid_t file_id) {
  MMG5_pMesh mesh;
  hid_t attr_id;

  mesh = parmesh->listgrp[0].mesh;

  attr_id = H5Aopen(file_id, "MeshVersionFormatted", H5P_DEFAULT);
  H5Aread(attr_id, H5T_NATIVE_INT, &mesh->ver);
  H5Aclose(attr_id);

  attr_id = H5Aopen(file_id, "Dimension", H5P_DEFAULT);
  H5Aread(attr_id, H5T_NATIVE_INT, &mesh->dim);
  H5Aclose(attr_id);

  if (mesh->dim != 3) {
    return 0;
  }

  return 1;
}

static int PMMG_loadCommunicators_hdf5(PMMG_pParMesh parmesh, hid_t grp_comm_id, hid_t dxpl_id, hsize_t *ncomm, int *color, int *nface) {
  hsize_t *ncommread, nprocsread;
  int *colorread, *part;
  int *nfaceread;
  int rank, nprocs, count;
  /* Metis variables */
  idx_t *xadj, *adjncy, *vwgt, *adjwgt;
  idx_t ncon = 1;
  hid_t dspace_mem_id, dspace_file_id;
  hid_t dset_id;

  ncommread = NULL; ncomm = NULL;
  colorread = NULL; color = NULL;
  nfaceread = NULL; nface = NULL;
  xadj = adjncy = vwgt = adjwgt = NULL;

  rank = parmesh->myrank;
  nprocs = parmesh->nprocs;

  /* Read the number of comms and the number of procs the mesh was saved with. */
  dset_id = H5Dopen(grp_comm_id, "NumberOfFaceCommunicators", H5P_DEFAULT);
  dspace_file_id = H5Dget_space(dset_id);
  H5Sget_simple_extent_dims(dspace_file_id, &nprocsread, NULL);
  H5Sclose(dspace_file_id);
  PMMG_MALLOC(parmesh, ncommread, nprocsread, hsize_t, "ncommread", return 0);
  H5Dread(dset_id, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL, dxpl_id, ncommread);
  H5Dclose(dset_id);

  /* If we try to read the mesh with less procs than the number of procs we saved
     the mesh with, we he have to merge some of the old partitions.
   */
  if (nprocs >= nprocsread) {
    PMMG_CALLOC(parmesh, xadj, nprocsread + 1, idx_t, "xadj", return 0);
    xadj[0] = 0;
    count = 0;
    for (int i = 0 ; i < nprocsread ; i++) {
      count += ncommread[i];
      xadj[i + 1] = count;
    }

    PMMG_CALLOC(parmesh, nfaceread, count, idx_t, "nfaceread", return 0);
    PMMG_CALLOC(parmesh, colorread, count, idx_t, "color", return 0);

    dset_id = H5Dopen(grp_comm_id, "Colors", H5P_DEFAULT);
    H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, dxpl_id, colorread);
    H5Dclose(dset_id);

    dset_id = H5Dopen(grp_comm_id, "NumberOfCommunicatorFaces", H5P_DEFAULT);
    H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, dxpl_id, nfaceread);
    H5Dclose(dset_id);

    /* Assign old partitions to each proc */
    PMMG_CALLOC(parmesh, part, nprocsread, idx_t, "part", return 0);
    METIS_PartGraphKway((idx_t*)&nprocsread, &ncon, xadj, colorread, NULL, NULL, NULL, &nprocs, NULL, NULL, NULL, NULL, part);

    PMMG_DEL_MEM(parmesh, xadj, idx_t, "xadj");
    PMMG_DEL_MEM(parmesh, nfaceread, idx_t, "nfaceread");
    PMMG_DEL_MEM(parmesh, colorread, idx_t, "colorread");
    PMMG_DEL_MEM(parmesh, part, idx_t, "part");
  }
  /* If we try to read the mesh with more procs than the number of procs we saved
     the mesh with, the excess procs do not work. */
  else {

  }

  return 1;
}

static int PMMG_loadMeshEntities_hdf5(PMMG_pParMesh parmesh, hid_t grp_entities_id, hid_t dxpl_id) {
  return 1;
}

static int PMMG_loadAllSols_hdf5(PMMG_pParMesh parmesh, hid_t grp_sols_id, hid_t dxpl_id) {
  return 1;
}

int PMMG_loadParmesh_hdf5(PMMG_pParMesh parmesh, const char *filename) {
  /* MMG variables */
  int ier = 1;
  int ntyp_entities = 20;
  hsize_t *nentities, *nentitiesl, *nentitiesg;

  /* Offset for parallel reading */
  hsize_t *offset;

  /* Communicators info */
  hsize_t *ncomm;
  int *color, *nface;

  /* HDF5 variables */
  hid_t file_id, grp_mesh_id, grp_comm_id, grp_entities_id, grp_sols_id; /* Objects */
  hid_t fapl_id, dxpl_id;                                                /* Property lists */
  herr_t status;

  /* MPI variables */
  MPI_Info info = MPI_INFO_NULL;
  MPI_Comm comm = parmesh->comm;
  int rank, root, nprocs;

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

  /* Set all buffers to NULL */
  nentities = NULL;
  nentitiesl = NULL;
  nentitiesg = NULL;
  ncomm = NULL;
  color = NULL;
  nface = NULL;
  /* Set MPI variables */
  nprocs = parmesh->nprocs;
  rank = parmesh->myrank;
  root = parmesh->info.root;

  /* Shut HDF5 error stack */
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);

  /* Create the property lists */
  fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  status = H5Pset_fapl_mpio(fapl_id, comm, info);    /* Parallel access to the file */
  status = H5Pset_all_coll_metadata_ops(fapl_id, 1); /* Collective metadata read */
  dxpl_id = H5Pcreate(H5P_DATASET_XFER);
  status = H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE); /* Collective dataset xfer operations */

  /* Open the HDF5 file */
  file_id = H5Fopen(filename, H5F_ACC_RDONLY, fapl_id);
  if (file_id < 0) {
    fprintf(stderr,"\n  ## Error: %s: Could not open the hdf5 file %s.\n",
            __func__, filename);
    return 0;
  }

  /* Load the header (version and dimension) */
  ier = PMMG_loadHeader_hdf5(parmesh, file_id);
  MPI_Allreduce(MPI_IN_PLACE, &ier, 1, MPI_INT, MPI_MIN, comm);
  if (ier == 0) {
    if (rank == root) {
      fprintf(stderr,"\n  ## Error: %s: Wrong mesh attributes in hdf5 file %s.\n",
              __func__, filename);
    }
    return 0;
  }

  /* Open the mesh group */
  grp_mesh_id = H5Gopen(file_id, "Mesh", H5P_DEFAULT);
  if (grp_mesh_id < 0) {
    fprintf(stderr,"\n  ## Error: %s: Could not open the /Mesh group in file %s.\n",
            __func__, filename);
    return 0;
  }

  /* Read the communicators and get the partitions */
  grp_comm_id = H5Gopen(grp_mesh_id, "FaceCommunicators", H5P_DEFAULT);
  if (grp_comm_id < 0) {
    fprintf(stderr,"\n  ## Error: %s: Could not open the /Mesh/FaceCommunicators group in file %s.\n",
            __func__, filename);
    return 0;
  }
  PMMG_loadCommunicators_hdf5(parmesh, grp_comm_id, dxpl_id, ncomm, color, nface);
  H5Gclose(grp_comm_id);

  /* Each proc reads the part of the mesh that is assigned to him */
  grp_entities_id = H5Gopen(grp_mesh_id, "MeshEntities", H5P_DEFAULT);
  if (grp_entities_id < 0) {
    fprintf(stderr,"\n  ## Error: %s: Could not open the /Mesh/MeshEntities group in file %s.\n",
            __func__, filename);
    return 0;
  }
  PMMG_loadMeshEntities_hdf5(parmesh, grp_entities_id, dxpl_id);
  H5Gclose(grp_entities_id);

  /* Close the mesh group */
  H5Gclose(grp_mesh_id);

  /* Each proc reads the part of the solutions/metric that is assigned to him */
  grp_sols_id = H5Gopen(file_id, "Solutions", H5P_DEFAULT);
  if (grp_sols_id < 0) {
    fprintf(stderr,"\n  ## Error: %s: Could not open the /Solutions group in file %s.\n",
            __func__, filename);
    return 0;
  }
  PMMG_loadAllSols_hdf5(parmesh, grp_sols_id, dxpl_id);
  H5Gclose(grp_sols_id);

  /*------------------------- RELEASE ALL HDF5 IDs -------------------------*/

  status = H5Fclose(file_id);
  status = H5Pclose(fapl_id);
  status = H5Pclose(dxpl_id);

  return 1;
}

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

#include "debug_pmmg.h"
#include <stdio.h> // fopen, fclose, fprintf
#include <stdlib.h> // abort
#ifdef __linux__
#include <malloc.h> // mallinfo
#endif

#define sw 4
#define sd 8

/**
 * \param name   filename to open
 * \param status file desriptor's desired mode
 *
 * \return file descriptor to newly opened file
 *
 * This is a fopen wrapper that checks the error code returned from fopen and
 * aborts on error. A file descriptor obtained from this function can safely
 * be used
 */
static FILE* PMMG_my_fopen( char *name, char *status )
{
   FILE *fp = fopen( name, status );
   if ( NULL == fp ) {
     fprintf( stderr, "error opening: %s. File: %s, line: %d \n", name, __FILE__, __LINE__);
     abort();
   }
   return fp;
}

/**
 * \param name filename to open
 * \param grp  pointer to list of mmg3d meshes
 * \param ngrp number of mmg3d meshes in grp
 *
 * This function writes all point and xpoint members of every mesh in grp to a
 * text file named "name"
 */
void PMMG_grplst_meshes_to_txt( char *name, PMMG_pGrp grp, int ngrp )
{
  int imsh,k;

  FILE *fp = PMMG_my_fopen( name, "w" );
  for ( imsh = 0; imsh < ngrp; ++imsh ) {
    fprintf( fp, "Points in mesh %d\n", imsh );
    for ( k = 0; k < grp[imsh].mesh->np + 2; k++ ) {
      fprintf( fp,
          "\tid#\t%10d\tcoords:\t(%9.6f,%9.6f,%9.6f )"
          "\tnormals:\t(%9.6f,%9.6f,%9.6f )"
          "\tref:\t%5d\t\txp\t%10d\ttmp\t%5d\ttag\t%5d",
          k,
          grp[ imsh ].mesh->point[ k ].c[ 0 ],
          grp[ imsh ].mesh->point[ k ].c[ 1 ],
          grp[ imsh ].mesh->point[ k ].c[ 2 ],
          grp[ imsh ].mesh->point[ k ].n[ 0 ],
          grp[ imsh ].mesh->point[ k ].n[ 1 ],
          grp[ imsh ].mesh->point[ k ].n[ 2 ],
          grp[ imsh ].mesh->point[ k ].ref,
          grp[ imsh ].mesh->point[ k ].xp,
          grp[ imsh ].mesh->point[ k ].tmp,
          grp[ imsh ].mesh->point[ k ].tag);
      if ( grp[ imsh ].mesh->point[ k ].xp != 0 )
        fprintf( fp, "\t\txp\t%10d\tn1\t(%9.6f,%9.6f,%9.6f )\tn2\t(%9.6f,%9.6f,%9.6f )",
            grp[ imsh ].mesh->point[ k ].xp,
            grp[ imsh ].mesh->xpoint[ grp[ imsh ].mesh->point[ k ].xp ].n1[ 0 ],
            grp[ imsh ].mesh->xpoint[ grp[ imsh ].mesh->point[ k ].xp ].n1[ 1 ],
            grp[ imsh ].mesh->xpoint[ grp[ imsh ].mesh->point[ k ].xp ].n1[ 2 ],
            grp[ imsh ].mesh->xpoint[ grp[ imsh ].mesh->point[ k ].xp ].n2[ 0 ],
            grp[ imsh ].mesh->xpoint[ grp[ imsh ].mesh->point[ k ].xp ].n2[ 1 ],
            grp[ imsh ].mesh->xpoint[ grp[ imsh ].mesh->point[ k ].xp ].n2[ 2 ]);
      fprintf( fp, "\n" );
    }
  }
  fclose(fp);
}

/**
 * \param name filename to open
 * \param mesh pointer to mmg3d mesh to write to file
 * \param num  number of tetras in mmg3d
 *
 * This function writes all the tetras' vertices of the given mmg3d mesh to a
 * text file named "name"
 */
void PMMG_tetras_of_mesh_to_txt( char *name, MMG5_pMesh mesh, int num )
{
  int k;

  FILE *fp = PMMG_my_fopen( name, "w" );
  fprintf( fp, "Tetras in  mesh %d.ne: %d, nei:%d\n", num, mesh->ne, mesh->nei );
  for ( k = 1; k < mesh->ne + 2; ++k )
    fprintf( fp,
        "%d %d %d %d\n",
        mesh->tetra->v[0],  mesh->tetra->v[1], mesh->tetra->v[2], mesh->tetra->v[3] );
  fclose(fp);
}

/**
 * \param name filename to open
 * \param grp  pointer to list of mmg3d meshes
 * \param nmsh number of mmg3d meshes in grp
 *
 * This function checks all the tetras in all meshes in grp for referencing
 * null vertice and outputs the ones detected to a text file named "name"
 */
void PMMG_find_tetras_referencing_null_points_to_txt( char *name,
                                                      PMMG_pGrp grp,
                                                      int nmsh )
{
  int imsh,tet,k;

  FILE *fp = PMMG_my_fopen( name, "w" );
  for ( imsh = 0; imsh < nmsh; ++imsh ) {
    for ( tet = 0; tet < grp[imsh].mesh->ne; ++tet ) {
      int check = 0;
      for ( k = 0; k < 4; ++k )
        if ( 0 == grp[imsh].mesh->tetra[tet].v[k] )
          ++check;
      if ( 3 == check )
        fprintf(fp, " mesh %d references point %d with all zero coordinates \n",imsh, tet );
    }
  }
  fclose(fp);
}

/**
 * \param mesh    pointer to mmg3d mesh
 * \param element tetra whose adjacent we want
 * \param face    face of tetra whose adjacent we want
 *
 * \return index in adjacency vector to adjacent tetra/face
 *
 * index in adjacency matrix of adjacent face for given tetra/face
 */
int PMMG_adja_idx_of_face( MMG5_pMesh mesh, int element, int face )
{
  int location = 4 * (element - 1) + 1 + face;
#ifndef NDEBUG
  int max_loc = 4 * (mesh->ne-1) + 5;
#endif

  assert( (face >= 0)    && (face < 4) && "There are only 4 faces per tetra" );
  assert( (location > 0) && (location < max_loc) && " adja out of bound "  );
  return mesh->adja[ location ];
}
/**
 * \param mesh    pointer to mmg3d mesh
 * \param element tetra whose adjacent tetra we are searching
 * \param face    face of tetra whose adjacent tetra we are searching
 *
 * \return base index of afjacent tetra to tetra/face
 *
 * find the idx in the adjacency vector of adjacent tetra to the given
 * tetra/face
 */
int PMMG_adja_tetra_to_face( MMG5_pMesh mesh, int element, int face )
{
  return PMMG_adja_idx_of_face( mesh, element, face ) / 4;
}
/**
 * \param mesh    pointer to mmg3d mesh
 * \param element tetra whose adjacent tetra we are searching
 * \param face    face of tetra whose adjacent tetra we are searching
 *
 * \return base index of afjacent tetra to tetra/face
 *
 * find the idx in the adjacency vector of the face of the adjacent tetra
 * to the given tetra/face
 */
int PMMG_adja_face_to_face( MMG5_pMesh mesh, int element, int face )
{
  return PMMG_adja_idx_of_face( mesh, element, face ) % 4;
}

/**
 * \param name filename to open
 * \param grp  pointer to list of mmg3d meshes
 * \param ngrp number of mmg3d meshes in grp
 *
 * create text file named "name" containing for every mesh's tetra/face
 * their adjacent tetra/face
 */
void PMMG_listgrp_meshes_adja_of_tetras_to_txt( char *name, PMMG_pGrp grp, int ngrp )
{
  int imsh,k,i;

  FILE *fp = PMMG_my_fopen( name, "w" );
  for ( imsh = 0; imsh < ngrp; ++imsh ) {
    fprintf( fp, "Mesh %d, ne= %d\n", imsh, grp[imsh].mesh->ne );
    for ( k = 1; k < grp[imsh].mesh->ne + 1; ++k ) {
      fprintf( fp, "tetra %d\t\t", k );
      for ( i = 0; i < 4; ++i ) {
        fprintf( fp, "adja[%d] %d", i, PMMG_adja_idx_of_face( grp[ imsh ].mesh, k, i ) );
        fprintf( fp, ", (tetra:%d, face:%1d)\t",
            PMMG_adja_tetra_to_face( grp[ imsh ].mesh, k, i ),
            PMMG_adja_face_to_face( grp[ imsh ].mesh, k, i ) );
      }
      fprintf( fp, "\n" );
    }
  }
  fclose(fp);
}

/**
 * \param mesh pointer toward the mesh structure
 * \param qual quality value.
 * \param inm pointer toward the solution file
 * \param bin 1 if binary file
 *
 * Write the quality value for a tetra in double precision.
 *
 */
void PMMG_writeDoubleQual3D(MMG5_pMesh mesh,double qual,FILE *inm,int bin) {

  if(!bin){
    fprintf(inm," %.15lg",qual);
  } else {
    fwrite((unsigned char*)&qual,sd,1,inm);
  }
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param filename name of file.
 * \param inm allocatable pointer toward the FILE structure.
 * \param ver file version (1=simple precision, 2=double).
 * \param bin 1 if the file is a binary.
 * \param np number of solutions of each type.
 * \param dim solution dimension.
 * \param nsols number of solutions of different types in the file.
 * \param type type of solutions.
 * \param size size of solutions.
 *
 * \return 0 if unable to open the file, 1 if success.
 *
 * Open the "filename" quality file and read the file header.
 *
 * \warning Binary version not tested.
 *
 */
int PMMG_saveQualHeader( MMG5_pMesh mesh,const char *filename,
                        FILE **inm,int ver,int *bin,int np,int dim,
                        int nsols,int *type,int *size) {
  MMG5_pTetra pt;
  int         binch,bpos;
  int         k;
  char        *ptr,*data,chaine[128];

  *bin = 0;

  MMG5_SAFE_CALLOC(data,strlen(filename)+6,char,return 0);
  strcpy(data,filename);
  ptr = strstr(data,".sol");
  if ( ptr ) {
    // filename contains the solution extension
    ptr = strstr(data,".solb");

    if ( ptr )  *bin = 1;

    if( !(*inm = fopen(data,"wb")) ) {
      fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
      MMG5_SAFE_FREE(data);
      return 0;
    }
  }
  else
  {
    // filename don't contains the solution extension
    ptr = strstr(data,".mesh");
    if ( ptr ) *ptr = '\0';

    strcat(data,".sol");
    if (!(*inm = fopen(data,"wb")) ) {
      ptr  = strstr(data,".solb");
      *ptr = '\0';
      strcat(data,".sol");
      if (!(*inm = fopen(data,"wb")) ) {
        fprintf(stderr,"  ** UNABLE TO OPEN %s.\n",data);
        MMG5_SAFE_FREE(data);
        return 0;
      }
      else *bin = 1;
    }
  }

  if ( mesh->info.imprim >= 0 )
    fprintf(stdout,"  %%%% %s OPENED\n",data);
  MMG5_SAFE_FREE(data);

  /*entete fichier*/
  binch=bpos=0;
  if(!*bin) {
    strcpy(&chaine[0],"MeshVersionFormatted\n");
    fprintf(*inm,"%s %d",chaine,ver);
    strcpy(&chaine[0],"\n\nDimension\n");
    fprintf(*inm,"%s %d",chaine,dim);
  } else {
    binch = 1; //MeshVersionFormatted
    fwrite(&binch,sw,1,*inm);
    binch = ver; //version
    fwrite(&binch,sw,1,*inm);
    binch = 3; //Dimension
    fwrite(&binch,sw,1,*inm);
    bpos = 20; //Pos
    fwrite(&bpos,sw,1,*inm);
    binch = dim;
    fwrite(&binch,sw,1,*inm);
  }

  np = 0;
  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( MG_EOK(pt) )  np++;
  }

  if(!*bin) {
    strcpy(&chaine[0],"\n\nSolAtTetrahedra\n");
    fprintf(*inm,"%s",chaine);
    fprintf(*inm,"%d\n",np);
    fprintf(*inm,"%d",nsols);
    for (k=0; k<nsols; ++k )
      fprintf(*inm," %d",type[k]);
    fprintf(*inm,"\n");
  } else {
    binch = 62; //Vertices
    fwrite(&binch,sw,1,*inm);
    bpos += 16;

    for (k=0; k<nsols; ++k )
      bpos += 4 + (size[k]*ver)*4*np; //Pos
    fwrite(&bpos,sw,1,*inm);

    fwrite(&np,sw,1,*inm);
    fwrite(&nsols,sw,1,*inm);
    for (k=0; k<nsols; ++k )
      fwrite(&type[k],sw,1,*inm);
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write isotropic or anisotropic metric.
 *
 */
int PMMG_saveMark(MMG5_pMesh mesh, const char *filename) {
  FILE*        inm;
  MMG5_pTetra  pt;
  int          ver,type,size;
  int          binch,bin,ier,k;

  ver = 2;
  type = MMG5_Scalar;
  size = 1;
  ier = PMMG_saveQualHeader( mesh,filename,&inm,ver,&bin,1,3,
                            1,&type,&size);

  if ( ier < 1 )  return ier;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    PMMG_writeDoubleQual3D(mesh,pt->mark,inm,bin);
    fprintf(inm,"\n");
  }

  /* End file */
  if(!bin) {
    fprintf(inm,"\n\nEnd\n");
  } else {
    binch = 54; //End
    fwrite(&binch,sw,1,inm);
  }
  fclose(inm);
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write isotropic or anisotropic metric.
 *
 */
int PMMG_saveQual(MMG5_pMesh mesh, const char *filename) {
  FILE*        inm;
  MMG5_pTetra  pt;
  int          ver,type,size;
  int          binch,bin,ier,k;

  ver = 2;
  type = MMG5_Scalar;
  size = 1;
  ier = PMMG_saveQualHeader( mesh,filename,&inm,ver,&bin,1,3,
                            1,&type,&size);

  if ( ier < 1 )  return ier;

  for (k=1; k<=mesh->ne; k++) {
    pt = &mesh->tetra[k];
    if ( !MG_EOK(pt) ) continue;

    PMMG_writeDoubleQual3D(mesh,pt->qual,inm,bin);
    fprintf(inm,"\n");
  }

  /* End file */
  if(!bin) {
    fprintf(inm,"\n\nEnd\n");
  } else {
    binch = 54; //End
    fwrite(&binch,sw,1,inm);
  }
  fclose(inm);
  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param grpId index of the group
 * \param basename filename prefix
 *
 * Write group mesh and metrics in medit format.
 *
 */
int PMMG_grp_to_saveMesh( PMMG_pParMesh parmesh, int grpId, char *basename ) {
  PMMG_pGrp  grp;
  MMG5_pMesh mesh;
  MMG5_pSol  met,field;
  char       name[ 2048 ];
  int        ier;

  grp  = &parmesh->listgrp[grpId];
  mesh = grp->mesh;
  met  = grp-> met;
  field= grp->field;

  assert( ( strlen( basename ) < 2048 - 14 ) && "filename too big" );
  sprintf( name, "%s-P%02d-%02d.mesh", basename, parmesh->myrank, grpId );

  PMMG_TRANSFER_AVMEM_FROM_PARMESH_TO_MESH(parmesh,mesh);

  /* Rebuild boundary */
  if ( !mesh->adja ) {
    ier = MMG3D_hashTetra( mesh, 0 );
  }
  MMG3D_bdryBuild( mesh ); //note: no error checking

  /* Save mesh */
  MMG3D_saveMesh( mesh, name );

  /* Destroy boundary */
  MMG5_DEL_MEM(mesh,mesh->tria);
  mesh->nt = 0;

  /* Save metrics */
  if ( met->m ) {
    sprintf( name, "%s-P%02d-%02d.sol", basename, parmesh->myrank, grpId );
    MMG3D_saveSol( mesh, met, name );
  }

  /* Save field */
  if ( mesh->nsols ) {
    assert ( field );
    sprintf( name, "%s-fields-P%02d-%02d.sol", basename, parmesh->myrank, grpId );
    MMG3D_saveAllSols( mesh, &field, name );
  }

  PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PARMESH(parmesh,mesh);

  return ier;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param grpId index of the group
 * \param basename filename prefix
 *
 * Write group mesh and tetra mark in medit format.
 *
 */
int PMMG_grp_mark_to_saveMesh( PMMG_pParMesh parmesh, int grpId, char *basename ) {
  PMMG_pGrp  grp;
  MMG5_pMesh mesh;
  char       name[ 2048 ];
  int        ier;

  grp  = &parmesh->listgrp[grpId];
  mesh = grp->mesh;

  assert( ( strlen( basename ) < 2048 - 14 ) && "filename too big" );
  sprintf( name, "%s-P%02d-%02d.mesh", basename, parmesh->myrank, grpId );

  PMMG_TRANSFER_AVMEM_FROM_PARMESH_TO_MESH(parmesh,mesh);
 
  ier = MMG3D_hashTetra( mesh, 0 );
  MMG3D_bdryBuild( mesh ); //note: no error checking
  /* Destroy boundary */
  MMG5_DEL_MEM(mesh,mesh->tria);
  mesh->nt = 0;

  MMG3D_saveMesh( mesh, name );
 
  sprintf( name, "%s-P%02d-%02d.sol", basename, parmesh->myrank, grpId );
  PMMG_saveMark( mesh, name );


  PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PARMESH(parmesh,mesh);

  return ier;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param grpId index of the group
 * \param basename filename prefix
 *
 * Write group mesh and quality in medit format.
 *
 */
int PMMG_grp_quality_to_saveMesh( PMMG_pParMesh parmesh, int grpId, char *basename ) {
  PMMG_pGrp  grp;
  MMG5_pMesh mesh;
  char       name[ 2048 ];
  int        ier;

  grp  = &parmesh->listgrp[grpId];
  mesh = grp->mesh;

  assert( ( strlen( basename ) < 2048 - 14 ) && "filename too big" );
  sprintf( name, "%s-P%02d-%02d.mesh", basename, parmesh->myrank, grpId );

  PMMG_TRANSFER_AVMEM_FROM_PARMESH_TO_MESH(parmesh,mesh);

  ier = MMG3D_hashTetra( mesh, 0 );
  MMG3D_bdryBuild( mesh ); //note: no error checking
  /* Destroy boundary */
  MMG5_DEL_MEM(mesh,mesh->tria);
  mesh->nt = 0;

  MMG3D_saveMesh( mesh, name );

  sprintf( name, "%s-P%02d-%02d.sol", basename, parmesh->myrank, grpId );
  PMMG_saveQual( mesh, name );


  PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PARMESH(parmesh,mesh);

  return ier;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param basename filenames prefix
 *
 * Write meshes and metrics of all groups in medit format.
 *
 */
int PMMG_listgrp_to_saveMesh( PMMG_pParMesh parmesh, char *basename ) {
  int grpId;

  for ( grpId = 0 ; grpId < parmesh->ngrp ; grpId++ )
    if( !PMMG_grp_to_saveMesh( parmesh, grpId, basename ) )
      return 0;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param basename filenames prefix
 *
 * Write meshes and tetra mark of all groups in medit format.
 *
 */
int PMMG_listgrp_mark_to_saveMesh( PMMG_pParMesh parmesh, char *basename ) {
  int grpId;

  for ( grpId = 0 ; grpId < parmesh->ngrp ; grpId++ )
    if( !PMMG_grp_mark_to_saveMesh( parmesh, grpId, basename ) )
      return 0;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param basename filenames prefix
 *
 * Write meshes and quality of all groups in medit format.
 *
 */
int PMMG_listgrp_quality_to_saveMesh( PMMG_pParMesh parmesh, char *basename ) {
  int grpId;

  for ( grpId = 0 ; grpId < parmesh->ngrp ; grpId++ )
    if( !PMMG_grp_quality_to_saveMesh( parmesh, grpId, basename ) )
      return 0;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param int_comm pointer toward the internal communicator (face or node one)
 * \param ext_comm pointer toward the external communicators (face or node ones)
 * \param next_comm number of external communicators
 *
 * For each item in the external communicators, print their index in the
 * internal communicator..
 *
 */
void PMMG_print_ext_comm( PMMG_pParMesh parmesh, PMMG_pInt_comm int_comm,
                          PMMG_pExt_comm ext_comm, int next_comm ) {
  PMMG_pExt_comm pext_comm;
  int iext_comm,iext;

  printf("myrank %d, int comm nitem %d\n",parmesh->myrank,int_comm->nitem);
  for( iext_comm = 0; iext_comm < next_comm; iext_comm++ ) {
    pext_comm = &ext_comm[iext_comm];
    for( iext = 0; iext < pext_comm->nitem; iext++ ) {
      printf("  myrank %d, ext comm %d, color out %d, item %d, idx %d\n",
             parmesh->myrank,iext_comm,pext_comm->color_out,iext,
             pext_comm->int_comm_index[iext]);
    }
  }

}

/**
 * \param msg custom message to print in the beggining of the report
 * \param id process number
 *
 * report memory usage as reported be glibc on linux systems
 */
void PMMG_dump_malloc_allocator_info( char *msg, int id )
{
  char name[ 16 ];
  FILE *fp;

  sprintf(name,"mem_info-%02d.txt", id );
  fp = PMMG_my_fopen( name, "a" );

#ifdef __linux__
  const int mb = 1024 * 1024;
  struct mallinfo me = mallinfo();

  fprintf( fp, "%s \n", msg );
  fprintf( fp, "** MALLOC ALLOCATOR INFO ***************************\n" );
  fprintf( fp, "* %4d \tNon-mmapped space allocated (mbytes)       *\n", me.arena / mb );
  fprintf( fp, "* %4d \tNumber of free chunks                      *\n", me.ordblks );
  fprintf( fp, "* %4d \tNumber of free fastbin blocks              *\n", me.smblks );
  fprintf( fp, "* %4d \tNumber of mmapped regions                  *\n", me.hblks );
  fprintf( fp, "* %4d \tSpace allocated in mmapped regions (mbytes)*\n", me.hblkhd / mb );
  fprintf( fp, "* %4d \tMaximum total allocated space (mbytes)     *\n", me.usmblks / mb );
  fprintf( fp, "* %4d \tSpace in freed fastbin blocks (mbytes)     *\n", me.fsmblks / mb );
  fprintf( fp, "* %4d \tTotal allocated space (mbytes)             *\n", me.uordblks / mb );
  fprintf( fp, "* %4d \tTotal free space (mbytes)                  *\n", me.fordblks / mb );
  fprintf( fp, "* %4d \tTop-most, releasable space (mbytes)        *\n", me.keepcost / mb );
  fprintf( fp, "****************************************************\n\n" );
#else
  fprintf( fp, "Extended information read directly from the malloc allocator is"
      "currently only implemented on linux\n" );
#endif

  fclose(fp);
}

/**
 * \param parmesh pointer to parmmg structure
 * \param msg     custom msg to include in the output messages
 *
 * check if any of the
 *   sum of memMax fields (in parmesh struct and in listgrp meshes)
 * or
 *   sum of memCur fields (in parmesh struct and in listgrp meshes)
 * exceed the memGloMax limit.
 */
void PMMG_check_mem_max_and_mem_cur( PMMG_pParMesh parmesh, const char *msg )
{
  size_t i,n_total = parmesh->memCur;
  const size_t mb = 1024 * 1024;

  for ( i = 0; i < parmesh->ngrp; ++i )
    n_total += parmesh->listgrp[ i ].mesh->memCur;

  if ( n_total > parmesh->memGloMax )
    fprintf( stderr,
             "%2d-%2d: %s: memCur check ERROR: memCur ( %8.2fMb ) > memGloMax ( %8.2fMb ) at %s %s %d\n",
             parmesh->myrank, parmesh->nprocs, msg,
             n_total / (float) mb, parmesh->memGloMax / (float) mb,
             __func__, __FILE__, __LINE__ );
//  else
//    fprintf( stderr,
//             "%2d-%2d: %s: memCur check OK: memCur = %8.2fMb - memGloMax = %8.2fMb \n",
//             parmesh->myrank, parmesh->nprocs, mesg,
//             n_total / (float) mb, parmesh->memGloMax / (float) mb );

  n_total = parmesh->memMax;
  for ( i = 0; i < parmesh->ngrp; ++i )
    n_total += parmesh->listgrp[ i ].mesh->memMax;
  if ( n_total > parmesh->memGloMax )
    fprintf( stderr,
             "%2d-%2d: %s: memMax check ERROR: memMax ( %8.2fMb ) > memGloMax ( %8.2fMb ) at %s %s %d\n",
             parmesh->myrank, parmesh->nprocs, msg,
             n_total / (float) mb, parmesh->memGloMax / (float) mb,
             __func__, __FILE__, __LINE__ );
//  else
//    fprintf( stderr,
//             "%2d-%2d: %s: memMax check OK: memMax = %8.2fMb - memGloMax = %8.2fMb \n",
//             parmesh->myrank, parmesh->nprocs, msg,
//             n_total / (float) mb, parmesh->memGloMax / (float) mb );
}

#include "debug_pmmg.h"
#include <stdio.h> // fopen, fclose, fprintf
#include <stdlib.h> // abort
#ifdef __linux__
#include <malloc.h> // mallinfo
#endif

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
  MMG5_pSol  met;
  size_t     memAv,oldMemMax;
  char       name[ 2048 ];
  int        ier;
  
  grp  = &parmesh->listgrp[grpId];
  mesh = grp->mesh;
  met  = grp-> met;

  assert( ( strlen( basename ) < 2048 - 14 ) && "filename too big" );
  sprintf( name, "%s-P%02d-%02d.mesh", basename, parmesh->myrank, grpId );
  
  oldMemMax = parmesh->memCur;
  memAv = parmesh->memMax-oldMemMax;
  PMMG_TRANSFER_AVMEM_FROM_PMESH_TO_MESH(parmesh,mesh,memAv,oldMemMax);
 
  ier = MMG3D_hashTetra( mesh, 1 );
  MMG3D_bdryBuild( mesh ); //note: no error checking
  MMG3D_saveMesh( mesh, name );
  if ( met->m ) {
    sprintf( name, "%s-P%02d-%02d.sol", basename, parmesh->myrank, grpId );
    MMG3D_saveSol( mesh, met, name );
  }

  PMMG_TRANSFER_AVMEM_FROM_MESH_TO_PMESH(parmesh,mesh,memAv,oldMemMax);

  return ier;
}


/**
 * \param parmesh filename to open
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

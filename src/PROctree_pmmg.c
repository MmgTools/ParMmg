/**
 * \file PROctree_pmmg.c
 * \brief Tools for local search around coordinates based on PROctree.
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \author Nikos Pattakos (Inria)
 * \author Luca Cirrottola (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */
#include "parmmg.h"

/**
 * \param mesh pointer toward the PROctree structure.
 * \param coord z-ordered coordinate in the octree
 * \param leaves double pointer to the leaf entities
 * \return 1 if ok 0 if fail
 *
 * Get the list of entities in the PROctree cell given a z-ordered coordinate.
 *
 */
int PMMG_getPROctree_leaves( MMG3D_pPROctree q,int64_t coord,int **leaves ) {
  MMG3D_PROctree_s *root,*root1;
  int64_t z;
  int nleaves;

  nleaves = 0;
  root = q->q0;
  z = coord;

  while( 1 )
    if( root->branches ) {
      /* Get the branch and right-shift z */
      root1 = &root->branches[ z % 8 ];
      if( root1->nbVer ) {
        root = root1;
        z = z >> 3;
      } else return 0;
    } else break;

  /* Get the leaves */
  *leaves = root->v;
  nleaves = root->nbVer;
  assert( *leaves );
  assert( nleaves );

  return nleaves;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param q pointer toward a pointer toward the global PROctree.
 *
 * Free the global PROctree structure.
 *
 */
void PMMG_freePROctrees(PMMG_pParMesh parmesh,MMG3D_pPROctree *q)
{
  int i;
  for( i = 0; i < parmesh->ngrp; i++ )
    MMG3D_freePROctree( parmesh->listgrp[i].mesh, &q[i] );
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward an PROctree cell.
 * \param ver vertex coordinates scaled such that the quadrant is [0;1]x[0;1]x[0;1]
 * \param no vertex index in the mesh.
 * \param nv maximum number of points in an PROctree cell.
 * \return 1 if ok 0 if memory saturated
 *
 * Add vertex in the suitable quadrant of the PROctree. This function is
 * recursively called until we reach the last one. At each step, the vertex
 * coordinates are scaled such as the quadrant is the [0;1]x[0;1]x[0;1] box.
 *
 */
int PMMG_addPROctreeRec(MMG5_pMesh mesh, MMG3D_PROctree_s* q, double* ver,
                        const int no, int nv)
{
  double   pt[3];
  int      dim, nbBitsInt,depthMax,i,j,k,ie,iloc;
  int      quadrant,sizBr;
  int      sizeRealloc;

  nbBitsInt = sizeof(int64_t)*8;
  dim       = mesh->dim;
  depthMax  = nbBitsInt/dim - 1; // maximum depth is to allow integer coordinates
  sizBr     = 1<<dim;

  if ( q->depth < depthMax ) // not at the maximum depth of the tree
  {
    if (q->nbVer < nv)  // not at the maximum number of vertice in the cell
    {

      if(q->nbVer == 0)  // first vertex list allocation
      {
        MMG5_ADD_MEM(mesh,sizeof(int),"PROctree vertice table", return 0);
        MMG5_SAFE_MALLOC(q->v,1,int,return 0);
      }
      else if(!(q->nbVer & (q->nbVer - 1))) //is a power of 2 -> reallocation of the vertex list
      {
        sizeRealloc = q->nbVer;
        sizeRealloc<<=1;
        MMG5_ADD_MEM(mesh,(sizeRealloc-sizeRealloc/2)*sizeof(int),"PROctree realloc",
                      return 0);
        MMG5_SAFE_REALLOC(q->v,q->nbVer,sizeRealloc,int,"PROctree",return 0);
      }

      q->v[q->nbVer] = no;
      q->nbVer++;
      return 1;
    }
    else if (q->nbVer == nv && q->branches==NULL)  //vertex list at maximum -> cell subdivision
    {
      /* creation of sub-branch and relocation of vertices in the sub-branches */
      MMG5_ADD_MEM(mesh,sizBr*sizeof(MMG3D_PROctree_s),"PROctree branches",
                    return 0);
      MMG5_SAFE_MALLOC(q->branches,sizBr,MMG3D_PROctree_s,return 0);

      for ( i = 0; i<sizBr; i++)
      {
        MMG3D_initPROctree_s(&(q->branches[i]));
        q->branches[i].depth = q->depth+1;
      }
      q->nbVer++;
      for (i = 0; i<nv; i++)
      {

        ie   = q->v[i] / 4;
        iloc = q->v[i] % 4;
        memcpy(&pt, mesh->point[mesh->tetra[ie].v[iloc]].c ,dim*sizeof(double));
        for ( j =0; j < q->depth; j++)
        {
          for (k = 0; k<dim; k++)
          {
            pt[k] -= ((double) (pt[k]>0.5))*0.5;
            pt[k] *= 2;
          }
        }
        if (!PMMG_addPROctreeRec(mesh, q, pt, q->v[i],nv))
          return 0;
        q->nbVer--;
      }
      if (!PMMG_addPROctreeRec(mesh, q, ver, no, nv))
        return 0;
      q->nbVer--;
      MMG5_DEL_MEM(mesh,q->v);

    }else // Recursive call in the corresponding sub cell
    {
      quadrant = 0;
      for ( i = 0; i<dim; i++)
      {
        quadrant += ((double) (ver[i]>0.5))*(1<<i);
        ver[i] -= ((double) (ver[i]>0.5))*0.5;
        ver[i] *= 2;
      }

      q->nbVer++;
      if (!PMMG_addPROctreeRec(mesh, &(q->branches[quadrant]), ver, no, nv))
        return 0;
    }
  }else // maximum PROctree depth reached
  {
    if (q->nbVer < nv)
    {
      if(q->nbVer == 0) // first allocation
      {
        MMG5_ADD_MEM(mesh,sizeof(int),"PROctree vertices table",
                      return 0);
        MMG5_SAFE_MALLOC(q->v,1,int,return 0);
      }
      else if(!(q->nbVer & (q->nbVer - 1))) //is a power of 2 -> normal reallocation
      {
        sizeRealloc = q->nbVer;
        sizeRealloc <<= 1;
        MMG5_ADD_MEM(mesh,(sizeRealloc-sizeRealloc/2)*sizeof(int),"PROctree realloc",
                      return 0);
        MMG5_SAFE_REALLOC(q->v,q->nbVer,sizeRealloc,int,"PROctree",return 0);
      }
    }
    else if (q->nbVer%nv == 0) // special reallocation of the vertex list because it is at maximum depth
    {
      MMG5_ADD_MEM(mesh,nv*sizeof(int),"PROctree realloc",
                    return 0);
      MMG5_SAFE_REALLOC(q->v,q->nbVer,q->nbVer+nv,int,"PROctree",return 0);
    }

    q->v[q->nbVer] = no;
    q->nbVer++;
  }

  return 1;
}


/**
 * \param mesh pointer toward the mesh structure.
 * \param q pointer toward the global PROctree
 * \param nv maximum number of vertices in each cell before subdivision
 * \return 1 if ok 0 if memory saturated
 *
 * Initialisation of the PROctree cell.
 *
 */
int PMMG_initPROctree(MMG5_pMesh mesh,MMG3D_pPROctree* q, int nv)
{
  MMG5_pTetra pt;
  double coor[3];
  int ie,iloc;

  MMG5_ADD_MEM(mesh,sizeof(MMG3D_PROctree),"PROctree structure",
                return 0);
  MMG5_SAFE_MALLOC(*q,1, MMG3D_PROctree, return 0);


  // set nv to the next power of 2
  nv--;
  nv |= nv >> 1;
  nv |= nv >> 2;
  nv |= nv >> 4;
  nv |= nv >> 8;
  nv |= nv >> 16;
  nv++;
  (*q)->nv = nv;

  // Number maximum of cells listed for the zone search
  (*q)->nc = MG_MAX(2048/nv,16);

  MMG5_ADD_MEM(mesh,sizeof(MMG3D_PROctree_s),"initial PROctree cell",
                return 0);

  MMG5_SAFE_MALLOC((*q)->q0,1, MMG3D_PROctree_s, return 0);
  MMG3D_initPROctree_s((*q)->q0);

  for( ie = 1; ie <= mesh->ne; ie++ ) {
    pt = &mesh->tetra[ie];
    for( iloc = 0; iloc < 4; iloc++ ) {
      memcpy(&coor,mesh->point[pt->v[iloc]].c,mesh->dim*sizeof(double));
      if( !PMMG_addPROctreeRec(mesh,(*q)->q0, coor, 4*ie+iloc, (*q)->nv) )
        return 0;
    }
  }

  return 1;
}

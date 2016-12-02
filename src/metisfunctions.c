/**
 * \file metisfunctions.c
 * \brief Partionning and repartitionning functions (with metis calls).
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */


#include "libparmmg.h"
#include "metis.h"

/**
 * \param parmesh pointer toward the mesh structure.
 * \param part pointer of an array containing the partitions (at the end).
 * \return 0 if the file is not found, -1 if we detect mismatch parameters,
 * 1 otherwise.
 *
 * Call metis to part the mesh.
 *
 */
int PMMG_metispartitioning(PMMG_pParMesh parmesh,int *part) {
  PMMG_pGrp  grp;
  MMG5_pMesh mesh;
  int        *xadj,*adjncy,*adja;
  int        j,k,iadr,jel,count,totCount,nelt,nproc,nbAdj,ncon,ier,objval;

  grp = parmesh->listgrp;
  mesh = grp[0].mesh;

  /* create tetra adjacency */
  if ( !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"  ## PMMG Hashing problem (1). Exit program.\n");
    return(0);
  }

  /*mesh -> graph*/
  _MMG5_SAFE_CALLOC(xadj,mesh->ne + 1,int);

  /*count neighboors*/
  xadj[0]  = 0;
  totCount = 0;
  for(k=1 ; k<=mesh->ne ; k++) {
    nbAdj = 0;
    iadr = 4*(k-1) + 1;
    adja = &mesh->adja[iadr];
    for(j=0 ; j<4 ; j++) {
      if(adja[j])   nbAdj++;
    }
    assert(nbAdj<=4);
    totCount += nbAdj;
    xadj[k] = totCount;
  }


  _MMG5_SAFE_CALLOC(adjncy,totCount + 1,int);

  count = 0;
  for(k=1 ; k<=mesh->ne ; k++) {
    iadr = 4*(k-1) + 1;
    adja = &mesh->adja[iadr];
    for(j=0 ; j<4 ; j++) {
      jel = adja[j]/4;

      if ( !jel ) continue;

      adjncy[count++] = jel-1;
    }

    if(count!=xadj[k]) printf("count %d %d %d\n",k,count,xadj[k]);
    assert(count==(xadj[k]));
  }

  /*call metis*/
  ncon       = 1;/*number of balancing constraint*/
  nelt       = mesh->ne;
  nproc      = parmesh->nprocs;

  ier =  METIS_PartGraphKway(&nelt,&ncon,xadj,adjncy,NULL/*vwgt*/,NULL/*vsize*/,
                             NULL/*adjwgt*/,&nproc,NULL/*tpwgts*/,
                             NULL/*ubvec*/,NULL/*options*/,&objval,part);
  /* ier =  METIS_PartGraphRecursive(&nelt,&ncon,xadj,adjncy,NULL/\*vwgt*\/,NULL/\*vsize*\/, */
  /*                         NULL/\*adjwgt*\/,&nproc,NULL/\*tpwgts*\/, */
  /*                         NULL/\*ubvec*\/,NULL/\*options*\/,&objval,part); */

  /**/

  _MMG5_SAFE_FREE(xadj);
  _MMG5_SAFE_FREE(adjncy);
  return(1);

}

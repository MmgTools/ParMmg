/**
 * \file metisfunctions.c
 * \brief Partionning and repartitionning functions (with metis calls).
 * \author Cécile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria/UBordeaux)
 * \version 5
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */


#include "parmmg.h"
#include "mpitypes.h"
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param xadj idx_t pointer.
 * \param adjncy idx_t pointer.
 * \return 1,
 *
 * Allocate and fill the xadj and adjncy arrays using for metis calls
 *
 */
int PMMG_mesh2metis(PMMG_pParMesh parmesh,idx_t **xadj,idx_t **adjncy) {
  PMMG_pGrp  grp;
  MMG5_pMesh mesh;
  int        *adja;
  int        j,k,iadr,jel,count,totCount,nbAdj;
  
  grp = parmesh->listgrp;
  mesh = grp[0].mesh;
  
  /* create tetra adjacency */
  if ( (!mesh->adja) && !MMG3D_hashTetra(mesh,1) ) {
    fprintf(stderr,"  ## PMMG Hashing problem (1). Exit program.\n");
    return(0);
  }
  
  /*mesh -> graph*/
  _MMG5_SAFE_CALLOC(*xadj,mesh->ne + 1,idx_t);
  
  /*count neighboors*/
  (*xadj)[0]  = 0;
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
    (*xadj)[k] = totCount;
  }
  
  
  _MMG5_SAFE_CALLOC(*adjncy,totCount + 1,idx_t);
  
  count = 0;
  for(k=1 ; k<=mesh->ne ; k++) {
    iadr = 4*(k-1) + 1;
    adja = &mesh->adja[iadr];
    for(j=0 ; j<4 ; j++) {
      jel = adja[j]/4;
      
      if ( !jel ) continue;
      
      (*adjncy)[count++] = jel-1;
    }
    
    if(count!=(*xadj)[k]) printf("count %d %d %d\n",k,count,(*xadj)[k]);
    assert(count==((*xadj)[k]));
  }
  
  return(1);
}
/**
 * \param parmesh pointer toward the mesh structure.
 * \param part pointer of an array containing the partitions (at the end).
 * \return 0 if the file is not found, -1 if we detect mismatch parameters,
 * 1 otherwise.
 *
 * Proc 0 call metis to part the mesh. Then broadcast the partition.
 *
 */
int PMMG_metispartitioning(PMMG_pParMesh parmesh,idx_t *part) {
  PMMG_pGrp  grp;
  MMG5_pMesh mesh;
  idx_t      *xadj,*adjncy;
  idx_t      nelt,ncon,nproc,objval, ier;

  if ( !parmesh->myrank ) {
    grp = parmesh->listgrp;
    mesh = grp[0].mesh;

    PMMG_mesh2metis(parmesh,&xadj,&adjncy);
 
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
	if ( ier != METIS_OK ) {
		fprintf(stderr, "Metis returned error value: " );
		switch ( ier ) {
			case METIS_ERROR_INPUT:
				fprintf(stderr, "METIS_ERROR_INPUT: input data error\n" );
			break;
			case METIS_ERROR_MEMORY:
				fprintf(stderr, "METIS_ERROR_MEMORY: could not allocate memory error\n" );
			break;
			case METIS_ERROR:
				fprintf(stderr, "METIS_ERROR: generic error\n" );
			break;
		}
		exit( EXIT_FAILURE );
	}

    /**/

    _MMG5_SAFE_FREE(xadj);
    _MMG5_SAFE_FREE(adjncy);
  }

  /* Send at each proc the graph */
  MPI_Bcast(&part[0],(parmesh->listgrp[0].mesh)->ne,MPI_INT,0,parmesh->comm);

  return(1);
}

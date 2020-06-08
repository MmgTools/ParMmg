#include "parmmgexterns.h"
#include "parmmg.h"

int (*PMMG_interp4bar)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTetra pt,int,double *phi)=NULL;
int (*PMMG_interp3bar)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTria ptr,int,double *phi)=NULL;

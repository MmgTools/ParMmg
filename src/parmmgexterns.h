#include "parmmg.h"

extern int (*PMMG_interp4bar)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTetra pt,int,PMMG_barycoord*);
extern int (*PMMG_interp3bar)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTria ptr,int,PMMG_barycoord*);
extern int (*PMMG_interp2bar)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTria ptr,int ip,int l,PMMG_barycoord *barycoord);

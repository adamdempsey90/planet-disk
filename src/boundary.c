#include "planetdisk.h"

void set_bc(Field *fld) {
/* Set the boundary conditions for the right boundary */
	int i,j, indx;
	double sig0 = (fld->Params->sig0);
	

/* Zero B.C */	
	for(i=0;i<NG;i++) {
		for(j=0;j<NC;j++) {
			indx = j+NC*(NG-1 + NG-i);
			fld->u[CINDX] = -conj(fld->u[indx]);
			fld->v[CINDX] = -conj(fld->v[indx]);
			fld->sig[CINDX] = conj(fld->sig[indx]);
		}
	}	
	
	for(i=Nx+NG;i<Nx+2*NG;i++) {
		for(j=0;j<NC;j++) {
			fld->u[CINDX] = 0;
			fld->v[CINDX] = 0;
			
			if(j==0) fld->sig[CINDX] = sig0;
			else	fld->sig[CINDX] = 0;
			

		}
	}

	return;
}

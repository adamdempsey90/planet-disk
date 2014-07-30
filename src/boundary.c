#include "planetdisk.h"

void set_bc(Field *fld) {
/* Set the boundary conditions for the right boundary */
	int i;
	
	for(i=0;i<NC;i++) {
		if(i==0) {
			fld->sigbc[i] = fld->params->sig0;
		}
		else {
			fld->sigbc[i] = 0;
		}
		fld->ubc[i] = 0;
		fld->vbc[i] = 0;
		fld->Tens->Pixxbc[i] = -(fld->c)*(fld->c)*(fld->sigbc[i]);
		fld->Tens->Pixybc[i] = -(fld->nu)*(fld->q)*(fld->omega)*(fld->sigbc[i]);
		fld->Tens->Piyybc[i] = -(fld->c)*(fld->c)*(fld->sigbc[i]);
		
	}

	return;
}



#include "planetdisk.h"

void set_bc(Field *fld) {
/* Set the boundary conditions for the right boundary */
	int i;
	double c = (fld->Params->c)*(fld->Params->c);
	double nuqom = (fld->Params->nu)*(fld->Params->omega)*(fld->Params->q);
	double sig0 = (fld->Params->sig0);
	for(i=0;i<NC;i++) {
		if(i==0) {
			fld->sigbc[i] = sig0;
		}
		else {
			fld->sigbc[i] = 0;
		}
		fld->ubc[i] = 0;
		fld->vbc[i] = 0;
		fld->Tens->Pixxbc[i] = -c*(fld->sigbc[i]);
		fld->Tens->Pixybc[i] = -nuqom*(fld->sigbc[i]);
		fld->Tens->Piyybc[i] = -c*(fld->sigbc[i]);
		
	}

	return;
}



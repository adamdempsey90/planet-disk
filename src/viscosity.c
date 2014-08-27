#include "planetdisk.h"

void visc_tens(Field *fld) {
/* Calculates the stress tensor including viscosity and pressure */

	int i;
	double twoth = (2.0/3.0);
	double qom = (fld->Params->q)*(fld->Params->omega);
	double nu = (fld->Params->nu);
	double complex divv;
	
/* Forms Navier-Stokes stress tensor T=grad(v) + grad(v)^T - 2/3 div(v) I */


	for(i=istart;i<iend;i++) {
	
		divv = fld->dxu[i] + fld->dyv[i];
		divv *= twoth;
		fld->Tens->Txx[i] = nu*(2*(fld->dxu[i])-divv);
		fld->Tens->Txy[i] = nu*(fld->dxv[i] + fld->dyu[i]) - nu*qom;
		fld->Tens->Tyy[i] = nu*(2*(fld->dyv[i]) - divv);
	
	}

/* Set Tensor B.C's */

	set_tens_bc(fld);

	return;
}

void add_visc(Field *fld) {

	visc_tens(fld);
	
/* Calculate div(Pi) */	
	calc_deriv(fld->Tens->Txx,fld->dtu,NULL,fld->Params->dx,fld->kk);
	calc_deriv(fld->Tens->Txy,fld->dtv,fld->dtu,fld->Params->dx,fld->kk);
	calc_deriv(fld->Tens->Tyy,NULL,fld->dtv,fld->Params->dx,fld->kk);
	

/* Convolve T \dot log(Sigma) */

 	convolve(&fld->dxsig[istart],&fld->Tens->Txx[istart],fld->dtu,1);
 	convolve(&fld->dysig[istart],&fld->Tens->Txy[istart],fld->dtu,1);
  	convolve(&fld->dxsig[istart],&fld->Tens->Txy[istart],fld->dtv,1);
 	convolve(&fld->dysig[istart],&fld->Tens->Tyy[istart],fld->dtv,1);
	
	return;
}

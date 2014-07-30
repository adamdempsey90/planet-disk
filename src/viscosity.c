#include "meanwave.h"

void visc_tens(Field *fld) {
	int i;
	double twoth = (2.0/3.0);
	double qom = (fld->q)*(fld->omega);
	double c2 = (fld->c)*(fld->c);
	for(i=0;i<NTOTC;i++) {
		divv = fld->dxu[i] + I*(fld->k[i])*(fld->v[i]);
		divv *= twoth;
		fld->Tens->Txx[i] = 2*(fld->dxu[i])-divv;
		fld->Tens->Txy[i] = fld->dxv[i] + I*(fld->k[i])*(fld->u[i]);
		fld->Tens->Tyy[i] = 2*I*(fld->k[i])*(fld->v[i]) - divv;
	}
	convolve(fld->Tens->Txx,fld->sig,fld->Tens->Pixx,fld->nu);
	convolve(fld->Tens->Txy,fld->sig,fld->Tens->Pixy,fld->nu);
	convolve(fld->Tens->Tyy,fld->sig,fld->Tens->Piyy,fld->nu);
	for(i=0;i<NTOTC;i++) {
		fld->Tens->Pixx[i] -= c*(fld->sig[i]);
		fld->Tens->Pixy[i] -= (fld->nu)*qom*(fld->sig[i]);
		fld->Tens->Piyy[i] -= c*(fld->sig[i]);
	}

	return;
}

void add_visc(Field *fld) {
	int i;

/* Calculate div(Pi) */	
calc_deriv(fld->u,fld->dxu,fld->params,"odd",fld->ubc,fld->dx);
	calc_deriv(fld->Tens->Pixx,fld->Tens->divPix,NULL,"even",fld->Pixxbc,fld->dx);
	calc_deriv(fld->Tens->Pixy,fld->Tens->divPiy,NULL,"even",fld->Pixybc,fld->dx);
	
	for(i=0;i<NTOTC;i++) {
		fld->Tens->divPix[i] += I*(fld->k[i])*(fld->Tens->Pixy[i]);
		fld->Tens->divPiy[i] += I*(fld->k[i])*(fld->Tens->Piyy[i]);
	}

/* Convolve with 1/Sigma */

	convolve_inv(fld->sig,fld->Tens->divPix,fld->dtu);
	convolve_inv(fld->sig,fld->Tens->divPiy,fld->dtv);
	
	return;
}
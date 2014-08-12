#include "planetdisk.h"

void visc_tens(Field *fld) {
/* Calculates the stress tensor including viscosity and pressure */

	int i,j,indx;
	double twoth = (2.0/3.0);
	double qom = (fld->Params->q)*(fld->Params->omega);
	double c = (fld->Params->c)*(fld->Params->c);
	double nu = (fld->Params->nu);
	double complex divv;
	
/* Forms Navier-Stokes stress tensor T=grad(v) + grad(v)^T - 2/3 div(v) I */

	for(i=istart;i<iend;i++) {
		divv = fld->dxu[i] + fld->dyv[i];
		divv *= twoth;
		fld->Tens->Txx[i] = 2*(fld->dxu[i])-divv;
		fld->Tens->Txy[i] = fld->dxv[i] + fld->dyu[i];
		fld->Tens->Tyy[i] = 2*(fld->dyv[i]) - divv;
		fld->Tens->Pixx[i] = -c*(fld->sig[i]);
		fld->Tens->Pixy[i] = -nu*qom*(fld->sig[i]);
		fld->Tens->Piyy[i] = -c*(fld->sig[i]);	
	}
	convolve(&fld->Tens->Txx[istart],&fld->sig[istart],&fld->Tens->Pixx[istart],nu);
	convolve(&fld->Tens->Txy[istart],&fld->sig[istart],&fld->Tens->Pixy[istart],nu);
	convolve(&fld->Tens->Tyy[istart],&fld->sig[istart],&fld->Tens->Piyy[istart],nu);

/* Set Tensor B.C's */

	for(i=0;i<NG;i++) {
		for(j=0;j<NC;j++) {
			indx = j+NC*(NG-1 + NG-i);
			fld->Tens->Pixx[CINDX] = conj(fld->Tens->Pixx[indx]);
			fld->Tens->Pixy[CINDX] = conj(fld->Tens->Pixy[indx]);
			fld->Tens->Piyy[CINDX] = conj(fld->Tens->Piyy[indx]);
		}
	}
	for(i=Nx+NG;i<Nx+2*NG;i++) {
		for(j=0;j<NC;j++) {
			fld->Tens->Pixx[CINDX] = -c*(fld->sig[CINDX]);
			fld->Tens->Pixy[CINDX] = -nu*qom*(fld->sig[CINDX]);
			fld->Tens->Piyy[CINDX] = -c*(fld->sig[CINDX]);
		}
	}
			

	return;
}

void add_visc(Field *fld) {

	visc_tens(fld);
	
/* Calculate div(Pi) */	
	calc_deriv(fld->Tens->Pixx,fld->Tens->divPix,NULL,fld->Params->dx,fld->k);
	calc_deriv(fld->Tens->Pixy,fld->Tens->divPiy,fld->Tens->divPix,fld->Params->dx,fld->k);
	calc_deriv(fld->Tens->Piyy,NULL,fld->Tens->divPiy,fld->Params->dx,fld->k);
	

/* Convolve with 1/Sigma */

 	convolve_inv(&fld->sig[istart],&fld->Tens->divPix[istart],fld->dtu,1);
 	convolve_inv(&fld->sig[istart],&fld->Tens->divPiy[istart],fld->dtv,1);
	
	return;
}

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
	convolve(fld->Tens->Txx,fld->sig,fld->Tens->Pixx,nu);
	convolve(fld->Tens->Txy,fld->sig,fld->Tens->Pixy,nu);
	convolve(fld->Tens->Tyy,fld->sig,fld->Tens->Piyy,nu);

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
	int i,j;

	visc_tens(fld);
	
/* Calculate div(Pi) */	
	calc_deriv(fld->Tens->Pixx,fld->Tens->divPix,NULL,
						fld->Params->dx,fld->k);
	calc_deriv(fld->Tens->Pixy,fld->Tens->divPiy,NULL,
						fld->Params->dx,fld->k);
	
	for(i=NG;i<Nx+NG;i++) {
		for(j=0;j<NC;j++) {
			fld->Tens->divPix[i] += I*(fld->k[j])*(fld->Tens->Pixy[i]);
			fld->Tens->divPiy[i] += I*(fld->k[j])*(fld->Tens->Piyy[i]);
		}
	}

/* Convolve with 1/Sigma */

	convolve_inv(fld->sig,fld->Tens->divPix,fld->dtu,1);
	convolve_inv(fld->sig,fld->Tens->divPiy,fld->dtv,1);
	
	return;
}

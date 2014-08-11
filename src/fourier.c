#include "planetdisk.h"
#include <fftw3.h>


double complex *wc1, *wc2, *wc3;
double *wr1, *wr2, *wr3;
double *mask;

fftw_plan r2c1, c2r1,r2c2, c2r2,r2c3, c2r3;

void convolve(double complex *q1, double complex *q2, double complex *res, double complex mult) {
/*		Convolution function using previously declared de-aliasing mask
	Inputs: q1 & q2 are the two complex arrays being convolved 
			res is an initialized array to which the convolution will be added 
			mult is a constant multiplication factor 
*/

	int i;

/* De-alias with 2/3 truncation rule */	
	for(i=0;i<NC*Nx;i++) {
		wc1[i] = q1[i]*mask[i]*mult; 
		wc2[i] = q2[i]*mask[i];
	}

/* FFT these to real space */	
	fftw_execute(c2r1);
	fftw_execute(c2r2);

/* Form product in real space */
	for(i=0;i<Nx*NR;i++) wr3[i] = wr1[i]*wr2[i]/(Ny*Ny);

/* FFT back to complex space */
	fftw_execute(r2c3);
	
/* add output */
	for(i=0;i<NC*Nx;i++) {
		res[i] += wc3[i];
	}	

	return;
}
void convolve_inv(double complex *q1, double complex *q2, double complex *res, double complex mult) {
/*		Convolution function using previously declared de-aliasing mask
		This function uses the inverse of q1, so make sure q1 != 0 anywhere.
	Inputs: q1 & q2 are the two complex arrays being convolved 
			res is an initialized array to which the convolution will be added 
			mult is a constant multiplication factor 
*/

	int i;

/* De-alias with 2/3 truncation rule */	
	for(i=0;i<NC*Nx;i++) {
		wc1[i] = mask[i]*q1[i]; 
		wc2[i] = q2[i]*mask[i]*mult;
	}

/* FFT these to real space */	
	fftw_execute(c2r1);
	fftw_execute(c2r2);

/* Form product in real space */
	for(i=0;i<Nx*NR;i++) wr3[i] = wr2[i]/(wr1[i]*Ny*Ny);

/* FFT back to complex space */
	fftw_execute(r2c3);
	
/* add output */
	for(i=0;i<Nx*NC;i++) {
		res[i] += wc3[i];
	}	

	return;
}


void init_fft(void) {
	int i,j,Nmax;
	wc1 = (double complex *)malloc(sizeof(double complex)*Nx*NC);
	wr1 = (double *)wc1;
	wc2 = (double complex *)malloc(sizeof(double complex)*Nx*NC);
	wr2 = (double *)wc2;
	wc3 = (double complex *)malloc(sizeof(double complex)*Nx*NC);
	wr3 = (double *)wc3;
	mask = (double *)malloc(sizeof(double)*Nx*NC);
	
	Nmax = (2./3)*(NC-1);
	
	for(i=0;i<Nx;i++) {
		for(j=0;j<NC;j++) {
			if (j < Nmax) mask[CINDX] = 1;
			else mask[CINDX] = 0;
		}
	}
	
	c2r1=fftw_plan_many_dft_c2r(1,&Nx,Ny,wc1,&Nx,1,NC,wr1,&Nx,1,NR,FFTW_ESTIMATE);
	if (c2r1 == NULL) printf("Problem with c2r1\n");
	r2c1=fftw_plan_many_dft_r2c(1,&Nx,Ny,wr1,&Nx,1,NR,wc1,&Nx,1,NC,FFTW_ESTIMATE);
	if (r2c1 == NULL) printf("Problem with r2c1\n");
	
	c2r2=fftw_plan_many_dft_c2r(1,&Nx,Ny,wc2,&Nx,1,NC,wr2,&Nx,1,NR,FFTW_ESTIMATE);
	if (c2r2 == NULL) printf("Problem with c2r2\n");
	r2c2=fftw_plan_many_dft_r2c(1,&Nx,Ny,wr2,&Nx,1,NR,wc2,&Nx,1,NC,FFTW_ESTIMATE);
	if (r2c2 == NULL) printf("Problem with r2c2\n");
	
	c2r3=fftw_plan_many_dft_c2r(1,&Nx,Ny,wc3,&Nx,1,NC,wr3,&Nx,1,NR,FFTW_ESTIMATE);
	if (c2r3 == NULL) printf("Problem with c2r3\n");
	r2c3=fftw_plan_many_dft_r2c(1,&Nx,Ny,wr3,&Nx,1,NR,wc3,&Nx,1,NC,FFTW_ESTIMATE);
	if (r2c3 == NULL) printf("Problem with r2c3\n");
	
	
	return;

}

void fft_free(void) {

	fftw_destroy_plan(r2c1);
	fftw_destroy_plan(c2r1);
	fftw_destroy_plan(r2c2);
	fftw_destroy_plan(c2r2);
	fftw_destroy_plan(r2c3);
	fftw_destroy_plan(c2r3);	
	
	free(wc1);free(wc2);free(wc3);
	free(mask);
	return;


}

void fft_phi(double *rphi, double complex *cphi) {
	memcpy(wr1,&rphi[NG*NR],sizeof(double)*Nx*NR);
	fftw_execute(r2c1);
	memcpy(&cphi[NG*NC],wc1,sizeof(double complex)*Nx*NC);
	return;
}
void fft_dxphi(double *rdxphi, double complex *cdxphi) {
	memcpy(wr1,&rdxphi[NG*NR],sizeof(double)*Nx*NR);
	fftw_execute(r2c1);
	memcpy(&cdxphi[NG*NC],wc1,sizeof(double complex)*Nx*NC);
	return;
}



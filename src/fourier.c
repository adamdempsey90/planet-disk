#include "planetdisk.h"
#include <fftw3.h>


int fNC, fNR;
double complex *wc1, *wc2, *wc3;
double *wr1, *wr2, *wr3;
double *mask;

fftw_plan r2c1, c2r1,r2c2, c2r2,r2c3, c2r3;

/* This file contains all FFFTW function calls.
	We are using real to complex transforms with no default normalization.
	We normalize by Ny when doing a r2c transform
*/


void convolve(double complex *q1, double complex *q2, double complex *res, double complex mult) {
/*		Convolution function using previously declared de-aliasing mask
	Inputs: q1 & q2 are the two complex arrays being convolved 
			res is an initialized array to which the convolution will be added 
			mult is a constant multiplication factor 
*/

	int i,j;

#ifndef LINEAR
/* De-alias with 2/3 truncation rule */	
	for(i=0;i<NC*Nx;i++) {
		wc1[i] = q1[i]*mask[i]*mult; 
		wc2[i] = q2[i]*mask[i];
	}

/* FFT these to real space */	
	fftw_execute(c2r1);
	fftw_execute(c2r2);

/* Form product in real space */

	for(i=0;i<Nx*NR;i++) wr3[i] = wr1[i]*wr2[i];

/* FFT back to complex space */
	fftw_execute(r2c3);
	
/* add output */

	for(i=0;i<Nx;i++) {
		for(j=0;j<Nmax;j++) {
			res[CINDX] += wc3[CINDX]/Ny;
		}
	}	
#else
/* No interaction b/w modes except with zero mode */
	for(i=0;i<Nx;i++) {
		res[i*NC] += q1[i*NC]*q2[i*NC]*mult;
		for(j=1;j<NC;j++) {
			res[i*NC] += 2*creal(conj(q1[CINDX])*q2[CINDX]*mult);
			res[CINDX] += (q1[i*NC]*q2[CINDX] + q1[CINDX]*q2[i*NC])*mult;
		}
	}
#endif
	return;
}
void convolve_inv(double complex *q1, double complex *q2, double complex *res, double complex mult) {
/*		Convolution function using previously declared de-aliasing mask
		This function uses the inverse of q1, so make sure q1 != 0 everywhere.
	Inputs: q1 & q2 are the two complex arrays being convolved 
			res is an initialized array to which the convolution will be added 
			mult is a constant multiplication factor 
*/

	int i,j;

#ifndef LINEAR
/* De-alias with 2/3 truncation rule */	
	
	for(i=0;i<Nx*NC;i++) {
		wc1[i] = mask[i]*q1[i]; 
		wc2[i] = q2[i]*mask[i]*mult;
//		printf("%lg %lg\n",cabs(wc1[i]), cabs(wc2[i]));
	}

/* FFT these to real space */	
	fftw_execute(c2r1);
	fftw_execute(c2r2);

/* Form product in real space */

	for(i=0;i<Nx;i++) {
		for(j=0;j<NR;j++) {
			if(wr1[RINDX]!=0) wr3[RINDX] = wr2[RINDX]/wr1[RINDX];
			else wr3[RINDX]= 0;
		}
	}
/* FFT back to complex space */
	fftw_execute(r2c3);
	
/* add output */
	
	for(i=0;i<Nx;i++) {
		for(j=0;j<Nmax;j++) {
			res[CINDX] += wc3[CINDX]/Ny;
		}
	}	
#else
/* Expand out to second order in q1 w.r.t zero mode q1 */
	for(i=0;i<Nx;i++) {
		res[i*NC] += mult*q2[i*NC]/q1[i*NC];
		for(j=1;j<NC;j++) {
			res[i*NC] += 2*creal(pow(q1[i*NC],-3.0)*(-q1[i*NC]*conj(q1[CINDX])*q2[CINDX]
								+ conj(q1[CINDX])*q1[CINDX]*q2[i*NC])*mult);
			res[CINDX] += pow(q1[i*NC],-2.0)*(q1[i*NC]*q2[CINDX]-q1[CINDX]*q2[i*NC])*mult;
		}
	}

#endif
	return;
}


void init_fft(void) {
	int i,j;
	
	int n[] = {Ny};
	int howmany = Nx;
	int idist=NC, odist=NR;
	int istride=1, ostride=1;
	int *inembed=NULL, *onembed=NULL;

#ifdef LINEAR
	Nmax = NC-1;
#else
	Nmax = floor(2./3*NC);
#endif	
	wc1 = (double complex *)malloc(sizeof(double complex)*Nx*NC);
	wr1 = (double *)wc1;
	wc2 = (double complex *)malloc(sizeof(double complex)*Nx*NC);
	wr2 = (double *)wc2;
	wc3 = (double complex *)malloc(sizeof(double complex)*Nx*NC);
	wr3 = (double *)wc3;

	mask = (double *)malloc(sizeof(double)*Nx*NC);
	
	
	for(i=0;i<Nx;i++) {
		for(j=0;j<NC;j++) {
			if (j < Nmax) {
				mask[CINDX] = 1;
#ifdef FILTER
				mask[CINDX] *= .5*(1+cos(j*M_PI/Nmax));
				
//				if (j!=0) mask[CINDX] *= sin(j*M_PI/Nmax)/(j*M_PI/Nmax);
#endif	
			}
			else mask[CINDX] = 0;
		}
	}
	
// //	c2r1=fftw_plan_many_dft_c2r(1,&Nx,Ny,wc1,&Nx,1,NC,wr1,&Nx,1,NR,FFTW_ESTIMATE);
// 	if (c2r1 == NULL) printf("Problem with c2r1\n");
// //	r2c1=fftw_plan_many_dft_r2c(1,&Nx,Ny,wr1,&Nx,1,NR,wc1,&Nx,1,NC,FFTW_ESTIMATE);
// 	if (r2c1 == NULL) printf("Problem with r2c1\n");
// 	
// //	c2r2=fftw_plan_many_dft_c2r(1,&Nx,Ny,wc2,&Nx,1,NC,wr2,&Nx,1,NR,FFTW_ESTIMATE);
// 	if (c2r2 == NULL) printf("Problem with c2r2\n");
// //	r2c2=fftw_plan_many_dft_r2c(1,&Nx,Ny,wr2,&Nx,1,NR,wc2,&Nx,1,NC,FFTW_ESTIMATE);
// 	if (r2c2 == NULL) printf("Problem with r2c2\n");
// 	
// //	c2r3=fftw_plan_many_dft_c2r(1,&Nx,Ny,wc3,&Nx,1,NC,wr3,&Nx,1,NR,FFTW_ESTIMATE);
// 	if (c2r3 == NULL) printf("Problem with c2r3\n");
// //	r2c3=fftw_plan_many_dft_r2c(1,&Nx,Ny,wr3,&Nx,1,NR,wc3,&Nx,1,NC,FFTW_ESTIMATE);
// 	if (r2c3 == NULL) printf("Problem with r2c3\n");
// 	
	c2r2=fftw_plan_many_dft_c2r(1,n,howmany,wc2, inembed, istride, idist, wr2, onembed, ostride, odist,FFTW_ESTIMATE);
	if (c2r2 == NULL) printf("Problem with c2r1 on processor %d\n",rank);
	r2c2=fftw_plan_many_dft_r2c(1,n,howmany,wr2, onembed, ostride, odist, wc2, inembed, istride, idist,FFTW_ESTIMATE);
	if (r2c2 == NULL) printf("Problem with c2r1 on processor %d\n",rank);

	c2r1=fftw_plan_many_dft_c2r(1,n,howmany,wc1, inembed, istride, idist, wr1, onembed, ostride, odist,FFTW_ESTIMATE);
	if (c2r1 == NULL) printf("Problem with c2r1 on processor %d\n",rank);

	r2c1=fftw_plan_many_dft_r2c(1,n,howmany,wr1, onembed, ostride, odist, wc1, inembed, istride, idist,FFTW_ESTIMATE);
	if (r2c1 == NULL) printf("Problem with c2r1 on processor %d\n",rank);


	c2r3=fftw_plan_many_dft_c2r(1,n,howmany,wc3, inembed, istride, idist, wr3, onembed, ostride, odist,FFTW_ESTIMATE);
	if (c2r3 == NULL) printf("Problem with c2r1 on processor %d\n",rank);
	
	r2c3=fftw_plan_many_dft_r2c(1,n,howmany,wr3, onembed, ostride, odist, wc3, inembed, istride, idist,FFTW_ESTIMATE);
	if (r2c3 == NULL) printf("Problem with c2r1 on processor %d\n",rank);

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
	int i;
	memcpy(wr1,rphi,sizeof(double)*Nx*NR);
	fftw_execute(r2c1);
	for(i=istart;i<iend;i++) cphi[i] = wc1[i-istart]/Ny;
		
//	memcpy(&cphi[istart],wc1,sizeof(double complex)*Nx*NC);
	return;
}
void transform(Field *fld) {
	
	memcpy(wc1,&(fld->u[istart]),sizeof(double complex)*Nx*NC);
	memcpy(wc2,&(fld->v[istart]),sizeof(double complex)*Nx*NC);
	memcpy(wc3,&(fld->sig[istart]),sizeof(double complex)*Nx*NC);
	fftw_execute(c2r1);
	fftw_execute(c2r2);
	fftw_execute(c2r3);
	memcpy(fld->vx,wr1,sizeof(double)*Nx*NR);
	memcpy(fld->vy,wr2,sizeof(double)*Nx*NR);
	memcpy(fld->dens,wr3,sizeof(double)*Nx*NR);
	return;

}



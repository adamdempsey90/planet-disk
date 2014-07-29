#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define NR 2*(Ny/2+1)
#define NC (Ny/2+1)
#define RINDX  j+NR*i
#define CINDX	j+NC*i
//#define OPENMP
fftw_plan r2c, c2r;
int Nx,Ny;
double complex *cphi;
double *rphi, *tphi;
double *wri;
double complex *wco;
void fft_init(void);
void fft_free(void);
void fft_r2c(int start);
void fft_c2r(int start);
int main(void) {
	int i,j;
	FILE *rfp, *cfp, *xfp, *yfp;
	Nx = 256; Ny=256;	
	cphi = (double complex *)malloc(sizeof(double complex)*Nx*(Ny/2+1));
	rphi = (double *)cphi;
	double *x,*y;
	x = (double *)malloc(sizeof(double)*Nx);
	y = (double *)malloc(sizeof(double)*Ny);


	for(i=0;i<Nx;i++) x[i]=-15 + i*(30.0/Nx);
	for(i=0;i<Ny;i++) y[i]=-15 + i*30.0/Ny;
	
	rfp = fopen("phi.dat","w");
	xfp = fopen("x.dat","w");
	yfp = fopen("y.dat","w");

	for(i=0;i<Nx;i++) {
		for(j=0;j< Ny;j++) {
			rphi[RINDX] = 	-1.0/sqrt(x[i]*x[i]+y[j]*y[j]+.6*.6);		
			fprintf(rfp,"%lg ", rphi[RINDX]);
			fprintf(xfp,"%lg ",x[i]);
			fprintf(yfp,"%lg ",y[j]);
		}
		fprintf(rfp,"\n"); 
		fprintf(xfp,"\n");
		fprintf(yfp,"\n");
	}
	fclose(rfp); fclose(xfp); fclose(yfp);
	
	rfp = fopen("rphi.dat","w");
	cfp = fopen("cphi.dat","w");
	printf("Inializing...\n");

	c2r=fftw_plan_many_dft_c2r(1,&Nx,Ny,cphi,&Nx,1,NC,rphi,&Nx,1,NR,FFTW_ESTIMATE);
	if (c2r == NULL) printf("Problem with c2r\n");
	r2c=fftw_plan_many_dft_r2c(1,&Nx,Ny,rphi,&Nx,1,NR,cphi,&Nx,1,NC,FFTW_ESTIMATE);
	if (r2c == NULL) printf("Problem with r2c\n");

//	fft_init();
	printf("Executing...\n");

	fftw_execute(r2c);
	for(i=0;i<Nx*NC;i++) cphi[i] /= Ny;
	printf("Outputting...\n");
	for(i=0;i<Nx;i++) {
		for(j=0;j<NC;j++) {
			fprintf(rfp,"%lg ",creal(cphi[CINDX]));
			fprintf(cfp,"%lg ",cimag(cphi[CINDX]));
		}
		fprintf(rfp,"\n");
		fprintf(cfp,"\n");
	}
	fclose(rfp); fclose(cfp);
	fftw_execute(c2r);
 	rfp = fopen("phi2.dat","w");
 	for(i=0;i<Nx;i++) {
 		for(j=0;j< Ny;j++) {
 			fprintf(rfp,"%lg ", rphi[RINDX]);
 		}
 		fprintf(rfp,"\n");
 	}
 	fclose(rfp);
	
	free(x); free(y);
	free(cphi); 
	free(wco);
	fft_free();
	return;
}

void fft_free(void) {
	fftw_destroy_plan(r2c);
	fftw_destroy_plan(c2r);
	return;
}
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>


//#define OPENMP
fftw_plan r2c, c2r;
int Nx,Ny;
double complex *cphi;
double *rphi;

void fft_init(void);
void fft_free(void);
int main(void) {
	int i,j;
	FILE *fp, *ofp;
	Ny = 16; 
	cphi = (double complex *)malloc(sizeof(double complex)*Nx*(Ny/2+1));
//	rphi = (double *)malloc(sizeof(double)*2*(Ny/2+1));
	rphi = (double *)cphi;
	double *y;
	y = (double *)malloc(sizeof(double)*Ny);
	for (i=0;i<Ny;i++) y[i] = -M_PI + i*2*M_PI/Ny;
	
	fp = fopen("rphi.dat","w");
	for(j=0;j<Ny;j++) {
			rphi[j] = sin(y[j]);
			fprintf(fp,"%lg \t %lg \n ", y[j],rphi[j]);
	}

	
	fclose(fp);
	
	
	printf("Inializing...\n");
	fft_init();
	
	printf("Executing...\n");
	fftw_execute(r2c);
	printf("Outputting...\n");
	fp = fopen("cphi.dat","w");

	for(j=0;j<Ny/2+1;j++) {
			cphi[j] /= Ny;
			fprintf(fp,"%lg \t %lg \n", creal(cphi[j]),cimag(cphi[j]));
	}

	
	
	fclose(fp);
	
	fftw_execute(c2r);
	fp = fopen("rphi2.dat","w");
	
	for(j=0;j<Ny;j++) fprintf(fp,"%lg %lg\n",y[j],rphi[j]);
	fclose(fp);
	free(y);
	free(cphi); 
	//free(rphi);
	fft_free();
	return;
}
void fft_init(void) {


	c2r = fftw_plan_dft_c2r_1d(Ny, cphi,rphi,FFTW_ESTIMATE);
	if (c2r == NULL) printf("Problem with c2r\n");
	r2c = fftw_plan_dft_r2c_1d(Ny, rphi,cphi,FFTW_ESTIMATE);
	if (r2c == NULL) printf("Problem with r2c\n");


	return;
}

void fft_free(void) {
	fftw_destroy_plan(r2c);
	fftw_destroy_plan(c2r);
	return;
}
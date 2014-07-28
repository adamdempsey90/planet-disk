#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define RINDX  i+Nx*j
#define CINDX	j+(Ny/2+1)*i
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
	wco = (double complex *)malloc(sizeof(double complex)*(Ny/2+1));
	wri = (double *)wco;
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
	fft_init();
	printf("Executing...\n");
	for(i=0;i<Nx;i++) {
		fft_r2c(i*2*(Ny/2+1));
	}
	printf("Outputting...\n");
	for(i=0;i<Nx;i++) {
		for(j=0;j<Ny/2+1;j++) {
			fprintf(rfp,"%lg ",creal(cphi[CINDX]));
			fprintf(cfp,"%lg ",cimag(cphi[CINDX]));
		}
		fprintf(rfp,"\n");
		fprintf(cfp,"\n");
	}
	fclose(rfp); fclose(cfp);
	for(i=0;i<Nx;i++)	fft_c2r(i*(Ny/2+1));
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
void fft_r2c(int j) {
	int i;
	for(i=0;i<2*(Ny/2+1);i++) {
		wri[i] = rphi[2*(Ny/2+1)*j+i];
	}
	fftw_execute(r2c);
	for(i=0;i<Ny/2+1;i++) {
		cphi[(Ny/2+1)*j+i] = wco[i]/Ny;
	}
	
	return;

}
void fft_c2r(int j) {
	int i;
	for(i=0;i<Ny/2+1;i++) {
		wco[i] = cphi[j*(Ny/2+1)+i];
	}
	fftw_execute(c2r);
	for(i=0;i<2*(Ny/2+1);i++) {
		rphi[2*(Ny/2+1)*j+i] = wri[i];
	}	
	return;

}
void fft_init(void) {

// 	int rank = 1; /* not 2: we are computing 1d transforms */
//     int n[] = {Nx}; /* 1d transforms of length 10 */
// 	int howmany = Ny;
//     int idist = 1,odist = 1; 
//     int rstride = Ny , cstride=Ny; /* distance between two elements in
//                                       the same column */
//     int *inembed = n, *onembed = n;

	
// 	
// 	c2r = fftw_plan_many_dft_c2r(rank, n, howmany,
//                                   cphi, inembed,
//                                   cstride, idist,
//                                   rphi, onembed,
//                                   rstride, odist,
//                                   FFTW_ESTIMATE);
// 	if (c2r == NULL) printf("Problem with c2r\n");
// 	r2c = fftw_plan_many_dft_r2c(rank, n, howmany,
//                                   rphi, inembed,
//                                   rstride, idist,
//                                   cphi, onembed,
//                                   cstride, odist,
//                                   FFTW_ESTIMATE);
//     if (r2c == NULL) printf("Problem with r2c\n");

	c2r = fftw_plan_dft_c2r_1d(Ny, wco,wri,FFTW_ESTIMATE);
	if (c2r == NULL) printf("Problem with c2r\n");
	r2c = fftw_plan_dft_r2c_1d(Ny, wri,wco,FFTW_ESTIMATE);
	if (r2c == NULL) printf("Problem with r2c\n");
	
	
	return;
}


void fft_free(void) {
	fftw_destroy_plan(r2c);
	fftw_destroy_plan(c2r);
	return;
}
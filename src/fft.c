#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define RINDX  j+2*(Ny/2+1)*i
#define CINDX	j+(Ny/2+1)*i
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
	Nx = 32; Ny=16;	
	rphi = (double *)malloc(sizeof(double)*Nx*2*(Ny/2+1));
	cphi = (double complex *)rphi;
	double *x,*y;
	x = (double *)malloc(sizeof(double)*Nx);
	y = (double *)malloc(sizeof(double)*2*(Ny/2+1));
	
	for(i=0;i<Nx;i++) x[i]=-15 + i*(30.0/Nx);
	for(i=0;i<2*(Ny/2+1);i++) y[i]=-15 + i*30.0/(2*(Ny/2+1));
	
	fp = fopen("rphi.dat","w");
	for(i=0;i<Nx;i++) {
		for(j=0;j<2*(Ny/2+1);j++) {
			rphi[RINDX] = 	-1.0/sqrt(x[i]*x[i]+y[j]*y[j]+.6*.6);	
			fprintf(fp,"%lg ", rphi[RINDX]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	fp = fopen("cphi.dat","w");
	
	printf("Inializing...\n");
	fft_init();
	printf("Executing...\n");
	fftw_execute(r2c);
	printf("Outputting...\n");
	for(i=0;i<Nx;i++) {
		for(j=0;j<Ny/2+1;j++) {
			fprintf(fp,"%lg+I%lg ", creal(cphi[CINDX]),cimag(cphi[CINDX]));
		}
		fprintf(fp,"\n");
	}
	
	fclose(fp);
	free(x); free(y);
	free(cphi); free(rphi);
	fft_free();
	return;
}
void fft_init(void) {

	int rank = 1; /* not 2: we are computing 1d transforms */
    int n[] = {Nx}; /* 1d transforms of length 10 */
	int howmany = Ny;
    int idist = 1,odist = 1; 
    int rstride = 2*(Ny/2+1) , cstride = (Ny/2+1); /* distance between two elements in
                                      the same column */
    int *inembed = n, *onembed = n;
#ifdef OPENMP	
	fftw_plan_with_nthreads( 1 );
#endif
	
	
	c2r = fftw_plan_many_dft_c2r(rank, n, howmany,
                                  cphi, inembed,
                                  cstride, idist,
                                  rphi, onembed,
                                  rstride, odist,
                                  FFTW_BACKWARD);

	r2c = fftw_plan_many_dft_r2c(rank, n, howmany,
                                  rphi, inembed,
                                  rstride, idist,
                                  cphi, onembed,
                                  cstride, odist,
                                  FFTW_FORWARD);
	return;
}

void fft_free(void) {
	fftw_destroy_plan(r2c);
	fftw_destroy_plan(c2r);
	return;
}
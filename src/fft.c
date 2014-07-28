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
int main(void) {
	int i,j;
	Nx = 32; Ny=16;	
	cphi = (double complex *)malloc(sizeof(double complex)*Nx*Ny);
	rphi = (double *)cphi;
	double *x,*y;
	
	for(i=0;i<Nx;i++) x[i]=-15 + i*(30.0/Nx);
	for(i=0;i<Ny;i++) y[i]=-15 + i*(30.0/Ny);
	
	for(i=0;i<Nx;i++) {
		for(j=0;j<Ny;j++) {
			rphi[j+Ny*i] = 	-1.0/sqrt(x[i]*x[i]+y[j]*y[j]+.6*.6);	
		}
	}
	
	for(i=0;i<Nx;i++) {
		for(j=0;j<Ny/2+1;j++) {
			printf("%lg ", cphi[j+Ny*i]);
		}
	}
	
	fftw_execute(r2c);

	free(x); free(y);
	free(cphi); free(rphi);
	return;
}
void fft_init(void) {

	int rank = 1; /* not 2: we are computing 1d transforms */
    int n[] = {Nx}; /* 1d transforms of length 10 */
	int howmany = Ny;
    int idist = 1,odist = 1; 
    int istride =Ny , ostride = Ny; /* distance between two elements in
                                      the same column */
    int *inembed = n, *onembed = n;
#ifdef OPENMP	
	fftw_plan_with_nthreads( 1 );
#endif
	
	
	c2r = fftw_plan_many_dft_c2r(rank, n, howmany,
                                  cphi, inembed,
                                  istride, idist,
                                  rphi, onembed,
                                  ostride, odist,
                                  FFTW_BACKWARD);

	r2c = fftw_plan_many_dft_r2c(rank, n, howmany,
                                  rphi, inembed,
                                  istride, idist,
                                  cphi, onembed,
                                  ostride, odist,
                                  FFTW_FORWARD);
	return;
}

void fft_free(void) {
	fftw_destroy_plan(r2c);
	fftw_destroy_plan(c2r);
	return;
}
#include "planetdisk.h"





void global_c2r(double *y, Field *fld) {
/* Used to pack the complex u,v,sig arrays into the real gsl y array  */
	
	memcpy(&y[0],(double *)&fld->u[istart],sizeof(double complex)*Nx*NC);
	memcpy(&y[Nx*NR],(double *)&fld->v[istart],sizeof(double complex)*Nx*NC);
	memcpy(&y[2*Nx*NR],(double *)&fld->sig[istart],sizeof(double complex)*Nx*NC);
	
	return;


}
void global_c2r_dt(double *y, Field *fld) {
/* Used to pack the complex u,v,sig arrays into the real gsl y array  */
	
	memcpy(&y[0],(double *)fld->dtu,sizeof(double complex)*Nx*NC);
	memcpy(&y[Nx*NR],(double *)fld->dtv,sizeof(double complex)*Nx*NC);
	memcpy(&y[2*Nx*NR],(double *)fld->dtsig,sizeof(double complex)*Nx*NC);
	
	return;


}
void global_r2c(const double *y, Field *fld) {
/* Used to unpack the real gsl y array into the complex u,v,sig arrays */

	memcpy(&fld->u[istart],(double complex *)&y[0],sizeof(double complex)*Nx*NC);
	memcpy(&fld->v[istart],(double complex *)&y[Nx*NR],sizeof(double complex)*Nx*NC);
	memcpy(&fld->sig[istart],(double complex *)&y[2*Nx*NR],sizeof(double complex)*Nx*NC);
	return;
}
void calc_deriv(double complex *in, double complex *dxout, double complex *dyout
					, double dx, double *k) {
/* Calculate the derivatives of in; both x and "y"
	bc contains the boundary condition for the left boundary
		an odd b.c means that the real part is antisymmetric about 0 and the imaginary
		part is even about 0, i.e u(-x) = -conj(u(x))
		opposite for an even b.c, i.e u(-x) = conj(u(x))
	If one of the output arrays is given as NULL, then it is ignored.
	
		The result is added to dxout and dyout arrays.
*/
	int i,d;
	double complex temp;
#ifdef OPENMP 
	#pragma omp parallel private(i,temp,d) shared(k,in,dxout,dyout,dx) num_threads(NUMTHREADS)
	#pragma omp for schedule(static)
#endif	
	for(i=istart;i<iend;i++) {
			if (dyout != NULL) dyout[i] += I*k[i-istart]*in[i];
			if (dxout != NULL) {
				temp = 0;
				for(d=0;d<2*NG;d++) {
					temp += in[i+NC*deriv.ind[d]]*deriv.coeffs[d];

				}						
				dxout[i] += temp/dx;
			}
	}				

	
	return;	
		
}

void print_time(double t) {
	int hr, min;	
	hr = (int)floor(t/(60.*60.)); 
	t -= hr*60*60;	
	min = (int)floor(t/60);
	t -= min*60;
	
	
	if (hr==0) {
		if (min == 0) {
			printf("Total Runtime:\t%.3lgs\n",t);
			
		}
		else {
			printf("Total Runtime:\t%dm%.3lgs\n",min,t);	
		}
	}
	else {
		printf("Total Runtime:\t%dh%dm%.3lgs\n",hr,min,t);
	}
	return;
}

void apply_filter(Field *fld) {
/* Apply a raised cosine filter to the solution */
	int i,j;
	double sigk;
	double fac=2*M_PI*Ny;
	
	for(i=NG;i<Nx+NG;i++) {
		for(j=0;j<Nmax;j++) {
			sigk = .5*(1+cos(j*fac));
			fld->u[CINDX] *= sigk;
			fld->v[CINDX] *= sigk;
			fld->sig[CINDX] *= sigk;	
		}
	}

	return;
}
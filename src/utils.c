#include "planetdisk.h"


void global_c2r(double *y, Field *fld) {
/* Used to pack the complex u,v,sig arrays into the real gsl y array  */
	
	memcpy(&y[0],(double *)fld->u,sizeof(double complex)*NTOTC);
	memcpy(&y[NTOTR],(double *)fld->v,sizeof(double complex)*NTOTC);
	memcpy(&y[2*NTOTR],(double *)fld->sig,sizeof(double complex)*NTOTC);
	
	return;


}
void global_c2r_dt(double *y, Field *fld) {
/* Used to pack the complex u,v,sig arrays into the real gsl y array  */
	
	memcpy(&y[0],(double *)fld->dtu,sizeof(double complex)*NTOTC);
	memcpy(&y[NTOTR],(double *)fld->dtv,sizeof(double complex)*NTOTC);
	memcpy(&y[2*NTOTR],(double *)fld->dtsig,sizeof(double complex)*NTOTC);
	
	return;


}
void global_r2c(const double *y, Field *fld) {
/* Used to unpack the real gsl y array into the complex u,v,sig arrays */

	memcpy(fld->u,(double complex *)&y[0],sizeof(double complex)*NTOTC);
	memcpy(fld->v,(double complex *)&y[NTOTR],sizeof(double complex)*NTOTC);
	memcpy(fld->sig,(double complex *)&y[2*NTOTR],sizeof(double complex)*NTOTC);
	return;
}
void calc_deriv(double complex *in, double complex *dxout, double complex *dyout
					, double dx, double *k, char *lbc, double complex *rbc) {
/* Calculate the derivatives of in; both x and "y"
	bc contains the boundary condition for the left boundary
		an odd b.c means that the real part is antisymmetric about 0 and the imaginary
		part is even about 0
		opposite for an even b.c
	If one of the output arrays is given as NULL, then it is ignored.
	
		The result is added to dxout and dyout arrays.
*/
	int i,j,d,ind;
	
	if (dyout != NULL) {
		for(i=0;i<NTOTC;i++) dyout[i]+= I*k[i]*in[i];
	}
	if (dxout != NULL) {
		for(i=0;i<Nx;i++) {
			for(j=0;j<NC;j++) { 
				for(d=0;d<2*NG;d++) {
					ind = i + deriv.ind[d];
					if (ind < 0) {		// Left B.C 
						if (strcmp("odd",lbc)==0)
							dxout[CINDX] += -conj(in[j-NC*(1+ind)])*deriv.coeffs[d];
						else 
							dxout[CINDX] += conj(in[j-NC*(1+ind)])*deriv.coeffs[d];
					}
					else if (ind > Nx-1) {		// Right B.C
						dxout[CINDX] += rbc[j]*deriv.coeffs[d];
					}
					else {
						dxout[CINDX] += in[j+NC*ind]*deriv.coeffs[d];
					}
	
				}						
				dxout[CINDX] /= dx;
			}
		}				
	}
	
	return;	
		
}

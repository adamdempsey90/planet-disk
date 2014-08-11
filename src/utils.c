#include "planetdisk.h"


void global_c2r(double *y, Field *fld) {
/* Used to pack the complex u,v,sig arrays into the real gsl y array  */
	
	memcpy(&y[0],(double *)&fld->u[NG*NC],sizeof(double complex)*Nx*NC);
	memcpy(&y[Nx*NR],(double *)&fld->v[NG*NC],sizeof(double complex)*Nx*NC);
	memcpy(&y[2*Nx*NR],(double *)&fld->sig[NG*NC],sizeof(double complex)*Nx*NC);
	
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

	memcpy(&fld->u[NG*NC],(double complex *)&y[0],sizeof(double complex)*Nx*NC);
	memcpy(&fld->v[NG*NC],(double complex *)&y[Nx*NR],sizeof(double complex)*Nx*NC);
	memcpy(&fld->sig[NG*NC],(double complex *)&y[2*Nx*NR],sizeof(double complex)*Nx*NC);
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
	int i,j,d;
	double complex temp;

	for(i=NG;i<Nx+NG;i++) {
		for(j=0;j<NC;j++) { 
			if (dyout != NULL) dyout[CINDX] += I*k[j]*in[CINDX];
			if (dxout != NULL) {
				temp = 0;
				for(d=0;d<2*NG;d++) {
					temp += in[j+NC*(i+deriv.ind[d])]*deriv.coeffs[d];

				}						
				dxout[CINDX] += temp/dx;
			}
		}
	}				

	
	return;	
		
}

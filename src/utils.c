#include "meanwave.h"
void r_2_c(const double *y, double complex *cy, const int Nx) {
	int i,j,m;
#ifdef OPENMP
	#pragma omp parallel private(i,j,m) shared(cy,y) num_threads(NUMTHREADS)
	#pragma omp for schedule(static)
#endif	
	for (j=0;j<Nx+2*NG;j++) {
		for(m=0;m<NK;m++) {
			for(i=0;i<3;i++) {
				if(m==0)	cy[CINDX(i,0)] = y[i+3*m+3*(2*NK-1)*j];
				else	cy[CINDX(i,m)]=y[i+3*(2*m-1)+3*(2*NK-1)*j]
										+ I*y[i+3*2*m+3*(2*NK-1)*j];
			}
		}
	}
	return;
}
void c_2_r(double *y, const double complex *cy, const int Nx) {

	int i,j,m;

#ifdef OPENMP
	#pragma omp parallel private(i,j,m) shared(cy,y) num_threads(NUMTHREADS)
	#pragma omp for schedule(static)
#endif
	for (j=0;j<Nx+2*NG;j++) {
		for(m=0;m<NK;m++) {
			for(i=0;i<3;i++) {
				if (m==0) y[i+3*m+3*(2*NK-1)*j]=creal(cy[CINDX(i,m)]);
				else {
					y[i+3*(2*m-1)+3*(2*NK-1)*j]=creal(cy[CINDX(i,m)]);
					y[i+3*2*m+3*(2*NK-1)*j]=cimag(cy[CINDX(i,m)]);
				}
			}
		}
	}
	return;
}

void calc_derivs(parameters *p, int flag) {

	int i,j,m,d,n;

#ifdef OPENMP
	#pragma omp parallel private(j,i,m,d,n) shared(p) num_threads(NUMTHREADS)
	#pragma omp for schedule(static)
#endif
	for(j=NG;j<p->Nx+NG;j++) {
	for(m=0;m<NK;m++) {
	if (flag==0) {      // Just do cy
		for(i=0;i<3;i++) {
			p->dxcy[CINDX(i,m)]=0;
			for(d=0;d<2*NG;d++) {
				 p->dxcy[CINDX(i,m)]+=(p->cy[i+3*m+3*NK*(j+deriv_ind[d])])*deriv_coeffs[d];
			}
			p->dxcy[CINDX(i,m)] /= (p->dx);
		}	
	}
	if (flag==1) {      // Just do stress and mass flux
		for(i=0;i<2;i++) {
			for(n=i;n<2;n++) {
				DXPIP(i,n,m) =  0;
				DXPP(i,n,m) = 0;
				for(d=0;d<2*NG;d++) {
					DXPIP(i,n,m) += 
								  (p->VT->p[i][n][m+NK*(j+deriv_ind[d])])*deriv_coeffs[d];  
					DXPP(i,n,m) += 
										(p->VT->pp[i][n][m+NK*(j+deriv_ind[d])])*deriv_coeffs[d];  
				
				}
				DXPIP(i,n,m) /= (p->dx);
				DXPP(i,n,m) /= (p->dx);
			}
		}	
		DXMF(m) = 0; 
		for(d=0;d<2*NG;d++) {
			DXMF(m) += (p->MF->mf[m+NK*(j+deriv_ind[d])])*deriv_coeffs[d];
		}
		DXMF(m) /= (p->dx); 
	}
	}
	}
	

	return;

}

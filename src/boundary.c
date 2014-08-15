#include "planetdisk.h"

void set_bc(Field *fld) {
/* Set the boundary conditions for the right boundary */
	int i,indx;
//	double sig0 = (fld->Params->sig0);
	

/* Zero Gradient B.C */	

	
	for(i=0;i<NG;i++) {
		indx = NC*(2*NG-(1+i));
		memcpy(&(fld->u[i*NC]),&(fld->u[indx]),sizeof(double complex)*NC);
		memcpy(&(fld->v[i*NC]),&(fld->v[indx]),sizeof(double complex)*NC);
		memcpy(&(fld->sig[i*NC]),&(fld->sig[indx]),sizeof(double complex)*NC);

		indx = NC*(2*(NG+Nx)-(i+Nx+NG+1));
		memcpy(&(fld->u[(i+Nx+NG)*NC]),&(fld->u[indx]),sizeof(double complex)*NC);
		memcpy(&(fld->v[(i+Nx+NG)*NC]),&(fld->v[indx]),sizeof(double complex)*NC);
		memcpy(&(fld->sig[(i+Nx+NG)*NC]),&(fld->sig[indx]),sizeof(double complex)*NC);
	}	
		
//	 for(i=0;i<NG;i++) {
// 		for(j=0;j<NC;j++) {
// 			fld->u[CINDX] = 0;
// 			fld->v[CINDX] = 0;
// 			
// 			if(j==0) fld->sig[CINDX] = sig0;
// 			else	fld->sig[CINDX] = 0;
// 			
// 
// 		}
// 	}

	return;
}

void wavekillbc(Field *fld,double dt)
{
	int i;
	double R,tau,x;
	double x_inf = -(fld->Params->Lx)*.5*.95;
	double x_sup = (fld->Params->Lx)*.5*0.95;
// #ifdef OPENMP 
// 	#pragma omp parallel private(i,x,R,tau) shared(fld) num_threads(NUMTHREADS)
// 	#pragma omp for schedule(static)
// #endif	
	for(i=istart;i<iend;i++) {
		x = fld->xx[i];
		R=0;
		if (x > x_sup) R = (x-x_sup)/(fld->x[Nx+NG-1] - x_sup);
		if (x < x_inf) R = (x_inf - x)/(x_inf - fld->x[NG]);

		R *= R;
		tau = 2*M_PI/(30*fld->Params->omega);

		if (R>0.0) {
			tau /= R;
			if (fld->kk[i-istart]==0) {
#ifdef BACKEVOLVE	
				fld->u[i] = (fld->u[i])/(1+dt);
				fld->v[i] = (fld->v[i])/(1+dt);
				fld->sig[i] = ((fld->sig[i])*tau + (fld->Params->sig0)*dt)/(dt+tau);
#endif
			}
			else {
#ifdef WAVEEVOLVE	
					fld->u[i] = (fld->u[i])/(1+dt);
					fld->v[i] = (fld->v[i])/(1+dt);
					fld->sig[i]=(fld->sig[i])/(1+dt);
#endif
			}
		}
	
	}
	
	return;
}
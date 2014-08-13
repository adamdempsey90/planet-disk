#include "planetdisk.h"

void set_bc(Field *fld) {
/* Set the boundary conditions for the right boundary */
	int i,j, indx;
	double sig0 = (fld->Params->sig0);
	

/* Zero B.C */	
	for(i=0;i<NG;i++) {
		for(j=0;j<NC;j++) {
			indx = j+NC*(NG-1 + NG-i);
			fld->u[CINDX] = -conj(fld->u[indx]);
			fld->v[CINDX] = -conj(fld->v[indx]);
			fld->sig[CINDX] = conj(fld->sig[indx]);
		}
	}	
	
	for(i=Nx+NG;i<Nx+2*NG;i++) {
		for(j=0;j<NC;j++) {
			fld->u[CINDX] = 0;
			fld->v[CINDX] = 0;
			
			if(j==0) fld->sig[CINDX] = sig0;
			else	fld->sig[CINDX] = 0;
			

		}
	}

	return;
}

void wavekillbc(Field *fld,double dt)
{
	int i,j;
	double R,tau,x;
	double x_sup = fld->x[Nx+NG-1] - (fld->Params->Lx)*0.05;
#ifdef OPENMP 
	#pragma omp parallel private(i,j,x,R,tau) shared(fld) num_threads(NUMTHREADS)
	#pragma omp for schedule(static)
#endif	
	for(i=NG;i<(Nx+NG);i++) {
		x = fld->x[i];
		for(j=0;j<NC;j++) { 
			R=0;
			if (x > x_sup) R = (x-x_sup)/(fld->x[Nx+NG-1] - x_sup);
			R *= R;
			tau = 2*M_PI/(10*fld->Params->omega);
			if (R>0.0) {
				tau /= R;
				if (j==0) {
#ifdef BACKEVOLVE	
				fld->u[CINDX] = (fld->u[CINDX])/(1+dt/tau);
				fld->v[CINDX] = (fld->v[CINDX])/(1+dt/tau);
				fld->sig[CINDX] = ((fld->sig[CINDX])*tau + 1.0*dt)/(dt+tau);
#endif
				}
				else {
#ifdef WAVEEVOLVE	
					fld->u[CINDX] = (fld->u[CINDX])/(1+dt/tau);
					fld->v[CINDX] = (fld->v[CINDX])/(1+dt/tau);
					fld->sig[CINDX]=(fld->sig[CINDX])/(1+dt/tau);
#endif
				}
			}
		}
	}
	
	return;
}
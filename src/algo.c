#include "planetdisk.h"


void func(double t, double complex *y, double complex *f,Field *fld) {

	func_calls++;
//	printf("\t\t # Function calls = %d \n",func_calls);
//	printf("Function Called at time = %lg...\n",t);

/*	Get the new u,v,sig */	
//	printf("Copy Data...\n");
	
	y_2_fld(fld,y);
	
	zero_derivs(fld);
	set_bc(fld);
/* Calculate the RHS */	
//	printf("Fill RHS...\n");

	fill_rhs(fld,t);

	
	memcpy(&f[0],fld->dtu,sizeof(double complex)*Nx*NC);
	memcpy(&f[Nx*NC],fld->dtv,sizeof(double complex)*Nx*NC);
	memcpy(&f[2*Nx*NC],fld->dtsig,sizeof(double complex)*Nx*NC);
//	output_rhs(fld); output_pi(fld);
/* Copy complex data into real array with u,v,sig all combined for gsl */
//	printf("Copy Data to GSL...\n");

// 	memcpy(&f[0],(double *)fld->dtu,sizeof(double complex)*NTOTC);
// 	memcpy(&f[NTOTR],(double *)fld->dtv,sizeof(double complex)*NTOTC);
// 	memcpy(&f[2*NTOTR],(double *)fld->dtsig,sizeof(double complex)*NTOTC);
	
//	printf("Finished Function Call...\n");

	return; 


}
void fill_rhs(Field *fld, double t) {
/* 		Fill the RHS of the EOM */
	int i;
	double qom = (fld->Params->q)*(fld->Params->omega);
	double om2 = 2*(fld->Params->omega);
	double k;
	double complex phi, dxphi;
	
/* Fill the derivative arrays */

	calc_deriv(fld->u,fld->dxu,fld->dyu,fld->Params->dx,fld->kk);
	calc_deriv(fld->v,fld->dxv,fld->dyv,fld->Params->dx,fld->kk);
	calc_deriv(fld->sig,fld->dxsig,fld->dysig,fld->Params->dx,fld->kk);

/*	Fill RHS arrays with any non-convolution terms*/	

#ifdef OPENMP 
	#pragma omp parallel private(i,k,phi,dxphi) shared(fld) num_threads(NUMTHREADS)
	#pragma omp for schedule(static)
#endif
	for(i=istart;i<iend;i++) {
			k = fld->kk[i-istart];			
#ifndef BACKEVOLVE
			if( k != 0) {
#endif
#ifndef WAVEEVOLVE 
			if( k == 0 ) {
#endif
			if (k < kmax) {
				phi = calc_pot(fld->phi[i],t,fld->Params->tau);
				dxphi = calc_pot(fld->dxphi[i],t,fld->Params->tau);

				fld->dtu[i-istart] =  om2*(fld->v[i]) - dxphi;

				fld->dtv[i-istart] = (qom-om2)*(fld->u[i]) -I*k*phi;
				fld->dtsig[i-istart] = 0;
				fld->dtu[i-istart] += qom*(fld->dyu[i])*(fld->xx[i-istart]);
				fld->dtv[i-istart] += qom*(fld->dyv[i])*(fld->xx[i-istart]);
				fld->dtsig[i-istart] += qom*(fld->dysig[i])*(fld->xx[i-istart]);
			}
			else {
				fld->dtu[i-istart] = 0;
				fld->dtv[i-istart] = 0;
				fld->dtsig[i-istart]=0;
			}
#ifndef BACKEVOLVE
			}
#endif
#ifndef WAVEEVOLVE 
			}
#endif
	
	}
	
/* Start adding the convolutions. 
   Start with advection and mass flux terms */

	convolve(&fld->u[istart],&fld->dxu[istart],fld->dtu,-1);
	convolve(&fld->v[istart],&fld->dyu[istart],fld->dtu,-1);
	
	convolve(&fld->u[istart],&fld->dxv[istart],fld->dtv,-1);
	convolve(&fld->v[istart],&fld->dyv[istart],fld->dtv,-1);

	convolve(&fld->sig[istart],&fld->dxu[istart],fld->dtsig,-1);
	convolve(&fld->dxsig[istart],&fld->u[istart],fld->dtsig,-1);
	convolve(&fld->sig[istart],&fld->dyv[istart],fld->dtsig,-1);
	convolve(&fld->dysig[istart],&fld->v[istart],fld->dtsig,-1);
	
/* Add viscosity and pressure */
	
	add_visc(fld);

/* Make sure everything is zeroed that's supposed to be */	
#ifndef BACKEVOLVE 

	for(i=0;i<Nx;i++) {
		fld->dtu[NC*i] = 0;
		fld->dtv[NC*i] = 0;
		fld->dtsig[NC*i] = 0;
	}
#endif
#ifndef WAVEEVOLVE

	for(i=NC;i<Nx*NC;i++) {
		for(j=1;j<NC;j++) {	
			fld->dtu[CINDX] = 0;
			fld->dtv[CINDX] = 0;
			fld->dtsig[CINDX] = 0;
		}
	}

#endif

// 	for(i=0;i<Nx*NC;i++) {
// 		if (fld->kk[i+istart] >= kmax) {
// 			fld->dtu[i] = 0;
// 			fld->dtv[i] = 0;
// 			fld->dtsig[i] = 0;
// 		}
// 	}
	return;	
}	


double complex calc_pot(double complex phi,double t, double tau) {
	
	if (tau==0) return phi;
	else return (1-exp(-t/tau))*phi;
	
}
void zero_derivs(Field *fld) {
	int i;

#ifdef OPENMP 
	#pragma omp parallel private(i) shared(fld) num_threads(NUMTHREADS)
	#pragma omp for schedule(static)
#endif	
	for(i=0;i<NTOTC;i++) {
		fld->dxu[i] = 0;
		fld->dxv[i]= 0;
		fld->dxsig[i] = 0;
		fld->dyu[i] = 0;
		fld->dyv[i]= 0;
		fld->dysig[i] = 0;
		fld->Tens->divPix[i] = 0;
		fld->Tens->divPiy[i] = 0;
		if (i<Nx*NC) {
			fld->dtu[i] = 0;
			fld->dtv[i] = 0;
			fld->dtsig[i] = 0;
		}
	}
	return;
}
void shear_advection(Field *fld,double dt) {
	int i;
	double qom = (fld->Params->q)*(fld->Params->omega);
	dt /= 2;
#ifdef OPENMP 
	#pragma omp parallel private(i) shared(fld) num_threads(NUMTHREADS)
	#pragma omp for schedule(static)
#endif	
	for(i=istart;i<iend;i++) {
		fld->u[i] *= (1+I*(fld->kk[i])*qom*(fld->xx[i-istart])*dt)	
						/(1-I*(fld->kk[i-istart])*qom*(fld->xx[i-istart])*dt);
		fld->v[i] *= (1+I*(fld->kk[i])*qom*(fld->xx[i-istart])*dt)
						/(1-I*(fld->kk[i-istart])*qom*(fld->xx[i-istart])*dt);
		fld->sig[i] *= (1+I*(fld->kk[i])*qom*(fld->xx[i-istart])*dt)
						/(1-I*(fld->kk[i-istart])*qom*(fld->xx[i-istart])*dt);
// 		fld->u[i] *= (1+I*(fld->kk[i-istart])*qom*(fld->xx[i-istart])*dt);
// 		fld->v[i] *= (1+I*(fld->kk[i-istart])*qom*(fld->xx[i-istart])*dt);
// 		fld->sig[i] *= (1+I*(fld->kk[i-istart])*qom*(fld->xx[i-istart])*dt);
	}

	

	return;
}

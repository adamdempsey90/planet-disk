#include "planetdisk.h"


int func (double t, const double y[], double f[],void *params) {

//	printf("Function Called at time = %lg...\n",t);
	Field *fld = (Field *)params;

/*	Get the new u,v,sig */	
//	printf("Copy Data...\n");

	global_r2c(y, fld);

/* Calculate the RHS */	
//	printf("Fill RHS...\n");

	fill_rhs(fld,t);

/* Copy complex data into real array with u,v,sig all combined for gsl */
//	printf("Copy Data to GSL...\n");

// 	memcpy(&f[0],(double *)fld->dtu,sizeof(double complex)*NTOTC);
// 	memcpy(&f[NTOTR],(double *)fld->dtv,sizeof(double complex)*NTOTC);
// 	memcpy(&f[2*NTOTR],(double *)fld->dtsig,sizeof(double complex)*NTOTC);

	global_c2r_dt(f,fld);
	
//	printf("Finished Function Call...\n");

	return  GSL_SUCCESS; 


}
void fill_rhs(Field *fld, double t) {
/* 		Fill the RHS of the EOM */
	int i;
	double qom = (fld->Params->q)*(fld->Params->omega);
	double om2 = 2*(fld->Params->omega);
	double k;

/* Fill the derivative arrays */
//	printf("\tCalculating Derivatives...\n");

	calc_deriv(fld->u,fld->dxu,fld->dyu,fld->Params->dx,fld->k,"odd",fld->ubc);
	calc_deriv(fld->v,fld->dxv,fld->dyv,fld->Params->dx,fld->k,"odd",fld->vbc);
	calc_deriv(fld->sig,fld->dxsig,fld->dysig,fld->Params->dx,fld->k,"even",fld->sigbc);

/*	Fill RHS arrays with any non-convolution terms*/	
//	printf("\tAdding Non-Convolution Terms...\n");

	for(i=0;i<NTOTC;i++) {
	
#ifndef BACKEVOLVE
		if( fmod(i/(float)NC,1) != 0) {
#endif
#ifndef WAVEEVOLVE 
		if( fmod(i/(float)NC,1) == 0 ) {
#endif

		k = fld->k[i];
		fld->dtu[i] = qom*I*k*(fld->u[i])*(fld->x[i]); + om2*(fld->v[i]) 
						- calc_pot(fld->dxphi[i],t,fld->Params->tau);
		fld->dtv[i] = qom*I*k*(fld->v[i])*(fld->x[i]) +(qom-om2)*(fld->u[i])
						-I*k*calc_pot(fld->phi[i],t,fld->Params->tau);
		fld->dtsig[i] = qom*I*k*(fld->x[i])*(fld->sig[i]);
		
#ifndef BACKEVOLVE
		}
#endif
#ifndef WAVEEVOLVE 
		}
#endif
	}
	
/* Start adding the convolutions. */
//	printf("\tAdding Convolution Terms...\n");

	convolve(fld->u,fld->dxu,fld->dtu,-1);
	convolve(fld->v,fld->dyu,fld->dtu,-1);
	
	convolve(fld->u,fld->dxv,fld->dtv,-1);
	convolve(fld->v,fld->dyv,fld->dtv,-1);

	convolve(fld->sig,fld->dxu,fld->dtsig,-1);
	convolve(fld->dxsig,fld->u,fld->dtsig,-1);
	convolve(fld->sig,fld->dyv,fld->dtsig,-1);
	convolve(fld->dysig,fld->v,fld->dtsig,-1);
	
/* Add viscosity and pressure */
//	printf("\tAdding Viscosity...\n");
	
	add_visc(fld);
	

	
	return;	
}	


double complex calc_pot(double complex phi,double t, double tau) {
	
	if (tau==0) return phi;
	else return (1-exp(-t/tau))*phi;
	
}


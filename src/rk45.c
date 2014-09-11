#include "rk45.h"

double complex *y, *yerr, *f, *oldy;


int rk45_step_apply(Field *fld, double *t, double *h) {
	int status;
	double tol = fld->Params->tol;
	double oldh;
	fld_2_y(fld,oldy);
	do {
	
		memcpy(y,oldy,sizeof(double complex)*rk_size);
		oldh = *h;
		rk45_step(y, yerr, f, *t, *h, fld);
		status = new_h(yerr,h,tol);
	
		if (*h < MIN_STEP) return -1;
		
	} while (status != 0);
	
	*t += oldh;
	y_2_fld(fld,y);

	
	
//	printf("new h = %lg\n", *h);

	return 1;

}

int new_h(double complex *yerr, double *h, double tol) {
	int i;
	double r;
	double peps= DBL_MIN;
	double eps;

	
	for(i=0;i<rk_size;i++) {
		peps = fmax(peps,fabs(yerr[i]));
	}

	MPI_Allreduce(&peps,&eps,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

	eps /= tol;
	if (eps>1.1) {
		r = SAFETY*pow(eps,-1.0/(rk_order));
		if (r < 0.2) r=0.2;
		*h *= r;
//	printf("%d\trel_eps = %lg,  r = %lg,   new h = %.12e\n",rank,eps,r,*h);

		return 1;
	}
	else if (eps < .5) {
		r = SAFETY*pow(eps,-1.0/(rk_order+1.0));
//		if (r > 5) r=5;
		if (r < 1) r=1;
		*h *= r;
//	printf("%d\trel_eps = %lg,  r = %lg,   new h = %.12e\n",rank,eps,r,*h);

		return 0;
	}
	else {
		r = 1;
//	printf("%d\trel_eps = %lg,  r = %lg,   new h = %.12e\n",rank,eps,r,*h);

		return 0;
	}
	
}



void fld_2_y(Field *fld, double complex *q) {

	memcpy(&q[0],&(fld->u[istart]),sizeof(double complex)*Nx*NC);
	memcpy(&q[Nx*NC],&(fld->v[istart]),sizeof(double complex)*Nx*NC);
	memcpy(&q[2*Nx*NC],&(fld->sig[istart]),sizeof(double complex)*Nx*NC);
	
	return;
}
void y_2_fld(Field *fld, double complex *q) {

	memcpy(&(fld->u[istart]),&q[0],sizeof(double complex)*Nx*NC);
	memcpy(&(fld->v[istart]),&q[Nx*NC],sizeof(double complex)*Nx*NC);
	memcpy(&(fld->sig[istart]),&q[2*Nx*NC],sizeof(double complex)*Nx*NC);
	
	return;
}


void init_rk45(void) {
	
	rk_size = Nx*NC*3;
	
	y = (double complex *)malloc(sizeof(double complex)*rk_size);
	f = (double complex *)malloc(sizeof(double complex)*rk_size);
	yerr = (double complex *)malloc(sizeof(double complex)*rk_size);
 	oldy = (double complex *)malloc(sizeof(double complex)*rk_size);
 	
	rk45_step_init();
	
	return;

}

void free_rk45(void) {
	free(yerr); free(oldy); free(y); free(f);
	rk45_step_free();
	return;
}
#include "planetdisk.h"

#define MIN_STEP 1e-7
#define SAFETY .9

#ifdef RKCK
#define NRK 6
#define RKORDER 4
static double a[NRK]={0,0.2, .3, .6, 1., .875};
static double b[NRK][NRK-1] = {
			{ 0,0,0,0,0},
			{ .2, 0,0,0,0},
			{3./40,9./40,0,0,0},
			{.3,-.9,1.2,0,0},
			{-11./54,2.5,-70./27,35./27,0},
			{1631./55296,175./512,575./13824,44275./110592,253./4096} };	
static double c5[NRK]={37./378,0,250./621,125./594,0,512./1771};
static double c4[NRK]={2825./27648,0,18575./48384,13525./55296,277./14366,.25};
#endif


#ifdef RKF
#define NRK 6
#define RKORDER 4
static double a[NRK]={0,0.25, .375, 12./13, 1., .5};
static double b[NRK][NRK-1] = {
			{ 0,0,0,0,0},
			{ .25, 0,0,0,0},
			{3./32,9./32,0,0,0},
			{1932./2197,-7200./2197,7296./2197,0,0},
			{439./216,-8.,3680./513,-845./4104,0},
			{-8./27,2.,-3544./2565,1859./4104,-11./40} };	
static double c5[NRK]={16./135,0,6656./12825, 28561./56340,-9./50,2./55};
static double c4[NRK]={25./216,0,1408./2565,2197./4104,-.2,0};
#endif


#ifdef RKDP
#define NRK 7
#define RKORDER 5
static double a[NRK]={0,.2,.3,.8,8./9,1.,1.};
static double b[NRK][NRK-1] = {
			{ 0,0,0,0,0,0},
			{ .2, 0,0,0,0,0},
			{3./40,9./40,0,0,0,0},
			{44./45,-56./15,32./9,0,0,0},
			{19372./6561,-25360./2187,64448./6561,-212./729,0,0},
			{9017./3168,-355./33,46732./5247,49./176,-5103./18656,0},
			{35./384,0,500./1113,125./192,-2187./6784,11./84} };	
static double c5[NRK]={5179./57600,0,7571./16695,393./640,-92097./339200,187./2100,1./40};
static double c4[NRK]={35./384,0,500./1113,125./192,-2187./6784,11./84,0};
#endif


double complex *d[NRK];
double complex *y4, *oldy, *y, *f;
int total;

int new_h(double *h, double eps);
double rk45_step(Field *fld, double t, double h); 

// void f2y(Field *fld, double complex *y);
// void y2f(Field *fld, double complex *y);


int rk45_step_apply(Field *fld, double *t, double *h) {
	int status;
	double eps;
	double tol = fld->Params->tol;
	double oldh;
	fld_2_y(fld,oldy);
	do {
		oldh = *h;
		eps = rk45_step(fld,*t,*h);
		printf("eps=%.12e \n", eps);
		status = new_h(h,eps/tol);

		if (*h < MIN_STEP) return -1;
	} while (status != 0);
	
//	printf("EXIT LOOP\n");
	*t += oldh;
	y_2_fld(fld,y4);
//	printf("new h = %lg\n", *h);

	return 1;

}

int new_h(double *h, double eps) {
	double r;
	

	if (eps>1.1) {
		r = SAFETY *pow(eps,-1.0/(RKORDER));
		if (r < 0.2) r=0.2;
		*h *= r;
	printf("rel_eps = %.12e,  r = %lg,   new h = %.12e\n",eps,r,*h);

		return 1;
	}
	else if (eps < .5) {
		r = SAFETY*pow(eps,-1.0/(RKORDER+1.0));
		if (r > 5) r=5;
		if (r < 1) r=1;
		*h *= r;
	printf("rel_eps = %.12e,  r = %lg,   new h = %.12e\n",eps,r,*h);

		return 0;
	}
	else {
		r = 1;
	printf("rel_eps = %.12e,  r = %lg,   new h = %.12e\n",eps,r,*h);

		return 0;
	}
	
}


// double rk45_step(Field *fld, double t, double h) {
// 	int i,j,m, indx;
// 	int inc = Nx*NC;
// 	double eps = 0;
// 	int total = 0;
// 	memcpy(y5,oldy,sizeof(double complex)*Nx*NC*3);
// 
// 	for(j=0;j<6;j++) {
// 		func(t+a[j]*h,fld);
// 		
// 		for(i=0;i<Nx*NC;i++) {
// 		
// 			indx = i;
// 			d[j][indx] = fld->dtu[i]*h;
// 			y5[indx] += d[j][indx]*c5[j];
// 			eps += cabs((c5[j]-c4[j])*d[j][indx]);
// 			total++;
// 			if (j<5) {
// 				fld->u[i+istart] = oldy[indx];
// 				for (m=0;m<=j;m++) {
// 					fld->u[i+istart] += b[j+1][m]*d[m][indx];
// 				}
// 			}
// 			else {
// 				fld->u[i+istart] = y5[indx];
// 			}
// 			
// 			indx += inc;
// 			d[j][indx] = fld->dtv[i]*h;
// 			y5[indx] += d[j][indx]*c5[j];
// 			eps += cabs((c5[j]-c4[j])*d[j][indx]);
// 			total++;
// 			if (j<5) {
// 				fld->v[i+istart] = oldy[indx];
// 				for (m=0;m<=j;m++) {
// 					fld->v[i+istart] += b[j+1][m]*d[m][indx];
// 				}
// 			}
// 			else {
// 				fld->v[i+istart] = y5[indx];
// 			}
// 			
// 			indx += inc;
// 			d[j][indx] = fld->dtsig[i]*h;
// 			y5[indx] += d[j][indx]*c5[j];
// 			eps += cabs((c5[j]-c4[j])*d[j][indx]);
// 			total++;
// 			if (j<5) {
// 				fld->sig[i+istart] = oldy[indx];
// 				for (m=0;m<=j;m++) {
// 					fld->sig[i+istart] += b[j+1][m]*d[m][indx];
// 				}
// 			}	
// 			else {
// 				fld->sig[i+istart] = y5[indx];
// 			}		
// 				
// 
// 	
// 		}
// 	
// 	}
// 
// 	return eps/total;
// }

double rk45_step(Field *fld, double t, double h) {
	int i,j,m;
	double eps = DBL_MIN;
	double complex err;
	memcpy(y4,oldy,sizeof(double complex)*total);
	memcpy(y,oldy,sizeof(double complex)*total);

	for(j=0;j<NRK;j++) {
		func(t+a[j]*h,fld,y,f);

// #ifdef OPENMP 
// 	#pragma omp parallel private(i) shared(fld,y,y4,d) num_threads(NUMTHREADS)
// 	#pragma omp for schedule(static) reduction(+:eps)
// #endif

		for(i=0;i<total;i++) {
			
			d[j][i] = f[i]*h;
			y4[i] += d[j][i]*c4[j];
			err = (c5[j]-c4[j])*d[j][i];
			eps = fmax(eps,fmax(fabs(creal(err)),fabs(cimag(err))));
			if (j<NRK-1) {
				y[i] = oldy[i];
				for (m=0;m<=j;m++) {
					y[i] += b[j+1][m]*d[m][i];
				}
			}
		}
	
	}

	return eps;
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

// void f2y(Field *fld, double complex *y) {
// 	
// 	memcpy(&y[0], &(fld->u[istart]),sizeof(double complex)*Nx*NC);
// 	memcpy(&y[Nx*NC],&(fld->v[istart]),sizeof(double complex)*Nx*NC);
// 	memcpy(&y[2*Nx*NC],&(fld->sig[istart]),sizeof(double complex)*Nx*NC);
// 	
// 	return;
// 
// }
// 
// void y2f(Field *fld, double complex *y) {
// 
// 	memcpy(&(fld->u[istart]),&y[0], sizeof(double complex)*Nx*NC);
// 	memcpy(&(fld->v[istart]),&y[Nx*NC], sizeof(double complex)*Nx*NC);
// 	memcpy(&(fld->sig[istart]),&y[2*Nx*NC], sizeof(double complex)*Nx*NC);
// 	
// 	return;
// }

void init_rk45(void) {
	int i;
	total = Nx*NC*3;
	for(i=0;i<NRK;i++) d[i] = (double complex *)malloc(sizeof(double complex)*total);
	
	y = (double complex *)malloc(sizeof(double complex)*total);
	f = (double complex *)malloc(sizeof(double complex)*total);
	y4 = (double complex *)malloc(sizeof(double complex)*total);
 	oldy = (double complex *)malloc(sizeof(double complex)*total);

	return;

}

void free_rk45(void) {
	int i;
	free(y4); free(oldy); free(y); free(f);
	for(i=0;i<NRK;i++) free(d[i]);
	return;

}
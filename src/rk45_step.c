#include "rk45.h"

static const double a[7] = {0,.2,.3,.8,8./9,1.,1.};
static const double b1 = .2;
static const double b2[6]	= {3./40,9./40,0,0,0,0};
static const double b3[6]	= {44./45,-56./15,32./9,0,0,0};
static const double b4[6]	= {19372./6561,-25360./2187,64448./6561,-212./729,0,0};
static const double b5[6]	= {9017./3168,-355./33,46732./5247,49./176,-5103./18656,0};	
static const double b6[6]	= {35./384,0,500./1113,125./192,-2187./6784,11./84};	

static double c[7] = {5179./57600,0,7571./16695,393./640,-92097./339200,187./2100,1./40};

static double ec[7] = {-(71./57600), 0, 71./16695, -(71./1920), 17253./339200, -(22./525), 1./40};

double complex *k1,*k2,*k3,*k4,*k5,*k6;
double complex *ytmp;

void rk45_step( double complex *y, double complex *yerr, double complex *f, 
									double t, double h, Field *fld) 
{
/*	Runge-Kutta Dormand-Prince Method */

	int i;
	
	
	func(t,y,f,fld);



	for(i=0;i<rk_size;i++) {	
		k1[i] = f[i];
		ytmp[i] = y[i] + h*k1[i]*b1;
	}
	
	func(t+a[1]*h,ytmp,f,fld);

	for(i=0;i<rk_size;i++) {
		k2[i] = f[i];
		ytmp[i] =  y[i] + h*(k1[i]*b2[0]+k2[i]*b2[1]);
	}
	
	func(t+a[2]*h,ytmp,f,fld);


	for(i=0;i<rk_size;i++) {
		k3[i] = f[i];
		ytmp[i] = y[i] + h*(k1[i]*b3[0] +k2[i]*b3[1] + k3[i]*b3[2]);
	}
	
	func(t+a[3]*h,ytmp,f,fld);

	for(i=0;i<rk_size;i++) {
		k4[i] = f[i];
		ytmp[i] = y[i] + h*(k1[i]*b4[0] + k2[i]*b4[1] + k3[i]*b4[2] + k4[i]*b4[3]);
	}

	func(t+a[4]*h,ytmp,f,fld);
	

	for(i=0;i<rk_size;i++) {
		k5[i] = f[i];
		ytmp[i] = y[i] + h*(k1[i]*b5[0] + k2[i]*b5[1] + k3[i]*b5[2] 
											+ k4[i]*b5[3] + k5[i]*b5[4]);
	}
	
	func(t+a[5]*h,ytmp,f,fld);

	for(i=0;i<rk_size;i++) {
		k6[i] = f[i];
		ytmp[i] = y[i] + h*(k1[i]*b6[0] + k2[i]*b6[1] + k3[i]*b6[2] 
											+ k4[i]*b6[3] + k5[i]*b6[4] + k6[i]*b6[5]);
	}
	
	func(t+a[6]*h,ytmp,f,fld);

	for(i=0;i<rk_size;i++) {									
		y[i] += h*(k1[i]*c[0] + k2[i]*c[1] + k3[i]*c[2] 
									+ k4[i]*c[3] + k5[i]*c[4] + k6[i]*c[5]+f[i]*c[6]);
		yerr[i] = h*(k1[i]*ec[0] + k2[i]*ec[1] + k3[i]*ec[2] 
									+ k4[i]*ec[3] + k5[i]*ec[4] + k6[i]*ec[5]+f[i]*c[6]);
	}

	return;
	

}

void rk45_step_init(void) {

	k1 = (double complex *)malloc(sizeof(double complex)*rk_size);
	k2 = (double complex *)malloc(sizeof(double complex)*rk_size);
	k3 = (double complex *)malloc(sizeof(double complex)*rk_size);
	k4 = (double complex *)malloc(sizeof(double complex)*rk_size);
	k5 = (double complex *)malloc(sizeof(double complex)*rk_size);
	k6 = (double complex *)malloc(sizeof(double complex)*rk_size);
	ytmp = (double complex *)malloc(sizeof(double complex)*rk_size);
	rk_order = 5;
	return;
}

void rk45_step_free(void) {
	free(k1); free(k2); free(k3); 
	free(k4); free(k5);	free(k6);
	free(ytmp);
	return;
}
#include "rk45.h"


static const double a[5] = {0.25, .375, 12./13, 1., .5};
static const double b1 = .25;
static const double b2[2]	= {3./32,9./32};
static const double b3[3]	= {1932./2197,-7200./2197,7296./2197};
static const double b4[4]	= {439./216,-8.,3680./513,-845./4104};
static const double b5[5]	= {-8./27,2.,-3544./2565,1859./4104,-11./40};	

static double c[6] = {16./135, 0, 6656./12825, 28561./56430, -9./50, 2./55};

static double ec[6] = {1./360, 0, -(128./4275), -(2197./75240), 1./50, 2./55};

double complex *k1,*k2,*k3,*k4,*k5;
double complex *ytmp;

void rk45_step( double complex *y, double complex *yerr, double complex *f, 
									double t, double h, Field *fld) 
{
/*	Runge-Kutta-Fehlberg Method */

	int i;
	
	
	func(t,y,f,fld);
	for(i=0;i<rk_size;i++) {	
		k1[i] = f[i];
		ytmp[i] = y[i] + h*k1[i]*b1;
	}
	
	func(t+a[0]*h,ytmp,f,fld);
	
	for(i=0;i<rk_size;i++) {
		k2[i] = f[i];
		ytmp[i] =  y[i] + h*(k1[i]*b2[0]+k2[i]*b2[1]);
	}
	
	func(t+a[1]*h,ytmp,f,fld);
	
	for(i=0;i<rk_size;i++) {
		k3[i] = f[i];
		ytmp[i] = y[i] + h*(k1[i]*b3[0] +k2[i]*b3[1] + k3[i]*b3[2]);
	}
	
	func(t+a[2]*h,ytmp,f,fld);
	
	for(i=0;i<rk_size;i++) {
		k4[i] = f[i];
		ytmp[i] = y[i] + h*(k1[i]*b4[0] + k2[i]*b4[1] + k3[i]*b4[2] + k4[i]*b4[3]);
	}

	func(t+a[3]*h,ytmp,f,fld);
	
	for(i=0;i<rk_size;i++) {
		k5[i] = f[i];
		ytmp[i] = y[i] + h*(k1[i]*b5[0] + k2[i]*b5[1] + k3[i]*b5[2] 
											+ k4[i]*b5[3] + k5[i]*b5[4]);
	}
	
	func(t+a[4]*h,ytmp,f,fld);
	
	for(i=0;i<rk_size;i++) {
		y[i] += h*(k1[i]*c[0] + k2[i]*c[1] + k3[i]*c[2] 
									+ k4[i]*c[3] + k5[i]*c[4] + f[i]*c[5]);
		yerr[i] = h*(k1[i]*ec[0] + k2[i]*ec[1] + k3[i]*ec[2] 
									+ k4[i]*ec[3] + k5[i]*ec[4] + f[i]*ec[5]);
	}

	return;
	

}

void rk45_step_init(void) {

	k1 = (double complex *)malloc(sizeof(double complex)*rk_size);
	k2 = (double complex *)malloc(sizeof(double complex)*rk_size);
	k3 = (double complex *)malloc(sizeof(double complex)*rk_size);
	k4 = (double complex *)malloc(sizeof(double complex)*rk_size);
	k5 = (double complex *)malloc(sizeof(double complex)*rk_size);
	ytmp = (double complex *)malloc(sizeof(double complex)*rk_size);
	rk_order = 5;
	return;
}

void rk45_step_free(void) {
	free(k1); free(k2); free(k3); 
	free(k4); free(k5);	
	free(ytmp);
	return;
}
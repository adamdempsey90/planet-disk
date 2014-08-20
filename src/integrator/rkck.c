#include "rk45.h"

static const double a[5] = {0.2, .3, .6, 1., .875};
static const double b1 = .2;
static const double b2[2]	= {3./40,9./40};
static const double b3[3]	= {.3,-.9,1.2};
static const double b4[4]	= {-11./54,2.5,-70./27,35./27};
static const double b5[5]	= {1631./55296,175./512,575./13824,44275./110592,253./4096};	

static double c[6] = {37./378,0,250./621,125./594,0,512./1771};

static double ec[6] = { -(277./64512), 0, 6925./370944, -(6925./202752), -(277./14366), 277./7084};

double complex *k1,*k2,*k3,*k4,*k5;
double complex *ytmp;

void rk45_step( double complex *y, double complex *yerr, double complex *f, 
									double t, double h, Field *fld) 
{

/*	Runge-Kutta Cash-Karp Method */

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
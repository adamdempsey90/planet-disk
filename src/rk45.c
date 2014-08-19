#include "planetdisk.h"

#define MIN_STEP 1e-7

double a[6]={0,0.2, .3, .6, 1, .875};
double b[6][5] = {
			{ 0,0,0,0,0},
			{ .2, 0,0,0,0},
			{3./40,9./40,0,0,0},
			{.3,-.9,1.2,0,0},
			{-11./54,2.5,-70./27,35./27,0},
			{1631./55296,175./512,575./13824,44275./110592,253./4096} };
			
// double b1 = .2;			
// double b2[2] =	{3./40,9./40};
// double b3[3] = {.3,-.9,1.2};
// double b4[4] =	{-11./54,2.5,-70./27,35./27};
// double b5[5] =	{1631./55296,175./512,575./13824,44275./110592,253./4096};			
			
double c5[6]={37./378,0,250./621,125./594,0,512./1771};
double c4[6]={2825./27648,0,18575./48384,13525./55296,277./14366,.25};


double complex *d[6];
double complex *y5, *oldy;


void new_h(double *h, double eps);
double rk45_step(Field *fld, double t, double h); 
void f2y(Field *fld, double complex *y);
void y2f(Field *fld, double complex *y);
void init_rk45(void);
void free_rk45(void);

int rk45_apply_step(Field *fld, double *t, double *h) {
	double eps;
	double tol = fld->Params->tol;
	
	do {
	
		eps = rk45_step(fld,*t,*h);
		
		if (eps > tol) {
			new_h(h,eps/tol);
			y2f(fld,oldy);
		}
		if (*h < MIN_STEP) return -1;
	} while (eps > tol);
	
	
	*t += *h;
	new_h(h,tol/eps);

	return 1;

}
void func(double t, Field *fld) {
	return;
}

void new_h(double *h, double eps) {

	if (eps >= 1)	*h *= pow(eps,0.2);
	else	*h *= pow(eps,0.25);

	return;
}


double rk45_step(Field *fld, double t, double h) {
	int i,j,m;
	int inc = Nx*NC;
	double eps = 0;
	int total = Nx*NC*3;
	f2y(fld,oldy);
	memcpy(y5,oldy,sizeof(double complex)*Nx*NC*3);

	for(j=0;j<5;j++) {
		func(t+a[j]*h,fld);
		
		for(i=0;i<Nx*NC;i++) {

			d[j][i] = fld->dtu[i]*h;
			y5[i] += d[j][i]*c5[j];
			eps += cabs((c5[j]-c4[j])*d[j][i]);
			fld->u[i+istart] = oldy[i];
			for(m=0;m<j;m++) fld->u[i+istart] += b[j+1][m]*d[m][i];
			
			d[j][i+inc] = fld->dtv[i]*h;
			y5[i+inc] += d[j][i+inc]*c5[j];
			eps += cabs((c5[j]-c4[j])*d[j][i+inc]);
			fld->v[i+istart] = oldy[i+inc];
			for(m=0;m<j;m++) fld->v[i+istart] += b[j+1][m]*d[m][i+inc];
			
			d[j][i+2*inc] = fld->dtsig[i]*h;
			y5[i+2*inc] += d[j][i+2*inc]*c5[j];
			eps += cabs((c5[j]-c4[j])*d[j][i+2*inc]);
			fld->sig[i+istart] = oldy[i+2*inc];
			for(m=0;m<j;m++) fld->sig[i+istart] += b[j+1][m]*d[m][i+2*inc];
	
		}
	
	}

	return eps;
}


void f2y(Field *fld, double complex *y) {
	
	memcpy(&y[0], &(fld->u[istart]),sizeof(double complex)*Nx*NC);
	memcpy(&y[Nx*NC],&(fld->v[istart]),sizeof(double complex)*Nx*NC);
	memcpy(&y[2*Nx*NC],&(fld->sig[istart]),sizeof(double complex)*Nx*NC);
	
	return;

}

void y2f(Field *fld, double complex *y) {

	memcpy(&(fld->u[istart]),&y[0], sizeof(double complex)*Nx*NC);
	memcpy(&(fld->v[istart]),&y[Nx*NC], sizeof(double complex)*Nx*NC);
	memcpy(&(fld->sig[istart]),&y[2*Nx*NC], sizeof(double complex)*Nx*NC);
	
	return;
}

void init_rk45(void) {
	int i;
	
	for(i=0;i<6;i++) d[i] = (double complex *)malloc(sizeof(double complex)*Nx*NC);
	
	
	y5 = (double complex *)malloc(sizeof(double complex)*Nx*NC*3);
 	oldy = (double complex *)malloc(sizeof(double complex)*Nx*NC*3);

	return;

}

void free_rk45(void) {
	int i;
	free(y5); free(oldy); 
	for(i=0;i<6;i++) free(d[i]);
	return;

}
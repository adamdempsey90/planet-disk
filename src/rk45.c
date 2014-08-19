#include "planetdisk.h"

double a[6]={0,0.2, .3, .6, 1, .875};
// double b[6][5] = {
// 			{ 0,0,0,0,0},
// 			{ .2, 0,0,0,0},
// 			{3./40,9./40,0,0,0},
// 			{.3,-.9,1.2,0,0},
// 			{-11./54,2.5,-70./27,35./27,0},
// 			{1631./55296,175./512,575./13824,44275./110592,253./4096} };
			
double b1 = .2;			
double b2[2] =	{3./40,9./40};
double b3[3] = {.3,-.9,1.2};
double b4[4] =	{-11./54,2.5,-70./27,35./27};
double b5[5] =	{1631./55296,175./512,575./13824,44275./110592,253./4096};			
			
double c5[6]={37./378,0,250./621,125./594,0,512./1771};
double c4[6]={2825./27648,0,18575./48384,13525./55296,277./14366,.25};


double complex *d0,*d1,*d2,*d3,*d4,*d5;
double complex *ku, *kv, *ks;
double complex *y4, *y5, *y45, *y0;
double newt;






void rk45_step(Field *fld,double *t, double *h) {
/* Runge-Kutta Cash-Karp 4/5 adaptive step method.
	Butcher Tableau taken from Numerical Recipes 2nd Ed.
*/
	int i;
	double dt=h;
	memcpy(&y0[0],fld->u,sizeof(double complex)*NTOTC);
	memcpy(&y0[NTOTC],fld->v,sizeof(double complex)*NTOTC);
	memcpy(&y0[2*NTOTC],fld->sig,sizeof(double complex)*NTOTC);
	
	func(*t,fld);
	for(i=0;i<Nx*NC;i++) {
	
		indx = i
		d0[indx] = fld->dtu[i]*h;
		y4[indx] = y0[i+istart]+d0[indx]*c4[0];
		y5[indx] = y0[i+istart]+d0[indx]*c5[0];
		y45[indx] = cabs(y5[indx]-y4[indx]);
		fld->u[i+istart] = y0[i+istart] + b1*d0[indx];
	
		indx += NTOTC
		d0[indx] = fld->dtv[i]*h;
		y4[indx] = y0[indx+istart]*d0[indx]*c4[0];
		y5[indx] = y0[indx+istart]*d0[indx]*c5[0];
		y45[indx] = y5[indx]-y4[indx];
		fld->v[i] = y0[indx+istart] + b1*d0[indx];

		indx += NTOTC
		d0[indx] = fld->dtsig[i]*h;
		y4[indx] = y0[indx+istart]*d0[indx]*c4[0];
		y5[indx] = y0[indx+istart]*d0[indx]*c5[0];
		y45[indx] = y5[indx]-y4[indx];
		fld->sig[i] = y0[indx+istart] + b1*d0[indx];


	}
	
	func(*t + h*a[1],fld);
	
	for(i=0;i<Nx*NC;i++) {
	
		indx = i
		d1[indx] = fld->dtu[i]*h;
		y4[indx] += d1[indx]*c4[1];
		y5[indx] += d1[indx]*c5[1];
		y45[indx] += cabs(y5[indx]-y4[indx]);
		fld->u[i+istart] += y0[i+istart] + b1*d0[indx];
	
		indx += NTOTC
		d0[indx] = fld->dtv[i]*h;
		y4[indx] = y0[indx+istart]*d0[indx]*c4[0];
		y5[indx] = y0[indx+istart]*d0[indx]*c5[0];
		y45[indx] = y5[indx]-y4[indx];
		fld->v[i] = y0[indx+istart] + b1*d0[indx];

		indx += NTOTC
		d0[indx] = fld->dtsig[i]*h;
		y4[indx] = y0[indx+istart]*d0[indx]*c4[0];
		y5[indx] = y0[indx+istart]*d0[indx]*c5[0];
		y45[indx] = y5[indx]-y4[indx];
		fld->sig[i] = y0[indx+istart] + b1*d0[indx];


	}
	
}
















double rk45_step(Field *fld,double *t,double h) {
	int i,j;
	double t0 = *t;
	indx = j+Nx*NC*i
	for(i=istart;i<iend;i++) {
		for(j=0;j<6;j++) {
			indx = j+Nx*NC*(i-istart);
			add_y(j,fld)
			func(t0+a[j]*h,yn);
			
		}
		
	}

}
void add_y(int j, Field *fld) {
	int i,m;
	
	for(i=0;i<Nx*NC;i++) {
		yn[i] = fld->u[i+istart];
		yn[i+Nx*NC] = fld->v[i+istart];
		yn[i+2*Nx*NC] = fld->sig[i+istart];
		for(m=1;m<j;m++) {
			yn[i] += b[j][m-1]*k[m-1 + Nx*NC*i];
			yn[i+Nx*NC] += b[j][m-1]*k[m-1 + Nx*NC*(i+Nx*NC)];
			yn[i+2*Nx*NC] += b[j][m-1]*k[m-1 + Nx*NC*(i+2*Nx*NC)];
		}	
	}
	return;
}
void init_rk45(void) {

	d1 = (double complex *)malloc(sizeof(double complex)*Nx*NC*3);
	kv = (double complex *)malloc(sizeof(double complex)*Nx*NC*6);
	ks = (double complex *)malloc(sizeof(double complex)*Nx*NC*6);
	
	yn = (double complex *)malloc(sizeof(double complex)*Nx*NC*3);
	
	return;

}

void free_rk45(void) {

	free(ku); free(kv); free(ks);
	free(yn);
	return;

}
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


int rk45_step_apply(Field *fld, double *t, double *h) {
	double eps;
	double tol = fld->Params->tol;
	
//	printf("%lg %lg \n", *t, *h);
	do {
		f2y(fld,oldy);
		eps = rk45_step(fld,*t,*h);
//		printf("eps=%.12e \n", eps);
		if (eps > tol) {
			new_h(h,tol/eps);
//			printf("new h = %lg \n",*h);
			y2f(fld,oldy);
		}
		if (*h < MIN_STEP) return -1;
	} while (eps > tol);
	
//	printf("EXIT LOOP\n");
	*t += *h;
	new_h(h,tol/eps);
//	printf("new h = %lg\n", *h);

	return 1;

}

void new_h(double *h, double eps) {
	double S = .98;
	if (eps >= 1)	*h *= (2-S)*pow(eps,0.2);
	else	*h *= S*pow(eps,0.25);

	return;
}


double rk45_step(Field *fld, double t, double h) {
	int i,j,m, indx;
	int inc = Nx*NC;
	double eps = 0;
	int total = 0;
	memcpy(y5,oldy,sizeof(double complex)*Nx*NC*3);

	for(j=0;j<6;j++) {
		func(t+a[j]*h,fld);
		
		for(i=0;i<Nx*NC;i++) {
		
			indx = i;
			d[j][indx] = fld->dtu[i]*h;
			y5[indx] += d[j][indx]*c5[j];
			eps += cabs((c5[j]-c4[j])*d[j][indx]);
			total++;
			if (j<5) {
				fld->u[i+istart] = oldy[indx];
				for (m=0;m<=j;m++) {
					fld->u[i+istart] += b[j+1][m]*d[m][indx];
				}
			}
			else {
				fld->u[i+istart] = y5[indx];
			}
			
			indx += inc;
			d[j][indx] = fld->dtv[i]*h;
			y5[indx] += d[j][indx]*c5[j];
			eps += cabs((c5[j]-c4[j])*d[j][indx]);
			total++;
			if (j<5) {
				fld->v[i+istart] = oldy[indx];
				for (m=0;m<=j;m++) {
					fld->v[i+istart] += b[j+1][m]*d[m][indx];
				}
			}
			else {
				fld->v[i+istart] = y5[indx];
			}
			
			indx += inc;
			d[j][indx] = fld->dtsig[i]*h;
			y5[indx] += d[j][indx]*c5[j];
			eps += cabs((c5[j]-c4[j])*d[j][indx]);
			total++;
			if (j<5) {
				fld->sig[i+istart] = oldy[indx];
				for (m=0;m<=j;m++) {
					fld->sig[i+istart] += b[j+1][m]*d[m][indx];
				}
			}	
			else {
				fld->sig[i+istart] = y5[indx];
			}		
				

	
		}
	
	}

	return eps/total;
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
	
	for(i=0;i<6;i++) d[i] = (double complex *)malloc(sizeof(double complex)*Nx*NC*3);
	
	
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <fftw3.h>

fftw_plan r2c, c2r;

#define NC (Ny/2+1)
#define NR 2*NC
#define NTOTR Nx*NR
#define NTOTC Nx*NC

#define CINDX (j+i*NC)
#define RINDX (j+i*NR)

typedef struct Field {
	
	double *vx, *vy, *x, *y, *k;
	double complex *u, *v;

	int Nx, Ny, Nk;
	
} Field;

void free_field(Field *fld);
void fft_free(void);
void init_fft(int Nx, int Ny);
void init_Field(Field *fld);
void convolve(double complex *q1, double complex *q2, double complex *res, int Nx, int Ny);

double complex *wc1, *wc2, *wc3;
double *wr1, *wr2, *wr3;
double *mask;

fftw_plan r2c1, c2r1,r2c2, c2r2,r2c3, c2r3;

int main(void) {
	int i,j;
	int Nx,Ny,Nk;
	FILE *rfp1, *cfp1, *rfp2, *cfp2;
	Field *fld = (Field *)malloc(sizeof(Field));
	fld->Nx=256;
	fld->Ny =256;
	fld->Nk = (fld->Ny)/2 + 1;
	Nx = fld->Nx;
	Ny = fld->Ny;
	Nk = fld->Nk;

	double complex *uv = (double complex *)malloc(sizeof(double complex)*NTOTC);

	
	
	printf("Initializing Field...\n");
	init_Field(fld);
	printf("Initializing FFT...\n");
	init_fft(fld->Nx,fld->Ny);
	
	printf("Outputting Field...\n");	
	rfp1 = fopen("fftout/reu.dat","w"); cfp1 = fopen("fftout/imu.dat","w");
	rfp2 = fopen("fftout/rev.dat","w"); cfp2 = fopen("fftout/imv.dat","w");
	
	for(i=0;i<Nx;i++) {
		for(j=0;j<NC;j++) {
			fprintf(rfp1,"%lg ",creal(fld->u[CINDX]));
			fprintf(cfp1,"%lg ",cimag(fld->u[CINDX]));
			fprintf(rfp2,"%lg ",creal(fld->v[CINDX]));
			fprintf(cfp2,"%lg ",cimag(fld->v[CINDX]));
		}
		fprintf(rfp1,"\n"); fprintf(cfp1,"\n");
		fprintf(rfp2,"\n"); fprintf(cfp2,"\n");
	}
	fclose(rfp1); fclose(cfp1);
	fclose(rfp2); fclose(cfp2);
	
	printf("Convolving...\n");
	convolve(fld->u,fld->v, uv, fld->Nx, fld->Ny);
	
	printf("Outputting Result...\n");
	
	rfp1 = fopen("fftout/rres.dat","w"); cfp1 = fopen("fftout/imres.dat","w");
	for(i=0;i<Nx;i++) {
		for(j=0;j<NC;j++) {
			fprintf(rfp1,"%lg ",creal(uv[CINDX]));
			fprintf(cfp1,"%lg ",cimag(uv[CINDX]));
			
		}
		fprintf(rfp1,"\n"); fprintf(cfp1,"\n");
	}
	fclose(rfp1); fclose(cfp1);	
	
	
	free(uv);
	fft_free();
	free_field(fld);
	return;

}
void convolve(double complex *q1, double complex *q2, double complex *res, int Nx, int Ny) {
	int i,j;

/* De-alias with 2/3 truncation rule */	
	for(i=0;i<NTOTC;i++) {wc1[i] = q1[i]*mask[i]; wc2[i] = q2[i]*mask[i];}

/* FFT these to real space */	
	fftw_execute(c2r1);
	fftw_execute(c2r2);

/* Form product in real space */
	for(i=0;i<NTOTR;i++) wr3[i] = wr1[i]*wr2[i]/(Ny*Ny);

/* FFT back to complex space */
	fftw_execute(r2c3);
/* Copy output */
	memcpy(res,wc3,sizeof(double complex)*NTOTC);	

	return;
}


void init_Field(Field *fld) {
	int i,j;
	int Nx,Ny;
	Nx = fld->Nx;
	Ny = fld->Ny;
	FILE *fpx, *fpy, *fpk;
	double ld, x, k;
	fld->x = (double *)malloc(sizeof(double)*fld->Nx);
	fld->y = (double *)malloc(sizeof(double)*fld->Ny);
	fld->k = (double *)malloc(sizeof(double)*NC);
	fld->u = (double complex *)malloc(sizeof(double complex)*NTOTC);
	fld->v = (double complex *)malloc(sizeof(double complex)*NTOTC);
	fld->vx = (double *)fld->u;
	fld->vy = (double *)fld->v;
	
	for(i=0;i<fld->Nx;i++) fld->x[i] = -15 + i*30.0/fld->Nx;
	for(i=0;i<fld->Ny;i++) fld->y[i] = -15 + i*30.0/fld->Ny;
	for(i=0;i<NC;i++) fld->k[i] = i*2*M_PI/30.0;

	for(i=0;i<Nx;i++) {
		x = fld->x[i];
		for(j=0;j<NC;j++) {
			k = fld->k[j];
			fld->u[CINDX] = cexp(I*k*k*x*x);
			fld->v[CINDX] = cexp(-I*k*k*x*x);
		}
	}

	fpx = fopen("fftout/x.dat","w"); 
	fpy = fopen("fftout/y.dat","w");
	fpk = fopen("fftout/k.dat","w");
	for(i=0;i<Nx;i++) {
		for(j=0;j<Ny;j++) {
			fprintf(fpx,"%lg ",fld->x[i]);
			fprintf(fpy,"%lg ",fld->y[j]);
			if (j<NC)	fprintf(fpk, "%lg ",fld->k[j]);
		}
		fprintf(fpx,"\n");
		fprintf(fpy,"\n");
		fprintf(fpk,"\n");
	}
	fclose(fpx); fclose(fpy); fclose(fpk);
	return;
}
void init_fft(int Nx, int Ny) {
	int i,j,Nmax;
	wc1 = (double complex *)malloc(sizeof(double complex)*NTOTC);
	wr1 = (double *)wc1;
	wc2 = (double complex *)malloc(sizeof(double complex)*NTOTC);
	wr2 = (double *)wc2;
	wc3 = (double complex *)malloc(sizeof(double complex)*NTOTC);
	wr3 = (double *)wc3;
	mask = (double *)malloc(sizeof(double)*NTOTC);
	
	Nmax = (2./3)*(NC-1);
	
	for(i=0;i<Nx;i++) {
		for(j=0;j<NC;j++) {
			if (j < Nmax) mask[CINDX] = 1;
			else mask[CINDX] = 0;
		}
	}
	
	c2r1=fftw_plan_many_dft_c2r(1,&Nx,Ny,wc1,&Nx,1,NC,wr1,&Nx,1,NR,FFTW_ESTIMATE);
	if (c2r1 == NULL) printf("Problem with c2r1\n");
	r2c1=fftw_plan_many_dft_r2c(1,&Nx,Ny,wr1,&Nx,1,NR,wc1,&Nx,1,NC,FFTW_ESTIMATE);
	if (r2c1 == NULL) printf("Problem with r2c1\n");
	
	c2r2=fftw_plan_many_dft_c2r(1,&Nx,Ny,wc2,&Nx,1,NC,wr2,&Nx,1,NR,FFTW_ESTIMATE);
	if (c2r2 == NULL) printf("Problem with c2r2\n");
	r2c2=fftw_plan_many_dft_r2c(1,&Nx,Ny,wr2,&Nx,1,NR,wc2,&Nx,1,NC,FFTW_ESTIMATE);
	if (r2c2 == NULL) printf("Problem with r2c2\n");
	
	c2r3=fftw_plan_many_dft_c2r(1,&Nx,Ny,wc3,&Nx,1,NC,wr3,&Nx,1,NR,FFTW_ESTIMATE);
	if (c2r3 == NULL) printf("Problem with c2r3\n");
	r2c3=fftw_plan_many_dft_r2c(1,&Nx,Ny,wr3,&Nx,1,NR,wc3,&Nx,1,NC,FFTW_ESTIMATE);
	if (r2c3 == NULL) printf("Problem with r2c3\n");
	
	
	return;

}

void fft_free(void) {

	fftw_destroy_plan(r2c);
	fftw_destroy_plan(c2r);
	
	free(wc1);free(wc2);free(wc3);
	free(mask);
	return;


}
void free_field(Field *fld) {
	
	free(fld->u); free(fld->v);
	free(fld->x); free(fld->k); free(fld->y);
	free(fld);
	return;

}

#include "defines.h"
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_bessel.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <string.h>
#include <complex.h>
#include <sys/stat.h>

#define NREAL 9
#define NCOMPLEX 6



#ifdef OPENMP
#define NUMTHREADS 8
#endif


#ifdef SECDERIV
	#define NG 1
	const double deriv_coeffs[2] = {-0.5,0.5};
	const int deriv_ind[2] = {-1,1};
#endif

#ifdef FOURTHDERIV
	#define NG 2
	const double deriv_coeffs[4] = {1.0/12, -2.0/3, 2.0/3, -1.0/12};
	const int deriv_ind[4] = {-2, -1, 1, 2};
#endif

#ifdef SIXTHDERIV
	#define NG 3
	const double deriv_coeffs[6] = {-1.0/60, 3.0/20, -0.75, 0.75, -3.0/20, 1.0/60};
	const int deriv_ind[6] = {-3,-2, -1, 1, 2,3};
#endif

#ifdef EIGTHDERIV
	#define NG 4
	const double deriv_coeffs[8] = {1.0/280,-4.0/105, 1.0/5, 
								-4.0/5,4.0/5,-1.0/5,4.0/105,-1.0/280};
	const int deriv_ind[8] = {-4,-3,-2, -1, 1, 2,3,4};
#endif

#ifdef TENTHDERIV
	#define NG 5
	const double deriv_coeffs[10] = {-2.0/2520, 25.0/2520, -150.0/2520, 600.0/2520, -2100.0/2520,  
							2100.0/2520, -600.0/2520,150.0/2520, -25.0/2520, 2.0/2520};
	const int deriv_ind[10] = {-5, -4,-3,-2, -1, 1, 2, 3, 4, 5};
#endif



#define RINDX(i) i+NREAL*j
#define RINDXP(i) i+NREAL*(j+1)
#define RINDXM(i) i+NREAL*(j-1)
#define RINDXM2(i) i+NREAL*(j-2)
#define RINDXP2(i) i+NREAL*(j+2)

#define CINDX(i,m) i+3*m+3*NK*j
#define CINDXP(i,m) i+3*m+3*NK*(j+1)
#define CINDXM(i,m) i+3*m+3*NK*(j-1)
#define CINDXM2(i,m) i+3*m+3*NK*(j-2)
#define CINDXP2(i,m) i+3*m+3*NK*(j+2)


#define SIGNUM(Z)	((Z > 0) - (Z < 0))
#define ISZERO(j) (j>1)

#define UVEL(m)  (p->cy[CINDX(0,m)])
#define VVEL(m) (p->cy[CINDX(1,m)])
#define SIG(m) (p->cy[CINDX(2,m)])
#define VXBAR (p->cy[CINDX(0,0)])
#define VYBAR (p->cy[CINDX(1,0)])
#define DBAR (p->cy[CINDX(2,0)])
#define XCORD (p->x[j])
#define PIP(i,n,m) (p->VT->p[i][n][m+NK*j])
#define PIPP(i,n,m) (p->VT->pp[i][n][m+NK*j])
#define DXUVEL(m)  (p->dxcy[CINDX(0,m)])
#define DXVVEL(m) (p->dxcy[CINDX(1,m)])
#define DXSIG(m) (p->dxcy[CINDX(2,m)])
#define DXVXBAR (p->dxcy[CINDX(0,0)])
#define DXVYBAR (p->dxcy[CINDX(1,0)])
#define DXDBAR (p->dxcy[CINDX(2,0)])
#define DXPIP(i,n,m) (p->VT->dxp[i][n][m+NK*j])
#define DXPP(i,n,m) (p->VT->dxpp[i][n][m+NK*j])
#define DXMF(m) (p->MF->dxmf[m+NK*j])
#define MF(m) (p->MF->mf[m+NK*j])
#define PHI(m) (p->phi[m+NK*(j-NG)])
#define DXPHI(m) (p->dxphi[m+NK*(j-NG)])

typedef struct stress {
	double complex *Txx, *Txy, *Tyy;
	double complex *Pixx, *Pixy, *Piyy;
	double complex *divPix, *divPiy;
	

	
} stress;

		
typedef struct parameters {

	int Nx,Ntot,Ny;
	double Lx,dx,xsoft,c,Mp,nu,q,omega,Ly,kfac;
	double t0,endt,tau;
	int numf;

	double *k;
	double *x, *y;
	double complex *rhs;
	double complex *u, *v, *sig, *dxu, *dxv, *dxsig;
	double complex *dtu, *dtv, *dtsig;
	double complex *phi,*dxphi;
	double *t;
	double lastt;
	
	char restartfname[50];
	stress *VT;
	
	massflux *MF;
	
	int conv_flag;
	double conv_tol;
	
	double xsplit;	
	
}	parameters;


int func (double t, const double y[], double f[],void *params);
void output_func(parameters *p, double ti);
void read_athena(parameters *p);
void read_input(parameters *p);
void c_2_r(double *y,const double complex *cy, const int Nx);
void r_2_c(const double *y,  double complex *cy, const int Nx);
void output_params(parameters *P);
void restart(parameters *p);
void read_params(parameters *P); 
void output_amf(parameters *p, double ti);
void free_params(parameters *p);
void calc_stress( parameters *p );
void init(parameters *p, double *y);
void allocate_params(parameters *p);
void calc_derivs(parameters *p, int flag);
void bounds( parameters *p, double t);
void wavekillbc(parameters *p,double dt);
void calcTwd(parameters *p, double ti);
void initphi(parameters *p);
double complex calc_pot(double complex phi,double t, double tau);
void write_phi(parameters *p);

int func_count;
int NK;
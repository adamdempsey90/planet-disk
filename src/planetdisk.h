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


#define NC (Ny/2 + 1)
#define NR 2*NC
#define NTOTR  Nx*NR
#define NTOTC Nx*NC

#define CINDX j+i*NC
#define RINDX j+i*NR




typedef struct Field {


	Parameters *Params;
	
	Stress *Tens;

	double *x, *y, *k;
	
	double complex *u, *v, *sig;
	double complex *dxu, *dxv, *dxsig;
	double complex *dyu, *dyv, *dysig;
	double complex *dtu, *dtv, *dtsig;
	double complex *phi, *dxphi;
	
	

} Field;
typedef struct Parameters {
	
	int Nx, Ny, Nk;
	double Lx, Ly, dx;
	double c, omega, xs, nu, q, Mp, sig0;
	double t0, endt, tau;
	int numf;
	char restartfname[50];
	
} Parameters;

typedef struct Stress {
	double complex *Txx, *Txy, *Tyy;
	double complex *Pixx, *Pixy, *Piyy;
	double complex *divPix, *divPiy;

} Stress;

		
int Nx, Ny;
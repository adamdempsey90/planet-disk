#include "defines.h"
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_bessel.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <string.h>
#include <complex.h>
#include <sys/stat.h>



#ifdef SECDERIV
	#define NG 1
#endif

#ifdef FOURTHDERIV
	#define NG 2
#endif

#ifdef SIXTHDERIV
	#define NG 3
#endif

#ifdef EIGTHDERIV
	#define NG 4
#endif

#ifdef TENTHDERIV
	#define NG 5
#endif





#define NC (Ny/2 + 1)
#define NR 2*NC
#define NTOTR  (Nx+2*NG)*NR
#define NTOTC (Nx+2*NG)*NC
#define istart NG*NC
#define iend (NG+Nx)*NC

#define CINDX j+i*NC
#define RINDX j+i*NR



typedef struct Parameters {
	
	int Nx, Ny, Nk;
	double Lx, Ly, dx;
	double h,c, omega, xs, nu, q, Mp, sig0;
	double t0, endt, tau;
	int numf;
	char restartfname[50];
	
} Parameters;

typedef struct Stress {
	double complex *Txx, *Txy, *Tyy;
	double complex *Pixx, *Pixy, *Piyy;
	double complex *divPix, *divPiy;

} Stress;

typedef struct Field {


	Parameters *Params;
	
	Stress *Tens;

	double *x, *y, *k, *kk, *xx;
	double complex *u, *v, *sig;
	double complex *dxu, *dxv, *dxsig;
	double complex *dyu, *dyv, *dysig;
	double complex *dtu, *dtv, *dtsig;
	double complex *phi, *dxphi;
	double *vx, *vy, *dens;
		

} Field;
typedef struct Cons {
	
	double complex *Mx, *My;
	double *rMx, *rMy;

} Cons;

typedef struct Derivative {

	double coeffs[2*NG];
	int ind[2*NG];

} Derivative;
		
int Nx, Ny, outnum, dxoutnum, func_calls, dtoutnum,pioutnum;
Derivative deriv; 

int func (double t, const double y[], double f[],void *params);
void fill_rhs(Field *fld, double t);
double complex calc_pot(double complex phi,double t, double tau);
void convolve(double complex *q1, double complex *q2, double complex *res, double complex mult);
void convolve_inv(double complex *q1, double complex *q2, double complex *res, double complex mult);
void init_fft(void);
void fft_free(void);
void fft_phi(double *rphi, double complex *cphi);
void init(Field *fld);
void allocate_field(Field *fld);
void free_field(Field *fld);
void initphi(Field *fld);
void output(Field *fld);
void output_coords(Field *fld);
void read_input(Field *fld);
void global_c2r(double *y, Field *fld);
void global_r2c(const double *y, Field *fld);
void global_c2r_dt(double *y, Field *fld);
void calc_deriv(double complex *in, double complex *dxout, double complex *dyout
					, double dx, double *k);
void visc_tens(Field *fld);
void add_visc(Field *fld);
void init_derivs(void);
void output_derivs(Field *fld);
void zero_derivs(Field *fld);
void set_bc(Field *fld);
void output_rhs(Field *fld);
void output_pi(Field *fld);
void transform(Field *fld);
void output_reals(Field *fld);
void wavekillbc(Field *fld,double dt);
void output_defines(void);
void restart(Field *fld);
void shear_advection(Field *fld,double dt);
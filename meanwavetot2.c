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

//#define RESTART
//#define COMPVISC
//#define BACKEVOLVE
#define WAVEEVOLVE
#define NONLINEAR
#define NLVISC
//#define OUTGHOST
//#define OUTPHI

//#define SIMPLEVISC

#define OUTBC
//#define PERIODICBC
#define WAVEKILLBC


//#define OUTDERIV

//#define SECDERIV
//#define FOURTHDERIV
//#define SIXTHDERIV
//#define EIGTHDERIV
#define TENTHDERIV

//#define INTFAC

//#define HIGHLEVEL
#define LOWLEVEL

#define NK 25
#define NREAL 9
#define NCOMPLEX 6

#define OPENMP


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
#define DTUVEL(m)  (p->rhs[CINDX(0,m)])
#define DTVVEL(m) (p->rhs[CINDX(1,m)])
#define DTSIG(m) (p->rhs[CINDX(2,m)])
#define DTVXBAR (p->rhs[CINDX(0,0)])
#define DTVYBAR (p->rhs[CINDX(1,0)])
#define DTDBAR (p->rhs[CINDX(2,0)])
#define DXPIP(i,n,m) (p->VT->dxp[i][n][m+NK*j])
#define DXPP(i,n,m) (p->VT->dxpp[i][n][m+NK*j])
#define DXMF(m) (p->MF->dxmf[m+NK*j])
#define MF(m) (p->MF->mf[m+NK*j])
#define PHI(m) (p->phi[m+NK*(j-NG)])
#define DXPHI(m) (p->dxphi[m+NK*(j-NG)])

//#define RESTART
//#define EVOLVEBACK

typedef struct stress {

	double complex *p[2][2], *pp[2][2];
	double complex *dxp[2][2], *dxpp[2][2];
	
} stress;

typedef struct massflux {

	double complex *mf, *dxmf;
} massflux;
		
typedef struct parameters {

	int Nx,Ntot,Ny;
	double Lx,dx,xsoft,c,Mp,nu,q,omega,Ly,kfac;
	double t0,endt,tau;
	int numf;

	double *k;
	double *x, *realrhs;
	double complex *cy, *rhs, *dxcy, *phi, *dxphi;
	double *t;
	double lastt;
	
	char restartfname[50];
	stress *VT;
	
	massflux *MF;
	
	int conv_flag;
	double conv_tol;
	
	double xsplit;	
	
}	parameters;


int func_count;



double phi(double z,double xsoft, double Mp, double k, double t, double tin);
double phip(double z,double xsoft, double Mp, double k, double t, double tin);
double phib(double z, double xsoft, double Mp, double Ly, double t, double tin);
void output_func(parameters *p, double ti);
double calc_max_diff(double *y1, double *y2);
double calc_avg_diff(double *y1, double *y2);
double calc_bc_diff(double *y);
double calc_dtdx(double x, double sig);
double calc_fh(double *y);
void output_rhs(double *x, double *y, double *f, double ti);
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
void remove_shear(parameters *p, double t);
void integrate(double *in, double *out, double dx, int jstart, int jend);
void bounds( parameters *p, double t);
double calc_nu(double x, double Lx);
void write_timestep(double t, double dt);
void advection(double *y, double dt,parameters *p);
void check_convergence(double *intTh, double *intTwb, double *Fp, double Fp0, parameters *p);
void calc_split(parameters *p);
void output_derivs(parameters *p);
void calc_deriv_single(double dx, int Nx, const double *in, double *out);
void wavekillbc(parameters *p,double dt);
void calcTwd(parameters *p, double ti);
void initphi(parameters *p);
double complex calc_pot(double complex phi,double t, double tau);
void write_phi(parameters *p);
void convolution(parameters *p, int j);

int func (double t, const double y[], double f[],void *params)
{

 	parameters *p = (parameters *)params;
	int i,j,m,n;


 	double k;
  	double c=p->c;
//  	double nu=p->nu;
 	double q=p->q;

 	double om=p->omega;
 	double Mp=p->Mp;
 	double xsoft=p->xsoft;
 	double dt = t-p->lastt;
 	
 	
// Copy real array into complex array
 	func_count++; 
 	r_2_c(y,p->cy,p->Nx);
	bounds(p,t);
	

	
	calc_stress(p);
	calc_derivs(p,1);
		
// Do work with complex
// Re-zero RHS array 

	for(i=0;i<(3*NK*(p->Ntot));i++) 	p->rhs[i]=0;
		
#ifdef OPENMP 
	#pragma omp parallel private(j,m,k) shared(p) num_threads(NUMTHREADS)
	#pragma omp for schedule(static)
#endif
 	for(j=NG;j<(p->Nx+NG);j++) { 
	
#ifdef BACKEVOLVE

/* vxbar */
		DTVXBAR	+= -VXBAR*DXVXBAR+2*om*(VYBAR)
				   + DXPIP(0,0,0)/DBAR;
								
//	    DTVXBAR -= calc_pot(DXPHI(m),t,p->tau);
	
/* vybar */
		DTVYBAR += -VXBAR*(DXVYBAR-q*om)-2*om*VXBAR
				   +(1/DBAR)*DXPIP(0,1,0);

/* dbar */
		DTDBAR += -DXDBAR*VXBAR - DBAR*DXVBAR;

#else
		DTVXBAR = 0;
		DTVYBAR = 0;
		DTDBAR = 0;
#endif
		
		for(m=1;m<NK;m++) {
 			k = p->k[m];
#ifdef WAVEEVOLVE

/* u */ 	
			DTUVEL(m) +=-VXBAR*DXUVEL(m)
						-I*k*(VYBAR-q*om*XCORD)*UVEL(m)
						-UVEL(m)*DXVXBAR
						+2*om*VVEL(m);
						
			DTUVEL(m) -= calc_pot(DXPHI(m),t,p->tau);
			
			DTUVEL(m) += (1/DBAR)*(DXPIP(0,0,m)+I*k*PIP(0,1,m)) 
						  - (SIG(m)/(DBAR*DBAR))*DXPIP(0,0,0);
						  
/* v */						  
								
			DTVVEL(m) += -VXBAR*DXVVEL(m)
						 -(VYBAR-q*om*XCORD)*I*k*VVEL(m)
						 -UVEL(m)*(DXVYBAR+(2-q)*om);
						 
			DTVVEL(m) -= I*k*calc_pot(PHI(m),t,p->tau);
			
			DTVVEL(m) += (1/DBAR)*(DXPIP(0,1,m)+I*k*PIP(1,1,m))
						 -(SIG(m)/(DBAR*DBAR))*DXPIP(0,1,0);
/* sigma */

			DTSIG(m) -=( DXDBAR*UVEL(m)+SIG(m)*DXVXBAR+DBAR*DXUVEL(m)+DXSIG(m)*VXBAR
						+ I*k*(SIG(m)*(VYBAR-q*om*XCORD)+DBAR*VVEL(m)));
						
#ifdef BACKEVOLVE

			DTVXBAR += 2*creal(-conj(UVEL(m))*DXUVEL(m)
					   -I*k*conj(VVEL(m))*UVEL(m) 
					   +(conj(SIG(m))*SIG(m)/(DBAR*DBAR*DBAR))*DXPIP(0,0,0)
					   +(DXPP(0,0,m)/DBAR)
					   -(conj(SIG(m))/(DBAR*DBAR))*(DXPIP(0,0,m)
					   +I*k*PIP(0,1,m)));
					   
			DTVYBAR += 2*creal(-conj(UVEL(m))*DXVVEL(m)
						+(conj(SIG(m))*SIG(m)/(DBAR*DBAR*DBAR))*DXPIP(0,1,0)
						+(DXPP(0,1,m)/DBAR)
						-(conj(SIG(m))/(DBAR*DBAR))*(DXPIP(0,1,m)+I*k*PIP(1,1,m)));
					
			DTDBAR -= 2*creal(conj(DXSIG(m))*UVEL(m) + conj(SIG(m))*DXUVEL(m));

#endif

#else
		DTUVEL(m) = 0;
		DTVVEL(m) = 0;
		DTSIG(m) = 0;
#endif	
		}
		
#ifdef WAVEEVOLVE
#ifdef NONLINEAR
		convolution(p,j);				
#endif
#endif
	}
			

		
// Go back to reals 

	c_2_r(f,p->rhs,p->Nx);
 
// 	memcpy((void *)f,(void *)p->realrhs,sizeof(double)*((p->Ntot)*NREAL));	
// 	for(i=0;i<(p->Nx+2)*NREAL;i++) f[i]=p->realrhs[i];
 
 // Ghost Zones
//  	for(i=0;i<n;i++) {
//  		j=0; f[INDX] = 0;
//  		j=Nx+1; f[INDX]= 0;
//  	}
 
  return GSL_SUCCESS;
}


int main (void)
{
	int i, cutind;
	double time=0;
 
  double *y;
  
  func_count=0;
  parameters *P = (parameters *)malloc(sizeof(parameters));
  
  read_input(P);
  printf("Allocating structs...\n");
  allocate_params(P);


  

  y = (double *)malloc((P->Ntot)*(6*NK-3)*sizeof(double));
  
  
  printf("Initialzing...\n");
  init(P,y);
#ifdef RESTART
	restart(P);
#endif
//  read_athena(P);
	printf("Appying BC's...\n");
  bounds(P,0);
	 

  c_2_r(y,P->cy,P->Nx);	
  printf("\tOUTPUT 0\n");
  output_func(P,0);
  
  gsl_odeiv2_system sys = {func, NULL, (6*NK-3)*(P->Ntot), P};


 	
//	gsl_odeiv2_driver_set_hmin(d,0.01);
  
	
	
	printf("Nx=%d, c=%g, Mp=%g, nu=%g, Nk=%d, Lx=%g dx=%g, xsoft=%g, t_end=%g \n",
			P->Nx, P->c, P->Mp, P->nu, NK, P->Lx, 
			P->dx, P->xsoft,P->endt);
	output_params(P);
	

  
	
	printf("Starting loop\n");
//	Evolve with high-level wrapper

#ifdef HIGHLEVEL
	 gsl_odeiv2_driver * d = 
 			 gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
				  1e-8, 1e-4, 0.0);
   for (i = 1; i <= P->numf; i++)
    {
      if (P->conf_flag == 1) {
      	printf("CONVERGENCE DETECTED \n");
      	break;
      }
      P->t[i-1] = P->t0 + i * (P->endt) / ((double) P->numf);
      P->lastt = time;
      int status = gsl_odeiv2_driver_apply (d, &time, P->t[i-1], y);
		printf("finished step %d, at t=%g\n", i, P->t[i-1]);
		printf("\t\t\t deriv count = %d\n", func_count);
		r_2_c(y,P->cy,P->Nx);
#ifdef INTFAC
		remove_shear(P,time-P->lastt);
#endif
		c_2_r(y,P->cy,P->Nx);
		bounds(P,time);
		
//		output_func(y,P,P->t[i-1]);
//		output_amf(P,P->t[i-1]);
      if (status != GSL_SUCCESS)
		{
		  printf ("error, return value=%d\n", status);
		  break;
		}
    }
    gsl_odeiv2_driver_free (d);
#endif

// Evolve with low-level wrapper to see step sizes.
	
	
#ifdef LOWLEVEL
	
	const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;

  	gsl_odeiv2_step * s  = gsl_odeiv2_step_alloc (T, (6*NK-3)*(P->Ntot));
 	 gsl_odeiv2_control * c  = gsl_odeiv2_control_y_new (1e-8, 1e-8);
 	 gsl_odeiv2_evolve * e  = gsl_odeiv2_evolve_alloc ((6*NK-3)*(P->Ntot));
 
 
  int status;
  double	h = .1;
  double 	t=P->t0;
  double 	t1 = P->endt;
  double dt;
  i=1;
  
  
//	output_derivs(P);

  while (t < t1)
    {
      
      dt = t;
      P->lastt = t;
      status = gsl_odeiv2_evolve_apply (e, c, s,
                                           &sys, 
                                           &t, t1,
                                           &h, y);
       dt = t-dt;

#ifdef BACKEVOLVE
#ifdef WAVEKILLBC
       r_2_c(y,P->cy,P->Nx);
	   wavekillbc(P,dt);
	   c_2_r(y,P->cy,P->Nx);
#endif
#endif
#ifdef INTFAC
       r_2_c(y,P->cy,P->Nx);
	   remove_shear(P,dt);
	   c_2_r(y,P->cy,P->Nx);
#endif
	  
//       calc_split(P);	  
	   printf ("step size = %.5e, suggested step size = %.5e, at t=%.5e \n", dt, h,t);
	   printf("\t\t\t deriv count = %d\n", func_count); 
    	printf("\t\t\t CUT at X = %lg \n", P->xsplit);
//	   write_timestep(t,dt);
      if (status != GSL_SUCCESS)
          break;

     
      
   	if( t >= P->t0 + i * (P->endt) / ((double) P->numf)) { 
   		 printf ("\t OUTPUT step size = %.5e, at t=%.5e \n", h,t);
		r_2_c(y,P->cy,P->Nx);

//		c_2_r(y,P->cy,P->Nx);
		calcTwd(P, floor(t));
		output_func(P,floor(t));
//		output_amf(P,floor(t));
      	i++;
     }
    }
  if (status != GSL_SUCCESS) printf("\n\n\tIntegration failed to converge...\n");
  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
	
#endif	
	
  free(y); 
  free_params(P); 
  return 0;
}


void convolution(parameters *p, int j) {
	int m,n;
	double c = (p->c)*(p->c);
	for(m=1;m<NK;m++) {	// k
		for(n=1;n<NK;n++) {	// k'
			if(n!=m) {
				if(abs(n-m)<NK) {
					if (m > n) {
						DTUVEL(m) += -UVEL(n)*DXUVEL(m-n)
									-VVEL(n)*I*(p->k[m-n])*UVEL(m-n)
									+(c/(DBAR*DBAR))*SIG(n)*DXSIG(m-n)
									-(c*DXDBAR/(DBAR*DBAR*DBAR))*SIG(n)*SIG(m-n);
	//							-SIG(n)*(DXPIP(0,0,m-n)+I*(p->k[m-n])*PIP(0,1,m-n))/(DBAR*DBAR);
									
						DTVVEL(m) += -UVEL(n)*DXVVEL(m-n)
									-I*(p->k[m]/2)*VVEL(n)*VVEL(m-n)
									+(c/(DBAR*DBAR))*I*(p->k[m]/2)*SIG(n)*SIG(m-n);
	//							-SIG(n)*(DXPIP(0,1,m-n)+I*(p->k[m-n])*PIP(1,1,m-n))/(DBAR*DBAR);
						DTSIG(m) += -SIG(n)*DXUVEL(m-n)
									-DXSIG(n)*UVEL(m-n)
									-I*(p->k[m])*SIG(n)*VVEL(m-n);
	
					}
					else {
						DTUVEL(m) += -UVEL(n)*conj(DXUVEL(n-m))
									-VVEL(n)*conj(I*(p->k[n-m])*UVEL(n-m))
								+(c/(DBAR*DBAR))*SIG(n)*conj(DXSIG(n-m))
									-(c*DXDBAR/(DBAR*DBAR*DBAR))*SIG(n)*conj(SIG(n-m));
//									-SIG(n)*conj((DXPIP(0,0,n-m)+I*(p->k[n-m])*PIP(0,1,n-m)))/(DBAR*DBAR);	
						DTVVEL(m) += -UVEL(n)*conj(DXVVEL(n-m))
									 -I*(p->k[m]/2)*VVEL(n)*conj(VVEL(n-m))
									 +(c/(DBAR*DBAR))*I*(p->k[m]/2)*SIG(n)*conj(SIG(n-m));
	//							-SIG(n)*conj((DXPIP(0,1,n-m)+I*(p->k[n-m])*PIP(1,1,n-m)))/(DBAR*DBAR);
						DTSIG(m) += -SIG(n)*conj(DXUVEL(n-m))
									-DXSIG(n)*conj(UVEL(n-m))
									-I*(p->k[m])*SIG(n)*conj(VVEL(n-m));	
					}
				}
				if(m+n <NK) {
					DTUVEL(m) += -conj(UVEL(n))*DXUVEL(m+n)
								-conj(VVEL(n))*I*(p->k[m+n])*UVEL(m+n)
								+(c/(DBAR*DBAR))*conj(SIG(n))*DXSIG(m+n)
								-(c*DXDBAR/(DBAR*DBAR*DBAR))*conj(SIG(n))*SIG(m+n);
//							-conj(SIG(n))*(DXPIP(0,0,m+n)+I*(p->k[m+n])*PIP(0,1,m+n))/(DBAR*DBAR);								
					
					DTVVEL(m) += -conj(UVEL(n))*DXVVEL(m+n)
								 -I*(p->k[m]/2)*conj(VVEL(n))*VVEL(m+n)
							 +(c/(DBAR*DBAR))*I*(p->k[m]/2)*conj(SIG(n))*SIG(m+n);
//							-conj(SIG(n))*(DXPIP(0,1,m+n)+I*(p->k[m+n])*PIP(1,1,m+n))/(DBAR*DBAR);					
					DTSIG(m) += -conj(SIG(n))*DXUVEL(m+n)
								-conj(DXSIG(n))*UVEL(m+n)
								-I*(p->k[m])*conj(SIG(n))*VVEL(m+n);

				}
			}
		}	
	}
					
							
	return;									
							
}

double phi(double z,double xsoft, double Mp, double k, double t, double tau) {

	double	xs = sqrt(xsoft*xsoft + z*z);
	
	if (tau == 0) return (-Mp/M_PI) * gsl_sf_bessel_K0(fabs(xs*k));
	else	return (1-exp(-t/tau))*(-Mp/M_PI) * gsl_sf_bessel_K0(fabs(xs*k));

//	return -Mp*cos(z);
}
double phip(double z, double xsoft, double Mp, double k, double t, double tau) {

	double	xs = sqrt(xsoft*xsoft + z*z);
	
	if (tau ==0) return (z/xs)*(Mp/M_PI)*k*gsl_sf_bessel_K1(fabs(xs*k));
	else return (1-exp(-t/tau))*(z/xs)*(Mp/M_PI)*k*gsl_sf_bessel_K1(fabs(xs*k));

//	return Mp*sin(z);
}
double phib(double z, double xsoft, double Mp, double Ly, double t, double tau) {

	double	xs = sqrt(xsoft*xsoft + z*z);
	
	return (1-exp(-t/tau))*(2.0*Mp*z/(xs*xs))/sqrt(4.0*xs*xs+Ly*Ly);

}
void read_input(parameters *p) {
	int i;
	printf("Reading inputs...\n");


	p->Nx=8192; 
	p->Ny = 1024;
	p->Ntot = p->Nx+2*NG;
	p->Lx=120.0; 
	p->Ly=30.0;
	p->dx=(p->Lx)/(p->Nx); 
	p->xsoft = .6;
	p->c=1; 
	p->Mp=.5; 
	p->nu=0.0024; 
	p->q=1.5; 
	p->omega=1.0;
 
	p->kfac = 2.0;
//  	p->k *= p->kfac;

	p->t0 = 0;
	p->tau = 10;
	p->endt = 200;
	p->numf = 200;	

//	sprintf(p->restartfname, "newm0.5nu0.0024_lx60_means.dat");
	sprintf(p->restartfname, "newmeans.dat");

	p->conv_flag = 0;
	
	p->conv_tol = 1e-4;

//	read_params(p);	

	
	
	printf("Finished reading inputs...\n");
	return;
}
void init(parameters *p, double *y) {
	printf("hi\n");
	int i,j,n,m;
	
	printf("k=\n");
	for(i=0;i<NK;i++) {
		p->k[i] = i*(2*M_PI)/(p->Ly);
		printf("\t %lg",p->k[i]);
	}
	printf("\n");
#ifdef OPENMP
	#pragma omp parallel private(i,j,n,m) shared(p,y) num_threads(NUMTHREADS)
{	

	#pragma omp for schedule(static)
#endif
	for(i=0;i<(6*NK*(p->Ntot)-3*(p->Ntot));i++) {
		if (i < 3*NK*(p->Ntot) ) {
			p->rhs[i]=0;
 		}
 		y[i]=0; 
 		p->realrhs[i]=0;
 	 }

#ifdef OPENMP
	#pragma omp for schedule(static)
#endif  
	for(j=0;j<p->Ntot;j++) {
		p->x[j] = -(p->Lx)/2 + (p->dx)*(0.5+j-NG);
		for(m=0;m<NK;m++) {
			if (m==0) {
				VXBAR = 0;
				VYBAR = 0;
				DBAR = 1.0;
				DXVXBAR = 0;
				DXVYBAR = 0;
				DXDBAR = 0;
				PIP(0,0,0) =-(p->c)*(p->c);
				PIP(1,1,0)=-(p->c)*(p->c);
				PIP(0,1,0)= -(p->nu)*(p->q)*(p->omega);
			}
			else {
				UVEL(m) = 0;
				VVEL(m) = 0;
				SIG(m) = 0;
				DXUVEL(m) = 0;
				DXVVEL(m) = 0;
				DXSIG(m) = 0;
				PIP(0,0,m) = 0;
				PIP(1,1,m)= 0;
				PIP(0,1,m)= 0;					
				DXPIP(0,0,m) = 0;
				DXPIP(0,1,m) = 0;
				DXPIP(0,0,m) = 0;				
			}
		
			for(i=0;i<2;i++) {
				p->MF->mf[m+NK*j] = 0;
				DXMF(m) = 0;
				for(n=i;n<2;n++) {
					PIPP(i,n,m) = 0;
					DXPP(i,n,m) = 0;
				}
			}
		}
	
	}		
#ifdef OPENMP
}
#endif
	
	char dir[10];
	for(m=0;m<NK;m++) {
		sprintf(dir,"k%d",m);
		mkdir(dir,0777);
	}
	
	initphi(p);
	return;
}
void r_2_c(const double *y, double complex *cy, const int Nx) {
	int i,j,m;
#ifdef OPENMP
	#pragma omp parallel private(i,j,m) shared(cy,y) num_threads(NUMTHREADS)
	#pragma omp for schedule(static)
#endif	
	for (j=0;j<Nx+2*NG;j++) {
		for(m=0;m<NK;m++) {
			for(i=0;i<3;i++) {
				if(m==0)	cy[CINDX(i,0)] = y[i+3*m+3*(2*NK-1)*j];
				else	cy[CINDX(i,m)]=y[i+3*(2*m-1)+3*(2*NK-1)*j]
										+ I*y[i+3*2*m+3*(2*NK-1)*j];
			}
		}
	}
	return;
}
void c_2_r(double *y, const double complex *cy, const int Nx) {

	int i,j,m;

#ifdef OPENMP
	#pragma omp parallel private(i,j,m) shared(cy,y) num_threads(NUMTHREADS)
	#pragma omp for schedule(static)
#endif
	for (j=0;j<Nx+2*NG;j++) {
		for(m=0;m<NK;m++) {
			for(i=0;i<3;i++) {
				if (m==0) y[i+3*m+3*(2*NK-1)*j]=creal(cy[CINDX(i,m)]);
				else {
					y[i+3*(2*m-1)+3*(2*NK-1)*j]=creal(cy[CINDX(i,m)]);
					y[i+3*2*m+3*(2*NK-1)*j]=cimag(cy[CINDX(i,m)]);
				}
			}
		}
	}
	return;
}
void output_params(parameters *P) {
	FILE *fp;
	fp=fopen("run_params.txt","w");
	
	
	fprintf(fp,"Nx = %d \n", P->Nx);
	fprintf(fp,"Mp = %g \n", P->Mp);
	fprintf(fp,"nu = %g \n", P->nu);
	fprintf(fp,"Lx = %g \n", P->Lx);
	fprintf(fp,"dx = %g \n", P->dx);
	fprintf(fp,"xsoft = %g \n", P->xsoft);
	fprintf(fp,"Nk = %d \n", NK);
	fprintf(fp,"q = %g \n", P->q);
	fprintf(fp,"c = %g \n", P->c);
	fprintf(fp,"omega = %g \n", P->omega);
	fprintf(fp,"Ly = %g \n", P->Ly);
	fprintf(fp,"end time = %g \n", P->endt);
	fprintf(fp,"number of files = %d \n", P->numf);
	
	
	fclose(fp);
	
	return;


}

void output_func(parameters *p, double ti) {

	FILE *fp;
	char fname[50];
	int j,m;

	for(m=0;m<NK;m++) {
		sprintf(fname,"k%d/output-k%d_%g.txt",m,m,ti);
		fp=fopen(fname,"w");

#ifdef OUTGHOST
		for(j=0;j<p->Ntot;j++) 
#else
		for(j=NG;j<p->Nx+NG;j++) 
#endif
		{
			if(m==0)
				fprintf(fp,"%12.8e \t %12.8e \t %12.8e \t %12.8e \n",
					XCORD,creal(VXBAR),creal(VYBAR),creal(DBAR));
			else
				fprintf(fp,"%12.8e \t %12.8e \t %12.8e \t %12.8e \t %12.8e \t %12.8e \t %12.8e \n",
					XCORD,creal(UVEL(m)),cimag(UVEL(m)),creal(VVEL(m)),cimag(VVEL(m)),
					creal(SIG(m)),cimag(SIG(m)));
		}
	
		fclose(fp);
	}
	return;

}
void restart(parameters *p) {
	FILE *fp;
	double ru,iu,rv,iv,rd,id,vx,vy,d,x;
	int j;

//	fp=fopen("restart.txt","r");
	fp=fopen(p->restartfname,"r");
	j=NG;
	while (!feof(fp)) {
//   		if (fscanf(fp, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &x, &ru, &iu, &rv, 
//   					&iv, &rd, &id, &vx, &vy, &d) != 10)   break;
//		if (fscanf(fp,"%lg",&d) != 1) break;
		if (fscanf(fp,"%lg %lg %lg %lg",&x,&d,&vx,&vy) != 4) break;
//		XCORD = x;
//		UVEL = ru + I*iu;
//		VVEL = rv + I*iv;
//		SIG = rd + I*id;
		VXBAR = vx;
		VYBAR = vy;
		DBAR = d;
		j++;
	}
	
	fclose(fp);
	return;
}


void calcTwd(parameters *p, double ti) {
	int i,j,m;
	FILE *fp;
	char fname[50];
	double *Twd, *dxFb, *dxFp, *Th;
	double om = p->omega; 
	double nu = p->nu;
	double q = p->q;
	
	Twd = (double *)malloc(sizeof(double)*(p->Nx));
	dxFp = (double *)malloc(sizeof(double)*(p->Nx));
	dxFb = (double *)malloc(sizeof(double)*(p->Nx));
	Th = (double *)malloc(sizeof(double)*(p->Nx));

	for(i=0;i<p->Nx;i++) {Twd[i] = 0; dxFp[i]=0; dxFb[i]=0; Th[i]=0;}
	
	calc_stress(p);
	calc_derivs(p,1);

#ifdef OPENMP 
	#pragma omp parallel private(j,m) shared(p) num_threads(NUMTHREADS)
	#pragma omp for schedule(static)
#endif
	for(j=NG;j<p->Nx+NG;j++) {
		for(m=1;m<NK;m++) {
			
			Twd[j-NG] += 2*creal(-conj(UVEL(m))*DXVVEL(m)*DBAR
					 + (DXVYBAR+(2*om-q))*conj(SIG(m))*UVEL(m)
					 + conj(SIG(m))*SIG(m)*DXPIP(0,1,0)/(DBAR*DBAR)
					 - (conj(SIG(m))/DBAR)*(DXPIP(0,1,m)+I*(p->k[m])*PIP(1,1,m)));
			
			dxFp[j-NG] += 2*creal(DXMF(m)*conj(VVEL(m))+MF(m)*conj(DXVVEL(m)));
			dxFb[j-NG] += 2*creal((VYBAR+(2*om-q)*XCORD)*(conj(SIG(m))*DXUVEL(m)+conj(DXSIG(m))*UVEL(m))
						+conj(SIG(m))*UVEL(m)*(DXVYBAR+2*om-q)
						-DXPIP(0,1,m));	
			Th[j-NG] += 2*creal(conj(SIG(m))*I*(p->k[m])*PHI(m));	
		
		}
		dxFb[j-NG] += creal(DXMF(0)*(VYBAR+(2*om-q)*XCORD)+MF(0)*(DXVYBAR+2*om-q)-DXPIP(0,1,0));
	}

	sprintf(fname,"k0/amf_%g.txt",ti);
	fp=fopen(fname,"w");


	for(j=NG;j<p->Nx+NG;j++) {
		fprintf(fp,"%12.8e \t %12.8e \t %12.8e \t %12.8e \t %12.8e \n",
				XCORD,Twd[j-NG],dxFp[j-NG],Th[j-NG],dxFb[j-NG]);
	}

	fclose(fp);
	
	free(Twd); free(dxFp); free(dxFb); free(Th);

	return;
}
// void read_params(parameters *P) {
// 	FILE *fp;
// 	fp=fopen("input_params.txt","r");
// 	
// 	
// 	fscanf(fp,"Nx = %d \n", P->Nx);
// 	fscanf(fp,"Mp = %g \n", P->Mp);
// 	fscanf(fp,"nu = %g \n", P->nu);
// 	fscanf(fp,"Lx = %g \n", P->Lx);
// 	fscanf(fp,"dx = %g \n", P->dx);
// 	fscanf(fp,"xsoft = %g \n", P->xsoft);
// 	fscanf(fp,"k = %g \n", P->k);
// 	fscanf(fp,"q = %g \n", P->q);
// 	fscanf(fp,"c = %g \n", P->c);
// 	fscanf(fp,"omega = %g \n", P->omega);
// 	fscanf(fp,"Ly = %g \n", P->Ly);
// 	fscanf(fp,"start time = %g \n", P->t0);
// 	fscanf(fp,"end time = %g \n", P->endt);
// 	fscanf(fp,"number of files = %d \n", P->numf);
// 	
// 	
// 	fclose(fp);
// 	
// 	return;
// 
// 
// }
void allocate_params(parameters *p) {
	int i,j;
	
	p->k = (double *)malloc(NK*sizeof(double));
	p->x = (double *)malloc((p->Ntot)*sizeof(double));
  	p->phi = (double complex *)malloc((p->Nx)*NK*sizeof(double complex));
  	p->dxphi = (double complex *)malloc((p->Nx)*NK*sizeof(double complex));

  	p->cy = (double complex *)malloc((p->Ntot)*3*NK*sizeof(double complex));
  	p->dxcy = (double complex *)malloc((p->Ntot)*3*NK*sizeof(double complex));
  	p->t = (double *)malloc((p->numf)*sizeof(double));
  	
	p->rhs = (double complex *)malloc((p->Ntot)*(6*NK-3)*sizeof(double complex));
	p->realrhs = (double *)malloc((p->Ntot)*(6*NK-3)*sizeof(double));

	p->VT = (stress *)malloc(sizeof(stress));
	p->MF = (massflux *)malloc(sizeof(massflux));
	p->MF->mf = (double complex *)malloc((p->Ntot)*NK*sizeof(double complex));
	p->MF->dxmf = (double complex *)malloc((p->Ntot)*NK*sizeof(double complex));

	for(i=0;i<2;i++) {
		for(j=i;j<2;j++) {
			p->VT->p[i][j]=(double complex *)malloc((p->Ntot)*NK*sizeof(double complex));
			p->VT->pp[i][j]=(double complex *)malloc((p->Ntot)*NK*sizeof(double complex));
			p->VT->dxp[i][j]=(double complex *)malloc((p->Ntot)*NK*sizeof(double complex));
			p->VT->dxpp[i][j]=(double complex *)malloc((p->Ntot)*NK*sizeof(double complex));
		}
	}
	return;
}
void free_params(parameters *p) {
	int i,j;
	free(p->MF->mf);
	free(p->MF->dxmf);	
	
	for(i=0;i<2;i++) {	
		for(j=i;j<2;j++) {
			free(p->VT->p[i][j]);	
			free(p->VT->pp[i][j]);	
			free(p->VT->dxp[i][j]);	
			free(p->VT->dxpp[i][j]);	
		}
	}
	free(p->VT);
	free(p->MF);
	free(p->cy); free(p->x); free(p->rhs);
	free(p->dxcy); free(p->k);
	free(p->phi); free(p->dxphi);
	free(p->realrhs); free(p->t);
	
	free(p);
	
	return;

}

void calc_stress( parameters *p ) 
{
	
	int i,j,n,m;
	
	double c  = p->c;
	double nu = p->nu;
	double k;
	double q = p->q;
	double om = p->omega;
	double complex divvk =0;
	double divvb = 0;
	calc_derivs(p,0);

#ifdef OPENMP
	#pragma omp parallel private(j,i,n,m,k,divvb,divvk) shared(p) num_threads(NUMTHREADS)
	#pragma omp for schedule(static)
#endif
	for(j=NG;j<(p->Nx+NG);j++) {
		divvb = DXVXBAR*(2.0/3.0);
		for(m=0;m<NK;m++) {
			k = p->k[m];
			divvk = (2.0/3.0)*(DXUVEL(m)+I*k*VVEL(m));
		if (m==0) {
			p->MF->mf[m+NK*j] = DBAR*VXBAR;
			
#ifdef SIMPLEVISC		 
			p->VT->p[0][0][m+NK*j] = -c*c*DBAR + nu*DXVXBAR;
			p->VT->p[0][1][m+NK*j] = nu*(DXVYBAR-q*om);
			p->VT->p[1][1][m+NK*j] = -c*c*DBAR-nu*DXVXBAR;
#else
			p->VT->p[0][0][m+NK*j] = -c*c*DBAR + nu*DBAR*(2*DXVXBAR-divvb);
			p->VT->p[0][1][m+NK*j] = nu*DBAR*(DXVYBAR-q*om);
			p->VT->p[1][1][m+NK*j] = -c*c*DBAR - nu*DBAR*divvb;
#endif
			
		
		}
		else {
			p->MF->mf[m+NK*j] = SIG(m)*VXBAR + DBAR*UVEL(m);
	
#ifdef SIMPLEVISC		 
		p->VT->p[0][0][m+NK*j] = -c*c*SIG(m)+nu*(DXUVEL(m)-I*k*VVEL(m));
		p->VT->p[0][1][m+NK*j] = nu*(DXVVEL(m)+I*k*UVEL(m));
		p->VT->p[1][1][m+NK*j] = -c*c*SIG(m)+nu*(I*k*VVEL(m)-DXUVEL(m));
#else
		p->VT->p[0][0][m+NK*j] = -c*c*SIG(m)+nu*(DBAR*(2*DXUVEL(m)-divvk)+SIG(m)*(2*DXVXBAR-divvb));
		p->VT->p[0][1][m+NK*j] = nu*DBAR*(DXVVEL(m)+I*k*UVEL(m))+nu*SIG(m)*(DXVYBAR-q*om);
		p->VT->p[1][1][m+NK*j] = -c*c*SIG(m)+nu*(DBAR*(2*I*k*VVEL(m) - divvk)+SIG(m)*(-divvb));
	
		p->VT->pp[0][0][m+NK*j] = nu*conj(SIG(m))*(2*DXUVEL(m)-divvk);
		p->VT->pp[0][1][m+NK*j] = nu*conj(SIG(m))*(DXVVEL(m)+I*k*UVEL(m));
		p->VT->pp[1][1][m+NK*j] = nu*conj(SIG(m))*(2*I*k*VVEL(m)-divvk);
		 
#endif		
		}
		}

	}
	
	return;
}
void calc_derivs(parameters *p, int flag) {

	int i,j,m,d,n;

#ifdef OPENMP
	#pragma omp parallel private(j,i,m,d,n) shared(p) num_threads(NUMTHREADS)
	#pragma omp for schedule(static)
#endif
	for(j=NG;j<p->Nx+NG;j++) {
	for(m=0;m<NK;m++) {
	if (flag==0) {      // Just do cy
		for(i=0;i<3;i++) {
			p->dxcy[CINDX(i,m)]=0;
			for(d=0;d<2*NG;d++) {
				 p->dxcy[CINDX(i,m)]+=(p->cy[i+3*m+3*NK*(j+deriv_ind[d])])*deriv_coeffs[d];
			}
			p->dxcy[CINDX(i,m)] /= (p->dx);
		}	
	}
	if (flag==1) {      // Just do stress and mass flux
		for(i=0;i<2;i++) {
			for(n=i;n<2;n++) {
				DXPIP(i,n,m) =  0;
				DXPP(i,n,m) = 0;
				for(d=0;d<2*NG;d++) {
					DXPIP(i,n,m) += 
								  (p->VT->p[i][n][m+NK*(j+deriv_ind[d])])*deriv_coeffs[d];  
					DXPP(i,n,m) += 
										(p->VT->pp[i][n][m+NK*(j+deriv_ind[d])])*deriv_coeffs[d];  
				
				}
				DXPIP(i,n,m) /= (p->dx);
				DXPP(i,n,m) /= (p->dx);
			}
		}	
		DXMF(m) = 0; 
		for(d=0;d<2*NG;d++) {
			DXMF(m) += (p->MF->mf[m+NK*(j+deriv_ind[d])])*deriv_coeffs[d];
		}
		DXMF(m) /= (p->dx); 
	}
	}
	}
	

	return;

}
void wavekillbc(parameters *p,double dt)
{
	int j,m;
	double R,tau;
	double x_inf = p->x[0] + (p->Lx)*0.05;
	double x_sup = p->x[p->Nx-1] - (p->Lx)*0.05;
	 
	for(j=NG;j<(p->Nx+NG);j++) { 
		
		if (XCORD > x_sup) R = (XCORD-x_sup)/(p->x[p->Nx-1] - x_sup);
		if (XCORD < x_inf) R = (x_inf - XCORD)/(x_inf - p->x[0]);
		R *= R;
 		tau = 2*M_PI/30;
 		if (R>0.0) {
 			tau /= R;	
 			p->cy[CINDX(0,0)] = (p->cy[CINDX(0,0)])/(1+dt);
			p->cy[CINDX(1,0)] = (p->cy[CINDX(1,0)])/(1+dt);
			p->cy[CINDX(2,0)] = ((p->cy[CINDX(2,0)])*tau + 1.0*dt)/(dt+tau); 
		}
	}
	for(j=0;j<NG;j++) { 
		
		if (XCORD > x_sup) R = (XCORD-x_sup)/(p->x[p->Nx-1] - x_sup);
		if (XCORD < x_inf) R = (x_inf - XCORD)/(x_inf - p->x[0]);
		R *= R;
 		tau = 2*M_PI/30;
 		if (R>0.0) {
 			tau /= R;	
 			p->cy[CINDX(0,0)] = (p->cy[CINDX(0,0)])/(1+dt);
			p->cy[CINDX(1,0)] = (p->cy[CINDX(1,0)])/(1+dt);
			p->cy[CINDX(2,0)] = ((p->cy[CINDX(2,0)])*tau + 1.0*dt)/(dt+tau); 
		}
	}
	
	return;
}

void bounds( parameters *p, double t) {
	int i,j,m;
											
	for(j=0;j<NG;j++) {
		 for(m=0;m<NK;m++) {
			 for (i=0;i<3;i++) {
			 
#ifdef OUTBC
			p->cy[CINDX(i,m)] = p->cy[i+3*m+3*NK*(2*NG-j)];
#endif	

			}
		}
	}

		
	for(j=p->Nx+NG;j<p->Ntot;j++) {
		for(m=0;m<NK;m++) {	
			for(i=0;i<3;i++) {
#ifdef OUTBC
				p->cy[CINDX(i,m)] = p->cy[i+3*m+3*NK*(2*(p->Nx+NG-1)-j)];
#endif
			}
		}
	}
	
  

}

void initphi(parameters *p) {
	int j,m;
	double re, im, dxre,dxim;
	FILE *fp,*fpx;
	char pystr[50];
	
	sprintf(pystr,"python fftpot.py %lg %lg %lg %d %d %d -1",p->xsoft,p->Lx,p->Ly,p->Nx,p->Ny,NK);
	system(pystr);
	
	fp = fopen("ftpot.dat","r");
	fpx = fopen("dxftpot.dat","r");
	
	for(j=NG;j<p->Nx+NG;j++) {
		for(m=0;m<NK;m++) {
			fscanf(fp,"%lg %lg",&re,&im);
			fscanf(fpx,"%lg %lg",&dxre,&dxim);
			PHI(m) = p->Mp*(re+I*im);
			DXPHI(m) = p->Mp*(dxre+I*dxim);
		}
	}
	fclose(fp); fclose(fpx);

#ifdef OUTPHI
	write_phi(p);
#endif	
	
	return;

}
double complex calc_pot(double complex phi,double t, double tau) {
	
	if (tau==0) return phi;
	else return (1-exp(-t/tau))*phi;
	
}
void write_phi(parameters *p) {
	int j,m;
	FILE *fp, *dxfp;
	
	fp = fopen("planet.dat","w");
	dxfp = fopen("dxplanet.dat","w");
	
	for(j=NG;j<p->Nx+NG;j++) {
		for(m=0;m<NK;m++) {
			fprintf(fp,"%lg \t %lg \t",creal(PHI(m)),cimag(PHI(m)));
			fprintf(dxfp,"%lg \t %lg \t",creal(DXPHI(m)),cimag(DXPHI(m)));
		}
		fprintf(fp,"\n");
		fprintf(dxfp,"\n");
	}

	fclose(fp); fclose(dxfp);
	
	return;
}

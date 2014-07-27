#include "meanwave.h"

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

double complex calc_pot(double complex phi,double t, double tau) {
	
	if (tau==0) return phi;
	else return (1-exp(-t/tau))*phi;
	
}

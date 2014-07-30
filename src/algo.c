#include "meanwave.h"


int func (double t, const double y[], double f[],void *params) {

	Field *fld = (Field *)params;
	
	fill_rhs(fld,t);

/* Copy complex data into real array with u,v,sig all combined for gsl */
	for(i=0;i<NTOTC;i++) {
		y[i] = creal(



}
void fill_rhs(Field *fld, double t) {
/* 		Fill the RHS of the EOM */
	int i;
	double qom = (fld->q)*(fld->omega);
	double om2 = 2*(fld->omega);
	double k;
/* Fill the derivative arrays */
	calc_deriv(fld->u,fld->dxu);
	calc_deriv(fld->v,fld->dxv);
	calc_deriv(fld->sig,fld->dxsig);

/*	Fill RHS arrays with any non-convolution terms*/	
	for(i=0;i<NTOTC;i++) {
		k = fld->k[i];
		fld->dtu[i] = qom*I*k*(fld->u[i])*(fld->x[i]); + om2*(fld->v[i]) 
						- calc_pot(fld->dxphi[i],t);
		fld->dtv[i] = qom*I*k*(fld->v[i])*(fld->x[i]) +(qom-om2)*(fld->u[i])
						-I*k*calc_pot(fld->phi[i],t);
		fld->dtsig[i] = qom*I*k*(fld->x[i])*(fld->sig[i]);
	}
	
/* Start adding the convolutions. */
	
	convolve(fld->u,fld->dxu,fld->dtu,-1);
	convolve(fld->v,fld->dyu,fld->dtu,-1);
	
	convolve(fld->u,fld->dxv,fld->dtv,-1);
	convolve(fld->v,fld->dyv,fld->dtv,-1);

	convolve(fld->sig,fld->dxu,fld->dtsig,-1);
	convolve(fld->dxsig,fld->u.fld->dtsig,-1);
	convolve(fld->sig,fld->dyv,fld->dtsig,-1);
	convolve(fld->dysig,fld->v,fld->dtsig,-1);
	
/* Add viscosity and pressure */
	
	add_viscosity(fld);
	

	
	return;	
}	

void add_advec(Field *fld) {
	int i;

	convolve(fld->u,fld->dxu,fld->dtu,-1);
	convolve(fld->v,fld->dyu,fld->dtu,-1);
	
	convolve(fld->u,fld->dxv,fld->dtv,-1);
	convolve(fld->v,fld->dyv,fld->dtv,-1);

	convolve(fld->sig,fld->dxu,fld->dtsig,-1);
	convolve(fld->dxsig,fld->u.fld->dtsig,-1);
	convolve(fld->sig,fld->dyv,fld->dtsig,-1);
	convolve(fld->dysig,fld->v,fld->dtsig,-1);
	
}

void add_mass(Field *fld) {
	int i;

}

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
double complex calc_pot(double complex phi,double t, double tau) {
	
	if (tau==0) return phi;
	else return (1-exp(-t/tau))*phi;
	
}

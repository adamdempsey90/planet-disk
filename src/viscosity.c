#include "meanwave.h"
void calc_stress( parameters *p ) {
	
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
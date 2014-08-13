#include "planetdisk.h"

double complex *Twd, *dxFb, *dxFp, *Th;

void amf(Field *fld) {
	int i,j;
	
	calc_deriv(fld->u,fld->dxu,fld->dyu,fld->Params->dx,fld->k);
	calc_deriv(fld->v,fld->dxv,fld->dyv,fld->Params->dx,fld->k);
	calc_deriv(fld->sig,fld->dxsig,fld->dysig,fld->Params->dx,fld->k);
	visc_tens(fld);

	calc_deriv(fld->Tens->Pixy,fld->Tens->divPiy,NULL,fld->Params->dx,fld->k);
	calc_deriv(fld->Tens->Piyy,NULL,fld->Tens->divPiy,fld->Params->dx,fld->k);

 	convolve_inv(&fld->sig[istart],&fld->Tens->divPiy[istart],Twd,1);



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

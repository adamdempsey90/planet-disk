#include "planetdisk.h"

Cons *cns; 
double *Twd, *dxFb, *dxFp, *Th;
double *Fb, *Fp, *Ft;
double complex *twdpi;

void amf(Field *fld) {
	int i,j,indx;
	double om2 = 2*(fld->Params->omega);
	
	
	for(i=0;i<Nx*NC;i++) twdpi[i] = 0;
	
	calc_deriv(fld->u,fld->dxu,fld->dyu,fld->Params->dx,fld->k);
	calc_deriv(fld->v,fld->dxv,fld->dyv,fld->Params->dx,fld->k);
	calc_deriv(fld->sig,fld->dxsig,fld->dysig,fld->Params->dx,fld->k);
	visc_tens(fld);

	calc_deriv(fld->Tens->Pixy,fld->Tens->divPiy,NULL,fld->Params->dx,fld->k);
	calc_deriv(fld->Tens->Piyy,NULL,fld->Tens->divPiy,fld->Params->dx,fld->k);

	convolve(&fld->u[istart],&fld->sig[istart],&cns->Mx[istart],1);
	convolve(&fld->v[istart],&fld->sig[istart],&cns->My[istart],1);

	convolve_inv(&fld->sig[istart],&fld->Tens->divPiy[istart],twdpi,1);


	for(i=0;i<Nx+2*NG;i++) {
		Twd=[i] = twdpi[i*NC];
		Fb[i] = (fld->vy[i*NC]+om2*(fld->x[i]))*(fld->sig[i*NC])*(fld->u[i*NC])
				- fld->Tens->Pixy[i*NC];
		for(j=0;j<NC;j++) {
			
			if (i<Nx) {
				indx = j+(i+NG)*NC;
				Twd[i] += -(fld->sig[istart+i*NC])*conj(fld->u[indx])*(fld->dxv[indx])
					+(fld->dxv[istart+i*NC]+om2*(fld->x[i+NG]))*conj(fld->sig[indx])*(fld->u[indx]);
			}
			Fb[i] += (fld->vy[i*NC]+om2*(fld->x[i]))*conj(fld->sig[CINDX])*(fld->u[CINDX]);

		}
	
	calc_deriv(Fb,dxFb,NULL,fld->Params->dx,NULL);
	

	for(i=0;i<Nx+2*NG;i++) {
		Fb
		for(j=1;j<NC;j++) {
			Fb[i] += 2*creal((fld->sig[CINDX])*conj(fld->u[CINDX])*(fld->v[NC*i]+om2*(fld->x[i]))
					- fld->Tens->Pixy[NC*i]);
			


	return;
}



void init_amf(void) {

	cns = (Cons *)malloc(sizeof(Cons));
	cns->Mx = (double complex *)malloc(sizeof(double complex)*NTOTC);
	cns->My = (double complex *)malloc(sizeof(double complex)*NTOTC);
	cns->rMx = (double *)(cns->Mx);
	cns->rMy = (double *)(cns->My); 

	Twd = (double *)malloc(sizeof(double)*Nx);
	dxFb = (double *)malloc(sizeof(double)*Nx);
	dxFp = (double *)malloc(sizeof(double)*Nx);
	Fb = (double *)malloc(sizeof(double)*(Nx+2*NG));
	Fp = (double *)malloc(sizeof(double)*(Nx+2*NG));
	Ft = (double *)malloc(sizeof(double)*(Nx+2*NG));
	Th = (double *)malloc(sizeof(double)*Nx);
	
	twdpi = (double complex *)malloc(sizeof(double complex)*Nx*NC);
	
	return;
}

void free_amf(void) {

	free(cns->Mx); free(cns->My);
	free(cns);
	free(Twd); free(dxFb); free(dxFp);
	free(Fb); free(Fp); free(Ft);
	free(Th);

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

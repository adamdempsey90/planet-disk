#include "planetdisk.h"

void init(Field *fld) {
	int i,j;
	double Lx = fld->Params->Lx;
	double Ly = fld->Params->Ly;
	
	for(i=0;i<Nx;i++) fld->x[i] = -.5*Lx + (i+.5)*(Lx/Nx);
	for(i=0;i<Ny;i++) fld->y[i] = -.5*Ly + (i+.5)*(Ly/Ny);
	for(i=0;i<NC;i++) fld->k[i]= i*(M_PI*2/Ly);
	
	initphi(fld);
	
	for(i=0;i<Nx;i++) {
		for(j=0;j<NC;j++) {
			fld->u[CINDX] = 0;
			fld->v[CINDX] = 0;
			if(j==0)	fld->sig[CINDX] = fld->Params->sig0;
			
			else	fld->sig[CINDX] = 0;
			fld->dxu[CINDX] = 0;
			fld->dxv[CINDX] = 0;
			fld->dxsig[CINDX] = 0;
			
			fld->Tens->Txx[CINDX] = 0;
			fld->Tens->Txy[CINDX] = 0;
			fld->Tens->Tyy[CINDX] = 0;
			
			fld->Tens->Pixx[CINDX] = 0;
			fld->Tens->Pixy[CINDX] = 0;
			fld->Tens->Piyy[CINDX] = 0;

			fld->Tens->divPix[CINDX] = 0;
			fld->Tens->divPiy[CINDX] = 0;

		}
	}
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





void allocate_field(Field *fld) {

	
	fld->Tens = (Stress *)malloc(sizeof(Stress));
		
	fld->x = (double *)malloc(sizeof(double)*Nx);
	fld->y = (double *)malloc(sizeof(double)*Ny);
	fld->k = (double *)malloc(sizeof(double)*NC);
	
	fld->u = (double complex *)malloc(sizeof(double complex)*NTOTC);
	fld->v = (double complex *)malloc(sizeof(double complex)*NTOTC);
	fld->sig = (double complex *)malloc(sizeof(double complex)*NTOTC);

	fld->dxu = (double complex *)malloc(sizeof(double complex)*NTOTC);
	fld->dxv = (double complex *)malloc(sizeof(double complex)*NTOTC);
	fld->dxsig = (double complex *)malloc(sizeof(double complex)*NTOTC);
	
	fld->dyu = (double complex *)malloc(sizeof(double complex)*NTOTC);
	fld->dyv = (double complex *)malloc(sizeof(double complex)*NTOTC);
	fld->dysig = (double complex *)malloc(sizeof(double complex)*NTOTC);		

	fld->dtu = (double complex *)malloc(sizeof(double complex)*NTOTC);
	fld->dtv = (double complex *)malloc(sizeof(double complex)*NTOTC);
	fld->dtsig = (double complex *)malloc(sizeof(double complex)*NTOTC);	
	
	fld->phi = (double complex *)malloc(sizeof(double complex)*NTOTC);
	fld->dxphi = (double complex *)malloc(sizeof(double complex)*NTOTC);
	
	fld->Tens->Txx = (double complex *)malloc(sizeof(double complex)*NTOTC);
	fld->Tens->Txy = (double complex *)malloc(sizeof(double complex)*NTOTC);
	fld->Tens->Tyy = (double complex *)malloc(sizeof(double complex)*NTOTC);

	fld->Tens->Pixx = (double complex *)malloc(sizeof(double complex)*NTOTC);
	fld->Tens->Pixy = (double complex *)malloc(sizeof(double complex)*NTOTC);
	fld->Tens->Piyy = (double complex *)malloc(sizeof(double complex)*NTOTC);

	fld->Tens->divPix = (double complex *)malloc(sizeof(double complex)*NTOTC);
	fld->Tens->divPiy = (double complex *)malloc(sizeof(double complex)*NTOTC);

	return;
}

void free_field(Field *fld) {

	free(fld->x);
	free(fld->y);
	free(fld->k);
	
	free(fld->u);
	free(fld->v);
	free(fld->sig);

	free(fld->dxu);
	free(fld->dxv);
	free(fld->dxsig);
	
	free(fld->dyu);
	free(fld->dyv);
	free(fld->dysig);		

	free(fld->dtu);
	free(fld->dtv);
	free(fld->dtsig);	
	
	free(fld->phi);
	free(fld->dxphi);
	
	free(fld->Tens->Txx);
	free(fld->Tens->Txy);
	free(fld->Tens->Tyy);

	free(fld->Tens->Pixx);
	free(fld->Tens->Pixy);
	free(fld->Tens->Piyy);

	free(fld->Tens->divPix);
	free(fld->Tens->divPiy);

	free(fld->Tens);
	free(fld->Params);
	free(fld);
	
	return;

}


void initphi(Field *fld) {
	int i,j;
	double rad, xs;
	double complex *cphi = (double complex *)malloc(sizeof(double complex)*NTOC);
	double *rphi = (double *)cphi;
	double complex *cdxphi = (double complex *)malloc(sizeof(double complex)*NTOC);
	double *rdxphi = (double *)cdxphi;
	
	for(i=0;i<Nx;i++) {
		for(j=0;j<NR;j++) {
			if (j<Ny) {
				xs = x[i]*x[i] + (fld->Params->xs)*(fld->Params->xs);
				rad = y[j]*y[j]+xs;
				rphi[j+i*NR] = -(fld->Params->Mp)/sqrt(rad);
				rdxphi[j+i*NR] = (fld->Params->Mp)*xs*pow(rad,-1.5);
			}
			else {
				rphi[j+i*NR]=0;
				rdxphi[j+i*NR] = 0;
			}
		}
	}
	fft_phi(rphi,cphi);
	fft_dxphi(rdxphi,cdxphi);
	
	
	memcpy(fld->phi,cphi,sizeof(double complex)*NTOTC);
	memcpy(fld->dxphi,cdxphi,sizeof(double complex)*NTOTC);

	free(cphi); free(cdxphi);
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
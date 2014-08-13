#include "planetdisk.h"

void init(Field *fld) {
	int i,j;
	double Lx = fld->Params->Lx;
	double Ly = fld->Params->Ly;
	
	printf("\t Initializing Coordinates...\n");
	for(i=0;i<(Nx+2*NG);i++) fld->x[i] = (i-NG+.5)*(Lx/Nx);
	for(i=0;i<Ny;i++) fld->y[i] = -.5*Ly + (i+.5)*(Ly/Ny);
	printf("\t Initializing Boundary Conditions...\n");
	for(i=0;i<NC;i++) {
		fld->k[i]= i*(M_PI*2/Ly);
	}
	printf("\t Initializing Gravitational Potential...\n");
	initphi(fld);

/* Set initial conditions here */
	for(i=NG;i<Nx+NG;i++) {
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
			
// 			fld->Tens->Pixx[CINDX] = -(fld->Params->c)(fld->Params->c)*(fld->Params->sig0);
// 			fld->Tens->Pixy[CINDX] = -(fld->Params->nu)*(fld->Params->q)*(fld->Params->omega)*(fld->Params->sig0);
// 			fld->Tens->Piyy[CINDX] = -(fld->Params->c)(fld->Params->c)*(fld->Params->sig0);

			fld->Tens->Pixx[CINDX] = 0;
			fld->Tens->Pixy[CINDX] = 0;
			fld->Tens->Piyy[CINDX] = 0;

			fld->Tens->divPix[CINDX] = 0;
			fld->Tens->divPiy[CINDX] = 0;

		}
	}
	for(i=0;i<istart;i++) {
		fld->dxu[i] = 0;
		fld->dxv[i] = 0;
		fld->dxsig[i] = 0;
		fld->phi[i] = 0;
		fld->dxphi[i] = 0;
		fld->Tens->divPix[i] = 0;
		fld->Tens->divPiy[i] = 0;
		fld->Tens->Txx[i] = 0;
		fld->Tens->Txy[i] = 0;
		fld->Tens->Tyy[i] = 0;
	}
	for(i=iend;i<NTOTC;i++) {
		fld->dxu[i] = 0;
		fld->dxv[i] = 0;
		fld->dxsig[i] = 0;
		fld->phi[i] = 0;
		fld->dxphi[i] = 0;
		fld->Tens->divPix[i] = 0;
		fld->Tens->divPiy[i] = 0;
		fld->Tens->Txx[i] = 0;
		fld->Tens->Txy[i] = 0;
		fld->Tens->Tyy[i] = 0;
	}
	return;
}

// void restart(parameters *p) {
// 	FILE *fp;
// 	double ru,iu,rv,iv,rd,id,vx,vy,d,x;
// 	int j;
// 
// //	fp=fopen("restart.txt","r");
// 	fp=fopen(p->restartfname,"r");
// 	j=NG;
// 	while (!feof(fp)) {
// //   		if (fscanf(fp, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &x, &ru, &iu, &rv, 
// //   					&iv, &rd, &id, &vx, &vy, &d) != 10)   break;
// //		if (fscanf(fp,"%lg",&d) != 1) break;
// 		if (fscanf(fp,"%lg %lg %lg %lg",&x,&d,&vx,&vy) != 4) break;
// //		XCORD = x;
// //		UVEL = ru + I*iu;
// //		VVEL = rv + I*iv;
// //		SIG = rd + I*id;
// 		VXBAR = vx;
// 		VYBAR = vy;
// 		DBAR = d;
// 		j++;
// 	}
// 	
// 	fclose(fp);
// 	return;
// }




void init_derivs(void) {

	int i, d;
	for(i=-NG, d=0;i<=NG;i++) {
		if (i!=0) {
			deriv.ind[d] = i;
			d++;
		}
	}
	
#ifdef SECDERIV	
	deriv.coeffs[0] = -0.5;
	deriv.coeffs[1] = 0.5; 
#endif

#ifdef FOURTHDERIV
	deriv.coeffs[0] = 1.0/12;
	deriv.coeffs[1] = -2.0/3;
	deriv.coeffs[2] = 2.0/3;
	deriv.coeffs[3] = -1.0/12;	
#endif

#ifdef SIXTHDERIV
	deriv.coeffs[0] = -1.0/60;
	deriv.coeffs[1] = 3.0/20;
	deriv.coeffs[2] = -0.75;
	deriv.coeffs[3] = 0.75;
	deriv.coeffs[4] = -3.0/20;
	deriv.coeffs[5] = 1.0/60;
#endif

#ifdef EIGTHDERIV
	deriv.coeffs[0] = 1.0/280;
	deriv.coeffs[1] = -4.0/105;
	deriv.coeffs[2] = 1.0/5;
	deriv.coeffs[3] = -4.0/5;
	deriv.coeffs[4] = 4.0/5;
	deriv.coeffs[5] = -1.0/5;
	deriv.coeffs[6] = 4.0/105;
	deriv.coeffs[7] = -1.0/280;
#endif

#ifdef TENTHDERIV
	deriv.coeffs[0] = -2.0/2520;
	deriv.coeffs[1] = 25.0/2520;
	deriv.coeffs[2] = -150.0/2520;
	deriv.coeffs[3] = 600.0/2520;
	deriv.coeffs[4] = -2100.0/2520;
	deriv.coeffs[5] = 2100.0/2520;
	deriv.coeffs[6] = -600.0/2520;
	deriv.coeffs[7] = 150.0/2520;
	deriv.coeffs[8] = -25.0/2520;
	deriv.coeffs[9] = 2.0/2520;
#endif

	return;
}

void allocate_field(Field *fld) {

	
	fld->Tens = (Stress *)malloc(sizeof(Stress));
		
	fld->x = (double *)malloc(sizeof(double)*(Nx+2*NG));
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

	fld->dtu = (double complex *)malloc(sizeof(double complex)*Nx*NC);
	fld->dtv = (double complex *)malloc(sizeof(double complex)*Nx*NC);
	fld->dtsig = (double complex *)malloc(sizeof(double complex)*Nx*NC);	

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

	fld->vx = (double *)malloc(sizeof(double)*Nx*NR);
	fld->vy =  (double *)malloc(sizeof(double)*Nx*NR);
	fld->dens =  (double *)malloc(sizeof(double)*Nx*NR);

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

	free(fld->vx); 
	free(fld->vy);
	free(fld->dens);
	
	free(fld->Tens);
	free(fld->Params);
	free(fld);
	
	return;

}


void initphi(Field *fld) {
	int i,j;
	double rad, xs;
	double *rphi = (double *)malloc(sizeof(double)*Nx*NR);
	double *rdxphi = (double *)malloc(sizeof(double)*Nx*NR);


	for(i=0;i<Nx;i++) {
		for(j=0;j<NR;j++) {
			if (j<Ny) {
				xs = (fld->x[i+NG])*(fld->x[i+NG]) + (fld->Params->xs)*(fld->Params->xs);
				rad = (fld->y[j])*(fld->y[j])+xs;
				rphi[RINDX] = -(fld->Params->Mp)/sqrt(rad);
				rdxphi[RINDX] = (fld->Params->Mp)*(fld->x[i])*pow(rad,-1.5);
			}
			else {
				rphi[RINDX]=0;
				rdxphi[RINDX] = 0;
			}
		}
	}
	
	fft_phi(rphi,fld->phi);
	fft_phi(rdxphi,fld->dxphi);


	free(rphi); free(rdxphi);
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



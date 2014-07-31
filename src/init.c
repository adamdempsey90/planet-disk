#include "planetdisk.h"

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
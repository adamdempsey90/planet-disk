#include "meanwave.h"
// void read_input(parameters *p) {
// 	int i;
// 	printf("Reading inputs...\n");
// 
// 
// 	p->Nx=2048; 
// 	p->Ny = 512;
// 	p->Ntot = p->Nx+2*NG;
// 	p->Lx=120.0; 
// 	p->Ly=30.0;
// 	p->dx=(p->Lx)/(p->Nx); 
// 	p->xsoft = .6;
// 	p->c=1; 
// 	p->Mp=.5; 
// 	p->nu=0.0024; 
// 	p->q=1.5; 
// 	p->omega=1.0;
//  
// 	p->kfac = 2.0;
// //  	p->k *= p->kfac;
// 
// 	p->t0 = 0;
// 	p->tau = 0;
// 	p->endt = 100;
// 	p->numf = 100;	
// 
// 	sprintf(p->restartfname, "newm0.5nu0.0024_lx60_means.dat");
// 	p->conv_flag = 0;
// 	
// 	p->conv_tol = 1e-4;
// 
// //	read_params(p);	
// 
// 	
// 	
// 	printf("Finished reading inputs...\n");
// 	return;
// }

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
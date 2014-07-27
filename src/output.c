#include "meanwave.h"
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

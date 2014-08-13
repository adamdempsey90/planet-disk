#include "planetdisk.h"

/* Option to write output in real space or complex space */

void output(Field *fld) {
	FILE *fu, *fv, *fs;
	char fnameu[50], fnamev[50], fnames[50];


	sprintf(fnameu,"outputs/vx_%d.dat",outnum);
	sprintf(fnamev,"outputs/vy_%d.dat",outnum);
	sprintf(fnames,"outputs/dens_%d.dat",outnum);


	fu = fopen(fnameu,"wb");
	if (fu == NULL) printf("ERROR: Couldn't open vx file\n");
	fv = fopen(fnamev,"wb");
	if (fv == NULL) printf("ERROR: Couldn't open vy file\n");
	fs = fopen(fnames,"wb");
	if (fs == NULL) printf("ERROR: Couldn't open dens file\n");



#ifdef OUTGHOST
	fwrite((double *)fld->u,sizeof(double),NTOTR,fu);
	fwrite((double *)fld->v,sizeof(double),NTOTR,fv);
	fwrite((double *)fld->sig,sizeof(double),NTOTR,fs);
#else
	fwrite((double *)&fld->u[istart],sizeof(double),Nx*NR,fu);
	fwrite((double *)&fld->v[istart],sizeof(double),Nx*NR,fv);
	fwrite((double *)&fld->sig[istart],sizeof(double),Nx*NR,fs);

#endif	

	fclose(fu); fclose(fv); fclose(fs);	
	
	outnum++;
	
	return;

}

void output_coords(Field *fld) {
	FILE *f;
	f = fopen("outputs/coords.dat","wb");
	fwrite(&fld->x[NG],sizeof(double),Nx,f);
	fwrite(fld->k,sizeof(double),NC,f);
	fwrite(fld->y,sizeof(double),Ny,f);
	fclose(f);
	return;
}

void output_derivs(Field *fld) {
	FILE *fu, *fv, *fs;
	char fnameu[50], fnamev[50], fnames[50];
	
	sprintf(fnameu,"outputs/dxvx_%d.dat",dxoutnum);
	sprintf(fnamev,"outputs/dxvy_%d.dat",dxoutnum);
	sprintf(fnames,"outputs/dxdens_%d.dat",dxoutnum);

	fu = fopen(fnameu,"wb");
	if (fu == NULL) printf("ERROR: Couldn't open dxvx file\n");
	fv = fopen(fnamev,"wb");
	if (fv == NULL) printf("ERROR: Couldn't open dxvy file\n");
	fs = fopen(fnames,"wb");
	if (fs == NULL) printf("ERROR: Couldn't open dxdens file\n");

#ifdef OUTGHOST	
	fwrite((double *)fld->dxu,sizeof(double),NTOTR,fu);
	fwrite((double *)fld->dxv,sizeof(double),NTOTR,fv);
	fwrite((double *)fld->dxsig,sizeof(double),NTOTR,fs);
#else
	fwrite((double *)&fld->dxu[istart],sizeof(double),Nx*NR,fu);
	fwrite((double *)&fld->dxv[istart],sizeof(double),Nx*NR,fv);
	fwrite((double *)&fld->dxsig[istart],sizeof(double),Nx*NR,fs);
#endif
	fclose(fu); fclose(fv); fclose(fs);	
	dxoutnum++;
	return;
}
void output_rhs(Field *fld) {
	
	FILE *fu, *fv, *fs;
	char fnameu[50], fnamev[50], fnames[50];


	sprintf(fnameu,"outputs/dtvx_%d.dat",dtoutnum);
	sprintf(fnamev,"outputs/dtvy_%d.dat",dtoutnum);
	sprintf(fnames,"outputs/dtdens_%d.dat",dtoutnum);

	fu = fopen(fnameu,"wb");
	if (fu == NULL) printf("ERROR: Couldn't open dtvx file\n");
	fv = fopen(fnamev,"wb");
	if (fv == NULL) printf("ERROR: Couldn't open dtvy file\n");
	fs = fopen(fnames,"wb");
	if (fs == NULL) printf("ERROR: Couldn't open dtdens file\n");

	fwrite((double *)fld->dtu,sizeof(double),Nx*NR,fu);
	fwrite((double *)fld->dtv,sizeof(double),Nx*NR,fv);
	fwrite((double *)fld->dtsig,sizeof(double),Nx*NR,fs);

	fclose(fu); fclose(fv); fclose(fs);	
	dtoutnum++;
	return;
}
void output_pi(Field *fld) {
	
	FILE *fxx, *fxy, *fyy, *fdx, *fdy;
	char fnamexx[50], fnameyy[50], fnamexy[50], fnamedx[50], fnamedy[50];


	sprintf(fnamexx,"outputs/pixx_%d.dat",pioutnum);
	sprintf(fnamexy,"outputs/pixy_%d.dat",pioutnum);
	sprintf(fnameyy,"outputs/piyy_%d.dat",pioutnum);
	sprintf(fnamedy,"outputs/divpix_%d.dat",pioutnum);
	sprintf(fnamedx,"outputs/divpiy_%d.dat",pioutnum);
	
	fxx = fopen(fnamexx,"wb");
	if (fxx == NULL) printf("ERROR: Couldn't open pixx file\n");
	fxy = fopen(fnamexy,"wb");
	if (fxy == NULL) printf("ERROR: Couldn't open pixy file\n");
	fyy = fopen(fnameyy,"wb");
	if (fyy == NULL) printf("ERROR: Couldn't open piyy file\n");
	fdx = fopen(fnamedx,"wb");
	if (fdx == NULL) printf("ERROR: Couldn't open divpix file\n");
	fdy = fopen(fnamedy,"wb");
	if (fdy == NULL) printf("ERROR: Couldn't open divpiy file\n");
	
#ifdef OUTGHOST	
	fwrite((double *)fld->Tens->Pixx,sizeof(double),NTOTR,fxx);
	fwrite((double *)fld->Tens->Pixy,sizeof(double),NTOTR,fxy);
	fwrite((double *)fld->Tens->Piyy,sizeof(double),NTOTR,fyy);
	fwrite((double *)fld->Tens->divPix,sizeof(double),NTOTR,fdx);
	fwrite((double *)fld->Tens->divPiy,sizeof(double),NTOTR,fdy);
#else
	fwrite((double *)&fld->Tens->Pixx[istart],sizeof(double),Nx*NR,fxx);
	fwrite((double *)&fld->Tens->Pixy[istart],sizeof(double),Nx*NR,fxy);
	fwrite((double *)&fld->Tens->Piyy[istart],sizeof(double),Nx*NR,fyy);
	fwrite((double *)&fld->Tens->divPix[istart],sizeof(double),Nx*NR,fdx);
	fwrite((double *)&fld->Tens->divPiy[istart],sizeof(double),Nx*NR,fdy);
#endif

	fclose(fxx); fclose(fxy); fclose(fyy); fclose(fdx); fclose(fdy);	
	pioutnum++;
	return;
}
void output_reals(Field *fld) {
	FILE *fu, *fv, *fs;
	char fnameu[50], fnamev[50], fnames[50];
	
	transform(fld);
	sprintf(fnameu,"outputs/rvx_%d.dat",dxoutnum);
	sprintf(fnamev,"outputs/rvy_%d.dat",dxoutnum);
	sprintf(fnames,"outputs/rdens_%d.dat",dxoutnum);

	fu = fopen(fnameu,"wb");
	if (fu == NULL) printf("ERROR: Couldn't open dxvx file\n");
	fv = fopen(fnamev,"wb");
	if (fv == NULL) printf("ERROR: Couldn't open dxvy file\n");
	fs = fopen(fnames,"wb");
	if (fs == NULL) printf("ERROR: Couldn't open dxdens file\n");


	fwrite(fld->vx,sizeof(double),Nx*NR,fu);
	fwrite(fld->vy,sizeof(double),Nx*NR,fv);
	fwrite(fld->dens,sizeof(double),Nx*NR,fs);

	fclose(fu); fclose(fv); fclose(fs);	
	dxoutnum++;
	return;
}
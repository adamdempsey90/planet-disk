#include "planetdisk.h"

#define STRLEN 100
/* Option to write output in real space or complex space */

void write_header(FILE *f);
void write_cheader(FILE *f);



void output(Field *fld) {
	FILE *fu, *fv, *fs;
	char fnameu[STRLEN], fnamev[STRLEN], fnames[STRLEN];
	char uname[STRLEN], vname[STRLEN], sname[STRLEN];
	
	strcpy(fnameu,fld->Params->outdir); 
	strcpy(fnamev,fld->Params->outdir); 
	strcpy(fnames,fld->Params->outdir);
	sprintf(uname,"vx_%d.dat",outnum);
	sprintf(vname,"vy_%d.dat",outnum);
	sprintf(sname,"dens_%d.dat",outnum);
	
	strcat(fnameu,uname); 
	strcat(fnamev,vname);
	strcat(fnames,sname);

// 	sprintf(fnameu,"outputs/id%d/vx_%d.dat",rank,outnum);
// 	sprintf(fnamev,"outputs/id%d/vy_%d.dat",rank,outnum);
// 	sprintf(fnames,"outputs/id%d/dens_%d.dat",rank,outnum);


	fu = fopen(fnameu,"wb");
	if (fu == NULL) printf("ERROR: Couldn't open vx file\n");
	fv = fopen(fnamev,"wb");
	if (fv == NULL) printf("ERROR: Couldn't open vy file\n");
	fs = fopen(fnames,"wb");
	if (fs == NULL) printf("ERROR: Couldn't open dens file\n");


	write_cheader(fu); write_cheader(fv); write_cheader(fs);
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
	char fname[STRLEN];
	
	strcpy(fname,fld->Params->outdir);
	strcat(fname,"coords.dat");
	
//	sprintf(fname,"outputs/id%d/coords.dat",rank);

	f = fopen(fname,"wb");
	write_header(f);
#ifdef OUTGHOST
	fwrite(fld->x,sizeof(double),Nx+2*NG,f);

#else
	
	fwrite(&fld->x[NG],sizeof(double),Nx,f);
#endif
	fwrite(fld->k,sizeof(double),NC,f);
	fwrite(fld->y,sizeof(double),Ny,f);
	fclose(f);
	return;
}

void output_derivs(Field *fld) {
	FILE *fu, *fv, *fs;
	char fnameu[STRLEN], fnamev[STRLEN], fnames[STRLEN];
	char uname[STRLEN], vname[STRLEN], sname[STRLEN];
	
	strcpy(fnameu,fld->Params->outdir); 
	strcpy(fnamev,fld->Params->outdir); 
	strcpy(fnames,fld->Params->outdir);
	sprintf(uname,"dxvx_%d.dat",dxoutnum);
	sprintf(vname,"dxvy_%d.dat",dxoutnum);
	sprintf(sname,"dxdens_%d.dat",dxoutnum);
	
	strcat(fnameu,uname); 
	strcat(fnamev,vname);
	strcat(fnames,sname);
	
// 	sprintf(fnameu,"outputs/id%d/dxvx_%d.dat",rank,dxoutnum);
// 	sprintf(fnamev,"outputs/id%d/dxvy_%d.dat",rank,dxoutnum);
// 	sprintf(fnames,"outputs/id%d/dxdens_%d.dat",rank,dxoutnum);

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
	char fnameu[STRLEN], fnamev[STRLEN], fnames[STRLEN];
	char uname[STRLEN], vname[STRLEN], sname[STRLEN];
	
	strcpy(fnameu,fld->Params->outdir); 
	strcpy(fnamev,fld->Params->outdir); 
	strcpy(fnames,fld->Params->outdir);
	sprintf(uname,"dtvx_%d.dat",dtoutnum);
	sprintf(vname,"dtvy_%d.dat",dtoutnum);
	sprintf(sname,"dtdens_%d.dat",dtoutnum);
	
	strcat(fnameu,uname); 
	strcat(fnamev,vname);
	strcat(fnames,sname);
	

// 	sprintf(fnameu,"outputs/id%d/dtvx_%d.dat",rank,dtoutnum);
// 	sprintf(fnamev,"outputs/id%d/dtvy_%d.dat",rank,dtoutnum);
// 	sprintf(fnames,"outputs/id%d/dtdens_%d.dat",rank,dtoutnum);

	fu = fopen(fnameu,"wb");
	if (fu == NULL) printf("ERROR: Couldn't open dtvx file\n");
	fv = fopen(fnamev,"wb");
	if (fv == NULL) printf("ERROR: Couldn't open dtvy file\n");
	fs = fopen(fnames,"wb");
	if (fs == NULL) printf("ERROR: Couldn't open dtdens file\n");

	write_cheader(fu); write_cheader(fv); write_cheader(fs);
	fwrite((double *)fld->dtu,sizeof(double),Nx*NR,fu);
	fwrite((double *)fld->dtv,sizeof(double),Nx*NR,fv);
	fwrite((double *)fld->dtsig,sizeof(double),Nx*NR,fs);

	fclose(fu); fclose(fv); fclose(fs);	
	dtoutnum++;
	return;
}
void output_pi(Field *fld) {
	
	FILE *fxx, *fxy, *fyy;
	char fnamexx[STRLEN], fnameyy[STRLEN], fnamexy[STRLEN];
	char xxname[STRLEN], yyname[STRLEN], xyname[STRLEN];


	strcpy(fnamexx,fld->Params->outdir); 
	strcpy(fnameyy,fld->Params->outdir); 
	strcpy(fnamexy,fld->Params->outdir);

	sprintf(xxname,"pixx_%d.dat",pioutnum);
	sprintf(xyname,"pixy_%d.dat",pioutnum);
	sprintf(yyname,"piyy_%d.dat",pioutnum);
	
	strcat(fnamexx,xxname); 
	strcat(fnamexy,xyname);
	strcat(fnameyy,yyname);
	
	

// 	sprintf(fnamexx,"outputs/id%d/pixx_%d.dat",rank,pioutnum);
// 	sprintf(fnamexy,"outputs/id%d/pixy_%d.dat",rank,pioutnum);
// 	sprintf(fnameyy,"outputs/id%d/piyy_%d.dat",rank,pioutnum);

	
	fxx = fopen(fnamexx,"wb");
	if (fxx == NULL) printf("ERROR: Couldn't open pixx file\n");
	fxy = fopen(fnamexy,"wb");
	if (fxy == NULL) printf("ERROR: Couldn't open pixy file\n");
	fyy = fopen(fnameyy,"wb");
	if (fyy == NULL) printf("ERROR: Couldn't open piyy file\n");

	
#ifdef OUTGHOST	
	fwrite((double *)fld->Tens->Pixx,sizeof(double),NTOTR,fxx);
	fwrite((double *)fld->Tens->Pixy,sizeof(double),NTOTR,fxy);
	fwrite((double *)fld->Tens->Piyy,sizeof(double),NTOTR,fyy);
#else
	fwrite((double *)&fld->Tens->Pixx[istart],sizeof(double),Nx*NR,fxx);
	fwrite((double *)&fld->Tens->Pixy[istart],sizeof(double),Nx*NR,fxy);
	fwrite((double *)&fld->Tens->Piyy[istart],sizeof(double),Nx*NR,fyy);
#endif

	fclose(fxx); fclose(fxy); fclose(fyy); 	
	pioutnum++;
	return;
}
void output_reals(Field *fld) {
	FILE *fu, *fv, *fs;
	char fnameu[STRLEN], fnamev[STRLEN], fnames[STRLEN];
	char vxname[STRLEN], vyname[STRLEN], densname[STRLEN];
	
	strcpy(fnameu,fld->Params->outdir); 
	strcpy(fnamev,fld->Params->outdir); 
	strcpy(fnames,fld->Params->outdir);
	sprintf(vxname,"rvx_%d.dat",dxoutnum);
	sprintf(vyname,"rvy_%d.dat",dxoutnum);
	sprintf(densname,"rdens_%d.dat",dxoutnum);
	
	strcat(fnameu,vxname); 
	strcat(fnamev,vyname);
	strcat(fnames,densname);

// 	sprintf(fnameu,"outputs/id%d/rvx_%d.dat",rank,dxoutnum);
// 	sprintf(fnamev,"outputs/id%d/rvy_%d.dat",rank,dxoutnum);
// 	sprintf(fnames,"outputs/id%d/rdens_%d.dat",rank,dxoutnum);

	fu = fopen(fnameu,"wb");
	if (fu == NULL) printf("ERROR: Couldn't open dxvx file\n");
	fv = fopen(fnamev,"wb");
	if (fv == NULL) printf("ERROR: Couldn't open dxvy file\n");
	fs = fopen(fnames,"wb");
	if (fs == NULL) printf("ERROR: Couldn't open dxdens file\n");
	
	transform(fld);
	write_header(fu); write_header(fv); write_header(fs);
	fwrite(fld->vx,sizeof(double),Nx*NR,fu);
	fwrite(fld->vy,sizeof(double),Nx*NR,fv);
	fwrite(fld->dens,sizeof(double),Nx*NR,fs);

	fclose(fu); fclose(fv); fclose(fs);	
	dxoutnum++;
	return;
}
void output_phi(Field *fld) {
	FILE *fp, *fdp;
	char fnamep[STRLEN], fnamedp[STRLEN];
	
	strcpy(fnamep,fld->Params->outdir); strcpy(fnamedp,fld->Params->outdir);
	strcat(fnamep,"phi.dat"); strcat(fnamedp,"dxphi.dat");
	
// 	sprintf(fnamep,"outputs/id%d/phi.dat",rank);
// 	sprintf(fnamedp,"outputs/id%d/dxphi.dat",rank);

	fp = fopen(fnamep,"wb");
	if (fp == NULL) printf("ERROR: Couldn't open phi file\n");
	fdp = fopen(fnamedp,"wb");
	if (fdp == NULL) printf("ERROR: Couldn't open dxphi file\n");



	write_cheader(fp); write_cheader(fdp);
	fwrite((double *)&fld->phi[istart],sizeof(double),Nx*NR,fp);
	fwrite((double *)&fld->dxphi[istart],sizeof(double),Nx*NR,fdp);

	fclose(fp); fclose(fdp); 	
	return;
}
void output_defines(char *dir) {
	char fname[STRLEN];
	strcpy(fname,dir);
	strcat(fname,"params.txt");
	FILE *f = fopen(fname,"w");
	fprintf(f,"Configuration parameters:\n");
#ifdef RESTART
	printf("\t Restarting from file\n");
	fprintf(f,"\tRestarting from file\n");
#endif
#ifdef WAVEEVOLVE
	printf("\t Evolving wave componenets\n");
	fprintf(f,"\tEvolving wave componenets\n");
#endif
#ifdef BACKEVOLVE
	printf("\t Evolving mean componenets\n");
	fprintf(f,"\tEvolving mean componenets\n");
#endif
	printf("\t Finite difference order is %d \n",2*NG);
	fprintf(f,"\tFinite difference order is %d \n",2*NG);
#ifdef SIMPLEVISC
	printf("\t Using the incompressible viscous stress tensor\n");
	fprintf(f,"\tUsing the incompressible viscous stress tensor\n");
#endif
#ifdef WAVEKILLBC
	printf("\t Enforcing wave killing boundary conditions\n");
	fprintf(f,"\tEnforcing wave killing boundary conditions\n");
#endif
#ifdef OPENMP
	printf("\t Using OpenMp with %d threads\n", NUMTHREADS);
	fprintf(f,"\tUsing OpenMp with %d threads\n", NUMTHREADS);
#endif

#ifdef OUTGHOST
	printf("\t Outputting the ghost zones\n");
	fprintf(f,"\tOutputting the ghost zones\n");
#endif
	fclose(f);
	return;
}
void output_params(Field *fld) {
	int i,Nxtot;
	char fname[STRLEN];
	
	for(i=0, Nxtot=0;i<np;i++) Nxtot+=Nxproc[i];
	
	strcpy(fname,fld->Params->outdir);
	strcat(fname,"params.txt");
	FILE *f = fopen(fname,"a");
	fprintf(f,"Input Parameters: \n \
		Nx = %d\n \
		Ny = %d\n \
		Lx = %lg\n \
		Ly =  %lg\n \
		dx = %lg\n \
		xsoft =  %lg\n \
		h =  %lg\n \
		Mp =  %lg\n \
		nu =  %lg\n \
		q =  %lg\n \
		omega =  %lg\n \
		sig0 = %lg\n \
		Time Parameters \n \
		t0 =  %lg\n \
		tau =  %lg\n \
		endt =  %lg\n \
		numf =  %d\n \
		tol =  %lg\n",
			  Nxtot, fld->Params->Ny, fld->Params->Lx, fld->Params->Ly, fld->Params->dx, fld->Params->xs, fld->Params->h, 
			  fld->Params->Mp, fld->Params->nu, fld->Params->q, fld->Params->omega, fld->Params->sig0,
			  fld->Params->t0, fld->Params->tau,   fld->Params->endt, fld->Params->numf,fld->Params->tol);
	fclose(f);


	return;

}
void write_header(FILE *f) {

#ifdef OUTGHOST
	double dNx = (double)(Nx+2*NG);
#else
	double dNx = (double)Nx;
#endif
	double dNy = (double)Ny;
	fwrite(&dNx,sizeof(double),1,f);
	fwrite(&dNy,sizeof(double),1,f);

	return;
}
void write_cheader(FILE *f) {

#ifdef OUTGHOST
	double complex dNx = (double complex)(Nx+2*NG);
#else
	double complex dNx = (double complex)Nx;
#endif
	double complex dNy = (double complex)Ny;
	fwrite((double *)&dNx,sizeof(double complex),1,f);
	fwrite((double *)&dNy,sizeof(double complex),1,f);

	return;
}
void init_output(char *dir) {
	char idstr[50];	
	size_t len = strlen(dir);
	if (dir[len-1] != '/') dir[len] = '/';
	mkdir(dir,0777);
	sprintf(idstr,"id%d/",rank);
	strcat(dir,idstr);
	mkdir(dir,0777);

	return;
}
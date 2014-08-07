#include "planetdisk.h"

/* Option to write output in real space or complex space */

void output(Field *fld, double t) {
	FILE *fu, *fv, *fs;
	char fnameu[50], fnamev[50], fnames[50];
	
	sprintf(fnameu,"outputs/vx_%d.dat",(int)t);
	sprintf(fnameu,"outputs/vy_%d.dat",(int)t);
	sprintf(fnameu,"outputs/dens_%d.dat",(int)t);


	fu = fopen(fnameu,"wb");
	fv = fopen(fnamev,"wb");
	fs = fopen(fnames,"wb");
	
	
#ifdef REALOUT

	fftw_execute_plan(fld->u,c2r1);
	fftw_execute_plan(fld->v,c2r1);
	fftw_execute_plan(fld->sig,c2r1);


	fwrite(fld->vx,sizeof(double),NTOTR,fu);
	fwrite(fld->vy,sizeof(double),NTOTR,fv);
	fwrite(fld->dens,sizeof(double),NTOTR,fs);
#else
	fwrite((double *)fld->y,sizeof(double),NTOTR,fu);
	fwrite((double *)fld->v,sizeof(double),NTOTR,fv);
	fwrite((double *)fld->sig,sizeof(double),NTOTR,fs);

#endif	

	fclose(fu); fclose(fv); fclose(fs);	
	

	
	return;

}

void output_coords(Field *fld) {
	int i;
	FILE *f;
	f = fopen("outputs/coords.dat","w");
	for(i=0;i<fld->Params->Nx;i++) fprintf(f,"%lg \t",fld->x[i]);
	fprintf(f,"\n");
	for(i=0;i<fld->Params->Nk;i++) fprintf(f,"%lg \t",fld->k[i]);


#ifdef REALOUT
	fprintf("\n");
	for(i=0;i<fld->Params->Ny;i++) fprintf(f,"%lg \t",fld->y[i]);
#endif

	return;
}


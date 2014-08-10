#include "planetdisk.h" 

void read_input(Field *fld) {
	FILE *f;
	char garbage[100];
	
	f=fopen("params.in","r");
	fgets(garbage,sizeof(garbage),f);
	fscanf(f,"Nx = %d \n",&fld->Params->Nx);
	fscanf(f,"Ny = %d \n",&fld->Params->Ny); 
	fscanf(f,"Lx = %lg \n",&fld->Params->Lx);
	fscanf(f,"Ly =  %lg \n",&fld->Params->Ly);
	fscanf(f,"xsoft =  %lg \n",&fld->Params->xs);
	fscanf(f,"c =  %lg \n",&fld->Params->c);
	fscanf(f,"Mp =  %lg \n",&fld->Params->Mp);
	fscanf(f,"nu =  %lg \n",&fld->Params->nu);
	fscanf(f,"q =  %lg \n",&fld->Params->q);
	fscanf(f,"omega =  %lg \n",&fld->Params->omega);
	fscanf(f,"sig0 =  %lg \n",&fld->Params->sig0);
	fgets(garbage,sizeof(garbage),f);
	fscanf(f,"t0 =  %lg \n",&fld->Params->t0);
	fscanf(f,"tau =  %lg \n",&fld->Params->tau);
	fscanf(f,"endt =  %lg \n",&fld->Params->endt);
	fscanf(f,"numf =  %d \n",&fld->Params->numf);

	fclose(f);

  	printf("\tInput Parameters \n \
  	  Nx = %d\n \
  	  Ny = %d\n \
  	  Lx = %lg\n \
  	  Ly =  %lg\n \
  	  xsoft =  %lg\n \
  	  c =  %lg\n \
  	  Mp =  %lg\n \
  	  nu =  %lg\n \
  	  q =  %lg\n \
  	  omega =  %lg\n \
  	  sig0 = %lg\n \
  	  Time Parameters \n \
  	  t0 =  %lg\n \
  	  tau =  %lg\n \
  	  endt =  %lg\n \
  	  numf =  %d\n", 
  	  fld->Params->Nx, fld->Params->Ny, fld->Params->Lx, fld->Params->Ly, fld->Params->xs, fld->Params->c, 
  	  fld->Params->Mp, fld->Params->nu, fld->Params->q, fld->Params->omega, fld->Params->sig0,
  	  fld->Params->t0, fld->Params->tau,   fld->Params->endt, fld->Params->numf);
  	
  	Nx = fld->Params->Nx;
  	Ny = fld->Params->Ny;
  	fld->Params->dx = (fld->Params->Lx/Nx);
  	
	return;


}

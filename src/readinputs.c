#include "planetdisk.h" 

void read_input(char *inputdir, Field *fld) {
	FILE *f;
	char garbage[100], outdir[100], ifname[100];
	double h,Lx,Ly,xs,Mp,nu,q,omega,sig0,t0,tau,endt,tol;
	int numf;
	
	size_t len = strlen(inputdir);
	if (inputdir[len-1] != '/') inputdir[len] = '/';
	strcpy(ifname,inputdir);
	strcat(ifname,"params.in");
	f=fopen(ifname,"r");
	if (f==NULL) printf("\n\nERROR Can't Find Input File!\n\n");
	fgets(garbage,sizeof(garbage),f);
	fscanf(f,"Nx = %d \n",&Nx);
	fscanf(f,"Ny = %d \n",&Ny); 
	fscanf(f,"Lx = %lg \n",&Lx);
	fscanf(f,"Ly =  %lg \n",&Ly);
	fscanf(f,"xsoft =  %lg \n",&xs);
	fscanf(f,"h =  %lg \n",&h);
	fscanf(f,"Mp =  %lg \n",&Mp);
	fscanf(f,"nu =  %lg \n",&nu);
	fscanf(f,"q =  %lg \n",&q);
	fscanf(f,"omega =  %lg \n",&omega);
	fscanf(f,"sig0 =  %lg \n",&sig0);
	fgets(garbage,sizeof(garbage),f);
	fscanf(f,"t0 =  %lg \n",&t0);
	fscanf(f,"tau =  %lg \n",&tau);
	fscanf(f,"endt =  %lg \n",&endt);
	fscanf(f,"numf =  %d \n",&numf);
	fscanf(f,"tol =  %lg \n",&tol);
	fscanf(f,"outputdir = %s \n",outdir);
		
	fclose(f);
	
	fld->Params->Nx = Nx;
	fld->Params->Ny= Ny;
	fld->Params->Lx = Lx*h; 
	fld->Params->Ly = Ly*h; 
	fld->Params->xs = xs*h;
	fld->Params->h = h;
	fld->Params->c = h*omega;
	fld->Params->Mp = Mp;
	fld->Params->nu = nu*h*h*omega;
	fld->Params->q = q;
	fld->Params->omega = omega;
	fld->Params->sig0 = sig0;
	fld->Params->t0 = t0;
	fld->Params->tau = tau;
	fld->Params->endt = endt;
	fld->Params->numf = numf;	
	fld->Params->tol = tol;		
  	fld->Params->dx = (fld->Params->Lx/Nx);
  	strcpy(fld->Params->outdir,outdir);
   
  	mkdir(outdir,0777);
	
  	MPI_Printf("Input Parameters \n \
  	 \t Nx = %d\n \
  	 \t Ny = %d\n \
  	 \t Lx = %lg\n \
  	 \t Ly =  %lg\n \
  	 \t dx = %lg\n \
  	 \t xsoft =  %lg\n \
  	 \t h =  %lg\n \
  	 \t Mp =  %lg\n \
  	 \t nu =  %lg\n \
  	 \t q =  %lg\n \
  	 \t omega =  %lg\n \
  	 \t sig0 = %lg\n \
  	 \t Time Parameters \n \
  	 \t t0 =  %lg\n \
  	 \t tau =  %lg\n \
  	 \t endt =  %lg\n \
  	 \t numf =  %d\n \
  	 \t tol =  %lg\n \
  	 \t outputdir = %s\n", 
  	  fld->Params->Nx, fld->Params->Ny, fld->Params->Lx, fld->Params->Ly, fld->Params->dx, fld->Params->xs, fld->Params->h, 
  	  fld->Params->Mp, fld->Params->nu, fld->Params->q, fld->Params->omega, fld->Params->sig0,
  	  fld->Params->t0, fld->Params->tau,   fld->Params->endt, fld->Params->numf,fld->Params->tol,fld->Params->outdir);
  
//   	if (rank==0) {
// 		f = fopen("outputs/params.txt","a");
// 		fprintf(f,"Input Parameters: \n \
// 		Nx = %d\n \
// 		Ny = %d\n \
// 		Lx = %lg\n \
// 		Ly =  %lg\n \
// 		dx = %lg\n \
// 		xsoft =  %lg\n \
// 		h =  %lg\n \
// 		Mp =  %lg\n \
// 		nu =  %lg\n \
// 		q =  %lg\n \
// 		omega =  %lg\n \
// 		sig0 = %lg\n \
// 		Time Parameters \n \
// 		t0 =  %lg\n \
// 		tau =  %lg\n \
// 		endt =  %lg\n \
// 		numf =  %d\n \
// 		tol =  %lg\n \
// 		outputdir = %s\n", 
// 			  fld->Params->Nx, fld->Params->Ny, fld->Params->Lx, fld->Params->Ly, fld->Params->dx, fld->Params->xs, fld->Params->h, 
// 			  fld->Params->Mp, fld->Params->nu, fld->Params->q, fld->Params->omega, fld->Params->sig0,
// 			  fld->Params->t0, fld->Params->tau,   fld->Params->endt, fld->Params->numf,fld->Params->tol,fld->Params->outdir);
// 		fclose(f);
// 	}
// 	
	
/* Send out data to rest of processors */
	
	int i;
	if (rank==0) {	
		for(i=0;i<np;i++) Nxproc[i] = Nx/np;
		Nxproc[np-1] += Nx % np;
	}
	

	return;


}

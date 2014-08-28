#include "planetdisk.h"
#include <time.h>

int main (int argc, char *argv[]) {
	int i;
	clock_t tic, toc;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&np);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if (rank==0)  tic= clock();
	Nxproc = (int *)malloc(sizeof(int)*np);
	MPI_Printf("Welcome to the planet disk code...\n");


	Field *fld = (Field *)malloc(sizeof(Field));
	fld->Params = (Parameters *)malloc(sizeof(Parameters));

	func_calls=0;
	
	
	init_derivs();	
	
	MPI_Printf("Reading Inputs...\n");
	
	if (rank==0) read_input(fld);

	MPI_Printf("Broadcasting Parameters...\n");
	MPI_Bcast(fld->Params,sizeof(Parameters),MPI_BYTE,0,MPI_COMM_WORLD);
	MPI_Bcast(Nxproc,np,MPI_INT,0,MPI_COMM_WORLD);
	
	fld->Params->Nx=Nxproc[rank];
	Nx = fld->Params->Nx;
	Ny = fld->Params->Ny;
	
	init_output(fld->Params->outdir);
	
	MPI_Printf("Code compiled with...\n");
	if (rank==0) {
		output_defines(fld->Params->outdir);
		output_params(fld);
	}
	printf("%d %s\n", rank, fld->Params->outdir);
  	  
	MPI_Printf("Initializing data structures...\n");

	allocate_field(fld);

	MPI_Printf("Initializing FFTW...\n");

	init_fft();
	init_rk45();
	init_buff();
	MPI_Printf("Initializing Field...\n");

	init(fld);

#ifdef OUTPHI
	MPI_Printf("Outputting Potential...\n");
	output_phi(fld);
#endif
	set_bc(fld);
	MPI_Printf("Outputting Coordinates...\n");

	
	output_coords(fld);
	outnum=0; dxoutnum=0; dtoutnum=0; pioutnum=0;
	
	MPI_Printf("Outputting Initial Conditions...\n");

	output(fld);



//	output_reals(fld);
 

  	double	h = .1;
  	double 	t=fld->Params->t0;
 	double 	t1 = fld->Params->endt;
  	double dt;
  	i=1;
  	func_calls = 0;
  
  	MPI_Printf("Max mode #%d at k=%lg\n",Nmax, kmax);
  	MPI_Printf("Effective y resolution of dy = %lg\n\n",2*M_PI/kmax);
	MPI_Printf("Starting the Time Loop...\n");
	MPI_Printf("\t Starting Time = %lg \t Ending Time = %lg \n",t,t1);
	
	numstep=0;
  while (t < t1)
    {
      
    dt = t;

    int status = rk45_step_apply(fld,&t,&h); 
    numstep++;
     if (status == -1) {
     	MPI_Printf("ERROR With Step...\nTerminating Run...\n");
        break;
    }
    dt = t-dt;
    
#ifdef SHEARSPLIT
	shear_advection(fld,dt);
#endif
	MPI_Printf ("\t step #%d, step size = %.5e, at t=%.5e \n", numstep,dt, t);
   

#ifdef WAVEKILLBC
	if (rank==0 || rank==np-1)	wavekillbc(fld,dt);
#endif
     
      
   	if( t >= fld->Params->t0 + i * (fld->Params->endt) / ((double) fld->Params->numf)) { 
  		 MPI_Printf ("\t\t OUTPUT %d, step size = %.5e, at t=%.5e \n", outnum,h,t);
		
		output(fld);
      	i++;
     }
    }


  
 	free(Nxproc);
  	free_field(fld);
  	free_buff();
  	fft_free(); 
  	free_rk45();
  
	if (rank==0) toc = clock();
	MPI_Printf("Elapsed time: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);	
	MPI_Printf("# steps per second: %f\n", numstep /((double)(toc - tic) / CLOCKS_PER_SEC));
  	int mpi_status = MPI_Finalize();
 	return mpi_status;
}

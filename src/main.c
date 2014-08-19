#include "planetdisk.h"

int main (void) {
	int i;

	printf("Welcome to the planet disk code...\n");
	printf("Code compiled with...\n");
	output_defines();

	Field *fld = (Field *)malloc(sizeof(Field));
	fld->Params = (Parameters *)malloc(sizeof(Parameters));

	func_calls=0;
	init_derivs();	
	
	printf("Reading Inputs...\n");
	
	read_input(fld);

	printf("Initializing data structures...\n");

	allocate_field(fld);

	printf("Initializing FFTW...\n");

	init_fft();
	init_rk45();
	printf("Initializing Field...\n");

	init(fld);
	set_bc(fld);
	printf("Outputting Coordinates...\n");

	output_coords(fld);
	outnum=0; dxoutnum=0; dtoutnum=0; pioutnum=0;
	
	printf("Outputting Initial Conditions...\n");

	output(fld);

	printf("Defining the ODE System...\n");


//	output_reals(fld);



 

  	double	h = .1;
  	double 	t=fld->Params->t0;
 	double 	t1 = fld->Params->endt;
  	double dt;
  	i=1;
  	func_calls = 0;
  
  	printf("Max mode #%d at k=%lg\n",Nmax, kmax);
  	printf("Effective y resolution of dy = %lg\n\n",2*M_PI/kmax);
	printf("Starting the Time Loop...\n");
	printf("\t Starting Time = %lg \t Ending Time = %lg \n",t,t1);
	
	
  while (t < t1)
    {
      
    dt = t;

    int status = rk45_step_apply(fld,&t,&h); 
    
     if (status == -1) {
     	printf("ERROR With Step...\nTerminating Run...\n");
        break;
    }
    dt = t-dt;
    
#ifdef SHEARSPLIT
	shear_advection(fld,dt);
#endif
	printf ("\t step size = %.5e, at t=%.5e \n", dt, t);
   

#ifdef WAVEKILLBC
	wavekillbc(fld,dt);
#endif
     
      
   	if( t >= fld->Params->t0 + i * (fld->Params->endt) / ((double) fld->Params->numf)) { 
   		 printf ("\t\t OUTPUT %d, step size = %.5e, at t=%.5e \n", outnum,h,t);
		
		output(fld);
      	i++;
     }
    }

	
  free_field(fld);
  fft_free(); 
  free_rk45();
  return 0;
}

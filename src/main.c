#include "planetdisk.h"

int main (void) {
	int i;
 	double *y;
	int gsl_size;
	printf("Welcome to the planet disk code...\n");
	printf("Reading Inputs...\n");
	printf("Initializing data structures...\n");
	Field *fld = (Field *)malloc(sizeof(Field));
	fld->Params = (Parameters *)malloc(sizeof(Parameters));

	func_calls=0;
	init_derivs();	
	
	printf("Reading Inputs...\n");
	
	read_input(fld);
	allocate_field(fld);

	printf("Initializing FFTW...\n");

	init_fft();

	printf("Initializing Field...\n");

	init(fld);
	set_bc(fld);
	printf("Outputting Coordinates...\n");

	output_coords(fld);
	outnum=0; dxoutnum=0; dtoutnum=0; pioutnum=0;
	
	printf("Outputting Initial Conditions...\n");

	output(fld);

	printf("Defining the ODE System...\n");

	gsl_size = 3*Nx*NR;

//	output_reals(fld);

	y = (double *)malloc(sizeof(double)*gsl_size);
	global_c2r(y,fld);

	gsl_odeiv2_system sys = {func, NULL, gsl_size, fld};


	
	const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;

  	gsl_odeiv2_step * s  = gsl_odeiv2_step_alloc (T, gsl_size);
 	gsl_odeiv2_control * c  = gsl_odeiv2_control_y_new (1e-8, 1e-8);
 	gsl_odeiv2_evolve * e  = gsl_odeiv2_evolve_alloc (gsl_size);
 
 

  	double	h = .1;
  	double 	t=fld->Params->t0;
 	double 	t1 = fld->Params->endt;
  	double dt;
  	i=1;
  	func_calls = 0;
  
	printf("Starting the Time Loop...\n");
	printf("\t Starting Time = %lg \t Ending Time = %lg \n",t,t1);
	
	
  while (t < t1)
    {
      
    dt = t;

    int status = gsl_odeiv2_evolve_apply (e, c, s,
                                           &sys, 
                                           &t, t1,
                                           &h, y);
    dt = t-dt;

	printf ("step size = %.5e, at t=%.5e \n", dt, t);
    if (status != GSL_SUCCESS) {
     	printf("ERROR With Step...\nTerminating Run...\n");
        break;
    }

     
      
   	if( t >= fld->Params->t0 + i * (fld->Params->endt) / ((double) fld->Params->numf)) { 
   		 printf ("\t OUTPUT %d, step size = %.5e, at t=%.5e \n", outnum,h,t);

		global_r2c(y,fld);
		
		output(fld);
      	i++;
     }
    }

  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
	
	
  free(y); 
  free_field(fld);
  fft_free(); 
  return 0;
}

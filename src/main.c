#include "planetdisk.h"

int main (void) {
	int i;
	double time=0;
 	double *y;

	printf("Welcome to the planet disk code...\n");
	printf("Reading Inputs...\n");
	printf("Initializing data structures...\n");
	Field *fld = (Field *)malloc(sizeof(Field));
	fld->Params = (Parameters *)malloc(sizeof(Parameters));

	init_derivs();	
	read_input(fld);
	allocate_field(fld);
	printf("Initializing FFTW...\n");
	init_fft();
	printf("Initializing Field...\n");
	init(fld);
	printf("Outputting Coordinates...\n");
//	output_coords(fld);
	
	printf("Defining the ODE System...\n");
	y = (double *)malloc(sizeof(double)*NTOTR);

	gsl_odeiv2_system sys = {func, NULL, NTOTR, fld};


	
	const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;

  	gsl_odeiv2_step * s  = gsl_odeiv2_step_alloc (T, NTOTR);
 	gsl_odeiv2_control * c  = gsl_odeiv2_control_y_new (1e-8, 1e-8);
 	gsl_odeiv2_evolve * e  = gsl_odeiv2_evolve_alloc (NTOTR);
 
 

  double	h = .1;
  double 	t=fld->Params->t0;
  double 	t1 = fld->Params->endt;
  double dt;
  i=1;
  
  
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
   		 printf ("\t OUTPUT step size = %.5e, at t=%.5e \n", h,t);

		global_r2c(y,fld);
		
//		output(fld,floor(t));
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

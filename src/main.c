#include "planetdisk.h"
int main (void) {
	int i;
	double time=0;
 	double *y;
  
	(Field *fld) = (Field *)malloc(sizeof(Field));
	fld->Params = (Parameters *)malloc(sizeof(Parameters));
	
	read_input(fld);
	allocate_field(fld);
	init(fld);

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
  
  
//	output_derivs(P);

  while (t < t1)
    {
      
      dt = t;

      int status = gsl_odeiv2_evolve_apply (e, c, s,
                                           &sys, 
                                           &t, t1,
                                           &h, y);
       dt = t-dt;

#ifdef WAVEKILLBC
       r_2_c(y,P->cy,P->Nx);
	   wavekillbc(P,dt);
	   c_2_r(y,P->cy,P->Nx);
#endif

	   printf ("step size = %.5e, at t=%.5e \n", dt, t);
      if (status != GSL_SUCCESS)
          break;

     
      
   	if( t >= fld->Params->t0 + i * (fld->Params->endt) / ((double) fld->Params->numf)) { 
   		 printf ("\t OUTPUT step size = %.5e, at t=%.5e \n", h,t);

		global_r2c(y,fld);
		
		calcTwd(P, floor(t));
		output_func(P,floor(t));
      	i++;
     }
    }

  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
	
#endif	
	
  free(y); 
  free_field(fld);
  free_fft(); 
  return 0;
}

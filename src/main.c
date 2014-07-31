#include "meanwave.h"
int main (void) {
	int i, cutind;
	double time=0;
 	double *y;
  
	(Field *fld) = (Field *)malloc(sizeof(Field));
	fld->Params = (Parameters *)malloc(sizeof(Parameters));
	
	read_input(fld);
	allocate_field(fld);






  gsl_odeiv2_system sys = {func, NULL, (6*NK-3)*(P->Ntot), P};



	
#ifdef LOWLEVEL
	
	const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;

  	gsl_odeiv2_step * s  = gsl_odeiv2_step_alloc (T, (6*NK-3)*(P->Ntot));
 	 gsl_odeiv2_control * c  = gsl_odeiv2_control_y_new (1e-8, 1e-8);
 	 gsl_odeiv2_evolve * e  = gsl_odeiv2_evolve_alloc ((6*NK-3)*(P->Ntot));
 
 

  double	h = .1;
  double 	t=P->t0;
  double 	t1 = P->endt;
  double dt;
  i=1;
  
  
//	output_derivs(P);

  while (t < t1)
    {
      
      dt = t;
      P->lastt = t;
      int status = gsl_odeiv2_evolve_apply (e, c, s,
                                           &sys, 
                                           &t, t1,
                                           &h, y);
       dt = t-dt;

#ifdef BACKEVOLVE
#ifdef WAVEKILLBC
       r_2_c(y,P->cy,P->Nx);
	   wavekillbc(P,dt);
	   c_2_r(y,P->cy,P->Nx);
#endif
#endif
#ifdef INTFAC
       r_2_c(y,P->cy,P->Nx);
	   remove_shear(P,dt);
	   c_2_r(y,P->cy,P->Nx);
#endif
	  
//       calc_split(P);	  
	   printf ("step size = %.5e, suggested step size = %.5e, at t=%.5e \n", dt, h,t);
	   printf("\t\t\t deriv count = %d\n", func_count); 
    	printf("\t\t\t CUT at X = %lg \n", P->xsplit);
//	   write_timestep(t,dt);
      if (status != GSL_SUCCESS)
          break;

     
      
   	if( t >= P->t0 + i * (P->endt) / ((double) P->numf)) { 
   		 printf ("\t OUTPUT step size = %.5e, at t=%.5e \n", h,t);
		r_2_c(y,P->cy,P->Nx);

//		c_2_r(y,P->cy,P->Nx);
		calcTwd(P, floor(t));
		output_func(P,floor(t));
//		output_amf(P,floor(t));
      	i++;
     }
    }

  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
	
#endif	
	
  free(y); 
  free_params(P); 
  return 0;
}

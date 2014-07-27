#include "meanwave.h"
int main (void)
{
	int i, cutind;
	double time=0;
 
  double *y;
  
  func_count=0;
  parameters *P = (parameters *)malloc(sizeof(parameters));
  
  read_input(P);
  printf("Allocating structs...\n");
  allocate_params(P);


  

  y = (double *)malloc((P->Ntot)*(6*NK-3)*sizeof(double));
  
  
  printf("Initialzing...\n");
  init(P,y);
#ifdef RESTART
	restart(P);
#endif
//  read_athena(P);
	printf("Appying BC's...\n");
  bounds(P,0);
	 

  c_2_r(y,P->cy,P->Nx);	
  printf("\tOUTPUT 0\n");
  output_func(P,0);
  
  gsl_odeiv2_system sys = {func, NULL, (6*NK-3)*(P->Ntot), P};


 	
//	gsl_odeiv2_driver_set_hmin(d,0.01);
  
	
	
	printf("Nx=%d, c=%g, Mp=%g, nu=%g, Nk=%d, Lx=%g dx=%g, xsoft=%g, t_end=%g \n",
			P->Nx, P->c, P->Mp, P->nu, NK, P->Lx, 
			P->dx, P->xsoft,P->endt);
	output_params(P);
	

  
	
	printf("Starting loop\n");
//	Evolve with high-level wrapper

#ifdef HIGHLEVEL
	 gsl_odeiv2_driver * d = 
 			 gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
				  1e-8, 1e-4, 0.0);
   for (i = 1; i <= P->numf; i++)
    {
      if (P->conf_flag == 1) {
      	printf("CONVERGENCE DETECTED \n");
      	break;
      }
      P->t[i-1] = P->t0 + i * (P->endt) / ((double) P->numf);
      P->lastt = time;
      int status = gsl_odeiv2_driver_apply (d, &time, P->t[i-1], y);
		printf("finished step %d, at t=%g\n", i, P->t[i-1]);
		printf("\t\t\t deriv count = %d\n", func_count);
		r_2_c(y,P->cy,P->Nx);
#ifdef INTFAC
		remove_shear(P,time-P->lastt);
#endif
		c_2_r(y,P->cy,P->Nx);
		bounds(P,time);
		
//		output_func(y,P,P->t[i-1]);
//		output_amf(P,P->t[i-1]);
      if (status != GSL_SUCCESS)
		{
		  printf ("error, return value=%d\n", status);
		  break;
		}
    }
    gsl_odeiv2_driver_free (d);
#endif

// Evolve with low-level wrapper to see step sizes.
	
	
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

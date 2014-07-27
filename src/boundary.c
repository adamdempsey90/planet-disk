#include "meanwave.h"
void wavekillbc(parameters *p,double dt)
{
	int j,m;
	double R,tau;
	double x_inf = p->x[0] + (p->Lx)*0.05;
	double x_sup = p->x[p->Nx-1] - (p->Lx)*0.05;
	 
	for(j=NG;j<(p->Nx+NG);j++) { 
		
		if (XCORD > x_sup) R = (XCORD-x_sup)/(p->x[p->Nx-1] - x_sup);
		if (XCORD < x_inf) R = (x_inf - XCORD)/(x_inf - p->x[0]);
		R *= R;
 		tau = 2*M_PI/30;
 		if (R>0.0) {
 			tau /= R;	
 			p->cy[CINDX(0,0)] = (p->cy[CINDX(0,0)])/(1+dt);
			p->cy[CINDX(1,0)] = (p->cy[CINDX(1,0)])/(1+dt);
			p->cy[CINDX(2,0)] = ((p->cy[CINDX(2,0)])*tau + 1.0*dt)/(dt+tau); 
		}
	}
	for(j=0;j<NG;j++) { 
		
		if (XCORD > x_sup) R = (XCORD-x_sup)/(p->x[p->Nx-1] - x_sup);
		if (XCORD < x_inf) R = (x_inf - XCORD)/(x_inf - p->x[0]);
		R *= R;
 		tau = 2*M_PI/30;
 		if (R>0.0) {
 			tau /= R;	
 			p->cy[CINDX(0,0)] = (p->cy[CINDX(0,0)])/(1+dt);
			p->cy[CINDX(1,0)] = (p->cy[CINDX(1,0)])/(1+dt);
			p->cy[CINDX(2,0)] = ((p->cy[CINDX(2,0)])*tau + 1.0*dt)/(dt+tau); 
		}
	}
	
	return;
}

void bounds( parameters *p, double t) {
	int i,j,m;
											
	for(j=0;j<NG;j++) {
		 for(m=0;m<NK;m++) {
			 for (i=0;i<3;i++) {
			 
#ifdef OUTBC
			p->cy[CINDX(i,m)] = p->cy[i+3*m+3*NK*(2*NG-j)];
#endif	

			}
		}
	}

		
	for(j=p->Nx+NG;j<p->Ntot;j++) {
		for(m=0;m<NK;m++) {	
			for(i=0;i<3;i++) {
#ifdef OUTBC
				p->cy[CINDX(i,m)] = p->cy[i+3*m+3*NK*(2*(p->Nx+NG-1)-j)];
#endif
			}
		}
	}
	
  

}

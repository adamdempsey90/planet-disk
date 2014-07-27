#include "meanwave.h" 
void read_input(parameters *p) { 
	 p->Nx = 2048;
	 NK = 25;
	 p->Ny = 512;
	 p->Lx = 120.0;
	 p->Ly = 30.0;
	 p->xsoft = .6;
	 p->c = 1;
	 p->Mp = .5;
	 p->nu = 0.0024;
	 p->q = 1.5;
	 p->omega = 1.0;
	 p->t0 = 0;
	 p->tau = 0;
	 p->endt = 100;
	 p->numf = 100;
	 sprintf(p->restartfname,"newm0.5nu0.0024_lx60_means.dat");
	 p->dx=(p->Lx)/(p->Nx);
	 p->Ntot = p->Nx+2*NG;
	 return; 
}

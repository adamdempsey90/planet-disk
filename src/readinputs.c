#include "planetdisk.h" 

void read_input(Field *fld) { 
	fld->Params->Nx = 1024;
	fld->Params->Ny = 16;
	Nx = fld->Params->Nx;
	Ny = fld->Params->Ny;
	fld->Params->Lx = 120.0;
	fld->Params->Ly = 30.0;
	fld->Params->xs = .6;
	fld->Params->c = 1;
	fld->Params->Mp = .5;
	fld->Params->nu = 0.0024;
	fld->Params->q = 1.5;
	fld->Params->omega = 1.0;
	fld->Params->sig0 = 1.0;
	fld->Params->t0 = 0;
	fld->Params->tau = 0;
	fld->Params->endt = 1e-8;
	fld->Params->numf = 1;
	sprintf(fld->Params->restartfname,"newm0.5nu0.0024_lx60_means.dat");
	fld->Params->dx=(fld->Params->Lx)/(fld->Params->Nx);
	return; 
}

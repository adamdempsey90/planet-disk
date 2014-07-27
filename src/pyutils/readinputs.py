#!/usr/bin/python

with open("../params.in","r") as f:
	lines = f.readlines()
	
lines = [x.split('\n')[0] for x in lines if x.split() != [] and '#' not in x]	

with open("../readinputs.c","w") as f:
	f.write('#include "meanwave.h" \n')
	f.write('void read_input(parameters *p) { \n')
	for line in lines:
		if 'restartfile' not in line:
			if 'NK' not in line:
				f.write('\t p->'+ line.split()[0]+' = '+line.split()[-1] + ';\n')
			else:
				f.write('\t '+ line.split()[0]+' = '+line.split()[-1] + ';\n')
		else:
			f.write('\t sprintf(p->restartfname,"'+line.split('=')[-1].split()[0]+'");\n')
	
	f.write('\t p->dx=(p->Lx)/(p->Nx);\n')
	f.write('\t p->Ntot = p->Nx+2*NG;\n')
	f.write('\t return; \n')
	f.write('}\n')
	

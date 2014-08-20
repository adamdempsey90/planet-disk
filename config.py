#!/usr/bin/python

from subprocess import call

with open("inputs/params.opt","r") as f:
	lines = f.readlines()

defines = [x.split('+')[-1].split('\n')[0] for x in lines if x.split() != [] and '#' not in x]	

print 'Planet Disk Code Configured With:'
for x in defines:
	print '\t'+x
	
if 'RKF' in defines:
	call(["cp", "src/integrator/rk45.c", "src/rk45.c"])
	call(["cp","src/integrator/rk45.h", "src/rk45.h"])
	call(["cp", "src/integrator/rkf.c", "src/rk45_step.c"])
	defines.remove('RKF')
elif 'RKCK' in defines:
	call(["cp", "src/integrator/rk45.c", "src/rk45.c"])
	call(["cp","src/integrator/rk45.h", "src/rk45.h"])
	call(["cp", "src/integrator/rkck.c", "src/rk45_step.c"])
	defines.remove('RKCK')

elif 'RKDP' in defines:
	call(["cp", "src/integrator/rk45.c", "src/rk45.c"])
	call(["cp","src/integrator/rk45.h", "src/rk45.h"])
	call(["cp", "src/integrator/rkdp.c", "src/rk45_step.c"])
	defines.remove('RKDP')

else:
	print 'No time integrator found!'
	
with open("src/defines.h","w") as g:
	for x in defines:
		g.write('#define ' + x + '\n')
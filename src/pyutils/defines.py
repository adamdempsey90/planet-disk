#!/usr/bin/python

with open("../params.opt","r") as f:
	lines = f.readlines()

lines = [x.split('+')[-1] for x in lines if x.split() != [] and '#' not in x]	
with open("../defines.h","w") as g:
	for x in lines:
		g.write('#define ' + x)
		
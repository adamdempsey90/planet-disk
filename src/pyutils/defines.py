#!/usr/bin/python

with open("inputs/params.opt","r") as f:
	lines = f.readlines()

lines = [x.split('+')[-1] for x in lines if x.split() != [] and '#' not in x]	
with open("src/defines.h","w") as g:
	for x in lines:
		g.write('#define ' + x)
		
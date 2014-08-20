planet-disk
===========

planet-disk interaction code


Set configuration options in inputs/params.opt
To compile the code run

./configure && make

To compile the executable planetdisk.

Set runtime prarameters in inputs/params.in and the execute the code to begin.

Outputs are in the outputs/ directory and are in binary C-major order for each variable. The python script pyutils/read.py has functions to read these output files. 

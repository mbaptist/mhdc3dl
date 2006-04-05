#!/usr/bin/python

#import mhdc3dl_python module
from mhdc3dl_python import *

#import the defaults from the file "defaults.cfg"
#into the __builtin__ module 
import imp
default="default.cfg"
imp.load_source('__builtin__',default)

#Main run cycle
#change as needed
for seed in range(100,200):
	for visc in range(0.5,0.01,-0.01):
		for diff in range(0.5,0.01,-0.01):
			for tcond in range(0.5,0.01,-0.01):
				lss_run("__builtin__")



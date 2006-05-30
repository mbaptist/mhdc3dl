#!/usr/bin/python

###############################################
# Run script and Input parameters for mhdc3dl #
###############################################

#This is part of mhdc3dl code
#Manuel Baptista 2004,2005,2006

### Description ###

#The input is done via this python file
#which is parsed by the c++ input class.
#Only variables declared and parsed in the 
#input class are passed to the code.

#This file can be loaded by c++, by running the
#mhdc3dl program, or invoked as a script from the
#command line or python interpreter. 

#This file is divided in two big sections:
#  default values - where all the input parameters are defined
#                    (all the variables parsed in the c++ input 
#                       class must e defined here)
#  launch script - this part runs only when this file is invoked 
#                    as a script; it calls the functions exposed
#                    in the module mhdc3dl_python.


print __name__


##################
# Default values #
##################

from math import pi,sqrt


### Run's Parameters ###

#Run's Name
runsname="default"



#Numerical grid
n1=32
n2=32
n3=16

#Physical sizes of system
l1=2*pi
l2=2*pi
l3=pi

#Physical extent for plan forms
ell2=4.
ell1=ell2*2.
#l1=2*pi/ell1
#l2=2*pi/ell2


#Physical parameters
visc=1.
omegaz=0.
compress=-1.
g=-1.
diff=1.
deltat=-1.
econd=0.
tcond=1.



### Basic Fields' Parameters ###

#Initialisation mode (load,random,plan,expression)
basic_mode='random'

#load mode parameters
basic_vel_fname=runsname+"_basic_vel.dat"
basic_mag_fname=runsname+"_basic_mag.dat"
basic_temp_fname=runsname+"_basic_temp.dat"

#random mode parameters
br_spectrum="power"
br_seed=102
br_ki=0
br_kf=5
br_alpha=4.  #ramdom power spectrum mode
br_rms_norm=1
br_kind=0
br_sym=1

#plan mode parameters


#expression mode parameters



### Large Scale Linear Stability ###
#linear solver parameters
ls_eps=1e-16
qq=.75
qq_adj=1.5
kk=5 
small=.01
small_adj=.005
#auxiliary problems
apbfname="aux"
#file for appending growth rate
grfname=""

### Short Scale Linear Stability ###
#vzdeigen parameters
ep=1e-8
thr=0.1
mp=2
sc=3000.
nseq=50
#initial fields
sss_ifname=''
sss_int_ofbname='sss_intermediate.dat'
sym_sub=1 #make both in code
sss_seed=666

### Time evolution ###



#################
# Launch script #
#################

# Invokes the functions exposed by the mhdc3dl_python interface
# For example, eddy=lss_run("__main__") runs the large scale
# stability code. You may redefine parameters as you please
# and make cycle for various sets of parameters


if (__name__=="__main__"):
  #import mhdc3dl_python module
  from mhdc3dl_python import *

  ofile=open("output.dat","w")

  visc=1.
  diff=1.
  tcond=1.

  #Main run cycle

#  for s in range(100,200):
#    br_seed=s
  eddy=lss_run("__main__") #this line runs the lss code for the parameters defined so far
  output_line=str(br_seed)+" "+str(visc)+" "+str(diff)+" "+str(tcond)+" "+str(eddy[0])+" "+str(eddy[1])+"\n"
  print(output_line)
#    ofile.write(line)
#    ofile.flush()



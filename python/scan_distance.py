#!/usr/bin/python
"""
A python script to launch PATRIC jobs 
"""

from math import *
from odict import * 
from socket import gethostname
import os
import shutil
import sys  
import re  
import time    

# constants
e_kin = 11.582
mp = 1.6605e-27
qe = 1.6022e-19
clight = 2.988e8
circum = 216.72 

#machines=['lxir054','lxir055','lxir056','lxir058']
machines=['lxir061','lxir060','lxir058','lxir059']    #49, 53, 54, 55, 56, 58, 59, 60, 61, 66 
machines=machines*8   

# Flexibity
bumpI = 1
#Bumpf = 160e-6  # korretur 
#BumpA = 85e-3
#Chopf = 65e-6
#ChopV = 15e-6# 55e-6-50e-6  # korretur  

#gamma0 = 1.0 + (e_kin*1e6*qe)/(mp*clight*clight) 
#beta0  = sqrt((gamma0*gamma0-1.0)/(gamma0*gamma0)) 
#t_turn = circum/(beta0*clight)
mtj=80#int(floor(Chopf/t_turn))    

# paret in Flexibity
#Amp0 = 77.8e-3
#delAmp = 4.57e-3
#BumpA = Amp0+ChopV/t_turn*delAmp
#Bumpf = BumpA*t_turn/delAmp      

#Amp0=BumpA-BumpA/Bumpf*ChopV      
#print Amp0


#emitanz
ex=10
ey=16

Qvec=3.30    
Q_hor=4.17    
Bumpf = 360e-6  # korretur 
BumpA = 60e-3
Chopf = 65e-6
ChopV = 0e-6# 55e-6-50e-6  # korretur
Amp0=BumpA-BumpA/Bumpf*ChopV     
gamma0 = 1.0 + (e_kin*1e6*qe)/(mp*clight*clight) 
beta0  = sqrt((gamma0*gamma0-1.0)/(gamma0*gamma0)) 
t_turn = circum/(beta0*clight)
delAmp =(BumpA/Bumpf*t_turn)
# Qhor
#q=4.02;  dq=0.01; qmax=4.02  
x0=0.07; dx0=0.001; x0max=0.09
path='/d/bhs01/appel/patric'  # path for output
expath='/u/sappel/codes/patric_mti'  # base directory of code
subdir = '/position_modi_360/'   

# clear target directory 
if(os.path.exists(path+subdir)):
    print "exists"
else:
	os.mkdir(path+subdir)   

#assert os.path.exists(path+subdir),  os.mkdir(path+subdir) 
scan=str(x0)
patric_dict=OrderedDict()

patric_dict['NPIC']=10000  # particles per beamlett; SP
patric_dict['NX']=128
patric_dict['NY']=128
patric_dict['NZ']=256
patric_dict['NZ_bunch']=64
patric_dict['cells']=18
patric_dict['lossTol']=1.2  # tolerable relative losses; SP
patric_dict['e_kin']=e_kin 
patric_dict['Z']=7.0
patric_dict['A']=14.0
patric_dict['current']=21e-6 # [A]
patric_dict['piperadius']=0.1
patric_dict['coll_halfgap']=0.07  # distance of septum
patric_dict['image_x']=0.1*1.0  # horizontal boundary for Green's function i.e. image currents; 0 means open boundary; SP
patric_dict['image_y']=0.07*1.0
patric_dict['circum']=circum 
patric_dict['gamma_t']=4.79
patric_dict['CF_advance_h']=97.3*pi/180.0
patric_dict['CF_advance_v']=97.3*pi/180.0 # 129.3*pi/180.0
patric_dict['CF_R']=0.0
patric_dict['CF_length']=18.0
patric_dict['NCF']=16
patric_dict['koct']=8.0
patric_dict['dQxm']=-0.0054*1.0
patric_dict['dQym']=-0.0054*1.0
patric_dict['dqx_detune']=-0.0025/6.0  # rms
patric_dict['dqy_detune']=-0.0025/6.0  # rms
patric_dict['pic_subset']=1  # 10000
patric_dict['init_pic_xy']=1  # 0 (WB), 1 (KV), 2 (SG), 3 (GS)
patric_dict['init_pic_z']=2  # 0 (elliptic, coast), 1 (elliptic, bunch), 2 (Gauss, coast)  # was 4 !
patric_dict['momentum_spread']=0.5e-3  # rms
patric_dict['rms_emittance_x0']=ex/4.  # rms, mm mrad
patric_dict['rms_emittance_y0']=ey/4.  # rms, mm mrad
patric_dict['mismatch_x']=1.0  # initial beam size does not match beta function (1 = match)
patric_dict['mismatch_y']=1.0
patric_dict['x_septum']=0.07  # distance of septum from nominal orbit
patric_dict['offcenter']=x0   #BumpA-BumpA/Bumpf*ChopV
patric_dict['inj_angle']=7.25e-3  # injection angle in rad; added by SP
patric_dict['max_inj']=mtj;
patric_dict['bumpI']=bumpI;   # 0  (version SP) , 1 (flexibility)  ; SA
patric_dict['amp0']=Amp0    # amplitude  bump if bumpI=1; added by SA
patric_dict['ampp0']=6e-3  # amplitude prime bump if bumpI=1; added by SA  
patric_dict['delAmp']=delAmp #(BumpA/Bumpf*t_turn)   # ramp rate if bumpI=1; added by SA
patric_dict['bunchfactor']=1.0
patric_dict['dqci']=0.0
patric_dict['dqcr']=-0.0
patric_dict['Rs']=0.0 # 8.5e6
patric_dict['nres']=10.0
patric_dict['Qs']=2.0
patric_dict['leit']=1.0e5
patric_dict['Zimage']=-7.0e7 # -3.0e7/2.0
patric_dict['madx_input_file']=1  # 0 (no, i.e. use CF), 1 (yes)
patric_dict['space_charge']=0  #   0 (off), 1 (self-consistent), 2 (linear), 3 (nonlinear)
patric_dict['imp_kick']=0       #   0 (off), 1 (on)
patric_dict['sliced']=0
patric_dict['cavity']=0   # 0 (off), 1 (rf), 2 (barrier)
patric_dict['octupole_kick']=0
patric_dict['ampdetun_kick']=0  # commented out in Main.cpp, as incompatibel with modifications
patric_dict['chroma']=1  # 0 off, 1 on
patric_dict['bc_end']=1  # 0 (open), 1 (periodic)
patric_dict['print_cell']=1  # with MADX input one cell is the whole lattice
patric_dict['footprint']=0
patric_dict['btf']=0
patric_dict['btf_harmonic']=0
patric_dict['ausgabe']=path+subdir+scan   
patric_dict['input']=expath                                      # wo liegt patric

# machine counter	
machine_id=0 
x0min=x0
while x0min <= x0max:
	
   	runid=1            # run identification number (only out.dat)      
 
	#--------------------- MPI options ------------------------
	num_procs=1            # 1 (scan), 2 (no scan) 

	#---------------------------scan parameter-----------------------------
	scan=str(x0min)
	patric_dict['ausgabe']=path+subdir+scan		

	#-------------------------- run ------------------------------   
   	if(os.path.exists(path+subdir+scan)):
	    print "exists"
	else: 
		os.mkdir(path+subdir+scan)
  
    # produce twiss and sectormap file     
	os.chdir(expath+"/mad/")         
	cmd='sed -ie "s/Q_hor = .*/Q_hor = %1.2f;/" sis18_inj.mad;' %(Q_hor) 
	os.system(cmd)
	cmd='sed -ie "s/Q_vec = .*/Q_vec = %1.2f;/" sis18_inj.mad;' %(Qvec) 
	os.system(cmd)   
	os.system("/u/sappel/.madx/madx < sis18_inj.mad >out.madx.dat")
	shutil.copy("/u/sappel/codes/patric_mti/mad/sis18_inj.mad", path+subdir+scan+"/")
	shutil.copy("/u/sappel/codes/patric_mti/mad/out.madx.dat", path+subdir+scan+"/")
	os.chdir("../python/")
                                
	patric_dict['offcenter']=x0min 
	# create new PATRIC configuration file
	outfile=open('%s/patric.cfg' % (path+subdir+scan),'w')
	for x,y in patric_dict.iteritems():
	   outfile.write('%s:  %s\n' %(x,y))
	outfile.close()


	# mpi machine configuation file:
	machinefile=open('machines.tmp','w')
	machinefile.write('%s:4\n' %(machines[machine_id]))
	machinefile.close()

 	# start PATRIC/mpi job:       # mpiexec path: /d/bhs01/local/mpich2-install/bin   
	thishost=gethostname()
	if thishost == machines[machine_id]:
		cmd="nohup mpiexec -n %d -machinefile machines.tmp %s/bin/patric %s/patric > %s/out.%s.dat &" %(num_procs,expath,path+subdir+scan,path+subdir+scan,runid)
	else:
		cmd="nohup ssh %s /d/bhs01/local/mpich2-install/bin/mpiexec -n %d -machinefile /u/sappel/codes/patric_mti/python/machines.tmp %s/bin/patric %s/patric > %s/out.%s.dat &" %(machines[machine_id],num_procs,expath,path+subdir+scan,path+subdir+scan,runid)   
   	os.system(cmd)
	os.system("echo ' '")  # release prompt (optically) 
	# next machine
	time.sleep(1) 
	machine_id+=1
	x0min  +=dx0 
	

#-------------------------------------------------------------

import os
import string
import sys
from numpy import *
from scipy.fftpack import fft, fftfreq    
import matplotlib.pyplot as plt   

#12   
scan=12
current_scan=array([26,30,34,38,42,46,49,53,56,60,64,68])  
NPIC_scan=zeros(scan,float)  
Int_scan=zeros(scan,float)
Int_mtj_scan=zeros(scan,float)
Loss_scan=zeros(scan,float)
emitx_scan=zeros(scan,float)
emity_scan=zeros(scan,float)   
dqx_scan=zeros(scan,float)
dqy_scan=zeros(scan,float)

# input     
run="proton_"+str(current_scan[0])+"/"      
inputdir="/d/bhs01/appel/patric/"+run

# read idl.dat
datafile=open('%sidl.dat' %(inputdir),'r')
lines=datafile.readlines()        
datafile.close()
e_kin=float(lines[1])       
current=float(lines[4]) 
circum=float(lines[5])       
cell_length=float(lines[7])  

# phy. values
qe=1.6022e-19 
rp=1.55e-18 
mp=1.6605e-27
clight=2.99e8
gamma0 = 1.0 + (e_kin*1e6*qe)/(mp*clight*clight) 
beta0  = sqrt((gamma0*gamma0-1.0)/(gamma0*gamma0))   
freq0=beta0*clight/circum
gammat=5.4   # !!!!!!
eta0=1.0/gammat**2-1.0/gamma0**2
N0=current*circum/(qe*beta0*clight)  

ex=150e-6*2/3
ey=50e-6*2/3
N=0.5*0.35*(pi*beta0**2*gamma0**3)*(ey+sqrt(ex*ey))/rp/1.4

s =0   

while s <= scan-1: 
	# read idl.dat
	datafile=open('%sidl.dat' %(inputdir),'r')
	lines=datafile.readlines()        
	datafile.close()  
	current=float(lines[4]) 
	N0=current*circum/(qe*beta0*clight)  
 
	   
	# read patric.dat 
	run="proton_"+str(current_scan[s])+"/"      
	inputdir="/d/bhs01/appel/patric/"+run
	datafile=open('%spatric.dat' %(inputdir),'r')
	lines=datafile.readlines()
	datafile.close()
	Np=len(lines)
    
    
	for j in xrange(Np):
		words=string.split(lines[j])
		emitx_scan[s]=float(words[1])
		emity_scan[s]=float(words[2]) 
		if j == 0:
			NPIC_scan[s]=N0/float(words[10])
		Int_scan[s]=float(words[10])  
		Loss_scan[s]=float(words[18])
		Int_mtj_scan[s]=float(words[19]) 
	

	dqx_scan=Int_scan*NPIC_scan*rp*1.4/(pi*beta0**2*gamma0**3)*1/(emitx_scan+sqrt(emitx_scan*emity_scan))/4       
	dqy_scan=Int_scan*NPIC_scan*rp*1.4/(pi*beta0**2*gamma0**3)*1/(emity_scan+sqrt(emitx_scan*emity_scan))/4   
	s +=1 
	
    


# output per current    
fig = plt.figure()      
fig.subplots_adjust(wspace=0.5, hspace=0.3)
ax = fig.add_subplot(222)
#loss_plot=plt.plot(s/cell_length,100*(1-Int/Int_mtj))
plt.plot(current_scan,dqy_scan, '-o',label='y')  
plt.plot(current_scan,dqx_scan, '-o', label='x')
plt.ylim(0,0.4) 
plt.axhline(0.5*0.36)   
ax.legend(loc=0)        
plt.xlabel("currents")
plt.ylabel(r"$\Delta$ Q")   


ax=fig.add_subplot(223) 
#lossloc_plot=plt.plot(s/cell_length,loss*NPIC)  
loss_plot=plt.plot(current_scan,NPIC_scan*(Int_mtj_scan-Int_scan), '-o')     
plt.xlabel("currents")
plt.ylabel("Particle loss")   

ax=fig.add_subplot(224) 
emitx_plot=plt.plot(current_scan,emitx_scan*4/1e-6, '-o',label="hor")
emitx_plot=plt.plot(current_scan,emity_scan*4/1e-6,'-o', label="ver")      
plt.axhline(100)
plt.axhline(33)    
#plt.ylim(0,150)
ax.legend(loc=0)
plt.xlabel("currents")
plt.ylabel(r"$\epsilon$ [mm mrad]")   

 
ax=fig.add_subplot(221)     
inj_plot=plt.plot(current_scan,Int_mtj_scan*NPIC_scan,'-o', label="injected")
store_plot=plt.plot(current_scan,Int_scan*NPIC_scan, '-o',label="store")   
plt.xlabel("currents")
plt.ylabel("Intensity")  
plt.axhline(N)   
ax.legend(loc=0)
plt.ylim([0,6e12])


print N
plt.show()     







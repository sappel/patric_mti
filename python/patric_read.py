import os
import string
import sys
from numpy import *
from scipy.fftpack import fft, fftfreq    
import matplotlib.pyplot as plt  

# input     
run="proton_56/"     
inputdir="/d/bhs01/appel/patric/"+run
inputdir='/d/bhs01/appel/patric/tmp/' 
# read idl.dat

datafile=open('%sidl.dat' %(inputdir),'r')
lines=datafile.readlines()        
datafile.close()

numprocs=int(lines[0])
e_kin=float(lines[1])   
Zb=float(lines[2])      
Ab=float(lines[3])      
current=float(lines[4]) 
circum=float(lines[5])       
Nelements=float(lines[6])   
cell_length=float(lines[7])  
pipe_radius=float(lines[9]) 
NX=int(lines[10]) 
NY=int(lines[11]) 
NZ=int(lines[12]) 
tunex=float(lines[13])
tuney=float(lines[14])
fsyn=float(lines[15])

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


# read patric.dat

datafile=open('%spatric.dat' %(inputdir),'r')
lines=datafile.readlines()
datafile.close()

Np=len(lines)

s=zeros(Np,float)
emitx=zeros(Np,float) # rms emittance x
emity=zeros(Np,float) 
x_max=zeros(Np,float)
Int=zeros(Np,float)
Int_mtj=zeros(Np,float)
loss=zeros(Np,float)

for j in xrange(Np):
	words=string.split(lines[j])
	s[j]=float(words[0])
	emitx[j]=float(words[1])
	emity[j]=float(words[2])
	x_max[j]=float(words[3])
	Int[j]=float(words[10])
	loss[j]=float(words[18])
	Int_mtj[j]=float(words[19])   

NPIC=N0/Int[0]   
                             

dqx=Int*NPIC*rp*1.4/(pi*beta0**2*gamma0**3)*1/(emitx+sqrt(emitx*emity))/4       
dqy=Int*NPIC*rp*1.4/(pi*beta0**2*gamma0**3)*1/(emity+sqrt(emitx*emity))/4                  

# output per turn    
fig = plt.figure()      
fig.subplots_adjust(wspace=0.5, hspace=0.3)
ax = fig.add_subplot(222)
#loss_plot=plt.plot(s/cell_length,100*(1-Int/Int_mtj))
#plt.plot(s/cell_length,dqy, label='y')  
#plt.plot(s/cell_length,dqx, label='x') 
#plt.axhline(0.5*0.36)  
#plt.axhline(0.5)   
plt.plot(s/cell_length,x_max, label="ver") 
ax.legend(loc=0)        
plt.xlabel("turns")
plt.ylabel(r"$\Delta$ Q")   


ax=fig.add_subplot(223) 
plt.plot(s/cell_length,loss*NPIC, label="on septum")  
plt.plot(s/cell_length,NPIC*(Int_mtj-Int), label="total")     
ax.legend(loc=0) 
plt.xlabel("turns")
plt.ylabel("Particle loss")   

ax=fig.add_subplot(224) 
emitx_plot=plt.plot(s/cell_length,emitx*4/1e-6, label="hor")
emitx_plot=plt.plot(s/cell_length,emity*4/1e-6, label="ver")      
#plt.plot(s/cell_length,x_max, label="ver") 
plt.axhline(100)
plt.axhline(33)    
plt.ylim(0,300)
ax.legend(loc=0)
plt.xlabel("turns")
plt.ylabel(r"$\epsilon$ [mm mrad]")   
ax=fig.add_subplot(221)     
inj_plot=plt.plot(s/cell_length,Int_mtj*NPIC, label="injected")
store_plot=plt.plot(s/cell_length,Int*NPIC, label="store")   
plt.xlabel("turns")
plt.ylabel("Intensity")  
plt.axhline(3.5e12)   
plt.title(run)   
ax.legend(loc=0)
#plt.ylim([0,6e12])


plt.show()     


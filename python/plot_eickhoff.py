#!/usr/bin/python
"""
A python script to create a tune diagramm    
"""      
  

import sys, math
from math import *    
from numpy import *   
import matplotlib.pyplot as mplplot
import numpy as np    
#import numarray 
from matplotlib import colors, ticker, cm, mpl   
from matplotlib.font_manager import FontProperties 


mplplot.matplotlib.rc('font', size = 20) 
fig = mplplot.figure()   
ax = fig.add_subplot(111)   
                                 

# constants
e_kin = 11.582         #11.4
mp = 1.6605e-27
qe = 1.6022e-19
clight = 2.988e8
circum = 216.72    
current=2e-3      
Z=7
gamma0 = 1.0 + (e_kin*1e6*qe)/(mp*clight*clight) 
beta0  = sqrt((gamma0*gamma0-1.0)/(gamma0*gamma0))
Ni = current*circum/(beta0*clight*qe*Z)
Ni = 8.25133e9

# loss sim. low 
Bumpfmin=20; dBumpf=20; Bumpfmax=360     # Madx: Q_hor; step, end   
 
dir='/d/bhs01/appel/patric/eickhoff/neu_nsp/'#'/ex_1/'                
scan=str(Bumpfmin) 
path=dir+'ex_6/'  
filename=path+scan+'/idl.dat' 
file = open(filename, 'r') 
N = []   
for line in file:
	 pair = line.split() 
	 N.append(float(pair[0]))   
file.close()     


cells=N[23] 


#path='/d/bhs01/appel/patric/uran_flex/' 

Int =  []
loss = []
pos = []
Bumpf=Bumpfmin
while Bumpf <= Bumpfmax:
	scan=str(Bumpf) 
	filename=path+scan+'/mti.dat' 
	#print filename
	file = open(filename, 'r')
	k=0
	while k < cells-1:                       
		file.readline() 
		k +=1
	for line in file:
		pair = line.split() 
		Int.append(float(pair[1])) 
		loss.append((float(pair[2])))
		pos.append(Bumpf) 
	file.close() 
	Bumpf +=dBumpf
                   
for i in range(len(loss)):
	loss[i]= (1-loss[i])*100     
	
for i in range(len(loss)):
	pos[i]= (pos[i])

                   
for i in range(len(loss)):
	Int[i]= (Int[i]/Ni)


mplplot.plot(pos,Int,'h-', label=r'$\epsilon_{x}=$ 6 mm mrad')         

path='/d/bhs01/appel/patric/eickhoff/neu_sp/ex_6/'
Int =  []
loss = []
pos = []
Bumpf=Bumpfmin
while Bumpf <= Bumpfmax:
	scan=str(Bumpf) 
	filename=path+scan+'/mti.dat' 
	#print filename
	file = open(filename, 'r')
	k=0
	while k < cells-1:                       
		file.readline() 
		k +=1
	for line in file:
		pair = line.split() 
		Int.append(float(pair[1])) 
		loss.append((float(pair[2])))
		pos.append(Bumpf) 
	file.close() 
	Bumpf +=dBumpf
                   
for i in range(len(loss)):
	loss[i]= (1-loss[i])*100     
	
for i in range(len(loss)):
	pos[i]= (pos[i])

                   
for i in range(len(loss)):
	Int[i]= (Int[i]/Ni)


mplplot.plot(pos,Int,'D--', label=r'$\epsilon_{x}=$ 6 mm mrad sp')
 
#dir='/d/bhs01/appel/patric/eickhoff/neu_nsp/'   
path='/d/bhs01/appel/patric/eickhoff/neu_sp/ex_1/'  #dir#+'ex_6/' 
Int =  []
loss = []
pos = []
Bumpf=Bumpfmin
while Bumpf <= Bumpfmax:
	scan=str(Bumpf) 
	filename=path+scan+'/mti.dat' 
	#print filename
	file = open(filename, 'r')
	k=0
	while k < cells-1:                       
		file.readline() 
		k +=1
	for line in file:
		pair = line.split() 
		Int.append(float(pair[1])) 
		loss.append((float(pair[2])))
		pos.append(Bumpf) 
	file.close() 
	Bumpf +=dBumpf
                   
for i in range(len(loss)):
	loss[i]= (1-loss[i])*100     
	
for i in range(len(loss)):
	pos[i]= (pos[i])

                   
for i in range(len(loss)):
	Int[i]= (Int[i]/Ni)


mplplot.plot(pos,Int,'h--', label=r'$\epsilon_{x}=$ 1 mm mrad sp')      

         

path='/d/bhs01/appel/patric/eickhoff/neu_nsp/ex_2/'#dir+'ex_10/' 
Int =  []
loss = []
pos = []
Bumpf=Bumpfmin
while Bumpf <= Bumpfmax:
	scan=str(Bumpf) 
	filename=path+scan+'/mti.dat' 
	#print filename
	file = open(filename, 'r')
	k=0
	while k < cells-1:                       
		file.readline() 
		k +=1
	for line in file:
		pair = line.split() 
		Int.append(float(pair[1])) 
		loss.append((float(pair[2])))
		pos.append(Bumpf) 
	file.close() 
	Bumpf +=dBumpf
                   
for i in range(len(loss)):
	loss[i]= (1-loss[i])*100     
	
for i in range(len(loss)):
	pos[i]= (pos[i])

                   
for i in range(len(loss)):
	Int[i]= (Int[i]/Ni)


mplplot.plot(pos,Int,'D-', label=r'$\epsilon_{x}=$ 1 mm mrad ')



#18 dy=20 dx=100 m=0.2
theorie=arange(20,380,20)  


mplplot.plot(theorie,theorie*0.22,'ks-', label=r'no loss')      

#theory=3*np.array(pos)  
#mplplot.plot(pos,theory)     


fontP = FontProperties()
fontP.set_size('small')  
ax.legend(loc=0, prop = fontP)	  
#ax.get_frame().set_visible(False) 
mplplot.axis([0, 365, 0,50])       
#mplplot.title('without space charge effects')  
#mplplot.annotate(r'$\epsilon_{x,1}$=1 mm mrad',xy=(200,5))     
#mplplot.annotate(r'$\epsilon_{x,6}$=6 mm mrad',xy=(200,1))   
#mplplot.annotate("(2,2,1)",xy=(qx1,(p-l*qx1)/m+0.01), rotation=theta,fontsize=12, color='w')    
mplplot.xlabel(r"Bump ramp rate $(\mu s)$" )
mplplot.ylabel(r"MTI efficiency ($N/N_0$)")       

#mplplot.savefig('test.eps', transparent=True)
mplplot.show()
























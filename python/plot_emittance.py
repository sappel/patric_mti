#!/usr/bin/python
"""
A python script to create a tune diagramm, need ABLASS output     
"""      
  

import sys, math
from math import *    
from numpy import *    
import matplotlib.pyplot as mplplot
import numpy as np    
import numarray 
from matplotlib import colors, ticker, cm, mpl 

mplplot.matplotlib.rc('font', size = 18) 
fig = mplplot.figure()   
ax = fig.add_subplot(111)   


# loss sim. low 
exmin=2; dxex=1; exmax=10     # Madx: Q_hor; step, end  
                
scan=str(exmin) 
path='/d/bhs01/appel/patric/emittance/low/' 
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
ex=exmin
while ex <= exmax:
	scan=str(ex) 
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
		pos.append(ex) 
	file.close() 
	ex +=dxex
                   
for i in range(len(loss)):
	loss[i]= (1-loss[i])*100     
	
for i in range(len(loss)):
	pos[i]= (pos[i])


mplplot.plot(pos,loss,'s-', label='low')  



path='/d/bhs01/appel/patric/emittance/high/' 

Int =  []
loss = []
pos = []
ex=exmin   
while ex <= exmax:
	scan=str(ex) 
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
		pos.append(ex) 
	file.close() 
	ex +=dxex
                   
for i in range(len(loss)):
	loss[i]= (1-loss[i])*100     
	
for i in range(len(loss)):
	pos[i]= (pos[i])




mplplot.plot(pos,loss,'ro--', label='high')   

#theory=3*np.array(pos)  
#mplplot.plot(pos,theory)
ax.legend(loc=0	)  
#mplplot.title('with orbit deformation')
mplplot.xlabel(r"hor. emittance / mm mrad")
mplplot.ylabel(r"Loss / %")       

#mplplot.savefig('test.eps', transparent=True)
mplplot.show()
























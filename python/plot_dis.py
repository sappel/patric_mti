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
x0min=0.07;  dx0=0.002; x0max=0.09    # Madx: Q_hor; step, end  
                
scan=str(x0min) 
path='/d/bhs01/appel/patric/'#position/emittance10/'
folder='position_modi_360/'
path=path+folder 
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
x0=x0min
while x0 <= x0max:
	scan=str(x0) 
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
		pos.append(x0) 
	file.close() 
	x0 +=dx0
                   
for i in range(len(loss)):
	loss[i]= (1-loss[i])*100     
	
for i in range(len(loss)):
	pos[i]= (pos[i])*1e3


mplplot.plot(pos,loss,'s-', label='360$\mu s$')

"""
path='/d/bhs01/appel/patric/position/emittance10/'
folder='position_modi_160/'
path=path+folder

Int =  []
loss = []
pos = []
x0=x0min
while x0 <= x0max:
	scan=str(x0) 
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
		pos.append(x0) 
	file.close() 
	x0 +=dx0
                   
for i in range(len(loss)):
	loss[i]= (1-loss[i])*100     
	
for i in range(len(loss)):
	pos[i]= (pos[i])*1e3

mplplot.plot(pos,loss,'ro--', label='160$\mu s$')   

path='/d/bhs01/appel/patric/position/emittance10/'
folder='position_modi_140/'
path=path+folder

Int =  []
loss = []
pos = []
x0=x0min
while x0 <= x0max:
	scan=str(x0) 
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
		pos.append(x0) 
	file.close() 
	x0 +=dx0
                   
for i in range(len(loss)):
	loss[i]= (1-loss[i])*100     
	
for i in range(len(loss)):
	pos[i]= (pos[i])*1e3

mplplot.plot(pos,loss,'mo--', label='140$\mu s$')

path='/d/bhs01/appel/patric/position/emittance10/'
folder='position_modi_120/'
path=path+folder

Int =  []
loss = []
pos = []
x0=x0min
while x0 <= x0max:
	scan=str(x0) 
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
		pos.append(x0) 
	file.close() 
	x0 +=dx0
                   
for i in range(len(loss)):
	loss[i]= (1-loss[i])*100     
	
for i in range(len(loss)):
	pos[i]= (pos[i])*1e3

mplplot.plot(pos,loss,'co--', label='120$\mu s$')


path='/d/bhs01/appel/patric/position/emittance5/'
folder='position_modi_80/'
path=path+folder

Int =  []
loss = []
pos = []
x0=x0min
while x0 <= x0max:
	scan=str(x0) 
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
		pos.append(x0) 
	file.close() 
	x0 +=dx0
                   
for i in range(len(loss)):
	loss[i]= (1-loss[i])*100     
	
for i in range(len(loss)):
	pos[i]= (pos[i])*1e3


mplplot.plot(pos,loss,'gD-.', label='80$\mu s$')  
"""




mplplot.axis([70, 90, 0,100])                          
ax.legend(loc=0	)  
#mplplot.title('with orbit deformation')
mplplot.xlabel(r"inj beam position / mm")
mplplot.ylabel(r"Loss / %")   


mplplot.show()	


""" 
numproc=N[0]
e_kin=N[1]
qb=N[2] 
mb=N[3]
current=N[4]
C=N[5]
Nelements=N[6] 
cell_length=N[7]    
Np=[8]
pipe_radius=N[9] 
NX=N[10]
NY=N[11]
NZ=N[12]
print_cel=N[13]
Q_hor=N[14]
Q_vec=N[15]
betax=N[16]
alpx=N[17]
x0=N[18]
x0prime=N[19]   
amp0=N[20]    
ampp0=N[21]
delAmp=N[22]  
cells=N[23] 
max_inj=N[24]
filename=path+'mti.dat' 
file = open(filename, 'r') 
turn = []  
Int =  []
loss = []
for line in file:
	 pair = line.split() 
	 turn.append(float(pair[0])) 
	 Int.append(float(pair[1])) 
	 loss.append(float(pair[2]))   
file.close()
"""



























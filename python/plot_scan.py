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

mplplot.matplotlib.rc('font', size = 25) 
fig = mplplot.figure()   
ax = fig.add_subplot(111)   
Qh_shift=0.04

# tune exp.
tune_exp=np.array([4.02, 4.04, 4.06, 4.08, 4.1, 4.12, 4.13, 4.14, 4.15, 4.16, 4.17, 4.18, 4.2, 4.22, 4.24, 4.26, 4.28, 4.29, 4.31, 4.33])

#loss exp. low
loss_low=np.array([67.16, 54.65, 38.79, 27.06, 24.54, 16.60, 20.60, 22.01, 22.36, 13.74, 20.19, 20.62, 25.10, 15.49, 17.99, 23.99, 30.70, 31.49, 32.71, 28.98])   

#mplplot.axvline(x=4.24, color='k') 
#mplplot.axvline(x=4.33, color='k') 

#mplplot.plot(tune_exp, loss_low,'go-', label='Exp.') 
#mplplot.errorbar(tune_exp, loss_low, xerr=0, yerr=5, fmt=None, ecolor='g')

mplplot.plot(tune_exp+Qh_shift, loss_low,'kD-', label='Exp + $\Delta Q$') 
mplplot.errorbar(tune_exp+Qh_shift, loss_low, xerr=0, yerr=5, fmt=None, ecolor='k')
                                                                                 
# loss exp. high
loss_high=100-np.array([69.34, 60.91, 56.30, 43.38, 35.91, 34.78, 33.26, 31.24, 31.72, 31.54, 30.27, 29.18, 32.69, 34.63, 31.18, 31.19, 37.97, 31.46, 46.43, 29.00])
 
#mplplot.plot(tune_exp+Qh_shift, loss_high,'rs--', label='high')   
#mplplot.errorbar(tune_exp+Qh_shift, loss_high, xerr=0, yerr=5, fmt=None, ecolor='r') 



# loss sim. low              q=4.02;  dq=0.025; qmax=4.5 
qmin=4.02;  dq=0.02; qmax=4.33    # Madx: Q_hor; step, end  
               
scan=str(qmin) 
path='/d/bhs01/appel/patric/orbit_test/meas_11/'  
#path='/d/bhs01/appel/patric/exp/high/meas_20/'
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
turn = []
q=qmin
while q <= qmax:
	scan=str(q) 
	filename=path+scan+'/mti.dat' 
	#print filename
	file = open(filename, 'r')
	k=0
	while k < cells-1:                      # line after 3 S09DT_ML 
		file.readline() 
		k +=1
	for line in file:
		pair = line.split() 
		Int.append(float(pair[1])) 
		loss.append((float(pair[2])))
		turn.append(q) 
	file.close() 
	q +=dq
                   
for i in range(len(loss)):
	loss[i]= (1-loss[i])*100#-30
                             

mplplot.plot(turn, loss,'rD--', label='x=0') 


# loss sim. high
#path='/d/bhs01/appel/patric/uran_sp/'
#path='/d/bhs01/appel/patric/neon_high2_no/'  
 
scan=str(qmin)
path='/d/bhs01/appel/patric/orbit_test/meas_12/' 
#path='/d/bhs01/appel/patric/exp/high/meas_1/'         
filename=path+scan+'/idl.dat' 
file = open(filename, 'r') 
N = []   
for line in file:
	 pair = line.split() 
	 N.append(float(pair[0]))   
file.close()     


cells=N[23]


q=qmin

Int =  []
loss = []
turn = []
while q <= qmax:
	scan=str(q) 
	filename=path+scan+'/mti.dat' 
	file = open(filename, 'r') 
	k=0
	while k < cells-1:                      # line after 3 S09DT_ML 
		file.readline() 
		k +=1
	for line in file:
		pair = line.split() 
		Int.append(float(pair[1])) 
		loss.append((float(pair[2]))) 
		turn.append(q)   
	file.close() 
	q +=dq
                   
for i in range(len(loss)):
      loss[i]= (1-loss[i])*100      

 
mplplot.plot(turn,loss,'bs--', label='x=+5mm') 

scan=str(qmin)
path='/d/bhs01/appel/patric/orbit_test/meas_13/' 
#path='/d/bhs01/appel/patric/exp/high/meas_1/'         
filename=path+scan+'/idl.dat' 
file = open(filename, 'r') 
N = []   
for line in file:
	 pair = line.split() 
	 N.append(float(pair[0]))   
file.close()     


cells=N[23]


q=qmin

Int =  []
loss = []
turn = []
while q <= qmax:
	scan=str(q) 
	filename=path+scan+'/mti.dat' 
	file = open(filename, 'r') 
	k=0
	while k < cells-1:                      # line after 3 S09DT_ML 
		file.readline() 
		k +=1
	for line in file:
		pair = line.split() 
		Int.append(float(pair[1])) 
		loss.append((float(pair[2]))) 
		turn.append(q)   
	file.close() 
	q +=dq
                   
for i in range(len(loss)):
      loss[i]= (1-loss[i])*100      

 
mplplot.plot(turn,loss,'gs--', label='x=-5mm')





mplplot.axis([3.98, 4.4, 0,100])                          
ax.legend(loc=0)  
#mplplot.title(folder) 
ax.set_xticks([4.00,4.1,4.2,4.3,4.4])
mplplot.xlabel(r"Horizontal tune")
mplplot.ylabel(r"Loss in %")   

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



























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

mplplot.matplotlib.rc('font', size = 20) 
fig = mplplot.figure()   
ax = fig.add_subplot(111)


# rampe paret
                  

path='/d/bhs01/appel/patric/uran/'
scan='4.11' 
filename=path+scan+'/idl.dat' 
file = open(filename, 'r') 
N = []   
for line in file:
	 pair = line.split() 
	 N.append(float(pair[0]))   
file.close()     

amp0=N[20]    
ampp0=N[21]
delAmp=N[22]
cells=N[23]

x=arange(0,25,1)
rampe=-x*delAmp+amp0
mplplot.plot(x,rampe/1e-3, label='rampe paret')   

# rampe flexibility
path='/d/bhs01/appel/patric/uran_flex/'
scan='4.11' 
filename=path+scan+'/idl.dat' 
file = open(filename, 'r') 
N = []   
for line in file:
	 pair = line.split() 
	 N.append(float(pair[0]))   
file.close()
   

amp0=N[20]    
ampp0=N[21]
delAmp=N[22]
cells=N[23] 

rampe=-x*delAmp+amp0 
mplplot.plot(x,rampe/1e-3, 'r--', label='rampe flex')
mplplot.axis([0, 25, 0,100])
ax.legend() 
mplplot.ylabel(r"Amplitude / mm")
mplplot.xlabel(r"Turn")
mplplot.show()	






























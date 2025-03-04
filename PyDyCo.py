# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 14:14:16 2022

@author: Mahdi-Amine MAMOUNI
#test

"""
# -*- coding: utf-8 -*-
"""
The first pydco was created on : 
Created on Fri May 31 15:55:48 2019
by: 
@author: physique
"""

"""++++""""""+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++""""""++++"""
"""++++""""""+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++""""""++++"""
# ==================================================Library importation  ===============================================
"""++++""""""+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++""""""++++"""
"""++++""""""+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++""""""++++"""
from PySpice.Spice.Netlist import SubCircuitFactory, SubCircuit
from scipy import signal
from math import pi
from scipy import interpolate
import math
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from scipy.signal import *
from numpy import *
import scipy.fftpack
from scipy.integrate import quad
from scipy.signal.signaltools import correlate
import statistics
import matplotlib.gridspec as gridspec
from decimal import Decimal
from PySpice.Unit import *
from PySpice.Spice.Netlist import Circuit
import numpy as np
import matplotlib.pyplot as plt
import PySpice.Logging.Logging as Logging
logger = Logging.setup_logging()
from foncdyco.foncdycobj  import *
import time 

"""++++""""""+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++""""""++++"""
"""++++""""""+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++""""""++++"""
start = time.time()
    
#========================Boundatry conditions =================================
#Define limit condition: Neuman/Dirichlet with Voltage and Current 
BC='Voltage'
#BC='Current'                                                                 #
#Amplitude:                                                                   #
A=  0.01                                                                      #
# Define frequency by uncommenting:                                           #
#fqc=0.001                              #For a specific value: works with 'transient' simulation 
#fqcrange=np.geomspace(0.0001,0.1,50)   #For a specific range of values: works with 'transient' simulations                                    #
fqc='dc'                                #For DC: works with 'transient' and 'operating point' simulations                               #
#Define temperature 
T1= 300                                                                       #
T2= 300                                                                       #
T0=300              #Initial temperature                                                          #
                                                                              #
#==============================================================================


#=======================Type of simulation ====================================
#For transient, it's possible to have both transient DC and transient AC. You can define the transient condition in Sim, line 120.

#Sim ='TRANSIENT'        #if transient then uncomment this line

Sim = 'Operating Point'  #if Op.point  then uncomment this line
#==============================================================================

#If you work in a range of frequency domain, uncomment the next line and check the indentation.

#for fqc in fqcrange:    #DO NOT FORGET THE IDENTATION!!!!!!
circuit = Circuit("thermoélectrical Network")

#=======================Type of materials in matériaux.txt======================================
#circuit.subcircuit(OnsaNod('Name', S[V/K], Sigma[Sm^-1], Kappa[Wk^-1m^-1] ,Cp[J/(Km^3)],TO [K] , L[m], A[m²], rho[kg/m^3]))
#you can add as much materias as you want 
circuit.subcircuit( OnsaNod('m',0, 0.001,   0.2, 1, T0,1,1,1))
circuit.subcircuit( OnsaNod('p',200e-6, 10, 0.2, 1, T0,1,1,1))
circuit.subcircuit( OnsaNod('n',-200e-6, 10, 0.2, 1, T0,1,1,1))
"""
circuit.subcircuit( OnsaNod('e',20e-6, 1e6, 57.3, 3.7e6, T0,5e-6,1e-12,8500) )
circuit.subcircuit( OnsaNod('s',0, 1e6, 7.3, 3.7e6, T0,5e-6,1e-12,8500) )
circuit.subcircuit( OnsaNod('k',0, 1e6, 100000007.3, 3.7e6, T0,5e-6,1e-12,8500) )
circuit.subcircuit( OnsaNod('i',0, 1e6, 0.0000001, 3.7e6, T0,5e-3,1e-12,8500) )
circuit.subcircuit( OnsaNod('r',0, 2.55e-4, 148,  1.631e6, T0,2.5e-7,1e-12,2330) )
#circuit.subcircuit( OnsaNod('M',0, 1e6, 1.7, 3.7e6, T0,5e-10,1e-12,8500) )
#circuit.subcircuit( OnsaNod('P',0, 1e6, 1.3, 3.7e6, T0,5e-10,1e-12,8500) )
circuit.subcircuit( OnsaNod('E',0, 1e6, 57.3, 3.7e6, T0,5e-10,1e-12,8500) )
circuit.subcircuit( OnsaNod('S',0, 1e6, 0.03, 3.7e6, T0,5e-10,1e-12,8500) )
circuit.subcircuit( OnsaNod('K',0, 1e6, 107.3, 3.7e6, T0,5e-10,1e-12,8500) )
circuit.subcircuit( OnsaNod('I',0, 1e6, 0.000001, 3.7e6, T0,5e-10,1e-8,8500) )

#circuit.subcircuit( OnsaNod('c',0e-6, 1e6, 17.3, 3.7e6, T0,0.2e-2,1e-8,8500) )

#circuit.subcircuit( OnsaNod('c',35e-6, 1e6, 22.7, 3.5e6, T0,0.2e-2,1e-8,8910) ) 
circuit.subcircuit( OnsaNod('c',35e-6, 1e6,  22.7, 3.5e6, T0,0.2e-2,1e-8,8500) )
circuit.subcircuit(OnsaNod('s',0, 1e-6, 10000.7, 3.5e6, T0,1e-2,1e-8,8910))
"""
 
#==============================================================================



m=[]
ma= open("matériaux.txt","r",encoding="utf-16")
print(ma)
Ntwork(ma,circuit)
print("Network creation")
print(ma)
ma= open("matériaux.txt","r",encoding="utf-16")
boundaries(ma,BC,A,fqc,Sim,T1,T2,circuit) 
print("Set boudary conditions ")
ma= open("matériaux.txt","r",encoding="utf-16")
simulator = circuit.simulator(temperature=25,nominal_temperature=25)

if Sim == 'TRANSIENT' : 
  if fqc=='dc' : 
     #Define step_time and end_time  
     analysis = simulator.transient(step_time=0.000000001, end_time=0.000001)
  else : 
      #Define step_time and end_time 
      analysis = simulator.transient(step_time=1/(200*fqc), end_time=50+(2/fqc))
elif  Sim =='Operating Point':
    analysis  = simulator.operating_point()
    print("Operating point analysis ongoing ")
#MV,MT= ExtractP(ma,analysis,Sim)
#Extraction of temperature and potential in the case of operating point
if Sim =='Operating Point':  
 print("Potential extraction : Start")    
 MV,MT= ExtractP(ma,analysis,Sim)
 

 print("Potential extraction: Finished")
 vtk_volt(MV)
 vtk_temp(MT)
 print("Volt/temp VtK creation ")
 print("Current extraction : Start")  
 ma= open("matériaux.txt","r",encoding="utf-16")       
 I,IE=ExtractC(ma,analysis,Sim)
 print("Current extraction : Finished") 
 ma= open("matériaux.txt","r",encoding="utf-16")    
 Delta_V,V_moy=Delta_Moy_P(ma,MV,I)
 ma= open("matériaux.txt","r",encoding="utf-16") 
 Delta_T,T_moy=Delta_Moy_P(ma,MT,I)
 #Heat energy 
 IQ= IE-V_moy*I
 IS = np.divide(IQ, T_moy, out=np.zeros_like(IQ, dtype=float), where=T_moy != 0)
 #IS=np.where(T_moy != 0, IQ/T_moy, 0)
 print("Volt/temp VtK creation")
 ma= open("matériaux.txt","r",encoding="utf-16")
 vtk_IN_vc(I,ma)
 ma= open("matériaux.txt","r",encoding="utf-16")
 vtk_IE_vc(IE,ma)
 
 ma= open("matériaux.txt","r",encoding="utf-16")
 vtk_IQ_vc(IQ,ma)
 
 ma= open("matériaux.txt","r",encoding="utf-16")
 vtk_IS_vc(IS,ma)
end = time.time()
elapsed = end - start
print(f'Temps d\'exécution : {elapsed:.2}')

# -*- coding: utf-8 -*-
"""
Created on Fri May 31 12:40:52 2019

@author: Mahdi-Amine MAMOUNI
"""
"""++++""""""+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++""""""++++"""
"""++++""""""+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++""""""++++"""
#############################################Importation des libraires ##########################################
"""++++""""""+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++""""""++++"""
"""++++""""""+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++""""""++++"""
from tqdm import tqdm                                                                         
import time                                                                                   
from math import exp,cos,sin                                                               
import PySpice.Logging.Logging as Logging                                           
logger = Logging.setup_logging()                                                    
import matplotlib.pyplot as plt	                                              
import numpy as np                                                                  
from PySpice.Spice.Netlist import Circuit                                           
from PySpice.Unit import *                                                        
from matplotlib.animation import FuncAnimation                                    
import matplotlib.animation as animation                                          
from decimal import Decimal                                                       
import math      

from scipy import interpolate 
from scipy.integrate import quad 
from math import pi  
import statistics
from scipy import interpolate 
from scipy.signal.signaltools import correlate
from scipy.integrate import quad
import scipy.fftpack
import numpy as np 
import matplotlib.pyplot as plt 
from math import pi 
from numpy import *
from scipy.signal import *
from scipy import signal 
from PySpice.Spice.Netlist import SubCircuitFactory,SubCircuit

                                                                 
"""++++""""""+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++""""""++++"""
"""++++""""""+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++""""""++++"""
#Class of type of material for eatch dipole 

class OnsaNod(SubCircuit):
    # inE entrée circuit électrique
    # outE sortie circuit électrique
    # inT température haute
    # outT température basse
    NODES = ('inE', 'outE', 'inT', 'ouT')
    ##############################################
    def __init__(self ,name,S=200e-6, Sigma= 0.1e-3, Kappa=0.1, Cp= 7.4e-4,T0=300,L=0.01,A=10e-8,rho=8400 ):
        SubCircuit.__init__(self,name,*self.NODES)
        # passage des paramètres
        self._alpha = S
        self._re = 1/(Sigma*A/L)
        self._rt =1/( Kappa*A/L)
        self._ct = Cp*L*A
        self._K = self._alpha / self._re
        self._T0=T0
        # Partie électrique avec contre réaction
        #       |--BE--|
        #     -----RE-----
        REle=(self.R('E', 'inE', 'outE', self._re))
        REle.minus.add_current_probe(self)
        #Vc=self.VCCS('E', 'inE', 'outE', 'inT', 'ouT', self._alpha/self._re)
        Vc=self.B('E', 'inE', 'outE',i=' (-V(ouT)+V(inT))*{:5.8f}'.format(self._K))
        Vc.minus.add_current_probe(self)
        # Partie thermique avec contre réaction
        #      |---BT---|
        #    ----R-|-R----
        #          C
        #          |
        #        E(T0)
        #          |
        #         gnd
        RTa=self.R('Ta', 'inT', 'a', self._rt/2)
        RTb=self.R('Tb', 'a', 'ouT', self._rt/2)
        RTa.minus.add_current_probe(self)
        RTb.minus.add_current_probe(self)
        capacity=self.C('T', 'a', 'b', self._ct)
        capacity.minus.add_current_probe(self)
        self.V('T', 'b', 0, self._T0)
 #       self.V('T', 'b', 0, 300)
        self.my_value = ('({:5.8f}*((V(inT)+V(ouT))/2)+(V(inE)+V(outE))/{:5.8f})*(-V(outE)+V(inE))+({:5.8f}*((V(inE)+V(outE))/2)+ {:5.8f}*{:5.8f}*(V(inT)+V(ouT))/2)*(-V(ouT)+V(inT)))'.format(self._K, 2*self._re, self._K, self._alpha, self._K))
        B=self.B('T', 'inT', 'ouT', i=self.my_value)
        B.minus.add_current_probe(self)        
 


def Ntwork(mat,circuit):
 
 Mat=[]
 tup=()
 for line in mat:    
     line.strip('\n') 
     r=line.split() 
     Mat+=r  
    
 l=len(Mat) # nombre de ligne 
 c=len(Mat[1]) # nomre de colonne
 m=0
 k=0
 q=1 
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$    
 if Mat[0][-1]=='0': # cela voufrai dire que la preière ligne est composé de éléments horizontaux 

      for r in range(0,(l+1)//2):
              
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""          
          if r%2 !=1 :#""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""              
              #========================================================================= 
              if k < (l) :  
               for i in range(1+r*c,c*(1+r)) : 
                  X=(circuit.X(str(q),Mat[k][m], str(i)+'e',str(i+1)+'e',str(i)+'t',str(i+1)+'t'),)
                  tup+=X
                  q+=1
                  m+=1
               k+=1
               #checked !
              else :
                  break
              #========================================================================= 
              
              
              
               #========================================================================= 
              if k < (l) :
                  
                  for i, j in zip ( range (c*(1+r),r*c,-1) , range (c+1+c*r, 2*c+c*r+2 )):
                   
                   X=(circuit.X(str(q),Mat[k][m], str(i)+'e',str(j)+'e',str(i)+'t',str(j)+'t'),)
                   tup+=X
                   q+=1
                   m-=1
                  k+=1
                  m+=1
              else :
                    break#(col+l*col,1+l*col,-1)
               #=========================================================================   
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""                    
          else :#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""              
              #=========================================================================            
              if k < (l) :
                  for j  in  range(c*(1+r),1+r*c,-1):
                      X=(circuit.X(str(q),Mat[k][m],str(j)+'e',str(j-1)+'e',str(j)+'t',str(j-1)+'t'),)
                      tup+=X
                      q+=1
                      m+=1
                  k+=1
                  
              else :
                  break
              #=========================================================================   
              
              #=========================================================================   
              if k < (l-1) :
         
                for i, j in zip ( range (1+r*c,c+r*c+1) , range ( 2*c+c*r,c+c*r, -1 )):
                    X=(circuit.X(str(q),str(Mat[k][m]), str(i)+'e',str(j)+'e',str(i)+'t',str(j)+'t'),)
                    tup+=X
                    q+=1
                    m-=1
                k+=1
                m+=1
              else :
                  break
              #=========================================================================          
        
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$      
 else : # cela voufrai dire que la première ligne est composé de éléments verticaux

      
      for r in range(0,(l+1)//2):
          
          
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""      
          if r%2 !=1 :#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""          
         
              #========================================================================= 
              if k < (l) :  
               for i, j in zip (range (1+r*c,c*(1+r)+1,+1) , range ((c)*(2+r),c+c*r,-1  )): 
                  
                  X=(circuit.X(str(q),str(Mat[k][m]), str(i)+'e',str(j)+'e',str(i)+'t',str(j)+'t'),)
                  tup+=X
                  q+=1
                  m+=1
               k+=1
               m-=2
              else : 
                  break 
               #========================================================================= 
 

               
               #========================================================================= 
                  
              if k < (l) :
                  for j  in  range((1+r)*c+1,(c)*(2+r),1 ): 
                      
                      X=(circuit.X(str(q),str(Mat[k][m]),str(j)+'e',str(j+1)+'e',str(j)+'t',str(j+1)+'t'),)
                      tup+=X
                      q+=1
                      m-=1
                  k+=1
              
              else :
                  break
                  
               #=========================================================================     

#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""            
          else : #"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""              
              
  
               #========================================================================= 
 
              if k < (l) :
                    m+=1
                    for i, j in zip ( range ((c)*(1+r),c*r-1,-1  ) , range ((1+r)*c+1,c*(r+2)+1,+1) ):
                        
                        X=(circuit.X(str(q),str(Mat[k][m]), str(i)+'e',str(j)+'e',str(i)+'t',str(j)+'t'),)
                        tup+=X
                        q+=1
                        m+=1
                    k+=1
                    m-=2
              else :
                      break
                #=========================================================================            
              
              #========================================================================= 
              if k < (l) :
                  for j  in  range ((2+r)*c,(c)*(1+r)+1,-1 ):
                      
                      X=(circuit.X(str(q),str(Mat[k][m]),str(j)+'e',str(j-1)+'e',str(j)+'t',str(j-1)+'t'),)
                      tup+=X
                      q+=1
                      m-=1
                  k+=1
                  m+=1
                  
 #for i in tup: 
     #i.minus.add_current_probe(circuit)


def boundaries(mat,BC,A,fqc,Sim,T1,T2,circuit): 
    Mat=[]
    source=()
    
    for line in mat:    
        line.strip('\n')
        r=line.split() 
        Mat+=r   
    l=len(Mat) # nombre de ligne 
    #print(l)
    c=len(Mat[1]) # nomre de colonne
    
    #m=0
    #k=0
    #q=1
    
    
    
      #bottom boundary 2  

    if Mat[0][-1]=='0':
        s=(l-1)//2
        #print(s)
        if s%2 ==1 :
                
                for  ss in range (c*(s+1),c*(s),-1) : 
                    
                    source+=(circuit.V(str(ss)+'V2',str(ss)+'e', 0,0 ),)
                    source+=(circuit.V(str(ss)+'T2', str(ss)+'t',0 ,T2),)
                    
        else: 
               for ss in range (c*s+1,c*(s+1)+1,1): 
                   
                   source+=(circuit.V(str(ss)+'V2',str(ss)+'e', 0,0 ),)
                   source+=(circuit.V(str(ss)+'T2', str(ss)+'t',0 ,T2),)
                      
    else :
        s=(l+1)//2
        if s%2 ==1 :
                
                for  ss in range (c*(s+1),c*(s),-1) : 
                    
                    source+=(circuit.V(str(ss)+'V2',str(ss)+'e', 0,0 ),)
                    source+=(circuit.V(str(ss)+'T2', str(ss)+'t',0 ,T2),)
                    
        else: 
               for ss in range (c*(s)+1,c*(s+1)+1):
                   
                   source+=(circuit.V(str(ss)+'V2',str(ss)+'e', 0,0 ),)
                   source+=(circuit.V(str(ss)+'T2', str(ss)+'t',0 ,T2),)
        
        
    #Top boundary +simulation condition     
        
    if Sim == 'TRANSIENT' : 
        if fqc!='dc':    
            if BC == 'Current': 
                for kk in range (1,c+1):
                   source+=(circuit.SinusoidalCurrentSource(str(kk)+'I1',str(kk)+'e',0 , amplitude=A, frequency=fqc),)
                   source+=(circuit.V(str(kk)+'T1', str(kk)+'t',0 ,T1),)
           
            elif BC == 'Voltage':    
                for kk in range (1,c+1):
                   source+=(circuit.SinusoidalVoltageSource(str(kk)+'V1',str(kk)+'e',0 , amplitude=A, frequency=fqc),)
                   source+=(circuit.V(str(kk)+'T1', str(kk)+'t',0 ,T1),)
            for i in source :
                 i.minus.add_current_probe(circuit)
            
        else :         
            if BC == 'Current': 
                for kk in range (1,c+1):
                   source+=(circuit.I(str(kk)+'I1',str(kk)+'e',0 , A),)
                   #source+=(circuit.V(str(kk)+'T1', str(kk)+'t',0 ,T1),)
                   source+=(circuit.PulseVoltageSource(str(kk)+'T1', str(kk)+'t',0, initial_value=300, pulsed_value=T1, pulse_width=40, period=40),)
                
            elif BC == 'Voltage':    
                 
                for kk in range (1,c+1):
                   source+=(circuit.PulseVoltageSource(str(kk)+'T1', str(kk)+'t',0, initial_value=300, pulsed_value=T1, pulse_width=40, period=40),) 
                   source+=(circuit.V(str(kk)+'V1',str(kk)+'e',0 , A),)
                   #source+=(circuit.V(str(kk)+'T1', str(kk)+'t',0 ,T1),)
                   
            for i in source :
                 i.minus.add_current_probe(circuit)
            
        
    elif  Sim =='Operating Point': 
            
            if BC == 'Current': 
                for kk in range (1,c+1):
                   source+=(circuit.I(str(kk)+'I1',str(kk)+'e',0 , A),)
                   source+=(circuit.V(str(kk)+'T1', str(kk)+'t',0 ,T1),)
                
            elif BC == 'Voltage':    
                 
                for kk in range (1,c+1):
                   source+=(circuit.V(str(kk)+'V1',str(kk)+'e',0 , A),)
                   source+=(circuit.V(str(kk)+'T1', str(kk)+'t',0 ,T1),)
                   
            for i in source :
                 i.minus.add_current_probe(circuit)




def Delta_Moy_P(mat, MatVolt, I):
  #  Rappelle M[i,j] tel que i est lindexe des lignes et j est lindexe des colonnes 
    x = I.shape[1]  # Nombre de colonnes de I
    y = I.shape[0]  # Nombre de lignes de I
    xx = MatVolt.shape[1]  # Colonnes de MatVolt
    yy = MatVolt.shape[0]  # Lignes de MatVolt
    
    moy_mat = np.zeros((y, x))  
    delta_mat= np.zeros((y, x))
    Mat = []
    for line in mat:
        line.strip('\n')
        r = line.split()
        Mat += r

    if Mat[0][-1] == '0':
        
        j = 0
        m = 0

        while j < y:
            n = 0
            for i in range(x):
                if n+1 < xx and m < yy:  # Vérifier avant d'accéder à MatVolt[n+1, m]
                    moy_mat[j, i] = (MatVolt[m, n] + MatVolt[m, n+1]) / 2
                    delta_mat[j, i] = (MatVolt[m, n] + MatVolt[m, n+1]) / 2
                    
                    
                else : break
                n += 1  # Correction
                if n >= xx or i >= x:
                    break
            
            j += 2
            m += 1
            if j >= y or m > yy:
                break

        
        j = 1
        m = 0

        while j < y:
            n = 0
            for i in range(x):
                if j+1 < y or m+1 < xx:  # Vérification pour éviter l'erreur
                    moy_mat[j, i] = (MatVolt[m, n] + MatVolt[m+1, n]) / 2
                    delta_mat[j, i] = (MatVolt[m, n] + MatVolt[m+1, n]) / 2
                    
                else:
                    break

                n += 1
                if n >= xx:
                    break
             
            j += 2
            m += 1
            if j > y or m >= yy:
                break

    return delta_mat, moy_mat


                   
#This function extract the temperature and the potential from eatch node             
def ExtractP(mat,analys,simul):
   #nodict={} 
   Mat=[]
   
   for line in mat:    
       line.strip('\n')
       r=line.split() 
       Mat+=r   
   lin=len(Mat) # nombre de ligne 
   #print(l)
   col=len(Mat[1]) # nomre de colonne
   #print(Mat,lin,col )
   fn=[]
   sn=[]
   for i in analys.nodes.values() : 
       if simul == 'Operating Point': 
        fn.append(float(i))        
        sn.append(str(i))
       else :
        fn.append(np.array(i))
        sn.append(str(i))
   #print(fn, sn)        
   fn.reverse()
   sn.reverse()
   nft=[]
   nst=[]
   nfv=[]
   nsv=[]
   for i in (range(len(fn)//2)): 
       for j,k in zip(sn,fn) : 
         #print(str(i+1)+'t',j)  
         if str(i+1)+'t' == j :        
           nst.append(j)
           nft.append(k)
         elif str(i+1)+'e' == j :        
           nsv.append(j)
           nfv.append(k) 
   #for ii, i, j , k in zip (nsv,nst,nft, nfv): 
   #    print(ii,i,j,k)            
   j =0
   k=0
   if Mat[0][-1]=='0': 
    MT=np.zeros((((lin+2)//2),col))
    MV=np.zeros((((lin+2)//2),col))
   else : 
       MT=np.zeros((((lin+3)//2),col))
       MV=np.zeros((((lin+3)//2),col))
       
   while k<len(nft) :    
       for i in range(col) : 
           MT[j,i] = nfv[k]
           MV[j,i] = nft[k]
           k+=1
       j+=1
       if  k==len(nfv) :
           break 
       for i in range(col-1, -1,-1 ) :   
           MT[j,i] = nfv[k]
           MV[j,i] = nft[k]
           k+=1
       j+=1
       if  k==len(nfv) :
           break
   return MT, MV             
#========================================================================= 
#Moyenne Arithmetique  
#Transform the Temperature and the électric potential to VTK file 
def vtk_temp(mat):
  x=mat.shape[1]
  y=mat.shape[0]
  z=1
  c = open("temp.vtk", "w")	
  c.write('# vtk DataFile Version 3.0\nCube example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(y)+ ' ' + str(z)+'\nPOINTS  '+  str(x*y*z)+'  float')
  c.close()
  c = open("temp.vtk", "a")
  for k in range (z) :
     for j in range (y) :
        for i in range (x): 
             c.write('\n\t'+str(float(i))+'\t'+str(float(-j*x/y))+'\t'+str(float(k)))	
  c.write('\nPOINT_DATA \t' +str(x*y*z))
  c.write(' \n \nSCALARS Temp.(K) float\nLOOKUP_TABLE default')
  for i in range (y): 
      for j in range (x):
          c.write('\n\t'+str(mat[i,j]))
  c.close()
  return 
def vtk_volt(mat):
  x=mat.shape[1]
  y=mat.shape[0]
  z=1
  f = open("volt.vtk", "w")	
  f.write('# vtk DataFile Version 3.0\nCube example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(y)+ ' ' + str(z)+
'\nPOINTS  '+  str(x*y*z)+'  float')
  f.close()
  f = open("volt.vtk", "a")
  for k in range (z) :
     for j in range (y):
        for i in range (x): 
             f.write('\n\t'+str(float(i))+'\t'+str(float(-j*x/y))+'\t'+str(float(k)))	
  f.write('\nPOINT_DATA \t' +str(x*y*z))
  f.write(' \n \nSCALARS Pot.(V) float\nLOOKUP_TABLE default')
  for i in range (y): 
        for j in range (x):
          f.write('\n\t'+str(mat[i,j]))
  f.close()
  return


def ExtractC(mat,analys,simul):
   #J'extrait les caractéristiques de la matrice 
   Mat=[]
   for line in mat:    
       line.strip('\n')
       r=line.split() 
       Mat+=r   
   lin=len(Mat) # nombre de ligne 
   #print(l)
   col=len(Mat[1]) # nomre de colonne
   #print(Mat,lin,col )
   fn=[]
   sn=[]
   for i in analys.branches.values() : 
       if simul == 'Operating Point': 
        fn.append(float(i))        
        sn.append(str(i))
       else :
        fn.append(np.array(i))
        sn.append(str(i))
   #print(sn, fn) 
   #Valeur du courant ohmique passant par la resistance electrique        
   CuOhmValue=[]
   #nom de chaque noeud 
   CuOhmString=[]
   #Valeur du courant fourrier passant par la resistance thermique d'entrée a
   CuFourValue_a=[]
   #nom de chaque noeud
   CuFourString_a=[]
   #valeur du courant de la contreréaction électrique  
   TE_eValue=[]
   #nom e la contre réaction électrique  
   TE_eString=[]
   #valeur du courant de la contreréaction électrique
   TE_tValue=[]
   #nom de la contre réaction thermique 
   TE_tString=[]
   for ii in (range(len(fn)//2)): 
       #print(ii,'voici i')
       for j,k in zip(sn,fn) : 
         #print(str(i+1),j) 
         if 'v.x'+str(ii+1)+'.vre_minus' == j :
           CuOhmString.append(j)
           CuOhmValue.append(k)
         elif 'v.x'+str(ii+1)+'.vbe_minus' == j :
           TE_eString.append(j)
           TE_eValue.append(k)
         elif 'v.x'+str(ii+1)+'.vrta_minus' == j :        
           CuFourString_a.append(j)
           CuFourValue_a.append(k)
         elif 'v.x'+str(ii+1)+'.vbt_minus' == j : 
           TE_tString.append(j)
           TE_tValue.append(k)  
   #print(CuOhmValue,CuOhmString,CuFourValue_a,CuFourString_a, TE_eValue,TE_eString,TE_tValue,TE_tString)            
   #print(CuOhmValue, TE_eValue)
   current_elec=np.array(CuOhmValue)+np.array(TE_eValue)            
   current_energ=np.array(CuFourValue_a)+np.array(TE_tValue)
   j =0
   k=0 
   MEC=np.zeros((lin,col)) #matrix electric current
   MTC=np.zeros((lin,col)) #matrix energy current
   #print(Mat,"begingin of mat")
   if Mat[0][-1] =='0':
       while k<len(current_elec)-1 :   
           for i in range(0,col-1) :
               
               
               MTC[j,i] = current_energ[k]
               MEC[j,i] = current_elec[k]
               k+=1
           j+=1
           if  k==len(CuOhmValue) :
               break 
           for i in range(col-1, -1,-1 ) :   
               #print(k,'hereitis...0')
               MEC[j,i] = current_elec[k]
               MTC[j,i] = current_energ[k]
               k+=1
           j+=1
           if  k==len(current_elec) :
               break
   else : 
       while k<len(current_elec)-1 :
           for i in range(col) : 
               #print(k,'hereitis')
               MEC[j,i] = current_elec[k]
               MTC[j,i] = current_energ[k]
               k+=1
           j+=1
           
           if k==len(current_elec) :
               
               break 
           for i in range(col-2, -1,-1 ) :
               #print(k,'hereitis')
               MEC[j,i] = current_elec[k]
               MTC[j,i] = current_energ[k]
               k+=1
           j+=1
           
           if  (k)==(len(current_elec)) :
               
               break
                       
   return MEC ,MTC          
    


def vtk_IN_vc(I,mat):
  Mat=[]
  source=()
    
  for line in mat:    
        line.strip('\n')
        r=line.split() 
        Mat+=r   
  l=len(Mat) # nombre de ligne 
    #print(l)
  c=len(Mat[1]) # nomre de colonne
       
  x=I.shape[1]
  #print(x)
  y=I.shape[0]
  #print(y)
  fin=I.shape[0]
  #print(I[1,1])
  z=1
  
 #----------------------------------------------------------------------------------------------------- 
 #----------------------------------------------------------------------------------------------------- 
 #Si on Commence par une ligne horizontale 
  if Mat[0][-1] =="0":
    
    #-----------------------------------------------------------------------------------------------------
    #Et termine par une ligne horizontale   
    if fin%2 ==1 :   
      f = open("IN_vc.vtk", "w")	
      f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(int((y+1)/2))+ ' ' + str(z)+
    '\nPOINTS  '+  str(int(x*((y+1)/2)*z))+'  float')  
      f.close()
      
      f = open("IN_vc.vtk", "a")
      j=0.0
      for k in range (z) :
         while j<(y+1)/2:
            for i in range (x): 
              #f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
              f.write('\n\t'+str(float(((i))))+'\t'+str(float((-j)*(2*x/(y+1))))+'\t'+str(float(k))) 
            j=j+1
            if j==(y+1)/2:
               break                	
      f.write('\nPOINT_DATA \t' +str(int(x*((y+1)/2)*z)))
      f.write(' \n \nVECTORS IN float\n')
      i=0 
      for j in range (x):
              f.write('\n\t\t'+str(-I[i,j])+'\t\t0'+'\t\t0')
      i=i+1
      while i< (y):
          for j in range (x):####A revoire !!!
              f.write('\n\t\t'+str(-I[i+1,j])+'\t\t'+str(I[i,j])+'\t\t0')
          i=i+2
          if i==fin:
             break 
            #for j in range (x):
             # f.write('\n\t\t0'+'\t\t'+str(I[i,j])+'\t\t0')
              #f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0') 
            #i=i+1
      f.close()
    #-----------------------------------------------------------------------------------------------------  
    else : #si ici si ça se termine par une ligne verticale et non pas horizontale 
        f = open("IN_vc.vtk", "w")	
        f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(int((y)/2))+ ' ' + str(z)+
      '\nPOINTS  '+  str(int(x*((y)/2)*z))+'  float')  
        f.close()
        f = open("IN_vc.vtk", "a")
        j=0.0
        for k in range (z) :
           while j<(y)/2:
              for i in range (x): 
                f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j*(2*x/y)))+'\t'+str(float(k))) 
                #f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
              j=j+1
              if j==(y)/2:
                 break
              #for i in range (x): 
              #     f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k)))
              #j=j+1                  	
        f.write('\nPOINT_DATA \t' +str(int(x*(y/2)*z)))
        f.write(' \n \nVECTORS IN float\n')
        i=0 
        #for j in range (x):
        #        f.write('\n\t\t'+str(-I[i,j])+'\t\t0'+'\t\t0')
        #i=i+1
        while i< (fin-1):
            for j in range (x):
                f.write('\n\t\t'+str(-I[i,j])+'\t\t'+str(I[i,j+1    ])+'\t\t0')
            i=i+2
            if i==fin-1:
               break   
        f.close()
  #-----------------------------------------------------------------------------------------------------#-----------------------------------------------------------------------------------------------------
  #-----------------------------------------------------------------------------------------------------      
  else : 
      if fin%2 ==0 :   
        f = open("IN_vc.vtk", "w")	
        f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(int((y)/2))+ ' ' + str(z)+
      '\nPOINTS  '+  str(int(x*((y)/2)*z))+'  float')  
        f.close()
        f = open("IN_vc.vtk", "a")
        j=0.0
        for k in range (z) :
           while j<(y/2):
              for i in range (x): 
                #f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
                f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j*(2*x/y)))+'\t'+str(float(k))) 
              j=j+1
              if j==y:
                 break                  	
        f.write('\nPOINT_DATA \t' +str(int(x*(y/2)*z)))
        f.write(' \n \nVECTORS IN float\n')
        i=0 
        while i< (fin):
            for j in range (x):
                f.write('\n\t\t'+str(-I[i+1,j])+'\t\t'+str(I[i,j])+'\t\t0')
            i=i+2
            if i==fin:
               break
        f.close()
      else :
          f = open("IN_vc.vtk", "w")	
          f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(int((y+1)/2))+ ' ' + str(z)+
        '\nPOINTS  '+  str(x*(int((y+1)/2)*z))+'  float')  
          f.close()
          f = open("IN_vc.vtk", "a")
          j=0.0
          for k in range (z) :
             while j<(y/2):
                for i in range (x): 
                  #f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
                  f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j*(2*x/y)))+'\t'+str(float(k))) 
                j=j+1
                if j==y/2:
                   break                 	
          f.write('\nPOINT_DATA \t' +str(x*(int((y+1)/2)*z)))
          f.write(' \n \nVECTORS IN float\n')
          i=0 
          
          while i< (fin-1):
              for j in range (x):
                  f.write('\n\t\t'+str(-I[i+1,j])+'\t\t'+str(I[i,j])+'\t\t0')
              i=i+2
              
              if i==fin-1:
                 break
          for j in range (x):
                     f.write('\n\t\t0'+'\t\t'+str(I[i,j])+'\t\t0')
          f.close()
  f.close()
  return
   
def vtk_IE_vc(Ie,mat):
  Mat=[]
  source=()
    
  for line in mat:    
        line.strip('\n')
        r=line.split() 
        Mat+=r   
  l=len(Mat) # nombre de ligne 
    #print(l)
  c=len(Mat[1]) # nomre de colonne
       
  x=Ie.shape[1]
  #print(x)
  y=Ie.shape[0]
  #print(y)
  fin=Ie.shape[0]
  #print(I[1,1])
  z=1
  
 #----------------------------------------------------------------------------------------------------- 
 #----------------------------------------------------------------------------------------------------- 
 #Si on Commence par une ligne horizontale 
  if Mat[0][-1] =="0":
    
    #-----------------------------------------------------------------------------------------------------
    #Et termine par une ligne horizontale   
    if fin%2 ==1 :   
      f = open("IE_vc.vtk", "w")	
      f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(int((y+1)/2))+ ' ' + str(z)+
    '\nPOINTS  '+  str(int(x*((y+1)/2)*z))+'  float')  
      f.close()
      
      f = open("IE_vc.vtk", "a")
      j=0.0
      for k in range (z) :
         while j<(y+1)/2:
            for i in range (x): 
              #f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
              f.write('\n\t'+str(float(((i))))+'\t'+str(float((-j)*(2*x/(y+1))))+'\t'+str(float(k))) 
            j=j+1
            if j==(y+1)/2:
               break                	
      f.write('\nPOINT_DATA \t' +str(int(x*((y+1)/2)*z)))
      f.write(' \n \nVECTORS IN float\n')
      i=0 
      for j in range (x):
              f.write('\n\t\t'+str(-Ie[i,j])+'\t\t0'+'\t\t0')
      i=i+1
      while i< (y):
          for j in range (x):####A revoire !!!
              f.write('\n\t\t'+str(-Ie[i+1,j])+'\t\t'+str(Ie[i,j])+'\t\t0')
          i=i+2
          if i==fin:
             break 
            #for j in range (x):
             # f.write('\n\t\t0'+'\t\t'+str(I[i,j])+'\t\t0')
              #f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0') 
            #i=i+1
      f.close()
    #-----------------------------------------------------------------------------------------------------  
    else : #si ici si ça se termine par une ligne verticale et non pas horizontale 
        f = open("IE_vc.vtk", "w")	
        f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(int((y)/2))+ ' ' + str(z)+
      '\nPOINTS  '+  str(int(x*((y)/2)*z))+'  float')  
        f.close()
        f = open("IE_vc.vtk", "a")
        j=0.0
        for k in range (z) :
           while j<(y)/2:
              for i in range (x): 
                f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j*(2*x/y)))+'\t'+str(float(k))) 
                #f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
              j=j+1
              if j==(y)/2:
                 break
              #for i in range (x): 
              #     f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k)))
              #j=j+1                  	
        f.write('\nPOINT_DATA \t' +str(int(x*(y/2)*z)))
        f.write(' \n \nVECTORS IN float\n')
        i=0 
        #for j in range (x):
        #        f.write('\n\t\t'+str(-I[i,j])+'\t\t0'+'\t\t0')
        #i=i+1
        while i< (fin-1):
            for j in range (x):
                f.write('\n\t\t'+str(-Ie[i,j])+'\t\t'+str(Ie[i,j+1    ])+'\t\t0')
            i=i+2
            if i==fin-1:
               break   
        f.close()
  #-----------------------------------------------------------------------------------------------------#-----------------------------------------------------------------------------------------------------
  #-----------------------------------------------------------------------------------------------------      
  else : 
      if fin%2 ==0 :   
        f = open("IE_vc.vtk", "w")	
        f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(int((y)/2))+ ' ' + str(z)+
      '\nPOINTS  '+  str(int(x*((y)/2)*z))+'  float')  
        f.close()
        f = open("IN_vc.vtk", "a")
        j=0.0
        for k in range (z) :
           while j<(y/2):
              for i in range (x): 
                #f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
                f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j*(2*x/y)))+'\t'+str(float(k))) 
              j=j+1
              if j==y:
                 break                  	
        f.write('\nPOINT_DATA \t' +str(int(x*(y/2)*z)))
        f.write(' \n \nVECTORS IN float\n')
        i=0 
        while i< (fin):
            for j in range (x):
                f.write('\n\t\t'+str(-Ie[i+1,j])+'\t\t'+str(Ie[i,j])+'\t\t0')
            i=i+2
            if i==fin:
               break
        f.close()
      else :
          f = open("IE_vc.vtk", "w")	
          f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(int((y+1)/2))+ ' ' + str(z)+
        '\nPOINTS  '+  str(x*(int((y+1)/2)*z))+'  float')  
          f.close()
          f = open("IE_vc.vtk", "a")
          j=0.0
          for k in range (z) :
             while j<(y/2):
                for i in range (x): 
                  #f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
                  f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j*(2*x/y)))+'\t'+str(float(k))) 
                j=j+1
                if j==y/2:
                   break                 	
          f.write('\nPOINT_DATA \t' +str(x*(int((y+1)/2)*z)))
          f.write(' \n \nVECTORS IN float\n')
          i=0 
          
          while i< (fin-1):
              for j in range (x):
                  f.write('\n\t\t'+str(-Ie[i+1,j])+'\t\t'+str(Ie[i,j])+'\t\t0')
              i=i+2
              
              if i==fin-1:
                 break
          for j in range (x):
                     f.write('\n\t\t0'+'\t\t'+str(Ie[i,j])+'\t\t0')
          f.close()
  f.close()
  return

def vtk_IQ_vc(I,mat):
  Mat=[]
  source=()
    
  for line in mat:    
        line.strip('\n')
        r=line.split() 
        Mat+=r   
  l=len(Mat) # nombre de ligne 
    #print(l)
  c=len(Mat[1]) # nomre de colonne
       
  x=I.shape[1]
  #print(x)
  y=I.shape[0]
  #print(y)
  fin=I.shape[0]
  #print(I[1,1])
  z=1
  
 #----------------------------------------------------------------------------------------------------- 
 #----------------------------------------------------------------------------------------------------- 
 #Si on Commence par une ligne horizontale 
  if Mat[0][-1] =="0":
    
    #-----------------------------------------------------------------------------------------------------
    #Et termine par une ligne horizontale   
    if fin%2 ==1 :   
      f = open("IQ_vc.vtk", "w")	
      f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(int((y+1)/2))+ ' ' + str(z)+
    '\nPOINTS  '+  str(int(x*((y+1)/2)*z))+'  float')  
      f.close()
      
      f = open("IQ_vc.vtk", "a")
      j=0.0
      for k in range (z) :
         while j<(y+1)/2:
            for i in range (x): 
              #f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
              f.write('\n\t'+str(float(((i))))+'\t'+str(float((-j)*(2*x/(y+1))))+'\t'+str(float(k))) 
            j=j+1
            if j==(y+1)/2:
               break                	
      f.write('\nPOINT_DATA \t' +str(int(x*((y+1)/2)*z)))
      f.write(' \n \nVECTORS IN float\n')
      i=0 
      for j in range (x):
              f.write('\n\t\t'+str(-I[i,j])+'\t\t0'+'\t\t0')
      i=i+1
      while i< (y):
          for j in range (x):####A revoire !!!
              f.write('\n\t\t'+str(-I[i+1,j])+'\t\t'+str(I[i,j])+'\t\t0')
          i=i+2
          if i==fin:
             break 
            #for j in range (x):
             # f.write('\n\t\t0'+'\t\t'+str(I[i,j])+'\t\t0')
              #f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0') 
            #i=i+1
      f.close()
    #-----------------------------------------------------------------------------------------------------  
    else : #si ici si ça se termine par une ligne verticale et non pas horizontale 
        f = open("IQ_vc.vtk", "w")	
        f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(int((y)/2))+ ' ' + str(z)+
      '\nPOINTS  '+  str(int(x*((y)/2)*z))+'  float')  
        f.close()
        f = open("IQ_vc.vtk", "a")
        j=0.0
        for k in range (z) :
           while j<(y)/2:
              for i in range (x): 
                f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j*(2*x/y)))+'\t'+str(float(k))) 
                #f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
              j=j+1
              if j==(y)/2:
                 break
              #for i in range (x): 
              #     f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k)))
              #j=j+1                  	
        f.write('\nPOINT_DATA \t' +str(int(x*(y/2)*z)))
        f.write(' \n \nVECTORS IN float\n')
        i=0 
        #for j in range (x):
        #        f.write('\n\t\t'+str(-I[i,j])+'\t\t0'+'\t\t0')
        #i=i+1
        while i< (fin-1):
            for j in range (x):
                f.write('\n\t\t'+str(-I[i,j])+'\t\t'+str(I[i,j+1    ])+'\t\t0')
            i=i+2
            if i==fin-1:
               break   
        f.close()
  #-----------------------------------------------------------------------------------------------------#-----------------------------------------------------------------------------------------------------
  #-----------------------------------------------------------------------------------------------------      
  else : 
      if fin%2 ==0 :   
        f = open("IQ_vc.vtk", "w")	
        f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(int((y)/2))+ ' ' + str(z)+
      '\nPOINTS  '+  str(int(x*((y)/2)*z))+'  float')  
        f.close()
        f = open("IQ_vc.vtk", "a")
        j=0.0
        for k in range (z) :
           while j<(y/2):
              for i in range (x): 
                #f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
                f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j*(2*x/y)))+'\t'+str(float(k))) 
              j=j+1
              if j==y:
                 break                  	
        f.write('\nPOINT_DATA \t' +str(int(x*(y/2)*z)))
        f.write(' \n \nVECTORS IN float\n')
        i=0 
        while i< (fin):
            for j in range (x):
                f.write('\n\t\t'+str(-I[i+1,j])+'\t\t'+str(I[i,j])+'\t\t0')
            i=i+2
            if i==fin:
               break
        f.close()
      else :
          f = open("IQ_vc.vtk", "w")	
          f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(int((y+1)/2))+ ' ' + str(z)+
        '\nPOINTS  '+  str(x*(int((y+1)/2)*z))+'  float')  
          f.close()
          f = open("IQ_vc.vtk", "a")
          j=0.0
          for k in range (z) :
             while j<(y/2):
                for i in range (x): 
                  #f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
                  f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j*(2*x/y)))+'\t'+str(float(k))) 
                j=j+1
                if j==y/2:
                   break                 	
          f.write('\nPOINT_DATA \t' +str(x*(int((y+1)/2)*z)))
          f.write(' \n \nVECTORS IN float\n')
          i=0 
          
          while i< (fin-1):
              for j in range (x):
                  f.write('\n\t\t'+str(-I[i+1,j])+'\t\t'+str(I[i,j])+'\t\t0')
              i=i+2
              
              if i==fin-1:
                 break
          for j in range (x):
                     f.write('\n\t\t0'+'\t\t'+str(I[i,j])+'\t\t0')
          f.close()
  f.close()
  return


def vtk_IS_vc(I,mat):
  Mat=[]
  source=()
    
  for line in mat:    
        line.strip('\n')
        r=line.split() 
        Mat+=r   
  l=len(Mat) # nombre de ligne 
    #print(l)
  c=len(Mat[1]) # nomre de colonne
       
  x=I.shape[1]
  #print(x)
  y=I.shape[0]
  #print(y)
  fin=I.shape[0]
  #print(I[1,1])
  z=1
  
 #----------------------------------------------------------------------------------------------------- 
 #----------------------------------------------------------------------------------------------------- 
 #Si on Commence par une ligne horizontale 
  if Mat[0][-1] =="0":
    
    #-----------------------------------------------------------------------------------------------------
    #Et termine par une ligne horizontale   
    if fin%2 ==1 :   
      f = open("IS_vc.vtk", "w")	
      f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(int((y+1)/2))+ ' ' + str(z)+
    '\nPOINTS  '+  str(int(x*((y+1)/2)*z))+'  float')  
      f.close()
      
      f = open("IS_vc.vtk", "a")
      j=0.0
      for k in range (z) :
         while j<(y+1)/2:
            for i in range (x): 
              #f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
              f.write('\n\t'+str(float(((i))))+'\t'+str(float((-j)*(2*x/(y+1))))+'\t'+str(float(k))) 
            j=j+1
            if j==(y+1)/2:
               break                	
      f.write('\nPOINT_DATA \t' +str(int(x*((y+1)/2)*z)))
      f.write(' \n \nVECTORS IN float\n')
      i=0 
      for j in range (x):
              f.write('\n\t\t'+str(-I[i,j])+'\t\t0'+'\t\t0')
      i=i+1
      while i< (y):
          for j in range (x):####A revoire !!!
              f.write('\n\t\t'+str(-I[i+1,j])+'\t\t'+str(I[i,j])+'\t\t0')
          i=i+2
          if i==fin:
             break 
            #for j in range (x):
             # f.write('\n\t\t0'+'\t\t'+str(I[i,j])+'\t\t0')
              #f.write('\n\t\t'+str(I[i,j])+'\t\t0'+'\t\t0') 
            #i=i+1
      f.close()
    #-----------------------------------------------------------------------------------------------------  
    else : #si ici si ça se termine par une ligne verticale et non pas horizontale 
        f = open("IS_vc.vtk", "w")	
        f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(int((y)/2))+ ' ' + str(z)+
      '\nPOINTS  '+  str(int(x*((y)/2)*z))+'  float')  
        f.close()
        f = open("IS_vc.vtk", "a")
        j=0.0
        for k in range (z) :
           while j<(y)/2:
              for i in range (x): 
                f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j*(2*x/y)))+'\t'+str(float(k))) 
                #f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
              j=j+1
              if j==(y)/2:
                 break
              #for i in range (x): 
              #     f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k)))
              #j=j+1                  	
        f.write('\nPOINT_DATA \t' +str(int(x*(y/2)*z)))
        f.write(' \n \nVECTORS IN float\n')
        i=0 
        #for j in range (x):
        #        f.write('\n\t\t'+str(-I[i,j])+'\t\t0'+'\t\t0')
        #i=i+1
        while i< (fin-1):
            for j in range (x):
                f.write('\n\t\t'+str(-I[i,j])+'\t\t'+str(I[i,j+1    ])+'\t\t0')
            i=i+2
            if i==fin-1:
               break   
        f.close()
  #-----------------------------------------------------------------------------------------------------#-----------------------------------------------------------------------------------------------------
  #-----------------------------------------------------------------------------------------------------      
  else : 
      if fin%2 ==0 :   
        f = open("IS_vc.vtk", "w")	
        f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(int((y)/2))+ ' ' + str(z)+
      '\nPOINTS  '+  str(int(x*((y)/2)*z))+'  float')  
        f.close()
        f = open("IS_vc.vtk", "a")
        j=0.0
        for k in range (z) :
           while j<(y/2):
              for i in range (x): 
                #f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
                f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j*(2*x/y)))+'\t'+str(float(k))) 
              j=j+1
              if j==y:
                 break                  	
        f.write('\nPOINT_DATA \t' +str(int(x*(y/2)*z)))
        f.write(' \n \nVECTORS IN float\n')
        i=0 
        while i< (fin):
            for j in range (x):
                f.write('\n\t\t'+str(-I[i+1,j])+'\t\t'+str(I[i,j])+'\t\t0')
            i=i+2
            if i==fin:
               break
        f.close()
      else :
          f = open("IS_vc.vtk", "w")	
          f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(int((y+1)/2))+ ' ' + str(z)+
        '\nPOINTS  '+  str(x*(int((y+1)/2)*z))+'  float')  
          f.close()
          f = open("IS_vc.vtk", "a")
          j=0.0
          for k in range (z) :
             while j<(y/2):
                for i in range (x): 
                  #f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
                  f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j*(2*x/y)))+'\t'+str(float(k))) 
                j=j+1
                if j==y/2:
                   break                 	
          f.write('\nPOINT_DATA \t' +str(x*(int((y+1)/2)*z)))
          f.write(' \n \nVECTORS IN float\n')
          i=0 
          
          while i< (fin-1):
              for j in range (x):
                  f.write('\n\t\t'+str(-I[i+1,j])+'\t\t'+str(I[i,j])+'\t\t0')
              i=i+2
              
              if i==fin-1:
                 break
          for j in range (x):
                     f.write('\n\t\t0'+'\t\t'+str(I[i,j])+'\t\t0')
          f.close()
  f.close()
  return


"""revoir et corriger cette fonction
def Magnefile(I,length):
       
    x=I.shape[1]  #nombre de colonnes au vecteur I  
    y=I.shape[0]  # nombre de lignes au vecteur Y 
    fin=I.shape[0] #
    z=0*length #la hauteur a laquelle on va faire le calcule 
    k=z
   #----------------------------------------------------------------------------------------------------- 
   #----------------------------------------------------------------------------------------------------- 
   #Commence par une ligne horizontale 
    if I[0][-1] ==0.0:
      #Et termine par une ligne horizontale   
      #-----------------------------------------------------------------------------------------------------
      if fin%2 ==1 : # il se termine par une ligne horizontale 
        # La première ligne va etre considéré comme 
        f = open("Magnetfile.txt", "w")
        #Les intitulées de mes elements   
        f.write('x\t y\t z\t  Ix\t Iy\t Iz ')
        #courdonnées des premiers elements + les valeurs pour j==0  
        for i in range (x) :
         f.write('\n' + str(float(i*length))+'\t0\t 0 \t'+str(float(I[0][0]))+'\t0'+'\t'+'0')
        f.close()
        #ensuite je referme le fichier 
        #j'ouvre le fichier encre une fois et je mets les coordonées et valeurs des restes des elements   
        f = open("Magnetfile.txt", "a")
        j=1        
        while j<(y)/2:
              for i in range (x): 
                f.write('\n'+str(float(i*length))+'\t'+str(float((-j)*(length)))+'\t'+str(float(k))+'\t'+str(float(I[i][j+1]))+'\t'+str(float(I[i][j]))+'\t'+'0') 
              j+=1
              if j==(y)/2:
                 break                	
        f.close()
      #-----------------------------------------------------------------------------------------------------  
      else : #ici si ça se termine par une ligne verticale et non pas horizontale 
          # La première ligne va etre considéré comme 
          f = open("Magnetfile.txt", "w")
          #Les intitulées de mes elements   
          f.write('x\t y\t z\t  Ix\t Iy\t Iz  ')
          #courdonnées des premiers elements + les valeurs pour j==0  
          for i in range (x) :
           f.write('\n' + str(float(i*length))+'\t0\t 0 \t'+str(float(I[0][0]))+'\t0'+'\t'+'0')
          f.close()
          #ensuite je referme le fichier 
          #j'ouvre le fichier encre une fois et je mets les coordonées et valeurs des restes des elements   
          f = open("Magnetfile.txt", "a")
          j=1        
          while j<(((y)/2)):
                for i in range (x): 
                  f.write('\n'+str(float(i*length))+'\t'+str(float((-j)*(length)))+'\t'+str(float(k))+'\t'+str(float(I[i][j+1]))+'\t'+str(float(I[i][j]))+'\t'+'0') 
                j+=1
                if j==(((y)/2)):
                   break
          #j+=1     
          for i in range(x) : 
             f.write('\n'+str(float(i*length))+'\t'+str(float((-j)*(length)))+'\t'+str(float(k))+'\t'+str(float(0))+'\t'+str(float(I[i][j]))+'\t'+'0')   
          f.close()
    #-----------------------------------------------------------------------------------------------------#-----------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------      
    else : #commence par une ligne verticale
        if fin%2 ==0 :  #se termine par une ligne horizontale ?  
          # La première ligne va etre considéré comme 
          f = open("Magnetfile.txt", "w")
          #Les intitulées de mes elements   
          f.write('x\t y\t z\t  Ix\t Iy\t Iz ')
          f.close()
          #j'ouvre le fichier encre une fois et je mets les coordonées et valeurs des restes des elements   
          f = open("Magnetfile.txt", "a")
          j=0    
          print(y)
          while j<(y)/2:
                for i in range (x): 
                  f.write('\n'+str(float(i*length))+'\t'+str(float((-j)*(length)))+'\t'+str(float(k))+'\t'+str(float(I[i][j+1]))+'\t'+str(float(I[i][j]))+'\t'+'0') 
                j+=1
                if j==(y)/2:
                   break                	
          f.close()
   #-------------------------------------------------------------------------------------------------
        else :#se termine par une ligne verticale 
            f = open("IN_vc.vtk", "w")	
            f.write('# vtk DataFile Version 3.0\nVector example\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS ' +str(x)+' ' + str(int((y+1)/2))+ ' ' + str(z)+
          '\nPOINTS  '+  str(x*(int((y+1)/2)*z))+'  float')  
            f.close()
            f = open("IN_vc.vtk", "a")
            j=0.0
            for k in range (z) :
               while j<(y/2):
                  for i in range (x): 
                    #f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j)*(x-1)/y)+'\t'+str(float(k))) 
                    f.write('\n\t'+str(float(((i))))+'\t'+str(float(-j*(2*x/y)))+'\t'+str(float(k))) 
                  j=j+1
                  if j==y/2:
                     break                 	
            f.write('\nPOINT_DATA \t' +str(x*(int((y+1)/2)*z)))
            f.write(' \n \nVECTORS IN float\n')
            i=0 
            
            while i< (fin-1):
                for j in range (x):
                    f.write('\n\t\t'+str(-I[i+1,j])+'\t\t'+str(I[i,j])+'\t\t0')
                i=i+2
                
                if i==fin-1:
                   break
            for j in range (x):
                       f.write('\n\t\t0'+'\t\t'+str(I[i,j])+'\t\t0')
            f.close()
    f.close()
"""


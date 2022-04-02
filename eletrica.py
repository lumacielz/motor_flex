#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 16 12:57:13 2022

@author: luiza.maciel
"""


import numpy as np 
import sympy as sp
import scipy.integrate as integrate
from scipy.misc import derivative
from math import radians,exp

A=np.array([[complex(10,-40),complex(20,10),-10],[-10,0,20],[1,-1,0]])
b=np.array([[0],[5],[complex(9.3969,3.4202)]])
#
x=np.linalg.solve(A,b)
#
#teta0=-164
#teta_comb=-4
#delta_teta=37
#tetaf=146
#
#tk=np.linspace(teta0,tetaf,500)
#def find_nearest(array, value):
#    idx = (np.abs(array - value)).argmin()
#    if array[idx]<=value:
#        return array[idx]
#    else:
#        return array[idx-1]
#    
#find_nearest(tk,teta_comb)
#a,m=2,1
#teta=sp.symbols('teta')
#w=1-sp.exp(-a*((teta_comb+delta_teta-teta_comb)/teta)**(m+1))



#teta,R,D,r,L=sp.symbols("t,R,D,r,L")
#s=R*sp.cos(teta)+sp.sqrt(L**2-R**2*(sp.sin(teta))**2)  #deslocamento
#Ap=np.pi*D*(D/2+L+R-s+(2*R/(r-1)))  #area das paredes
#Vo=(np.pi*D**2/4)*(L+R-s+2*R/(r-1)) #volume 
#a0,a1,a2,a3,a3,a4,a5=sp.symbols("a0,a1,a2,a3,a3,a4,a5")
#Ru,AC,concentracaoCO2,concentracaoH2O,concentracaoN2=sp.symbols("Ru,AC,concentracaoCO2,concentracaoH2O,concentracaoN2")
#Ts=sp.symbols('Ts')
#Cp=a0+a1*sp.log(Ts)+a2*sp.log(Ts)**2+a3*sp.log(Ts)**3+a4*sp.log(Ts)**4+a5*sp.log(Ts)**5
#CpO=10228.3426-7184.92333*sp.log(Ts)+2010.86808*sp.log(Ts)**2-279.69496*sp.log(Ts)**3+19.34823*sp.log(Ts)**4-0.53257*sp.log(Ts)**5
#CpN=0.54497*sp.log(Ts)**5-18.69984*sp.log(Ts)**4+254.29554*sp.log(Ts)**3-1712.17390*sp.log(Ts)**2+5708.38047*sp.log(Ts)-7513.3642
#CpCO2=0.20754*sp.log(Ts)**5-6.43522*sp.log(Ts)**4+77.54809*sp.log(Ts)**3-452.81197*sp.log(Ts)**2+1288.4677*sp.log(Ts)-1412.36785
#CpH2O=0.64541*sp.log(Ts)**5-23.54277*sp.log(Ts)**4+339.33662*sp.log(Ts)**3-2414.77575*sp.log(Ts)**2+8490.5218*sp.log(Ts)-11780.76495
#KReag=((Cp/(Cp-Ru))+AC*(CpO/(CpO-Ru))+3.76*AC*(CpN/(CpN-Ru)))/(1+AC+3.76*AC)
#KProd=((CpCO2/(CpCO2-Ru))*concentracaoCO2+(CpH2O/(CpH2O-Ru))*concentracaoH2O+(CpN/(CpN-Ru))*AC*3.76)/(concentracaoCO2+concentracaoH2O+concentracaoN2)
#a,m=sp.symbols("a,m")
#teta_comb,delta_teta=radians(-4),radians(37.6)
#a,m=2,1
#X_teta=1-sp.exp(-a*((teta-teta_comb)/delta_teta)**(m+1))
#dXdt=sp.diff(X_teta,teta)
#def dx(teta):
#    #return a*((teta - teta_comb)/delta_teta)**(m + 1)*(m + 1)*exp(-a*((teta - teta_comb)/delta_teta)**(m + 1))/(teta - teta_comb)
#    return -exp(-a*((teta-teta_comb)/delta_teta)**(m+1))*(-a*(m+1)*((teta-teta_comb)/delta_teta)**m)*1/delta_teta
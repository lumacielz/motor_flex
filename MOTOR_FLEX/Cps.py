#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 20:23:28 2022

@author: luiza.maciel
"""


import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from math import log
import numpy as np

oc = 0.493  
et = 0.507
h2o=0
def Cps(Ts):
    CpH2O = 0.14904476*log(Ts)**5 -5.60452380*log(Ts)**4 + 81.54220310*log(Ts)**3 -573.43988184*log(Ts)**2 + 1954.64102748*log(Ts) -2558.72544088
    CpC2H5OH = 1.40634*log(Ts)**5 -45.60008*log(Ts)**4 + 578.32500*log(Ts)**3 -3579.01508*log(Ts)**2+ 10834.12151 *log(Ts)-12834.47773
    CpC8H18 = 4.9695987*log(Ts)**5 -157.61132*log(Ts)**4 +1966.430702 *log(Ts)**3 -12036.09466*log(Ts)**2 +36241.19904*log(Ts)-43029.69896

    return CpC8H18*oc+et*CpC2H5OH+h2o*CpH2O

def CpRef(Ts):
    a0=-8102.9675849653
    a1=7963.2048889506
    a2=-2923.24238668569
    a3=509.144026930886
    a4=-42.2764373650608
    a5=1.35175803975684
    Cpgasolina =a0+a1*log(Ts)+a2*log(Ts)**2+a3*log(Ts)**3+a4*log(Ts)**4+a5*log(Ts)**5
    return Cpgasolina

c1 = [(Cps(t)/(Cps(t)-8.314)) for t in range(300,3000)]
c2 = [(CpRef(t)/(CpRef(t)-8.314)) for t in range(300,3000)]

#plt.plot(range(300,3000),c1)
plt.plot(range(300,3000),c2)

def fitLog(Ts, a0,a1,a2,a3,a4,a5):
    return a0+a1*np.log(Ts)+a2*np.log(Ts)**2+a3*np.log(Ts)**3+a4*np.log(Ts)**4+a5*np.log(Ts)**5

td = np.linspace(300,3000)
c2 = [(CpRef(t)/(CpRef(t)-8.314)) for t in td]
a0,a1,a2,a3,a4,a5=curve_fit(fitLog,td, c2)[0]

#5 a0,a1,a2,a3,a4,a5= res[0]
#def cpF(Ts):
#    Cpgasolina =a0+a1*log(Ts)+a2*log(Ts)**2+a3*log(Ts)**3+a4*log(Ts)**4+a5*log(Ts)**5
#    return Cpgasolina
#c3 = [(cpF(t)/(cpF(t)-8.314)) for t in td]
#
#plt.plot(td,fitLog(td,*res[0]))
from sympy import log,symbols,diff
a0,a1,a2,a3,a4,a5,Ts,OC,Ru=symbols("a0,a1,a2,a3,a4,a5,Ts,OC,Ru")
Cp=a0+a1*log(Ts)+a2*log(Ts)**2+a3*log(Ts)**3+a4*log(Ts)**4+a5*log(Ts)**5
CpO=0.06180721*log(Ts)**5-1.91808290*log(Ts)**4+23.18329753*log(Ts)**3-135.04796311*log(Ts)**2+377.22447890*log(Ts)-373.49246047
CpN=0.21448719*log(Ts)**5 -7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 -645.21170218*log(Ts)**2 +2087.25937009*log(Ts) -2624.96917520
Cpt =Cp/(1+OC+OC*3.76) + OC/(1+OC+OC*3.76) * CpO + CpN*3.76*OC/(1+OC+OC*3.76)
KReag= ((Cp/(Cp-Ru))+OC*(CpO/(CpO-Ru))+OC*3.76*(CpN/(CpN-Ru)))/(1+OC+OC*3.76)
#KReag = Cpt/(Cpt-Ru)
print(diff(KReag,Ts))
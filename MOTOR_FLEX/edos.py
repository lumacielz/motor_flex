#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 11:15:36 2022

@author: luiza.maciel
"""
from math import radians,degrees,sin,cos,sqrt,log,exp
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp,odeint
from constants import *
from equilibrium import *

#dados lanzafame
def KReag(Ts):
    Cp=a0+a1*log(Ts)+a2*log(Ts)**2+a3*log(Ts)**3+a4*log(Ts)**4+a5*log(Ts)**5
    CpO=0.06180721*log(Ts)**5-1.91808290*log(Ts)**4+23.18329753*log(Ts)**3-135.04796311*log(Ts)**2+377.22447890*log(Ts)-373.49246047
    CpN=0.21448719*log(Ts)**5 -7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 -645.21170218*log(Ts)**2 +2087.25937009*log(Ts) -2624.96917520
    Cpt =Cp/(1+OC+3.772*OC) + OC/(1+OC+3.773*OC) * CpO + CpN*OC*3.772/(1+OC+3.773*OC)
    KReag= ((Cp/(Cp-Ru))+OC*(CpO/(CpO-Ru))+stoichometricN2*(CpN/(CpN-Ru)))/(1+OC+stoichometricN2)
    #KReag = Cpt/(Cpt-Ru)
    return KReag
    
def KProd(Ts):
    CpN=0.21448719*log(Ts)**5 -7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 -645.21170218*log(Ts)**2 +2087.25937009*log(Ts) -2624.96917520
    CpCO2=0.09995132*log(Ts)**5-2.84698619*log(Ts)**4+29.97175824*log(Ts)**3-139.36548820*log(Ts)**2+262.46393332*log(Ts)-77.60476812
    CpH2O=0.14904476*log(Ts)**5 -5.60452380*log(Ts)**4 + 81.54220310*log(Ts)**3 -573.43988184*log(Ts)**2 + 1954.64102748*log(Ts) -2558.72544088
    Cpt=(CpCO2*stoichometricCO2+CpH2O*stoichometricH2O+CpN*stoichometricN2)/(stoichometricCO2+stoichometricH2O+stoichometricN2)
    KProd=((CpCO2/(CpCO2-Ru))*stoichometricCO2+(CpH2O/(CpH2O-Ru))*stoichometricH2O+(CpN/(CpN-Ru))*stoichometricN2)/(stoichometricCO2+stoichometricH2O+stoichometricN2)
    #KProd = Cpt/(Cpt-Ru)
    return KProd

def dkreag(Ts):
    #usando ks
    dkr= (7.0274904719501e-6*OC*(-0.30903605*log(Ts)**4/Ts + 7.6723316*log(Ts)**3/Ts - 69.54989259*log(Ts)**2/Ts + 270.09592622*log(Ts)/Ts - 377.2244789/Ts)*(0.06180721*log(Ts)**5 - 1.9180829*log(Ts)**4 + 23.18329753*log(Ts)**3 - 135.04796311*log(Ts)**2 + 377.2244789*log(Ts) - 373.49246047)/(-0.00265094143125609*Ru + 0.000163847293739346*log(Ts)**5 - 0.00508472542819384*log(Ts)**4 + 0.0614575639354141*log(Ts)**3 - 0.358004240615044*log(Ts)**2 + log(Ts) - 0.990106637721702)**2 + OC*(0.30903605*log(Ts)**4/Ts - 7.6723316*log(Ts)**3/Ts + 69.54989259*log(Ts)**2/Ts - 270.09592622*log(Ts)/Ts + 377.2244789/Ts)/(-Ru + 0.06180721*log(Ts)**5 - 1.9180829*log(Ts)**4 + 23.18329753*log(Ts)**3 - 135.04796311*log(Ts)**2 + 377.2244789*log(Ts) - 373.49246047) + 1.45128124948838e-7*stoichometricN2*(-1.07243595*log(Ts)**4/Ts + 29.22307272*log(Ts)**3/Ts - 293.86871481*log(Ts)**2/Ts + 1290.42340436*log(Ts)/Ts - 2087.25937009/Ts)*(0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752)/(-0.000380956854445275*Ru + 8.1710365221206e-5*log(Ts)**5 - 0.00278318246515918*log(Ts)**4 + 0.0373171004046311*log(Ts)**3 - 0.245797820513774*log(Ts)**2 + 0.795155764040913*log(Ts) - 1)**2 + stoichometricN2*(1.07243595*log(Ts)**4/Ts - 29.22307272*log(Ts)**3/Ts + 293.86871481*log(Ts)**2/Ts - 1290.42340436*log(Ts)/Ts + 2087.25937009/Ts)/(-Ru + 0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752) + (-a1/Ts - 2*a2*log(Ts)/Ts - 3*a3*log(Ts)**2/Ts - 4*a4*log(Ts)**3/Ts - 5*a5*log(Ts)**4/Ts)*(a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)/(-Ru + a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)**2 + (a1/Ts + 2*a2*log(Ts)/Ts + 3*a3*log(Ts)**2/Ts + 4*a4*log(Ts)**3/Ts + 5*a5*log(Ts)**4/Ts)/(-Ru + a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5))/(OC + stoichometricN2 + 1)
    #usando cpt
    #dkr =0.0702840261838923*(-OC*(0.30903605*log(Ts)**4/Ts - 7.6723316*log(Ts)**3/Ts + 69.54989259*log(Ts)**2/Ts - 270.09592622*log(Ts)/Ts + 377.2244789/Ts)/(4.772*OC + 1) - 3.772*OC*(1.07243595*log(Ts)**4/Ts - 29.22307272*log(Ts)**3/Ts + 293.86871481*log(Ts)**2/Ts - 1290.42340436*log(Ts)/Ts + 2087.25937009/Ts)/(4.772*OC + 1) - (a1/Ts + 2*a2*log(Ts)/Ts + 3*a3*log(Ts)**2/Ts + 4*a4*log(Ts)**3/Ts + 5*a5*log(Ts)**4/Ts)/(4.772*OC + 1))*(OC*(0.06180721*log(Ts)**5 - 1.9180829*log(Ts)**4 + 23.18329753*log(Ts)**3 - 135.04796311*log(Ts)**2 + 377.2244789*log(Ts) - 373.49246047)/(4.772*OC + 1) + 3.772*OC*(0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752)/(4.772*OC + 1) + (a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)/(4.772*OC + 1))/(0.265111346765642*OC*(0.06180721*log(Ts)**5 - 1.9180829*log(Ts)**4 + 23.18329753*log(Ts)**3 - 135.04796311*log(Ts)**2 + 377.2244789*log(Ts) - 373.49246047)/(4.772*OC + 1) + OC*(0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752)/(4.772*OC + 1) - 0.265111346765642*Ru + 0.265111346765642*(a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)/(4.772*OC + 1))**2 + (OC*(0.30903605*log(Ts)**4/Ts - 7.6723316*log(Ts)**3/Ts + 69.54989259*log(Ts)**2/Ts - 270.09592622*log(Ts)/Ts + 377.2244789/Ts)/(4.772*OC + 1) + 3.772*OC*(1.07243595*log(Ts)**4/Ts - 29.22307272*log(Ts)**3/Ts + 293.86871481*log(Ts)**2/Ts - 1290.42340436*log(Ts)/Ts + 2087.25937009/Ts)/(4.772*OC + 1) + (a1/Ts + 2*a2*log(Ts)/Ts + 3*a3*log(Ts)**2/Ts + 4*a4*log(Ts)**3/Ts + 5*a5*log(Ts)**4/Ts)/(4.772*OC + 1))/(OC*(0.06180721*log(Ts)**5 - 1.9180829*log(Ts)**4 + 23.18329753*log(Ts)**3 - 135.04796311*log(Ts)**2 + 377.2244789*log(Ts) - 373.49246047)/(4.772*OC + 1) + 3.772*OC*(0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752)/(4.772*OC + 1) - Ru + (a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)/(4.772*OC + 1))
    return dkr

def dkprod(Ts):
    #usando ks
    dkp= (1.45164604139169e-5*stoichometricCO2*(-0.4997566*log(Ts)**4/Ts + 11.38794476*log(Ts)**3/Ts - 89.91527472*log(Ts)**2/Ts + 278.7309764*log(Ts)/Ts - 262.46393332/Ts)*(0.09995132*log(Ts)**5 - 2.84698619*log(Ts)**4 + 29.97175824*log(Ts)**3 - 139.3654882*log(Ts)**2 + 262.46393332*log(Ts) - 77.60476812)/(-0.00381004729812071*Ru + 0.000380819256709598*log(Ts)**5 - 0.0108471520409965*log(Ts)**4 + 0.114193816502239*log(Ts)**3 - 0.530989101767684*log(Ts)**2 + log(Ts) - 0.29567783709689)**2 + stoichometricCO2*(0.4997566*log(Ts)**4/Ts - 11.38794476*log(Ts)**3/Ts + 89.91527472*log(Ts)**2/Ts - 278.7309764*log(Ts)/Ts + 262.46393332/Ts)/(-Ru + 0.09995132*log(Ts)**5 - 2.84698619*log(Ts)**4 + 29.97175824*log(Ts)**3 - 139.3654882*log(Ts)**2 + 262.46393332*log(Ts) - 77.60476812) + 1.52739943457638e-7*stoichometricH2O*(-0.7452238*log(Ts)**4/Ts + 22.4180952*log(Ts)**3/Ts - 244.6266093*log(Ts)**2/Ts + 1146.87976368*log(Ts)/Ts - 1954.64102748/Ts)*(0.14904476*log(Ts)**5 - 5.6045238*log(Ts)**4 + 81.5422031*log(Ts)**3 - 573.43988184*log(Ts)**2 + 1954.64102748*log(Ts) - 2558.72544088)/(-0.000390819579163632*Ru + 5.82496103797445e-5*log(Ts)**5 - 0.00219035763292856*log(Ts)**4 + 0.0318682894996174*log(Ts)**3 - 0.224111533296352*log(Ts)**2 + 0.763911983775703*log(Ts) - 1)**2 + stoichometricH2O*(0.7452238*log(Ts)**4/Ts - 22.4180952*log(Ts)**3/Ts + 244.6266093*log(Ts)**2/Ts - 1146.87976368*log(Ts)/Ts + 1954.64102748/Ts)/(-Ru + 0.14904476*log(Ts)**5 - 5.6045238*log(Ts)**4 + 81.5422031*log(Ts)**3 - 573.43988184*log(Ts)**2 + 1954.64102748*log(Ts) - 2558.72544088) + 1.45128124948838e-7*stoichometricN2*(-1.07243595*log(Ts)**4/Ts + 29.22307272*log(Ts)**3/Ts - 293.86871481*log(Ts)**2/Ts + 1290.42340436*log(Ts)/Ts - 2087.25937009/Ts)*(0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752)/(-0.000380956854445275*Ru + 8.1710365221206e-5*log(Ts)**5 - 0.00278318246515918*log(Ts)**4 + 0.0373171004046311*log(Ts)**3 - 0.245797820513774*log(Ts)**2 + 0.795155764040913*log(Ts) - 1)**2 + stoichometricN2*(1.07243595*log(Ts)**4/Ts - 29.22307272*log(Ts)**3/Ts + 293.86871481*log(Ts)**2/Ts - 1290.42340436*log(Ts)/Ts + 2087.25937009/Ts)/(-Ru + 0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752))/(stoichometricCO2 + stoichometricH2O + stoichometricN2)
    #usando cpt
    #dkp =(stoichometricCO2*(0.4997566*log(Ts)**4/Ts - 11.38794476*log(Ts)**3/Ts + 89.91527472*log(Ts)**2/Ts - 278.7309764*log(Ts)/Ts + 262.46393332/Ts) + stoichometricH2O*(0.7452238*log(Ts)**4/Ts - 22.4180952*log(Ts)**3/Ts + 244.6266093*log(Ts)**2/Ts - 1146.87976368*log(Ts)/Ts + 1954.64102748/Ts) + stoichometricN2*(1.07243595*log(Ts)**4/Ts - 29.22307272*log(Ts)**3/Ts + 293.86871481*log(Ts)**2/Ts - 1290.42340436*log(Ts)/Ts + 2087.25937009/Ts))/((-Ru + (stoichometricCO2*(0.09995132*log(Ts)**5 - 2.84698619*log(Ts)**4 + 29.97175824*log(Ts)**3 - 139.3654882*log(Ts)**2 + 262.46393332*log(Ts) - 77.60476812) + stoichometricH2O*(0.14904476*log(Ts)**5 - 5.6045238*log(Ts)**4 + 81.5422031*log(Ts)**3 - 573.43988184*log(Ts)**2 + 1954.64102748*log(Ts) - 2558.72544088) + stoichometricN2*(0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752))/(stoichometricCO2 + stoichometricH2O + stoichometricN2))*(stoichometricCO2 + stoichometricH2O + stoichometricN2)) - (stoichometricCO2*(0.4997566*log(Ts)**4/Ts - 11.38794476*log(Ts)**3/Ts + 89.91527472*log(Ts)**2/Ts - 278.7309764*log(Ts)/Ts + 262.46393332/Ts) + stoichometricH2O*(0.7452238*log(Ts)**4/Ts - 22.4180952*log(Ts)**3/Ts + 244.6266093*log(Ts)**2/Ts - 1146.87976368*log(Ts)/Ts + 1954.64102748/Ts) + stoichometricN2*(1.07243595*log(Ts)**4/Ts - 29.22307272*log(Ts)**3/Ts + 293.86871481*log(Ts)**2/Ts - 1290.42340436*log(Ts)/Ts + 2087.25937009/Ts))*(stoichometricCO2*(0.09995132*log(Ts)**5 - 2.84698619*log(Ts)**4 + 29.97175824*log(Ts)**3 - 139.3654882*log(Ts)**2 + 262.46393332*log(Ts) - 77.60476812) + stoichometricH2O*(0.14904476*log(Ts)**5 - 5.6045238*log(Ts)**4 + 81.5422031*log(Ts)**3 - 573.43988184*log(Ts)**2 + 1954.64102748*log(Ts) - 2558.72544088) + stoichometricN2*(0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752))/((-Ru + (stoichometricCO2*(0.09995132*log(Ts)**5 - 2.84698619*log(Ts)**4 + 29.97175824*log(Ts)**3 - 139.3654882*log(Ts)**2 + 262.46393332*log(Ts) - 77.60476812) + stoichometricH2O*(0.14904476*log(Ts)**5 - 5.6045238*log(Ts)**4 + 81.5422031*log(Ts)**3 - 573.43988184*log(Ts)**2 + 1954.64102748*log(Ts) - 2558.72544088) + stoichometricN2*(0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752))/(stoichometricCO2 + stoichometricH2O + stoichometricN2))**2*(stoichometricCO2 + stoichometricH2O + stoichometricN2)**2)
    return dkp

#dados tadeu
#def KReag(Ts):
#    Cp=a0+a1*log(Ts)+a2*log(Ts)**2+a3*log(Ts)**3+a4*log(Ts)**4+a5*log(Ts)**5
#    CpO=10228.3426-7184.92333*log(Ts)+2010.86808*log(Ts)**2-279.69496*log(Ts)**3+19.34823*log(Ts)**4-0.53257*log(Ts)**5
#    CpN=0.54497*log(Ts)**5-18.69984*log(Ts)**4+254.29554*log(Ts)**3-1712.17390*log(Ts)**2+5708.38047*log(Ts)-7513.3642
#    #Cpt =Cp/(1+OC+3.772*OC) + OC/(1+OC+3.772*OC) * CpO + CpN*OC*3.772/(1+OC+3.772*OC)
#    KReag=((Cp/(Cp-Ru))+OC*(CpO/(CpO-Ru))+stoichometricN2*(CpN/(CpN-Ru)))/(1+OC+stoichometricN2)
#    return KReag
#    
#def KProd(Ts):
#    CpN=0.54497*log(Ts)**5-18.69984*log(Ts)**4+254.29554*log(Ts)**3-1712.17390*log(Ts)**2+5708.38047*log(Ts)-7513.3642
#    CpCO2=0.20754*log(Ts)**5-6.43522*log(Ts)**4+77.54809*log(Ts)**3-452.81197*log(Ts)**2+1288.4677*log(Ts)-1412.36785
#    CpH2O=0.64541*log(Ts)**5-23.54277*log(Ts)**4+339.33662*log(Ts)**3-2414.77575*log(Ts)**2+8490.5218*log(Ts)-11780.76495
#    #Cpt=(CpCO2*stoichometricCO2+CpH2O*stoichometricH2O+CpN*stoichometricN2)/(stoichometricCO2+stoichometricH2O+stoichometricN2)
#    KProd=((CpCO2/(CpCO2-Ru))*stoichometricCO2+(CpH2O/(CpH2O-Ru))*stoichometricH2O+(CpN/(CpN-Ru))*stoichometricN2)/(stoichometricCO2+stoichometricH2O+stoichometricN2)
#    return KProd
#
#def dkreag(Ts):
#    dkr= (OC*(-2.66285*log(Ts)**4/Ts + 77.39292*log(Ts)**3/Ts - 839.08488*log(Ts)**2/Ts + 4021.73616*log(Ts)/Ts - 7184.92333/Ts)/(-Ru - 0.53257*log(Ts)**5 + 19.34823*log(Ts)**4 - 279.69496*log(Ts)**3 + 2010.86808*log(Ts)**2 - 7184.92333*log(Ts) + 10228.3426) + 9.55849389871466e-9*OC*(2.66285*log(Ts)**4/Ts - 77.39292*log(Ts)**3/Ts + 839.08488*log(Ts)**2/Ts - 4021.73616*log(Ts)/Ts + 7184.92333/Ts)*(-0.53257*log(Ts)**5 + 19.34823*log(Ts)**4 - 279.69496*log(Ts)**3 + 2010.86808*log(Ts)**2 - 7184.92333*log(Ts) + 10228.3426)/(-9.77675503360632e-5*Ru - 5.20680642824772e-5*log(Ts)**5 + 0.00189162905043873*log(Ts)**4 - 0.0273450910805432*log(Ts)**3 + 0.196597646230583*log(Ts)**2 - 0.70245235332653*log(Ts) + 1)**2 + 1.77145905099399e-8*stoichometricN2*(-2.72485*log(Ts)**4/Ts + 74.79936*log(Ts)**3/Ts - 762.88662*log(Ts)**2/Ts + 3424.3478*log(Ts)/Ts - 5708.38047/Ts)*(0.54497*log(Ts)**5 - 18.69984*log(Ts)**4 + 254.29554*log(Ts)**3 - 1712.1739*log(Ts)**2 + 5708.38047*log(Ts) - 7513.3642)/(-0.000133096170155042*Ru + 7.25334198493932e-5*log(Ts)**5 - 0.00248887708651206*log(Ts)**4 + 0.0338457624615083*log(Ts)**3 - 0.227883788729422*log(Ts)**2 + 0.759763578344838*log(Ts) - 1)**2 + stoichometricN2*(2.72485*log(Ts)**4/Ts - 74.79936*log(Ts)**3/Ts + 762.88662*log(Ts)**2/Ts - 3424.3478*log(Ts)/Ts + 5708.38047/Ts)/(-Ru + 0.54497*log(Ts)**5 - 18.69984*log(Ts)**4 + 254.29554*log(Ts)**3 - 1712.1739*log(Ts)**2 + 5708.38047*log(Ts) - 7513.3642) + (-a1/Ts - 2*a2*log(Ts)/Ts - 3*a3*log(Ts)**2/Ts - 4*a4*log(Ts)**3/Ts - 5*a5*log(Ts)**4/Ts)*(a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)/(-Ru + a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)**2 + (a1/Ts + 2*a2*log(Ts)/Ts + 3*a3*log(Ts)**2/Ts + 4*a4*log(Ts)**3/Ts + 5*a5*log(Ts)**4/Ts)/(-Ru + a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5))/(OC + stoichometricN2 + 1)
#    return dkr
#
#def dkprod(Ts):
#    dkp= (5.01307675179101e-7*stoichometricCO2*(-1.0377*log(Ts)**4/Ts + 25.74088*log(Ts)**3/Ts - 232.64427*log(Ts)**2/Ts + 905.62394*log(Ts)/Ts - 1288.4677/Ts)*(0.20754*log(Ts)**5 - 6.43522*log(Ts)**4 + 77.54809*log(Ts)**3 - 452.81197*log(Ts)**2 + 1288.4677*log(Ts) - 1412.36785)/(-0.000708030843381205*Ru + 0.000146944721235335*log(Ts)**5 - 0.0045563342439436*log(Ts)**4 + 0.0549064395653016*log(Ts)**3 - 0.320604841012205*log(Ts)**2 + 0.912274872300442*log(Ts) - 1)**2 + stoichometricCO2*(1.0377*log(Ts)**4/Ts - 25.74088*log(Ts)**3/Ts + 232.64427*log(Ts)**2/Ts - 905.62394*log(Ts)/Ts + 1288.4677/Ts)/(-Ru + 0.20754*log(Ts)**5 - 6.43522*log(Ts)**4 + 77.54809*log(Ts)**3 - 452.81197*log(Ts)**2 + 1288.4677*log(Ts) - 1412.36785) + 7.20531576341265e-9*stoichometricH2O*(-3.22705*log(Ts)**4/Ts + 94.17108*log(Ts)**3/Ts - 1018.00986*log(Ts)**2/Ts + 4829.5515*log(Ts)/Ts - 8490.5218/Ts)*(0.64541*log(Ts)**5 - 23.54277*log(Ts)**4 + 339.33662*log(Ts)**3 - 2414.77575*log(Ts)**2 + 8490.5218*log(Ts) - 11780.76495)/(-8.48841313992942e-5*Ru + 5.47850672464185e-5*log(Ts)**5 - 0.00199840758218336*log(Ts)**4 + 0.0288042942406724*log(Ts)**3 - 0.204976142062829*log(Ts)**2 + 0.720710568119772*log(Ts) - 1)**2 + stoichometricH2O*(3.22705*log(Ts)**4/Ts - 94.17108*log(Ts)**3/Ts + 1018.00986*log(Ts)**2/Ts - 4829.5515*log(Ts)/Ts + 8490.5218/Ts)/(-Ru + 0.64541*log(Ts)**5 - 23.54277*log(Ts)**4 + 339.33662*log(Ts)**3 - 2414.77575*log(Ts)**2 + 8490.5218*log(Ts) - 11780.76495) + 1.77145905099399e-8*stoichometricN2*(-2.72485*log(Ts)**4/Ts + 74.79936*log(Ts)**3/Ts - 762.88662*log(Ts)**2/Ts + 3424.3478*log(Ts)/Ts - 5708.38047/Ts)*(0.54497*log(Ts)**5 - 18.69984*log(Ts)**4 + 254.29554*log(Ts)**3 - 1712.1739*log(Ts)**2 + 5708.38047*log(Ts) - 7513.3642)/(-0.000133096170155042*Ru + 7.25334198493932e-5*log(Ts)**5 - 0.00248887708651206*log(Ts)**4 + 0.0338457624615083*log(Ts)**3 - 0.227883788729422*log(Ts)**2 + 0.759763578344838*log(Ts) - 1)**2 + stoichometricN2*(2.72485*log(Ts)**4/Ts - 74.79936*log(Ts)**3/Ts + 762.88662*log(Ts)**2/Ts - 3424.3478*log(Ts)/Ts + 5708.38047/Ts)/(-Ru + 0.54497*log(Ts)**5 - 18.69984*log(Ts)**4 + 254.29554*log(Ts)**3 - 1712.1739*log(Ts)**2 + 5708.38047*log(Ts) - 7513.3642))/(stoichometricCO2 + stoichometricH2O + stoichometricN2)
#    return dkp

def K(Ts,composicao):
    
    Ncomb=composicao[0]
    NO2=composicao[1]
    NH2O=composicao[2]
    NCO2=composicao[3]
    NCO=composicao[4]
    NO=composicao[5]
    NN=composicao[6]
    NNO2=composicao[7]
    NNO=composicao[8]
    NOH=composicao[9]
    NH2=composicao[10]
    NH=composicao[11]
    NN2=composicao[12]
    
    total_mols = sum(composicao)
    
    Cp=a0+a1*log(Ts)+a2*log(Ts)**2+a3*log(Ts)**3+a4*log(Ts)**4+a5*log(Ts)**5
    
    CpO2=0.06180721*log(Ts)**5-1.91808290*log(Ts)**4+23.18329753*log(Ts)**3-135.04796311*log(Ts)**2+377.22447890*log(Ts)-373.49246047
    CpO=0.07027490*log(Ts)**5-2.31335083*log(Ts)**4 + 30.19191662*log(Ts)**3 -194.76523964*log(Ts)**2 +618.57117075*log(Ts) -748.19500164
    
    CpN2=0.21448719*log(Ts)**5 -7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 -645.21170218*log(Ts)**2 +2087.25937009*log(Ts) -2624.96917520
    CpN=0.14239267*log(Ts)**5 -4.53379847*log(Ts)**4 + 57.24529880*log(Ts)**3 -358.12331558*log(Ts)**2 + 1109.54271743*log(Ts) -1340.57218387
    
    CpH2=-0.1571988*log(Ts)**5 + 5.07583945*log(Ts)**4 -64.62974617*log(Ts)**3 + 406.63764375*log(Ts)**2 -1265.40995607*log(Ts) + 1586.40614067
    CpH=20.78600000
    
    CpH2O=0.14904476*log(Ts)**5 -5.60452380*log(Ts)**4 + 81.54220310*log(Ts)**3 -573.43988184*log(Ts)**2 + 1954.64102748*log(Ts) -2558.72544088
    CpOH=0.03478011*log(Ts)**5 -1.49678745*log(Ts)**4 + 23.95819609 *log(Ts)**3-180.33447646*log(Ts)**2 + 643.43907081*log(Ts) -844.81094064
    
    CpCO2=0.09995132*log(Ts)**5-2.84698619*log(Ts)**4+29.97175824*log(Ts)**3-139.36548820*log(Ts)**2+262.46393332*log(Ts)-77.60476812
    CpCO=0.22518502*log(Ts)**5 -7.60385762*log(Ts)**4 + 101.08195929*log(Ts)**3 -660.23411744*log(Ts)**2 + 2118.61214474*log(Ts) -2644.11606414
    
    CpNO=0.22148616*log(Ts)**5 -7.36223109*log(Ts)**4 +  96.15006747*log(Ts)**3 -615.33775005*log(Ts)**2 + 1927.66396650*log(Ts) -2333.02836259

    Cps = [Cp, CpO2, CpH2O, CpCO2, CpCO, CpO, CpN, 0, CpNO, CpOH, CpH2, CpH, CpN2]
    Kt = 0
    for n, c in zip(composicao, Cps):
        Kt += (c/(c-Ru))*n/total_mols
       
    dK = 1.43033896777609e-7*NCO*(-1.1259251*log(Ts)**4/Ts + 30.41543048*log(Ts)**3/Ts - 303.24587787*log(Ts)**2/Ts + 1320.46823488*log(Ts)/Ts - 2118.61214474/Ts)*(0.22518502*log(Ts)**5 - 7.60385762*log(Ts)**4 + 101.08195929*log(Ts)**3 - 660.23411744*log(Ts)**2 + 2118.61214474*log(Ts) - 2644.11606414)/(total_mols*(-0.000378198224186218*Ru + 8.5164574677338e-5*log(Ts)**5 - 0.00287576544884884*log(Ts)**4 + 0.0382290175007416*log(Ts)**3 - 0.249699370762963*log(Ts)**2 + 0.801255350880023*log(Ts) - 1)**2) + NCO*(1.1259251*log(Ts)**4/Ts - 30.41543048*log(Ts)**3/Ts + 303.24587787*log(Ts)**2/Ts - 1320.46823488*log(Ts)/Ts + 2118.61214474/Ts)/(total_mols*(-Ru + 0.22518502*log(Ts)**5 - 7.60385762*log(Ts)**4 + 101.08195929*log(Ts)**3 - 660.23411744*log(Ts)**2 + 2118.61214474*log(Ts) - 2644.11606414)) + 1.45164604139169e-5*NCO2*(-0.4997566*log(Ts)**4/Ts + 11.38794476*log(Ts)**3/Ts - 89.91527472*log(Ts)**2/Ts + 278.7309764*log(Ts)/Ts - 262.46393332/Ts)*(0.09995132*log(Ts)**5 - 2.84698619*log(Ts)**4 + 29.97175824*log(Ts)**3 - 139.3654882*log(Ts)**2 + 262.46393332*log(Ts) - 77.60476812)/(total_mols*(-0.00381004729812071*Ru + 0.000380819256709598*log(Ts)**5 - 0.0108471520409965*log(Ts)**4 + 0.114193816502239*log(Ts)**3 - 0.530989101767684*log(Ts)**2 + log(Ts) - 0.29567783709689)**2) + NCO2*(0.4997566*log(Ts)**4/Ts - 11.38794476*log(Ts)**3/Ts + 89.91527472*log(Ts)**2/Ts - 278.7309764*log(Ts)/Ts + 262.46393332/Ts)/(total_mols*(-Ru + 0.09995132*log(Ts)**5 - 2.84698619*log(Ts)**4 + 29.97175824*log(Ts)**3 - 139.3654882*log(Ts)**2 + 262.46393332*log(Ts) - 77.60476812)) + NH2*(-0.785994*log(Ts)**4/Ts + 20.3033578*log(Ts)**3/Ts - 193.88923851*log(Ts)**2/Ts + 813.2752875*log(Ts)/Ts - 1265.40995607/Ts)/(total_mols*(-Ru - 0.1571988*log(Ts)**5 + 5.07583945*log(Ts)**4 - 64.62974617*log(Ts)**3 + 406.63764375*log(Ts)**2 - 1265.40995607*log(Ts) + 1586.40614067)) + 3.97348186706385e-7*NH2*(0.785994*log(Ts)**4/Ts - 20.3033578*log(Ts)**3/Ts + 193.88923851*log(Ts)**2/Ts - 813.2752875*log(Ts)/Ts + 1265.40995607/Ts)*(-0.1571988*log(Ts)**5 + 5.07583945*log(Ts)**4 - 64.62974617*log(Ts)**3 + 406.63764375*log(Ts)**2 - 1265.40995607*log(Ts) + 1586.40614067)/(total_mols*(-0.000630355603375099*Ru - 9.90911444238415e-5*log(Ts)**5 + 0.00319958383913988*log(Ts)**4 - 0.0407397226429698*log(Ts)**3 + 0.25632631728106*log(Ts)**2 - 0.797658256375362*log(Ts) + 1)**2) + 1.52739943457638e-7*NH2O*(-0.7452238*log(Ts)**4/Ts + 22.4180952*log(Ts)**3/Ts - 244.6266093*log(Ts)**2/Ts + 1146.87976368*log(Ts)/Ts - 1954.64102748/Ts)*(0.14904476*log(Ts)**5 - 5.6045238*log(Ts)**4 + 81.5422031*log(Ts)**3 - 573.43988184*log(Ts)**2 + 1954.64102748*log(Ts) - 2558.72544088)/(total_mols*(-0.000390819579163632*Ru + 5.82496103797445e-5*log(Ts)**5 - 0.00219035763292856*log(Ts)**4 + 0.0318682894996174*log(Ts)**3 - 0.224111533296352*log(Ts)**2 + 0.763911983775703*log(Ts) - 1)**2) + NH2O*(0.7452238*log(Ts)**4/Ts - 22.4180952*log(Ts)**3/Ts + 244.6266093*log(Ts)**2/Ts - 1146.87976368*log(Ts)/Ts + 1954.64102748/Ts)/(total_mols*(-Ru + 0.14904476*log(Ts)**5 - 5.6045238*log(Ts)**4 + 81.5422031*log(Ts)**3 - 573.43988184*log(Ts)**2 + 1954.64102748*log(Ts) - 2558.72544088)) + 5.56441602198163e-7*NN*(-0.71196335*log(Ts)**4/Ts + 18.13519388*log(Ts)**3/Ts - 171.7358964*log(Ts)**2/Ts + 716.24663116*log(Ts)/Ts - 1109.54271743/Ts)*(0.14239267*log(Ts)**5 - 4.53379847*log(Ts)**4 + 57.2452988*log(Ts)**3 - 358.12331558*log(Ts)**2 + 1109.54271743*log(Ts) - 1340.57218387)/(total_mols*(-0.000745950133854913*Ru + 0.000106217831246459*log(Ts)**5 - 0.0033819875755677*log(Ts)**4 + 0.0427021383024245*log(Ts)**3 - 0.267142135193466*log(Ts)**2 + 0.827663538584653*log(Ts) - 1)**2) + NN*(0.71196335*log(Ts)**4/Ts - 18.13519388*log(Ts)**3/Ts + 171.7358964*log(Ts)**2/Ts - 716.24663116*log(Ts)/Ts + 1109.54271743/Ts)/(total_mols*(-Ru + 0.14239267*log(Ts)**5 - 4.53379847*log(Ts)**4 + 57.2452988*log(Ts)**3 - 358.12331558*log(Ts)**2 + 1109.54271743*log(Ts) - 1340.57218387)) + 1.45128124948838e-7*NN2*(-1.07243595*log(Ts)**4/Ts + 29.22307272*log(Ts)**3/Ts - 293.86871481*log(Ts)**2/Ts + 1290.42340436*log(Ts)/Ts - 2087.25937009/Ts)*(0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752)/(total_mols*(-0.000380956854445275*Ru + 8.1710365221206e-5*log(Ts)**5 - 0.00278318246515918*log(Ts)**4 + 0.0373171004046311*log(Ts)**3 - 0.245797820513774*log(Ts)**2 + 0.795155764040913*log(Ts) - 1)**2) + NN2*(1.07243595*log(Ts)**4/Ts - 29.22307272*log(Ts)**3/Ts + 293.86871481*log(Ts)**2/Ts - 1290.42340436*log(Ts)/Ts + 2087.25937009/Ts)/(total_mols*(-Ru + 0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752)) + 1.83721491689154e-7*NNO*(-1.1074308*log(Ts)**4/Ts + 29.44892436*log(Ts)**3/Ts - 288.45020241*log(Ts)**2/Ts + 1230.6755001*log(Ts)/Ts - 1927.6639665/Ts)*(0.22148616*log(Ts)**5 - 7.36223109*log(Ts)**4 + 96.15006747*log(Ts)**3 - 615.33775005*log(Ts)**2 + 1927.6639665*log(Ts) - 2333.02836259)/(total_mols*(-0.000428627450928138*Ru + 9.49350481766618e-5*log(Ts)**5 - 0.00315565434525059*log(Ts)**4 + 0.0412125583262346*log(Ts)**3 - 0.263750651263787*log(Ts)**2 + 0.826249692206919*log(Ts) - 1)**2) + NNO*(1.1074308*log(Ts)**4/Ts - 29.44892436*log(Ts)**3/Ts + 288.45020241*log(Ts)**2/Ts - 1230.6755001*log(Ts)/Ts + 1927.6639665/Ts)/(total_mols*(-Ru + 0.22148616*log(Ts)**5 - 7.36223109*log(Ts)**4 + 96.15006747*log(Ts)**3 - 615.33775005*log(Ts)**2 + 1927.6639665*log(Ts) - 2333.02836259)) + 1.78636579731184e-6*NO*(-0.3513745*log(Ts)**4/Ts + 9.25340332*log(Ts)**3/Ts - 90.57574986*log(Ts)**2/Ts + 389.53047928*log(Ts)/Ts - 618.57117075/Ts)*(0.0702749*log(Ts)**5 - 2.31335083*log(Ts)**4 + 30.19191662*log(Ts)**3 - 194.76523964*log(Ts)**2 + 618.57117075*log(Ts) - 748.19500164)/(total_mols*(-0.00133654996064937*Ru + 9.39259148296387e-5*log(Ts)**5 - 0.0030919089608047*log(Ts)**4 + 0.0403530049703902*log(Ts)**3 - 0.260313473376708*log(Ts)**2 + 0.826751273924749*log(Ts) - 1)**2) + NO*(0.3513745*log(Ts)**4/Ts - 9.25340332*log(Ts)**3/Ts + 90.57574986*log(Ts)**2/Ts - 389.53047928*log(Ts)/Ts + 618.57117075/Ts)/(total_mols*(-Ru + 0.0702749*log(Ts)**5 - 2.31335083*log(Ts)**4 + 30.19191662*log(Ts)**3 - 194.76523964*log(Ts)**2 + 618.57117075*log(Ts) - 748.19500164)) + 7.0274904719501e-6*NO2*(-0.30903605*log(Ts)**4/Ts + 7.6723316*log(Ts)**3/Ts - 69.54989259*log(Ts)**2/Ts + 270.09592622*log(Ts)/Ts - 377.2244789/Ts)*(0.06180721*log(Ts)**5 - 1.9180829*log(Ts)**4 + 23.18329753*log(Ts)**3 - 135.04796311*log(Ts)**2 + 377.2244789*log(Ts) - 373.49246047)/(total_mols*(-0.00265094143125609*Ru + 0.000163847293739346*log(Ts)**5 - 0.00508472542819384*log(Ts)**4 + 0.0614575639354141*log(Ts)**3 - 0.358004240615044*log(Ts)**2 + log(Ts) - 0.990106637721702)**2) + NO2*(0.30903605*log(Ts)**4/Ts - 7.6723316*log(Ts)**3/Ts + 69.54989259*log(Ts)**2/Ts - 270.09592622*log(Ts)/Ts + 377.2244789/Ts)/(total_mols*(-Ru + 0.06180721*log(Ts)**5 - 1.9180829*log(Ts)**4 + 23.18329753*log(Ts)**3 - 135.04796311*log(Ts)**2 + 377.2244789*log(Ts) - 373.49246047)) + 1.40113809460065e-6*NOH*(-0.17390055*log(Ts)**4/Ts + 5.9871498*log(Ts)**3/Ts - 71.87458827*log(Ts)**2/Ts + 360.66895292*log(Ts)/Ts - 643.43907081/Ts)*(0.03478011*log(Ts)**5 - 1.49678745*log(Ts)**4 + 23.95819609*log(Ts)**3 - 180.33447646*log(Ts)**2 + 643.43907081*log(Ts) - 844.81094064)/(total_mols*(-0.00118369679166611*Ru + 4.11691046207945e-5*log(Ts)**5 - 0.0017717425023711*log(Ts)**4 + 0.0283592398458406*log(Ts)**3 - 0.21346134121249*log(Ts)**2 + 0.761636763750423*log(Ts) - 1)**2) + NOH*(0.17390055*log(Ts)**4/Ts - 5.9871498*log(Ts)**3/Ts + 71.87458827*log(Ts)**2/Ts - 360.66895292*log(Ts)/Ts + 643.43907081/Ts)/(total_mols*(-Ru + 0.03478011*log(Ts)**5 - 1.49678745*log(Ts)**4 + 23.95819609*log(Ts)**3 - 180.33447646*log(Ts)**2 + 643.43907081*log(Ts) - 844.81094064)) + Ncomb*(-a1/Ts - 2*a2*log(Ts)/Ts - 3*a3*log(Ts)**2/Ts - 4*a4*log(Ts)**3/Ts - 5*a5*log(Ts)**4/Ts)*(a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)/(total_mols*(-Ru + a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)**2) + Ncomb*(a1/Ts + 2*a2*log(Ts)/Ts + 3*a3*log(Ts)**2/Ts + 4*a4*log(Ts)**3/Ts + 5*a5*log(Ts)**4/Ts)/(total_mols*(-Ru + a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5))
    
    return Kt,dK


def wiebe(teta):
    return 1-exp(-a*((teta-teta_comb)/delta_teta)**(m+1))

def dx(teta):
    return -exp(-a*((teta-teta_comb)/delta_teta)**(m+1))*(-a*(m+1)*((teta-teta_comb)/delta_teta)**m)*1/delta_teta

def dc (Concentrations,t,T):
    C=Concentrations[0]
    O2=Concentrations[1]
    CO=Concentrations[2]
    CO2=Concentrations[3]
    H2O=Concentrations[4]
    
    ka=Ae*exp(-Ea/T)*np.sign(C)*abs(C)**me*O2**ne
    kbf=10**14.6*exp(-20131/T)*CO*H2O**0.5*O2**0.25
    kbr=5*10**8*exp(-20131/T)*CO2
    
    dC=-ka/N
    dO2=(-(c/2+h/4)*ka-0.5*c*kbf+c*0.5*kbr)/N
    dCO=(-kbf+ka+kbr)*c/N
    dCO2=(kbf-kbr)*c/N
    dH2O=h/2*ka/N
    return dC,dO2,dCO,dCO2,dH2O

def twoStepKinects (Concentrations0,ts,Temperature,Pressure):
    concentration = odeint(dc,Concentrations0,ts,args=(Temperature,))[-1]
    ro = Pressure * 10**-3 * M_m /(Ru * Temperature)
    dx = (-dc(concentration,ts[-1],Temperature)[0]*M_c *10**3/ ro)/Y0
    mass_burned = (Y0 - concentration[0] * M_c *10**3/ro)/Y0
    return mass_burned,dx,concentration



def woschni(P,Pm,T):
    vg=2.28*vp+0.00324*(P-Pm)*Vd*T0/(P0*V0)
    h=3.26*D**-0.2*(P*10**-3)**0.8*T**-0.55*vg**0.8
    return h

def volume(teta):
    s=R*cos(teta)+sqrt(L**2-R**2*(sin(teta))**2)  #deslocamento
    V=(np.pi*D**2/4)*(L+R-s+2*R/(r-1)) #volume 
    return V

def dvdt(teta):
    return 0.785398163397448*D**2*(R**2*sin(teta)*cos(teta)/sqrt(L**2 - R**2*sin(teta)**2) + R*sin(teta))

def area(teta):
    s=R*cos(teta)+sqrt(L**2-R**2*(sin(teta))**2)  #deslocamento
    A=np.pi*D*(D/2+L+R-s+(2*R/(r-1)))  #area das paredes
    return A
    

def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    if array[idx]<=value:
        return idx
    else:
        return idx-1


V0=volume(teta0)
T0=P0*10**-3*V0/(m_m*Rg)  

def motored(t,x,spark):
    P=x[0]
    T=x[1]
    Qa=x[2]
    Qp=x[3]
    W=x[4]
    
    A=area(t)
    V=volume(t)
    dVdt=dvdt(t)  #printar dvdt e ver os intervalos em que e desprezivel
    
    k=KReag(T)
    dkdT=dkreag(T)
    
    dxdt=0
    h=0
    
    dWdt=P*dVdt  #joule
    dQpdt=0 #rads por s
    dTdt=(((1/(P*V))*(Qtot*dxdt-fcor*dQpdt+spark)-dVdt/V)*(k-1))/((1/T)-(1/(k-1))*dkdT)
    dPdt=((((Qtot*dxdt-fcor*dQpdt+spark)-P*dVdt)*(k-1))-P*dVdt+(P*V*dkdT*dTdt/(k-1)))/V
    dQadt=P*dVdt+(1/(k-1))*(V*dPdt+P*dVdt-(P*V*dkdT*dTdt/(k-1)))
    
    return(dPdt,dTdt,dQadt,dQpdt,dWdt)  

    
def motor(t,x,spark):
    P=x[0]
    T=x[1]
    Qa=x[2]
    Qp=x[3]
    W=x[4]
    
    A=area(t)
    V=volume(t)
    dVdt=dvdt(t)
    
    dxdt=0
    
    k=KReag(T)
    dkdT=dkreag(T)
    
    vg=2.28*vp
    h=3.26*D**(-0.2)*(P*10**-3)**0.8*T**(-0.53)*vg**0.8
    
    dWdt=P*dVdt  #joule
    dQpdt=(h*A*(T-Tp)/(N*fcor)) #rads por s
    dTdt=(((1/(P*V))*(Qtot*dxdt-fcor*dQpdt+spark)-dVdt/V)*(k-1))/((1/T)-(1/(k-1))*dkdT)
    dPdt=((((Qtot*dxdt-fcor*dQpdt + spark)-P*dVdt)*(k-1))-P*dVdt+(P*V*dkdT*dTdt/(k-1)))/V
    dQadt=P*dVdt+(1/(k-1))*(V*dPdt+P*dVdt-(P*V*dkdT*dTdt/(k-1)))
    
    return(dPdt,dTdt,dQadt,dQpdt,dWdt)  

def motorComb(t,x,Pm,xt,dxdt):
    P=x[0]
    T=x[1]
    Qa=x[2]
    Qp=x[3]
    W=x[4]
    
    A=area(t)
    V=volume(t)
    dVdt=dvdt(t)
    
    k=xt*KProd(T)+(1-xt)*KReag(T)
    dkdT=xt*dkprod(T)+(1-xt)*dkreag(T)
    
    #dxdt=dx(t)

    vg=2.28*vp+0.00324*(P-Pm)*Vd*T0/(P0*V0)
    h=3.26*D**-0.2*(P*10**-3)**0.8*T**(-0.53)*vg**0.8
  
    dWdt=P*dVdt  #joule
    dQpdt=(h*A*(T-Tp)/(N*fcor)) #rads por s
    dTdt=(((1/(P*V))*(Qtot*dxdt-fcor*dQpdt)-dVdt/V)*(k-1))/((1/T)-(1/(k-1))*dkdT)
    dPdt=((((Qtot*dxdt-fcor*dQpdt)-P*dVdt)*(k-1))-P*dVdt+(P*V*dkdT*dTdt/(k-1)))/V
    dQadt=P*dVdt+(1/(k-1))*(V*dPdt+P*dVdt-(P*V*dkdT*dTdt/(k-1)))
    
    return(dPdt,dTdt,dQadt,dQpdt,dWdt)      

def motoredEq(t,x,reactants_composition):
    P=x[0]
    T=x[1]
    Qa=x[2]
    Qp=x[3]
    W=x[4]
    
    A=area(t)
    V=volume(t)
    dVdt=dvdt(t)  #printar dvdt e ver os intervalos em que e desprezivel
    
    k, dkdT=K(T, reactants_composition)
    
    dxdt=0
    h=0
    
    dWdt=P*dVdt  #joule
    dQpdt=0 #rads por s
    dTdt=(((1/(P*V))*(Qtot*dxdt-fcor*dQpdt)-dVdt/V)*(k-1))/((1/T)-(1/(k-1))*dkdT)
    dPdt=((((Qtot*dxdt-fcor*dQpdt)-P*dVdt)*(k-1))-P*dVdt+(P*V*dkdT*dTdt/(k-1)))/V
    dQadt=P*dVdt+(1/(k-1))*(V*dPdt+P*dVdt-(P*V*dkdT*dTdt/(k-1)))
    
    return(dPdt,dTdt,dQadt,dQpdt,dWdt)  

    
def motorEq(t,x,reactants_composition):
    P=x[0]
    T=x[1]
    Qa=x[2]
    Qp=x[3]
    W=x[4]
    
    A=area(t)
    V=volume(t)
    dVdt=dvdt(t)
    
    dxdt=0
    
    k, dkdT=K(T, reactants_composition)
    
    vg=2.28*vp
    h=3.26*D**(-0.2)*(P*10**-3)**0.8*T**(-0.55)*vg**0.8
    
    dWdt=P*dVdt  #joule
    dQpdt=(h*A*(T-Tp)/(N*fcor)) #rads por s
    dTdt=(((1/(P*V))*(Qtot*dxdt-fcor*dQpdt)-dVdt/V)*(k-1))/((1/T)-(1/(k-1))*dkdT)
    dPdt=((((Qtot*dxdt-fcor*dQpdt)-P*dVdt)*(k-1))-P*dVdt+(P*V*dkdT*dTdt/(k-1)))/V
    dQadt=P*dVdt+(1/(k-1))*(V*dPdt+P*dVdt-(P*V*dkdT*dTdt/(k-1)))
    
    return(dPdt,dTdt,dQadt,dQpdt,dWdt)  

def motorCombEquilibrium(t,x,Pm,composition0,composition1,dt):
    P=x[0]
    T=x[1]
    Qa=x[2]
    Qp=x[3]
    W=x[4]
    
    A=area(t)
    V=volume(t)
    dVdt=dvdt(t)
    
    vg=2.28*vp+0.00324*(P-Pm)*Vd*T0/(P0*V0)
    h=3.26*D**-0.2*(P*10**-3)**0.8*T**-0.55*vg**0.8
    
    dQ = heatRelease(composition0,composition1,T)/dt
    
    k,dkdT= K(T,composition0)
    
    dWdt=P*dVdt  #joule
    dQpdt=(h*A*(T-Tp)/(N*fcor)) #rads por s
    dTdt=(((1/(P*V))*(dQ-fcor*dQpdt)-dVdt/V)*(k-1))/((1/T)-(1/(k-1))*dkdT)
    dPdt=((((dQ-fcor*dQpdt)-P*dVdt)*(k-1))-P*dVdt+(P*V*dkdT*dTdt/(k-1)))/V
    dQadt=P*dVdt+(1/(k-1))*(V*dPdt+P*dVdt-(P*V*dkdT*dTdt/(k-1)))
    
    return(dPdt,dTdt,dQadt,dQpdt,dWdt)           

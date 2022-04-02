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

def KReag(Ts):
    Cp=a0+a1*log(Ts)+a2*log(Ts)**2+a3*log(Ts)**3+a4*log(Ts)**4+a5*log(Ts)**5
    CpO=10228.3426-7184.92333*log(Ts)+2010.86808*log(Ts)**2-279.69496*log(Ts)**3+19.34823*log(Ts)**4-0.53257*log(Ts)**5
    CpN=0.54497*log(Ts)**5-18.69984*log(Ts)**4+254.29554*log(Ts)**3-1712.17390*log(Ts)**2+5708.38047*log(Ts)-7513.3642
    KReag=((Cp/(Cp-Ru))+OC*(CpO/(CpO-Ru))+3.76*OC*(CpN/(CpN-Ru)))/(1+OC+3.76*OC)
    return KReag
    
def KProd(Ts):
    CpN=0.54497*log(Ts)**5-18.69984*log(Ts)**4+254.29554*log(Ts)**3-1712.17390*log(Ts)**2+5708.38047*log(Ts)-7513.3642
    CpCO2=0.20754*log(Ts)**5-6.43522*log(Ts)**4+77.54809*log(Ts)**3-452.81197*log(Ts)**2+1288.4677*log(Ts)-1412.36785
    CpH2O=0.64541*log(Ts)**5-23.54277*log(Ts)**4+339.33662*log(Ts)**3-2414.77575*log(Ts)**2+8490.5218*log(Ts)-11780.76495
    KProd=((CpCO2/(CpCO2-Ru))*stoichometricCO2+(CpH2O/(CpH2O-Ru))*stoichometricH2O+(CpN/(CpN-Ru))*OC*3.76)/(stoichometricCO2+stoichometricH2O+stoichometricN2)
    return KProd

def dkreag(Ts):
    dkr=(6.66068603173739e-8*OC*(-2.72485*log(Ts)**4/Ts + 74.79936*log(Ts)**3/Ts - 762.88662*log(Ts)**2/Ts + 3424.3478*log(Ts)/Ts - 5708.38047/Ts)*(0.54497*log(Ts)**5 - 18.69984*log(Ts)**4 + 254.29554*log(Ts)**3 - 1712.1739*log(Ts)**2 + 5708.38047*log(Ts) - 7513.3642)/(-0.000133096170155042*Ru + 7.25334198493932e-5*log(Ts)**5 - 0.00248887708651206*log(Ts)**4 + 0.0338457624615083*log(Ts)**3 - 0.227883788729422*log(Ts)**2 + 0.759763578344838*log(Ts) - 1)**2 + OC*(-2.66285*log(Ts)**4/Ts + 77.39292*log(Ts)**3/Ts - 839.08488*log(Ts)**2/Ts + 4021.73616*log(Ts)/Ts - 7184.92333/Ts)/(-Ru - 0.53257*log(Ts)**5 + 19.34823*log(Ts)**4 - 279.69496*log(Ts)**3 + 2010.86808*log(Ts)**2 - 7184.92333*log(Ts) + 10228.3426) + 9.55849389871466e-9*OC*(2.66285*log(Ts)**4/Ts - 77.39292*log(Ts)**3/Ts + 839.08488*log(Ts)**2/Ts - 4021.73616*log(Ts)/Ts + 7184.92333/Ts)*(-0.53257*log(Ts)**5 + 19.34823*log(Ts)**4 - 279.69496*log(Ts)**3 + 2010.86808*log(Ts)**2 - 7184.92333*log(Ts) + 10228.3426)/(-9.77675503360632e-5*Ru - 5.20680642824772e-5*log(Ts)**5 + 0.00189162905043873*log(Ts)**4 - 0.0273450910805432*log(Ts)**3 + 0.196597646230583*log(Ts)**2 - 0.70245235332653*log(Ts) + 1)**2 + 3.76*OC*(2.72485*log(Ts)**4/Ts - 74.79936*log(Ts)**3/Ts + 762.88662*log(Ts)**2/Ts - 3424.3478*log(Ts)/Ts + 5708.38047/Ts)/(-Ru + 0.54497*log(Ts)**5 - 18.69984*log(Ts)**4 + 254.29554*log(Ts)**3 - 1712.1739*log(Ts)**2 + 5708.38047*log(Ts) - 7513.3642) + (-a1/Ts - 2*a2*log(Ts)/Ts - 3*a3*log(Ts)**2/Ts - 4*a4*log(Ts)**3/Ts - 5*a5*log(Ts)**4/Ts)*(a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)/(-Ru + a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)**2 + (a1/Ts + 2*a2*log(Ts)/Ts + 3*a3*log(Ts)**2/Ts + 4*a4*log(Ts)**3/Ts + 5*a5*log(Ts)**4/Ts)/(-Ru + a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5))/(4.76*OC + 1)
    return dkr

def dkprod(Ts):
    dkp=(6.66068603173739e-8*OC*(-2.72485*log(Ts)**4/Ts + 74.79936*log(Ts)**3/Ts - 762.88662*log(Ts)**2/Ts + 3424.3478*log(Ts)/Ts - 5708.38047/Ts)*(0.54497*log(Ts)**5 - 18.69984*log(Ts)**4 + 254.29554*log(Ts)**3 - 1712.1739*log(Ts)**2 + 5708.38047*log(Ts) - 7513.3642)/(-0.000133096170155042*Ru + 7.25334198493932e-5*log(Ts)**5 - 0.00248887708651206*log(Ts)**4 + 0.0338457624615083*log(Ts)**3 - 0.227883788729422*log(Ts)**2 + 0.759763578344838*log(Ts) - 1)**2 + 3.76*OC*(2.72485*log(Ts)**4/Ts - 74.79936*log(Ts)**3/Ts + 762.88662*log(Ts)**2/Ts - 3424.3478*log(Ts)/Ts + 5708.38047/Ts)/(-Ru + 0.54497*log(Ts)**5 - 18.69984*log(Ts)**4 + 254.29554*log(Ts)**3 - 1712.1739*log(Ts)**2 + 5708.38047*log(Ts) - 7513.3642) + 5.01307675179101e-7*stoichometricCO2*(-1.0377*log(Ts)**4/Ts + 25.74088*log(Ts)**3/Ts - 232.64427*log(Ts)**2/Ts + 905.62394*log(Ts)/Ts - 1288.4677/Ts)*(0.20754*log(Ts)**5 - 6.43522*log(Ts)**4 + 77.54809*log(Ts)**3 - 452.81197*log(Ts)**2 + 1288.4677*log(Ts) - 1412.36785)/(-0.000708030843381205*Ru + 0.000146944721235335*log(Ts)**5 - 0.0045563342439436*log(Ts)**4 + 0.0549064395653016*log(Ts)**3 - 0.320604841012205*log(Ts)**2 + 0.912274872300442*log(Ts) - 1)**2 + stoichometricCO2*(1.0377*log(Ts)**4/Ts - 25.74088*log(Ts)**3/Ts + 232.64427*log(Ts)**2/Ts - 905.62394*log(Ts)/Ts + 1288.4677/Ts)/(-Ru + 0.20754*log(Ts)**5 - 6.43522*log(Ts)**4 + 77.54809*log(Ts)**3 - 452.81197*log(Ts)**2 + 1288.4677*log(Ts) - 1412.36785) + 7.20531576341265e-9*stoichometricH2O*(-3.22705*log(Ts)**4/Ts + 94.17108*log(Ts)**3/Ts - 1018.00986*log(Ts)**2/Ts + 4829.5515*log(Ts)/Ts - 8490.5218/Ts)*(0.64541*log(Ts)**5 - 23.54277*log(Ts)**4 + 339.33662*log(Ts)**3 - 2414.77575*log(Ts)**2 + 8490.5218*log(Ts) - 11780.76495)/(-8.48841313992942e-5*Ru + 5.47850672464185e-5*log(Ts)**5 - 0.00199840758218336*log(Ts)**4 + 0.0288042942406724*log(Ts)**3 - 0.204976142062829*log(Ts)**2 + 0.720710568119772*log(Ts) - 1)**2 + stoichometricH2O*(3.22705*log(Ts)**4/Ts - 94.17108*log(Ts)**3/Ts + 1018.00986*log(Ts)**2/Ts - 4829.5515*log(Ts)/Ts + 8490.5218/Ts)/(-Ru + 0.64541*log(Ts)**5 - 23.54277*log(Ts)**4 + 339.33662*log(Ts)**3 - 2414.77575*log(Ts)**2 + 8490.5218*log(Ts) - 11780.76495))/(stoichometricCO2 + stoichometricH2O + stoichometricN2)
    return dkp

def K(Ts,composition):
    comb = composition[0]
    O2 = composition[1]
    CO2 = composition[2]
    H2O =composition[3]
    N2 = composition[4]
    
    Cp=a0+a1*log(Ts)+a2*log(Ts)**2+a3*log(Ts)**3+a4*log(Ts)**4+a5*log(Ts)**5
    CpO=10228.3426-7184.92333*log(Ts)+2010.86808*log(Ts)**2-279.69496*log(Ts)**3+19.34823*log(Ts)**4-0.53257*log(Ts)**5
    CpN=0.54497*log(Ts)**5-18.69984*log(Ts)**4+254.29554*log(Ts)**3-1712.17390*log(Ts)**2+5708.38047*log(Ts)-7513.3642
    CpCO2=0.20754*log(Ts)**5-6.43522*log(Ts)**4+77.54809*log(Ts)**3-452.81197*log(Ts)**2+1288.4677*log(Ts)-1412.36785
    CpH2O=0.64541*log(Ts)**5-23.54277*log(Ts)**4+339.33662*log(Ts)**3-2414.77575*log(Ts)**2+8490.5218*log(Ts)-11780.76495
    Kt = (comb*(Cp/(Cp-Ru)) + O2 * (CpO/(CpO-Ru)) + CO2 * (CpCO2/(CpCO2-Ru)) + H2O * (CpH2O/(CpH2O-Ru)) + N2 * (CpN/(CpN-Ru)))/(comb+O2+CO2+H2O+N2)
    
    dK = (5.01307675179101e-7*CO2*(-1.0377*log(Ts)**4/Ts + 25.74088*log(Ts)**3/Ts - 232.64427*log(Ts)**2/Ts + 905.62394*log(Ts)/Ts - 1288.4677/Ts)* \
          (0.20754*log(Ts)**5 - 6.43522*log(Ts)**4 + 77.54809*log(Ts)**3 - 452.81197*log(Ts)**2 + 1288.4677*log(Ts) - 1412.36785)/ \
          (-0.000708030843381205*Ru + 0.000146944721235335*log(Ts)**5 - 0.0045563342439436*log(Ts)**4 + 0.0549064395653016*log(Ts)**3 - \
           0.320604841012205*log(Ts)**2 + 0.912274872300442*log(Ts) - 1)**2 + CO2*(1.0377*log(Ts)**4/Ts - 25.74088*log(Ts)**3/Ts + 232.64427*log(Ts)**2/Ts - 905.62394*log(Ts)/Ts + 1288.4677/Ts)/ \
           (-Ru + 0.20754*log(Ts)**5 - 6.43522*log(Ts)**4 + 77.54809*log(Ts)**3 - 452.81197*log(Ts)**2 + 1288.4677*log(Ts) - 1412.36785) + \
           7.20531576341265e-9*H2O*(-3.22705*log(Ts)**4/Ts + 94.17108*log(Ts)**3/Ts - 1018.00986*log(Ts)**2/Ts + 4829.5515*log(Ts)/Ts - 8490.5218/Ts)* \
           (0.64541*log(Ts)**5 - 23.54277*log(Ts)**4 + 339.33662*log(Ts)**3 - 2414.77575*log(Ts)**2 + 8490.5218*log(Ts) - 11780.76495)/(-8.48841313992942e-5*Ru + \
           5.47850672464185e-5*log(Ts)**5 - 0.00199840758218336*log(Ts)**4 + 0.0288042942406724*log(Ts)**3 - 0.204976142062829*log(Ts)**2 + 0.720710568119772*log(Ts) - 1)**2 + \
           H2O*(3.22705*log(Ts)**4/Ts - 94.17108*log(Ts)**3/Ts + 1018.00986*log(Ts)**2/Ts - 4829.5515*log(Ts)/Ts + 8490.5218/Ts)/(-Ru + 0.64541*log(Ts)**5 - 23.54277*log(Ts)**4 + \
          339.33662*log(Ts)**3 - 2414.77575*log(Ts)**2 + 8490.5218*log(Ts) - 11780.76495) + 1.77145905099399e-8*N2*(-2.72485*log(Ts)**4/Ts + 74.79936*log(Ts)**3/Ts - \
        762.88662*log(Ts)**2/Ts + 3424.3478*log(Ts)/Ts - 5708.38047/Ts)*(0.54497*log(Ts)**5 - 18.69984*log(Ts)**4 + 254.29554*log(Ts)**3 - 1712.1739*log(Ts)**2 + 5708.38047*log(Ts) - \
            7513.3642)/(-0.000133096170155042*Ru + 7.25334198493932e-5*log(Ts)**5 - 0.00248887708651206*log(Ts)**4 + 0.0338457624615083*log(Ts)**3 - 0.227883788729422*log(Ts)**2 + \
            0.759763578344838*log(Ts) - 1)**2 + N2*(2.72485*log(Ts)**4/Ts - 74.79936*log(Ts)**3/Ts + 762.88662*log(Ts)**2/Ts - 3424.3478*log(Ts)/Ts + 5708.38047/Ts)/(-Ru + 0.54497*log(Ts)**5 - \
        18.69984*log(Ts)**4 + 254.29554*log(Ts)**3 - 1712.1739*log(Ts)**2 + 5708.38047*log(Ts) - 7513.3642) + O2*(-2.66285*log(Ts)**4/Ts + 77.39292*log(Ts)**3/Ts - 839.08488*log(Ts)**2/Ts + 4021.73616*log(Ts)/Ts - 7184.92333/Ts)/ \
        (-Ru - 0.53257*log(Ts)**5 + 19.34823*log(Ts)**4 - 279.69496*log(Ts)**3 + 2010.86808*log(Ts)**2 - 7184.92333*log(Ts) + 10228.3426) + 9.55849389871466e-9*O2*(2.66285*log(Ts)**4/Ts - 77.39292*log(Ts)**3/Ts + 839.08488*log(Ts)**2/Ts - 4021.73616*log(Ts)/Ts + 7184.92333/Ts)* \
        (-0.53257*log(Ts)**5 + 19.34823*log(Ts)**4 - 279.69496*log(Ts)**3 + 2010.86808*log(Ts)**2 - 7184.92333*log(Ts) + 10228.3426)/ \
        (-9.77675503360632e-5*Ru - 5.20680642824772e-5*log(Ts)**5 + 0.00189162905043873*log(Ts)**4 - 0.0273450910805432*log(Ts)**3 + 0.196597646230583*log(Ts)**2 - 0.70245235332653*log(Ts) + 1)**2 + \
        comb*(-a1/Ts - 2*a2*log(Ts)/Ts - 3*a3*log(Ts)**2/Ts - 4*a4*log(Ts)**3/Ts - 5*a5*log(Ts)**4/Ts)*(a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)/ \
        (-Ru + a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)**2 + \
        comb*(a1/Ts + 2*a2*log(Ts)/Ts + 3*a3*log(Ts)**2/Ts + 4*a4*log(Ts)**3/Ts + 5*a5*log(Ts)**4/Ts)/ \
        (-Ru + a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5))/ \
        (CO2 + H2O + N2 + O2 + comb)
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
    h=3.26*D**-0.2*(P*10**-3)**0.8*T**-0.53*vg**0.8
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
    dQpdt=(h*A*(T-Tp)/(N*fcor)) #rads por s
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
    h=3.26*D**-0.2*(P*10**-3)**0.8*T**-0.53*vg**0.8
  
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
    h=3.26*D**-0.2*(P*10**-3)**0.8*T**-0.53*vg**0.8
    
    dQ = heatRelease(composition0,composition1,T)/dt
    
    k,dkdT= K(T,[composition0[0],composition0[1],composition0[3],composition0[2],composition0[12]])
    
    dWdt=P*dVdt  #joule
    dQpdt=(h*A*(T-Tp)/(N*fcor)) #rads por s
    dTdt=(((1/(P*V))*(dQ-fcor*dQpdt)-dVdt/V)*(k-1))/((1/T)-(1/(k-1))*dkdT)
    dPdt=((((dQ-fcor*dQpdt)-P*dVdt)*(k-1))-P*dVdt+(P*V*dkdT*dTdt/(k-1)))/V
    dQadt=P*dVdt+(1/(k-1))*(V*dPdt+P*dVdt-(P*V*dkdT*dTdt/(k-1)))
    
    return(dPdt,dTdt,dQadt,dQpdt,dWdt)           

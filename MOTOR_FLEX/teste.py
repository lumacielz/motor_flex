#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 12 08:23:38 2022

@author: luiza.maciel
"""


from constants import *
from equilibrium import *
from scipy.integrate import solve_ivp,odeint

def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    if array[idx]<=value:
        return idx
    else:
        return idx-1
def KReag(Ts):
    Cp=a0+a1*log(Ts)+a2*log(Ts)**2+a3*log(Ts)**3+a4*log(Ts)**4+a5*log(Ts)**5
    CpO=0.06180721*log(Ts)**5-1.91808290*log(Ts)**4+23.18329753*log(Ts)**3-135.04796311*log(Ts)**2+377.22447890*log(Ts)-373.49246047
    CpN=0.21448719*log(Ts)**5 -7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 -645.21170218*log(Ts)**2 +2087.25937009*log(Ts) -2624.96917520
    Cpt =Cp/(1+OC+OC*3.76) + OC/(1+OC+OC*3.76) * CpO + CpN*OC*3.76/(1+OC+OC*3.76)
    KReag= ((Cp/(Cp-Ru))+OC*(CpO/(CpO-Ru))+OC*3.76*(CpN/(CpN-Ru)))/(1+OC+OC*3.76)
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

def volume(teta):
    s=R*cos(teta)+sqrt(L**2-R**2*(sin(teta))**2)  #deslocamento
    V=(np.pi*D**2/4)*(L+R-s+2*R/(r-1)) #volume 
    return V

def dvdt(teta):
    return np.pi*D**2*(R**2*sin(teta)*cos(teta)/sqrt(L**2 - R**2*sin(teta)**2) + R*sin(teta))/4

def area(teta):
    s=R*cos(teta)+sqrt(L**2-R**2*(sin(teta))**2)  #deslocamento
    A=np.pi*D*(D/2+L+R-s+(2*R/(r-1)))  #area das paredes
    return A
    
def dkreag(Ts):
    #usando ks
    dkr=(5.45681749807633e-7*OC*(-1.07243595*log(Ts)**4/Ts + 29.22307272*log(Ts)**3/Ts - 293.86871481*log(Ts)**2/Ts + 1290.42340436*log(Ts)/Ts - 2087.25937009/Ts)*(0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752)/(-0.000380956854445275*Ru + 8.1710365221206e-5*log(Ts)**5 - 0.00278318246515918*log(Ts)**4 + 0.0373171004046311*log(Ts)**3 - 0.245797820513774*log(Ts)**2 + 0.795155764040913*log(Ts) - 1)**2 + 7.0274904719501e-6*OC*(-0.30903605*log(Ts)**4/Ts + 7.6723316*log(Ts)**3/Ts - 69.54989259*log(Ts)**2/Ts + 270.09592622*log(Ts)/Ts - 377.2244789/Ts)*(0.06180721*log(Ts)**5 - 1.9180829*log(Ts)**4 + 23.18329753*log(Ts)**3 - 135.04796311*log(Ts)**2 + 377.2244789*log(Ts) - 373.49246047)/(-0.00265094143125609*Ru + 0.000163847293739346*log(Ts)**5 - 0.00508472542819384*log(Ts)**4 + 0.0614575639354141*log(Ts)**3 - 0.358004240615044*log(Ts)**2 + log(Ts) - 0.990106637721702)**2 + OC*(0.30903605*log(Ts)**4/Ts - 7.6723316*log(Ts)**3/Ts + 69.54989259*log(Ts)**2/Ts - 270.09592622*log(Ts)/Ts + 377.2244789/Ts)/(-Ru + 0.06180721*log(Ts)**5 - 1.9180829*log(Ts)**4 + 23.18329753*log(Ts)**3 - 135.04796311*log(Ts)**2 + 377.2244789*log(Ts) - 373.49246047) + 3.76*OC*(1.07243595*log(Ts)**4/Ts - 29.22307272*log(Ts)**3/Ts + 293.86871481*log(Ts)**2/Ts - 1290.42340436*log(Ts)/Ts + 2087.25937009/Ts)/(-Ru + 0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752) + (-a1/Ts - 2*a2*log(Ts)/Ts - 3*a3*log(Ts)**2/Ts - 4*a4*log(Ts)**3/Ts - 5*a5*log(Ts)**4/Ts)*(a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)/(-Ru + a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)**2 + (a1/Ts + 2*a2*log(Ts)/Ts + 3*a3*log(Ts)**2/Ts + 4*a4*log(Ts)**3/Ts + 5*a5*log(Ts)**4/Ts)/(-Ru + a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5))/(4.76*OC + 1)
    #usando cpt
    #dkr =(-OC*(0.30903605*log(Ts)**4/Ts - 7.6723316*log(Ts)**3/Ts + 69.54989259*log(Ts)**2/Ts - 270.09592622*log(Ts)/Ts + 377.2244789/Ts)/(OC + stoichometricN2 + 1) - stoichometricN2*(1.07243595*log(Ts)**4/Ts - 29.22307272*log(Ts)**3/Ts + 293.86871481*log(Ts)**2/Ts - 1290.42340436*log(Ts)/Ts + 2087.25937009/Ts)/(OC + stoichometricN2 + 1) - (a1/Ts + 2*a2*log(Ts)/Ts + 3*a3*log(Ts)**2/Ts + 4*a4*log(Ts)**3/Ts + 5*a5*log(Ts)**4/Ts)/(OC + stoichometricN2 + 1))*(OC*(0.06180721*log(Ts)**5 - 1.9180829*log(Ts)**4 + 23.18329753*log(Ts)**3 - 135.04796311*log(Ts)**2 + 377.2244789*log(Ts) - 373.49246047)/(OC + stoichometricN2 + 1) + stoichometricN2*(0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752)/(OC + stoichometricN2 + 1) + (a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)/(OC + stoichometricN2 + 1))/(OC*(0.06180721*log(Ts)**5 - 1.9180829*log(Ts)**4 + 23.18329753*log(Ts)**3 - 135.04796311*log(Ts)**2 + 377.2244789*log(Ts) - 373.49246047)/(OC + stoichometricN2 + 1) - Ru + stoichometricN2*(0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752)/(OC + stoichometricN2 + 1) + (a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)/(OC + stoichometricN2 + 1))**2 + (OC*(0.30903605*log(Ts)**4/Ts - 7.6723316*log(Ts)**3/Ts + 69.54989259*log(Ts)**2/Ts - 270.09592622*log(Ts)/Ts + 377.2244789/Ts)/(OC + stoichometricN2 + 1) + stoichometricN2*(1.07243595*log(Ts)**4/Ts - 29.22307272*log(Ts)**3/Ts + 293.86871481*log(Ts)**2/Ts - 1290.42340436*log(Ts)/Ts + 2087.25937009/Ts)/(OC + stoichometricN2 + 1) + (a1/Ts + 2*a2*log(Ts)/Ts + 3*a3*log(Ts)**2/Ts + 4*a4*log(Ts)**3/Ts + 5*a5*log(Ts)**4/Ts)/(OC + stoichometricN2 + 1))/(OC*(0.06180721*log(Ts)**5 - 1.9180829*log(Ts)**4 + 23.18329753*log(Ts)**3 - 135.04796311*log(Ts)**2 + 377.2244789*log(Ts) - 373.49246047)/(OC + stoichometricN2 + 1) - Ru + stoichometricN2*(0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752)/(OC + stoichometricN2 + 1) + (a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)/(OC + stoichometricN2 + 1))
    return dkr

def dkprod(Ts):
    #usando ks
    dkp= (1.45164604139169e-5*stoichometricCO2*(-0.4997566*log(Ts)**4/Ts + 11.38794476*log(Ts)**3/Ts - 89.91527472*log(Ts)**2/Ts + 278.7309764*log(Ts)/Ts - 262.46393332/Ts)*(0.09995132*log(Ts)**5 - 2.84698619*log(Ts)**4 + 29.97175824*log(Ts)**3 - 139.3654882*log(Ts)**2 + 262.46393332*log(Ts) - 77.60476812)/(-0.00381004729812071*Ru + 0.000380819256709598*log(Ts)**5 - 0.0108471520409965*log(Ts)**4 + 0.114193816502239*log(Ts)**3 - 0.530989101767684*log(Ts)**2 + log(Ts) - 0.29567783709689)**2 + stoichometricCO2*(0.4997566*log(Ts)**4/Ts - 11.38794476*log(Ts)**3/Ts + 89.91527472*log(Ts)**2/Ts - 278.7309764*log(Ts)/Ts + 262.46393332/Ts)/(-Ru + 0.09995132*log(Ts)**5 - 2.84698619*log(Ts)**4 + 29.97175824*log(Ts)**3 - 139.3654882*log(Ts)**2 + 262.46393332*log(Ts) - 77.60476812) + 1.52739943457638e-7*stoichometricH2O*(-0.7452238*log(Ts)**4/Ts + 22.4180952*log(Ts)**3/Ts - 244.6266093*log(Ts)**2/Ts + 1146.87976368*log(Ts)/Ts - 1954.64102748/Ts)*(0.14904476*log(Ts)**5 - 5.6045238*log(Ts)**4 + 81.5422031*log(Ts)**3 - 573.43988184*log(Ts)**2 + 1954.64102748*log(Ts) - 2558.72544088)/(-0.000390819579163632*Ru + 5.82496103797445e-5*log(Ts)**5 - 0.00219035763292856*log(Ts)**4 + 0.0318682894996174*log(Ts)**3 - 0.224111533296352*log(Ts)**2 + 0.763911983775703*log(Ts) - 1)**2 + stoichometricH2O*(0.7452238*log(Ts)**4/Ts - 22.4180952*log(Ts)**3/Ts + 244.6266093*log(Ts)**2/Ts - 1146.87976368*log(Ts)/Ts + 1954.64102748/Ts)/(-Ru + 0.14904476*log(Ts)**5 - 5.6045238*log(Ts)**4 + 81.5422031*log(Ts)**3 - 573.43988184*log(Ts)**2 + 1954.64102748*log(Ts) - 2558.72544088) + 1.45128124948838e-7*stoichometricN2*(-1.07243595*log(Ts)**4/Ts + 29.22307272*log(Ts)**3/Ts - 293.86871481*log(Ts)**2/Ts + 1290.42340436*log(Ts)/Ts - 2087.25937009/Ts)*(0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752)/(-0.000380956854445275*Ru + 8.1710365221206e-5*log(Ts)**5 - 0.00278318246515918*log(Ts)**4 + 0.0373171004046311*log(Ts)**3 - 0.245797820513774*log(Ts)**2 + 0.795155764040913*log(Ts) - 1)**2 + stoichometricN2*(1.07243595*log(Ts)**4/Ts - 29.22307272*log(Ts)**3/Ts + 293.86871481*log(Ts)**2/Ts - 1290.42340436*log(Ts)/Ts + 2087.25937009/Ts)/(-Ru + 0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752))/(stoichometricCO2 + stoichometricH2O + stoichometricN2)
    #usando cpt
    #dkp =(stoichometricCO2*(0.4997566*log(Ts)**4/Ts - 11.38794476*log(Ts)**3/Ts + 89.91527472*log(Ts)**2/Ts - 278.7309764*log(Ts)/Ts + 262.46393332/Ts) + stoichometricH2O*(0.7452238*log(Ts)**4/Ts - 22.4180952*log(Ts)**3/Ts + 244.6266093*log(Ts)**2/Ts - 1146.87976368*log(Ts)/Ts + 1954.64102748/Ts) + stoichometricN2*(1.07243595*log(Ts)**4/Ts - 29.22307272*log(Ts)**3/Ts + 293.86871481*log(Ts)**2/Ts - 1290.42340436*log(Ts)/Ts + 2087.25937009/Ts))/((-Ru + (stoichometricCO2*(0.09995132*log(Ts)**5 - 2.84698619*log(Ts)**4 + 29.97175824*log(Ts)**3 - 139.3654882*log(Ts)**2 + 262.46393332*log(Ts) - 77.60476812) + stoichometricH2O*(0.14904476*log(Ts)**5 - 5.6045238*log(Ts)**4 + 81.5422031*log(Ts)**3 - 573.43988184*log(Ts)**2 + 1954.64102748*log(Ts) - 2558.72544088) + stoichometricN2*(0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752))/(stoichometricCO2 + stoichometricH2O + stoichometricN2))*(stoichometricCO2 + stoichometricH2O + stoichometricN2)) - (stoichometricCO2*(0.4997566*log(Ts)**4/Ts - 11.38794476*log(Ts)**3/Ts + 89.91527472*log(Ts)**2/Ts - 278.7309764*log(Ts)/Ts + 262.46393332/Ts) + stoichometricH2O*(0.7452238*log(Ts)**4/Ts - 22.4180952*log(Ts)**3/Ts + 244.6266093*log(Ts)**2/Ts - 1146.87976368*log(Ts)/Ts + 1954.64102748/Ts) + stoichometricN2*(1.07243595*log(Ts)**4/Ts - 29.22307272*log(Ts)**3/Ts + 293.86871481*log(Ts)**2/Ts - 1290.42340436*log(Ts)/Ts + 2087.25937009/Ts))*(stoichometricCO2*(0.09995132*log(Ts)**5 - 2.84698619*log(Ts)**4 + 29.97175824*log(Ts)**3 - 139.3654882*log(Ts)**2 + 262.46393332*log(Ts) - 77.60476812) + stoichometricH2O*(0.14904476*log(Ts)**5 - 5.6045238*log(Ts)**4 + 81.5422031*log(Ts)**3 - 573.43988184*log(Ts)**2 + 1954.64102748*log(Ts) - 2558.72544088) + stoichometricN2*(0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752))/((-Ru + (stoichometricCO2*(0.09995132*log(Ts)**5 - 2.84698619*log(Ts)**4 + 29.97175824*log(Ts)**3 - 139.3654882*log(Ts)**2 + 262.46393332*log(Ts) - 77.60476812) + stoichometricH2O*(0.14904476*log(Ts)**5 - 5.6045238*log(Ts)**4 + 81.5422031*log(Ts)**3 - 573.43988184*log(Ts)**2 + 1954.64102748*log(Ts) - 2558.72544088) + stoichometricN2*(0.21448719*log(Ts)**5 - 7.30576818*log(Ts)**4 + 97.95623827*log(Ts)**3 - 645.21170218*log(Ts)**2 + 2087.25937009*log(Ts) - 2624.9691752))/(stoichometricCO2 + stoichometricH2O + stoichometricN2))**2*(stoichometricCO2 + stoichometricH2O + stoichometricN2)**2)
    return dkp


def dc (t,Concentrations,T):
    C=Concentrations[0]
    O2=Concentrations[1]
    CO=Concentrations[2]
    CO2=Concentrations[3]
    H2O=Concentrations[4]
    N2=Concentrations[5]
    
    ka=Ae*np.exp(-Ea/T)*np.sign(C)*abs(C)**me*O2**ne
#    kbf=10**14.6*exp(-20131/T)*CO*H2O**0.5*O2**0.25
#    kbr=5*10**8*exp(-20131/T)*CO2
    kbf=3.98*10**8*np.exp(-10/(1.987*10**-3)/T)*CO*H2O**0.5*O2**0.25
    kbr=6.16*10**13*np.sign(T)*abs(T)**-0.97*np.exp(-78.4/(1.987*10**-3)/T)*CO2*H2O**0.5*O2**-0.25
    
    dC=-ka/N
    dO2=(-(c/2+h/4)*ka-0.5*c*kbf+c*0.5*kbr)/N
    dCO=(-kbf+ka+kbr)*c/N
    dCO2=(kbf-kbr)*c/N
    dH2O=h/2*ka/N
    dN2=0
    
    return dC,dO2,dCO,dCO2,dH2O,dN2 #gmol;cm3


def twoStepKinects (Concentrations0,ts,Temperature,Pressure):
    concentration = solve_ivp(dc,ts,Concentrations0,args=(Temperature,)).y
        
    C = concentration[0][-1]
    O2=concentration[1][-1]
    CO = concentration[2][-1]
    CO2=concentration[3][-1]
    H2O = concentration[4][-1]
    N2 = concentration[5][-1]
    
    next_concentration=[C,O2,CO,CO2,H2O,N2]
    
    return next_concentration

def motorComb(t,x,xt,dxdt):
    P=x[0]
    T=x[1]
    Qa=x[2]
    Qp=x[3]
    
    A=area(t)
    V=volume(t)
    dVdt=dvdt(t)
    
    k = KReag(T)*(1-xt)+KProd(T)*xt
    dkdT=xt*dkprod(T)+(1-xt)*dkreag(T)

    dWdt=P*dVdt  #joule

    dTdt=(((1/(P*V))*(Qtot*dxdt)-dVdt/V)*(k-1))/((1/T)-(1/(k-1))*dkdT)
    dPdt=((((Qtot*dxdt)-P*dVdt)*(k-1))-P*dVdt+(P*V*dkdT*dTdt/(k-1)))/V
    dQadt=P*dVdt+(1/(k-1))*(V*dPdt+P*dVdt-(P*V*dkdT*dTdt/(k-1)))
    
    return(dPdt,dTdt,dQadt,dWdt)   

time = np.linspace(radians(-164),radians(146),1000)
xt,dxdt=0,0
massFraction=[]
V0=volume(teta0)
T0=P0*10**-3*V0/(m_m*Rg) 

x0 = [P0,T0,0,0]
P=np.ones(len(time))*x0[0]
T=np.ones(len(time))*x0[1]
Qa=np.ones(len(time))*x0[2]
W=np.ones(len(time))*x0[3]
  


c0fuel=P0*X0/(Ru*T0)
c0O2=P0*X0*OC/(Ru*T0)
c0N2=P0*X0*OC*3.76/(Ru*T0)
concentration0 = [c0fuel*10**-6,c0O2*10**-6,0,0,0,c0N2*10**-6] #kmol/cm3
tig=find_nearest(time,radians(-4))
for t in range(len(time)-1):
    ts=[time[t],time[t+1]]
    
    x=solve_ivp(motorComb,ts,x0,args=(xt,dxdt),method='DOP853').y #T(t+1),P(t+1) com x(t)
    next_c = twoStepKinects(concentration0,ts,T[t],P[t]) #x(t+1)
    
    P[t+1]=x[0][-1]
    T[t+1]=x[1][-1]
    Qa[t+1]=x[2][-1]
    W[t+1]=x[3][-1]
    
    dcdt = -dc(ts[-1],next_c,T[t+1])[0]
    #molar_fractions = np.asarray(next_c)*Ru*T[t+1]*10**6/P[t+1]
    molar = np.asarray(next_c)*10**3*volume(time[t+1])
    
    xt=(m_c-molar[0]*M_c)/m_c
    dxdt = dcdt*10**3*volume(time[t+1])*M_c/m_c
    concentration0=next_c
    
#        c0fuel=mols_c/volume(time[t])
#        c0O2=mols_c/volume(time[t])
#        c0N2=mols_c/volume(time[t])
#        concentration0 = [c0fuel*10**-6,c0O2*10**-6,0,0,0,c0N2*10**-6] #kmol/cm3
#    M_m=0
#    for M,ys in zip(Molecular_Mass_Kinects,molar_fractions):
#        M_m += M*ys
#        
#    mass_fraction=molar_fractions[0]*(M_c/M_m)
#    xt=(Y0-mass_fraction)/Y0
#    dxdt = (dcdt*Ru*T[t+1]/P[t+1])*(M_c/M_m)*10**6/Y0
    print(xt,dxdt)
    
    
    x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1]]
    massFraction.append(xt)
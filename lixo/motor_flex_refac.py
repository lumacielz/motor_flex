#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 08:53:35 2022

@author: luiza.maciel
"""


# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 17:19:20 2021

@author: lumac
"""

from math import radians,degrees,sin,cos,sqrt,log,exp
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp,odeint
from scipy import signal

print("\n------------------------------ANALISE DE MOTOR FLEX------------------------------\n")
#dados_motor=input('Dados do Motor: ').split(',')
dados_motor=11,0.144,0.081,0.0864,4

r=float(dados_motor[0])
L=float(dados_motor[1])
D=float(dados_motor[2])
S=float(dados_motor[3])
cilindros=float(dados_motor[4])


#avanco=float(input("Ângulo de avanço: "))
avanco=21
if avanco<=45 and avanco>36:
    fcor=(avanco-36)*((1.77-1.67)/(45-36))+1.67
elif avanco<=36 and avanco >27:
    fcor=(avanco-27)*((1.67-1.68)/(36-27))+1.68
elif avanco<=27 and avanco >=18:
    fcor=(avanco-18)*((1.68-1.45)/(27-18))+1.45
print('fcor= '+ str(fcor))

#teta0d=float(input("ângulo de fechamento da válvula de admissão [em graus]: "))
#teta_combd=float(input("ângulo de início da combustão[em graus]: "))
#delta_tetad=float(input("duração da combustão[em graus]: "))
#tetafd=float(input("ângulo de abertura da válvula de escape[em graus] "))
teta0d=-164
teta_combd=-4
delta_tetad=37.6
tetafd=146
teta0=radians(teta0d)
tetaf=radians(tetafd)
teta_comb=radians(teta_combd)
delta_teta=radians(delta_tetad)


Tp=373
Var=89.25
P0=67.49*10**3
rotacao=2492.5

combustivel="gasolina"

while 1>0:
    #combustivel=input("Com qual combustível deseja operar?  ")
    if combustivel=="gasolina":
        AC=9.6 #mols de oxigenio 
        ACm=12.99 #massica
        c=6.67
        hid=12.8
        ox=0.533
        concentracaoCO2=6.67
        concentracaoH2O=6.4
        concentracaoN2=3.76*AC
        PCI=39249*10**3 #kJ/kg
        a0=-8102.9675849653
        a1=7963.204888950
        a2=-2923.24238668569
        a3=509.144026930886
        a4=-42.2764373650608
        a5=1.35175803975684
        M_c=(12.0107*c+hid*1.008+16*ox)
        M_m=(M_c+AC*(32+3.76*28.013))/(1+AC+AC*3.76)  #massa molecular da mistura
        #cinetica
        Ae=(5.7*10**11) #octano
        Ea=15098
        me=0.25
        ne=1.5
        break
    elif combustivel=="alcool":
        AC=3.2
        ACm=8.427
        c=2.15
        hid=6.62
        ox=1.23
        concentracaoCO2=2.15
        concentracaoH2O=3.31
        concentracaoN2=AC*3.76
        PCI=24804*10**3
        a0=-12482.8740213179
        a1=10263.4623453335
        a2=-3316.7850402598
        a3=526.309291795851
        a4=-40.8869367350809
        a5=1.24555084441151
        mol=((12*c+hid+16*ox)+AC*(32+3.76*28))/(1+AC+AC*3.76)
        break
    elif combustivel=="GNV":
        AC=2.112 #molar
        ACm=16.64
        met=0.92285
        #92.285 em volume
        et=5.455*10**-2
        prop=1.115*10**-2
        but=0.125*10**-2
        co2=0.255*10**-2
        r=0.765*10**-2
        concentracaoH2O=2.0722
        concentracaoCO2=1.0789
        concentracaoN2=AC*3.76+r
        PCI=48737*10**3
        a0=-12713.9307245771
        a1=10272.7734235011
        a2=-3251.2253041433
        a3=504.476942937525
        a4=-38.3629908265519
        a5=1.14655193729902
        mol=((met*16+et*30+prop*44+but*58+co2*44+28*r)+AC*(32+3.76*28))/(1+AC+AC*3.76)
        break
    else:
        print("digite um combustível válido!!!")
        
#coefwiebe=input('Coeficientes da equação de Wiebe m , a: ').split(',')
coefwiebe=1,2
m=1
a=2
M = [M_c, 32, 12.0107+16, 12.0107+32, 2*1.008*16]
M_a = 28.962

def KReag(Ts):
    Cp=a0+a1*log(Ts)+a2*log(Ts)**2+a3*log(Ts)**3+a4*log(Ts)**4+a5*log(Ts)**5
    CpO=10228.3426-7184.92333*log(Ts)+2010.86808*log(Ts)**2-279.69496*log(Ts)**3+19.34823*log(Ts)**4-0.53257*log(Ts)**5
    CpN=0.54497*log(Ts)**5-18.69984*log(Ts)**4+254.29554*log(Ts)**3-1712.17390*log(Ts)**2+5708.38047*log(Ts)-7513.3642
    KReag=((Cp/(Cp-Ru))+AC*(CpO/(CpO-Ru))+3.76*AC*(CpN/(CpN-Ru)))/(1+AC+3.76*AC)
    return KReag
    
def KProd(Ts):
    CpN=0.54497*log(Ts)**5-18.69984*log(Ts)**4+254.29554*log(Ts)**3-1712.17390*log(Ts)**2+5708.38047*log(Ts)-7513.3642
    CpCO2=0.20754*log(Ts)**5-6.43522*log(Ts)**4+77.54809*log(Ts)**3-452.81197*log(Ts)**2+1288.4677*log(Ts)-1412.36785
    CpH2O=0.64541*log(Ts)**5-23.54277*log(Ts)**4+339.33662*log(Ts)**3-2414.77575*log(Ts)**2+8490.5218*log(Ts)-11780.76495
    KProd=((CpCO2/(CpCO2-Ru))*concentracaoCO2+(CpH2O/(CpH2O-Ru))*concentracaoH2O+(CpN/(CpN-Ru))*AC*3.76)/(concentracaoCO2+concentracaoH2O+concentracaoN2)
    return KProd

def dkreag(Ts):
    dkr=(6.66068603173739e-8*AC*(-2.72485*log(Ts)**4/Ts + 74.79936*log(Ts)**3/Ts - 762.88662*log(Ts)**2/Ts + 3424.3478*log(Ts)/Ts - 5708.38047/Ts)*(0.54497*log(Ts)**5 - 18.69984*log(Ts)**4 + 254.29554*log(Ts)**3 - 1712.1739*log(Ts)**2 + 5708.38047*log(Ts) - 7513.3642)/(-0.000133096170155042*Ru + 7.25334198493932e-5*log(Ts)**5 - 0.00248887708651206*log(Ts)**4 + 0.0338457624615083*log(Ts)**3 - 0.227883788729422*log(Ts)**2 + 0.759763578344838*log(Ts) - 1)**2 + AC*(-2.66285*log(Ts)**4/Ts + 77.39292*log(Ts)**3/Ts - 839.08488*log(Ts)**2/Ts + 4021.73616*log(Ts)/Ts - 7184.92333/Ts)/(-Ru - 0.53257*log(Ts)**5 + 19.34823*log(Ts)**4 - 279.69496*log(Ts)**3 + 2010.86808*log(Ts)**2 - 7184.92333*log(Ts) + 10228.3426) + 9.55849389871466e-9*AC*(2.66285*log(Ts)**4/Ts - 77.39292*log(Ts)**3/Ts + 839.08488*log(Ts)**2/Ts - 4021.73616*log(Ts)/Ts + 7184.92333/Ts)*(-0.53257*log(Ts)**5 + 19.34823*log(Ts)**4 - 279.69496*log(Ts)**3 + 2010.86808*log(Ts)**2 - 7184.92333*log(Ts) + 10228.3426)/(-9.77675503360632e-5*Ru - 5.20680642824772e-5*log(Ts)**5 + 0.00189162905043873*log(Ts)**4 - 0.0273450910805432*log(Ts)**3 + 0.196597646230583*log(Ts)**2 - 0.70245235332653*log(Ts) + 1)**2 + 3.76*AC*(2.72485*log(Ts)**4/Ts - 74.79936*log(Ts)**3/Ts + 762.88662*log(Ts)**2/Ts - 3424.3478*log(Ts)/Ts + 5708.38047/Ts)/(-Ru + 0.54497*log(Ts)**5 - 18.69984*log(Ts)**4 + 254.29554*log(Ts)**3 - 1712.1739*log(Ts)**2 + 5708.38047*log(Ts) - 7513.3642) + (-a1/Ts - 2*a2*log(Ts)/Ts - 3*a3*log(Ts)**2/Ts - 4*a4*log(Ts)**3/Ts - 5*a5*log(Ts)**4/Ts)*(a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)/(-Ru + a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5)**2 + (a1/Ts + 2*a2*log(Ts)/Ts + 3*a3*log(Ts)**2/Ts + 4*a4*log(Ts)**3/Ts + 5*a5*log(Ts)**4/Ts)/(-Ru + a0 + a1*log(Ts) + a2*log(Ts)**2 + a3*log(Ts)**3 + a4*log(Ts)**4 + a5*log(Ts)**5))/(4.76*AC + 1)
    return dkr

def dkprod(Ts):
    dkp=(6.66068603173739e-8*AC*(-2.72485*log(Ts)**4/Ts + 74.79936*log(Ts)**3/Ts - 762.88662*log(Ts)**2/Ts + 3424.3478*log(Ts)/Ts - 5708.38047/Ts)*(0.54497*log(Ts)**5 - 18.69984*log(Ts)**4 + 254.29554*log(Ts)**3 - 1712.1739*log(Ts)**2 + 5708.38047*log(Ts) - 7513.3642)/(-0.000133096170155042*Ru + 7.25334198493932e-5*log(Ts)**5 - 0.00248887708651206*log(Ts)**4 + 0.0338457624615083*log(Ts)**3 - 0.227883788729422*log(Ts)**2 + 0.759763578344838*log(Ts) - 1)**2 + 3.76*AC*(2.72485*log(Ts)**4/Ts - 74.79936*log(Ts)**3/Ts + 762.88662*log(Ts)**2/Ts - 3424.3478*log(Ts)/Ts + 5708.38047/Ts)/(-Ru + 0.54497*log(Ts)**5 - 18.69984*log(Ts)**4 + 254.29554*log(Ts)**3 - 1712.1739*log(Ts)**2 + 5708.38047*log(Ts) - 7513.3642) + 5.01307675179101e-7*concentracaoCO2*(-1.0377*log(Ts)**4/Ts + 25.74088*log(Ts)**3/Ts - 232.64427*log(Ts)**2/Ts + 905.62394*log(Ts)/Ts - 1288.4677/Ts)*(0.20754*log(Ts)**5 - 6.43522*log(Ts)**4 + 77.54809*log(Ts)**3 - 452.81197*log(Ts)**2 + 1288.4677*log(Ts) - 1412.36785)/(-0.000708030843381205*Ru + 0.000146944721235335*log(Ts)**5 - 0.0045563342439436*log(Ts)**4 + 0.0549064395653016*log(Ts)**3 - 0.320604841012205*log(Ts)**2 + 0.912274872300442*log(Ts) - 1)**2 + concentracaoCO2*(1.0377*log(Ts)**4/Ts - 25.74088*log(Ts)**3/Ts + 232.64427*log(Ts)**2/Ts - 905.62394*log(Ts)/Ts + 1288.4677/Ts)/(-Ru + 0.20754*log(Ts)**5 - 6.43522*log(Ts)**4 + 77.54809*log(Ts)**3 - 452.81197*log(Ts)**2 + 1288.4677*log(Ts) - 1412.36785) + 7.20531576341265e-9*concentracaoH2O*(-3.22705*log(Ts)**4/Ts + 94.17108*log(Ts)**3/Ts - 1018.00986*log(Ts)**2/Ts + 4829.5515*log(Ts)/Ts - 8490.5218/Ts)*(0.64541*log(Ts)**5 - 23.54277*log(Ts)**4 + 339.33662*log(Ts)**3 - 2414.77575*log(Ts)**2 + 8490.5218*log(Ts) - 11780.76495)/(-8.48841313992942e-5*Ru + 5.47850672464185e-5*log(Ts)**5 - 0.00199840758218336*log(Ts)**4 + 0.0288042942406724*log(Ts)**3 - 0.204976142062829*log(Ts)**2 + 0.720710568119772*log(Ts) - 1)**2 + concentracaoH2O*(3.22705*log(Ts)**4/Ts - 94.17108*log(Ts)**3/Ts + 1018.00986*log(Ts)**2/Ts - 4829.5515*log(Ts)/Ts + 8490.5218/Ts)/(-Ru + 0.64541*log(Ts)**5 - 23.54277*log(Ts)**4 + 339.33662*log(Ts)**3 - 2414.77575*log(Ts)**2 + 8490.5218*log(Ts) - 11780.76495))/(concentracaoCO2 + concentracaoH2O + concentracaoN2)
    return dkp

def wiebe(teta):
    return 1-exp(-a*((teta-teta_comb)/delta_teta)**(m+1))

#def dx(teta):
#    return -exp(-a*((teta-teta_comb)/delta_teta)**(m+1))*(-a*(m+1)*((teta-teta_comb)/delta_teta)**m)*1/delta_teta

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
    dO2=(-(c/2+hid/4)*ka-0.5*c*kbf+c*0.5*kbr)/N
    dCO=(-kbf+ka+kbr)*c/N
    dCO2=(kbf-kbr)*c/N
    dH2O=hid/2*ka/N
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
    

#ti=np.linspace(teta0,teta_comb,100)[:-1]
#tj=np.linspace(teta_comb,teta_comb+delta_teta,100)[:-1]
#tf=np.linspace(teta_comb+delta_teta,tetaf,100)
tk=np.linspace(teta0,tetaf,1000)
#tk[:len(ti)]=ti
#tk[len(ti):len(ti)+len(tj)] = tj
#tk[len(tj)+len(ti):]=tf
ti=tk[:find_nearest(tk,teta_comb)]
tj=tk[find_nearest(tk,teta_comb):find_nearest(tk,teta_comb+delta_teta)]
tf=tk[find_nearest(tk,teta_comb+delta_teta):]

R=S/2 #raio do virabrequim
N=rotacao*2*np.pi/60 #rad/s

m_ar=Var/(0.5*3600*cilindros*(N/(2*np.pi))) #rot por s
m_m=m_ar/ACm+m_ar #massa da mistura admitida
m_c=(m_m/(1+ACm)) 
print(m_c)
 
Y0 = m_c/(m_m)
   
ef=0.87
Qtot=ef*m_c*PCI #kJ
Vd=np.pi*D**2*R/2 #volume deslocado m^3
vp=2*S*N #velocidade do pistão m/s

Ru=8.314 #kJ/kmolK
V0=volume(teta0)
Rg=Ru/M_m #constante dos gases
T0=P0*10**-3*V0/(m_m*Rg)  
Qa0=0
Qp0=0
W0=0


x0=[P0,T0,Qa0,Qp0,W0] #condicoes iniciais
print(x0)

P=np.ones(len(tk))*P0
T=np.ones(len(tk))*T0
Qa=np.ones(len(tk))*Qa0
Qp=np.ones(len(tk))*Qp0
W=np.ones(len(tk))*W0
P1=np.ones(len(tk))*P0
T1=np.ones(len(tk))*T0

Spark = signal.unit_impulse(len(tk),find_nearest(tk,radians(-avanco)))*10000

def motored(t,x):
    P=x[0]
    T=x[1]
    Qa=x[2]
    Qp=x[3]
    W=x[4]
    
    A=area(t)
    V=volume(t)
    dVdt=dvdt(t)
    
    k=KReag(T)
    dkdT=dkreag(T)
    
    dxdt=0
    h=0
    
    dWdt=P*dVdt  #joule
    dQpdt=(h*A*(T-Tp)/(N*fcor)) #rads por s
    dTdt=(((1/(P*V))*(Qtot*dxdt-fcor*dQpdt)-dVdt/V)*(k-1))/((1/T)-(1/(k-1))*dkdT)
    dPdt=((((Qtot*dxdt-fcor*dQpdt)-P*dVdt)*(k-1))-P*dVdt+(P*V*dkdT*dTdt/(k-1)))/V
    dQadt=P*dVdt+(1/(k-1))*(V*dPdt+P*dVdt-(P*V*dkdT*dTdt/(k-1)))
    
    return(dPdt,dTdt,dQadt,dQpdt,dWdt)  

    
def motor(t,x):
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
    dTdt=(((1/(P*V))*(Qtot*dxdt-fcor*dQpdt+spark)-dVdt/V)*(k-1))/((1/T)-(1/(k-1))*dkdT)
    dPdt=((((Qtot*dxdt-fcor*dQpdt)-P*dVdt+spark)*(k-1))-P*dVdt+(P*V*dkdT*dTdt/(k-1)))/V
    dQadt=P*dVdt+(1/(k-1))*(V*dPdt+P*dVdt-(P*V*dkdT*dTdt/(k-1)))
    
    return(dPdt,dTdt,dQadt,dQpdt,dWdt)      
    

#curvaa de pressao sem combustao
for l in range(len(tk)-1):
    ts=[tk[l],tk[l+1]]
    spark = Spark[l]
    x=solve_ivp(motored,ts,x0,method='DOP853').y
    
    P1[l+1]=x[0][-1]
    T1[l+1]=x[1][-1]
    
    x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]

x0=[P0,T0,Qa0,Qp0,W0]
#fechamento da valvula ao inicio da combustao
for i in range(len(ti)-1): #range(0,98)
    ts=[ti[i],ti[i+1]]
    spark = Spark[i]
    
    x=solve_ivp(motor,ts,x0,method='DOP853').y
    
    P[i+1]=x[0][-1] #ultima linha,coluna1
    T[i+1]=x[1][-1]
    Qa[i+1]=x[2][-1]
    Qp[i+1]=x[3][-1]
    W[i+1]=x[4][-1]
    
    x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]

#concentracao inicial
ro_fuel = P[i+1]* 10**-3 * M_c /(Ru* T[i+1])
ro_o2 = P[i+1]* 10**-3 * 32 /(Ru * T[i+1])
ro =  P[i+1]* 10**-3 * M_m /(Ru* T[i+1])

c0fuel=(m_c/m_m)*ro*10**-3/M_c
c0O2=((m_ar*0.232)/m_m)*ro*10**-3/32
concentration0 = [c0fuel,c0O2,0,0,0]

xt=(Y0 - c0fuel * M_c *10**3/ro)/Y0
print(xt)
wiebeEv=[]
print(concentration0)

dxdt=(-dc(concentration0,ts[-1],T[i+1])[0]*M_c*10**3 / ro)/Y0

for j in range (len(ti)-1,len(ti)+len(tj)-1):
    ts=[tk[j],tk[j+1]]
    spark = Spark[j]
    x=solve_ivp(motorComb,ts,x0,args=(P1[j],xt,dxdt),method='DOP853').y
    
    P[j+1]=x[0][-1]
    T[j+1]=x[1][-1]
    Qa[j+1]=x[2][-1]
    Qp[j+1]=x[3][-1]
    W[j+1]=x[4][-1]
    
    xt, dxdt, concentration = twoStepKinects(concentration0,ts,T[j+1],P[j+1])
    concentration0=concentration
    print(xt,dxdt)
    x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
   
    wiebeEv.append(xt)

for f in range(len(tj)+len(ti)-1,len(tk)-1): #200 ao 298
    ts=[tk[f],tk[f+1]]
    spark=Spark[f]
    x=solve_ivp(motorComb,ts,x0,args=(P1[f],xt,dxdt),method='DOP853').y
    
    P[f+1]=x[0][-1]
    T[f+1]=x[1][-1]
    Qa[f+1]=x[2][-1]
    Qp[f+1]=x[3][-1]
    W[f+1]=x[4][-1]
    
    xt, dxdt, concentration = twoStepKinects(concentration0,ts,T[f+1],P[f+1])
    concentration0=concentration
    
    x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
    

    wiebeEv.append(xt)


tkd=list(map(lambda x:degrees(x),tk))
#print(P,T,Qa,Qp,W)
print('\n'+'Pmax: '+ str(max(P)*10**-6)+' MPa' )
print('Tmax: '+str(max(T))+' K' +str(degrees(tk[list(T).index(max(T))])))

fig,(ax1,ax2)=plt.subplots(2,3)
ax1[0].plot(tkd,P*10**(-6))
ax1[1].plot(tkd,T)
ax1[2].plot(tk[len(ti):],wiebeEv)
ax2[0].plot(tkd,Qa)
ax2[1].plot(tkd,Qp)
ax2[2].plot(tkd,W)
ax1[0].set_ylabel('P[MPa]')
ax1[1].set_ylabel('T[K]')
ax1[2].set_ylabel('X(teta)')
ax2[0].set_ylabel('Qa[J]')
ax2[1].set_ylabel('Qp[J]')
ax2[2].set_ylabel('W[J]')
plt.show()
#path='./dados.txt'
#with open(path,'w') as f:
#    f.write("fuel mass per cicle\n"+str(m_c)+"\n")
#    f.write("Temperatures:\n")
#    for d in range(len(T)):
#        f.write(str(T[d])+'\n')
#    f.write("Pressures:\n")
#    for d in range(len(P)):
#        f.writelines(str(P[d])+'\n')
#    f.close()

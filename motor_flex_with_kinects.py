# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 17:19:20 2021

@author: lumac
"""

from data import *
import sympy as sp
from scipy.integrate import odeint,solve_ivp

teta=sp.symbols("teta")
Ts=sp.symbols('Ts')

s=R*sp.cos(teta)+sp.sqrt(L**2-R**2*(sp.sin(teta))**2)  #deslocamento
Ap=sp.pi*D*(D/2+L+R-s+(2*R/(r-1)))  #area das paredes
Vo=(sp.pi*D**2/4)*(L+R-s+2*R/(r-1)) #volume 

A=Ap.subs(teta,teta0)
V0=Vo.subs(teta,teta0)
V=V0
dVdt=sp.diff(Vo,teta).subs(teta,teta0)

T0=float(P0*10**-3*V0/(m_m*Rg))  
Qa0=0
Qp0=0
W0=0
x0=[float(P0),float(T0),Qa0,Qp0,W0] #condicoes iniciais
print(x0)

ti=np.linspace(teta0,teta_comb,100) #fechamento da admissao ate inicio da combustao
tj=np.linspace(teta_comb,teta_comb+delta_teta,100) #combustao
tf=np.linspace(teta_comb+delta_teta,tetaf,100) #fim da combustao ate abertura do escapamento
tk=np.ones(len(ti)+len(tj)+len(tf)) 
tk[:len(ti)]=ti #0 a 99
tk[len(ti):len(tj)+len(ti)]=tj #100 ao 199
tk[len(tj)+len(ti):len(tf)+len(tj)+len(ti)]=tf #200 ao 299

P=np.ones(len(tk))*P0
T=np.ones(len(tk))*T0
Qa=np.ones(len(tk))*Qa0
Qp=np.ones(len(tk))*Qp0
W=np.ones(len(tk))*W0
P1=np.ones(len(tk))*P0
T1=np.ones(len(tk))*T0

x0=[float(P0),float(T0),Qa0,Qp0,W0]

def motor(t,x):
    P=x[0]
    T=x[1]
    Qa=x[2]
    Qp=x[3]
    W=x[4]
    dWdt=P*dVdt  #joule
    dQpdt=(h*A*(T-Tp)/(N*fcor)) #rads por s
    dTdt=((((1/(P*V))*(Qtot*dxdt-fcor*dQpdt)-dVdt/V)*(k-1))+dkdt/(k-1))*T
    dPdt=((((Qtot*dxdt-fcor*dQpdt)-P*dVdt)*(k-1))-P*dVdt+(P*V*dkdt/(k-1)))/V
    dQadt=P*dVdt+(1/(k-1))*(V*dPdt+P*dVdt-(P*V*dkdt/(k-1)))
    return(dPdt,dTdt,dQadt,dQpdt,dWdt)  

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

#def dc (C,t,Temp):
#    O2=12.5*C
#    dC=-(Ae*exp(-Ea/Temp)*np.sign(C)*abs(C)**me*(O2)**ne)/N
#    return dC
#
#def oneStep(C0,ts,Temperature,Pressure):
#    x=odeint(dc,C0,ts,args=(Temperature,))
#    next_c=x[-1]
#    mass_burned=(m0-(x[-1][0]*8315*Temperature/(Pressure*10**-3))*m_c)/m0
#    return next_c,mass_burned

def twoStepKinects (Concentrations0,ts,Temperature,Pressure):
    x=odeint(dc,Concentrations0,ts,args=(Temperature,),atol=10**-17,rtol=10**-17)
    mass_burned=(m0-(x[-1][0]*8315*Temperature/(Pressure*10**-3))*m_cm)/m0
    next_t=x[-1]
    return next_t,mass_burned

K=KReag
k=K.subs(Ts,T0)
dkdt=sp.diff(K,teta).subs([(Ts,T0),(teta,teta0)])

dxdt=0
h=0

#curva sem combustao
for l in range(len(tk)-1):
    ts=[tk[l],tk[l+1]]
    x=solve_ivp(motor,ts,x0,rtol=1e-10,atol=1e-10).y
    P1[l+1]=x[0][-1]
    T1[l+1]=x[1][-1]
    x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]

    V=Vo.subs(teta,tk[l+1])
    A=Ap.subs(teta,tk[l+1])
    dVdt=sp.diff(Vo,teta).subs(teta,ts[-1])
    k=K.subs(Ts,T1[l+1])
    dkdt=sp.diff(K,teta).subs([(Ts,T1[l+1]),(teta,ts[-1])])

#fechamento da admissao ao inicio da combustao 
for i in range(len(ti)-1): 
    ts=[ti[i],ti[i+1]]
    x=solve_ivp(motor,ts,x0,rtol=1e-10,atol=1e-10).y
    P[i+1]=x[0][-1] #ultima linha,coluna1
    T[i+1]=x[1][-1]
    Qa[i+1]=x[2][-1]
    Qp[i+1]=x[3][-1]
    W[i+1]=x[4][-1]
    
    x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
    
    V=Vo.subs(teta,ti[i+1])
    A=Ap.subs(teta,ti[i+1])
    dVdt=sp.diff(Vo,teta).subs(teta,ts[-1])
    
    k=K.subs(Ts,T[i+1])
    dkdt=sp.diff(K,teta).subs([(Ts,T[i+1]),(teta,ts[-1])])
    
    vg=2.28*vp
    h=0.013*D**(-0.2)*(P[i+1])**0.8*T[i+1]**(-0.53)*vg**0.8
    
    print(k,dkdt)

wiebe=[]
C0=(1/((1+AC)+3.76*AC))*(P[99]*10**-3/(8315*T[99]))
Concentrations0=[C0,AC*C0,0,0,0]

dTdt=motor(ts[1],x0)[1]
dPdt=motor(ts[1],x0)[0]
dCdt=dc(Concentrations0,ts[1],T[99])

#x_teta=(m0-(C0*8315*T[99]/(P[99]*10**-3))*m_cm)/m0
#dxdt=-dCdt*8315*T[99]*m_cm/(P[99]*10**-3*m0)
#dxdt=-(m_cm*Rg/m0)*(P[99]*(dCdt*T[99]+C*dTdt)-dPdt*C*T[99])/((P[99]*10**-3)**2)
print(float(dxdt))

#durante a combustao
for j in range(len(tk))[len(ti)-1:len(ti)+len(tj)-1]: 
    ts=[tk[j],tk[j+1]]
    x=solve_ivp(motor,ts,x0,rtol=1e-10,atol=1e-10).y
    P[j+1]=x[0][-1]
    T[j+1]=x[1][-1]
    Qa[j+1]=x[2][-1]
    Qp[j+1]=x[3][-1]  
    W[j+1]=x[4][-1]
    x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
    
    V=Vo.subs(teta,tk[j+1])
    A=Ap.subs(teta,tk[j+1])
    dVdt=sp.diff(Vo,teta).subs(teta,ts[-1])
    
    K=KReag+(KProd-KReag)*(teta-teta0)/delta_teta
    
    k=K.subs([(teta,ts[-1]),(Ts,T[j+1])])
    dkdt=sp.diff(K,teta).subs([(Ts,T[j+1]),(teta,ts[-1])])
    
    next_t,x_teta=twoStepKinects(Concentrations0,ts,T[j+1],P[j+1]) #composicao em t[j+1]
    
    dCdt=dc(next_c,ts[-1],T[j+1])[0]
    Concentrations0=next_t
    MolarFractions=list(map(lambda x:x*8315*T[j+1]/(P[j+1]*10**-3),next_t))
    print(MolarFractions[0])

    dxdt=-dCdt*8315*T[j+1]*m_cm/(P[j+1]*10**-3*m0)

    print(float(next_t[0]))
  
    wiebe.append(x_teta)
    
    vg=2.28*vp+0.00324*(P[j+1]-P1[j+1])*Vd*T0/(P0*V0)
    h=0.013*D**-0.2*(P[j+1])**0.8*T[j+1]**-0.53*vg**0.8
      
K=KProd
dxdt=0

#fim da combustao ate abertura do escape
for f in range(len(tk))[len(tj)+len(ti)-1:len(tk)-1]: #200 ao 298
    ts=[tk[f],tk[f+1]]
    x=solve_ivp(motor,ts,x0,rtol=1e-10,atol=1e-10).y
    P[f+1]=x[0][-1]
    T[f+1]=x[1][-1]
    Qa[f+1]=x[2][-1]
    Qp[f+1]=x[3][-1]
    W[f+1]=x[4][-1]
    x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
    
    V=Vo.subs(teta,tk[f+1])
    A=Ap.subs(teta,tk[f+1])
    dVdt=sp.diff(Vo,teta).subs(teta,ts[-1])
    
    k=K.subs(Ts,T[f+1])
    dkdt=sp.diff(K,teta).subs([(Ts,T[f+1]),(teta,ts[-1])])
    
    vg=2.28*vp+0.00324*(P[f+1]-P1[f+1])*Vd*T0/(P0*V0)
    h=0.013*D**-0.2*(P[f+1])**0.8*T[f+1]**-0.53*vg**0.8

#plots 
#print(P,T,Qa,Qp,W)
print('\n'+'Pmax: '+ str(max(P)*10**-6)+' MPa' )
print('Tmax: '+str(max(T))+' K')
print(W[199]-W[100])
print(Qa[199]-Qa[100])
print(Qp[199]-Qp[100])
fig,(ax1,ax2)=plt.subplots(2,3)
ax1[0].plot(tk,P*10**(-6))
ax1[1].plot(tk,T)
ax1[2].plot(tj,wiebe)
ax2[0].plot(tk,Qa)
ax2[1].plot(tk,Qp)
ax2[2].plot(tk,W)
ax1[0].set_ylabel('P[MPa]')
ax1[1].set_ylabel('T[K]')
ax1[2].set_ylabel('X(teta)')
ax2[0].set_ylabel('Qa[J]')
ax2[1].set_ylabel('Qp[J]')
ax2[2].set_ylabel('W[J]')
plt.show()
path='./dados.txt'
with open(path,'w') as f:
    f.write("fuel mass per cicle\n"+str(m_c)+"\n")
    f.write("Temperatures:\n")
    for d in range(len(T)):
        f.write(str(T[d])+'\n')
    f.write("Pressures:\n")
    for d in range(len(P)):
        f.writelines(str(P[d])+'\n')
    f.close()

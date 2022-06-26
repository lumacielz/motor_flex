#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 20:35:18 2021

@author: luiza.maciel
"""


from data import *
import sympy as sp
from scipy.integrate import odeint,solve_ivp
from scipy import signal

ti=np.linspace(teta0,teta_comb,100) #fechamento da admissao ate inicio da combustao
tj=np.linspace(teta_comb,teta_comb+delta_teta,100) #combustao
tf=np.linspace(teta_comb+delta_teta,tetaf,100) #fim da combustao ate abertura do escapamento
tk=np.ones(len(ti)+len(tj)+len(tf)) 
tk[:len(ti)]=ti #0 a 99
tk[len(ti):len(tj)+len(ti)]=tj #100 ao 199
tk[len(tj)+len(ti):len(tf)+len(tj)+len(ti)]=tf #200 ao 299

aprox=[abs(i-radians(-avanco))for i in ti]
disp=aprox.index(min(aprox))

#Spark=signal.unit_impulse(len(tk),disp)*5000

V0=Vo.subs(teta,teta0)
T0=float(P0*10**-3*V0/(m_m*Rg))  

P=np.ones(len(tk))*P0
T=np.ones(len(tk))*T0
Qa=np.ones(len(tk))*Qa0
Qp=np.ones(len(tk))*Qp0
W=np.ones(len(tk))*W0
P1=np.ones(len(tk))*P0
T1=np.ones(len(tk))*T0

x0=[float(P0),float(T0),Qa0,Qp0,W0]

x_teta=1-sp.exp(-a*((teta-teta_comb)/delta_teta)**(m+1))

wiebe=[]
def motor(t,x):
    P=x[0]
    T=x[1]
    Qa=x[2]
    Qp=x[3]
    W=x[4]
    
    
    A=Ap.subs(teta,t)  #area das paredes
    V=Vo.subs(Vo,t)#volume 
    k=K.subs([(Ts,T),(teta,t)])
    print(A,V,k)
    
    
    dkdt=sp.diff(K,teta)
    dkdT=sp.diff(dkdt,Ts).subs(Ts,T).subs(teta,t)
    dVdt=sp.diff(V,teta).subs(teta,t)
    print(dkdT,dVdt)
    
    h=0.013*D**-0.2*(P)**0.8*T**-0.53*vg**0.8
    
    dWdt=P*dVdt  #joule
    dQpdt=(h*A*(T-Tp)/(N*fcor)) #rads por s
    dTdt=(((1/(P*V))*(Qtot*dxdt-fcor*dQpdt)-dVdt/V)*(k-1))/((1/T)-(1/(k-1))*dkdT)
    dPdt=((((Qtot*dxdt-fcor*dQpdt)-P*dVdt)*(k-1))-P*dVdt+(P*V*dkdT*dTdt/(k-1)))/V
    dQadt=P*dVdt+(1/(k-1))*(V*dPdt+P*dVdt-(P*V*dkdT*dTdt/(k-1)))
    return(dPdt,dTdt,dQadt,dQpdt,dWdt)  

K=KReag

for l in range(len(tk)-1):
    ts=[tk[l],tk[l+1]]
    x=solve_ivp(motor,ts,x0,rtol=1e-10,atol=1e-10).y
    P1[l+1]=x[0][-1]
    T1[l+1]=x[1][-1]
    x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
    dxdt=0
    vg=0
    print(P1[l])
    

dxdt=0 
for i in range(len(ti)-1): #range(0,98)
    ts=[ti[i],ti[i+1]]
    x=solve_ivp(motor,ts,x0,rtol=1e-10,atol=1e-10).y
    P[i+1]=x[0][-1] #ultima linha,coluna1
    T[i+1]=x[1][-1]
    Qa[i+1]=x[2][-1]
    Qp[i+1]=x[3][-1]
    W[i+1]=x[4][-1]
    x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
   
    vg=2.28*vp
    print(P[i])

dxdt=sp.diff(x_teta,teta).subs(teta,t)

for j in range(len(tk))[len(ti)-1:len(ti)+len(tj)-1]: #199 ao 198
    ts=[tk[j],tk[j+1]]
    wiebe.append(x_teta.subs(teta,tk[j]))
    x=solve_ivp(motor,ts,x0,rtol=1e-10,atol=1e-10).y
    P[j+1]=x[0][-1]
    T[j+1]=x[1][-1]
    Qa[j+1]=x[2][-1]
    Qp[j+1]=x[3][-1]
    W[j+1]=x[4][-1]
    x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
    
    K=((teta-teta0)*(KProd)+(teta0+delta_teta-teta)*(KReag))/delta_teta
    
    
    vg=2.28*vp+0.00324*(P[j+1]-P1[j+1])*Vd*T0/(P0*V0)
    print(P[j])
      

K=KProd


for f in range(len(tk))[len(tj)+len(ti)-1:len(tk)-1]: #200 ao 298
    ts=[tk[f],tk[f+1]]
    x=solve_ivp(motor,ts,x0,rtol=1e-10,atol=1e-10).y
    P[f+1]=x[0][-1]
    T[f+1]=x[1][-1]
    Qa[f+1]=x[2][-1]
    Qp[f+1]=x[3][-1]
    W[f+1]=x[4][-1]
    x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
    
    vg=2.28*vp+0.00324*(P[f+1]-P1[f+1])*Vd*T0/(P0*V0)
    print(P[f])

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
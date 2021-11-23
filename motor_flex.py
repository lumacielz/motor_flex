# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 17:19:20 2021

@author: lumac
"""
__name__=="m"
from math import radians
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

print("\n------------------------------ANALISE DE MOTOR FLEX------------------------------\n")
dados_motor=input('Dados do Motor: ').split(',')
r=float(dados_motor[0])
L=float(dados_motor[1])
D=float(dados_motor[2])
S=float(dados_motor[3])
cilindros=float(dados_motor[4])
R=S/2 #raio do virabrequim
teta=sp.symbols("teta")
s=R*sp.cos(teta)+sp.sqrt(L**2-R**2*(sp.sin(teta))**2)  #deslocamento
Ap=3.14*D*(D/2+L+R-s+(2*R/(r-1)))  #area das paredes
Vo=(3.14*D**2/4)*(L+R-s+2*R/(r-1)) #volume 
print(Vo.subs(teta,0))
avanco=float(input("Ângulo de avanço: "))
if avanco<=45 and avanco>36:
    fcor=(avanco-36)*((1.77-1.67)/(45-36))+1.67
elif avanco<=36 and avanco >27:
    fcor=(avanco-27)*((1.67-1.68)/(36-27))+1.68
elif avanco<=27 and avanco >=18:
    fcor=(avanco-18)*((1.68-1.45)/(27-18))+1.45
print('fcor= '+ str(fcor))
teta0d=float(input("ângulo de fechamento da válvula de admissão [em graus]: "))
teta_combd=float(input("ângulo de início da combustão[em graus]: "))
delta_tetad=float(input("duração da combustão[em graus]: "))
tetafd=float(input("ângulo de abertura da válvula de escape[em graus] "))
teta0=radians(teta0d)
tetaf=radians(tetafd)
teta_comb=radians(teta_combd)
delta_teta=radians(delta_tetad)
A=Ap.subs(teta,teta0)
V0=Vo.subs(teta,teta0)
Var=float(input("Vazão de ar[kg/h]: "))
P0=float(input("Pressão no ângulo de fechamento: "))*10**3
Tp=float(input("Temperatura das paredes: ")) #fazer input
Ru=8.314 #kJ/kmolK
rotacao=float(input("rotação do motor: "))
N=rotacao*2*3.14/60 #rad/s

while 1>0:
    combustivel=input("Com qual combustível deseja operar?  ")
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
        mol=((12*c+hid+16*ox)+AC*(32+3.76*28))/(1+AC+AC*3.76)  #massa molecular da mistura
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
        
coefwiebe=input('Coeficientes da equação de Wiebe m , a: ').split(',')
m=float(coefwiebe[0])
a=float(coefwiebe[1])
m_ar=Var/(0.5*3600*cilindros*(N/(2*3.14))) #rot por s
m_m=m_ar/ACm+m_ar #massa da mistura admitida
m_c=(m_m/(1+ACm)) 
print(m_c)    
Ts=sp.symbols('Ts')
Cp=a0+a1*sp.log(Ts)+a2*sp.log(Ts)**2+a3*sp.log(Ts)**3+a4*sp.log(Ts)**4+a5*sp.log(Ts)**5
CpO=10228.3426-7184.92333*sp.log(Ts)+2010.86808*sp.log(Ts)**2-279.69496*sp.log(Ts)**3+19.34823*sp.log(Ts)**4-0.53257*sp.log(Ts)**5
CpN=0.54497*sp.log(Ts)**5-18.69984*sp.log(Ts)**4+254.29554*sp.log(Ts)**3-1712.17390*sp.log(Ts)**2+5708.38047*sp.log(Ts)-7513.3642
CpCO2=0.20754*sp.log(Ts)**5-6.43522*sp.log(Ts)**4+77.54809*sp.log(Ts)**3-452.81197*sp.log(Ts)**2+1288.4677*sp.log(Ts)-1412.36785
CpH2O=0.64541*sp.log(Ts)**5-23.54277*sp.log(Ts)**4+339.33662*sp.log(Ts)**3-2414.77575*sp.log(Ts)**2+8490.5218*sp.log(Ts)-11780.76495
KReag=((Cp/(Cp-Ru))+AC*(CpO/(CpO-Ru))+3.76*AC*(CpN/(CpN-Ru)))/(1+AC+3.76*AC)
KProd=((CpCO2/(CpCO2-Ru))*concentracaoCO2+(CpH2O/(CpH2O-Ru))*concentracaoH2O+(CpN/(CpN-Ru))*AC*3.76)/(concentracaoCO2+concentracaoH2O+concentracaoN2)

ef=0.87
Qtot=ef*m_c*PCI #kJ
Vd=3.14*D**2*R/2 #volume deslocado m^3
vp=2*S*N #velocidade do pistão m/s
h=0 
Rg=Ru/mol #constante dos gases
T0=P0*10**-3*V0/(m_m*Rg)  
Qa0=0
Qp0=0
W0=0
x0=[P0,T0,Qa0,Qp0,W0] #condicoes iniciais
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
step=(tetaf-teta0)/(len(tk))

def motor(t,x):
    P=x[0]
    T=x[1]
    Qa=x[2]
    Qp=x[3]
    W=x[4]
    dWdt=P*dVdt  #joule
    dQpdt=(h*A*(T-Tp)/(N*fcor)) #rads por s
    dTdt=(((1/(P*V))*(Qtot*dxdt-fcor*dQpdt)-dVdt/V)*(k-1))/((1/T)-(1/(k-1))*dkdT)
    dPdt=((((Qtot*dxdt-fcor*dQpdt)-P*dVdt)*(k-1))-P*dVdt+(P*V*dkdT*dTdt/(k-1)))/V
    dQadt=P*dVdt+(1/(k-1))*(V*dPdt+P*dVdt-(P*V*dkdT*dTdt/(k-1)))
    return(dPdt,dTdt,dQadt,dQpdt,dWdt)  
K=KReag
dxdt=0
DKdT=sp.diff(K,Ts)
dkdT=DKdT.subs(Ts,T0)
V=V0
DVdt=sp.diff(Vo,teta)
dVdt=DVdt.subs(teta,teta0)
K=KReag
k=K.subs(Ts,T0)
x0=[P0,T0,Qa0,Qp0,W0]
Trab=[]
Qap=[]
Qpa=[]
for l in range(len(tk)-1):
    ts=[tk[l],tk[l+1]]
    x=solve_ivp(motor,ts,x0,rtol=1e-10,atol=1e-10).y
    P1[l+1]=x[0][-1]
    T1[l+1]=x[1][-1]
    x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
    dxdt=0
    h=0
    V=Vo.subs(teta,tk[l+1])
    A=Ap.subs(teta,tk[l+1])
    dVdt=DVdt.subs(teta,tk[l+1])
    k=K.subs(Ts,T1[l+1])
    dkdT=DKdT.subs(Ts,T1[l+1])

    
for i in range(len(ti)-1): #range(0,98)
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
    dVdt=DVdt.subs(teta,ti[i+1])
    k=K.subs(Ts,T[i+1])
    dxdt=0
    dkdT=DKdT.subs(Ts,T[i+1])
    vg=2.28*vp
    h=0.013*D**(-0.2)*(P[i+1])**0.8*T[i+1]**(-0.53)*vg**0.8
    


X_teta=1-sp.exp(-a*((teta-teta_comb)/delta_teta)**(m+1))
dXdt=sp.diff(X_teta,teta)
x_teta=X_teta.subs(teta,teta_comb)
dxdt=dXdt.subs(teta,teta_comb)
wiebe=[]

for j in range(len(tk))[len(ti)-1:len(ti)+len(tj)-1]: #199 ao 198
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
    dVdt=DVdt.subs(teta,tk[j+1])
    #K=((teta-teta0)*(KProd)+(teta0+delta_teta-teta)*(KReag))/delta_teta
    #ks=K.subs(teta,tk[j+1])
    x_teta=X_teta.subs(teta,tk[j+1])
    K=(x_teta*KProd+(1-x_teta)*KReag)
    DKdT=sp.diff(K,Ts)
    
    dxdt=dXdt.subs(teta,tk[j+1])
    k=K.subs(Ts,T[j+1])
    dkdT=DKdT.subs(Ts,T[j+1])
    wiebe.append(x_teta)
    vg=2.28*vp+0.00324*(P[j+1]-P1[j+1])*Vd*T0/(P0*V0)
    h=0.013*D**-0.2*(P[j+1])**0.8*T[j+1]**-0.53*vg**0.8
      

K=KProd
DKdT=sp.diff(K,Ts)
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
    dVdt=DVdt.subs(teta,tk[f+1])
    k=K.subs(Ts,T[f+1])
    dxdt=0
    dkdT=DKdT.subs(Ts,T[f+1])
    vg=2.28*vp+0.00324*(P[f+1]-P1[f+1])*Vd*T0/(P0*V0)
    h=0.013*D**-0.2*(P[f+1])**0.8*T[f+1]**-0.53*vg**0.8


print(P,T,Qa,Qp,W)
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

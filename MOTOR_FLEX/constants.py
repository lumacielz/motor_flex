#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 11:34:38 2022

@author: luiza.maciel
"""


from math import radians,degrees,sin,cos,sqrt,log,exp
import numpy as np
import sympy as sp
from scipy import signal


#dados da geometria do motor
dados_motor=input('Dados do Motor: ').split(',')
#dados_motor=11,0.144,0.081,0.0864,4

r=float(dados_motor[0])
L=float(dados_motor[1])
D=float(dados_motor[2])
S=float(dados_motor[3])
cilindros=float(dados_motor[4])

avanco=float(input("Ângulo de avanço: "))
if avanco<=45 and avanco>36:
    fcor=(avanco-36)*((1.77-1.67)/(45-36))+1.67
elif avanco<=36 and avanco >27:
    fcor=(avanco-27)*((1.67-1.68)/(36-27))+1.68
elif avanco<=27 and avanco >=18:
    fcor=(avanco-18)*((1.68-1.45)/(27-18))+1.45


teta0d=float(input("ângulo de fechamento da válvula de admissão [em graus]: "))
teta_combd=float(input("ângulo de início da combustão[em graus]: "))
delta_tetad=float(input("duração da combustão[em graus]: "))
tetafd=float(input("ângulo de abertura da válvula de escape[em graus] "))

teta0=radians(teta0d)
tetaf=radians(tetafd)
teta_comb=radians(teta_combd)
delta_teta=radians(delta_tetad)

#condicoes de operacao
Var=float(input("Vazão de ar[kg/h]: "))
P0=float(input("Pressão no ângulo de fechamento: "))*10**3
Tp=float(input("Temperatura das paredes: "))
rotacao=float(input("rotação do motor: "))


R=S/2 #raio do virabrequim
N=rotacao*2*np.pi/60 #rad/s
Vd=np.pi*D**2*R/2 #volume deslocado m^3
vp=2*S*N #velocidade do pistão m/s

Ru=8.314 #kJ/kmolK
M_air = 28.962

#combustivel = "gasolina"
combustivel=input("Com qual combustível deseja operar?  ")
def stringToFloat(inputLista):
    lista=[]
    for i in inputLista.split(","):
        lista.append(float(i))
    return lista

def composicaoCombustivel(c):
    coef=c.split("C")[1].split("H")
    coef[1:]=(coef[1].split("O"))
    return list(map(lambda x: float(x),coef))

while True:
    if combustivel=="gasolina":
    #combustivel equivalente
    #    carb=6.67
    #    h=12.8
    #    o=0.533  
    #    n=0
    #    OC=9.6
    #    complete_combustion_H2O=6.4
    #    complete_combustion_CO2=6.67
        
        c,h,o,n=8*0.5462+2*0.4538,18*0.5462+0.4538*6,0.4538,0
        
        a0=-8102.9675849653
        a1=7963.204888950
        a2=-2923.24238668569
        a3=509.144026930886
        a4=-42.2764373650608
        a5=1.35175803975684
        
        h0_comb=0.5462*(-208450)+0.4538*(-235310)
        s0_comb=0.5462*(466.514)+0.4538*(282.444)
        mcomb=(12*c+h+16*o)
        
        PCI=39249*10**3 #kJ/kg
        
        Ae=(5.7*10**11) #octano
        Ea=15098
        me=0.25
        ne=1.5
        
        break
        
    elif combustivel=="alcool":
        #combustivel equivalente
    #    OC=3.2
    #    carb=2.15
    #    h=6.62
    #    o=1.23
    #    complete_combustion_CO2=2.15
    #    complete_combustion_H2O=3.31
        
        c,h,o,n=1.709,5.418801,0
    
        a0=-12482.8740213179
        a1=10263.4623453335
        a2=-3316.7850402598
        a3=526.309291795851
        a4=-40.8869367350809
        a5=1.24555084441151
       
        h0_comb=0.8547*(-235310)+0.1453*(-241826)
        s0_comb=0.8547*(282.444)+0.1453*(188.835)
        mcomb=(12*c+h+16*o)
        
        #PCI=
        
        Ae=(1.8*10**12) #etanol
        Ea=15098
        me=0.15
        ne=1.6
        
        break
    
    elif combustivel=="GNV":
        met=0.92285
        et=5.455*10**-2
        prop=1.115*10**-2
        but=0.125*10**-2
        co2=0.255*10**-2
        r=0.765*10**-2
        c=met+2*et+3*prop+4*but+co2
        h=met*4+et*6+prop*8+but*10
        o=co2*2
        n=r
    
        a0=-12713.9307245771
        a1=10272.7734235011
        a2=-3251.2253041433
        a3=504.476942937525
        a4=-38.3629908265519
        a5=1.14655193729902
    
        h0_met=-74873
        s0_met=186.251
        h0_et=-84740
        s0_et=229.597
        h0_prop=-103900
        s0_prop=269.917
        h0_but=-126200
        s0_but=306.647
        h0_co2=-393522
        s0_co2=213.795
        h0_n2=0
        s0_n2=191.609
        h0_comb=met*h0_met+et*h0_et+prop*h0_prop+but*h0_but+co2*h0_co2+r*h0_n2
        s0_comb=met*s0_met+et*s0_et+prop*s0_prop+but*s0_but+co2*s0_co2+r*s0_n2
        mcomb=(met*16+et*30+prop*44+but*58+co2*44+28*r)
        
        #PCI =
        break 
    
    elif combustivel == "H2":
        c,h,o,n=0,2,0,0
        mcomb=0
        break
    elif combustivel == "outro":
        molecularComposition = input("digite a composição do combustível na forma CnHnOn: ")
        c,h,o,n=composicaoCombustivel(molecularComposition),0
        mcomb=c*12+h+o*16
        h0_comb=float(input("entalpia de formação do combustivel: "))
        s0_comb=float(input("entropia de formação do combustivel: "))
        a0,a1,a2,a3,a4,a5=stringToFloat(input("coeficientes da equação logaritimica de Cp em funcao da temperatura: "))
     
    else:
        print("digite um combustível válido!!!")

#massa molecular do combustivel        
M_c=(12.0107*c+h*1.008+16*o+28*n)

#quantidades estequimetricas de produtos
stoichometricH2O = h/2
stoichometricCO2 = c
OC = (2*stoichometricCO2 + stoichometricH2O - o)/2 #mols de oxigenio
stoichometricN2 = 3.772*OC

AC = OC*(32 + 3.772 * 28.16) / M_c  #razao ar combustivel
M_m=(M_c+OC*(32+3.772*28.16))/(1+OC+OC*3.772)  #massa molecular da mistura

#massa de ar admitida
m_air=Var/(0.5*3600*cilindros*(N/(2*np.pi))) #rot por s
#massa de mistura admitida
m_m=m_air/AC+m_air 
#massa de combustivel admitida
m_c=m_air/AC
mols_c = m_c/M_c
#mols de oxigenio e nitrogenio admitidos
NI=stoichometricN2*mols_c
O=OC*mols_c

#fracao massica inicial
Y0 = m_c/m_m

Rg=Ru/M_m #constante dos gases

#DADOS PARA WIEBE!
coefwiebe=input('Coeficientes da equação de Wiebe m , a: ').split(',')
m=1
a=2
ef=0.87
Qtot=ef*m_c*PCI #kJ
#spark = 0

#composicao molar inicial
reactants_composition = np.asarray([mols_c,OC*mols_c,0,0,0,0,0,0,0,0,0,0,3.773*OC*mols_c])
estimative = np.asarray([mols_c,OC*mols_c,10**-15,10**-15,10**-15,10**-15,10**-15,10**-15,10**-15,10**-15,10**-15,10**-15,3.773*OC*mols_c])

Ts = sp.symbols('Ts')
Cp=a0+a1*sp.log(Ts)+a2*sp.log(Ts)**2+a3*sp.log(Ts)**3+a4*sp.log(Ts)**4+a5*sp.log(Ts)**5

comp=[reactants_composition[0]]

massFraction=[]

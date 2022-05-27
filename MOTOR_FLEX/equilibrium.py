from scipy.optimize import minimize,curve_fit
import sympy as sp
import pandas as pd
from decimal import  Decimal
from math import radians,log
import numpy as np
import sys

file='./tabelaHS.xlsx'
df2=pd.read_excel(file)

Ts=sp.symbols("Ts")

def stringToFloat(inputLista):
    lista=[]
    for i in inputLista.split(","):
        lista.append(float(i))
    return lista

def molarToMassFraction(moles,molecular_mass):
    mass_c = [M*nm for M,nm in zip(molecular_mass,moles)]
    mass_fractions = np.asarray(mass_c)/sum(mass_c)        
    return mass_fractions

def massToMolarFraction(mass,molecular_mass):
    molar_c = [mn/M for M,mn in zip(molecular_mass,mass)]
    molar_fractions = np.asarray(molar_c)/sum(molar_c)    
    return molar_fractions

def massToComposition(mass,molecular_mass):
    molar_c = [mn/M for M,mn in zip(molecular_mass,mass)]
    return np.asarray(molar_c)

def composicaoCombustivel(c):
    coef=c.split("C")[1].split("H")
    coef[1:]=(coef[1].split("O"))
    return list(map(lambda x: float(x),coef))

def Cps(Ts):
    CpH2O = 0.14904476*log(Ts)**5 -5.60452380*log(Ts)**4 + 81.54220310*log(Ts)**3 -573.43988184*log(Ts)**2 + 1954.64102748*log(Ts) -2558.72544088
    CpC2H5OH = 1.40634*log(Ts)**5 -45.60008*log(Ts)**4 + 578.32500*log(Ts)**3 -3579.01508*log(Ts)**2+ 10834.12151 *log(Ts)-12834.47773
    CpC8H18 = 4.9695987*log(Ts)**5 -157.61132*log(Ts)**4 +1966.430702 *log(Ts)**3 -12036.09466*log(Ts)**2 +36241.19904*log(Ts)-43029.69896
    return CpC8H18*oc+et*CpC2H5OH+h2o*CpH2O

def fitLog(Ts, a0,a1,a2,a3,a4,a5):
    return a0+a1*np.log(Ts)+a2*np.log(Ts)**2+a3*np.log(Ts)**3+a4*np.log(Ts)**4+a5*np.log(Ts)**5

tempRange = np.linspace(298.15,3000,1000)

Ru=8.314

#entalpias de formacao 
h0_o2=0
h0_h2o=-241820
h0_co2=-393520
h0_co=-110530
h0_o = 249190
h0_n = 472650
h0_no2 = 33100
h0_no = 90291
h0_oh = 39460
h0_h2 = 0
h0_h = 218000
h0_n2=0

def gibbs (composicoes,T,P):
   
    Ntot=sum(composicoes)
    
    deltaS=interpolacaoS(T)
    constant=Ru*T
    
    h_comb,h_o2, h_h2o,h_co2,h_co,h_o,h_n,h_no2,h_no,h_oh,h_h2,h_h,h_n2=interpolacaoH(T)
    
    gt_comb=(h_comb-T*deltaS["deltasComb"])/constant
    gt_o2=(h_o2-T*deltaS["deltasO2"])/constant
    gt_h2o=(h_h2o-T*deltaS["deltasH2O"])/constant
    gt_co2=(h_co2-T*deltaS["deltasCO2"])/constant
    gt_co=(h_co-T*deltaS["deltasCO"])/constant
    gt_o = (h_o - T * deltaS["deltasO"]) / constant
    gt_n = (h_n - T * deltaS["deltasN"]) / constant
    gt_no2 = (h_no2 - T * deltaS["deltasNO2"]) / constant
    gt_no = (h_no - T * deltaS["deltasNO"]) / constant
    gt_oh = (h_oh - T * deltaS["deltasOH"]) / constant
    gt_h2 = (h_h2 - T * deltaS["deltasH2"])/ constant
    gt_h = (h_h - T * deltaS["deltasH"]) / constant
    gt_n2=(h_n2-T*deltaS["deltasN2"])/constant

    gt=[gt_comb,gt_o2,gt_h2o,gt_co2,gt_co,gt_o,gt_n,gt_no2,gt_no,gt_oh,gt_h2,gt_h,gt_n2]
    deltaG=0
    for n,g in zip(composicoes,gt):
        if n!=0 and n>sys.float_info.min:
            deltaG+=n * (g + float((Decimal(n) * Decimal("%.15f" % P) / Decimal(Ntot)).ln()))
   
    return (deltaG)

def restricaoC(composicoes):
    return mols_c*c-composicoes[0]*c-composicoes[3]-composicoes[4]
def restricaoH(composicoes):
    return mols_c*h-composicoes[0]*h-composicoes[2]*2-composicoes[9]-composicoes[10]*2-composicoes[11]
def restricaoO(composicoes):
    return mols_c*o+2*O-composicoes[0]*o-2*composicoes[1]-composicoes[2]-composicoes[3]*2-composicoes[4]-composicoes[5]-composicoes[7]*2-composicoes[8]-composicoes[9]
def restricaoN(composicoes):
    return mols_c*n+2*NI-n*composicoes[0]-composicoes[12]*2-composicoes[6]-composicoes[7]-composicoes[8]

def interpolacaoH(i):
    for k in range(36):
        if i>=df2['T'][k] and i<=df2['T'][k+1]:
            deltahH2O=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hH2O'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hH2O'][k]
            deltahCO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hCO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hCO2'][k]
            deltahO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hO2'][k]
            deltahCO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hCO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hCO'][k]
            deltahN2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hN2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hN2'][k]
            deltahNO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hNO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hNO'][k]
            deltahNO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hNO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hNO2'][k]
            deltahOH=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hOH'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hOH'][k]
            deltahH2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hH2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hH2'][k]
            deltahH=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hH'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hH'][k]
            deltahN=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hN'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hN'][k]
            deltahO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hO'][k]
            
    
    h_o2 = h0_o2 + deltahO2    
    h_h2o = h0_h2o + deltahH2O
    h_co2 = h0_co2 + deltahCO2
    h_co = h0_co + deltahCO
    h_o = h0_o + deltahO
    h_n = h0_n + deltahN
    h_no2 = h0_no2 + deltahNO2
    h_no = h0_no + deltahNO
    h_oh = h0_oh + deltahOH
    h_h2 = h0_h2 + deltahH2
    h_h = h0_h + deltahH
    h_n2 = h0_n2 + deltahN2

    if combustivel=="H2":
        h_comb=h_h2
        h_h2=0
    else:
        h_comb = float(h0_comb + sp.integrate(Cp, (Ts, 298.15, i)))
    
    return h_comb,h_o2,h_h2o,h_co2,h_co,h_o,h_n,h_no2,h_no,h_oh,h_h2,h_h,h_n2

def interpolacaoS(i):
    for k in range(36):
            if i>=df2['T'][k] and i<=df2['T'][k+1]:
                deltasH2O=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sH2O'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sH2O'][k]
                deltasCO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sCO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sCO2'][k]
                deltasO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sO2'][k]
                deltasCO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sCO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sCO'][k]
                deltasN2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sN2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sN2'][k]
                deltasNO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sNO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sNO'][k]
                deltasNO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sNO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sNO2'][k]
                deltasOH=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sOH'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sOH'][k]
                deltasH2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sH2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sH2'][k]
                deltasH=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sH'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sH'][k]
                deltasN=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sN'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sN'][k]
                deltasO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sO'][k]
                if combustivel=="H2":
                    deltasComb=deltasH2
                    deltasH2=0
                else:
                    deltasComb=float(s0_comb + sp.integrate(Cp/Ts,(Ts,298.15,i)))
                
    return {"deltasComb":deltasComb,"deltasO2":deltasO2,"deltasH2O":deltasH2O,"deltasCO2":deltasCO2,"deltasCO":deltasCO,"deltasO":deltasO,"deltasN":deltasN,"deltasNO2":deltasNO2,"deltasNO":deltasNO,"deltasOH":deltasOH,"deltasH2":deltasH2,"deltasH":deltasH,"deltasN2":deltasN2}


def equilibrium(estimative,T,P):
    sol = minimize(gibbs, estimative,args=(T,P), method='SLSQP', bounds=bnds, constraints=cons)
    return sol.x
#kJ
def heatRelease(initialComposition,finalComposition,T):
    entalpy0,entalpy1=0,0
    for c,h in zip(initialComposition,interpolacaoH(T)):
        entalpy0 +=c*h
    for c2,h2 in zip(finalComposition,interpolacaoH(T)):
        entalpy1 +=c2*h2
    return -(entalpy1-entalpy0)*10**3
        
def temperaturaAdiabatica(composicao,T0,T1,entalpy_reactants):
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
    
    
    while True:
        h_comb,h_o2,h_h2o,h_co2,h_co,h_o,h_n,h_no2,h_no,h_oh,h_h2,h_h,h_n2=interpolacaoH(T0)
        entalpy_products0=Ncomb*h_comb+NO2*h_o2+NH2O*h_h2o+NCO2*h_co2+NCO*h_co+NO*h_o+NN*h_n+NNO2*h_no2+NNO*h_no+NOH*h_oh+NH2*h_h2+NH*h_h+NN2*h_n2
        zero0=entalpy_products0-entalpy_reactants
        
        h_comb,h_o2,h_h2o,h_co2,h_co,h_o,h_n,h_no2,h_no,h_oh,h_h2,h_h,h_n2=interpolacaoH(T1)
        entalpy_products1=Ncomb*h_comb+NO2*h_o2+NH2O*h_h2o+NCO2*h_co2+NCO*h_co+NO*h_o+NN*h_n+NNO2*h_no2+NNO*h_no+NOH*h_oh+NH2*h_h2+NH*h_h+NN2*h_n2
        zero1=entalpy_products1-entalpy_reactants
        
        nextT=T1-zero1*(T1-T0)/(zero1-zero0)
        
        if abs(zero1)<0.2:
            break
        else:
            T0=T1
            T1=nextT
    return nextT

def equilibrioAdiabatico(estimative,T,P,reactants_composition,T0):
    entalpy_reactants = 0
    for n, h in zip(reactants_composition, interpolacaoH(T0)):
        entalpy_reactants += n*h
    
    while True:
        sol = equilibrium(estimative, T,P)
        
        nextTemp=float(temperaturaAdiabatica(sol,T,T+2,entalpy_reactants))
        
        if abs(nextTemp-T)<0.1:
            break
        else:
            T=nextTemp
            estimative=sol
            
    total_mols = sum(i for i in sol)
    molar_fraction = [s/total_mols for s in sol]        
    mass_fraction = molarToMassFraction(sol,Molecular_Mass)       
    return {'Tad': nextTemp, 'moles': sol, 'molar fraction': molar_fraction, 'mass fraction': mass_fraction}

cons=[{'type': 'eq', 'fun': restricaoC},{'type': 'eq', 'fun': restricaoO},{'type': 'eq', 'fun':restricaoH},{'type': 'eq', 'fun': restricaoN}]

def equilibrioAdiabaticoVolumeConstante(estimative, T,P0,V,reactants_composition,T0):
    entalpy_reactants = 0
    P = P0
    for nt, ht in zip(reactants_composition, interpolacaoH(T0)):
        entalpy_reactants += nt*ht
    while True:
        sol = equilibrium(estimative, T,P)
        
        nextTemp=float(temperaturaAdiabatica(sol,T,T+2,entalpy_reactants))
        
        if abs(nextTemp-T)<0.5:
            break
        else:
            T=nextTemp
            total_mols = sum(i for i in sol)
            P = total_mols*Ru*T/(V*101.325)
            estimative=sol
    
    molar_fraction = [s/total_mols for s in sol]
    mass_fraction = molarToMassFraction(sol,Molecular_Mass)
    return {'Tad': nextTemp, 'moles': sol, 'molar fraction': molar_fraction, 'mass fraction': mass_fraction, 'pressure': P}


def bounds (composicao):
    bnds=[]
    for cp in composicao:
        if cp==0:
            bnds.append((0.0,sys.float_info.min))
        else:
            bnds.append((0,None))
    return tuple(bnds)


if __name__ == '__main__'   :
    T0=float(input("Temperatura dos reagentes: "))
    P0=float(input("Pressão[atm]:  "))
    estimateT=float(input("Estimativa da temperatura final: "))
    combustivel=input("combustível [gasolina,etanol,GNV,H2,outro]")
    while True:
        if combustivel=="gasolina":
            ocV = float(input("%V octano: "))
            etV = float(input("%V etanol: "))
            n_oc = ocV*703/(114.232)
            n_et = etV*789/(46.069)
            
            oc = n_oc/(n_oc+n_et)
            et = n_et/(n_oc+n_et)
            
            h2o=0
            c,h,o,n=8*oc+2*et,18*oc+et*6,et,0
            
            cps = [Cps(t) for t in tempRange]
            a0,a1,a2,a3,a4,a5=curve_fit(fitLog,tempRange,cps)[0]
            Cp=a0+a1*sp.log(Ts)+a2*sp.log(Ts)**2+a3*sp.log(Ts)**3+a4*sp.log(Ts)**4+a5*sp.log(Ts)**5  #J/molK
            
            h0_comb=oc*(-208450)+et*(-235310)
            s0_comb=oc*(466.514)+et*(282.444)
            mcomb=(12*c+h+16*o)
            
            break
            
        elif combustivel=="alcool":
            et = 0.938756707502278
            h2o = 0.06124329249772197
            oc=0
            c,h,o,n = et*2,6*et+2*h2o,h2o+et,0
            
            cps = [Cps(t) for t in tempRange]
            a0,a1,a2,a3,a4,a5=curve_fit(fitLog,tempRange,cps)[0]
            Cp=a0+a1*sp.log(Ts)+a2*sp.log(Ts)**2+a3*sp.log(Ts)**3+a4*sp.log(Ts)**4+a5*sp.log(Ts)**5  #J/molK
            
            h0_comb=et*(-235310)+h2o*(-241826)
            s0_comb=et*(282.444)+h2o*(188.835)
            mcomb=(12*c+h+16*o)
            break
        
        elif combustivel=="GNV":
            met=0.92285
            et=5.455*10**-2
            prop=1.115*10**-2
            but=0.125*10**-2
            co2=0.255*10**-2
            n=0.765*10**-2
            c=met+2*et+3*prop+4*but+co2
            h=met*4+et*6+prop*8+but*10
            o=co2*2
            
        
            a0=-12713.9307245771
            a1=10272.7734235011
            a2=-3251.2253041433
            a3=504.476942937525
            a4=-38.3629908265519
            a5=1.14655193729902
            Cp=a0+a1*sp.log(Ts)+a2*sp.log(Ts)**2+a3*sp.log(Ts)**3+a4*sp.log(Ts)**4+a5*sp.log(Ts)**5  #J/molK
            
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
            h0_comb=met*h0_met+et*h0_et+prop*h0_prop+but*h0_but+co2*h0_co2+n*h0_n2
            s0_comb=met*s0_met+et*s0_et+prop*s0_prop+but*s0_but+co2*s0_co2+n*s0_n2
            mcomb=(met*16+et*30+prop*44+but*58+co2*44+28*n)
            break 
        
        elif combustivel == "H2":
            c,h,o,n=0,2,0,0
            mcomb=0
            break
        elif combustivel == "outro":
            molecularComposition = input("digite a composição do combustível na forma CnHnOn: ")
            c,h,o=composicaoCombustivel(molecularComposition)
            n=0
            mcomb=c*12+h+o*16
            h0_comb=float(input("entalpia de formação do combustivel: "))
            s0_comb=float(input("entropia de formação do combustivel: "))
            a0,a1,a2,a3,a4,a5=stringToFloat(input("coeficientes da equação logaritimica de Cp em funcao da temperatura: "))
            Cp=a0+a1*sp.log(Ts)+a2*sp.log(Ts)**2+a3*sp.log(Ts)**3+a4*sp.log(Ts)**4+a5*sp.log(Ts)**5  #J/molK
            
        else:
            print("digite um combustível válido!!!")


    mols_c=float(input("quantidade de combustivel: "))
    oxid=input("[0] ar atmosferico [1] somente O2 ")
    O=(input("quantidade de oxigenio: "))
    if O == "estequiometrica":
        stoichometricH2O = h/2
        stoichometricCO2 = c
        O = (2*stoichometricCO2 + stoichometricH2O - o)/2 * mols_c
    else:
        O = float(O)
    NI=0
    if oxid == "0":
        NI=3.76*O
    M_c=12.0107*c+h*1.008+16*o+28*n
        #comb,o2,h2o,co2,co,o,n,no2,no,oh,h2,h,n2
    Molecular_Mass=[M_c,31.999,18.015,44.01,28.01,16,14.007,46.005,30.006,17.007,2.016,1.008,28.013]

    total_mass=mols_c*M_c+32*O+NI*28
        
    estimative=stringToFloat(input("estimativa inicial para a composicao de equilibrio [Combustivel,O2,H2O,CO2,CO,O,N,NO2,NO,OH,H2,H,N2]: "))
    bnds = bounds(estimative)
    entalpy_reactants=interpolacaoH(T0)[0]*mols_c+O*interpolacaoH(T0)[1]+NI*interpolacaoH(T0)[12]
    reactants_composition = [mols_c,O,0,0,0,0,0,0,0,0,0,0,NI]
    total_mols = mols_c + O + NI
    
    P=P0
    T=estimateT
    V = total_mols*Ru*T0/(P0* 101.325)
    names = ['Combustivel','O2','H2O','CO2','CO','O','N','NO2','NO','OH','H2','H','N2']
#1,9.6,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,1e-15,36.096   
#composicao final
    restricao = input("restrições [P-T] | [H-P] | [H-V]: ")
    if restricao == "[P-T]":
        sol = equilibrium(estimative,T,P)
        total_mols=sum(sol)
        
        print("\n ------------------COMPOSICAO EQUILIBRIO ------------------- \nMOLES: ")
        for i in range(len(names)):
            print(names[i] + ': ' + str(sol[i]))
        print('\n MASS FRACTION:')   
        for i in range(len(names)):
            print(names[i] + ': ' + str(sol[i]*Molecular_Mass[i]/total_mass))
        print('\n MOLAR FRACTION:')   
        for i in range(len(names)):
            print(names[i] + ': ' + str(sol[i]/total_mols))
    if restricao == "[H-P]":
        adiab=equilibrioAdiabatico(estimative,T,P,reactants_composition,T0)
        
        print('\n ------------------TEMPERATURA ADIABATICA----------------------')
        print("Tad: " + str(adiab['Tad']))
        for n in range(len(names)):
            print(names[n] + ': ' + str(adiab['molar fraction'][n]))
    if restricao == "[H-V]":
        adiab=equilibrioAdiabaticoVolumeConstante(estimative,T,P,V,reactants_composition,T0)
        
        print('\n ------------------TEMPERATURA ADIABATICA----------------------')
        print("Tad: " + str(adiab['Tad']))
        for n in range(len(names)):
            print(names[n] + ': ' + str(adiab['molar fraction'][n]))
    
else:
    from constants import *
    bnds=bounds(estimative)
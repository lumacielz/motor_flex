from scipy.integrate import solve_ivp
from scipy.misc import derivative
from math import exp
import numpy as np
import sympy as sp

def funcao (t,x):
    c=derivative(c1,t,dx=1e-10)
    print(c)
    e=x[0]
    dxdt=c*np.exp(-e)
    return dxdt
def c1(time):
    return time**3

s=solve_ivp(funcao,[1,2],[2])
  
a=sp.symbols('s')
c2=a**3

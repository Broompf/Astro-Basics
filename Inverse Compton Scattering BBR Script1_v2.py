import streamlit as st
import sympy as sym
import numpy as np
#Programme Conditionals---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------X
#constants
h = 6.626e-34
L = 2.99e8
k = 1.38
r = 2.817e-15
P = 10
Q = 10**5
m = 9e-31
s = 2.3676e12

#Webapp Layout------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------X

st.title("Inverse Compton Spectra for Single Scattering")
st.header("Black Body Radiation Condition")

st.sidebar.write("input values")
p = st.sidebar.number_input("Enter value for Momentum: ",value=2.5)
E = st.sidebar.number_input("Enter value for Energy of the Seed Photon: ",value=2.72)
T = st.sidebar.number_input("Enter value for Temperature : ",value=2.72)
B = st.sidebar.number_input("Enter value for Magnetic Field: ",value=2.72)
l = st.sidebar.number_input("Enter value for L: ",value=2.72)

lower_E=st.sidebar.number_input("Enter value for lower limit of epsilon: ",value=1)
upper_E=st.sidebar.number_input("Enter value for lower limit of epsilon: ",value=100)

#Programme Calculations---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------X
def rangeE(l,u):
  liste=[]
  for i in range(float(l),float(u)+1):
    liste.append(i)
  return liste
liste= rangeE(lower_E,upper_E)
st.info(f" length of e {len(liste)}, first term = {liste[0]}, last term = {liste[len(liste)]}")

def findc(l,B):
  y = 3-p
  G = (P**y) - (Q**y)
  v = m*(L**2)

  N = (l/(s*(B**2))*(v**y))*(y/((Q**y)-(P**y)))

  c = (N*(v**(-p)))
  return c
def vol_emmissivity(p,E,T):
    a = p + 3
    #print(f"a={a}")                                                        
    b = (p + 5) / 2
    #print(f"b={b}")  
    c = findc(l,B)
    x = (p - 1) / 2
    #print(f"x={x}")                                                        
    A = ((p**2) + 4*p + 11) / ((a**2)*(2*b)*(p + 1))
    #print(f"A={A}")
    g = sym.gamma(b)
    #print(f"g={g}")
    Z = sym.zeta(b)
    #print(f"Z={Z}")
    F = A*g*Z
    #print(f"F={F}")
    i = 3.141592653589
    d = (8*(i**2)*(r**2)) / ((h**3)*(L**2))
    #print(f"d={d}")
    o = (k*T)**b
    #print(f"o={o:e}")

    X = (c*d*o*F*(E*(-x)))*((1e-23)**b)
    exp = X *(1.6*(10**(-16)))
    return exp

#Results Tab--------------------------------------------------------------------------------------------------------------------------------------------------------X

value=vol_emmissivity(p,E,T)
st.info(f"Volume Emmissivity is {value:e} Js-1KeV-1K-1")

st.write(" ")
st.header("Power Law Condition")

#Importing required Modules-----------------------------------------------------------------------------------------------------------------------------------------------------------------
import streamlit as st
import sympy as sym
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#constants----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
L = 2.99e8               #Speed of Light in m/s
r = 2.817e-15            #Classical Radius of Electron
R = 10                   #Lorentz Function Gamma 1
Q = 10**5                #Lorentz Function Gamma 2
m = 9e-31                #Mass of Electron
s = 2.3676e12            #Constant C3
h = 6.626e-34            #Planck's Constant

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
st.sidebar.write("input values")
p = st.sidebar.number_input("Enter value for Momentum: ",value=2.5)
B = st.sidebar.number_input("Enter value for Magnetic Field: ",value=2.72)
q = st.sidebar.number_input("Enter value for L: ",value=2.72)
alpha = st.sidebar.number_input("Enter value for alpha:",value=2.72)
F = st.sidebar.number_input("Enter value for Fvo:",value=2.72)
V = st.sidebar.number_input("Enter value for vo:",value=2.72)
lower_E=st.sidebar.number_input("Enter value for lower limit of epsilon: ",value=1)
upper_E=st.sidebar.number_input("Enter value for lower limit of epsilon: ",value=100)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#User Defined Functions---------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Defining P1 constant-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Defining the variable c in terms of p

def findc(q,B,p):
  y = 3-p
  G = (R**y) - (Q**y)
  v = m*(L**2)

  N = (q/(s*(B**2))*(v**y))*(y/((Q**y)-(R**y)))

  c = (N*(v**(-p)))
  return c

def Const1 (p):
    a = p + 3
    #print(f"a={a}")                                                        
    b = (p + 5) / 2
    #print(f"b={b}")  
    c = findc(q,B,p)
    x = (p - 1) / 2
    #print(f"x={x}")                                                        
    A = ((p**2) + 4*p + 11) / ((a**2)*(2*b)*(p + 1))
    i = 3.141592653589
    P1 = i*L*(r**2)*c*A
    
    return P1

#Defining the Integral V(Epsilon)Epsilon^(p-1/2) i.e. Constant P2---------------------------------------------------------------------------------------------------------------------------



# Extracting data from the Uploaded file

def createdata(file):
    df = pd.read_csv(file, delim_whitespace=True)
    # Check if the required columns exist in the DataFrame
    if 'Epsilon' not in df.columns or 'V_Epsilon' not in df.columns:
        st.error("Required columns 'Epsilon' or 'V_Epsilon' are missing in the dataset.")
        return

    # Extract the relevant data
    E = df['Epsilon']
    V = df['V_Epsilon']
    
    return E,V

# File uploader

st.sidebar.write("Upload a text file containing data of curve of V(Epsilon) against Epsilon")
uploaded_file = st.sidebar.file_uploader("Choose a text file", type="txt")
st.sidebar.write("**columns should be named 'Epsilon', 'V_Epsilon'")
#data collection and unit correction
if uploaded_file is None:
    st.sidebar.write("using sample dataset. upload file and provide asked values to process other dataset")
    E,V = createdata('Sample.txt')
   

if uploaded_file is not None:
    E,V = createdata(uploaded_file)
  

# -- Simpson's 1/3rd Rule Function --
def simpsons_one_third(E, V, p):
    """
    Applies Simpson's 1/3 Rule to the integral of V(E)*E^((p-1)/2) over E.
    Assumes E is equally spaced.
    """
    n = len(E) - 1

    if n < 2:
        raise ValueError("Need at least 3 points for Simpson's rule.")
    if n % 2 != 0:
        raise ValueError("Simpson's rule requires an even number of intervals (odd number of points).")

   # st.write("Epsilon List:", E)
    #st.write("V_Epsilon List:", V)
   
    #st.success(f" E[1] = {E[1]} with type {type(E[1])}")
    
    h = E[1] - E[0]
    if not np.allclose(np.diff(E), h, rtol=1e-5, atol=1e-8):
        raise ValueError("Epsilon values must be equally spaced.")

    # Transform the function as V(E) * E^((p-1)/2)
    transformed = (V * E**((p - 1) / 2)).to_numpy()
  
    # Apply Simpson's Rule
    result = transformed[0] + transformed[-1] + \
             4 * sum(transformed[1:n:2]) + \
             2 * sum(transformed[2:n-1:2])

    P2 = (h / 3) * result
    
    return P2

def int2(E, p, alpha):
    Ei = E.iloc[0]
    Ef = E.iloc[-1]
   
    if alpha == (p-1)/2:
        if Ei > 0:
            diff = np.log(Ef) - np.log(Ei)
        else:
            raise ValueError("Logarithm undefined for non-positive values.")
      
    else:
        diff = ((Ef**(((p-1)/2)-alpha))/(((p-1)/2)-alpha)) - ((Ef**(((p-1)/2)-alpha))/(((p-1)/2)-alpha))
        
    return diff

#Defining for Epsilon1 list-----------------------------------------------------------------------------------------------------------------------------------------------------------------

def rangeE(l,u):
    liste=[]
#  st.write(f" l = {l} {type(l)}")
    for i in range(l,u + 1):
        liste.append(i)

    return liste
liste= rangeE(lower_E,upper_E)
#st.info(f" length of e {len(liste)}, first term = {liste[0]}, last term = {liste[len(liste)-1]}")

#Graph--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Blackbody Radiation Graph------------------------------------------------------------------------------------------------------------------------------------------------------------------
def plot_it(liste,final,x_label,y_label,title):
    plt.figure(figsize=(10, 6))
    
    plt.xlabel(x_label)  # Set the x-axis label
    plt.ylabel(y_label)  # Set the y-axis label
    
    plt.plot(liste, final, color='grey', alpha=0.5)  # Connect points with a line
    plt.title(title)
    logscale=st.toggle("Show Graph in logscale", value=True)
    if logscale==True:
        plt.xscale('log')
        plt.yscale('log')
    plt.grid(True)
    plt.legend()
    st.pyplot(plt)
#Power Law Graph----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def plot_it(liste,final2,x_label,y_label,title):
    plt.figure(figsize=(10, 6))
    
    plt.xlabel(x_label)  # Set the x-axis label
    plt.ylabel(y_label)  # Set the y-axis label
    
    plt.plot(liste, final2, color='grey', alpha=0.5)  # Connect points with a line
    plt.title(title)
    logscale=st.toggle("Show Graph in logscale", value=True)
    if logscale==True:
        plt.xscale('log')
        plt.yscale('log')
    plt.grid(True)
    plt.legend()
    st.pyplot(plt)

#Calculating Volume emissivity--------------------------------------------------------------------------------------------------------------------------------------------------------------
#Defining the constant terms as a single term

P1 = Const1(p)  # Call the function to get P1

#E = np.array(E, dtype=float)
#V = np.array(V, dtype=float)

P2 = simpsons_one_third(E, V, p)

diff = int2(E, p, alpha)

o = (F*(V**alpha)) / (L*(h**(1-alpha)))

#st.write(f"P1: {P1}")
#st.write(f"P2: {P2}")
         
Constt1 = (P1)*(P2)

Constt2 = (P1)*(o)*(diff)
#st.write(f"Constt={Constt}")

#Obtaining the final result for volume emissivity in Js^-1KeV^-1K^-1

final = [Constt1 * (i**(-((p-1)/2))) * 1.6e-16 for i in liste]
final2 = [Constt2 * (i**(-((p-1)/2))) * 1.6e-16 for i in liste]

#Table of Dataset---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

data = pd.DataFrame({'Epsilon': liste, 'Volume Emmissivity': final})
#data["Epsilon"] = data["Epsilon"].apply(lambda x: '{:.6e}'.format(x))
#data["Volume Emmissivity"] = data["Volume Emmissivity"].apply(lambda x: '{:.6e}'.format(x))

data2 = pd.DataFrame({'Epsilon': liste, 'Volume Emmissivity': final2})
#data2["Epsilon"] = data["Epsilon"].apply(lambda x: '{:.6e}'.format(x))
#data2["Volume Emmissivity"] = data["Volume Emmissivity"].apply(lambda x: '{:.6e}'.format(x))


#Streamlit app layout-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
st.title("Inverse Compton Spectra for Single Scattering")
st.header("Black Body Condition")

st.write ("")
st.dataframe(data, use_container_width=True)
plot_it(liste, final, 'Epsilon1', 'Volume Emissivity', 'Inverse Compton Result')

st.write ("")
st.header ("Power Law Distribution")
st.write ("")
st.dataframe (data2, use_container_width=True)
plot_it(liste, final2, 'Epsilon1', 'Volume Emissivity', 'Inverse Compton Result')





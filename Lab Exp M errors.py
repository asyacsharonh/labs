#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd

# Load data from CSV file
data = np.loadtxt('dataexp.csv', delimiter=',')

# VALUES, denoted as (variable) and d(variable) for the uncertainty
V0_exp = data[:,0] # from csv data collected in lab
dV0 = 0.1

d =  0.00756
dd = 0.00001

n_vals = data[:,5]  # NOTE: might change due to temp, this value is for 20°C
dn_vals = data[:,6]

v0 = data[:,1] # from csv data
dv0 = data[:,2]

p = 886
dp = 1

sigma_vals = data[:,7]
dsigma_vals = data[:,-1]
print(dsigma_vals)

g = 9.81174734
dg = 0.00000003

# Lists to store results
ne_calc = []
ne_error_calc = []

# Loop over data to calculate ne and its uncertainty
for i in range(len(V0_exp)):
    V0 = float(V0_exp[i])          
    v = float(v0[i])               
    dv = float(dv0[i])             
    n = float(n_vals[i])           
    dn = float(dn_vals[i])        
    sigma = float(sigma_vals[i])   
    dsigma = float(dsigma_vals[i])

    # Calculate ne using the current V0 and v values
    ne = - (9 * d * np.pi * np.sqrt(2) / V0) * np.sqrt((n**3 * v**3) / (g * (p - sigma)))
    
    # Calculate partial derivatives for the current values of V0 and v
    partial_d = -(9 * np.pi * np.sqrt(2) / V0) * np.sqrt((n**3 * v**3) / (g * (p - sigma)))
    partial_V0 = (9 * d * np.pi * np.sqrt(2) / V0**2) * np.sqrt((n**3 * v**3) / (g * (p - sigma)))
    partial_n = -(27/2 * d * np.pi * np.sqrt(2) / V0) * np.sqrt((n * v**3) / (g * (p - sigma)))
    partial_v0 = -(27/2 * d * np.pi * np.sqrt(2) / V0) * np.sqrt((n**3 * v) / (g * (p - sigma)))
    partial_g = (9/2 * d * np.pi * np.sqrt(2) / V0) * np.sqrt((n**3 * v**3) / (g**3 * (p - sigma)))
    partial_p = (9/2 * d * np.pi * np.sqrt(2) / V0) * np.sqrt((n**3 * v**3) / (g * (p - sigma)**3))
    partial_sigma = -(9/2 * d * np.pi * np.sqrt(2) / V0) * np.sqrt((n**3 * v**3) / (g * (p - sigma)**3))

    # Calculate total uncertainty in ne for the current values of V0 and v
    delta_ne = np.sqrt(
        (partial_V0 * dV0)**2 +
        (partial_d * dd)**2 +
        (partial_n * dn)**2 +
        (partial_v0 * dv)**2 +
        (partial_g * dg)**2 +
        (partial_p * dp)**2 +
        (partial_sigma * dsigma)**2
    )
    
    # Append the calculated ne and uncertainty to the lists
    ne_calc.append(ne)
    ne_error_calc.append(delta_ne)

# Convert lists to arrays for output
ne_calc = np.array(ne_calc)
ne_error_calc = np.array(ne_error_calc)

df = pd.DataFrame({
    'ne_calc': ne_calc,
    'ne_error_calc': ne_error_calc
})

df.to_csv('ne_results.csv', index=False)

# Results
print(df)

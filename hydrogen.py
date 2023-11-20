import numpy as np
import math as m
import scipy.integrate
import matplotlib.pyplot as plt
from thermopack.saftvrmie import saftvrmie
from thermopack.saftvrqmie import saftvrqmie
from thermopack.cubic import cubic
#from thermopack.cubic import cubic
#from thermopack.cpa import cpa
import pandas as pd
from tabulate import tabulate 


#eos = saftvrqmie('H2', minimum_temperature = 20) # SAFT-VRQ Mie EoS for pure H2
eos_mix = saftvrqmie('P-H2,O-H2', minimum_temperature=20)
srk = cubic('P-H2,O-H2', 'SRK')
# eos_saftvrmie = saftvrmie('H2')
x = z = n = [0.5,0.5] # Total Molar composition


#get envelope
T_m, p_m = eos_mix.get_envelope_twophase(1e5, z)
T_srk, p_srk = srk.get_envelope_twophase(1e5, n)
#plt.plot(p, T,'b', label="Pure hydrogen") # Tp-projection of phase envelope

figure, axis = plt.subplots(2, 2)
axis[0, 0].plot(T_srk, p_srk, 'b', label="Ortho/para hydrogen, SRK")
axis[0, 0].plot(T_m, p_m, 'y', label="O-H2 P-H2 mix, Saftvrqmie")
axis[0,0].set_ylabel("pressure [Pa]")
axis[0,0].set_xlabel("Temperature [K]")
axis[0,0].set_title('Saftvrqmie vs. SRK')



def f(x_o):
     T  =  30
     calj = 4.1840
     h_p = 2023.1*calj  #300K 
     h_o = 2040.87*calj
     s_p = 31.212*calj #ved 300K
     s_o = 33.446*calj
  
     n = [1-x_o, x_o]
     lnphi = eos_mix.thermo(T, 1e5, n, eos_mix.VAPPH)  

     phi_p, = m.exp(lnphi[0])  # Replace 0 with the index of the relevant value
     phi_o, = m.exp(lnphi[1])
     DeltaG_o = (h_o-h_p)-T*(s_o-s_p)
     R = 8.314 # kJ/molK
     
     eq = x_o/(1-x_o) - (phi_p/phi_o)*m.exp(-DeltaG_o/(R*T))
     return eq

initial_guess = 0.5

# Use scipy.optimize.fsolve to find the root of the equation
result = scipy.optimize.root(f, initial_guess)

print("Optimal x_o:", result[0])

f(0.75)
vg_values = []
#specific volume for the ortho/para mxture
for T, p in zip(T_m, p_m):
    vg, = eos_mix.specific_volume(T, p, z, eos_mix.VAPPH)
    vg_values.append(vg,)

table1 = zip(T_m, p_m, vg_values)
print(tabulate(table1, headers = ('T [K]', 'P [kPa]', 'Vg [L/kg]')))
print()

#for T, p, vg in zip(T_m, p_m, vg_values):
#    print(tabulate(T,p,vg), headers=['T [K]', 'P [kPa]', 'Vg [L/kg]'])
#    print(f"Temperature = {round(T,1)}K, Pressure = {(p*1e-3)} kPa, Specific volume (vapor phase) = {round(vg,5)} L/kg")


#print('Temperature:',T)
# srk vs kvantetilstand ved lavere og lavere temperaturer, når begynner de å avvike fra hverandre, #
# finne eksperimentell data for hydrogen ved lave temperaturer
#ved å bruke dannelsesentalpier for orto og para, og finne likevektsfunksjoner
# tetthet, specific volume


u_p_values = []
u_o_values =[]

phi_p_values = []
phi_o_values = []

fug_o_values = []
fug_p_values = []

calj = 4.1840 #hvor 1cal = 4.18400 joule

for T, p, vg in zip(T_m, p_m, vg_values):
    R = 8.314 #J/kmol
    h_p = 2023.1*calj  #300K 
    h_o = 2040.87*calj
    s_p = 31.212*calj #ved 300K
    s_o = 33.446*calj
    eos_mix.set_ideal_enthalpy_reference_value(1,h_p)
    eos_mix.set_ideal_enthalpy_reference_value(2,h_o)
    eos_mix.set_ideal_entropy_reference_value(1,s_p)
    eos_mix.set_ideal_entropy_reference_value(2,s_o)

    u_, = eos_mix.chemical_potential_tv(T, vg, n)
    phi, = eos_mix.thermo(T, 1e5, n, eos_mix.VAPPH)
    fug, = eos_mix.fugacity_tv(T, vg, n)

    fug_p_values.append(fug[0],)
    fug_o_values.append(fug[1],)

    u_p_values.append(u_[0]/1000)
    u_o_values.append(u_[1]/1000)

    phi_p_values.append(m.exp(phi[0]),)
    phi_o_values.append(m.exp(phi[1]),)

table2 = zip(T_m, u_p_values, u_o_values, phi_p_values, phi_o_values, fug_p_values, fug_o_values)
print(tabulate(table2, headers = ('T [K]', 'u_ [kJ/mol](P)', 'u_ [kJ/mol](O)', 'fug. coeff(P)', 'fug. coeff(O)', 'fugacities(P)','fugacities(O)')))


axis[1,0].plot(T_m, u_p_values,'b') #Para
axis[1,0].plot(T_m, u_o_values,'r') #Ortho
axis[1,0].set_title('Chemical potential')
axis[1,0].set_xlabel('temperature [K]')
axis[1,0].set_ylabel('chemical potential [kJ/mol]')

axis[1,1].plot(T_m, phi_p_values, 'b') #Para
axis[1,1].plot(T_m, phi_o_values, 'r') #Ortho
axis[1,1].set_title('fugacity coefficient')
axis[1,1].set_xlabel('temperature [K]')
axis[1,1].set_ylabel('fugacity coefficient')

axis[0,1].plot(T_m, fug_p_values, 'b') #Para
axis[0,1].plot(T_m, fug_o_values, 'r') #Ortho
axis[0,1].set_title('fugacities')
axis[0,1].set_xlabel('temperature [K]')
axis[0,1].set_ylabel('fugacities')

plt.tight_layout()
plt.show()









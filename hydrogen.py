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
from pykingas.MieKinGas import MieKinGas



x = z = n = [0.25,0.75] # Total Molar composition

comps= 'P-H2,O-H2'

eos = saftvrqmie(comps, minimum_temperature=20)
srk = cubic(comps, 'SRK')
mie = MieKinGas(comps, use_eos=eos)

#get envelope
T_m, p_m = eos.get_envelope_twophase(1e5, z)
T_srk, p_srk = srk.get_envelope_twophase(1e5, n)
#plt.plot(p, T,'b', label="Pure hydrogen") # Tp-projection of phase envelope

figure, axis = plt.subplots(2, 2)
axis[0,0].plot(T_srk, p_srk, 'b', label="Ortho/para hydrogen, SRK")
axis[0,0].plot(T_m, p_m, 'y', label="O-H2 P-H2 mix, Saftvrqmie")
axis[0,0].set_ylabel("pressure [Pa]")

axis[0,0].set_xlabel("Temperature [K]")
axis[0,0].set_title('Saftvrqmie vs. SRK')



def f(x_o):
     T  =  300
     calj = 4.1840
     h_p = 2023.1*calj #300K 
     h_o = 2040.87*calj
     s_p = 31.212*calj #300K
     s_o = 33.446*calj
  
     n_ = [1-x_o, x_o]
     lnphi, = eos.thermo(T, 1e5, n_, eos.VAPPH)  

     #phi_p = lnphi[0], # Replace 0 with the index of the relevant value
     #phi_o = lnphi[1],
     DeltaG_o = (h_p-h_o)-T*(s_p-s_o)
     R = 8.314 # kJ/molK
     
     eq = x_o/(1-x_o) - (lnphi[0]/lnphi[1])*m.exp(-DeltaG_o/(R*T))
     return eq

initial_guess = 0.04

# Use scipy.optimize.fsolve to find the root of the equation
result = scipy.optimize.root(f, initial_guess)
print("Composition ortho:", result.x)

f(0.75)

T_list = np.arange(20, 307, 18).tolist()

p_list = []
for i in range(len(T_list)): 
    p = 1e5
    p_list.append(p)


vg_values = []
#specific volume for the ortho/para mxture
for T, p in zip(T_list, p_list):
    Vg, = eos.specific_volume(T, p, z, eos.VAPPH)
    vg_values.append(Vg,)

table1 = zip(T_list, p_list, vg_values)
print(tabulate(table1, headers = ('T [K]', 'P [kPa]', 'Vg [L/kg]')))
print()


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

calj = 4.1840 # 1cal = 4.18400 joule



for T, p, Vg in zip(T_list, p_list, vg_values):
    R = 8.314          #J/kmol
    h_p = 2023.1*calj  #300K 
    h_o = 2040.87*calj
    s_p = 31.212*calj  #300K
    s_o = 33.446*calj
    eos.set_ideal_enthalpy_reference_value(1,h_p)
    eos.set_ideal_enthalpy_reference_value(2,h_o)
    eos.set_ideal_entropy_reference_value(1,s_p)
    eos.set_ideal_entropy_reference_value(2,s_o)

    u_, = eos.chemical_potential_tv(T, Vg, n)
    phi, = eos.thermo(T, 1e5, n, eos.VAPPH)
    fug, = eos.fugacity_tv(T, Vg, n)

    fug_p_values.append(fug[0],)
    fug_o_values.append(fug[1],)

    u_p_values.append(u_[0]/1000)
    u_o_values.append(u_[1]/1000)

    phi_p_values.append(m.exp(phi[0]),)
    phi_o_values.append(m.exp(phi[1]),)

table2 = zip(T_list, u_p_values, u_o_values, phi_p_values, phi_o_values, fug_p_values, fug_o_values)
print(tabulate(table2, headers = ('T [K]', 'u_ [kJ/mol](P)', 'u_ [kJ/mol](O)', 'fug. coeff(P)', 'fug. coeff(O)', 'fugacities(P)','fugacities(O)')))



cond_val = []
visc_val = []
D_val = []
D_con_val = []
TD_fac = []
TD_val = []


for T, p, Vg in zip(T_list, p_list, vg_values):

    cond = mie.thermal_conductivity(T, Vg, x, N=2) # Thermal conductivity [W / m K]
    visc = mie.viscosity(T, Vg, x, N=2) # Shear viscosity [Pa s] #originally (T, Vm, x, N=2)
    D = mie.interdiffusion(T, Vg, x, N=2) #Binary diffusion coefficient [m^2 / s]
    D_CoN = mie.interdiffusion(T, Vg, x, N=2, frame_of_reference='CoN') # Diffusion coefficient
    alpha = mie.thermal_diffusion_factor(T, Vg, x, N=2) # Thermal diffusion factors [dimensionless]
    DT = mie.thermal_diffusion_coeff(T, Vg, x, N=2) # Thermal diffusion coefficients in the CoN FoR [mol / m s]

    cond_val.append(cond)
    visc_val.append(visc)
    D_val.append(D)
    D_con_val.append(D_CoN)
    TD_fac.append(alpha[1],)
    TD_val.append(abs(DT[1]))

table3 = zip(T_list, p_list, cond_val, visc_val, D_val, D_con_val)
print(tabulate(table3, headers = ('T [K]','P [kPa]','Therm. cond [W/mK]', 'Visc. [Pa S]', 'Diff. coeff [m^2/s]', 'Alt. Diff coeff [m^2/s]')))

table4 = zip(T_list, p_list, TD_fac, TD_val)
print(tabulate(table4, headers= ('T [K]', 'P [kPa]', 'Therm. Diff. Coeff. Fac','Therm. Diff. Coeff')))


#PLOTTING CHEMICAL POTENTIAL
axis[1,0].plot(T_list, u_p_values,'b') #Para
axis[1,0].plot(T_list, u_o_values,'r') #Ortho
axis[1,0].set_title('Chemical potential')
axis[1,0].set_xlabel('temperature [K]')
axis[1,0].set_ylabel('chemical potential [kJ/mol]')

#PLOTTING FUGACITY COEFFICIENT
axis[1,1].plot(T_list, phi_p_values, 'b') #Para
axis[1,1].plot(T_list, phi_o_values, 'r') #Ortho
axis[1,1].set_title('fugacity coefficient')
axis[1,1].set_xlabel('temperature [K]')
axis[1,1].set_ylabel('fugacity coefficient')

#PLOTTING FUGACITY
axis[0,1].plot(T_list, fug_p_values, 'b') #Para
axis[0,1].plot(T_list, fug_o_values, 'r') #Ortho
axis[0,1].set_title('fugacities')
axis[0,1].set_xlabel('temperature [K]')
axis[0,1].set_ylabel('fugacities')

plt.tight_layout()
plt.show()








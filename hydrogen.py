import numpy as np
import math as m
import scipy.integrate
from scipy import optimize
import matplotlib.pyplot as plt
from thermopack.saftvrmie import saftvrmie
from thermopack.saftvrqmie import saftvrqmie
from thermopack.cubic import cubic
#from thermopack.cubic import cubic
#from thermopack.cpa import cpa
import pandas as pd
from tabulate import tabulate 
from pykingas.MieKinGas import MieKinGas
from scipy.interpolate import interp1d


data = pd.read_excel(r'Para_percent.xlsx')

t_data, x_para, x_ortho = np.array(data['T[K]'], dtype=float), np.array(data['x_p'], dtype=float), np.array(data['x_o'], dtype=float)

t_x_data = zip(t_data, x_para, x_ortho)
print(tabulate(t_x_data, headers= ('T [K]', 'x (P)', 'x (O)')))

#INTERPOLATING DATA FOR EVERY 5 DEGREES
interp_tx = interp1d(t_data, x_para, kind='linear')
x_p_interp = interp_tx(np.linspace(20, 300, 10))

n_list = []

for i in range(len(x_p_interp)): 
    xp = x_p_interp[i]
    xo = 1-x_p_interp[i]
    n = xp,xo
    n_list.append(n)


x = [0.5,0.5]
comps= 'P-H2,O-H2'


eos = saftvrqmie(comps, minimum_temperature=20)
srk = cubic(comps, 'SRK')
mie = MieKinGas(comps, use_eos=eos)

#get envelope
T_m, p_m = eos.get_envelope_twophase(1e5, x)
T_srk, p_srk = srk.get_envelope_twophase(1e5, x)
#plt.plot(p, T,'b', label="Pure hydrogen") # Tp-projection of phase envelope
T_ = list(range(30,301))

figure, axis = plt.subplots(2, 2)
axis[0,0].plot(T_srk, p_srk, 'b', label="Ortho/para hydrogen, SRK")
axis[0,0].plot(T_m, p_m, 'y', label="O-H2 P-H2 mix, Saftvrqmie")
axis[0,0].set_ylabel("pressure [Pa]")

axis[0,0].set_xlabel("Temperature [K]")
axis[0,0].set_title('Saftvrqmie vs. SRK')


"""
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

#Use scipy.optimize.fsolve to find the root of the equation
result = scipy.optimize.root(f, initial_guess)
print("Composition ortho:", result.x)

f(0.75)
"""

#CREATING ARRAY OF TEMPERATURE AND PRESSURE
T_list = np.linspace(20,300,len(x_p_interp))

p_list = []
for i in range(len(T_list)): 
    p = 1e5
    p_list.append(p)


vg_values = []
#SPECIFIC VOLUME FOR THE MIXTURE
for T, p, n in zip(T_list, p_list, n_list):
    Vg, = eos.specific_volume(T, p, n, eos.VAPPH)
    vg_values.append(Vg,)

table1 = zip(T_list, p_list, vg_values)
print(tabulate(table1, headers = ('T [K]', 'P [kPa]', 'Vg [L/kg]')))
print()


#CALCULATING CHEMICAL POTENTIAL, FUGACITY COEFFICIENT AND FUGACITY FOR THE COMPONENTS IN THE MIXTURE
u_p_values = []
u_o_values =[]

phi_p_values = []
phi_o_values = []

fug_o_values = []
fug_p_values = []

calj = 4.1840 # 1cal = 4.18400 joule


for n, T, p, Vg in zip(n_list, T_list, p_list, vg_values):
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

#TABLE OF CHEMICAL POTENTIAL, FUGACITY COEFF. AND FUGACITY FOR THE COMPONENTS IN THE MIXTURE
table2 = zip(T_list, u_p_values, u_o_values, phi_p_values, phi_o_values, fug_p_values, fug_o_values)
print(tabulate(table2, headers = ('T [K]', 'u_ [kJ/mol](P)', 'u_ [kJ/mol](O)', 'fug. coeff(P)', 'fug. coeff(O)', 'fugacities(P)','fugacities(O)')))


#CALCULATING THERMAL CONDUCTIVITY, VISCOSITY, DIFFUSION COEFFICIENTS, THERMAL DIFFUSION FACTORS AND THERMAL DIFFUSION COEFFICIENTS FOR THE MIXTURE
cond_val = []
visc_val = []
D_val = []
D_con_val = []
TD_fac = []
TD_val = []


for T, p, Vg in zip(T_list, p_list, vg_values):

    cond = mie.thermal_conductivity(T, Vg, n, N=2) # Thermal conductivity [W / m K]
    visc = mie.viscosity(T, Vg, n, N=2) # Shear viscosity [Pa s] #originally (T, Vm, x, N=2)
    D = mie.interdiffusion(T, Vg, n, N=2) #Binary diffusion coefficient [m^2 / s]
    D_CoN = mie.interdiffusion(T, Vg, n, N=2, frame_of_reference='CoN') # Diffusion coefficient
    alpha = mie.thermal_diffusion_factor(T, Vg, n, N=2) # Thermal diffusion factors [dimensionless]
    DT = mie.thermal_diffusion_coeff(T, Vg, n, N=2) # Thermal diffusion coefficients in the CoN FoR [mol / m s]

    cond_val.append(cond)
    visc_val.append(visc)
    D_val.append(D)
    D_con_val.append(D_CoN)
    TD_fac.append(alpha[1],)
    TD_val.append(abs(DT[1]))

TD_fac1 = []
for i in TD_fac: 
    TD_fac1.append(i[0])

#TABLE OF THERM. COND. VISCOSITY, DIFFUSION COEFF. AND ALT. DIFFUSION COEFF, THERM. DIFF. COEFF. FAC., THERM. DIFF. COEFF.
table3 = zip(T_list, p_list, cond_val, visc_val, D_val, D_con_val, TD_fac1, TD_val)
print(tabulate(table3, headers = ('T [K]','P [kPa]','Therm. cond [W/mK]', 'Visc. [Pa S]', 'Diff. coeff [m^2/s]', 'Alt. Diff coeff [m^2/s]','Therm. Diff. Coeff. Fac','Therm. Diff. Coeff')))


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

#plt.tight_layout()
#plt.show()


#PLOTTING 
fig, axs = plt.subplots(3, 2)
axs[0,0].plot(T_list, cond_val, 'b') 
#axis[0,1].plot(T_list, fug_o_values, 'r') #Ortho
axs[0,0].set_title('Therm. cond.')
#axis[0,1].set_xlabel('temperature [K]')
#axis[0,1].set_ylabel('fugacities')

axs[0,1].plot(T_list, visc_val, 'b') 
axs[0,1].set_title('Viscosity')

axs[1,0].plot(T_list, D_val, 'r') 
axs[1,0].plot(T_list, D_con_val, 'b')
axs[1,0].set_title('Diff. coeff')

axs[1,1].plot(T_list, TD_fac1, 'b') 
axs[1,1].set_title('Therm. Diff. Coeff. Fac.')

axs[2,0].plot(T_list, TD_val, 'b') 
axs[2,0].set_title('Therm. Diff. Coeff.')

#axs[2,1].plot(T_list, visc_val, 'b') 
#axs[2,1].set_title('Viscosity')

plt.tight_layout()
plt.show()
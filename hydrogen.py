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
from matplotlib import figure

fig = plt.figure()

data = pd.read_excel(r'Para_percent.xlsx')

t_data, x_para, x_ortho = np.array(data['T[K]'], dtype=float), np.array(data['x_p'], dtype=float), np.array(data['x_o'], dtype=float)

t_x_data = zip(t_data, x_para, x_ortho)
print(tabulate(t_x_data, headers= ('T [K]', 'x (P)', 'x (O)')))

#PLOTTING VISCOSITY DATA
data2 = pd.read_excel(r'visc_data.xlsx')

t_visc, visc_data = np.array(data2['T[K]'], dtype=float), np.array(data2['n'], dtype=float)

t_visc_data = zip(t_visc, visc_data)
print(tabulate(t_visc_data, headers= ('T[K]', 'Visc. coeff')))


#INTERPOLATING DATA FOR EVERY 5 DEGREES
interp_tx = interp1d(t_data, x_para, kind='linear')
x_p_interp = interp_tx(np.linspace(20, 300, 10))

n_list = []

for i in range(len(x_p_interp)): 
    xp = x_p_interp[i]
    xo = 1-x_p_interp[i]
    n = xp,xo
    n_list.append(n)

x = [0.25,0.75]
comps= 'P-H2,O-H2'


eos = saftvrqmie(comps, minimum_temperature=20)
#eos = cubic(comps,'SRK')
srk = cubic(comps, 'SRK')
mie = MieKinGas(comps, use_eos=eos)

#get envelope
T_m, p_m = eos.get_envelope_twophase(1e5, x)
T_srk, p_srk = srk.get_envelope_twophase(1e5, x)
T_ = list(range(30,301))

figsrk = fig.add_subplot(211)
figsrk.plot(T_srk, p_srk, 'b', label="Ortho/para hydrogen, SRK")
figsrk.plot(T_m, p_m, 'y', label="O-H2 P-H2 mix, Saftvrqmie")
figsrk.set_ylabel("pressure [Pa]")
figsrk.set_xlabel("Temperature [K]")
figsrk.set_title('Saftvrqmie vs. SRK')


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


#PLOTTING COMPOSITION
figtx = fig.add_subplot(212)
figtx.plot(T_list, x_p_interp, 'b', label='Para')
figtx.plot(T_list, 1-x_p_interp, 'r', label='Ortho')
figtx.set_xlabel('Temperature [K]')
figtx.set_ylabel('Fraction')


vg_values = []
#SPECIFIC VOLUME FOR THE MIXTURE
for T, p, n in zip(T_list, p_list, n_list):
    Vg, = eos.specific_volume(T, p, n, eos.VAPPH)
    vg_values.append(Vg,)

table1 = zip(T_list, p_list, vg_values)
print(tabulate(table1, headers = ('T [K]', 'P [Pa]', 'Vg [L/kg]')))
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


for T, p, Vg, n in zip(T_list, p_list, vg_values, n_list):

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
    TD_fac.append(alpha[0],)
    TD_val.append(abs(DT[1]))


TD_fac1 = []
for alpha in TD_fac: 
    TD_fac1.append(alpha[1])

#TABLE OF THERM. COND. VISCOSITY, DIFFUSION COEFF. AND ALT. DIFFUSION COEFF, THERM. DIFF. COEFF. FAC., THERM. DIFF. COEFF.
table3 = zip(T_list, p_list, cond_val, visc_val, D_val, D_con_val, TD_fac1, TD_val)
print(tabulate(table3, headers = ('T [K]','P [kPa]','Therm. cond [W/mK]', 'Visc. [Pa S]', 'Diff. coeff [m^2/s]', 'Alt. Diff coeff [m^2/s]','Therm. Diff. Coeff. Fac','Therm. Diff. Coeff')))

interdif_list = []
for T, Vg, n in zip(T_list, vg_values, n_list):
    D_CoN = mie.interdiffusion(T, Vg, n, N=2, frame_of_reference='solvent', solvent_idx=1) # Diffusion coefficient
    interdif_list.append(D_CoN)

Lqq_list =[]
def Lqq(therm_cond, absT):
    Lqq = therm_cond*(absT**2)
    return Lqq_list.append(Lqq)

for cond, T in zip(cond_val, T_list):
    Lqq(cond, T)

L11_list = []
def L11(temp, diff, interdif):
    L11 = temp*diff/interdif
    return L11_list.append(L11)

for T, D, D_CoN in zip(T_list, D_val, interdif_list):
    L11(T, D, D_CoN)

L1q_list = []
def L1q(x1, x2, alphaT, L11, interdif):
    L1q = x1*x2*alphaT*L11*interdif
    return L1q_list.append(L1q)

for xp, alpha, L11, interdif in zip(x_p_interp, TD_fac1, L11_list, interdif_list):
    L1q(xp, 1-xp, alpha, L11, interdif)


table4 = zip(T_list, cond_val, Lqq_list, L11_list, L1q_list)
print(tabulate(table4, headers =('T [K]','Therm. cond [W/mK]', 'Lqq', 'L11', 'L1q' )))

#------------------------------------------------------------------------------------------------
#SAME FOR NORMAL HYDROGEN!!!!!!

vg_values_normal = []
#SPECIFIC VOLUME FOR THE MIXTURE
for T, p in zip(T_list, p_list):
    Vg, = eos.specific_volume(T, p, x, eos.VAPPH)
    vg_values_normal.append(Vg,)

#CALCULATING CHEMICAL POTENTIAL, FUGACITY COEFFICIENT AND FUGACITY FOR THE COMPONENTS IN THE MIXTURE
u_p_values_normal = []
u_o_values_normal =[]

phi_p_values_normal = []
phi_o_values_normal = []

fug_o_values_normal = []
fug_p_values_normal = []

calj = 4.1840 # 1cal = 4.18400 joule


for T, p, Vg in zip(T_list, p_list, vg_values_normal):
    R = 8.314          #J/kmol
    h_p = 2023.1*calj  #300K 
    h_o = 2040.87*calj
    s_p = 31.212*calj  #300K
    s_o = 33.446*calj
    eos.set_ideal_enthalpy_reference_value(1,h_p)
    eos.set_ideal_enthalpy_reference_value(2,h_o)
    eos.set_ideal_entropy_reference_value(1,s_p)
    eos.set_ideal_entropy_reference_value(2,s_o)

    u_, = eos.chemical_potential_tv(T, Vg, x)
    phi, = eos.thermo(T, 1e5, x, eos.VAPPH)
    fug, = eos.fugacity_tv(T, Vg, x)

    fug_p_values_normal.append(fug[0],)
    fug_o_values_normal.append(fug[1],)

    u_p_values_normal.append(u_[0]/1000)
    u_o_values_normal.append(u_[1]/1000)

    phi_p_values_normal.append(m.exp(phi[0]),)
    phi_o_values_normal.append(m.exp(phi[1]),)

#CALCULATING THERMAL CONDUCTIVITY, VISCOSITY, DIFFUSION COEFFICIENTS, THERMAL DIFFUSION FACTORS AND THERMAL DIFFUSION COEFFICIENTS FOR THE MIXTURE
cond_val_normal = []
visc_val_normal = []
D_val_normal = []
D_con_val_normal = []
TD_fac_normal = []
TD_val_normal = []


for T, p, Vg in zip(T_list, p_list, vg_values_normal):

    cond = mie.thermal_conductivity(T, Vg, x, N=2) # Thermal conductivity [W / m K]
    visc = mie.viscosity(T, Vg, x, N=2) # Shear viscosity [Pa s] #originally (T, Vm, x, N=2)
    D = mie.interdiffusion(T, Vg, x, N=2) #Binary diffusion coefficient [m^2 / s]
    D_CoN = mie.interdiffusion(T, Vg, x, N=2, frame_of_reference='CoN') # Diffusion coefficient
    alpha = mie.thermal_diffusion_factor(T, Vg, x, N=2) # Thermal diffusion factors [dimensionless]
    DT = mie.thermal_diffusion_coeff(T, Vg, x, N=2) # Thermal diffusion coefficients in the CoN FoR [mol / m s]

    cond_val_normal.append(cond)
    visc_val_normal.append(visc)
    D_val_normal.append(D)
    D_con_val_normal.append(D_CoN)
    TD_fac_normal.append(alpha[0],)
    TD_val_normal.append(abs(DT[1]))


TD_fac1_normal = []
for alpha in TD_fac_normal: 
    TD_fac1_normal.append(alpha[1])

#TABLE OF THERM. COND. VISCOSITY, DIFFUSION COEFF. AND ALT. DIFFUSION COEFF, THERM. DIFF. COEFF. FAC., THERM. DIFF. COEFF.
#table3 = zip(T_list, p_list, cond_val, visc_val, D_val, D_con_val, TD_fac1, TD_val)
#print(tabulate(table3, headers = ('T [K]','P [kPa]','Therm. cond [W/mK]', 'Visc. [Pa S]', 'Diff. coeff [m^2/s]', 'Alt. Diff coeff [m^2/s]','Therm. Diff. Coeff. Fac','Therm. Diff. Coeff')))

interdif_list_normal = []
for T, Vg in zip(T_list, vg_values_normal):
    D_CoN = mie.interdiffusion(T, Vg, x, N=2, frame_of_reference='solvent', solvent_idx=1) # Diffusion coefficient
    interdif_list_normal.append(D_CoN)

Lqq_list_normal =[]
def Lqq(therm_cond, absT):
    Lqq = therm_cond*(absT**2)
    return Lqq_list_normal.append(Lqq)

for cond, T in zip(cond_val_normal, T_list):
    Lqq(cond, T)

L11_list_normal = []
def L11(temp, diff, interdif):
    L11 = temp*diff/interdif
    return L11_list_normal.append(L11)

for T, D, D_CoN in zip(T_list, D_val_normal, interdif_list_normal):
    L11(T, D, D_CoN)

L1q_list_normal = []
def L1q(x1, x2, alphaT, L11, interdif):
    L1q = x1*x2*alphaT*L11*interdif
    return L1q_list_normal.append(L1q)

x_para_normal = 0.25
for alpha, L11, interdif in zip(TD_fac1_normal, L11_list_normal, interdif_list_normal):
    L1q(x_para_normal, 1-x_para_normal, alpha, L11, interdif)

#---------------------------------------------------------------------------------------------------


#PLOTTING CHEMICAL POTENTIAL
figf, axis = plt.subplots(3)
axis[0].plot(T_list, u_p_values,'g', label= 'Para') #Para
axis[0].plot(T_list, u_o_values,'m', label= 'Ortho') #Ortho
axis[0].plot(T_list, u_p_values_normal,'b', label= 'Para') #Para
axis[0].plot(T_list, u_o_values_normal,'r', label= 'Ortho') #Ortho
axis[0].set_title('Chemical potential')
axis[0].set_xlabel('Temperature [K]')
axis[0].set_ylabel('[kJ/mol]')

#PLOTTING FUGACITY COEFFICIENT
axis[1].plot(T_list, phi_p_values, 'g', label= 'Para') #Para
axis[1].plot(T_list, phi_o_values, 'm', label= 'Ortho') #Ortho
axis[1].plot(T_list, phi_p_values_normal, 'b', label= 'Para') #Para
axis[1].plot(T_list, phi_o_values_normal, 'r', label= 'Ortho') #Ortho
axis[1].set_title('Fugacity coefficient')
axis[1].set_xlabel('Temperature [K]')
#axis[1].set_ylabel('Fugacity coefficient')

#PLOTTING FUGACITY
axis[2].plot(T_list, fug_p_values, 'g', label= 'Para') #Para
axis[2].plot(T_list, fug_o_values, 'm', label= 'Ortho') #Ortho
axis[2].plot(T_list, fug_p_values_normal, 'b', label= 'Para') #Para
axis[2].plot(T_list, fug_o_values_normal, 'r', label= 'Ortho') #Ortho
axis[2].set_title('Fugacities')
axis[2].set_xlabel('Temperature [K]')
#axis[2].set_ylabel('fugacities')

#print(tabulate(table3, headers = ('T [K]','P [kPa]','Therm. cond [W/mK]', 'Visc. [Pa S]', 'Diff. coeff [m^2/s]', 'Alt. Diff coeff [m^2/s]','Therm. Diff. Coeff. Fac','Therm. Diff. Coeff')))

#PLOTTING 
fig, axs = plt.subplots(3, 2)
axs[0,0].plot(T_list, vg_values, 'r')
axs[0,0].plot(T_list, vg_values_normal, 'b')
axs[0,0].set_title('Specific volume')
#axis[0,0].set_xlabel('Temperature [K]')
axs[0,0].set_ylabel('[m^3/mol]')

axs[2,1].plot(T_list, cond_val, 'r') 
axs[2,1].plot(T_list, cond_val_normal, 'b') 
axs[2,1].set_title('Therm. cond.')
#axis[2,1].set_xlabel('Temperature [K]')
axs[2,1].set_ylabel('[W/mK]')

axs[0,1].plot(T_list, visc_val, 'r') 
axs[0,1].plot(T_list, visc_val_normal, 'b') 
axs[0,1].set_title('Viscosity')
#axis[0,1].set_xlabel('Temperature [K]')
axs[0,1].set_ylabel('[Pa s]')

axs[1,0].plot(T_list, D_val, 'm') 
axs[1,0].plot(T_list, D_con_val, 'r')
axs[1,0].plot(T_list, D_val_normal, 'g') 
axs[1,0].plot(T_list, D_con_val_normal, 'b')
axs[1,0].set_title('Diff. coeff')
#axis[1,0].set_xlabel('Temperature [K]')
axs[1,0].set_ylabel('[m^2/s]')

axs[1,1].plot(T_list, TD_fac1, 'r') 
axs[1,1].plot(T_list, TD_fac1_normal, 'b') 
axs[1,1].set_title('Therm. Diff. Coeff. Fac.')

axs[2,0].plot(T_list, TD_val, 'r') 
axs[2,0].plot(T_list, TD_val_normal, 'b') 
axs[2,0].set_title('Therm. Diff. Coeff.')



#PLOTTING Lqq AND Lq1
figt, axst = plt.subplots(3)
#figt.suptitle('Transport coefficients')
axst[0].plot(T_list, Lqq_list, 'r',label='accurate')
axst[0].plot(T_list, Lqq_list_normal, 'b', label='normal hydrogen')
axst[0].set_title('Transport coefficient (Lqq)')

axst[1].plot(T_list, L11_list, 'r')
axst[1].plot(T_list, L11_list_normal, 'b')
axst[1].set_title('Transport coefficient (L11)')

axst[2].plot(T_list, L1q_list, 'r')
axst[2].plot(T_list, L1q_list_normal, 'b')
axst[2].set_title('Transport coefficient (L1q)')

figvs, axvs = plt.subplots(2)
axvs[0].plot(T_list, visc_val, 'r') 
axvs[0].plot(T_list, visc_val_normal, 'b') 
axvs[0].set_title('Viscosity')
axvs[1].plot(T_list, visc_val, 'r') 
axvs[1].plot(T_list, visc_val_normal, 'b') 
axvs[1].set_title('Viscosity')
plt.scatter(t_visc, visc_data, marker='o')
axvs[1].set_xlabel('Temperature [K]')
axvs[1].set_ylabel('[Pa s]')

#plt.plot(T_list, visc_val, 'r')
#plt.plot(T_list, visc_val_normal, 'b')


plt.tight_layout()
plt.show()
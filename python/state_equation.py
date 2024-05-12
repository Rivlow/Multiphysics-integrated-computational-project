import numpy as np
import matplotlib.pyplot as plt
import os
import sys

SMALL_SIZE = 8
MEDIUM_SIZE = 14
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)   # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)   # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)   # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('text', usetex=True)
plt.rc('font', family='lmodern')


c_0 = 30
rho_0 = 1000
R = 8.314
T = 273.15 + 25
gamma = 7
B = np.power(c_0,2)*rho_0/gamma
M = 18e-3


rho = np.linspace(1, 10000, 100)
p= np.zeros_like(rho)
c= np.zeros_like(rho)

#state_equation = "Ideal_gas_law"
state_equation = "Quasi_incompressible_fluid"
current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))



if (state_equation == "Ideal_gas_law"):
    for i in range(len(rho)):
        p[i] = (R*T/M)*(rho[i]/rho_0 - 1)
        c[i] = c_0
elif(state_equation == "Quasi_incompressible_fluid"):
    for i in range(len(rho)):
        p[i] = B*np.power(rho[i]/rho_0, 1/gamma)
        c[i] = c_0*np.power(rho[i]/rho_0, gamma-1)
    
fig, ax = plt.subplots(1, 2, figsize=(6.9, 4))
ax[0].plot(rho, p, label = f"{state_equation}")
ax[0].set_xlabel(r"Density $\rho$ [kg/$m^3$]")
ax[0].set_ylabel(r"Pressure p($\rho$) [N/$m^2$] ")
ax[0].grid(True, alpha=0.5)

ax[1].plot(rho, c, label = f"{state_equation}")
ax[1].set_xlabel(r"Density $\rho$ [kg/$m^3$]")
ax[1].set_ylabel(r"Speed of sound c($\rho$) [m/s]")
ax[1].grid(True, alpha=0.5)

plt.tight_layout()
#plt.grid(True)
plt.savefig(f"{current_directory}\Pictures\state_equation_{state_equation}.PDF")
plt.show()

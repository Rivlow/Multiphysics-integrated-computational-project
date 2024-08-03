import numpy as np
import matplotlib.pyplot as plt
import os
import sys
current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))

def isLatex(latex):
    if latex:
        
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

def plotStateEquation(plot, save):

    for i in range(len(rho)):
        p_ig[i] = (R*T/M)*(rho[i]/rho_0 - 1)
        c_ig[i] = c_0
        p_qi[i] = B*(np.power(rho[i]/rho_0, gamma) - 1)
        c_qi[i] = c_0*np.power(rho[i]/rho_0, (gamma-1)/2)

            
    if plot["IG"]: 

        fig, ax = plt.subplots(1, 2, figsize=(6.9, 4))
        ax[0].plot(rho, p_ig, label = "Ideal gas law")
        ax[0].set_xlabel(r"Density $\rho$ [kg/$m^3$]")
        ax[0].set_ylabel(r"Pressure p($\rho$) [N/$m^2$] ")
        ax[0].grid(True, alpha=0.5)

        ax[1].plot(rho, c_ig, label = "Ideal gas law")
        ax[1].set_xlabel(r"Density $\rho$ [kg/$m^3$]")
        ax[1].set_ylabel(r"Speed of sound c($\rho$) [m/s]")
        ax[1].grid(True, alpha=0.5)

        plt.tight_layout()
        plt.grid(True)
        if save:
            plt.savefig(rf"{current_directory}\Pictures\state_equation_IG.PDF")
        plt.show()

    if plot["QI"]: 

        fig, ax = plt.subplots(1, 2, figsize=(6.9, 4))
        ax[0].plot(rho, p_qi, label = "Quasi incompressible fluid")
        ax[0].set_xlabel(r"Density $\rho$ [kg/$m^3$]")
        ax[0].set_ylabel(r"Pressure p($\rho$) [N/$m^2$] ")
        ax[0].grid(True, alpha=0.5)

        ax[1].plot(rho, c_qi, label = "Quasi incompressible fluid")
        ax[1].set_xlabel(r"Density $\rho$ [kg/$m^3$]")
        ax[1].set_ylabel(r"Speed of sound c($\rho$) [m/s]")
        ax[1].grid(True, alpha=0.5)
        
        plt.tight_layout()
        plt.grid(True)
        if save:
            plt.savefig(rf"{current_directory}\Pictures\state_equation_IG.PDF")
        plt.show()

plot = {"IG":True, "QI":True}
save = False
isLatex(False)

# Thermodynamic parameters
c_0 = 30
rho_0 = 1000
R = 8.314
T = 273.15 + 25
gamma = 7
B = np.power(c_0,2)*rho_0/gamma
print(f"B = {B}")
M = 18e-3

rho = np.linspace(1000, 1200, 100)
p_ig = np.zeros_like(rho)
c_ig = np.zeros_like(rho)
p_qi = np.zeros_like(rho)
c_qi = np.zeros_like(rho)

plotStateEquation(plot, save)







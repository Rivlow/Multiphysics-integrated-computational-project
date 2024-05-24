import numpy as np
import matplotlib.pyplot as plt
import os 
import sys

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



def plot_state_equation(state):
    
    
    p_qi = B *(np.power(rho/rho_0, gamma) - 1)
    c_qi = c = c_0 * np.power(rho / rho_0, 0.5*(gamma - 1))

    c_ig = c_0*np.ones(len(rho))
    p_ig = rho*R*T/M
    
    if state == "Quasi_incompressible":
        
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

        ax1.plot(rho, p_qi)
        ax1.set_ylabel("Pressure [Pa]")
        ax1.grid(True)

        ax2.plot(rho, c_qi)
        ax2.set_xlabel(r"Density [kg/$m^3$]")
        ax2.set_ylabel("Speed of Sound [m/s]")
        ax2.grid(True)
        
        plt.tight_layout()
        plt.savefig(f"{current_directory}/Pictures/{state}.PDF")
        plt.show()
        
    if state == "Ideal_gas_law":
        
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

        ax1.plot(rho, p_ig)
        ax1.set_ylabel("Pressure [Pa]")
        ax1.grid(True)

        ax2.plot(rho, c_ig)
        ax2.set_xlabel(r"Density [kg/$m^3$]")
        ax2.set_ylabel("Speed of Sound [m/s]")
        ax2.grid(True)
        
        plt.tight_layout()
        plt.savefig(f"{current_directory}/Pictures/{state}.PDF")
        plt.show()


# Constantes
c_0 = 30
rho_0 = 1000
gamma = 7  # Exemple de valeur pour gamma
B = (c_0**2)*rho_0/gamma
T = 273.15+25
M = 18e-3
R = 8.314
rho = np.linspace(0, 2000, 100)

state = "Quasi_incompressible"
#state = "Ideal_gas_law"
current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
latex = True

isLatex(latex)
plot_state_equation(state)




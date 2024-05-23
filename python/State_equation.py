import numpy as np
import matplotlib.pyplot as plt

# Constantes
c_0 = 30
rho_0 = 1000
gamma = 7  # Exemple de valeur pour gamma
B = (c_0**2)*rho_0/gamma
T = 273.15+25
M = 18e-3
R = 8.314
print(B)

def plot_state_equation(state):
    
    rho = np.linspace(960, 1040, 100)
    p_qi = B * ((rho / rho_0) ** gamma - 1)
    c_qi = c = c_0 * np.power(rho / rho_0, 0.5*(gamma - 1))

    c_ig = c_0
    p_ig = rho*R*T/M
    if state["Quasi incompressible"]:
        
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

        # Premier plot : p_qi
        ax1.plot(rho, p_qi)
        ax1.set_ylabel("Pressure [Pa]")
        ax1.grid(True)

        # Deuxième plot : c_qi
        ax2.plot(rho, c_qi)
        ax2.set_xlabel(r"Density [kg/$m^3$]")
        ax2.set_ylabel("Speed of Sound [m/s]")
        ax2.grid(True)
        
        plt.tight_layout()
        plt.show()


state = {"Quasi incompressible":True, "Ideal gas law":False}
plot_state_equation(state)
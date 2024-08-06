import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import sys
import pandas as pd



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

        
def StateEquation():
    

    # Arbitrary take last iteration (seeking steady-state)
    p = np.array(pd.read_csv("output/500_pressure.csv", sep = ',', decimal='.', header=None))[-1] 
    rho = np.array(pd.read_csv("output/500_rho.csv", sep = ',', decimal='.', header=None))[-1]

    c_0 = 30
    rho_0 = 1000
    gamma = 7  
    B = (c_0**2)*rho_0/gamma
    T = 273.15+25
    M = 18e-3
    R = 8.314

    plt.scatter(p, rho, s = 10, label = "state_equation_sph")

    p_compare = np.linspace(48.7799, 12397, 100)
    rho_compare = rho_0 * ((p_compare + 1) / B) ** (1 / gamma)
    #plt.xscale("log")
    #plt.plot(p_val[11:], rho_val[11:], color = 'b', label = 'SPH')
    #plt.plot(p_compare, rho_compare, color = 'r', label = 'quasi incompressible (theorical)')
    plt.xlabel("Pressure [Pa]")
    plt.ylabel(r"Density [kg/$m^3$]")
    plt.grid(True)
    plt.legend(loc = "best")
    
    
    
def Hydrostatic():


    p = np.array(pd.read_csv("output/p.csv", sep = ',', decimal='.', header=None))[-1]
    rho_ref = 1023 # average density for ghost particles (using "plotSingleVariable")
    g = 9.81
    z = np.linspace(0.25, 0.9, len(p))
    p_th = rho_ref*g*z
    
    print(p_th[0])
    print(p[-1])

    plt.figure()
    plt.scatter(z, p[::-1], color = 'r', label = 'SPH pressure')
    plt.plot(z, p_th, color = 'b', label = 'Theoretical hydrostatic pressure')

    plt.ylabel("Pressure [Pa]")
    plt.xlabel(r"Height [m]")
    plt.grid(True)
    plt.legend(loc = "best")
    plt.tight_layout()
    #plt.savefig(f'{current_directory}/Pictures/hydrostatic_pressure.PDF')

def PoiseuilleFlow():


    u_x = np.array(pd.read_csv("output/u_x.csv", sep = ',', decimal='.', header=None))[-1]
    z = np.linspace(0.02, 1.24, len(u_x))

    u_sph= u_x + 0.5

    u_max = max(u_sph)
    h = z[-1] - z[0]
    u_analytic = u_max*((4/h)*z - (4/np.pow(h,2)*np.pow(z,2)))
    
    plt.figure()
    plt.scatter(z, u_sph, s = 10, label = "sph")
    plt.plot(z, u_sph, label = "sph", ls = "--")
    plt.plot(z+0.02, u_analytic, label = "analytic", ls = "--") # +0.04 because we begin at the ghost particle
    plt.legend(loc = "best")

        

def main():
    
    current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
    isLatex(False)

    # Put in comment what you don't want to plot

    #PoiseuilleFlow()
    #Hydrostatic()
    StateEquation()

    plt.show()
    

    

if __name__ == "__main__":
    main()




"""import vtk
#from vtk.util.numpy_support import vtk_to_numpy
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys
import pandas as pd

def fonction():
    file = "output/p.csv"
    ite = 99
    pressure = pd.read_csv(file, sep = ',', decimal='.', header=None)
    pressure = np.array(pressure)
    print(pressure[ite])
    print(len(pressure[ite]))
    x = np.arange(0.125, 1.875, 0.025)
    print(len(x))
    x1 = np.arange(0.1,1.85,0.025)

    plt.plot(x,pressure[ite][::-1])
    plt.plot(x,1000*9.81*x1)

    plt.show()

def mrua():
    file = "output/p.csv"
    pos = pd.read_csv("output/50_pos_z.csv",sep = '.', decimal=',', header=None)
    time = pd.read_csv("output/time.csv",sep = '.', decimal=',', header=None)
    pos =  np.transpose(np.asarray(pos))
    time =  np.transpose(np.asarray(time))
   
    print(time)
    plt.scatter(time[0],pos[0], color = "red", label = "particle")
    time = np.array(time)
    pos= np.array(pos)
    
    i = 0
    time_2 = np.zeros(len(time[0]))
    while(i<len(time_2)):
        time_2[i] = time[0][i] *   time[0][i]
        i =i + 1
    y = pos[0][0] - (9.81/2) * time_2

    
    plt.plot(time[0], y, label = "theorie")
    plt.legend()
    plt.xlabel("time")
    plt.ylabel("z position")
    plt.ylim(-0.01,1.3)
    plt.xlim(0.01,0.5)
    plt.show()

mrua() """
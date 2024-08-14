import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import sys
import pandas as pd

#fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(19, 7))
# plt.figure(figsize=(8.5, 5))
save = True
Lat = True
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
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/state_equation.PDF')
    plt.show()
    
    
    
def Hydrostatic():


    '''p = np.array(pd.read_csv("output/csv/Euler_3D_hydrostatic_p.csv", sep = ',', decimal='.', header=None))[-1]
    
    time =  np.array(pd.read_csv("output/csv/3D_hydrostatic_time.csv", decimal='.', header=None).astype(float)).T
    rho_ref = 1000 # average density for ghost particles (using "plotSingleVariable")
    g = 9.81
    z = np.linspace(0.05, 1.15, len(p))
    z_th = np.linspace(0,1.1,len(p))
    p_th = rho_ref*g*z_th
    print(time)
    print(p_th[0])
    print(p[::-1])

    plt.figure(figsize=(8.5, 5))
    plt.scatter(z, p[::-1], color = 'r', label = 'SPH pressure')
    plt.plot(z_th+0.28, p_th, color = 'b', label = 'Theoretical hydrostatic pressure')

    plt.ylabel("Pressure [Pa]")
    plt.xlabel(r"Height [m]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/hydrostatic.PDF')
    plt.show()'''
    p = np.array(pd.read_csv("output/csv/Euler_2D_hydrostatic_p.csv", sep = ',', decimal='.', header=None))[-1]
    
    time =  np.array(pd.read_csv("output/csv/2D_hydrostatic_time.csv", decimal='.', header=None).astype(float)).T
    rho_ref = 1000 # average density for ghost particles (using "plotSingleVariable")
    g = 9.81
    z = np.linspace(0.04, 0.82, len(p))
    z_th = np.linspace(0,0.8,len(p))
    p_th = rho_ref*g*z_th
   

    plt.figure(figsize=(8.5, 5))
    plt.scatter(z, p[::-1], color = 'r', label = 'SPH pressure')
    plt.plot(z_th+0.215, p_th, color = 'b', label = 'Theoretical hydrostatic pressure')

    plt.ylabel(r"Pressure $p(z)$ [Pa]")
    plt.xlabel(r"Depth $z$ [m]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/3D_hydrostatic.PDF')
    plt.show()

    s = 0.02
    print((z_th[-5]-0.215))
    plt.figure(figsize=(8.5, 5))
    p_part = np.array(pd.read_csv("output/csv/Euler_2D_hydrostatic_p.csv", sep = ',', decimal='.', header=None)).T[5]
    plt.axhline(y=1000*9.81*(z_th[-5]-0.215), color='r', linestyle='--')
    plt.plot(time[0], p_part, color = 'b', label = 'Theoretical hydrostatic pressure')

    plt.ylabel(r"Pressure $p(z)$ [Pa]")
    plt.xlabel(r"Depth $z$ [m]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/2D_hydrostatic.PDF')
    plt.show()

    

def PoiseuilleFlow():


    u_x = np.array(pd.read_csv("output/csv/Euler_2D_Poiseuille_u_x.csv", sep = ',', decimal='.', header=None))[200]
    z = np.linspace(0.05, 1.3, len(u_x))

    u_sph= u_x + 0.5

    u_max = max(u_sph)
    h = z[-1] - z[0]
    print(z)
    u_analytic = u_max*((4/h)*z - (4/np.power(h,2)*np.power(z,2)))
    
    plt.figure()
    plt.scatter(z[1:-1], u_sph[1:-1], s = 10, label = "sph")
    plt.plot(z[1:-1], u_sph[1:-1], label = "sph", ls = "--")
    plt.plot(z[0:-2]+0.05, u_analytic[0:-2], label = "analytic", ls = "--") # +0.04 because we begin at the ghost particle
    plt.legend(loc = "best")
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/Poiseuille.PDF')
    plt.show()

def dam_break(): 

    pos_x = np.array(pd.read_csv("output/csv/3D_dam_break_max_pos_x.csv", sep=',', decimal='.', header=None).astype(float))     
    time = np.array(pd.read_csv("output/csv/3D_dam_break_time.csv", sep=',', decimal='.', header=None).astype(float))

    L = 0.2
    s= 0.05
    time_ad = time * np.sqrt(2 * 9.81 / L)
    pos_ad = (pos_x -3*s/2)/ L

    plt.plot(time_ad, pos_ad , color = "red")
    plt.xlim(0,4)
    plt.ylabel(r"x/L [-]")
    plt.xlabel(r"t$\sqrt{2\,g/L}$ [-]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/3D_dam_break.PDF')
    plt.show()

    pos_x = np.array(pd.read_csv("output/csv/2D_dam_break_max_pos_x.csv", sep=',', decimal='.', header=None).astype(float))     
    time = np.array(pd.read_csv("output/csv/2D_dam_break_time.csv", sep=',', decimal='.', header=None).astype(float))

    L = 0.2
    s= 0.05
    time_ad = time * np.sqrt(2 * 9.81 / L)
    pos_ad = (pos_x -3*s/2)/ L

    plt.plot(time_ad, pos_ad , color = "red")
    plt.xlim(0,4)
    plt.ylabel(r"x/L [-]")
    plt.xlabel(r"t$\sqrt{2\,g/L}$ [-]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/2D_dam_break.PDF')
    plt.show()

def MRUA():
    pos_z = np.array(pd.read_csv("output/csv/3D_splash_1105_pos_z.csv", sep=',', decimal='.', header=None).astype(float)).T
    time =  np.array(pd.read_csv("output/csv/3D_splash_time.csv", decimal='.', header=None).astype(float)).T
    pos_z_RK2 = np.array(pd.read_csv("output/csv/3D_RK2_splash_1105_pos_z.csv", sep=',', decimal='.', header=None).astype(float)).T
    time_RK2 =  np.array(pd.read_csv("output/csv/3D_RK2_splash_time.csv", decimal='.', header=None).astype(float)).T


    print(len(pos_z[0]))
    print(len(time[0]))
    Z_th = pos_z[0][0] - 0.5*9.81*time[0]**2
    Z_th_RK2 = pos_z_RK2[0][0] - 0.5*9.81*time_RK2[0]**2
    
    print(Z_th)
    plt.figure(figsize=(8.5, 5))
    
    plt.plot(time[0], Z_th, "b", label = "Analytical solution", zorder = 1)
    plt.scatter(time[0], pos_z[0] , color = "red", marker = "o", label = "SPH simulation" , zorder = 2)
    plt.scatter(time_RK2[0], pos_z_RK2[0] , color = "green", marker = "o", label = "SPH simulation" , zorder = 2)
    plt.ylim(0,0.7)
    plt.xlim(0,0.5)
    plt.ylabel(r"Postion $z(t)$ [m]")
    plt.xlabel(r"Time $t$ [s]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/MRUA.PDF')
    plt.show()
    
    
    error_RK2= np.abs((pos_z_RK2[0] - Z_th_RK2)/Z_th_RK2)
    error = np.abs((pos_z[0] - Z_th)/Z_th)
    plt.plot(time[0], error)
    plt.plot(time_RK2[0],error_RK2)
    plt.yscale("log")
    plt.show()

def main():
    
    current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
    isLatex(Lat)

    # Put in comment what you don't want to plot

    #PoiseuilleFlow()
    Hydrostatic()
    #StateEquation()
    #dam_break()
    #MRUA()
    
    

    

if __name__ == "__main__":
    main()

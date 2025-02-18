import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import sys
import pandas as pd

#fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(19, 7))
# plt.figure(figsize=(8.5, 5))
save = True
Lat = False
def isLatex(latex):
    if latex:
        
        SMALL_SIZE = 8
        MEDIUM_SIZE = 14
        BIGGER_SIZE = 18
        plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
        plt.rc('axes', titlesize=BIGGER_SIZE)    # fontsize of the axes title
        plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
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
    plt.xlabel(r"Pressure [Pa]")
    plt.ylabel(r"Density [kg/$m^3$]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/state_equation.PDF')
    plt.show()
    
    
    
def Hydrostatic():

    #3D hydrostatic Euler

    rho_ref = 1000 
    g = 9.81
    dist = 0.155
    part = 10

    p = np.array(pd.read_csv("output/csv/Euler_3D_hydrostatic_p.csv", sep = ',', decimal='.', header=None))[-1]
    time =  np.array(pd.read_csv("output/csv/3D_hydrostatic_time.csv", decimal='.', header=None).astype(float)).T
    p_part = np.array(pd.read_csv("output/csv/Euler_3D_hydrostatic_p.csv", sep = ',', decimal='.', header=None)).T[part]
    z = np.linspace(0.05, 0.85, len(p))
    z_th = np.linspace(0,0.8,len(p))
  
    p_th = rho_ref*g*z_th
    plt.figure(figsize=(7, 5))
    
    plt.scatter(z-dist, p[::-1], color = 'r', label = 'SPH pressure')
    plt.plot(z_th, p_th, color = 'b', label = 'Theoretical hydrostatic pressure')
    plt.ylabel(r"Pressure $p(z)$ [Pa]")
    plt.xlabel(r"Depth $z$ [m]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/3D_hydrostatic.PDF')
    plt.show()

    
    plt.figure(figsize=(7, 5))

    plt.axhline(y=1000*9.81*(z_th[-part+1]-dist), color='blue', label = 'Theoretical hydrostatic pressure', zorder = 1)
    plt.plot(time[0], p_part, color = 'red', zorder = 2, label = "SPH pressure")
    
    plt.ylabel(r"Pressure $p(t)$ [Pa]")
    plt.xlabel(r"Time $t$ [s]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/3D_hydrostatic_convergence.PDF')
    plt.show()



    rho_ref = 1000 
    g = 9.81
    
    #2D hydrostatic Euler
    dist = 0.116
    part = 10

    p = np.array(pd.read_csv("output/csv/Euler_2D_hydrostatic_p.csv", sep = ',', decimal='.', header=None))[-1]
    time =  np.array(pd.read_csv("output/csv/2D_hydrostatic_time.csv", decimal='.', header=None).astype(float)).T
    p_part = np.array(pd.read_csv("output/csv/Euler_2D_hydrostatic_p.csv", sep = ',', decimal='.', header=None)).T[part]

    p_RK2 = np.array(pd.read_csv("output/csv/RK22_2D_RK2_hydrostatic_p.csv", sep = ',', decimal='.', header=None))[-1]
    time_RK2 =  np.array(pd.read_csv("output/csv/2D_RK2_hydrostatic_time.csv", decimal='.', header=None).astype(float)).T
    p_part_RK2 = np.array(pd.read_csv("output/csv/RK22_2D_RK2_hydrostatic_p.csv", sep = ',', decimal='.', header=None)).T[part]

    z = np.linspace(0.04, 0.82, len(p))
    
    z_th = np.linspace(0,0.78,len(p))
    
    p_th = rho_ref*g*z_th
    
    plt.figure(figsize=(7, 5))
    
    plt.scatter(z-dist, p[::-1], color = 'r',marker="o", label = 'SPH pressure Euler')
    plt.scatter(z-dist, p_RK2[::-1], color = 'g',marker = "+", label = 'SPH pressure RK22')
    plt.plot(z_th, p_th, color = 'b', label = 'Theoretical hydrostatic pressure')
    plt.ylabel(r"Pressure $p(z)$ [Pa]")
    plt.xlabel(r"Depth $z$ [m]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/2D_hydrostatic_and_RK2.PDF')
    plt.show()
    
    plt.figure(figsize=(7, 5))
    plt.axhline(y=1000*9.81*(z_th[-part+1]-dist), color='blue', label = 'Theoretical hydrostatic pressure', zorder = 1)
    plt.plot(time[0], p_part, color = 'red', zorder = 2, label = "SPH pressure Euler")
    plt.plot(time_RK2[0], p_part_RK2, color = 'g', linestyle="--", zorder = 2, label = "SPH pressure RK22")
    plt.ylabel(r"Pressure $p(t)$ [Pa]")
    plt.xlabel(r"Time $t$ [s]")
    
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/2D_hydrostatic_convergence_and_RK2.PDF')
    plt.show()

    


    

def PoiseuilleFlow():


    u_x = np.array(pd.read_csv("output/csv/Euler_2D_Poiseuille_u_x.csv", sep = ',', decimal='.', header=None))[-1]
    z = np.linspace(0.05, 1.3, len(u_x))

    u_sph= u_x + 0.5

    u_max = max(u_sph)
    h = z[-1] - z[0]
    print(z)
    u_analytic = u_max*((4/h)*z - (4/np.power(h,2)*np.power(z,2)))
    
    plt.figure(figsize=(8.5, 5))
    plt.scatter(z[1:-1], u_sph[1:-1], label = "sph", marker ="o", color = "red", zorder = 2)
    
    plt.plot(z[0:-2]+0.05, u_analytic[0:-2], label = "analytic", color="blue", zorder = 1) # +0.04 because we begin at the ghost particle
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/Poiseuille.PDF')
    plt.show()

def dam_break(): 

    pos_x = np.array(pd.read_csv("output/csv/3D_dam_break_max_pos_x.csv", sep=',', decimal='.', header=None).astype(float))     
    time = np.array(pd.read_csv("output/csv/3D_dam_break_time.csv", sep=',', decimal='.', header=None).astype(float))

    image_j_x = np.array(pd.read_csv("output/csv/Results.csv", sep=',', decimal='.', header=0).astype(float)).T[-2]
    image_j_y = np.array(pd.read_csv("output/csv/Results.csv", sep=',', decimal='.', header=0).astype(float)).T[-1]
    print(image_j_x)
    print(image_j_y)

    unity_x = np.abs(image_j_x[1]-image_j_x[0])
    unity_y = np.abs(image_j_y[2]-image_j_y[0])

    points_x = np.abs(image_j_x[3:]-image_j_x[0])
    points_y = np.abs(image_j_y[3:]-image_j_y[0])

    coord_x = points_x/unity_x
    coord_y = points_y/unity_y

    print(coord_x)
    print(coord_y)
    

    

    
    L = 0.5
    s= 0.01
    time_ad = time * np.sqrt(2 * 9.81 / L)
    pos_ad = (pos_x -3*s/2)/ L

    plt.plot(time_ad, pos_ad , color = "red", label= "SPH simulation", zorder = 1)
    plt.scatter(coord_x, coord_y+1, color = 'b', marker = "o", label = "Koshizuka and Oka 1996", zorder = 2)
    plt.xlim(0,3.3)
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

    pos_x_05 = np.array(pd.read_csv("output/csv/2D_dam_05_break_max_pos_x.csv", sep=',', decimal='.', header=None).astype(float))     
    time_05 = np.array(pd.read_csv("output/csv/2D_dam_05_break_time.csv", sep=',', decimal='.', header=None).astype(float))


    L = 0.5
    s= 0.01
    time_ad = time * np.sqrt(2 * 9.81 / L)
    pos_ad = (pos_x -3*s/2)/ L

    time_ad_05 = time_05 * np.sqrt(2 * 9.81 / L)
    pos_ad_05 = (pos_x_05 -3*s/2)/ L


    plt.plot(time_ad, pos_ad , color = "red", label= r"SPH simulation $\alpha = 0.1$", zorder = 1)
    plt.plot(time_ad_05, pos_ad_05 , color = "green", label= r"SPH simulation $\alpha = 0.5$", zorder = 1)
    
    plt.scatter(coord_x, coord_y+1, color = 'b', marker = "o", label = "Koshizuka and Oka 1996", zorder = 2)
    plt.xlim(0,3.3)
    plt.ylabel(r"x/L [-]")
    plt.xlabel(r"t$\sqrt{2\,g/L}$ [-]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/2D_dam_break_alpha.PDF')

    plt.show()



def MRUA():
    pos_z = np.array(pd.read_csv("output/csv/2D_splash_0_pos_z.csv", sep=',', decimal='.', header=None).astype(float)).T
    time =  np.array(pd.read_csv("output/csv/2D_splash_time.csv", decimal='.', header=None).astype(float)).T
    pos_z_RK2 = np.array(pd.read_csv("output/csv/2D_splash_RK_0_pos_z.csv", sep=',', decimal='.', header=None).astype(float)).T
    time_RK2 =  np.array(pd.read_csv("output/csv/2D_splash_RK_time.csv", decimal='.', header=None).astype(float)).T
    
    Z_th =  np.round(0.6 - 0.5*9.81*time[0]**2, decimals = 15)
    Z_th_RK2 =  np.round(0.6 - 0.5*9.81*time_RK2[0]**2, decimals = 15)
    print(Z_th_RK2-pos_z)
    
    plt.figure(figsize=(7, 5))
    plt.plot(time[0], Z_th, "b", label = "Analytical solution", zorder = 1)
    plt.scatter(time[0][::2], pos_z[0][::2] , color = "red", marker = "o", label = "SPH simulation Euler" , zorder = 2)
    plt.scatter(time_RK2[0][::2], pos_z_RK2[0][::2] , color = "green", marker = "+", label = "SPH simulation RK22" , zorder = 2)
    
    plt.ylim(0,0.8)
    plt.xlim(0,0.5)
    plt.ylabel(r"Postion $z(t)$ [m]")
    plt.xlabel(r"Time $t$ [s]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/2D_MRUA.PDF')
    plt.show()


    error = np.abs((pos_z-Z_th)/Z_th)
    error_RK2 = np.abs((pos_z_RK2-Z_th_RK2)/Z_th_RK2)
    plt.plot(time[0], error[0], color = "orange", label = "Relative error Euler")
    plt.plot(time_RK2[0], error_RK2[0], color = "green", label = "Relative error RK22")
    plt.legend(loc = "best")
    plt.ylabel(r"Relative error [-]")
    plt.xlabel(r"Time $t$ [s]")
    plt.yscale("log")
    plt.xlim(0,0.5)
    plt.grid(True, alpha = 0.5)
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/2D_error_MRUA.PDF')
    plt.show()


    pos_z = np.array(pd.read_csv("output/csv/3D_splash_1105_pos_z.csv", sep=',', decimal='.', header=None).astype(float)).T
    time =  np.array(pd.read_csv("output/csv/3D_splash_time.csv", decimal='.', header=None).astype(float)).T
    
    
    Z_th = pos_z[0][0] - 0.5*9.81*time[0]**2
    
    
    plt.figure(figsize=(7, 5))
    plt.plot(time[0], Z_th, "b", label = "Analytical solution", zorder = 1)
    plt.scatter(time[0][::2], pos_z[0][::2] , color = "red", marker = "o", label = "SPH simulation" , zorder = 2)
    
    plt.ylim(0,0.8)
    plt.xlim(0,0.5)
    plt.ylabel(r"Postion $z(t)$ [m]")
    plt.xlabel(r"Time $t$ [s]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/3D_MRUA.PDF')
    plt.show()


    error = np.abs((pos_z-Z_th)/Z_th)
    plt.plot(time[0], error[0], color = "orange", label = "Relative error")
    plt.legend(loc = "best")
    plt.ylabel(r"Relative error [-]")
    plt.xlabel(r"Time $t$ [s]")
    plt.yscale("log")
    plt.xlim(0,0.5)
    plt.grid(True, alpha = 0.5)
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/3D_error_MRUA.PDF')
    plt.show()

    
def linked_compa():
    nb_part = [9261,68921,132651,1030301]
    t_link = [0.007402,0.04134,0.074223,0.6910427]
    t_naive = [ 0.017371,0.941326,3.51582,244.755]   
    plt.figure(figsize=(8.5, 5))
    plt.scatter(nb_part,t_link, color ="red", marker="o", label= "Linked list algorithm")
    plt.scatter(nb_part, t_naive, color ="blue", marker="o", label= "Naive algorithm")
    plt.plot(nb_part,t_link, color ="red", linestyle ="--")
    plt.plot(nb_part, t_naive, color ="blue", linestyle ="--")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"Number of particles [-]")
    plt.ylabel(r"Time [s]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/compa_algo.PDF')
    plt.show()
    plt.show()

def surface_tension():
    '''p = np.array(pd.read_csv("output/csv/Euler_2D_surface_tension_p.csv", sep = ',', decimal='.', header=None))[-1]
    max_pos_x = np.array(pd.read_csv("output/csv/2D_surface_tension_max_pos_x.csv", sep = ',', decimal='.', header=None)).T[0]
    max_pos_z = np.array(pd.read_csv("output/csv/2D_surface_tension_max_pos_z.csv", sep = ',', decimal='.', header=None)).T[0]
    min_pos_x = np.array(pd.read_csv("output/csv/2D_surface_tension_min_pos_x.csv", sep = ',', decimal='.', header=None)).T[0]
    min_pos_z = np.array(pd.read_csv("output/csv/2D_surface_tension_min_pos_z.csv", sep = ',', decimal='.', header=None)).T[0]
    time =  np.array(pd.read_csv("output/csv/2D_surface_tension_time.csv", decimal='.', header=None).astype(float)).T[0]

    x = np.linspace(0,0.03,len(p))

    radius_x = max_pos_x - min_pos_x
    radius_z = max_pos_z - min_pos_z
    L=0.01
    R_th = L/((np.pi)**0.5)
    plt.plot(time, radius_x/2)
    plt.plot(time, radius_z/2)
    plt.axhline(y=R_th, color='r', linestyle='--')
    plt.show()
    
    plt.plot(x, p)
    plt.axhline(y=0.072/R_th, color='r', linestyle='--')
    plt.show()    '''

    p = np.array(pd.read_csv("output/csv/Euler_2D_surface_tension_p.csv", sep = ',', decimal='.', header=None))[-1]
    max_pos_x = np.array(pd.read_csv("output/csv/2D_surface_tension_max_pos_x.csv", sep = ',', decimal='.', header=None)).T[0]
    max_pos_z = np.array(pd.read_csv("output/csv/2D_surface_tension_max_pos_z.csv", sep = ',', decimal='.', header=None)).T[0]
    min_pos_x = np.array(pd.read_csv("output/csv/2D_surface_tension_min_pos_x.csv", sep = ',', decimal='.', header=None)).T[0]
    min_pos_z = np.array(pd.read_csv("output/csv/2D_surface_tension_min_pos_z.csv", sep = ',', decimal='.', header=None)).T[0]
    time =  np.array(pd.read_csv("output/csv/2D_surface_tension_time.csv", decimal='.', header=None).astype(float)).T[0]

    x = np.linspace(0,0.03,len(p))
    
    radius_x = max_pos_x - min_pos_x
    radius_z = max_pos_z - min_pos_z
    L = 0.01
    R_th = L/(np.sqrt(np.pi))
    rad = (radius_x + radius_z)/4
    print(rad[9])
    plt.figure(figsize=(7, 5))
    plt.plot(time,rad, label = "Average radius", color="red")

    plt.axhline(y=R_th, color='blue', linestyle='--', label="Theorical radius")

    plt.xlabel(r"Time $t$ [s]")
    plt.ylabel(r"Radius $R$ [m]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig("Pictures/radius_over_time.PDF")
    plt.show()

    plt.figure(figsize=(7, 5))
    plt.plot(x,p, label="Pressure ", color = "red")
    plt.axhline(y=0.072/R_th, color='b', linestyle='--', label="Theorical pressure ")
    plt.xlabel(r"Position $x$ [m]")
    plt.ylabel(r"Pressure $p$ [Pa]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig("Pictures/pressure_over_time.PDF")
    plt.show()
    
def main():
    
    current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
    isLatex(Lat)

    # Put in comment what you don't want to plot

    #PoiseuilleFlow()
    #Hydrostatic()
    #StateEquation()
    #dam_break()
    #MRUA()
    #linked_compa()
    surface_tension()

    

if __name__ == "__main__":
    main()

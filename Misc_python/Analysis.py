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
    z = np.linspace(0.04, 0.82, len(p))
    
    z_th = np.linspace(0,0.78,len(p))
    
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
        plt.savefig(f'{current_directory}/Pictures/2D_hydrostatic.PDF')
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
        plt.savefig(f'{current_directory}/Pictures/2D_hydrostatic_convergence.PDF')
    plt.show()

    

    plt.figure(figsize=(7, 5))
    plt.axhline(y=1000*9.81*(z_th[-part+1]-dist), color='blue', zorder = 1)
    plt.plot(time[0], p_part, color = 'red', label = 'Theoretical hydrostatic pressure', zorder = 2)
    plt.ylabel(r"Pressure $p(t)$ [Pa]")
    plt.xlabel(r"Time $t$ [s]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/2D_hydrostatic_convergence.PDF')
    plt.show()


    #2D hydrostatic RK22

    p_RK2 = np.array(pd.read_csv("output/csv/RK22_2D_RK2_hydrostatic_p.csv", sep = ',', decimal='.', header=None))[-1]
    time_RK2 =  np.array(pd.read_csv("output/csv/2D_RK2_hydrostatic_time.csv", decimal='.', header=None).astype(float)).T
    
    z_RK2 = np.linspace(0.04, 0.82, len(p_RK2))
    
    z_th_RK2 = np.linspace(0,0.78,len(p_RK2))
    
    p_th_RK2 = rho_ref*g*z_th_RK2
    dist = 0.116
    part = 10

    plt.figure(figsize=(7, 5))
    plt.scatter(z_RK2-dist, p_RK2[::-1], color = 'r', label = 'SPH pressure', marker="o")
    plt.plot(z_th_RK2, p_th_RK2, color = 'b', label = 'Theoretical hydrostatic pressure')
    plt.ylabel(r"Pressure $p(z)$ [Pa]")
    plt.xlabel(r"Depth $z$ [m]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/2D_RK22_hydrostatic.PDF')
    plt.show()

    p_part_RK2 = np.array(pd.read_csv("output/csv/RK22_2D_RK2_hydrostatic_p.csv", sep = ',', decimal='.', header=None)).T[part]
    

    plt.figure(figsize=(7, 5))
    plt.axhline(y=1000*9.81*(z_th_RK2[-part+1]-dist), color='blue', zorder = 1)
    plt.plot(time_RK2[0], p_part_RK2, color = 'red', label = 'Theoretical hydrostatic pressure', zorder = 2)
    plt.ylabel(r"Pressure $p(t)$ [Pa]")
    plt.xlabel(r"Time $t$ [s]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/2D_RK22_hydrostatic_convergence.PDF')
    plt.show()

    p = p[::-1]
    p = p[4:-2]
    z_th=z_th[4:-2]
    p_th = rho_ref*g*z_th
    p_th = p_th-p_th[0]
    z_th = z_th - z_th[0]
    error = np.abs((p[1:]-p_th[1:])/p_th[1:])

    p_RK2 = p_RK2[::-1]
    p_RK2 = p_RK2[4:-2]
    z_th_RK2=z_th_RK2[4:-2]
    p_th_RK2 = rho_ref*g*z_th_RK2
    p_th_RK2 = p_th_RK2 - p_th_RK2[0]
    z_th_RK2= z_th_RK2 - z_th_RK2[0]

    error_RK2 = np.abs((p_RK2[1:]-p_th_RK2[1:])/p_th_RK2[1:])
    plt.plot(z_RK2[4:-2]-dist, p_RK2 ,label ="SPH RK22", linestyle="--")
    plt.plot(z_th_RK2,p_th_RK2,label ="theo RK22")
    plt.plot(z[4:-2]-dist, p,label ="SPH Euler", linestyle="--")
    plt.plot(z_th,p_th,label ="theo Euler")
    plt.legend()
    plt.show()
    plt.yscale("log")
    plt.plot(z_th[1:],error, color = "red",label="Error Euler"  )
    plt.plot(z_th_RK2[1:],error_RK2, color = "blue",linestyle = "--",label="Error RK22"  )
    plt.ylabel(r"Relative error [-]")
    plt.xlabel(r"Depth $z$ [m]")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    if(save):
        current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
        plt.savefig(f'{current_directory}/Pictures/2D_Error_euler_RK22.PDF')
    plt.show()
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
    

    

    
    L = 0.2
    s= 0.05
    time_ad = time * np.sqrt(2 * 9.81 / L)
    pos_ad = (pos_x -3*s/2)/ L

    plt.plot(time_ad, pos_ad , color = "red", label= "SPH simulation", zorder = 1)
    plt.scatter(coord_x, coord_y+1, color = 'b', marker = "o", label = "Koshizuka and Oka 1996", zorder = 2)
    plt.xlim(0,3.5)
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
    s= 0.01
    time_ad = time * np.sqrt(2 * 9.81 / L)
    pos_ad = (pos_x -3*s/2)/ L

    plt.plot(time_ad, pos_ad , color = "red", label= "SPH simulation", zorder = 1)
    plt.scatter(coord_x, coord_y+1, color = 'b', marker = "o", label = "Koshizuka and Oka 1996", zorder = 2)
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
    pos_z = np.array(pd.read_csv("output/csv/2D_splash_78_pos_z.csv", sep=',', decimal='.', header=None).astype(float)).T
    time =  np.array(pd.read_csv("output/csv/2D_splash_time.csv", decimal='.', header=None).astype(float)).T
    
    
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
        plt.savefig(f'{current_directory}/Pictures/2D_MRUA.PDF')
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
    p = np.array(pd.read_csv("output/csv/Euler_2D_surface_tension_p.csv", sep = ',', decimal='.', header=None))[-1]
    max_pos_x = np.array(pd.read_csv("output/csv/2D_surface_tension_max_pos_x.csv", sep = ',', decimal='.', header=None)).T[0]
    max_pos_z = np.array(pd.read_csv("output/csv/2D_surface_tension_max_pos_z.csv", sep = ',', decimal='.', header=None)).T[0]
    min_pos_x = np.array(pd.read_csv("output/csv/2D_surface_tension_min_pos_x.csv", sep = ',', decimal='.', header=None)).T[0]
    min_pos_z = np.array(pd.read_csv("output/csv/2D_surface_tension_min_pos_z.csv", sep = ',', decimal='.', header=None)).T[0]
    time =  np.array(pd.read_csv("output/csv/2D_surface_tension_time.csv", decimal='.', header=None).astype(float)).T[0]

    x = np.linspace(0,0.03,len(p))

    print(p)
    radius_x = max_pos_x - min_pos_x
    radius_z = max_pos_z - min_pos_z

    
    plt.plot(time, radius_x/2)
    plt.plot(time, radius_z/2)
    plt.axhline(y=0.005641895835, color='r', linestyle='--')
    plt.show()
    
    plt.plot(x, p)
    plt.show()    


    max_pos_x = np.array(pd.read_csv("output/csv/3D_surface_tension_max_pos_x.csv", sep = ',', decimal='.', header=None)).T[0]
    max_pos_y = np.array(pd.read_csv("output/csv/3D_surface_tension_max_pos_y.csv", sep = ',', decimal='.', header=None)).T[0]
    max_pos_z = np.array(pd.read_csv("output/csv/3D_surface_tension_max_pos_z.csv", sep = ',', decimal='.', header=None)).T[0]
    min_pos_x = np.array(pd.read_csv("output/csv/3D_surface_tension_min_pos_x.csv", sep = ',', decimal='.', header=None)).T[0]
    min_pos_y = np.array(pd.read_csv("output/csv/3D_surface_tension_min_pos_y.csv", sep = ',', decimal='.', header=None)).T[0]
    min_pos_z = np.array(pd.read_csv("output/csv/3D_surface_tension_min_pos_z.csv", sep = ',', decimal='.', header=None)).T[0]
    time =  np.array(pd.read_csv("output/csv/3D_surface_tension_time.csv", decimal='.', header=None).astype(float)).T[0]

    x = np.linspace(0,0.03,len(p))

    print(p)
    radius_x = max_pos_x - min_pos_x
    radius_x = max_pos_y - min_pos_y
    radius_z = max_pos_z - min_pos_z

    
    plt.plot(time, radius_x/2)
    plt.plot(time, radius_z/2)
    plt.axhline(y=0.005641895835, color='r', linestyle='--')
    plt.show()

def main():
    
    current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
    isLatex(Lat)

    # Put in comment what you don't want to plot

    #PoiseuilleFlow()
    #Hydrostatic()
    #StateEquation()
    dam_break()
    #MRUA()
    #linked_compa()
    #surface_tension()

    

if __name__ == "__main__":
    main()

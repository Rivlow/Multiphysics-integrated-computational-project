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

s =  np.array([0.1, 0.05, 0.025, 0.01, 0.005])
nb_part = np.array([318, 1110, 4134, 24726, 97446])
OMP_time =  np.array([17.2, 31 , 76.5, 591,2239])
serial_time =  np.array([10, 40.3, 173.4, 1207, 4968])


def plot_tab(PLOT_time, PLOT_scale):
    if PLOT_time:

        plt.xscale('log')
        plt.semilogx(nb_part, OMP_time, color = "r")
        plt.scatter(nb_part, OMP_time, color = "r", label = "OMP")

        plt.scatter(nb_part, serial_time, color = "b", label = "Serial")
        plt.semilogx(nb_part, serial_time, color = "b")

        plt.xlabel("Number of particles [-]")
        plt.ylabel("Total computation time [s]")

        plt.grid(True)
        plt.legend(loc = "best")
        plt.tight_layout()
        plt.savefig(f"{current_directory}/Pictures/Speed_Up_OMP.PDF")
        plt.show()
        
    if PLOT_scale:
        plt.xscale('log')
        plt.scatter(nb_part, scaling, color = "r")
        plt.semilogx(nb_part, scaling, color = "r")
        plt.grid(True)
        plt.xlabel("Number of particles [-]")
        plt.ylabel("Speed-up [-]")
        plt.legend(loc = "best")
        plt.tight_layout()
        plt.savefig(f"{current_directory}/Pictures/Scaling_OMP.PDF")
        plt.show()
    
    
    
PLOT_time = False
PLOT_scale = True
current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))

scaling = serial_time/OMP_time
plot_tab(PLOT_time, PLOT_scale)




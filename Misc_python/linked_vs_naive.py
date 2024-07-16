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

def plotComparison(nb_part, linked, naive):
    
    plt.xscale("log")
    plt.yscale("log")

    plt.scatter(nb_part, linked, label="Linked-list", color = "red")
    plt.plot(nb_part, linked, color = "red", linestyle = "--")
    
    plt.scatter(nb_part, naive, label = "Naive", color = "blue")
    plt.plot(nb_part, naive, color = "blue", linestyle = "--")
    
    plt.ylabel('Average time of research per step [s]')
    plt.xlabel("Number of particles [-]")
    plt.legend(loc = "best")
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig(f"{current_directory}/Pictures/linked_naice_comparison.PDF")
    plt.show()

step = 100
nb_part = [24, 679, 2503, 14887, 23107, 30278, 58567]
linked = [0.00019067, 0.00028684, 0.00045728, 0.00232608, 0.003866, 0.00558104, 0.01220402]
naive = [0.00020428, 0.00052669, 0.00349627, 0.09740487, 0.26255906, 0.47741576, 2.50199889]

current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))


isLatex(True)
plotComparison(nb_part, linked, naive)

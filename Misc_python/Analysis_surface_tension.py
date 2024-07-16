import vtk
#from vtk.util.numpy_support import vtk_to_numpy
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys
import pandas as pd

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




def plotgraph(file, ite, s):
    pressure = pd.read_csv(file, sep = ',', decimal='.', header=None)

    x = np.arange(0.0, 3,s)

    plt.figure(figsize=(8.5,5))
    pressure = pressure.T
    print(len(x))
    print(len(pressure))
    plt.scatter(x,pressure[ite])
    plt.ylabel(r"Pressure [Pa]")
    plt.xlabel(r"Position [m]")
    plt.grid(alpha = 0.5)
    pressuremean = pressure[ite][21:40].mean()
    print(pressuremean)
    #plt.savefig(f"output/video/p_cube_to_sphere_{s}.PDF")
    plt.show(block = False)
    plt.close()
    return x, pressure 

file1 = "output/cube_to_sphere01/p.csv"
ite1 = 499
file2 ="output/cube_to_sphere0125/p.csv"
ite2 = 499
file3 = "output/cube_to_sphere005/p.csv"
ite3 = 799
x2 , p2 = plotgraph(file2,ite2,0.125)
x1 , p1 = plotgraph(file1,ite1,0.1)
x3 , p3 = plotgraph(file3,ite3,0.05)
plt.figure(figsize=(8.5,5))
plt.scatter(x2,p2[ite2], label = " s = 0.125")
plt.scatter(x1,p1[ite1], label = " s = 0.1")

plt.scatter(x3,p3[ite3], label = " s = 0.05")
plt.grid(alpha= 0.50)
plt.ylabel(r"Pressure [Pa]")
plt.xlabel(r"Position [m]")

plt.legend()
#☻plt.savefig(f"output/video/p_many_spacing.PDF")
plt.show()
plt.figure(figsize=(8.5,5))
file4 = "output/cube_to_sphere01_alpha50/p.csv"
ite4 = 499
x4, p4 = plotgraph(file4,ite4,0.1)
plt.scatter(x1,p1[ite1], label = r"$\alpha_{curv}$ = 10")
plt.scatter(x4,p4[ite4], label = r"$\alpha_{curv}$ = 50")
plt.grid(alpha= 0.50)
plt.ylabel(r"Pressure [Pa]")
plt.xlabel(r"Position [m]")
plt.legend()
plt.savefig(f"output/video/p_many_alpha_s01.PDF")
plt.show()

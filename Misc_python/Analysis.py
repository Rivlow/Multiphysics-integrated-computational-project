import vtk
#from vtk.util.numpy_support import vtk_to_numpy
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import sys
import pandas as pd


file = "output/p.csv"
ite = 100
pressure = pd.read_csv(file, sep = ',', decimal='.', header=None)
print(len(pressure))

x = np.arange(0.05, 0.7, 0.05)

x1 = np.arange(0,0.65,0.05)
pressure = pressure.T
print(x)
plt.plot(x,pressure[ite])
plt.plot(x,1000*9.81*x1)

plt.show()



'''
def read_vtp(path):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(path)
    reader.Update()
    data = reader.GetOutput().GetPointData()
    field_count = data.GetNumberOfArrays()
    return {data.GetArrayName(i): vtk_to_numpy(data.GetArray(i)) for i in range(field_count)}

def read_vtp_files_in_folder(folder_path):
    vtp_files = {}
    for filename in os.listdir(folder_path):
        if filename.endswith(".vtp"):
            file_path = os.path.join(folder_path, filename)
            file_data = read_vtp(file_path)
            for array_name, array_data in file_data.items():
                if array_name not in vtp_files:
                    vtp_files[array_name] = []
                vtp_files[array_name].append(array_data)
    return vtp_files

folder_path = "output/cube_bcp_particles_less"
vtp_files_data = read_vtp_files_in_folder(folder_path)


rho = vtp_files_data["rho_array"]
print((len(vtp_files_data)))


val_rho = []

for i, array_data in enumerate(rho):
    val_rho.append(array_data[2])
    print(f"Data from file {i + 1}: {array_data}")
    
plt.plot(val_rho)
plt.show()


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
    
def getData(directory):
    
    df, name = [], []
    
    # Search all csv. files
    if os.path.exists(directory) and os.path.isdir(directory):
        for filename in os.listdir(directory):
            if filename.endswith('.csv'):
                name.append(filename)
                csv_file = os.path.join(directory, filename)  
                if csv_file is not None:
                    df.append(np.loadtxt(csv_file, delimiter=",", dtype=float))
        
        return df, name
    else:
        print(f"No csv files in folder : {directory}.")

def plotSingleVariable(analysis_type, data, particle, iteration, name):
    
    if (analysis_type["Time"]):
        
        plt.plot(range(data.shape[0]), data[:, particle], label=f'test') 
        plt.legend(loc='best')
        plt.xlabel('Iterations')
        plt.ylabel('Values')
        plt.grid(True)
        plt.show()
        
    elif (analysis_type["Spacial"]):
        
        x = np.linspace(0.075, 1.15, len(data[iteration, :]))
        plt.plot(x, data[iteration, :])
        #plt.scatter(x, data[iteration, :],)
        plt.xlabel('Height [m]')
        plt.ylabel(f'{name}')
        plt.grid(True)
        plt.show()
        
    else:
        print("Error, choose one type of plot.")
        
def plotStateEquation(rho, p, iteration):
    
    p_val = p[iteration, :]
    rho_val = rho[iteration, :]
    x = np.linspace(0.01, 1.2, len(p_val))
    
    c_0 = 30
    rho_0 = 1000
    gamma = 7  # Exemple de valeur pour gamma
    B = (c_0**2)*rho_0/gamma
    T = 273.15+25
    M = 18e-3
    R = 8.314
  
    #plt.scatter(x[11:], p_val[11:])
    #plt.scatter(x[11:], rho_val[11:])
    p_compare = np.linspace(48.7799, 12397, 100)
    rho_compare = rho_0 * ((p_compare + 1) / B) ** (1 / gamma)
    plt.xscale("log")
    plt.plot(p_val[11:], rho_val[11:], color = 'b', label = 'SPH')
    #plt.plot(p_compare, rho_compare, color = 'r', label = 'quasi incompressible (theorical)')
    plt.xlabel("Pressure [Pa]")
    plt.ylabel(r"Density [kg/$m^3$]")
    plt.grid(True)
    plt.legend(loc = "best")
    plt.show()
    
    
    
def plotHydrostaticPressure(p, current_directory):

    rho_ref = 1023 # average density for ghost particles (using "plotSingleVariable")
    g = 9.81
    z = np.linspace(0.25, 0.9, len(p))
    p_th = rho_ref*g*z
    
    print(p_th[0])
    print(p[-1])
    
    print(((p[-1]-p[0])/(z[-1]-z[0]))/9.81)

    plt.figure()
    plt.scatter(z, p[::-1], color = 'r', label = 'SPH pressure')
    plt.plot(z, p_th, color = 'b', label = 'Theoretical hydrostatic pressure')

    plt.ylabel("Pressure [Pa]")
    plt.xlabel(r"Height [m]")
    plt.grid(True)
    plt.legend(loc = "best")
    plt.tight_layout()
    plt.savefig(f'{current_directory}/Pictures/hydrostatic_pressure.PDF')
    plt.show()



        

def main():
    
    current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
    outputFile = os.path.dirname(current_directory) + "\\output"  
    all_data = getData(outputFile)
    
    print(f"Data available : {all_data[1]}")

    analysis_type = {"Time":False, "Spacial":True} # chose only one as "True"
    particle = 0 # if Time plot desired, have to chose which particle to look at
    iteration = 999-1 # if Spacial plot desired, have to chose which iteration to look at
    
    rho = all_data[0][1]
    p = all_data[0][0]

    
    name = "Pressure"
    #name = "Rho"
    
    isLatex(True)
    #plotSingleVariable(analysis_type, rho, particle, iteration, name)
    #plotStateEquation(rho, p, iteration)
    plotHydrostaticPressure(p[iteration,:], current_directory)
    
    nstepT = 200000
    dt = 2.5e-5
    
    print(nstepT*dt)
    
    

 

    
if __name__ == "__main__":
    main()


'''


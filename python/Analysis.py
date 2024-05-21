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
ite = 800
pressure = pd.read_csv(file, sep = ',', decimal='.', header=None)
print(len(pressure))

x = np.arange(0.0, 0.03, 0.00125)


pressure = pressure.T
print(len(pressure))
plt.plot(x,pressure[ite])

plt.show()
maxp = max(pressure[ite])
print("max", maxp)
minp = min(pressure[ite])
print("min", minp)
deltaP = maxp - minp
print(deltaP)

"""
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

    
path = os.getcwd()


def plotData(all_data):
    
    data_tot = all_data[0]
    names = all_data[1]
    print(f"Data available : {names}")

    for data, name in zip(data_tot, names):
        
        t = np.arange(0, np.shape(data)[0], 1)
        var = np.zeros(len(t))
        
        print(data)

        for i, rows in enumerate(data):
            var[i] = rows[0]
        
        plt.figure()
        plt.scatter(t, var)
        plt.xlabel("time [s]")
        plt.ylabel(f"{name[:-4]}(t)")
        plt.title(f"{name[:-4]}")
        plt.tight_layout()
        plt.show()

        
    
    
def main():
    
    # Chemin absolu du répertoire du script actuel
    current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
    outputFile = os.path.dirname(current_directory) + "\\output"  
    all_data = getData(outputFile)
    
    analysis_type = {"Time":False, "Spacial":True} # chose only one as "True"
    
    print(all_data[1])
    data = all_data[0][1]
    # Création des subplots
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    # Plot 2D
    for i in range(data.shape[1]):
        axs[0].plot(range(data.shape[0]), data[:, i], label=f'part {i}')
    axs[0].legend(loc='best')
    axs[0].set_title('Plot 2D')
    axs[0].set_xlabel('Iterations')
    axs[0].set_ylabel('Values')

    # Plot 3D
    ax = fig.add_subplot(122, projection='3d')
    for i in range(data.shape[1]):
        ax.plot(range(data.shape[0]), data[:, i], zs=i)
    ax.set_title('Plot 3D')
    ax.set_xlabel('Iterations')
    ax.set_ylabel('Y')
    ax.set_zlabel('Ghost particles')

    plt.tight_layout()
    plt.show()
 
    
    
    
if __name__ == "__main__":
    main()
"""




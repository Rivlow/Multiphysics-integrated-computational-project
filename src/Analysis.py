import vtk
from vtk.util.numpy_support import vtk_to_numpy
import os
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

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
    """
    
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
        print(f"No csv. files in folder : {directory}.")

    
path = os.getcwd()


def plotData(all_data, particle):
    
    data_tot = all_data[0]
    names = all_data[1]
    print(f"Data available : {names}")

    for data, name in zip(data_tot, names):
        
        t = np.arange(0, np.shape(data)[0]/39, 1)
        var = np.zeros(len(t))

        for i, rows in enumerate(data/39):
            var[i] = rows[particle+i*39]
        
        plt.figure()
        plt.scatter(t, var)
        plt.xlabel("time [s]")
        plt.ylabel(f"{name[:-4]}(t)")
        plt.title(f"{name[:-4]}")
        plt.tight_layout()
        plt.show()
        
    
    
def main():
    
    outputFile = "output/"
    particle = 20
    
    all_data = getData(outputFile)
    plotData(all_data, particle)
    
    
    
if __name__ == "__main__":
    main()





import vtk
#from vtk.util.numpy_support import vtk_to_numpy
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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

        
    
    
def main(iteration):
    
    latex = False
    isLatex(latex)
    
    # Chemin absolu du répertoire du script actuel
    current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
    outputFile = os.path.dirname(current_directory) + "\\output"  
    all_data = getData(outputFile)
    
    analysis_type = {"Time":False, "Spacial":True} # chose only one as "True"
    
    print(all_data[1])
    name = all_data[1]
    data = all_data[0][0]
    p_val, rho_val, p_val_init, p_val_end = [], [], [], []
    
    #fig, axs = plt.subplots(1, 2, figsize=(12, 6))


    if analysis_type["Time"]:
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

    elif analysis_type["Spacial"]:
        '''
        print(len(all_data[0][0][iteration,:]))
        for i in range(data.shape[1]):
            p_val.append(all_data[0][0][iteration,i])
            rho_val.append(all_data[0][1][iteration,i])
            
        x_val = np.arange(0, len(p_val), 1)
        axs[0].set_title(f"At iteration {iteration}")
        axs[0].plot(x_val, p_val)
        axs[0].scatter(x_val, p_val)
        axs[0].legend(loc='best')
        axs[0].set_xlabel('Ghost particles')
        axs[0].set_ylabel('Pressure p [Pa]')
        
        axs[1].set_title(f"At iteration {iteration}")
        axs[1].plot(x_val, rho_val)
        axs[1].scatter(x_val, rho_val)
        axs[1].legend(loc='best')
        axs[1].set_xlabel('Ghost particles')
        axs[1].set_ylabel('Density ')
        '''
        
        rho = 1000
        g = 9.81
        
        for i in range(data.shape[1]):
            p_val_init.append(all_data[0][0][0,i])
            p_val_end.append(all_data[0][0][iteration,i])
            
        x_val = np.arange(0, len(p_val_end), 1)
        z = np.linspace(0, 1, len(x_val))
        p_hydro = rho*g*z
        plt.plot(x_val, p_hydro, label = "Theoretical pressure")
        plt.scatter(x_val, p_val_end, label = "SPH pressure", color = "red")
        plt.legend(loc='best')
        plt.xlabel('Ghost particles')
        plt.ylabel('Pressure p [Pa]')
        plt.grid(True)
        state_equation = "Quasi_incompressible_fluid"
        plt.savefig(f"{current_directory}\Pictures\hydrostatic_{state_equation}.PDF")

        
        
    else:
        print("No plot type chosen")

    plt.tight_layout()
    plt.show()
    
  
    
    
if __name__ == "__main__":
    iteration = 50
    main(iteration-1)





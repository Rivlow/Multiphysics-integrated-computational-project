import vtk
from vtk.util.numpy_support import vtk_to_numpy
import os
from matplotlib import pyplot as plt

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
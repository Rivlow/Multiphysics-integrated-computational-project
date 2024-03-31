import vtk

"""
# Charger le fichier .vtp
reader = vtk.vtkXMLPolyDataReader()
reader.SetFileName("sph.vtp")
reader.Update()

# Créer un mapper pour chaque donnée
mappers = []
for i in range(reader.GetNumberOfPointArrays()):
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(reader.GetOutput())
    mapper.SetScalarModeToUsePointData()
    mapper.SelectColorArray(i)
    mappers.append(mapper)

# Créer un acteur pour chaque mapper
actors = []
for mapper in mappers:
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actors.append(actor)

# Créer une scène
renderer = vtk.vtkRenderer()
for actor in actors:
    renderer.AddActor(actor)
renderer.SetBackground(1, 1, 1)  # Couleur de fond blanc

# Créer une fenêtre de rendu
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)

# Créer un interpréteur d'événements
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

# Démarrer le rendu
renderWindow.Render()
renderWindowInteractor.Start()
"""

import vtk
from vtk.util.numpy_support import vtk_to_numpy

def read_vtp(path):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(path)
    reader.Update()
    data = reader.GetOutput().GetPointData()
    field_count = data.GetNumberOfArrays()
    return {data.GetArrayName(i): vtk_to_numpy(data.GetArray(i)) for i in range(field_count)}

read_vtp("VLP-16_Dual_1.vtp")

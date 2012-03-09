import sys

import vtk
import ca_smoothing

stl = vtk.vtkSTLReader()
stl.SetFileName(sys.argv[1])
stl.Update()

normals = vtk.vtkPolyDataNormals()
normals.SetInput(stl.GetOutput())
normals.ComputeCellNormalsOn()
normals.Update()

clean = vtk.vtkCleanPolyData()
clean.SetInput(normals.GetOutput())
clean.Update()

pd = clean.GetOutput()
pd.BuildLinks()

tpd = ca_smoothing.ca_smoothing(pd, 0.7, 3, 0.2, 10)

ply = vtk.vtkPLYWriter()
ply.SetFileName(sys.argv[2])
ply.SetInput(tpd)
ply.Write()

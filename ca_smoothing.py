#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import argparse
import itertools
import numpy as np
import sys
import vtk

def find_staircase_artifacts(pd, stack_orientation=(0, 0, 1), T=0.7):
    """
    This function is used to find vertices at staircase artifacts, which are
    those vertices whose incident faces' orientation differences are
    aproximated 90°. The argument `pd' is a vtkPolydata, it have to be its
    faces orientation (normal) calculated before and make sure that
    BuildLinks() has been called.
    """
    output = []
    nv = pd.GetNumberOfPoints() # Number of vertices.
    scalars = vtk.vtkFloatArray()
    for vid in xrange(nv):
        idfaces = vtk.vtkIdList()
        pd.GetPointCells(vid, idfaces) # Getting faces connected to face vid.

        nf = idfaces.GetNumberOfIds()
        orientations = []
        # For each combination of connected faces
        for nid in xrange(nf):
            fid = idfaces.GetId(nid)
            ni = pd.GetCellData().GetArray("Normals").GetTuple(fid)

            of = 1 - np.dot(ni, stack_orientation)
            orientations.append(of)

        minimun = min(orientations)
        maximun = max(orientations)

        # Getting the ones which normals dot is 90°, its vertex is added to
        # output
        if maximun - minimun >= T:
            output.append(vid)
            scalars.InsertNextValue(1)
            print "aqui"
        else:
            scalars.InsertNextValue(0)


    pd.GetPointData().SetScalars(scalars)

    return output

def calc_artifacts_weight(pd, vertices_staircase, tmax=2.0):
    """
    Calculate the artifact weight based on distance of each vertex to its
    nearest staircase artifact vertex.
    pd - vtkPolydata;
    vertices_staircase - the identified staircase artifact vertices;
    tmax=2 - max distance the vertex must be to its nearest artifact vertex to
             considered to calculate the weight.
    """
    scalars = pd.GetPointData().GetScalars()
    for vid in vertices_staircase:
        vi = pd.GetPoint(vid)
        idfaces = vtk.vtkIdList()
        pd.GetPointCells(vid, idfaces) # Getting faces connected to face vid.

        nf = idfaces.GetNumberOfIds()
        # For each combination of connected faces
        for nid in xrange(nf):
            fid = idfaces.GetId(nid)
            face = pd.GetCell(fid)
            
            for i in xrange(3):
                vjid = face.GetPointId(i)
                #if vjid == vid:
                    #continue
                vj = pd.GetPoint(vjid)
                d = np.sqrt((vi[0]-vj[0])**2 + (vi[1]-vj[1])**2 + (vi[2]-vj[2])**2)
                if d <= tmax:
                    value = (1.0 - d/tmax)
                    scalars.SetValue(vjid, value)



def read_stl(filename):
    "It only reads a STL file, it can be binary or ascii format."

    stl_reader = vtk.vtkSTLReader()
    stl_reader.SetFileName(filename)
    stl_reader.Update()

    output = stl_reader.GetOutput()
    return output


def visualize(pd, vertices_staircase):
    lt = vtk.vtkLookupTable()
    lt.SetNumberOfColors(100)
    lt.SetValueRange(0, 1)
    #lt.Build()
    lt.SetTableValue(0, 0, 0, 1, 1)
    lt.SetTableValue(99, 1, 0, 0, 1)

    m = vtk.vtkPolyDataMapper()
    m.SetInput(pd)
    m.SetScalarRange(0.0, 1.0)
    m.SetLookupTable(lt)

    a = vtk.vtkActor()
    a.SetMapper(m)

    ren = vtk.vtkRenderer()
    ren.AddActor(a)

    renwin = vtk.vtkRenderWindow()
    renwin.AddRenderer(ren)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renwin)
    iren.Start()


def main():
    parser = argparse.ArgumentParser(prog='Context-aware mesh smoothing')
    parser.add_argument("file", help="A STL file", metavar="STL file")
    parser.add_argument("-t", "--threshold", help="Threshold to find staircase artifacts", 
                       type=float, default=0.7, dest='threshold')
    args = parser.parse_args()

    stl = read_stl(args.file)

    normals = vtk.vtkPolyDataNormals()
    normals.SetInput(stl)
    normals.ComputeCellNormalsOn()
    normals.Update()

    clean = vtk.vtkCleanPolyData()
    clean.SetInput(normals.GetOutput())
    clean.Update()
    
    pd = clean.GetOutput()
    pd.BuildLinks()
    vertices_staircase = find_staircase_artifacts(pd, T=args.threshold)
    calc_artifacts_weight(pd, vertices_staircase)
    visualize(pd, vertices_staircase)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInput(pd)
    writer.SetFileName("/tmp/test_cor.vtp")
    writer.Write()


if __name__ == '__main__':
    main()

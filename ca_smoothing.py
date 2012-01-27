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
        else:
            scalars.InsertNextValue(0)


    pd.GetPointData().SetScalars(scalars)

    return output

def get_near_vertices_to_v(pd, v, dmax):
    """
    Returns all vertices with distance at most "d" to the vertice "v" with
    their distances.
    pd - vtkPolydata
    v - the reference vertice
    dmax - the maximun distance.
    """
    visit = []
    vi = pd.GetPoint(v) # The position of vertice v
    n = 0
    while 1:
        idfaces = vtk.vtkIdList()
        pd.GetPointCells(v, idfaces) # Getting faces connected to face v
        nf = idfaces.GetNumberOfIds()
        # For each combination of connected faces
        for nid in xrange(nf):
            fid = idfaces.GetId(nid)
            face = pd.GetCell(fid)
            
            for i in xrange(3):
                vjid = face.GetPointId(i)
                if vjid not in visit:
                    vj = pd.GetPoint(vjid)
                    d = np.sqrt((vi[0]-vj[0])**2 + (vi[1]-vj[1])**2 + (vi[2]-vj[2])**2)
                    if d <= dmax:
                        yield vjid, d
                        visit.append(vjid)
        n += 1
        print n

        if n >= len(visit):
            break

        v = visit[n]

def calc_artifacts_weight(pd, vertices_staircase, tmax=5.0, bmin=0.1):
    """
    Calculate the artifact weight based on distance of each vertex to its
    nearest staircase artifact vertex.
    pd - vtkPolydata;
    vertices_staircase - the identified staircase artifact vertices;
    tmax=2 - max distance the vertex must be to its nearest artifact vertex to
             considered to calculate the weight;
    bmin=0.1 - The minimun weight.
    """
    weights = np.zeros(pd.GetNumberOfPoints(), dtype='float64')
    scalars = pd.GetPointData().GetScalars()
    for vi in vertices_staircase:
        for vj, d in get_near_vertices_to_v(pd, vi, tmax):
            value = (1.0 - d/tmax)
            if value > weights[vj]:
                weights[vj] = value
                scalars.SetValue(vj, value)
    weights = weights * (1 - bmin) + bmin
    return weights

def calculate_d(mesh, poly, pid):
    t = 0
    n = 0.0
    cell_ids = vtk.vtkIdList()
    p0 = np.array(poly.GetPoint(pid))
    mesh.GetPointCells(pid, cell_ids)
    for i in xrange(cell_ids.GetNumberOfIds()):
        point_ids = vtk.vtkIdList()
        mesh.GetCellPoints(cell_ids.GetId(i), point_ids)
        n += 1
        if point_ids.GetId(0) != pid:
            p1 = np.array(poly.GetPoint(point_ids.GetId(0)))
        else:
            p1 = np.array(poly.GetPoint(point_ids.GetId(1)))

        t = t + (p1 - p0)

    return t / n

def taubin_smooth(pd, weights, l, m, steps):
    edgesfilter = vtk.vtkExtractEdges()
    edgesfilter.SetInput(pd)
    edgesfilter.Update()

    edges = edgesfilter.GetOutput()

    new_pd = vtk.vtkPolyData()
    new_pd.DeepCopy(pd)

    points = new_pd.GetPoints()
    for s in xrange(steps):
        D = {}
        for i in xrange(edges.GetNumberOfPoints()):
            D[i] = calculate_d(edges, new_pd, i)
        for i in xrange(pd.GetNumberOfPoints()):
            p = np.array(points.GetPoint(i))
            pl = p + weights[i]*l*D[i]
            nx, ny, nz = pl
            points.SetPoint(i, nx, ny, nz)

        #D = {}
        #for i in xrange(edges.GetNumberOfPoints()):
            #D[i] = calculate_d(edges, new_pd, i)
        #for i in xrange(pd.GetNumberOfPoints()):
            #p = np.array(points.GetPoint(i))
            #pl = p + weights[i]*m*D[i]
            #nx, ny, nz = pl
            #points.SetPoint(i, nx, ny, nz)

        #D = {}
        #for i in xrange(edges.GetNumberOfPoints()):
            #D[i] = calculate_d(edges, i)
        #for i in xrange(pd.GetNumberOfPoints()):
            #x, y, z = points.GetPoint(i)
            #nx = x + m*D[i][0]
            #ny = y + m*D[i][1]
            #nz = z + m*D[i][2]
            #points.SetPoint(i, nx, ny, nz)

    #new_pd.SetPoints(points)
    return new_pd


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
    weights = calc_artifacts_weight(pd, vertices_staircase)
    print weights
    new_pd = taubin_smooth(pd, weights, 0.5, -0.53, 10)
    visualize(new_pd, vertices_staircase)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInput(new_pd)
    writer.SetFileName("/tmp/test_cor.vtp")
    writer.Write()


if __name__ == '__main__':
    main()

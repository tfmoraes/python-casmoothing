#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import vtk

def find_staircase_artifacts(pd):
    """
    This function is used to find vertices at staircase artifacts, which are
    those vertices whose incident faces' orientation differences are
    aproximated 90°. The argument `pd' is a vtkPolydata, it have to be its
    faces orientation (normal) calculated before.
    """
    pass

def read_stl(filename):
    "It only reads a STL file, it can be binary or ascii format."

    stl_reader = vtk.vtkSTLReader()
    stl_reader.SetFileName(filename)
    stl_reader.Update()

    output = stl_reader.GetOutput()
    return output

def main():
    pd = read_stl(sys.argv[1])


if __name__ == '__main__':
    main()

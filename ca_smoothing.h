//--------------------------------------------------------------------------
// Software:     Context Aware Smoothing
// Copyright:    (C) 2012  Centro de Pesquisas Renato Archer
// Homepage:     https://github.com/tfmoraes/context_aware_smoothing
// Contact:      tfmoraes@cti.gov.br
// License:      GNU - GPL 2 (LICENSE.txt/LICENCA.txt)
//--------------------------------------------------------------------------
//    Este programa e software livre; voce pode redistribui-lo e/ou
//    modifica-lo sob os termos da Licenca Publica Geral GNU, conforme
//    publicada pela Free Software Foundation; de acordo com a versao 2
//    da Licenca.
//
//    Este programa eh distribuido na expectativa de ser util, mas SEM
//    QUALQUER GARANTIA; sem mesmo a garantia implicita de
//    COMERCIALIZACAO ou de ADEQUACAO A QUALQUER PROPOSITO EM
//    PARTICULAR. Consulte a Licenca Publica Geral GNU para obter mais
//    detalhes.
//--------------------------------------------------------------------------

#include <math.h>

#ifdef HAVE_CXX0X
    #include<unordered_map>
    #define MAP std::unordered_map
#elif HAVE_TR1
    #include <tr1/unordered_map>
    #define MAP std::tr1::unordered_map
#else
    #include <map>
    #define MAP std::map
#endif

#include <queue>
#include <set>
#include <vector>
#include <stdlib.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSet.h>
#include <vtkSphereSource.h>
#include <vtkTriangleFilter.h>
#include <vtkExtractEdges.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkIdTypeArray.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkExtractSelection.h>
#include <vtkProperty.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkSTLReader.h>
#include <vtkPLYWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkCleanPolyData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkCell.h>

typedef struct _Point
{
  double x;
  double y;
  double z;
} Point;


vtkPolyData* ca_smoothing (vtkPolyData*, double, double, double, int);
vtkIdList* find_staircase_artifacts(vtkPolyData*, const double[3], double);
vtkIdList* get_near_vertices_to_v(vtkPolyData*, int, double);
vtkDoubleArray* calc_artifacts_weight(vtkPolyData*, vtkIdList*, double, double);
Point calc_d(vtkPolyData*, int);
vtkPolyData* taubin_smooth(vtkPolyData*, vtkDoubleArray*, double, double, int);
MAP<int, std::vector<int> > get_flat_areas(vtkPolyData*, double[3], double);

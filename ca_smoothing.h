#include <math.h>
#include <map>
#include<unordered_map>
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

vtkPolyData* read_stl(char*);
vtkIdList* find_staircase_artifacts(vtkPolyData*, const double[3], double);
vtkIdList* get_near_vertices_to_v(vtkPolyData*, int, double);
vtkDoubleArray* calc_artifacts_weight(vtkPolyData*, vtkIdList*, double, double);
Point calc_d(vtkPolyData*, int);
vtkPolyData* taubin_smooth(vtkPolyData*, vtkDoubleArray*, double, double, int);

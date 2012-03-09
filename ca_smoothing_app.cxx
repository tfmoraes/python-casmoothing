#include "ca_smoothing.h"

int main(int argc, char *argv[])
{
    const double stack_orientation[3] = { 0, 0, 1 };
    vtkPolyData *stl, *nm, *cl, *pd, *tpd;
    vtkPolyDataNormals *normals;
    vtkCleanPolyData *clean;
    vtkIdList *vertices_staircase;
    vtkDoubleArray* weights;
    
    printf("Reading STL\n");
    stl = read_stl(argv[1]);

    printf("Generating the normals\n");
    normals = vtkPolyDataNormals::New();
    normals->SetInput(stl);
    normals->ComputeCellNormalsOn();
    normals->Update();

    printf("Cleaning the polydata\n");
    clean = vtkCleanPolyData::New();
    clean->SetInput(normals->GetOutput());
    clean->Update();

    pd = clean->GetOutput();
    pd->BuildLinks();

    printf("Finding staircase artifacts\n");
    vertices_staircase = find_staircase_artifacts(pd, stack_orientation, atof(argv[2]));
    printf("Calculating the Weights\n");
    weights = calc_artifacts_weight(pd, vertices_staircase, atof(argv[3]), atof(argv[4]));
    printf("Taubin Smooth\n");
    tpd = taubin_smooth(pd, weights, 0.5, -0.53, atoi(argv[5]));

    vertices_staircase->Delete();
    weights->Delete();
    
    vtkXMLPolyDataWriter *writer = vtkXMLPolyDataWriter::New();
    writer->SetInput(pd);
    writer->SetFileName("saida.vtp");
    writer->Write();
    
    vtkPLYWriter *stl_writer = vtkPLYWriter::New();
    stl_writer->SetInput(tpd);
    stl_writer->SetFileName(argv[6]);
    stl_writer->Write();

    return 0;
}

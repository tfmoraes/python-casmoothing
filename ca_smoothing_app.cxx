#include "ca_smoothing.h"

vtkPolyData* read_stl(char* filename) {
    vtkSTLReader *stl_reader = vtkSTLReader::New();
    stl_reader->SetFileName(filename);
    stl_reader->Update();
    return stl_reader->GetOutput();
}

int main(int argc, char *argv[])
{
    vtkPolyData *stl, *nm, *cl, *pd, *tpd;
    vtkPolyDataNormals *normals;
    vtkCleanPolyData *clean;
    
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
    
    tpd = ca_smoothing(pd, atof(argv[2]), atof(argv[3]), atof(argv[4]),atoi(argv[5]));

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

#include <math.h>
#include <map>
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


int main(int argc, char *argv[])
{
    const double stack_orientation[3] = { 0, 0, 1 };
    vtkPolyData *stl, *nm, *cl, *pd, *tpd;
    vtkPolyDataNormals *normals;
    vtkCleanPolyData *clean;
    vtkIdList *vertices_staircase;
    vtkDoubleArray* weights;
    
    stl = read_stl(argv[1]);

    normals = vtkPolyDataNormals::New();
    normals->SetInput(stl);
    normals->ComputeCellNormalsOn();
    normals->Update();

    clean = vtkCleanPolyData::New();
    clean->SetInput(normals->GetOutput());
    clean->Update();

    pd = clean->GetOutput();
    pd->BuildLinks();

    vertices_staircase = find_staircase_artifacts(pd, stack_orientation, atof(argv[2]));
    weights = calc_artifacts_weight(pd, vertices_staircase, atof(argv[3]), atof(argv[4]));
    tpd = taubin_smooth(pd, weights, 0.5, -0.53, atoi(argv[5]));

    vertices_staircase->Delete();
    weights->Delete();
    
    vtkXMLPolyDataWriter *writer = vtkXMLPolyDataWriter::New();
    writer->SetInput(tpd);
    writer->SetFileName("saida.vtp");
    writer->Write();
    
    vtkPLYWriter *stl_writer = vtkPLYWriter::New();
    stl_writer->SetInput(tpd);
    stl_writer->SetFileName(argv[6]);
    stl_writer->Write();

    return 0;
}

vtkPolyData* read_stl(char* filename) {
    vtkSTLReader *stl_reader = vtkSTLReader::New();
    stl_reader->SetFileName(filename);
    stl_reader->Update();
    return stl_reader->GetOutput();
}

vtkIdList* find_staircase_artifacts(vtkPolyData* pd, const double stack_orientation[3], double T) {
    /*
    This function is used to find vertices at staircase artifacts, which are
    those vertices whose incident faces' orientation differences are
    aproximated 90°. The argument `pd' is a vtkPolydata, it have to be its
    faces orientation (normal) calculated before and make sure that
    BuildLinks() has been called.
    */
    int nv, nf, fid;
    double of, min, max;
    
    double *ni;
    vtkIdList *output = vtkIdList::New();
    vtkDoubleArray *scalars = vtkDoubleArray::New();
    vtkIdList *idfaces;//idfaces = vtk.vtkIdList()
    
    nv = pd->GetNumberOfPoints(); // Number of vertices.
    for (int vid=0; vid < nv; vid++){ //for vid in xrange(nv):
        idfaces = vtkIdList::New();//idfaces = vtk.vtkIdList()
    pd->GetPointCells(vid, idfaces); //pd.GetPointCells(vid, idfaces) # Getting faces connected to face vid.
        nf = idfaces->GetNumberOfIds();
    
    max = -1000;
    min = 1000;
    for (int nid=0; nid < nf; nid++) {
        fid = idfaces->GetId(nid);
        ni = pd->GetCellData()->GetArray("Normals")->GetTuple(fid);

        of = 1 - (ni[0]*stack_orientation[0] + ni[1]*stack_orientation[1] + ni[2]*stack_orientation[2]);

        if (of > max) max = of;
        if (of < min) min = of;
    }

        // Getting the ones which normals dot is 90°, its vertex is added to
        // output
    if (max - min >= T) {
        output->InsertNextId(vid);
        scalars->InsertNextValue(1);
    }
    else {
        scalars->InsertNextValue(0);
    }
    idfaces->Delete();
    }
    vtkPointData* pointData = pd->GetPointData();
    pointData->SetScalars(scalars);
    return output;
}

vtkIdList* get_near_vertices_to_v(vtkPolyData* pd, int v, double dmax){
    /*
    Returns all vertices with distance at most "d" to the vertice "v" with
    their distances.
    pd - vtkPolydata
    v - the reference vertice
    dmax - the maximun distance.
    */
    double vi[3], vj[3], d;
    int n=0, nf, fid;
    
    std::map <int, bool> status_v;
    std::queue <int> to_visit;

    vtkIdList* near_vertices = vtkIdList::New();
    vtkIdList* idfaces;

    pd->GetPoint(v, vi); // The position of vertex v

    while (1) {
        idfaces = vtkIdList::New();
        pd->GetPointCells(v, idfaces);
        nf = idfaces->GetNumberOfIds();
        for(int nid=0; nid < nf; nid++) {
            fid = idfaces->GetId(nid);
            vtkCell* face = pd->GetCell(fid);

            for(int i=0; i < 3; i++) {
                int vjid = face->GetPointId(i);
                if (status_v.find(vjid) == status_v.end() || !status_v[vjid]) {
                    pd->GetPoint(vjid, vj);
                    d = sqrt((vi[0] - vj[0]) * (vi[0] - vj[0])\
                            + (vi[1] - vj[1]) * (vi[1] - vj[1])\
                            + (vi[2] - vj[2]) * (vi[2] - vj[2]));
                    if (d <= dmax) {
                        near_vertices->InsertNextId(vjid);
                        to_visit.push(vjid);
                    }
                }
                status_v[vjid] = true;
            }
        }

        n++;

        if (to_visit.empty())
            break;
        v = to_visit.front();
        to_visit.pop();
        idfaces->Delete();
    }

    return near_vertices;
}

vtkDoubleArray* calc_artifacts_weight(vtkPolyData* pd, vtkIdList* vertices_staircase, double tmax, double bmin) {
    /*
    Calculate the artifact weight based on distance of each vertex to its
    nearest staircase artifact vertex.
    pd - vtkPolydata;
    vertices_staircase - the identified staircase artifact vertices;
    tmax=2 - max distance the vertex must be to its nearest artifact vertex to
             considered to calculate the weight;
    bmin=0.1 - The minimun weight.
    */
    double vi[3], vj[3], d, value;
    int viid, vjid, nnv;
    int nid = vertices_staircase->GetNumberOfIds();
    vtkIdList* near_vertices;
    vtkDoubleArray* weights = vtkDoubleArray::New();
    vtkDataArray* scalars = pd->GetPointData()->GetScalars();
    for (int i=0; i < pd->GetNumberOfPoints(); i++){
        weights->InsertNextValue(0);
    }
    for (int i=0; i < nid; i++) {
        viid = vertices_staircase->GetId(i);
        pd->GetPoint(viid, vi);
        near_vertices = get_near_vertices_to_v(pd, viid, tmax);
        nnv = near_vertices->GetNumberOfIds();
        for (int j=0; j < nnv; j++){
            vjid = near_vertices->GetId(j);
            pd->GetPoint(vjid, vj);
            d = sqrt((vi[0] - vj[0]) * (vi[0] - vj[0])\
                    + (vi[1] - vj[1]) * (vi[1] - vj[1])\
                    + (vi[2] - vj[2]) * (vi[2] - vj[2]));
            value = (1.0 - d/tmax) * (1 - bmin) + bmin;
            if (value > weights->GetValue(vjid)) {
                printf("%f\n", value);
                weights->SetValue(vjid, value);
                scalars->SetTuple1(vjid, value);
            }
        }
        near_vertices->Delete();
    }

    vtkPointData* pointData = pd->GetPointData();
    pointData->SetScalars(scalars);
    return weights;
}

Point calc_d(vtkPolyData* pd, int vid){
    Point D;
    int nf, fid, n=0;
    double vi[3], vj[3];
    std::set<int> vertices;
    std::set<int>::iterator it;
    vtkIdList* idfaces = vtkIdList::New();
    pd->GetPointCells(vid, idfaces);
    nf = idfaces->GetNumberOfIds();
    for (int nid=0; nid < nf; nid++) {
        fid = idfaces->GetId(nid);
        vtkCell* face = pd->GetCell(fid);
        for (int i=0; i < 3; i++) {
            int vjid = face->GetPointId(i);
            vertices.insert(vjid);
        }
    }
    D.x = 0;
    D.y = 0;
    D.z = 0;

    pd->GetPoint(vid, vi); // The position of vertex v 
    for (it=vertices.begin(); it!=vertices.end(); it++) {
        pd->GetPoint(*it, vj);
        D.x = D.x + (vi[0] - vj[0]);
        D.y = D.y + (vi[1] - vj[1]);
        D.z = D.z + (vi[2] - vj[2]);
        n++;
    }

    D.x = D.x / n;
    D.y = D.y / n;
    D.z = D.z / n;
    return D;
}

vtkPolyData* taubin_smooth(vtkPolyData* pd, vtkDoubleArray* weights, double l, double m, int steps){
    double vi[3];
    vtkPolyData* new_pd = vtkPolyData::New();
    new_pd->DeepCopy(pd);
    //std::vector<Point> D(pd->GetNumberOfPoints());
    //D.reserve(pd->GetNumberOfPoints());
    Point *D;
    D = (Point*) malloc(pd->GetNumberOfPoints() * sizeof(Point));
    
    for (int s=0; s < steps; s++) {
        for (int i=0; i < pd->GetNumberOfPoints(); i++) {
            D[i] = calc_d(new_pd, i);
        }

        for (int i=0; i < pd->GetNumberOfPoints(); i++) {
            new_pd->GetPoint(i, vi);
            vi[0] = vi[0] + weights->GetValue(i)*l*D[i].x;
            vi[1] = vi[1] + weights->GetValue(i)*l*D[i].y;
            vi[2] = vi[2] + weights->GetValue(i)*l*D[i].z;
            new_pd->GetPoints()->SetPoint(i, vi);
        }

        for (int i=0; i < pd->GetNumberOfPoints(); i++) {
            D[i] = calc_d(new_pd, i);
        }

        for (int i=0; i < pd->GetNumberOfPoints(); i++) {
            new_pd->GetPoint(i, vi);
            vi[0] = vi[0] + weights->GetValue(i)*m*D[i].x;
            vi[1] = vi[1] + weights->GetValue(i)*m*D[i].y;
            vi[2] = vi[2] + weights->GetValue(i)*m*D[i].z;
            new_pd->GetPoints()->SetPoint(i, vi);
        }
    }
    free(D);
    return new_pd;
}

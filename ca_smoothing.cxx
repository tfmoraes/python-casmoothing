#include "ca_smoothing.h"

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
    double of_z, of_y, of_x, min_z, max_z, min_y, max_y, min_x, max_x;
    
    double *ni;
    vtkIdList *output = vtkIdList::New();
    vtkDoubleArray *scalars = vtkDoubleArray::New();
    vtkIdList *idfaces;//idfaces = vtk.vtkIdList()
    
    nv = pd->GetNumberOfPoints(); // Number of vertices.
    for (int vid=0; vid < nv; vid++){ //for vid in xrange(nv):
        idfaces = vtkIdList::New();//idfaces = vtk.vtkIdList()
        pd->GetPointCells(vid, idfaces); //pd.GetPointCells(vid, idfaces) # Getting faces connected to face vid.
        nf = idfaces->GetNumberOfIds();
    
        max_z = -1000;
        min_z = 1000;
        max_y = -1000;
        min_y = 1000;
        max_x = -1000;
        min_x = 1000;
        for (int nid=0; nid < nf; nid++) {
            fid = idfaces->GetId(nid);
            ni = pd->GetCellData()->GetArray("Normals")->GetTuple(fid);

            of_z = 1 - abs(ni[0]*stack_orientation[0] + ni[1]*stack_orientation[1] + ni[2]*stack_orientation[2]);
            of_y = 1 - abs(ni[0]*0 + ni[1]*1 + ni[2]*0);
            of_x = 1 - abs(ni[0]*1 + ni[1]*0 + ni[2]*0);

            if (of_z > max_z) max_z = of_z;
            if (of_z < min_z) min_z = of_z;

            if (of_y > max_y) max_y = of_y;
            if (of_y < min_y) min_y = of_y;

            if (of_x > max_x) max_x = of_x;
            if (of_x < min_x) min_x = of_x;
        }

        // Getting the ones which normals dot is 90°, its vertex is added to
        // output
        if ((abs(max_z - min_z) >= T) || (abs(max_y - min_y) >= T) || (abs(max_x - min_x) >= T)) {
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
    
    std::unordered_map <int, bool> status_v;
    std::unordered_map <int, bool> status_f;
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
            if (status_f.find(fid) != status_f.end()) {
                if (to_visit.empty())
                    break;
                v = to_visit.front();
                to_visit.pop();
                continue;
            }
            status_f[fid] = true;
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
            if (vjid != vid)
                vertices.insert(vjid);
        }
    }
    idfaces->Delete();
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
    /*
     * Implementation of Taubin's smooth algorithm described in the paper "A
     * Signal Processing Approach To Fair Surface Design". His benefeat is it
     * avoids to 
     *
     */
    double vi[3];
    vtkPolyData* new_pd = vtkPolyData::New();
    new_pd->DeepCopy(pd);
    //std::vector<Point> D(pd->GetNumberOfPoints());
    //D.reserve(pd->GetNumberOfPoints());
    Point *D;
    D = (Point*) malloc(pd->GetNumberOfPoints() * sizeof(Point));
    
    for (int s=0; s < steps; s++) {
        printf("Step %d\n", s);
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

std::unordered_map<int, std::vector<int>> get_flat_areas(vtkPolyData* pd, double stack_orientation[3], double tmax, int msize) {
    /* Returns each flat area in the given *vtkPolydata pd*. To be considered a
     * flat area it must have normal bellow *tmax* related to the
     * *stack_orientation* and have at least *msize* triangles. The result is
     * returned as hash_map.
     */
    std::unordered_map<int, std::vector<int>> flat_areas;
    return flat_areas;
}

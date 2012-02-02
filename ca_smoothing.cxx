#include <math.h>
#include <map>
#include <queue>
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
#include <vtkPolyDataNormals.h>
#include <vtkCleanPolyData.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkCell.h>

vtkPolyData* read_stl(char*);
vtkIdList* find_staircase_artifacts(vtkPolyData*, const double[3], double);
vtkIdList* get_near_vertices_to_v(vtkPolyData*, int, double);
vtkDoubleArray* calc_artifacts_weight(vtkPolyData*, vtkIdList*, double, double);

int main(int argc, char *argv[])
{
	const double stack_orientation[3] = { 0, 0, 1 };
	vtkPolyData *stl, *nm, *cl, *pd;
	vtkPolyDataNormals *normals;
	vtkCleanPolyData *clean;
	vtkIdList *vertices_staircase;
	
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
    calc_artifacts_weight(pd, vertices_staircase, atof(argv[3]), 1);
	
	vtkXMLPolyDataWriter *writer = vtkXMLPolyDataWriter::New();
	writer->SetInput(pd);
	writer->SetFileName("/tmp/saida.vtp");
	writer->Write();
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
                    if (status_v.find(vjid) == status_v.end() or !status_v[vjid]) {
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
    }

    vtkPointData* pointData = pd->GetPointData();
    pointData->SetScalars(scalars);
    return weights;
}
/*
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

def calc_artifacts_weight(pd, vertices_staircase, tmax=5.0, bmin=1):
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

    new_pd.SetPoints(points)
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
*/

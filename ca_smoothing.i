%module ca_smoothing
/*%include "ca_smoothing.h"*/
%include exception.i

// Common stuff
#include "vtkObjectBase.h"
#include "vtkObject.h"
#include "vtkPolyData.h"

%include cpointer.i
%pointer_functions(int, intp);
%pointer_functions(float, floatp);
%pointer_functions  (double, doublep);
%pointer_functions  (bool, boolp);

/*%typemap(cstype) vtkPolyData * "vtkPolyData"*/
/*%typemap(csin) vtkPolyData * "$csinput.GetCppThis()"*/
/*%typemap(csout) (vtkPolyData*) {*/
  /*IntPtr rawCppThisSwig = $imcall;*/
  /*vtkPolyData data = new vtkPolyData( rawCppThisSwig, false, false );*/
  /*return data;*/
/*}*/

%{
    #include "vtkPythonUtil.h"
    #include "vtkPolyData.h"
%}

%typemap(out) vtkPolyData* {
  PyImport_ImportModule("vtk");
  $result = vtkPythonUtil::GetObjectFromPointer ( (vtkPolyData*)$1 );
}

%typemap(in) vtkPolyData* {
  /*$1 = NULL;*/
  $1 = (vtkPolyData*) vtkPythonUtil::GetPointerFromObject ( $input, "vtkPolyData" );
  if ( $1 == NULL ) { SWIG_fail; }
}

extern vtkPolyData* ca_smoothing (vtkPolyData*, double, double, double, int);
extern vtkPolyData* get_sphere (void);

%{
   #include "stdio.h"
   #include "math.h"
   #include "ca_smoothing.h"
%}

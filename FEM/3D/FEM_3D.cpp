#include "FEM_3D.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <memory>   
namespace py = pybind11;

/////////////////////////////////////////////////////////////////////////////////////
FEM_3D::FEM_3D() {}
/////////////////////////////////////////////////////////////////////////////////////
FEM_3D::~FEM_3D() {}
/////////////////////////////////////////////////////////////////////////////////////
void FEM_3D::main_func() {
    std::cout << "FEM_3D::main_func()" << std::endl;
    std::vector<std::shared_ptr<ELEMENT_HEXAHEDRAL_EDGE>> elements;
    std::vector<std::shared_ptr<NODE>> nodes;
    std::vector<std::shared_ptr<EDGE>> edges;

    
}
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
PYBIND11_MODULE(FEM_3D, m){
    m.doc() = "FEM_3D";
    py::class_<FEM_3D, std::shared_ptr<FEM_3D>>(m, "FEM_3D")
        .def(py::init<>())
        .def("main_func", &FEM_3D::main_func);
}
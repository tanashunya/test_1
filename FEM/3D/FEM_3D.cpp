#include "node.h"
#include "element.h"
#include <vector>
#include <memory>
#include <pybind11/pybind11.h>
#include <iostream>

namespace py = pybind11;

class FEM_3D{
    private:
    public:
        FEM_3D();
        ~FEM_3D();
        void main_func();
        std::vector<std::shared_ptr<NODE>> nodes;
        std::vector<std::shared_ptr<ELEMENT>> elements;
};
/////////////////////////////////////////////////////////////////////////////////////
FEM_3D::FEM_3D() {}
/////////////////////////////////////////////////////////////////////////////////////
FEM_3D::~FEM_3D() {}
/////////////////////////////////////////////////////////////////////////////////////
void FEM_3D::main_func() {
    std::cout << "FEM_3D" << std::endl;
}
/////////////////////////////////////////////////////////////////////////////////////
PYBIND11_MODULE(FEM_3D, m){
    m.doc() = "FEM_3D";
    py::class_<FEM_3D, std::shared_ptr<FEM_3D>>(m, "FEM_3D")
        .def(py::init<>())
        .def("main_func", &FEM_3D::main_func);
}
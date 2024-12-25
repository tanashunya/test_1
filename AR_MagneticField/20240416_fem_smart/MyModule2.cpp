#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include "magic.h"
#include "const.h"
#include "GenMatrix.h"
#include "MyFunction.h"
namespace py = pybind11;

class Cfunctions{
public:
    Cfunctions();
    ~Cfunctions();
    // void write_vtk(int num_row, int num_col, py::array_t<double> Bx, py::array_t<double> By,  std::string& DATA_PASS);
    void write_vtk(int num_row, int num_col, py::array_t<double> Bx, py::array_t<double> By,  std::string& path);
};

Cfunctions::Cfunctions() {

}

Cfunctions::~Cfunctions() {

}

void Cfunctions::write_vtk(int num_row, int num_col, py::array_t<double> Bx, py::array_t<double> By,  std::string& path){
    // std::cout << "Function: write_vtk" << std::endl;
    int num_cells = num_row*num_col;
    FILE *fp;
    fp = fopen(path.c_str(), "w");
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "B_vector\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    
    
    fprintf(fp,"DIMENSIONS %d %d 1\n", num_row+1, num_col+1);
    fprintf(fp,"ORIGIN %lf %lf %lf\n",0.0, 0.0, 0.0);
    fprintf(fp,"SPACING %lf %lf 1\n",8.0, 8.0);
    // Bx By 書き込み
    auto Bx_ptr = Bx.data(); // Bxのポインタを取得
    auto By_ptr = By.data(); // Byのポインタを取得
    fprintf(fp, "CELL_DATA %d\n", num_cells);
    fprintf(fp, "VECTORS B float\n");
    for(int i=0; i<num_cells; i++){
        fprintf(fp, "%.15lf %.15lf %.15lf\n", Bx_ptr[i], By_ptr[i], 0.0);
    }
    
    

    fclose(fp);
}


PYBIND11_MODULE(Cfunctions, m) {
    m.doc() = "This is a test module";
    py::class_<Cfunctions>(m, "Cfunctions")
		.def(py::init<>())
		.def("write_vtk", &Cfunctions::write_vtk);
}
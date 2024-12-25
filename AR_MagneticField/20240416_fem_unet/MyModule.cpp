#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include "magic.h"
#include "const.h"
#include "GenMatrix.h"
#include "MyFunction.h"
#include <vector>
#include <memory>
namespace py = pybind11;
//********************************************************************************
class Carray{
private:
    int n; // 磁石の数
    // double* theta_of_magnets; // 磁化方向の配列
    std::vector<double> theta_of_magnets;
public:
    Carray(double theta, int n, py::array_t<double> np_array, int width, int height, bool MAKE_DATASET, int data_index, const std::string& DATA_PASS, const std::string& MESH_PASS);
    ~Carray();
    int len();
    double sum();
    int read_mesh();
    int read_akimacsv();
    int analysis(bool nonliner);
    int make_data();
    int write_vtk(const std::string& VTK_TYPE);
    
    double theta;
    int width;
    int height;
    bool MAKE_DATASET;
    int data_index;
    std::string DATA_PASS;
    std::string MESH_PASS;
    int num_nodes;
    int num_elements;
    double akima_data[38][5];
    const char *akima_csv_name = "akima4.csv";
    std::vector<std::unique_ptr<ELEMENT>> elements; 
    std::unique_ptr<double[]> A;
    std::unique_ptr<double[]> del_A;
    std::unique_ptr<double[]> r;         
    std::unique_ptr<double[]> b;         
    std::vector<std::unique_ptr<NODE>> nodes;       
    int num_OnBoundNodes;
    
    int max_nrloop = 20;
    double M = 1.25;
    double ICCG_conv = 1E-6;
    double NR_conv = 1E-3;
    double ICCG_conv_MAX;
    // bool printinfo = false;
    bool printinfo = true;

    const char *name_vtk_vec = "vectors.vtk";
    const char *name_vtk_mat = "material.vtk";
    const char *name_vtk_con = "contour.vtk";
};
//********************************************************************************
Carray::Carray(double theta, int n, py::array_t<double> np_array, int width, int height, bool MAKE_DATASET, int data_index, const std::string& DATA_PASS, const std::string& MESH_PASS) {
    this->n = n;
    this->theta = theta;
    this->width = width;
    this->height = height;
    this->MAKE_DATASET = MAKE_DATASET;
    this->data_index = data_index;
    this->DATA_PASS = DATA_PASS;
    this->MESH_PASS = MESH_PASS;
    for(int i = 0; i < n; i++) {
        double tmp = np_array.mutable_at(i);
        this->theta_of_magnets.push_back(tmp);
    }
}
//********************************************************************************
Carray::~Carray() {
    
} 
//********************************************************************************
int Carray::read_mesh(){
    using namespace std;

    char buf[256];
    FILE *fp;
    const char *meshfilename = (MESH_PASS).c_str();

    fp = fopen(meshfilename, "r");
    if(fp == NULL)
    {
        printf("File open error.(meshfile)\n");
        return 1;
    }
    if(fscanf(fp, "%*[^\n]%*c") == 0){} //１行読み飛ばし
    if(fscanf(fp, "%d%*c", &num_nodes)==0){} //ノード数読み込み
    nodes.reserve(num_nodes);
    for (int i=0; i<num_nodes; i++) // nodeデータ格納
    {   
        int node_index;
        double node_x, node_y, node_z;
        if (fscanf(fp, "%d %lf %lf %lf", &node_index, &node_x, &node_y, &node_z) == 0){}
        
        auto node = std::make_unique<NODE>();
        node->setNode(node_index-1, node_x, node_y, node_z);
        nodes.push_back(std::move(node));
    }
    if (fgets(buf, sizeof(buf), fp) == 0); // ２行読み飛ばし?
    if (fgets(buf, sizeof(buf), fp) == 0); 
    if (fgets(buf, sizeof(buf), fp) == 0);
    if (fscanf(fp, "%d%*c", &num_elements) == 0){} // エレメント数読み込み  
    elements.reserve(num_elements);
    for (int i = 0; i < num_elements; i++) {
        int elemNum, matNum, node_idx1, node_idx2, node_idx3;
        fscanf(fp, "%d %*d %d %*d %*d %d %d %d", &elemNum, &matNum, &node_idx1, &node_idx2, &node_idx3);
        auto element = std::make_unique<ELEMENT>();
        element->setElement(elemNum-1, matNum, *nodes[node_idx1-1].get(), *nodes[node_idx2-1].get(), *nodes[node_idx3-1].get());
        elements.push_back(std::move(element));
    }
    // for (int i=0; i<num_elements; i++){
    //     cout << elements[i]->nodes[1].index << endl;
    // }
    fclose(fp);
    return 0;
}
//********************************************************************************
int Carray::read_akimacsv(){
    // cout << "This is Carray::read_akimacsv" << endl;
    FILE *fp;
    double tmp;

    fp = fopen(this->akima_csv_name, "r");
    if (fp==NULL)
    {
        cout << "File open error.(akima_data)" << endl;
    }
    for (int i=0; i<38; i++)
    {
        for (int j=0; j<5; j++)
        {   
            if (fscanf(fp, "%lf,", &akima_data[i][j])==1){}
        }
        // cout << i << ":  " <<  akima_data[i][0] << "   " << akima_data[i][1] << "   " << akima_data[i][2] << "   " << akima_data[i][3] << "   " << akima_data[i][4] << endl;
    }
    fclose(fp);
    return 0;
}
//********************************************************************************
int Carray::analysis(bool nonliner){
    int num_boundnodes = 0;
    int num_nonboundnodes = num_nodes - num_boundnodes;
    GenMatrix K(num_nodes, num_nodes);
    GenMatrix L(num_nonboundnodes, num_nonboundnodes);
    r = std::make_unique<double[]>(num_nonboundnodes);
    b = std::make_unique<double[]>(num_nodes); // 最終的な解ベクトル
    del_A = std::make_unique<double[]>(num_nodes);
    A = std::make_unique<double[]>(num_nodes);
    ICCG_conv_MAX = num_nodes;

    for (int n=0; n<num_nodes; n++){
        A[n] = 0.0;
    }
    
    // cとかdとかdeltaとか色々計算
    for (int e=0; e<num_elements; e++){
        elements[e]->cal_constants();
    }

    auto start1 = chrono::steady_clock::now();
    auto start2 = start1;

    // NR法
    for (int nrloop=0; nrloop<max_nrloop; nrloop++){
        // 配列初期化
        for (int n=0; n<num_nodes; n++){
            b[n] = 0.0;
            del_A[n] = 0.0;
            K.quickAllZero();
        }
        // K, b 作成
        for (int e=0; e<num_elements; e++){
            int mat = elements[e]->materialNumber;
            // U更新
            elements[e]->cal_U();
            for (int i=0; i<3; i++){
                for (int j=0; j<3; j++){
                    double tmp = elements[e]->nyu * elements[e]->S[i][j] + 2.0 * elements[e]->dnyu_dB2 * elements[e]->U[i] * elements[e]->U[j] / elements[e]->delta;
                    K.add(elements[e]->nodes[i].index, elements[e]->nodes[j].index, tmp);
                }
                b[elements[e]->nodes[i].index] += (elements[e]->J_0 * elements[e]->delta / 3.0) - (elements[e]->nyu * elements[e]->U[i]);
                if(4 <= mat && mat <= 10){
                    b[elements[e]->nodes[i].index] += elements[e]->nyu * (M * cos(theta_of_magnets[mat-4]) * elements[e]->d[i] - M * sin(theta_of_magnets[mat-4]) * elements[e]->c[i]) / 2.0;
                }
            }
        }

        // 固定境界の行番号を格納する
        // v1.push_back(番号);で入れていく
        vector<int> v1;

        K.setBoundaryCondition(b.get(), v1, 0.0);
        GenMatrix::icdcmp(K, L, 1.05);
        GenMatrix::iccgSolv(K, L, b.get(), del_A.get(), r.get(), ICCG_conv, ICCG_conv_MAX);
        
        for (int n=0; n<num_nodes; n++){
            A[n] += del_A[n];
            nodes[n]->add_A(del_A[n]);
        }
        for (int e=0; e<num_elements; e++){
            int index1 = elements[e]->nodes[0].index;
            int index2 = elements[e]->nodes[1].index;
            int index3 = elements[e]->nodes[2].index;
            elements[e]->update_node(*nodes[index1].get(), *nodes[index2].get(), *nodes[index3].get());
        }
        // 時間計測
        auto end = chrono::steady_clock::now();

        // 非線形性考慮しない場合はNRloop抜ける
        if (nonliner == false){
            for (int e=0; e<num_elements; e++){
                elements[e]->update_B();
            }
            if(printinfo == true){
                cout << "time[ms]: " << chrono::duration_cast<chrono::milliseconds>(end - start1).count() << endl;
            }
            break;
        }else{
            // delta_Aの大きさ計算
            double norm_delta_A = 0.0;
            for (int n=0; n<num_nodes; n++){
                norm_delta_A += del_A[n] * del_A[n];
            }
            norm_delta_A = sqrt(norm_delta_A);
            
            // 収束したらNRloop抜ける
            if (norm_delta_A < NR_conv){
                if(printinfo == true){
                    cout << "NR loop: " << nrloop + 1 
                         << "   norm_delta_A: " << setw(12) << left <<  norm_delta_A
                         << "   time[ms]: " << setw(4) << left <<  chrono::duration_cast<chrono::milliseconds>(end - start2).count()
                         << "   (終了)\n";
                }
                for (int e=0; e<num_elements; e++){
                    elements[e]->update_B();
                }
                break;
            }else{
                if(printinfo == true){
                    cout << "NR loop: " << nrloop + 1 
                         << "   norm_delta_A: " << setw(12) << left << norm_delta_A
                         << "   time[ms]: " << setw(4) << left <<  chrono::duration_cast<chrono::milliseconds>(end - start2).count()
                         << "   (継続)\n";
                }
                start2 = chrono::steady_clock::now();
                
                // B^2, nyu, dnyu/dB^2計算（要素ごと）
                for (int e=0; e<num_elements; e++){
                    elements[e]->update_B();
                    elements[e]->update_nyu(akima_data);
                }
            }
        }
    }
    // cout << "This is end of analysis()." << endl;

    // write_vtk("vector");

    return 0;
}
//********************************************************************************
int Carray::make_data(){
    // // 深層学習用出力側データ作成関数(20240328現在動作未確認)
    // if(this->MAKE_DATASET ==1){
    //     double x[3];
    //     double y[3];
    //     // ピクセルごとのBx, Byを求めcsvへ
    //     FILE *fp;
    //     const char *fout = (this->DATA_PASS).c_str();
    //     fp = fopen(fout, "w");
    //     if (fp == NULL){
    //         cout << fout << "ファイルが開けませんでした。" << endl;
    //         return -1;
    //     }
    //     for(int j=0; j<this->height; j++)
    //     {
    //         for(int k=0; k<this->width; k++)
    //         {   
    //             double px = k + 0.5; // ピクセル中心のx座標
    //             double py = this->height - ( j + 0.5 );
    //             // このピクセルがどの要素内にあるかを探索する
    //             for(int e=0; e<num_elements; e++)
    //             {
    //                 for (int i = 0; i < 3; i++)
    //                 {
    //                     x[i] = elements[e]->nodes[i].x;
    //                     y[i] = elements[e]->nodes[i].y;
    //                     // メモ：この要素は(x1, y1), (x2, y2), (x3, y3)に囲まれてる
    //                 }
    //                 double abXat = (x[1]-x[0])*(py-y[0])-(y[1]-y[0])*(px-x[0]);
    //                 double bcXbt = (x[2]-x[1])*(py-y[1])-(y[2]-y[1])*(px-x[1]);
    //                 double caXct = (x[0]-x[2])*(py-y[2])-(y[0]-y[2])*(px-x[2]);
    //                 if( abXat == 0.0 || bcXbt == 0.0 || caXct == 0.0){
    //                     fprintf(fp, "%.16lf,%.16lf,", elements[e]->Bx, elements[e]->By);
    //                     break;
    //                 }else if(( abXat > 0.0 && bcXbt > 0.0 && caXct > 0.0) || ( abXat < 0.0 && bcXbt < 0.0 && caXct < 0.0)){
    //                     fprintf(fp, "%.16lf,%.16lf,", elements[e]->Bx, elements[e]->By);
    //                     break;
    //                 }
    //                 if(e==num_elements-1){
    //                     cout << "error! : ピクセル（" << k + 1 << "," << j + 1 << "）は、どの要素にも属しません" << endl;
    //                     return 1;
    //                 }
    //             }
    //         }
    //         fprintf(fp, "\n");
    //     }
    //     fclose(fp);
    // }
    return 0;
}
//********************************************************************************
int Carray::write_vtk(const std::string& VTK_TYPE){
    const char *type = (VTK_TYPE).c_str();
    if(VTK_TYPE == "material"){
        write_vtk_material_2(name_vtk_mat, nodes, elements, num_nodes, num_elements);
    }else if (VTK_TYPE == "vector"){
        write_vtk_vector_2(name_vtk_vec, nodes, elements, num_nodes, num_elements);
    }else if (VTK_TYPE == "contour"){
        write_vtk_contour_2(name_vtk_con, nodes, elements, num_nodes, num_elements, A.get());
    }else{
        cout << "VTKファイルの書き込みができませんでした。" << endl;
        return 1;
    }
    return 0;
}
//********************************************************************************
int Carray::len() {
    return this->n;
}
//********************************************************************************
double Carray::sum() {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += theta_of_magnets[i];
    }
    return sum;
}
//********************************************************************************


PYBIND11_MODULE(Carray, m) {
    m.doc() = "This is a test module";
    py::class_<Carray>(m, "Carray")
		.def(py::init<double, int,py::array_t<double>, int, int, bool, int, const std::string&, const std::string&>())
		.def("len", &Carray::len)
        .def("sum", &Carray::sum)
        .def("read_mesh", &Carray::read_mesh)
        .def("read_akimacsv", &Carray::read_akimacsv)
        .def("analysis", &Carray::analysis)
        .def("write_vtk", &Carray::write_vtk);
}

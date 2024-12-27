#include <vector>
#include <iostream>
#include "MyFunction.h"
#include "magic.h"
#include "const.h"
using namespace std;
// *********************************************************************************************
ELEMENT::ELEMENT(){
    
}
void ELEMENT::setElement(int index, int matNum, std::shared_ptr<NODE> node1, std::shared_ptr<NODE> node2, std::shared_ptr<NODE> node3){
    this->index = index;
    this->materialNumber = matNum;
    this->nodes[0] = node1;
    this->nodes[1] = node2;
    this->nodes[2] = node3;
}
void ELEMENT::print_info(){
    // cout << d[0] << endl;
    cout << "要素番号: " << index << "  材料番号: " << materialNumber << "  node[0].index: " << this->nodes[0]->index<< "  node[0].x: " << this->nodes[0]->x << endl;
}
void ELEMENT::set_material_config(){
    // 材料番号ごとに設定

}
void ELEMENT::cal_constants(){
    dnyu_dB2 = 0.0;
    for (int i=0; i<3; i++){
        this->c[i] = nodes[(i+1)%3]->y - nodes[(i+2)%3]->y;
        this->d[i] = nodes[(i+2)%3]->x - nodes[(i+1)%3]->x;
    }
    this->delta = (nodes[0]->x * nodes[1]->y
                 + nodes[1]->x * nodes[2]->y
                 + nodes[2]->x * nodes[0]->y
                 - nodes[0]->y * nodes[1]->x
                 - nodes[1]->y * nodes[2]->x
                 - nodes[2]->y * nodes[0]->x) / 2.0;
    // this->delta = (nodes[0].x * nodes[1].y + nodes[1].x * nodes[2].y + nodes[2].x * nodes[0].y - nodes[0].y * nodes[1].x - nodes[1].y * nodes[2].x - nodes[2].y * nodes[0].x) / 2.0;
    // cout << "element index: " << index << 
    //         "  nodes[0].x: " << nodes[0].x << 
    //         "  delta: " << delta <<  endl;

    if(this->delta < 0){
        cout << "error!!! delta<0" << endl;
    }

    // S計算
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            this->S[i][j] = (c[i] * c[j] + d[i] * d[j]) / (4.0 * delta);
        }
    }
    
        
    // 電流条件付与
    if (materialNumber == 0 || materialNumber == 1 || materialNumber == 2) {
        // 空気、コア、磁石
        this->J_0 = 0.0;
    } else if (materialNumber == 31 || materialNumber == 32) {
        // U+
        this->J_0 = 0.0;  // 正の電流値を設定
    } else if (materialNumber == -31 || materialNumber == -32) {
        // U-
        this->J_0 = 0.0;  // 負の電流値を設定
    } else if (materialNumber == -41 || materialNumber == -42) {
        // V-
        this->J_0 = 0.0;  // 負の電流値を設定
    } else if (materialNumber == 51 || materialNumber == 52) {
        // W+
        this->J_0 = 0.0;  // 正の電流値を設定
    } else {
        cout << "J_0 エラー: 未定義の材料番号 " << materialNumber << endl;
    }

    // nyuの初期値設定
    if (materialNumber == 0){ 
        // 空気
        this->nyu = 1.0 / (myu0);
    }
    else if (materialNumber == 1){
        // コア
        this->nyu = 1.0 / (myu0 * myur);
    }
    else if (materialNumber == 2){
        // 永久磁石
        this->nyu = 1.0 / (myu0 * 1.05);
    }
    else{
        // コイル
        this->nyu = 1.0 / (myu0);
    }
}
void ELEMENT::update_B(){
    this->Bx =        (d[0] * nodes[0]->Az + d[1] * nodes[1]->Az + d[2] * nodes[2]->Az) / (2.0 * delta);
    this->By = -1.0 * (c[0] * nodes[0]->Az + c[1] * nodes[1]->Az + c[2] * nodes[2]->Az) / (2.0 * delta);
    this->B2 = Bx*Bx + By*By;
}
void ELEMENT::update_nyu(std::vector<std::vector<double>>& akima_data){

    // for (const auto& row : akima_data) {
    //     for (const auto& value : row) {
    //         cout << value << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;

    if (true){ 
        // コアなら
        double B2_dash, a0, a1, a2, a3, del;
        for (int i = 0; i < akima_data.size()-1; i++){
            if (akima_data[i][0] <= B2 && B2 < akima_data[i + 1][0]){ // 端から探索 効率悪くね？
                B2_dash = akima_data[i][0];
                a0 = akima_data[i][1];
                a1 = akima_data[i][2];
                a2 = akima_data[i][3];
                a3 = akima_data[i][4];
                del = B2 - B2_dash;
                // cout << "element index: " << index << "  nyu: " << nyu << endl;
                this->nyu = a0 + a1 * del + a2 * del * del + a3 * del * del * del;
                this->dnyu_dB2 = a1 + 2 * a2 * del + 3 * a3 * del * del;
                // cout << "element index: " << index << "  nyu: " << nyu << endl;
                // cout << "i: " << i << endl;
                break;
            }
            if (i == akima_data.size() - 1) {
                del = B2 - akima_data[i][0];
                this->nyu = akima_data[i][1] + akima_data[i][2] * del;
                this->dnyu_dB2 = akima_data[i][2];
                cout << "element index: " << index << "  磁気飽和しました " << endl;
            }
        }
    }
}
void ELEMENT::cal_U(){
    for (int i=0; i<3; i++){
        this->U[i] = S[i][0] * nodes[0]->Az + S[i][1] * nodes[1]->Az + S[i][2] * nodes[2]->Az;
    }
    // cout << index << ": " << U[0] << endl;;
}
void ELEMENT::update_node(std::shared_ptr<NODE> node1, std::shared_ptr<NODE> node2, std::shared_ptr<NODE> node3){
    this->nodes[0] = node1;
    this->nodes[1] = node2;
    this->nodes[2] = node3;
}

// *********************************************************************************************
NODE::NODE(){
    this->Az = 0.0;

    this->isSlideSlave = false;
    this->isPeriodicSlave = false;
    this->isSlideMaster = false;
    this->isPeriodicMaster = false;
    this->isOnDirichletBoundary = false;
    this->isRotor = false;

    this->C_index = -1; // エラーチェック用
    this->PeriodicMasterIndex  = -2; // マスターなし（補間されない、スレイブの節点じゃない）
    this->SlideMasterIndex_1  = -3; // マスターなし（補間されない、スレイブの節点じゃない）
    this->SlideMasterIndex_2  = -4; // マスターなし（補間されない、スレイブの節点じゃない)

}
void NODE::setNode(int node_index, double node_x, double node_y, double node_z){
    this->index = node_index;
    this->x = node_x;
    this->y = node_y;
    this->z = node_z;
}
void NODE::add_A(double delta_A){
    // if(this->isOnBoundary==false){
    //     this->Az += delta_A;
    // }
    this->Az += delta_A;

}
void NODE::cal_r_theta(){
    r = sqrt(x*x+y*y);
    theta = atan2(y, x);
}
void NODE::set_isOnSlideBoundary(bool TF){
    this->isSlideSlave = TF;
}
void NODE::set_isOnPeriodicBoundary(bool TF){
    this->isPeriodicSlave = TF;
}
void NODE::set_isOnDirichletBoundary(bool TF){
    this->isOnDirichletBoundary = TF;
}
void NODE::add_theta(double theta){
    this->theta += theta;
    this->x = r * cos(this->theta);
    this->y = r * sin(this->theta);
}

// *********************************************************************************************
void swap(int *a, int *b){
    int temp = *a;
    *a = *b;
    *b = temp;
}

// 未実装
void Readmesh(const char *filename){
    FILE *fp = fopen(filename, "r");
    if (fp == NULL){
        printf("File open error.\n");
    }

    // 最初の行を読み飛ばす
    if (fscanf(fp, "%*[^\n]%*c") == 0){}
    


    fclose(fp);
}

// paraview用書き出し(材料分布)
void write_vtk_material_2(const char *filename, std::vector<std::unique_ptr<NODE>>& nodes, std::vector<std::unique_ptr<ELEMENT>>& elements, int numnode, int numele)
{
    FILE *fp;
    fp = fopen(filename, "w");
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "material\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(fp, "POINTS %d float\n", numnode);
    for (int i = 0; i < numnode; i++)
    {
        // x, y, z座標
        fprintf(fp, "%.15lf %.15lf %.15lf\n", nodes[i]->x, nodes[i]->y, nodes[i]->z);
    }

    fprintf(fp, "CELLS %d %d\n", numele, 4 * numele);
    for (int i = 0; i < numele; i++)
    {   
        // 節点番号0、節点番号1、節点番号2
        fprintf(fp, "3 %d %d %d\n", elements[i]->nodes[0]->index, elements[i]->nodes[1]->index, elements[i]->nodes[2]->index);
    }

    fprintf(fp, "CELL_TYPES %d\n", numele);
    for (int i = 0; i < numele; i++)
    {
        fprintf(fp, "5\n");
    }

    fprintf(fp, "CELL_DATA %d\n", numele);
    fprintf(fp, "SCALARS cell_data int\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int i = 0; i < numele; i++)
    {
        // 材料番号
        if (fabs(elements[i]->materialNumber) > 30)
        {
            fprintf(fp, "%d\n", 3);
        }
        else
        {
            fprintf(fp, "%d\n", elements[i]->materialNumber);
        }
    }
    fclose(fp);
}

// paraview用書き出し(ベクトル)
void write_vtk_vector_2(const char *filename, std::vector<std::unique_ptr<NODE>>& nodes, std::vector<std::unique_ptr<ELEMENT>>& elements, int numnode, int numele)
{
    FILE *fp;
    fp = fopen(filename, "w");
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "material\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(fp, "POINTS %d float\n", numnode);
    for (int i = 0; i < numnode; i++)
    {
        fprintf(fp, "%.15lf %.15lf %.15lf\n", nodes[i]->x, nodes[i]->y, nodes[i]->z);
    }

    fprintf(fp, "CELLS %d %d\n", numele, 4 * numele);
    for (int i = 0; i < numele; i++)
    {
        fprintf(fp, "3 %d %d %d\n", elements[i]->nodes[0]->index, elements[i]->nodes[1]->index, elements[i]->nodes[2]->index);
    }

    fprintf(fp, "CELL_TYPES %d\n", numele);
    for (int i = 0; i < numele; i++)
    {
        fprintf(fp, "5\n");
    }

    fprintf(fp, "CELL_DATA %d\n", numele);
    fprintf(fp, "VECTORS cell_data float\n");
    for (int i = 0; i < numele; i++)
    {
        fprintf(fp, "%.15lf %.15lf %.15lf\n", elements[i]->Bx, elements[i]->By, 0.0);
    }
    fclose(fp);
}

// paraview用書き出し(contour)
void write_vtk_contour_2(const char *filename, std::vector<std::unique_ptr<NODE>>& nodes, std::vector<std::unique_ptr<ELEMENT>>& elements, int numnode, int numele, double *A)
{
    FILE *fp;
    fp = fopen(filename, "w");
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "material\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(fp, "POINTS %d float\n", numnode);
    for (int i = 0; i < numnode; i++)
    {
        // x, y, z
        fprintf(fp, "%.15lf %.15lf %.15lf\n", nodes[i]->x, nodes[i]->y, nodes[i]->z);
    }

    fprintf(fp, "CELLS %d %d\n", numele, 4 * numele);
    for (int i = 0; i < numele; i++)
    {
        // 3(角形), 節点１、節点２、節点３
        fprintf(fp, "3 %d %d %d\n", elements[i]->nodes[0]->index, elements[i]->nodes[1]->index, elements[i]->nodes[2]->index);
    }

    fprintf(fp, "CELL_TYPES %d\n", numele);
    for (int i = 0; i < numele; i++)
    {
        fprintf(fp, "5\n");
    }
    ////////////////////////////////////////////
    fprintf(fp, "POINT_DATA %d\n", numnode);
    fprintf(fp, "SCALARS A float\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int i=0; i<numnode; i++)
    {
        // |A|
        fprintf(fp, "%.15lf\n", A[i]);
    }
    ////////////////////////////////////////////
    fprintf(fp, "CELL_DATA %d\n", numele);
    fprintf(fp, "VECTORS cell_data float\n");
    for (int i = 0; i < numele; i++)
    {
        fprintf(fp, "%.15lf %.15lf %.15lf\n", elements[i]->Bx, elements[i]->By, 0.0);
    }
    fclose(fp);
}

// paraview用書き出し(材料分布)v3
void write_vtk_material_3(const char *filename, std::vector<std::shared_ptr<NODE>>& nodes, std::vector<std::shared_ptr<ELEMENT>>& elements, int numnode, int numele)
{
    FILE *fp;
    fp = fopen(filename, "w");
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "material\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(fp, "POINTS %d float\n", numnode);
    for (int i = 0; i < numnode; i++)
    {
        // x, y, z座標
        fprintf(fp, "%.15lf %.15lf %.15lf\n", nodes[i]->x, nodes[i]->y, nodes[i]->z);
    }

    fprintf(fp, "CELLS %d %d\n", numele, 4 * numele);
    for (int i = 0; i < numele; i++)
    {   
        // 節点番号0、節点番号1、節点番号2
        fprintf(fp, "3 %d %d %d\n", elements[i]->nodes[0]->C_index, elements[i]->nodes[1]->C_index, elements[i]->nodes[2]->C_index);
    }

    fprintf(fp, "CELL_TYPES %d\n", numele);
    for (int i = 0; i < numele; i++)
    {
        fprintf(fp, "5\n");
    }

    fprintf(fp, "CELL_DATA %d\n", numele);
    fprintf(fp, "SCALARS cell_data int\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int i = 0; i < numele; i++)
    {
        // 材料番号
        fprintf(fp, "%d\n", elements[i]->materialNumber);
    }
    fclose(fp);
}

// paraview用書き出し(ベクトル)v3
void write_vtk_vector_3(const char *filename, std::vector<std::shared_ptr<NODE>>& nodes, std::vector<std::shared_ptr<ELEMENT>>& elements, int numnode, int numele)
{
    FILE *fp;
    fp = fopen(filename, "w");
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "material\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(fp, "POINTS %d float\n", numnode);
    for (int i = 0; i < numnode; i++)
    {
        fprintf(fp, "%.15lf %.15lf %.15lf\n", nodes[i]->x, nodes[i]->y, nodes[i]->z);
    }

    fprintf(fp, "CELLS %d %d\n", numele, 4 * numele);
    for (int i = 0; i < numele; i++)
    {
        fprintf(fp, "3 %d %d %d\n", elements[i]->nodes[0]->C_index, elements[i]->nodes[1]->C_index, elements[i]->nodes[2]->C_index);
    }

    fprintf(fp, "CELL_TYPES %d\n", numele);
    for (int i = 0; i < numele; i++)
    {
        fprintf(fp, "5\n");
    }

    fprintf(fp, "CELL_DATA %d\n", numele);
    fprintf(fp, "VECTORS cell_data float\n");
    for (int i = 0; i < numele; i++)
    {
        fprintf(fp, "%.15lf %.15lf %.15lf\n", elements[i]->Bx, elements[i]->By, 0.0);
    }
    fclose(fp);
}

// paraview用書き出し(contour)v3
void write_vtk_contour_3(const char *filename, std::vector<std::shared_ptr<NODE>>& nodes, std::vector<std::shared_ptr<ELEMENT>>& elements, int numnode, int numele, double *A)
{
    FILE *fp;
    fp = fopen(filename, "w");
    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "material\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(fp, "POINTS %d float\n", numnode);
    for (int i = 0; i < numnode; i++)
    {
        // x, y, z
        fprintf(fp, "%.15lf %.15lf %.15lf\n", nodes[i]->x, nodes[i]->y, nodes[i]->z);
    }

    fprintf(fp, "CELLS %d %d\n", numele, 4 * numele);
    for (int i = 0; i < numele; i++)
    {
        // 3(角形), 節点１、節点２、節点３
        fprintf(fp, "3 %d %d %d\n", elements[i]->nodes[0]->C_index, elements[i]->nodes[1]->C_index, elements[i]->nodes[2]->C_index);
    }

    fprintf(fp, "CELL_TYPES %d\n", numele);
    for (int i = 0; i < numele; i++)
    {
        fprintf(fp, "5\n");
    }
    ////////////////////////////////////////////
    fprintf(fp, "POINT_DATA %d\n", numnode);
    fprintf(fp, "SCALARS A float\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int i=0; i<numnode; i++)
    {
        // |A|
        fprintf(fp, "%.15lf\n", A[i]);
    }
    ////////////////////////////////////////////
    fprintf(fp, "CELL_DATA %d\n", numele);
    fprintf(fp, "VECTORS cell_data float\n");
    for (int i = 0; i < numele; i++)
    {
        fprintf(fp, "%.15lf %.15lf %.15lf\n", elements[i]->Bx, elements[i]->By, 0.0);
    }
    fclose(fp);
}
//********************************************************************************
int count_csv_lines(const char* filename) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        cout << filename << " could not be opened." << endl;
        return -1;
    }

    int line_count = 0;
    char buf[256];
    while (fgets(buf, sizeof(buf), fp) != NULL) {
        line_count++;
    }

    fclose(fp);
    return line_count;
}
#include "MyFunction.h"
#include "magic.h"
#include "const.h"
using namespace std;
// *********************************************************************************************
ELEMENT::ELEMENT(){
    
}
void ELEMENT::setElement(int index, int matNum, NODE& node1, NODE& node2, NODE& node3){
    this->index = index;
    this->materialNumber = matNum;
    this->nodes[0] = node1;
    this->nodes[1] = node2;
    this->nodes[2] = node3;
}
void ELEMENT::print_info(){
    // cout << d[0] << endl;
    cout << "要素番号: " << index << "  材料番号: " << materialNumber << "  node[0].index: " << this->nodes[0].index<< "  node[0].x: " << this->nodes[0].x << endl;
}
void ELEMENT::set_material_config(){
    // 材料番号ごとに設定

}
void ELEMENT::cal_constants(){
    dnyu_dB2 = 0.0;
    for (int i=0; i<3; i++){
        this->c[i] = nodes[(i+1)%3].y - nodes[(i+2)%3].y;
        this->d[i] = nodes[(i+2)%3].x - nodes[(i+1)%3].x;
    }
    this->delta = (nodes[0].x * nodes[1].y
                 + nodes[1].x * nodes[2].y
                 + nodes[2].x * nodes[0].y
                 - nodes[0].y * nodes[1].x
                 - nodes[1].y * nodes[2].x
                 - nodes[2].y * nodes[0].x) / 2.0;
    // this->delta = (nodes[0].x * nodes[1].y + nodes[1].x * nodes[2].y + nodes[2].x * nodes[0].y - nodes[0].y * nodes[1].x - nodes[1].y * nodes[2].x - nodes[2].y * nodes[0].x) / 2.0;
    // cout << "element index: " << index << 
    //         "  nodes[0].x: " << nodes[0].x << 
    //         "  delta: " << delta <<  endl;

    if(delta < 0){
        cout << "error!!! delta<0" << endl;
    }

    // S計算
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            this->S[i][j] = (c[i] * c[j] + d[i] * d[j]) / (4.0 * delta);
        }
    }
    
        
    // 電流条件付与
    if (1 <= materialNumber){
        // 全部
        this->J_0 = 0.0;
    }else{
        cout << "J_0 error" << endl;
    }

    // nyuの初期値設定
    if (materialNumber == 1){ 
        // 空気
        this->nyu = 1.0 / (myu0);
    }
    else if (materialNumber == 2){
        // 鉄
        this->nyu = 1.0 / (myu0 * myur);
        // cout << nyu << endl;
    }
    else if (materialNumber == 3){
        // コイル
        this->nyu = 1.0 / (myu0);
    }
    else if (4 <= materialNumber && materialNumber <= 10){
        // 永久磁石
        this->nyu = 1.0 / (myu0 * 1.05);
    }
    else {
        cout << "材料番号 " << materialNumber << " は定義されていません。" << endl;
    }
}
void ELEMENT::update_B(){
    this->Bx =        (d[0] * nodes[0].Az + d[1] * nodes[1].Az + d[2] * nodes[2].Az) / (2.0 * delta);
    this->By = -1.0 * (c[0] * nodes[0].Az + c[1] * nodes[1].Az + c[2] * nodes[2].Az) / (2.0 * delta);
    this->B2 = Bx*Bx + By*By;
}
void ELEMENT::update_nyu(double akima_data[38][5]){
    if (materialNumber == 2){ 
        // 鉄なら
        double B2_dash, a0, a1, a2, a3, del;
        for (int i = 0; i < 40; i++){
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
                break;
            }
        }
    }
}
void ELEMENT::cal_U(){
    for (int i=0; i<3; i++){
        this->U[i] = S[i][0] * nodes[0].Az + S[i][1] * nodes[1].Az + S[i][2] * nodes[2].Az;
    }
    // cout << index << ": " << U[0] << endl;;
}
void ELEMENT::update_node(NODE& node1, NODE& node2, NODE& node3){
    this->nodes[0] = node1;
    this->nodes[1] = node2;
    this->nodes[2] = node3;
}
// *********************************************************************************************
NODE::NODE(){
    this->Az = 0.0;
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
// void write_vtk_material(const char *filename, double **node, int **element, int numnode, int numele)
// {
//     FILE *fp;
//     fp = fopen(filename, "w");
//     fprintf(fp, "# vtk DataFile Version 2.0\n");
//     fprintf(fp, "material\n");
//     fprintf(fp, "ASCII\n");
//     fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

//     fprintf(fp, "POINTS %d float\n", numnode);
//     for (int i = 0; i < numnode; i++)
//     {
//         fprintf(fp, "%.15lf %.15lf %.15lf\n", node[i][1], node[i][2], node[i][3]);
//     }

//     fprintf(fp, "CELLS %d %d\n", numele, 4 * numele);
//     for (int i = 0; i < numele; i++)
//     {
//         fprintf(fp, "3 %d %d %d\n", element[i][2], element[i][3], element[i][4]);
//     }

//     fprintf(fp, "CELL_TYPES %d\n", numele);
//     for (int i = 0; i < numele; i++)
//     {
//         fprintf(fp, "5\n");
//     }

//     fprintf(fp, "CELL_DATA %d\n", numele);
//     fprintf(fp, "SCALARS cell_data int\n");
//     fprintf(fp, "LOOKUP_TABLE default\n");
//     for (int i = 0; i < numele; i++)
//     {
//         if (fabs(element[i][1]) > 30)
//         {
//             fprintf(fp, "%d\n", 3);
//         }
//         else
//         {
//             fprintf(fp, "%d\n", element[i][1]);
//         }
//     }
//     fclose(fp);
// }

// // paraview用書き出し(ベクトル)
// void write_vtk_vector(const char *filename, double **node, int **element, double *B_x, double *B_y, int numnode, int numele)
// {
//     FILE *fp;
//     fp = fopen(filename, "w");
//     fprintf(fp, "# vtk DataFile Version 2.0\n");
//     fprintf(fp, "material\n");
//     fprintf(fp, "ASCII\n");
//     fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

//     fprintf(fp, "POINTS %d float\n", numnode);
//     for (int i = 0; i < numnode; i++)
//     {
//         fprintf(fp, "%.15lf %.15lf %.15lf\n", node[i][1], node[i][2], node[i][3]);
//     }

//     fprintf(fp, "CELLS %d %d\n", numele, 4 * numele);
//     for (int i = 0; i < numele; i++)
//     {
//         fprintf(fp, "3 %d %d %d\n", element[i][2], element[i][3], element[i][4]);
//     }

//     fprintf(fp, "CELL_TYPES %d\n", numele);
//     for (int i = 0; i < numele; i++)
//     {
//         fprintf(fp, "5\n");
//     }

//     fprintf(fp, "CELL_DATA %d\n", numele);
//     fprintf(fp, "VECTORS cell_data float\n");
//     for (int i = 0; i < numele; i++)
//     {
//         fprintf(fp, "%.15lf %.15lf %.15lf\n", B_x[i], B_y[i], 0.0);
//     }
//     fclose(fp);
// }

// // paraview用書き出し(contour)
// void write_vtk_contour(const char *filename, double **node, int **element, double *B_x, double *B_y, int numnode, int numele, double *A)
// {
//     FILE *fp;
//     fp = fopen(filename, "w");
//     fprintf(fp, "# vtk DataFile Version 2.0\n");
//     fprintf(fp, "material\n");
//     fprintf(fp, "ASCII\n");
//     fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

//     fprintf(fp, "POINTS %d float\n", numnode);
//     for (int i = 0; i < numnode; i++)
//     {
//         // x, y, z
//         fprintf(fp, "%.15lf %.15lf %.15lf\n", node[i][1], node[i][2], node[i][3]);
//     }

//     fprintf(fp, "CELLS %d %d\n", numele, 4 * numele);
//     for (int i = 0; i < numele; i++)
//     {
//         // 3(角形), 節点１、節点２、節点３
//         fprintf(fp, "3 %d %d %d\n", element[i][2], element[i][3], element[i][4]);
//     }

//     fprintf(fp, "CELL_TYPES %d\n", numele);
//     for (int i = 0; i < numele; i++)
//     {
//         fprintf(fp, "5\n");
//     }
//     ////////////////////////////////////////////
//     fprintf(fp, "POINT_DATA %d\n", numnode);
//     fprintf(fp, "SCALARS A float\n");
//     fprintf(fp, "LOOKUP_TABLE default\n");
//     for (int i=0; i<numnode; i++)
//     {
//         // |A|
//         fprintf(fp, "%.15lf\n", A[i]);
//     }
//     ////////////////////////////////////////////
//     fprintf(fp, "CELL_DATA %d\n", numele);
//     fprintf(fp, "VECTORS cell_data float\n");
//     for (int i = 0; i < numele; i++)
//     {
//         fprintf(fp, "%.15lf %.15lf %.15lf\n", B_x[i], B_y[i], 0.0);
//     }
//     fclose(fp);
// }


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
        fprintf(fp, "3 %d %d %d\n", elements[i]->nodes[0].index, elements[i]->nodes[1].index, elements[i]->nodes[2].index);
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
        fprintf(fp, "3 %d %d %d\n", elements[i]->nodes[0].index, elements[i]->nodes[1].index, elements[i]->nodes[2].index);
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
        fprintf(fp, "3 %d %d %d\n", elements[i]->nodes[0].index, elements[i]->nodes[1].index, elements[i]->nodes[2].index);
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
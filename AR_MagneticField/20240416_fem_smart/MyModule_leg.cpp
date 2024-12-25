#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include "magic.h"
#include "const.h"
#include "GenMatrix.h"
#include "MyFunction.h"
namespace py = pybind11;
using namespace std;
//********************************************************************************
class Carray{
private:
    int n; // 磁石の数
    // double* theta_of_magnets; // 磁化方向の配列
    vector<double> theta_of_magnets;
public:
    Carray(double theta, int n, py::array_t<double> np_array, int width, int height, bool MAKE_DATASET, int data_index, const std::string& DATA_PASS, const std::string& MESH_PASS);
    ~Carray();
    int len();
    double sum();
    int read_mesh();
    int read_akimacsv();
    int fem();
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
    ELEMENT *elements; 
    double *A;
    double *del_A;
    double *r;         
    double *b;         
    NODE *nodes;       
    int num_OnBoundNodes;
    
    int max_nrloop = 100;
    double M = 1.25;
    double ICCG_conv = 1E-6;
    double NR_conv = 1E-3;
    double ICCG_conv_MAX;
    bool printinfo = false;

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
    // std::cout << DATA_PASS << std::endl;
    // theta_of_magnets = new double[n]();
    theta_of_magnets.resize(n);
    for(int i = 0; i < n; i++) {
        double tmp = np_array.mutable_at(i);
        this->theta_of_magnets[i] = tmp;
    }
}
//********************************************************************************
Carray::~Carray() {
    
} 
//********************************************************************************
int Carray::read_mesh(){
    using namespace std;
    // cout << "This is read_mesh method" << endl;

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
    nodes = new NODE[num_nodes];
    for (int i=0; i<num_nodes; i++) // nodeデータ格納
    {   
        int node_index;
        double node_x, node_y, node_z;
        // [0]->節点番号, [1]->x座標, [2]->y座標, [3]->z座標
        if (fscanf(fp, "%d %lf %lf %lf", &node_index, &node_x, &node_y, &node_z) == 0){}
        nodes[i].setNode(node_index-1, node_x, node_y, node_z);
        /*境界条件判定*/
        if(1==0){ 
            nodes[i].isOnBoundary = true;
            nodes[i].Az = 0.0;
        }else{
            nodes[i].isOnBoundary = false;
        }
        // cout << i << "   " << nodes[i].index << "   " << nodes[i].x_coordinate << "   " << nodes[i].y_coordinate << "   " << nodes[i].y_coordinate << endl;
    }
    if (fgets(buf, sizeof(buf), fp) == 0); // ２行読み飛ばし?
    if (fgets(buf, sizeof(buf), fp) == 0); 
    if (fgets(buf, sizeof(buf), fp) == 0);
    if (fscanf(fp, "%d%*c", &num_elements) == 0){} // エレメント数読み込み  
    elements = new ELEMENT[num_elements];
    for (int i = 0; i < num_elements; i++) {
        int elemNum, matNum, node_idx1, node_idx2, node_idx3;
        // 要素番号、材料番号、節点番号を読み込む
        fscanf(fp, "%d %*d %d %*d %*d %d %d %d", &elemNum, &matNum, &node_idx1, &node_idx2, &node_idx3);
        // ELEMENTオブジェクトにデータを設定
        elements[i].setElement(elemNum, matNum, nodes[node_idx1-1], nodes[node_idx2-1], nodes[node_idx3-1]); // 節点番号を0ベースで調整
    }
    // elements[num_elements-1].print_info();
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
    r = new double[num_nonboundnodes];
    b = new double[num_nodes]; // 最終的な解ベクトル
    del_A = new double [num_nodes];
    A = new double[num_nodes];
    ICCG_conv_MAX = num_nodes;

    for (int n=0; n<num_nodes; n++){
        A[n] = 0.0;
    }
    
    // cとかdとかdeltaとか色々計算
    for (int e=0; e<num_elements; e++){
        elements[e].cal_constants();
    }

    auto start1 = chrono::steady_clock::now();
    auto start2 = start1;

    // NR法
    for (int nrloop=0; nrloop<max_nrloop; nrloop++){
        // 配列初期化
        for (int n=0; n<num_nodes; n++){
            b[n] = 0.0;
            del_A[n] = 0.0;
        }
        // K, b 作成
        for (int e=0; e<num_elements; e++){
            int mat = elements[e].materialNumber;
            // U更新
            elements[e].cal_U();
            for (int i=0; i<3; i++){
                for (int j=0; j<3; j++){
                    K.add(elements[e].nodes[i].index, elements[e].nodes[j].index, elements[e].nyu * elements[e].S[i][j] + 2.0 * elements[e].dnyu_dB2 * elements[e].U[i] * elements[e].U[j] / elements[e].delta);
                }
                b[elements[e].nodes[i].index] += (elements[e].J_0 * elements[e].delta / 3.0) - (elements[e].nyu * elements[e].U[i]);
                if(4 <= mat && mat <= 10){
                    // 永久磁石なら +=Jm
                    b[elements[e].nodes[i].index] += elements[e].nyu * (M * cos(theta_of_magnets[mat-4]) * elements[e].d[i] - M * sin(theta_of_magnets[mat-4]) * elements[e].c[i]) / 2.0;

                }
            }
            // cout << "nyu[e]: " << elements[e].nyu << endl;
        }

        // 固定境界の行番号を格納する
        // v1.push_back(番号);で入れていく
        vector<int> v1;

        // 境界条件の付与
        K.setBoundaryCondition(b, v1, 0.0);
        // 不完全コレスキー分解
        GenMatrix::icdcmp(K, L, 1.05);
        // ICCG本体
        GenMatrix::iccgSolv(K, L, b, del_A, r, ICCG_conv, ICCG_conv_MAX);
        
        // A更新
        for (int n=0; n<num_nodes; n++){
            A[n] += del_A[n];
            nodes[n].add_A(del_A[n]);
        }
        for (int e=0; e<num_elements; e++){
            int index1 = elements[e].nodes[0].index;
            int index2 = elements[e].nodes[1].index;
            int index3 = elements[e].nodes[2].index;
            elements[e].update_node(nodes[index1], nodes[index2], nodes[index3]);
        }
        // 時間計測
        auto end = chrono::steady_clock::now();

        // 非線形性考慮しない場合はNRloop抜ける
        if (nonliner == false){
            for (int e=0; e<num_elements; e++){
                elements[e].update_B();
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
                    elements[e].update_B();
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
                    elements[e].update_B();
                    elements[e].update_nyu(akima_data);
                    // elements[e].print_info();
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
    // 深層学習用出力側データ作成関数(20240328現在動作未確認)
    if(this->MAKE_DATASET ==1){
        double x[3];
        double y[3];
        // ピクセルごとのBx, Byを求めcsvへ
        FILE *fp;
        const char *fout = (this->DATA_PASS).c_str();
        fp = fopen(fout, "w");
        if (fp == NULL){
            cout << fout << "ファイルが開けませんでした。" << endl;
            return -1;
        }
        for(int j=0; j<this->height; j++)
        {
            for(int k=0; k<this->width; k++)
            {   
                double px = k + 0.5; // ピクセル中心のx座標
                double py = this->height - ( j + 0.5 );
                // このピクセルがどの要素内にあるかを探索する
                for(int e=0; e<num_elements; e++)
                {
                    for (int i = 0; i < 3; i++)
                    {
                        x[i] = elements[e].nodes[i].x;
                        y[i] = elements[e].nodes[i].y;
                        // メモ：この要素は(x1, y1), (x2, y2), (x3, y3)に囲まれてる
                    }
                    double abXat = (x[1]-x[0])*(py-y[0])-(y[1]-y[0])*(px-x[0]);
                    double bcXbt = (x[2]-x[1])*(py-y[1])-(y[2]-y[1])*(px-x[1]);
                    double caXct = (x[0]-x[2])*(py-y[2])-(y[0]-y[2])*(px-x[2]);
                    if( abXat == 0.0 || bcXbt == 0.0 || caXct == 0.0){
                        // cout << "辺と重なっています  ピクセル（" << k + 1 << "," << j + 1 << ")"<< endl;
                        // cout << abXat << "  "  << bcXbt << "  " << caXct << endl;
                        // どの要素かわかったらBx[e]格納していく
                        fprintf(fp, "%.16lf,%.16lf,", elements[e].Bx, elements[e].By);
                        break;
                    }else if(( abXat > 0.0 && bcXbt > 0.0 && caXct > 0.0) || ( abXat < 0.0 && bcXbt < 0.0 && caXct < 0.0)){
                        // cout << "ピクセル（" << k + 1 << "," << j + 1 << "）は、要素番号" << e << "内に存在し,Bx=" << Bx[e] << " By=" << By[e] << endl;
                        // どの要素かわかったらBx,By書き込んでいく
                        // fprintf(fp, "%.16lf,%.16lf,", Bx[e], By[e]);
                        fprintf(fp, "%.16lf,%.16lf,", elements[e].Bx, elements[e].By);
                        break;
                    }
                    if(e==num_elements-1){
                        cout << "error! : ピクセル（" << k + 1 << "," << j + 1 << "）は、どの要素にも属しません" << endl;
                        return 1;
                    }
                }
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
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
        write_vtk_contour_2(name_vtk_con, nodes, elements, num_nodes, num_elements, A);
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
int Carray::fem(){
    using namespace std;
    // cout << this->theta << endl;
    int aaa[10];
    int num_nodes;
    int num_elements;
    char ch;
    char buf[256];
    int num_boundnodes = 0;

    // 永久磁石パラメータ**************************
    double M = 1.25;
    // ********************************************

    // msh1ファイル読み取り*******************************************************************
    FILE *fppp;
    const char *fmesh = (this->MESH_PASS).c_str();
    fppp = fopen(fmesh, "r");
    if (fppp == NULL)
    {
        printf("File open error.\n");
        return 1;
    }
    // fscanfを使用して、ファイルから1行読み飛ばすには、"%*[^\n]%*c"という書式指定子を使用します。これは、改行文字（\n）が現れるまで、すべての文字を読み取り、そして改行文字自体も読み取り、読み取った文字を破棄することを示します。
    // この例では、ファイルの最初の行を読み飛ばし、2行目を読み取り、整数と改行文字を出力しています。%*[^\n]%*cを使用して、最初の行を読み飛ばしています。注意してほしいのは、読み飛ばしの後、次の行を読み込む前に、1文字だけ読み取る必要があるため、%d%cを使用して整数と改行文字を読み込んでいます。

    // 最初の行を読み飛ばす
    if (fscanf(fppp, "%*[^\n]%*c") == 0)
    {
    }

    // ノード数を読み込む
    if (fscanf(fppp, "%d%c", &num_nodes, &ch) == 0)
    {
    }

    // node配列確保
    double **node = new double *[num_nodes]; // -1
    for (int i = 0; i < num_nodes; i++)
    {
        node[i] = new double[4];
    }

    // nodeデータ格納
    for (int i = 0; i < num_nodes; i++)
    {
        if (fscanf(fppp, "%lf %lf %lf %lf", &node[i][0], &node[i][1], &node[i][2], &node[i][3]) == 0)
        {
        }
        // 確認用
        // cout << node[i][1] << " " << node[i][2] << endl;
    }

    // 二行読み飛ばし
    if (fgets(buf, sizeof(buf), fppp) == 0)
        ;
    if (fgets(buf, sizeof(buf), fppp) == 0)
        ;
    if (fgets(buf, sizeof(buf), fppp) == 0)
        ;

    // エレメント数を読み込む
    if (fscanf(fppp, "%d%c", &num_elements, &ch) == 0)
    {
    }

    // element配列確保
    int **element = new int *[num_elements];//-2
    for (int i = 0; i < num_elements; i++)
    {
        element[i] = new int[5];
    }

    // elementデータ格納
    for (int i = 0; i < num_elements; i++)
    {
        if (fscanf(fppp, "%d %d %d %d %d %d %d %d", &element[i][0], &aaa[0], &element[i][1], &aaa[1], &aaa[2], &element[i][2], &element[i][3], &element[i][4]) == 0)
        {
        }
        // 確認用
        element[i][2] -= 1;
        element[i][3] -= 1;
        element[i][4] -= 1;
        // cout << element[i][0] << " " << element[i][1] << " " << element[i][2] << " " << element[i][3] << " " << element[i][4] << endl;
    }
    fclose(fppp);
    // *****************************************************************************************
    // akimaデータ読み込み
    double Data[100][100];
    const char *fl = "akima4.csv"; // 間違ってるから使うな
    FILE *fpp;
    fpp = fopen(fl, "r");
    if (fpp == NULL)
    {
        printf("error\n");
    }
    else
    {
        double tmp;
        for (int i = 0; i < 40; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                if (fscanf(fpp, "%lf,", &Data[i][j]) == 1)
                {
                }
                // double data = tmp;
                // printf("data[%d][%d] = %f", i+1, j+1, Deta[i][j]);
            }
            // printf("\n");
        }
    }
    fclose(fpp);
    // ***************************************************

    // const char* fn = "test_element.csv";
    // const char* fm = "test_node.csv";

    // int **element = new int*[num_elements];//10
    // for (int i = 0; i < num_elements; i++) {
    //     element[i] = new int[5];
    // }

    // double **node = new double*[num_nodes];//11
    // for (int i = 0; i < num_nodes; i++) {
    //     node[i] = new double[4];
    // }

    // // エレメント格納**************************************************************

    // // // ファイルポインタfp
    // // FILE* fp;
    // // fp = fopen(fn, "r");
    // // if(fp == NULL) {
    // //   printf("error\n");
    // // }
    // // else {
    // //   int tmp;
    // //   for(int i = 0; i < num_elements; i++) {
    // //     for(int j = 0; j < 5; j++) {
    // //       if(fscanf(fp, "%d,", &tmp)==1){
    // //         element[i][j] = tmp;
    // //         // printf("data[%d][%d] = %d", i+ 1, j+1, element[i][j]);
    // //       }else{
    // //         cout << "要素の読み込みに失敗しました "<<endl;
    // //         return 1;
    // //       };

    // //     }
    // //     // printf("\n");
    // //   }
    // // }
    // // fclose(fp);

    // // ノード格納*******************************************************************
    // // fp = fopen(fm, "r");
    // // if(fp == NULL) {
    // //   printf("error\n");
    // // }
    // // else {
    // //   double tmp;
    // //   for(int i = 0; i < num_nodes; i++) {
    // //     for(int j = 0; j < 4; j++) {
    // //       if(fscanf(fp, "%lf,", &tmp)==1){
    // //         node[i][j] = tmp;
    // //         // printf("data[%d][%d] = %lf", i+ 1, j+1, node[i][j]);
    // //       }else{
    // //         return 1;
    // //       };
    // //     }
    // //     // printf("\n");
    // //   }
    // // }
    // // fclose(fp);

    // 変数定義***********************************************************************
    int num_nonboundnodes = num_nodes - num_boundnodes;
    int *vertex_node = new int[3]; // 1
    double *x = new double[3];     // 2
    double *y = new double[3];     // 3
    double *b = new double[3];     // 4
    double *c = new double[3];     // 5
    double *d = new double[3];     // 6
    // double *bb = new double[num_nodes]; // 7

    // double **K = new double *[num_nodes]; // 8
    // for (int i = 0; i < num_nodes; i++)
    // {
    //     K[i] = new double[num_nodes];
    // }

    int *known_then_1 = new int[num_nodes]; // 9
    double ep = 1.0e-5;
    int known_counter;
    double *a = new double[num_nodes]; // 13
    for (int i = 0; i < num_nodes; i++)
    {
        a[i] = 0.0;
    }
    double **xyzBxByBz = new double *[num_elements]; // 14
    for (int i = 0; i < num_elements; i++)
    {
        xyzBxByBz[i] = new double[6];
    }

    double ***S = new double **[3]; // 15
    for (int i = 0; i < 3; i++)
    {
        S[i] = new double *[3];
        for (int j = 0; j < 3; j++)
        {
            S[i][j] = new double[num_elements];
        }
    }
    // double S[3][3][num_elements];

    double U[3];
    // double *A = new double[num_nodes];           // 16
    double *nyu = new double[num_elements];      // 17
    double *dnyu_dB2 = new double[num_elements]; // 18
    double *delta = new double[num_elements];    // 19
    double *J_0 = new double[num_elements];      // 20
    double *delta_A = new double[num_nodes];     // 21
    double *B2 = new double[num_elements];       // 22
    double norm_delta_A;
    double *Bx = new double[num_elements];      // 23
    double *By = new double[num_elements];      // 24
    double *Jm = new double[num_elements];      // 25
    double *dJm_dAj = new double[num_elements]; // 26
    double *norm_B = new double[num_elements];  // 27
    double B22 = 0.0;
    double a0 = 0.0;
    double a1 = 0.0;
    double a2 = 0.0;
    double a3 = 0.0;
    double del;
    // for (int i = 0; i < num_nodes; i++)
    // {
    //     for (int j = 0; j < num_nodes; j++)
    //     {
    //         // K[i][j] = 0.0;
    //     }
    //     bb[i] = 0.0;
    //     A[i] = 0.0;
    // }
    GenMatrix K(num_nodes, num_nodes); //29
    // GenMatrix C(num_nodes, num_nonboundnodes);
    // GenMatrix Ct(num_nonboundnodes, num_nodes);
    // GenMatrix CtK(num_nonboundnodes, num_nodes);
    // GenMatrix CtKC(num_nonboundnodes, num_nonboundnodes);

    GenMatrix L(num_nonboundnodes, num_nonboundnodes); //30

    double *A = new double[num_nodes];//16
    // double *A = new double[num_nonboundnodes];
    double *r = new double[num_nonboundnodes]; // 28
    // double *Ctb = new double[num_nonboundnodes];

    double *bb = new double[num_nodes]; // 7 最終的な右辺ベクトル
    // double *CA = new double[num_nodes]; // 最終的な解ベクトル
    int mat;
    bool liner_or_nonliner = true;
    // double **pixel_B = new double *[this->width]; // 31
    // for (int i = 0; i < this->width; i++)
    // {
    //     pixel_B[i] = new double[this->width];
    // }

    // *****************************************************************************

    // cout << this->theta_of_magnets[0] << endl;

    // 固定境界上かどうかの判定
    known_counter = 0;
    for (int i = 0; i < num_nodes; i++)
    {
        // x=0またはy=1またはx=1
        if (0 == 1)
        {
            known_then_1[i] = 1;
            known_counter++;
            A[i] = 0.0;
        }
        else
        {
            // それ以外の点
            known_then_1[i] = 0;
        }
        // cout << "known_then_1[" << i << "] = " << known_then_1[i] << endl;
    }
    // cout << "固定境界上節点数 : " << known_counter << endl;

    // 要素ごとにSとデルタとJ_0の計算 **********************************************************
    for (int e = 0; e < num_elements; e++)
    {

        for (int i = 0; i < 3; i++)
        {
            vertex_node[i] = element[e][i + 2];
            x[i] = node[vertex_node[i]][1];
            y[i] = node[vertex_node[i]][2];
        }
        for (int i = 0; i < 3; i++)
        {
            // b[i] = x[(1+i)%3]*y[(2+i)%3]-x[(2+i)%3]*y[(1+i)%3];
            c[i] = y[(1 + i) % 3] - y[(2 + i) % 3];
            d[i] = x[(2 + i) % 3] - x[(1 + i) % 3];
        }

        // デルタ
        delta[e] = (x[0] * y[1] + x[1] * y[2] + x[2] * y[0] - y[0] * x[1] - y[1] * x[2] - y[2] * x[0]) / 2.0;
        // cout << e << ".  " << delta[e] << endl;

        // S（3×3×要素数）
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                S[i][j][e] = (c[i] * c[j] + d[i] * d[j]) / (4.0 * delta[e]);
            }
        }

        // 材料判定
        mat = element[e][1];
        if (mat == 1)
        {
            J_0[e] = 0.0;
        }
        else if (mat == 2)
        {
            J_0[e] = 0.0;
        }
        else if (mat == 3)
        {
            // コイルなら
            J_0[e] = 0.0;
            // J_0[e] = 1.0/S_coil;
        }
        else if (mat == 4)
        {
            // 永久磁石なら
            J_0[e] = 0.0;
        }
        else if (mat == 5)
        {
            // 永久磁石なら
            J_0[e] = 0.0;
        }
        else if (mat == 6)
        {
            // 永久磁石なら
            J_0[e] = 0.0;
        }
        else if (mat == 7)
        {
            // 永久磁石なら
            J_0[e] = 0.0;
        }
        else if (mat == 8)
        {
            // 永久磁石なら
            J_0[e] = 0.0;
        }
        else if (mat == 9)
        {
            // 永久磁石なら
            J_0[e] = 0.0;
        }
        else if (mat == 10)
        {
            // 永久磁石なら
            J_0[e] = 0.0;
        }
        else if (mat == 11)
        {
            // 永久磁石なら
            J_0[e] = 0.0;
        }
        else if (mat == 12)
        {
            // 永久磁石なら
            J_0[e] = 0.0;
        }
        else if (mat == 13)
        {
            // 永久磁石なら
            J_0[e] = 0.0;
        }
        else
        {
            cout << "定義されていない材料です" << endl;
            return 1;
        }
    }

    // A,nyu,dnyu_dB2の初期値設定***************************************************
    for (int n = 0; n < num_nodes; n++)
    {
        A[n] = 0.0;
    }
    for (int e = 0; e < num_elements; e++)
    {
        // 材料判定
        mat = element[e][1];
        if (mat == 1)
        { // 空気
            nyu[e] = 1.0 / (myu0);
            // nyu[e] = 1.0/(myu0*myur);
            dnyu_dB2[e] = 0.0;
            // cout << "空気！" << endl;
        }
        else if (mat == 2)
        { // 鉄
            nyu[e] = 1.0 / (myu0 * myur);
            dnyu_dB2[e] = 0.0;
            // cout << "鉄！" << endl;
        }
        else if (mat == 3)
        { // コイル
            nyu[e] = 1.0 / (myu0);
            dnyu_dB2[e] = 0.0;
        }
        else if (mat == 4)
        { // 永久磁石
            nyu[e] = 1.0 / (myu0 * 1.05);
            dnyu_dB2[e] = 0.0;
        }
        else if (mat == 5)
        { // 永久磁石
            nyu[e] = 1.0 / (myu0 * 1.05);
            dnyu_dB2[e] = 0.0;
        }
        else if (mat == 6)
        { // 永久磁石
            nyu[e] = 1.0 / (myu0 * 1.05);
            dnyu_dB2[e] = 0.0;
        }
        else if (mat == 7)
        { // 永久磁石
            nyu[e] = 1.0 / (myu0 * 1.05);
            dnyu_dB2[e] = 0.0;
        }
        else if (mat == 8)
        { // 永久磁石
            nyu[e] = 1.0 / (myu0 * 1.05);
            dnyu_dB2[e] = 0.0;
        }
        else if (mat == 8)
        { // 永久磁石
            nyu[e] = 1.0 / (myu0 * 1.05);
            dnyu_dB2[e] = 0.0;
        }
        else if (mat == 9)
        { // 永久磁石
            nyu[e] = 1.0 / (myu0 * 1.05);
            dnyu_dB2[e] = 0.0;
        }
        else if (mat == 10)
        { // 永久磁石
            nyu[e] = 1.0 / (myu0 * 1.05);
            dnyu_dB2[e] = 0.0;
        }
        else if (mat == 11)
        { // 永久磁石
            nyu[e] = 1.0 / (myu0 * 1.05);
            dnyu_dB2[e] = 0.0;
        }
        else if (mat == 12)
        { // 永久磁石
            nyu[e] = 1.0 / (myu0 * 1.05);
            dnyu_dB2[e] = 0.0;
        }
        else if (mat == 13)
        { // 永久磁石
            nyu[e] = 1.0 / (myu0 * 1.05);
            dnyu_dB2[e] = 0.0;
        }
        else
        {
            cout << "定義されていない材料です" << endl;
            return 1;
        }
    }

    // NR法***************************************************************
    for (int nrloop = 0; nrloop < 10; nrloop++)
    {
        auto start = chrono::steady_clock::now();
        // 初期化
        for (int i = 0; i < num_nodes; i++)
        {
            for (int j = 0; j < num_nodes; j++)
            {
                // K[i][j] = 0.0;
            }
            bb[i] = 0.0;
            delta_A[i] = 0.0;
        }
        // K,bb作成
        for (int e = 0; e < num_elements; e++)
        {
            // その要素の節点番号3つ調べる
            for (int i = 0; i < 3; i++)
            {
                vertex_node[i] = element[e][i + 2];
                x[i] = node[vertex_node[i]][1];
                y[i] = node[vertex_node[i]][2];
            }
            // cとdもとめる（3つずつ）
            for (int i = 0; i < 3; i++)
            {
                // b[i] = x[(1+i)%3]*y[(2+i)%3]-x[(2+i)%3]*y[(1+i)%3];
                c[i] = y[(1 + i) % 3] - y[(2 + i) % 3];
                d[i] = x[(2 + i) % 3] - x[(1 + i) % 3];
            }

            // Uの計算
            for (int i = 0; i < 3; i++)
            {
                U[i] = S[i][0][e] * A[vertex_node[0]] + S[i][1][e] * A[vertex_node[1]] + S[i][2][e] * A[vertex_node[2]];
            }


            // 全体Kに9か所bに3か所足しこむ
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {   
                    double tmp;
                    tmp = nyu[e] * S[i][j][e] + 2.0 * dnyu_dB2[e] * U[i] * U[j] / delta[e];
                    // cout << "Kval: " << tmp << endl;
                    K.add(vertex_node[i], vertex_node[j], nyu[e] * S[i][j][e] + 2.0 * dnyu_dB2[e] * U[i] * U[j] / delta[e]);
                    // K[vertex_node[i]][vertex_node[j]] += nyu[e] * S[i][j][e] + 2.0 * dnyu_dB2[e] * U[i] * U[j] / delta[e];
                    // cout << K[vertex_node[0]][vertex_node[2]] << endl;
                }
                bb[vertex_node[i]] += (J_0[e] * delta[e] / 3.0) - (nyu[e] * U[i]);
                // matが永久磁石だったらJmたす***************************************************
                mat = element[e][1];
                if (mat == 4)
                {
                    bb[vertex_node[i]] += nyu[e] * (M * cos(this->theta_of_magnets[0]) * d[i] - M * sin(this->theta_of_magnets[0]) * c[i]) / 2.0;
                    // cout << "M * cos(this->theta_of_magnets[0]) * d[i]: " << M * cos(this->theta_of_magnets[0]) * d[i] << endl;
                }
                else if (mat == 5)
                {
                    bb[vertex_node[i]] += nyu[e] * (M * cos(this->theta_of_magnets[1]) * d[i] - M * sin(this->theta_of_magnets[1]) * c[i]) / 2.0;
                }
                else if (mat == 6)
                {
                    bb[vertex_node[i]] += nyu[e] * (M * cos(this->theta_of_magnets[2]) * d[i] - M * sin(this->theta_of_magnets[2]) * c[i]) / 2.0;
                }
                // *************************************************************************
            }
        }

        // for (int n=0; n<num_nodes; n++){
        //     cout << "b[" << n << "]: " << bb[n] << endl;
        // }

        double norm = 0.0;
        for (int i = 0; i < num_nodes; i++)
        {
            // norm += Ctb[i] * Ctb[i];
        }
        double Bnorm = sqrt(norm);

        // 収束条件
        double ICCG_conv = 1E-6;

        // 最大試行回数
        double ICCG_conv_MAX = num_nodes;

        // boundノードの番号を格納するvector作成
        vector<int> v1;
        // 固定境界の行番号を格納する
        // v1.push_back(番号);で入れてく

        // 境界条件の付与
        K.setBoundaryCondition(bb, v1, 0.0);

        // 不完全コレスキー分解
        GenMatrix::icdcmp(K, L, 1.05);
        // ICCG本体
        GenMatrix::iccgSolv(K, L, bb, delta_A, r, ICCG_conv, ICCG_conv_MAX);

        // A更新
        for (int i = 0; i < num_nodes; i++)
        {
            if (known_then_1[i] != 1)
            {
                A[i] += delta_A[i];
            }
        }

        // 線形でいいときはこれで抜ける
        if(liner_or_nonliner==true){
            break;
        }else{
            // delta_Aの大きさの計算
            norm_delta_A = 0.0;
            for (int i = 0; i < num_nodes; i++)
            {
                norm_delta_A += delta_A[i] * delta_A[i];
            }
            norm_delta_A = sqrt(norm_delta_A);
        }

        // 収束したらNRloop抜ける
        if (norm_delta_A < 1.0e-3)
        {
            cout << "NR loop: " << nrloop + 1 <<  "   delta_Aのノルムは　" << norm_delta_A << "　よってNR法を終了します。\n";
            break;
        }
        else
        {
            // 収束しなければすべての要素のB^2,nyu,dnyu/dB^2計算
            // cout << "delta_Aのノルム: " << norm_delta_A << "  よって" << nrloop + 2 << "度目の計算に入ります。\n";
            

            auto end = chrono::steady_clock::now();
            // cout << "NR法の1ループにかかった時間: "
            //      << chrono::duration_cast<chrono::milliseconds>(end - start).count()
            //      << " ミリ秒" << endl;

            for (int e = 0; e < num_elements; e++)
            {
                // Bx
                Bx[e] = (d[0] * A[vertex_node[0]] + d[1] * A[vertex_node[1]] + d[2] * A[vertex_node[2]]) / (2.0 * delta[e]);

                // By
                By[e] = -1.0 * (c[0] * A[vertex_node[0]] + c[1] * A[vertex_node[1]] + c[2] * A[vertex_node[2]]) / (2.0 * delta[e]);

                // B^2
                B2[e] = Bx[e] * Bx[e] + By[e] * By[e];

                // 材料判定したあと　Bみてニュー変える
                if (element[e][1] == 2)
                { // 鉄なら
                    for (int i = 0; i < 40; i++)
                    {
                        if (Data[i][0] <= B2[e] && B2[e] < Data[i + 1][0])
                        {
                            B22 = Data[i][0];
                            a0 = Data[i][1];
                            a1 = Data[i][2];
                            a2 = Data[i][3];
                            a3 = Data[i][4];
                            break;
                        }
                    }

                    del = B2[e] - B22;

                    nyu[e] = a0 + a1 * del + a2 * del * del + a3 * del * del * del;
                    dnyu_dB2[e] = a1 + 2 * a2 * del + 3 * a3 * del * del;
                }
            }
        }
    }

    for (int i = 0; i < num_nodes; i++)
    {
        // cout << "A[" << i << "] = " << A[i] << endl;
    }

    // Bx,Byの計算*************************************************************
    // cout << "                     Bx               By" << endl;
    for (int e = 0; e < num_elements; e++)
    {
        for (int i = 0; i < 3; i++)
        {
            vertex_node[i] = element[e][i + 2];
            x[i] = node[vertex_node[i]][1];
            y[i] = node[vertex_node[i]][2];
        }
        for (int i = 0; i < 3; i++)
        {
            c[i] = y[(1 + i) % 3] - y[(2 + i) % 3];
            d[i] = x[(2 + i) % 3] - x[(1 + i) % 3];
        }
        // delta[e] = (x[0]*y[1]+x[1]*y[2]+x[2]*y[0]-y[0]*x[1]-y[1]*x[2]-y[2]*x[0])/2.0;

        // x
        xyzBxByBz[e][0] = (x[0] + x[1] + x[2]) / 3.0;
        // y
        xyzBxByBz[e][1] = (y[0] + y[1] + y[2]) / 3.0;
        // z
        xyzBxByBz[e][2] = 0.0;
        // Bx
        xyzBxByBz[e][3] = (d[0] * A[vertex_node[0]] + d[1] * A[vertex_node[1]] + d[2] * A[vertex_node[2]]) / (2.0 * delta[e]);
        Bx[e] = xyzBxByBz[e][3];
        // By
        xyzBxByBz[e][4] = -1.0 * (c[0] * A[vertex_node[0]] + c[1] * A[vertex_node[1]] + c[2] * A[vertex_node[2]]) / (2.0 * delta[e]);
        By[e] = xyzBxByBz[e][4];
        // Bz
        xyzBxByBz[e][5] = 0.0;
        // |B|
        norm_B[e] = sqrt(xyzBxByBz[e][3] * xyzBxByBz[e][3] + xyzBxByBz[e][4] * xyzBxByBz[e][4]);

        // cout << "要素 " << setw(4) << e + 1 << " : " << setw(15) << xyzBxByBz[e][3] << "  " << setw(15) << xyzBxByBz[e][4] << endl;
    }

    // BxBy.csvへの吐き出し(グラフR用）****************************************************************************
    // const char *fout = "BxBy.csv";
    // const char *fout = "vectors.csv";
    // const char *name_vtk_vec = "vectors.vtk";
    // const char *name_vtk_mat = "material.vtk";
    // const char *name_vtk_con = "contour.vtk";

    if(this->MAKE_DATASET ==1){
        // ピクセルごとのBx, Byを求めcsvへ
        FILE *fo;
        const char *fout = (this->DATA_PASS).c_str();
        fo = fopen(fout, "w");
        if (fo == NULL){
            cout << fout << "ファイルが開けませんでした。" << endl;
            return -1;
        }
        for(int j=0; j<this->height; j++)
        {
            for(int k=0; k<this->width; k++)
            {   
                double px = k + 0.5;
                double py = this->height - ( j + 0.5 );
                // このピクセルがどの要素内にあるかを探索する
                for(int e=0; e<num_elements; e++)
                {
                    for (int i = 0; i < 3; i++)
                    {
                        vertex_node[i] = element[e][i + 2];
                        x[i] = node[vertex_node[i]][1];
                        y[i] = node[vertex_node[i]][2];
                        // メモ：この要素は(x1, y1), (x2, y2), (x3, y3)に囲まれてる
                    }
                    double abXat = (x[1]-x[0])*(py-y[0])-(y[1]-y[0])*(px-x[0]);
                    double bcXbt = (x[2]-x[1])*(py-y[1])-(y[2]-y[1])*(px-x[1]);
                    double caXct = (x[0]-x[2])*(py-y[2])-(y[0]-y[2])*(px-x[2]);
                    if( abXat == 0.0 || bcXbt == 0.0 || caXct == 0.0){
                        // cout << "辺と重なっています  ピクセル（" << k + 1 << "," << j + 1 << ")"<< endl;
                        // cout << abXat << "  "  << bcXbt << "  " << caXct << endl;
                        // どの要素かわかったらBx[e]格納していく
                        fprintf(fo, "%.16lf,%.16lf,", Bx[e], By[e]);
                        break;
                    }else if(( abXat > 0.0 && bcXbt > 0.0 && caXct > 0.0) || ( abXat < 0.0 && bcXbt < 0.0 && caXct < 0.0)){
                        // cout << "ピクセル（" << k + 1 << "," << j + 1 << "）は、要素番号" << e << "内に存在し,Bx=" << Bx[e] << " By=" << By[e] << endl;
                        // どの要素かわかったらBx[e]格納していく
                        fprintf(fo, "%.16lf,%.16lf,", Bx[e], By[e]);
                        break;
                    }
                    if(e==num_elements-1){
                        cout << "error! : ピクセル（" << k + 1 << "," << j + 1 << "）は、どの要素にも属しません" << endl;
                        return 1;
                    }
                }
            }
            fprintf(fo, "\n");
        }
        fclose(fo);
    }






    
    

    // FILE *fo;
    // fo = fopen(fout, "w");

    // if (fo == NULL)
    // {
    //     cout << fout << "ファイルが開けませんでした。" << endl;
    //     return -1;
    // }

    // if (fout=="BxBy.csv"){
    //     // 最初の三行（グラフR用）
    //     fprintf(fo, "DataFormat,5,,,,\n");
    //     fprintf(fo, ",,,,,\n");
    //     fprintf(fo, ",,,,,\n");
    // }

    // // 行：要素　列：x,y,z,Bx,By,Bz
    // for (int i = 0; i < num_elements; i++)
    // {
    //     for (int j = 0; j < 6; j++)
    //     {
    //         double tmp;
    //         tmp = xyzBxByBz[i][j];
    //         fprintf(fo, "%.16lf,", tmp);
    //     }
    //     fprintf(fo, "\n");
    // }
    // fclose(fo);

    write_vtk_material(name_vtk_mat, node, element, num_nodes, num_elements);
    write_vtk_vector(name_vtk_vec, node, element, Bx, By, num_nodes, num_elements);
    write_vtk_contour(name_vtk_con, node, element, Bx, By, num_nodes, num_elements, A);
    // write_vtk_material(name_vtk_mat, node, element, )
    // ************************************************************************
    delete[] vertex_node; // 1
    delete[] x;           // 2
    delete[] y;           // 3
    delete[] b;           // 4
    delete[] c;           // 5
    delete[] d;           // 6
    // for (int i = 0; i < num_nodes; i++)
    // { // 8
    //     delete[] K[i];
    // }
    // delete[] K;

    delete[] known_then_1; // 9

    delete[] a; // 13

    for (int i = 0; i < num_elements; i++)
    { // 14
        delete[] xyzBxByBz[i];
    }
    delete[] xyzBxByBz;

    for (int i = 0; i < 3; i++)
    { // 15
        for (int j = 0; j < 3; j++)
        {
            delete[] S[i][j];
        }
        delete[] S[i];
    }
    delete[] S;

    delete[] A;        // 16
    delete[] nyu;      // 17
    delete[] dnyu_dB2; // 18
    delete[] delta;    // 19
    delete[] J_0;      // 20
    delete[] delta_A;  // 21
    delete[] B2;       // 22
    delete[] Bx;       // 23
    delete[] By;       // 24
    delete[] Jm;       // 25
    delete[] dJm_dAj;  // 26

    delete[] bb;          // 7
    
    for (int i = 0; i < num_nodes; i++)
    { // -1
        delete[] node[i];
    }
    delete[] node;

    for (int i = 0; i < num_elements; i++)
    { // -2
        delete[] element[i];
    }
    delete[] element;

    delete[] norm_B; // 27
    delete[] r;      // 28



    return 0;
}

PYBIND11_MODULE(Carray, m) {
    m.doc() = "This is a test module";
    py::class_<Carray>(m, "Carray")
		.def(py::init<double, int,py::array_t<double>, int, int, bool, int, const std::string&, const std::string&>())
		.def("len", &Carray::len)
        .def("sum", &Carray::sum)
        .def("read_mesh", &Carray::read_mesh)
        .def("read_akimacsv", &Carray::read_akimacsv)
        .def("analysis", &Carray::analysis)
        .def("write_vtk", &Carray::write_vtk)
        .def("fem", &Carray::fem);
}
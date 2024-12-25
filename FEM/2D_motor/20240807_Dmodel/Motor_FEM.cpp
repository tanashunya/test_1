#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include "magic.h"
#include "const.h"
#include "GenMatrix.h"
#include "MyFunction.h"
namespace py = pybind11;

using namespace std;

// D model解析用クラス
class MOTOR_FEM{
    private:
    public:
        MOTOR_FEM();
        ~MOTOR_FEM();
        void main_func();
        void count_lines();
        void print_mesh_info();
        int read_akima();
        int preprocess_periodic();
        void sort_nodes();
        void rotation(double angle); // deg
        int preprocess_for_make_C();
        int make_C();
        int NewtonRaphson_loop();
        int write_vtk();
        
        // メッシュ関連
        int read_rot_elementCSV();
        int read_rot_node_CSV();
        int read_sta_element_CSV();
        int read_sta_node_CSV();
        const char *name_rot_elementCSV = "mesh_Dmodel/new/RotElementData.csv";
        const char *name_rot_nodeCSV    = "mesh_Dmodel/new/RotNodeData.csv";
        const char *name_sta_elementCSV = "mesh_Dmodel/new/StaElementData.csv";
        const char *name_sta_nodeCSV    = "mesh_Dmodel/new/StaNodeData.csv";
        std::vector<std::shared_ptr<ELEMENT>> rot_elements;
        std::vector<std::shared_ptr<NODE>>    rot_nodes;
        std::vector<std::shared_ptr<ELEMENT>> sta_elements;
        std::vector<std::shared_ptr<NODE>>    sta_nodes;
        std::vector<std::shared_ptr<NODE>>    all_nodes; 
        std::vector<std::shared_ptr<ELEMENT>> all_elements;

        int num_rot_elements;
        int num_rot_nodes;
        int num_sta_elements;
        int num_sta_nodes;
        int TotalNodeNumber;
        int TotalnonBoundNodeNumber;
        int TotalElementNumber;

        int counter_slideslave;
        int counter_periodicslave;
        int counter_slidemaster;
        int counter_periodicmaster;
        int counter_dirichletboundary;
        int counter_slave;
        int counter_nonboundary;
        
        GenMatrix K;
        GenMatrix C;
        GenMatrix Ct;
        GenMatrix CtK;
        GenMatrix CtKC;
        GenMatrix L;
        std::unique_ptr<double[]> A;
        std::unique_ptr<double[]> del_A;
        std::unique_ptr<double[]> r;
        std::unique_ptr<double[]> Ctb;
        std::unique_ptr<double[]> b;
        std::unique_ptr<double[]> CA;

        // akima関連　
        const char *name_akima_coefCSV1 = "50A470/akima_coef.csv";
        const char *name_akima_coefCSV2 = "50A470/nyu_B^2.csv";
        int num_akima_coef; // 760?
        std::vector<std::vector<double>> akima_data; // B^2, p0, p1, p2, p3　double akima_data[num_akima_coef][5];と同義

        // 構造パラメータ
        double r_rotor  = 0.02775; // 回転子外形径[m]
        double r_stator = 0.07;    // 固定子外形径[m]

        // 解析条件
        double rotation_angle_deg ;  // 回転角 deg
        double M                  = 1.25; // 永久磁石磁化[T]
        double J_0                = 0.0; // 電流密度[A/mm^2]
        int shift;
        bool apply_current        = false; // 電流を負荷するかどうか
        bool nonliner             = false;
        double ICCG_conv          = 1E-6;
        double NR_conv            = 1E-3;
        double ICCG_conv_MAX;
        int max_nrloop            = 100;

        // 材料番号
        const int AIR = 0;        // 空気
        const int CORE = 1;       // コア
        const int MAGNET = 2;     // 磁石
        const int U_PLUS = 31;    // U+
        const int U_MINUS = -32;  // U-
        const int V_MINUS_1 = -41;// V- (1)
        const int V_MINUS_2 = -42;// V- (2)
        const int W_PLUS_1 = 51;  // W+ (1)
        const int W_PLUS_2 = 52;  // W+ (2)

        // vtkファイル名
        const char *name_vtk_vec = "./vtk/vectors.vtk";
        const char *name_vtk_mat = "./vtk/material.vtk";
        const char *name_vtk_con = "./vtk/contour.vtk";

        bool printinfo = true;
};
/////////////////////////////////////////////////////////////////////////////////////
void MOTOR_FEM::main_func(){
    cout << "MOTOR_FEM::main_func()" << endl;
    auto start = std::chrono::high_resolution_clock::now();

    count_lines();
    read_rot_node_CSV();
    read_rot_elementCSV();
    read_sta_node_CSV();
    read_sta_element_CSV();
    TotalnonBoundNodeNumber = TotalNodeNumber - counter_slave;
    read_akima();
    // print_mesh_info();

    preprocess_periodic();
    rotation(0.0); // deg
    sort_nodes();

    preprocess_for_make_C();
    make_C();

    NewtonRaphson_loop();
    
    write_vtk();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "処理時間: " << elapsed.count() << " 秒" << std::endl;
}
/////////////////////////////////////////////////////////////////////////////////////
MOTOR_FEM::MOTOR_FEM(){ // コンストラクタ
    counter_slideslave        = 0;
    counter_periodicslave     = 0;
    counter_slidemaster       = 0;
    counter_periodicmaster   = 0;
    counter_dirichletboundary = 0;
    counter_slave             = 0;
    counter_nonboundary       = 0;
    rotation_angle_deg        = 0.0;
}
/////////////////////////////////////////////////////////////////////////////////////
MOTOR_FEM::~MOTOR_FEM(){}
/////////////////////////////////////////////////////////////////////////////////////
void MOTOR_FEM::count_lines(){
    num_rot_nodes = count_csv_lines(name_rot_nodeCSV);
    num_rot_elements = count_csv_lines(name_rot_elementCSV);
    num_sta_nodes = count_csv_lines(name_sta_nodeCSV);
    num_sta_elements = count_csv_lines(name_sta_elementCSV);
    num_akima_coef = count_csv_lines(name_akima_coefCSV1);
    TotalNodeNumber  = num_rot_nodes + num_sta_nodes;
    TotalElementNumber = num_rot_elements + num_sta_elements;

    // cout << TotalNodeNumber << endl;

    // cout << "回転子節点CSVの行数: " << num_rot_nodes << endl;
    // cout << "回転子要素CSVの行数: " << num_rot_elements << endl;
    // cout << "固定子節点CSVの行数: " << num_sta_nodes << endl;
    // cout << "固定子要素CSVの行数: " << num_sta_elements << endl;
}
/////////////////////////////////////////////////////////////////////////////////////
void MOTOR_FEM::print_mesh_info(){
    using namespace std;
    cout << "MOTOR_FEM::print_mesh_info()" << endl;

    cout << "　　　　スレイブ節点数： " << counter_slave << endl;
    cout << "　　スライド境界節点数： " << counter_slideslave << endl;
    cout << "　　　　周期境界節点数： " << counter_periodicslave << endl;
    cout << "　　　　固定境界節点数： " << counter_dirichletboundary << endl;
    cout << "　　　　　　　総節点数： " << TotalNodeNumber << endl;
    cout << "スレイブじゃない節点数： " << TotalnonBoundNodeNumber << endl;
}
/////////////////////////////////////////////////////////////////////////////////////
int MOTOR_FEM::read_rot_node_CSV(){
    // 回転子節点CSV読み込み関数
    FILE *fp;
    char buf[256];

    fp = fopen(name_rot_nodeCSV, "r");
    if(fp == NULL){
        cout << name_rot_nodeCSV << "could not be opened." << endl;    
        return 1;
    }
    rot_nodes.reserve(num_rot_nodes);

    for(int i=0; i<num_rot_nodes; i++){
        bool isSlave = false;
        double x, y, z;

        if(fscanf(fp, "%lf,%lf,%lf", &x, &y, &z) == 0){}
        // cout << "x("<<i<<")= " << x << endl;    
        auto node = std::make_shared<NODE>();
        node->setNode(i, x, y, z);
        rot_nodes.push_back(node);

        rot_nodes[i]->cal_r_theta();
        rot_nodes[i]->isRotor = true;
        // printf("[%d]  x: %lf  y: %lf  r: %lf  theta: %lf\n",i , rot_nodes[i]->x, rot_nodes[i]->y, rot_nodes[i]->r, rot_nodes[i]->theta);

        // 周期境界slaveの節点（x=0）を探索（原点をのぞく）
        if (fabs(rot_nodes[i]->x - 0.0) <= eps && not(fabs(rot_nodes[i]->y - 0.0) <= eps)){
            // cout << "periodic! ";
            isSlave = true;
            rot_nodes[i]->isPeriodicSlave = true;
            counter_periodicslave++;
        }

        // 周期境界masterの節点（y=0）を探索（原点をのぞく）
        if (fabs(rot_nodes[i]->y - 0.0) <= eps && not(fabs(rot_nodes[i]->x - 0.0) <= eps)){
            rot_nodes[i]->isPeriodicMaster = true;
            counter_periodicmaster++;
        }

        // スライド境界のSlave節点を探索
        if (fabs(rot_nodes[i]->r - r_rotor) <= eps){ //r=27.75mmなら
            // cout << "slide!  ";
            isSlave = true;
            rot_nodes[i]->isSlideSlave = true;
            // cout << i << ": rot_nodes[i]->C_index = " << rot_nodes[i]->C_index << endl;
            counter_slideslave++;
        }

        // C_indexをつける　スレイブが後半になるように逆から入れる
        if (!rot_nodes[i]->isPeriodicSlave && !rot_nodes[i]->isSlideSlave){
            rot_nodes[i]->C_index = counter_nonboundary;
            // cout << i << ": rot_nodes[i]->C_index = " << rot_nodes[i]->C_index << endl;
            counter_nonboundary++;
        }else{
            rot_nodes[i]->C_index = TotalNodeNumber - 1 - counter_slave;
            counter_slave++;
        }
        // cout << counter_slave << "   " << rot_nodes[i]->index << ": rot_nodes[i]->C_index = " << rot_nodes[i]->C_index << endl;
        if ((!isSlave) && rot_nodes[i]->C_index > 6832) {
            cout << "Error: スレイブでない節点のC_indexが6832を超えています: " << rot_nodes[i]->C_index << endl;
        } else {
            // cout << "rot_nodes[i]->C_index: " << rot_nodes[i]->C_index << endl;
        }
    }
    fclose(fp);

    if (counter_periodicslave != counter_periodicmaster) {
        cout << "Error: counter_periodicslaveとcounter_periodicamasterの値が一致していません: " << "counter_periodicslave: " << counter_periodicslave << ", counter_periodicmaster: " << counter_periodicmaster << endl;
    }
    
    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////
int MOTOR_FEM::read_rot_elementCSV(){
    // 回転子要素CSV読み込み関数
    FILE *fp;
    char buf[256];

    fp = fopen(name_rot_elementCSV, "r");
    if(fp == NULL){
        cout << name_rot_elementCSV << "could not be opened." << endl;    
        return 1;
    }
    rot_elements.reserve(num_rot_elements);
    for(int i=0; i<num_rot_elements; i++){
        int material_number, node_idx1, node_idx2, node_idx3;

        if(fscanf(fp, "%d,%d,%d,%d", &material_number, &node_idx1, &node_idx2, &node_idx3) == 0){}
        
        // cout << "matnum("<<i<<"): " << material_number << endl;
        if(node_idx1==num_rot_nodes){cout << "error: 0オリジンにしてください" <<endl; return 2;}    
        auto element = std::make_shared<ELEMENT>();
        element->setElement(i,                             // 要素番号
                            material_number,               // 材料番号
                            rot_nodes[node_idx1],          // NODEクラスのオブジェクト1
                            rot_nodes[node_idx2],          // NODEクラスのオブジェクト2
                            rot_nodes[node_idx3]);         // NODEクラスのオブジェクト3
        rot_elements.push_back(element);
    }
    fclose(fp);
    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////
int MOTOR_FEM::read_sta_node_CSV(){
    // 固定子節点CSV読み込み関数
    FILE *fp;
    char buf[256];

    fp = fopen(name_sta_nodeCSV, "r");
    if(fp == NULL){
        cout << name_sta_nodeCSV << "could not be opened." << endl;    
        return 1;
    }
    sta_nodes.reserve(num_sta_nodes);

    for(int i=0; i<num_sta_nodes; i++){
        bool isSlave = false;
        double x, y, z;

        if(fscanf(fp, "%lf,%lf,%lf", &x, &y, &z) == 0){}
        // cout << "x("<<i<<")= " << x << endl;    
        auto node = std::make_shared<NODE>();
        node->setNode(i+num_rot_nodes, x, y, z); // i+回転子節点数を節点番号とする
        sta_nodes.push_back(node);

        sta_nodes[i]->cal_r_theta();
        sta_nodes[i]->isRotor = false;

        // 周期境界のSlave節点を探索
        if (fabs(sta_nodes[i]->x - 0.0) <= eps && !fabs(sta_nodes[i]->y - r_stator) <= eps){
            sta_nodes[i]->isPeriodicSlave = true;
            isSlave = true;
            sta_nodes[i]->C_index = TotalNodeNumber - 1 - counter_slave;
            // cout << i + num_rot_nodes << ": sta_nodes[i]->C_index = " << sta_nodes[i]->C_index << endl;
            counter_periodicslave++;
            counter_slave++;
        }else{
            //スレイブじゃない節点
            sta_nodes[i]->C_index = counter_nonboundary;
            counter_nonboundary++;
            // 固定境界の節点を探索
            if (fabs(sta_nodes[i]->r - r_stator) <= eps) { // r=70mmなら
                sta_nodes[i]->isOnDirichletBoundary = true;
                counter_dirichletboundary++;
            }
        }
        
        // スライド境界のMaster節点を探索
        if (fabs(sta_nodes[i]->r - r_rotor) <= eps){
            sta_nodes[i]->isSlideMaster = true;
            counter_slidemaster++;
        }

        // 周期境界のMaster節点を探索
        if (fabs(sta_nodes[i]->y - 0.0) <= eps && !fabs(sta_nodes[i]->x - r_stator) <= eps){
            sta_nodes[i]->isPeriodicMaster = true;
            counter_periodicmaster++;
        }

        // cout << i + num_rot_nodes << ": sta_nodes[i]->C_index = " << sta_nodes[i]->C_index << endl;

        if ((isSlave==false) && (sta_nodes[i]->C_index > 6832)) {
            cout << "Error: スレイブでない節点のC_indexが6832を超えています: " << sta_nodes[i]->C_index << endl;
        } else {
            // cout << "sta_nodes[i]->C_index: " << sta_nodes[i]->C_index << endl;
        }
    }
    fclose(fp);

    if (counter_periodicslave != counter_periodicmaster) {
        cout << "Error: counter_periodicslaveとcounter_periodicamasterの値が一致していません: " << "counter_periodicslave: " << counter_periodicslave << ", counter_periodicmaster: " << counter_periodicmaster << endl;
    }

    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////
int MOTOR_FEM::read_sta_element_CSV(){
    // 固定子要素CSV読み込み関数
    FILE *fp;
    char buf[256];

    fp = fopen(name_sta_elementCSV, "r");
    if(fp == NULL){
        cout << name_sta_elementCSV << "could not be opened." << endl;    
        return 1;
    }
    sta_elements.reserve(num_sta_elements);
    for(int i=0; i<num_sta_elements; i++){
        int material_number, node_idx1, node_idx2, node_idx3;

        if(fscanf(fp, "%d,%d,%d,%d", &material_number, &node_idx1, &node_idx2, &node_idx3) == 0){}
        
        // cout << "matnum("<<i<<"): " << material_number << endl;
        if(node_idx1==num_sta_nodes){cout << "error: 0オリジンにしてください" <<endl; return 2;}    
        auto element = std::make_shared<ELEMENT>();
        element->setElement(i,                             // 要素番号
                            material_number,               // 材料番号
                            sta_nodes[node_idx1],          // NODEクラスのオブジェクト1
                            sta_nodes[node_idx2],          // NODEクラスのオブジェクト2
                            sta_nodes[node_idx3]);         // NODEクラスのオブジェクト3
        sta_elements.push_back(element);
    }
    fclose(fp);
    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////
int MOTOR_FEM::read_akima(){

    akima_data.resize(num_akima_coef, std::vector<double>(5, 0.0));

    FILE *fp;

    fp = fopen(name_akima_coefCSV1, "r");
    if(fp == NULL){
        cout << name_akima_coefCSV1 << " could not be opened." << endl;    
        return 1;
    }
    for(int i=0; i<num_akima_coef; i++){
        if(fscanf(fp, "%lf,%lf,%lf,%lf", &akima_data[i][1], &akima_data[i][2], &akima_data[i][3], &akima_data[i][4]) == 0){}
    }
    fclose(fp);

    fp = fopen(name_akima_coefCSV2, "r");
    if(fp == NULL){
        cout << name_akima_coefCSV2 << " could not be opened." << endl;    
        return 1;
    }
    for(int i=0; i<num_akima_coef; i++){
        if(fscanf(fp, "%lf,%*lf", &akima_data[i][0]) == 0){}
    }

    // （確認用）akima_dataを出力
    // for(int i = 0; i < num_akima_coef; i++) {
    //     for(int j = 0; j < 5; j++) {
    //         printf("%10f  ", akima_data[i][j]);
    //     }
    //     printf("\n");
    // }

    fclose(fp);
    
    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////
int MOTOR_FEM::preprocess_periodic(){
    // cout << "MOTOR_FEM::preprocess_periodic()" << endl;
    // // マスタースレイブ関係認識（周期境界）
    // for(int i=0; i<rot_nodes.size(); i++){
    //     if(rot_nodes[i]->isPeriodicSlave){
    //         for(const auto& master_node : sta_nodes){
    //             if(fabs(rot_nodes[i]->y - master_node->x) < eps && master_node->isPeriodicMaster){
    //                 rot_nodes[i]->PeriodicMasterIndex = master_node->C_index;
    //                 break;
    //             }
    //         }
    //     }
    // }
    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////
void MOTOR_FEM::rotation(double angle) {
    for (auto& node : rot_nodes) {
        node->add_theta(angle*PI/180.0);
    }
    rotation_angle_deg += angle;
}
/////////////////////////////////////////////////////////////////////////////////////
void MOTOR_FEM::sort_nodes(){
    // sta_nodesを角度theta順に並べ替える　なおスライド境界を優先する
    std::sort(sta_nodes.begin(), sta_nodes.end(), [](const std::shared_ptr<NODE>& a, const std::shared_ptr<NODE>& b) {
        if (a->isSlideMaster && !b->isSlideMaster) {
            return true;
        } else if (!a->isSlideMaster && b->isSlideMaster) {
            return false;
        } else {
            return a->theta < b->theta;
        }
    });
    // rot_nodesを角度theta順に並べ替える　なおスライド境界を優先する
    std::sort(rot_nodes.begin(), rot_nodes.end(), [](const std::shared_ptr<NODE>& a, const std::shared_ptr<NODE>& b) {
        if (a->isSlideSlave && !b->isSlideSlave) {
            return true;
        } else if (!a->isSlideSlave && b->isSlideSlave) {
            return false;
        } else {
            return a->theta < b->theta;
        }
    });
    // for(const auto& node : sta_nodes){
    //     if(node->isSlideMaster){
    //        cout << node->index << ": " << node->theta << endl;
    //     }
    // }
    // for(const auto& node : sta_nodes){
    //     cout << node->theta << endl;
    // }
}
/////////////////////////////////////////////////////////////////////////////////////
int MOTOR_FEM::preprocess_for_make_C(){
    cout << "MOTOR_FEM::preprocess_for_make_C()" << endl;  

    C = GenMatrix(TotalNodeNumber, TotalnonBoundNodeNumber);	
    K = GenMatrix(TotalNodeNumber, TotalNodeNumber); // 全体行列
    Ct = GenMatrix(TotalnonBoundNodeNumber,TotalNodeNumber);
    CtK = GenMatrix(TotalnonBoundNodeNumber, TotalNodeNumber);
    CtKC = GenMatrix(TotalnonBoundNodeNumber, TotalnonBoundNodeNumber);//計算中に用いる全体行列
    L = GenMatrix(TotalnonBoundNodeNumber, TotalnonBoundNodeNumber);
    A = std::make_unique<double[]>(TotalnonBoundNodeNumber);//計算中に用いる解ベクトル
    r = std::make_unique<double[]>(TotalnonBoundNodeNumber);//残差
    Ctb = std::make_unique<double[]>(TotalnonBoundNodeNumber);//計算中に用いる右辺ベクトル
    b = std::make_unique<double[]>(TotalNodeNumber);//最終的な右辺ベクトル
    CA = std::make_unique<double[]>(TotalNodeNumber);
    del_A = std::make_unique<double[]>(TotalNodeNumber);//最終的な解ベクトル

    int counter_slide = 0;
    int counter_periodic = 0;
    int count = 0;

    // マスターとスレイブの関係認識
    for(int i = 0; i < rot_nodes.size(); i++){ //回転子ノードに対して

        // 周期境界slaveについて
        if(rot_nodes[i]->isPeriodicSlave){
            count++;
            // cout << "count: " << count << endl;
            for(const auto& master_node : rot_nodes){
                // 周期境界のマスターを探す(r_slave == r_master)
                if(fabs(rot_nodes[i]->r - master_node->r) < 1.0e-5 && master_node->isPeriodicMaster){ 
                    rot_nodes[i]->PeriodicMasterIndex = master_node->C_index;
                    // デバッグ出力
                    // cout << "rot_nodes[" << i << "]->PeriodicMasterIndex: " << rot_nodes[i]->PeriodicMasterIndex << endl;
                    // cout << "rot_nodes[" << i << "]->r: " << rot_nodes[i]->r << ", master_node->r: " << master_node->r << endl;
                    counter_periodic++;
                    // cout << "counter_periodic: " << counter_periodic << endl;
                    break;
                }
            }
        }

        // スライド境界について
        if(rot_nodes[i]->isSlideSlave){

            int idx1, idx2;
            // 最初のノードだけ内分点探索　shiftを求めておく
            if(i == 0){
                for (int j = 0; j < counter_slideslave; j++) { // 0~180
                    idx1 = j;
                    idx2 = (j + 1) % counter_slideslave;
                    if(idx2 == 0){
                        idx2 = counter_slideslave-1;
                    }
                    const auto& master_node1 = sta_nodes[idx1];
                    const auto& master_node2 = sta_nodes[idx2];
                    double t = (fmod(rot_nodes[i]->theta, PI/2.0) - master_node1->theta) / (master_node2->theta - master_node1->theta); // PI/2 超えたら，超えた分を使う fmod
                    if (t >= 0.0 && t < 1.0) {
                        shift = idx1;
                        cout << "shift: " << shift << ", t: " << t << endl;
                        // cout << "idx1: " << idx1 << ", idx2: " << idx2 << endl;
                        break;
                    }
                }
                rot_nodes[i]->SlideMasterIndex_1 = sta_nodes[idx1]->C_index;
                rot_nodes[i]->SlideMasterIndex_2 = sta_nodes[idx2]->C_index;
                // デバッグ出力
                // cout << "rot_nodes[" << i << "]->SlideMasterIndex_1: " << rot_nodes[i]->SlideMasterIndex_1 << ", SlideMasterIndex_2: " << rot_nodes[i]->SlideMasterIndex_2 << endl;
                counter_slide++;
            }else{
                // shiftを使ってMasterを認識する（固定子・回転子接触部分が同じ節点数のときという前提のもとで）
                idx1 = (shift + i) % (counter_slideslave-1);
                idx2 = (shift + i + 1) % (counter_slideslave-1);
                if(idx2 == 0){
                    idx2 = counter_slideslave-1;
                }
                // cout << "idx1: " << idx1 << ", idx2: " << idx2 << endl;
                if (sta_nodes[idx1]->isSlideMaster && sta_nodes[idx2]->isSlideMaster) {
                    rot_nodes[i]->SlideMasterIndex_1 = sta_nodes[idx1]->C_index;
                    rot_nodes[i]->SlideMasterIndex_2 = sta_nodes[idx2]->C_index;
                } else {
                    cout << "Error: スライドマスターが見つかりません: idx1: " << idx1 << ", idx2: " << idx2 << endl;
                }
                // デバッグ出力
                // cout << "rot_nodes[" << i << "]->SlideMasterIndex_1: " << rot_nodes[i]->SlideMasterIndex_1 << ", SlideMasterIndex_2: " << rot_nodes[i]->SlideMasterIndex_2 << endl;
                counter_slide++;
            }
        }
    }

    for(int i=0; i<sta_nodes.size(); i++){
        // 周期境界について
        if(sta_nodes[i]->isPeriodicSlave){
            for(const auto& master_node : sta_nodes){
                // 周期境界のマスターを探す(r_slave == r_master)
                if(fabs(sta_nodes[i]->r - master_node->r) < 1.0e-5 && master_node->isPeriodicMaster){ 
                    sta_nodes[i]->PeriodicMasterIndex = master_node->C_index;
                    counter_periodic++;
                    // cout << "counter_periodic: " << counter_periodic << endl;
                    break;
                }
            }
        }
    }


    if (counter_slide != counter_slideslave || counter_periodic != counter_periodicslave) {
        cout << "Error: counterの値が正しくありません: " << counter_slide << ", " << counter_periodic << endl;
    }

    // rot_nodesとsta_nodesを結合したall_nodesを作成する．その際，C_index順になるようにする．
    all_nodes.reserve(rot_nodes.size() + sta_nodes.size());
    for (const auto& node : rot_nodes) {
        all_nodes.push_back(node);
    }
    for (const auto& node : sta_nodes) {
        all_nodes.push_back(node);
    }
    std::sort(all_nodes.begin(), all_nodes.end(), [](const std::shared_ptr<NODE>& a, const std::shared_ptr<NODE>& b) {
        return a->C_index < b->C_index;
    });

    // rot_elementsとsta_elementsを結合したall_elementsを作成する
    all_elements.reserve(rot_elements.size() + sta_elements.size());
    for (const auto& element : rot_elements) {
        all_elements.push_back(element);
    }
    for (const auto& element : sta_elements) {
        all_elements.push_back(element);
    }

    // cout << "all_elementsの情報を出力します:" << endl;
    // for (const auto& element : all_elements) {
    //     // element->print_info();
    // }

    // for (const auto& node : all_nodes){
    //     if(node->SlideMasterIndex_1 < 0 || node->SlideMasterIndex_2 < 0 || node->PeriodicMasterIndex < 0){
            
    //         cout << "index: " << node->index 
    //              << ", C_index: " << node->C_index 
    //              << ", SlideMasterIndex_1: " << node->SlideMasterIndex_1 
    //              << ", SlideMasterIndex_2: " << node->SlideMasterIndex_2 
    //              << ", PeriodicMasterIndex: " << node->PeriodicMasterIndex << endl;
    //     }
    // }


    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////
int MOTOR_FEM::make_C(){
    // C行列作成
    cout << "MOTOR_FEM::make_C()" << endl;

    C.quickAllZero();
    for(int i=0; i<TotalnonBoundNodeNumber; i++){
        C.set(i, i , 1.0);
    }
    for(const auto& node : all_nodes){
        if(node->isSlideSlave){
            if (node->SlideMasterIndex_1 > TotalnonBoundNodeNumber || node->SlideMasterIndex_2 > TotalnonBoundNodeNumber) {
                // cout << "node->C_index: " << node->C_index << ", node->SlideMasterIndex_1: " << node->SlideMasterIndex_1 << ", node->SlideMasterIndex_2: " << node->SlideMasterIndex_2 << endl;
                int column1 = node->SlideMasterIndex_1;
                int column2 = all_nodes[node->SlideMasterIndex_2]->PeriodicMasterIndex;
                bool isAntiPeriodic = (fmod(node->theta, 2*PI) >= PI/2.0 
                                && fmod(node->theta, 2*PI) < PI) 
                                || (fmod(node->theta, 2*PI) >= 3*PI/2.0 
                                && fmod(node->theta, 2*PI) < 2*PI);
                node->virtual_theta = fmod(node->theta, PI/2.0);
                node->virtual_x = node->r * cos(node->virtual_theta);
                node->virtual_y = node->r * sin(node->virtual_theta);

                double d1 = sqrt(pow(node->virtual_x - all_nodes[node->SlideMasterIndex_1]->virtual_x, 2) + 
                                pow(node->virtual_y - all_nodes[node->SlideMasterIndex_1]->virtual_y, 2));
                double d2 = sqrt(pow(node->virtual_x - all_nodes[node->SlideMasterIndex_2]->virtual_x, 2) + 
                                pow(node->virtual_y - all_nodes[node->SlideMasterIndex_2]->virtual_y, 2));
                double weight1 = - d2 / (d1 + d2);
                double weight2 = - d1 / (d1 + d2);
                if(isAntiPeriodic){
                    weight1 *= -1.0;
                }
                C.set(node->C_index, column1, weight1);
                C.set(node->C_index, column2, weight2);
            }else{
                int column1 = node->SlideMasterIndex_1;
                int column2 = node->SlideMasterIndex_2;
                bool isAntiPeriodic = (fmod(all_nodes[node->SlideMasterIndex_1]->theta, 2*PI) >= PI/2.0 
                                && fmod(all_nodes[node->SlideMasterIndex_1]->theta, 2*PI) < PI) 
                                || (fmod(all_nodes[node->SlideMasterIndex_1]->theta, 2*PI) >= 3*PI/2.0 
                                && fmod(all_nodes[node->SlideMasterIndex_1]->theta, 2*PI) < 2*PI);
                node->virtual_theta = fmod(node->theta, PI/2.0);
                node->virtual_x = node->r * cos(node->virtual_theta);
                node->virtual_y = node->r * sin(node->virtual_theta);

                double d1 = sqrt(pow(node->virtual_x - all_nodes[node->SlideMasterIndex_1]->virtual_x, 2) + 
                                pow(node->virtual_y - all_nodes[node->SlideMasterIndex_1]->virtual_y, 2));
                double d2 = sqrt(pow(node->virtual_x - all_nodes[node->SlideMasterIndex_2]->virtual_x, 2) + 
                                pow(node->virtual_y - all_nodes[node->SlideMasterIndex_2]->virtual_y, 2));
                double weight1 = - d2 / (d1 + d2);
                double weight2 = - d1 / (d1 + d2);
                if(isAntiPeriodic){
                    weight1 *= -1.0;
                    weight2 *= -1.0;
                }
                C.set(node->C_index, column1, weight1);
                C.set(node->C_index, column2, weight2);
            }
        }
        else if(node->isPeriodicSlave){
            int column = node->PeriodicMasterIndex;
            if (column == -2) {
                cout << "node->C_index: " << node->C_index << ", node->PeriodicMasterIndex: " << node->PeriodicMasterIndex << endl;
            }
            C.set(node->C_index, column, -1.0);
        }
    }
    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////
int MOTOR_FEM::NewtonRaphson_loop(){
    // K行列作成
    cout << "MOTOR_FEM::NewtonRaphson_loop()" << endl;


    for(int i=0; i<TotalNodeNumber; i++){
        CA[i] = 0.0;
    }

    auto start = chrono::high_resolution_clock::now();

    // 全要素に対して定数計算を行う
    for(const auto& element : all_elements){
        element->cal_constants();
    }
    // NR法スタート
    for(int nrloop=0; nrloop<max_nrloop; nrloop++){

        // 配列初期化
        K.quickAllZero();
        CtK.quickAllZero();
        CtKC.quickAllZero();
        for (int i = 0; i < TotalnonBoundNodeNumber; i++) {
            Ctb[i] = 0.0;
            r[i] = 0.0;
        }
        for (int i = 0; i < TotalNodeNumber; i++) {
            b[i] = 0.0;
            del_A[i] = 0.0;
            // CA[i] = 0.0;
        }

        // K, b 作成
        for (const auto& e : all_elements) {
            int mat = e->materialNumber;
            // U更新
            e->cal_U();
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    double tmp = e->nyu * e->S[i][j] + 2.0 * e->dnyu_dB2 * e->U[i] * e->U[j] / e->delta;
                    
                    int row = e->nodes[i]->C_index; 
                    int col = e->nodes[j]->C_index;
                    K.add(row, col, tmp);
                    // cout << "K[" << row << "][" << col << "]: " << K.get(row, col) << endl;
                    // if(tmp != 0){
                    //     cout << "tmp: " << tmp << endl;
                    // }
                }
                b[e->nodes[i]->C_index] += (e->J_0 * e->delta / 3.0) - (e->nyu * e->U[i]);
                if(mat == MAGNET){ // 材料が磁石なら
                    // e->print_info();
                    b[e->nodes[i]->C_index] += e->nyu * (M * cos(rotation_angle_deg * PI / 180.0 + PI/4.0) * e->d[i] - M * sin(rotation_angle_deg * PI / 180.0 + PI/4.0) * e->c[i]) / 2.0;
                }
                // if (b[e->nodes[i]->C_index] != 0) {
                //     cout << "b[" << e->nodes[i]->C_index << "]: " << b[e->nodes[i]->C_index] << endl;
                // }
            }
            // if (mat == MAGNET) {
            //     cout << "e->nodes[0]->isRotor: " << e->nodes[0]->isRotor << endl;
            // }
        }


        GenMatrix::productGtG(C, K, CtK);
        GenMatrix::productGG(CtK, C, CtKC);//CtKC
        GenMatrix::productTranAX(C, b.get(), Ctb.get());//Ctb
        // 固定境界の行番号を格納する
        vector<int> v1;
        for(const auto& node : all_nodes){
            if(node->isOnDirichletBoundary){
                v1.push_back(node->C_index);      
            }
        }
        // 境界条件の設定
        K.setBoundaryCondition(Ctb.get(), v1, 0.0);

        double norm = 0.0;
        for (int i = 0; i<TotalnonBoundNodeNumber; i++) {
            norm += Ctb[i] * Ctb[i];
        }
        double Bnorm = sqrt(norm);
        // cout << "Bnorm: " << Bnorm << endl;

        // 不完全コレスキー分解
        GenMatrix::icdcmp(CtKC, L, 1.05);
        // ICCG本体
        GenMatrix::iccgSolv(CtKC, L, Ctb.get(), A.get(), r.get(), ICCG_conv*Bnorm, ICCG_conv_MAX=TotalnonBoundNodeNumber);
        // // 解ベクトル
        GenMatrix::productAX(C, A.get(), del_A.get());

        for(int i=0; i<TotalNodeNumber; i++){
            CA[i] += del_A[i];
            all_nodes[i]->add_A(del_A[i]);
            // cout << "del_A[" << i << "]: " << del_A[i] << ", CA[" << i << "]: " << CA[i] << endl;
        }

        for(const auto& element : all_elements){
            element->update_B();
        }

        auto end = chrono::high_resolution_clock::now();

        if(nonliner){
            double norm_del_A = 0.0;
            for(int i=0; i<TotalNodeNumber; i++){
                norm_del_A += del_A[i] * del_A[i];
            }
            norm_del_A = sqrt(norm_del_A);
            // cout << "norm_del_A: " << norm_del_A << endl;
            if(norm_del_A < NR_conv){
                if(printinfo == true){
                    cout << "NR loop: " << nrloop + 1 
                         << "   norm_delta_A: " << setw(12) << left <<  norm_del_A
                         << "   time[ms]: " << setw(4) << left <<  chrono::duration_cast<chrono::milliseconds>(end - start).count()
                         << "   (終了)\n";
                }
                break;
            }else{
                if(printinfo == true){
                    cout << "NR loop: " << nrloop + 1 
                         << "   norm_delta_A: " << setw(12) << left <<  norm_del_A
                         << "   time[ms]: " << setw(4) << left <<  chrono::duration_cast<chrono::milliseconds>(end - start).count()
                         << "   (続行)\n";
                }
                for(const auto& element : all_elements){
                    if(element->materialNumber == CORE){ // 鉄なら                        
                        element->update_nyu(akima_data);
                    }
                }
            }
        }else{
            break;
        }
        
    }
    // NR法終了



    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////
int MOTOR_FEM::write_vtk(){
    cout << "MOTOR_FEM::write_vtk()" << endl;
    write_vtk_material_3(name_vtk_mat, all_nodes, all_elements, TotalNodeNumber, TotalElementNumber);
    write_vtk_vector_3(name_vtk_vec, all_nodes, all_elements, TotalNodeNumber, TotalElementNumber);
    write_vtk_contour_3(name_vtk_con, all_nodes, all_elements, TotalNodeNumber, TotalElementNumber, CA.get());

    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////
PYBIND11_MODULE(MOTOR_FEM, m) {
    m.doc() = "This is a test module";
    py::class_<MOTOR_FEM>(m, "MOTOR_FEM")
		.def(py::init<>())
        .def("main_func", &MOTOR_FEM::main_func);
		
}
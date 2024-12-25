#include <vector>
#include <memory>
#include "GenMatrix.h" // Added this line to include the GenMatrix header

// 節点クラス
class NODE{
    private:
    public:
        NODE();
        void setNode(int node_index, double node_x, double node_y, double node_z);
        void add_A(double delta_A);
        void cal_r_theta();
        void set_isOnSlideBoundary(bool);
        void set_isOnPeriodicBoundary(bool);
        void set_isOnDirichletBoundary(bool);
        void add_theta(double theta);

        int index;         // rot と sta　足し合わせてる　0~（全節点数-1） 
        int C_index;       // C行列の何列目に入れるか
        double x;
        double y;
        double z;
        double r;
        double theta;
        double virtual_theta; // fmod(theta, PI/2)
        double virtual_x;
        double virtual_y;
        bool isOnBoundary;
        double Ax;
        double Ay;
        double Az; // ベクトルポテンシャルz成分
        double del_Ax;
        double del_Ay;
        double del_Az;

        bool isSlideSlave;
        bool isPeriodicSlave;
        bool isSlideMaster;
        bool isPeriodicMaster;
        bool isOnDirichletBoundary;
        bool isRotor;

        int PeriodicMasterIndex;
        int SlideMasterIndex_1;
        int SlideMasterIndex_2;

        
};
// 要素クラス
class ELEMENT{
    private:
    public:
        ELEMENT();
        void setElement(int index, int material_number, std::shared_ptr<NODE> node1, std::shared_ptr<NODE> node2, std::shared_ptr<NODE> node3);
        void print_info();
        void set_material_config();
        void cal_constants();
        void update_B();
        void update_nyu(std::vector<std::vector<double>>& akima_data);
        void cal_U();
        void update_node(std::shared_ptr<NODE> node1, std::shared_ptr<NODE> node2, std::shared_ptr<NODE> node3);
        
        int index;
        int materialNumber;
        int nodeNumbers[3];
        std::shared_ptr<NODE> nodes[3];         // 囲まれてる節点３つ
        double c[3];
        double d[3];
        double delta;
        double S[3][3];
        double U[3];
        bool isCurrentFlowing; // 電流流れてるか否か
        double J_0;            // 強制電流密度
        double Bx;             // 磁束密度x成分
        double By;             // 磁束密度y成分
        double B2;             // 磁束密度の大きさの２乗
        double nyu;            // 透磁率
        double dnyu_dB2;       // 透磁率をB^2で微分した値
};

void swap(int *a, int *b);
void Readmesh(const char *filename);
// void write_vtk_material(const char *filename, double **node, int **element, int numnode, int numele);
// void write_vtk_vector(const char *filename, double **node, int **element, double *B_x, double *B_y, int numnode, int numele);
// void write_vtk_contour(const char *filename, double **node, int **element, double *B_x, double *B_y, int numnode, int numele, double *A);
void write_vtk_material_2(const char *filename, std::vector<std::unique_ptr<NODE>>& nodes, std::vector<std::unique_ptr<ELEMENT>>& elements, int numnode, int numele);
void write_vtk_vector_2(const char *filename, std::vector<std::unique_ptr<NODE>>& nodes, std::vector<std::unique_ptr<ELEMENT>>& elements, int numnode, int numele);
void write_vtk_contour_2(const char *filename, std::vector<std::unique_ptr<NODE>>& nodes, std::vector<std::unique_ptr<ELEMENT>>& elements, int numnode, int numele, double *A);
void write_vtk_material_3(const char *filename, std::vector<std::shared_ptr<NODE>>& nodes, std::vector<std::shared_ptr<ELEMENT>>& elements, int numnode, int numele);
void write_vtk_vector_3(const char *filename, std::vector<std::shared_ptr<NODE>>& nodes, std::vector<std::shared_ptr<ELEMENT>>& elements, int numnode, int numele);
void write_vtk_contour_3(const char *filename, std::vector<std::shared_ptr<NODE>>& nodes, std::vector<std::shared_ptr<ELEMENT>>& elements, int numnode, int numele, double *A);
int count_csv_lines(const char* filename);
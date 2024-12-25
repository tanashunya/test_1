#include <vector>
#include <memory>

// 節点クラス
class NODE{
    private:
    public:
        NODE();
        void setNode(int node_index, double node_x, double node_y, double node_z);
        void add_A(double delta_A);
        int index;
        double x;
        double y;
        double z;
        bool isOnBoundary;
        double Ax;
        double Ay;
        double Az; // ベクトルポテンシャルz成分
        double del_Ax;
        double del_Ay;
        double del_Az;
};
// 要素クラス
class ELEMENT{
    private:
    public:
        ELEMENT();
        void setElement(int index, int matNum, NODE& node1, NODE& node2, NODE& node3);
        void print_info();
        void set_material_config();
        void cal_constants();
        void update_B();
        void update_nyu(double akima_data[38][5]);
        void cal_U();
        void update_node(NODE& node1, NODE& node2, NODE& node3);
        
        int index;
        int materialNumber;
        int nodeNumbers[3];
        NODE nodes[3];         // 囲まれてる節点３つ
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


#include "node.h"
#include "element.h"
#include <vector>
#include <memory>
#include <iostream>

class FEM_3D{
    private:
    public:
        FEM_3D();
        ~FEM_3D();
        void main_func();
        std::vector<std::shared_ptr<NODE>> nodes;
        std::vector<std::shared_ptr<ELEMENT_HEXAHEDRAL_EDGE>> elements;
};
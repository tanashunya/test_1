#include "node.h"
#include "edge.h"
#include <memory>
class ELEMENT_HEXAHEDRAL_EDGE{
    private:
    public:
        ELEMENT_HEXAHEDRAL_EDGE();
        ~ELEMENT_HEXAHEDRAL_EDGE();
        std::shared_ptr<EDGE> edges[12];
        std::shared_ptr<NODE> nodes[8];
};
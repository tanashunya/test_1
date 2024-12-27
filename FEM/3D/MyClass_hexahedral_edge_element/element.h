#pragma once

#include "node.h"
#include <memory>
class ELEMENT{
    private:
    public:
        ELEMENT();
        std::shared_ptr<NODE> nodes[8];
};
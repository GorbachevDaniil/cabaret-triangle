#ifndef Mesh_hpp
#define Mesh_hpp

#include "Cell.hpp"

#include <vector>

class Mesh {
public:
    std::vector<Cell> cells;
};

#endif
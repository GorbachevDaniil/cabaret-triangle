#ifndef Mesh_hpp
#define Mesh_hpp

#include "Node.hpp"
#include "Edge.hpp"
#include "Cell.hpp"

#include <vector>

class Mesh {
public:
    std::vector<Node> nodes;
    std::vector<Edge> edges;
    std::vector<Cell> cells;
};

#endif
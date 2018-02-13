#ifndef Edge_hpp
#define Edge_hpp

#include "Node.hpp"
#include "Cell.hpp"

#include <vector>

class Cell;

class Edge {
public:
    long id;

    std::vector<Node> nodes; // TODO we need order here
    std::vector<Cell> neigCells; // TODO we need order here
};

#endif
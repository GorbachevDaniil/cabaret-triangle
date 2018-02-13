#ifndef Edge_hpp
#define Edge_hpp

#include "Node.hpp"

#include <vector>

class Edge {
public:
    long id;

    Node node;

    std::vector<Cell> neigCells(2); // TODO we need order here
};

#endif
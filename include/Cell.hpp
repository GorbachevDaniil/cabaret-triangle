#ifndef Cell_hpp
#define Cell_hpp

#include "Node.hpp"
#include "Edge.hpp"

#include <vector>

class Cell {
public:
    long id;
    double volume;

    Node center;

    std::vector<Edge> neigEdges; // TODO we need order here for calculating div
};

#endif
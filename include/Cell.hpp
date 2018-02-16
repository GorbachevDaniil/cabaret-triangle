#ifndef Cell_hpp
#define Cell_hpp

#include "Node.hpp"
#include "Edge.hpp"

#include <vector>
#include <map>

class Cell {
public:
    long ID;
    long centerNodeID;
    double volume;

    std::vector<long> edgeIDs; // TODO do we need order here for calculating div?
    std::map<long, int> edgeIDNormalDirection;
};

#endif
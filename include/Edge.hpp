#ifndef Edge_hpp
#define Edge_hpp

#include <vector>

class Cell;

class Edge {
public:
    long ID;
    long centerNodeID;
    bool BoundEdge;

    std::vector<long> nodeIDs; // TODO do we need order here?
    std::vector<long> cellIDs; // TODO do we need order here?
};

#endif
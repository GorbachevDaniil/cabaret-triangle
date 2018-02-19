#ifndef Edge_hpp
#define Edge_hpp

#include <vector>

class Cell;

class Edge {
public:
    long ID;
    long centerNodeID;
    double length;
    bool boundEdge;

    std::vector<long> nodeIDs;
    std::vector<long> cellIDs;
};

#endif
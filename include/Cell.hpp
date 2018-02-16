#ifndef Cell_hpp
#define Cell_hpp

#include <vector>
#include <map>

class Cell {
public:
    long ID;
    long centerNodeID;
    double volume;

    std::vector<long> nodeIDs; // TODO we need order here for calculating normals
    std::vector<long> edgeIDs; // TODO do we need order here for calculating div?
    std::map<long, int> edgeIDNormalDirection;
};

#endif
#ifndef Cell_hpp
#define Cell_hpp

#include <vector>
#include <map>

class Cell {
public:
    long ID;
    long centerNodeID;
    double volume;

    std::vector<long> nodeIDs;
    std::vector<long> edgeIDs;
    std::map<long, int> edgeIDNormalDirection;

    long getNextNodeID(unsigned long nodeIDPos);
    long getPrevNodeID(unsigned long nodeIDPos);
    std::vector<long> getEdgeOrderedNodeIDs(std::vector<long> edgeUnorderedNodeIDs);
};

#endif
#ifndef Cell_hpp
#define Cell_hpp

#include "Vector.hpp"

#include <vector>
#include <map>

class Mesh;

class Cell {
public:
    long ID;
    long centerNodeID;
    double volume;
    double maxH;

    std::vector<long> nodeIDs;
    std::vector<long> edgeIDs;
    std::map<long, int> edgeToNormalDir;
    std::map<long, Vector> edgeToTransportDir;
    std::map<long, double> edgeToMedianLength;

    Cell() {};
    Cell(Mesh &mesh, long id, long node_id_1, long node_id_2, long node_id_3);
    double countVolume(Mesh &mesh, long node_id_1, long node_id_2, long node_id_3);
    long getNextNodeID(unsigned long nodeIDPos);
    long getPrevNodeID(unsigned long nodeIDPos);
    std::vector<long> getEdgeOrderedNodeIDs(std::vector<long> edgeUnorderedNodeIDs);
};

#endif

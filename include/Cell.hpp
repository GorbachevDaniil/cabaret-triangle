#ifndef Cell_hpp
#define Cell_hpp

#include "Vector.hpp"

#include <map>
#include <vector>
#include <armadillo>

class Mesh;

class Cell {
   public:
    long ID;
    long centerNodeID;
    double volume;

    std::vector<long> nodeIDs;
    std::vector<long> edgeIDs;
    std::map<long, int> edgeToNormalDir;
    std::map<long, Vector> nodeToTransferVector;
    std::map<long, long> nodeIDToOppositeNodeID;

    arma::mat interpolationMat;

    Cell(){};
    Cell(Mesh &mesh, long ID, long nodeID1, long nodeID2, long nodeID3);
    double countVolume(Mesh &mesh, long nodeID1, long nodeID2, long nodeID3);
    long getNextNodeID(unsigned long nodeIDPos);
    long getPrevNodeID(unsigned long nodeIDPos);
    std::vector<long> getEdgeOrderedNodeIDs(std::vector<long> edgeUnorderedNodeIDs);

   private:
    void assignOppositeNodeIDs(Mesh &mesh);
    void buildInterpolationMat(Mesh &mesh);
};

#endif

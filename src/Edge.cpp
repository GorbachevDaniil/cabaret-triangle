#include "Edge.hpp"

#include "Mesh.hpp"
#include "Node.hpp"

#include <cassert>
#include <cmath>

Edge::Edge(Mesh &mesh, long ID, long startNodeID, long endNodeID, bool boundEdge,
           int innerNodeNum) {
    this->ID = ID;
    this->boundEdge = boundEdge;
    endNodeIDs.push_back(startNodeID);
    endNodeIDs.push_back(endNodeID);

    double startNodeX, startNodeY, endNodeX, endNodeY;
    startNodeX = mesh.nodes[startNodeID].data.coords.x;
    startNodeY = mesh.nodes[startNodeID].data.coords.y;
    endNodeX = mesh.nodes[endNodeID].data.coords.x;
    endNodeY = mesh.nodes[endNodeID].data.coords.y;

    length = sqrt(pow((endNodeX - startNodeX), 2) + pow((endNodeY - startNodeY), 2));

    // In nodeIDs there will be always nodes order from start node to end node
    nodeIDs.push_back(startNodeID);
    Vector startNode(startNodeX, startNodeY);
    Vector edgeVector(endNodeX - startNodeX, endNodeY - startNodeY);
    for (int i = 1; i < innerNodeNum + 1; i++) {
        Vector nodeCoords = edgeVector / (innerNodeNum + 1) * i + startNode;
        Node *node = new Node(mesh, nodeCoords.x, nodeCoords.y, true, false);
        mesh.nodes.push_back(*node);
        nodeIDs.push_back(node->ID);
    }
    nodeIDs.push_back(endNodeID);

    std::pair<int, int> tempNodeIDs;
    tempNodeIDs = std::make_pair(startNodeID, endNodeID);
    mesh.mapNodesWithEdge.insert(std::pair<std::pair<int, int>, int>(tempNodeIDs, ID));
    tempNodeIDs = std::make_pair(endNodeID, startNodeID);
    mesh.mapNodesWithEdge.insert(std::pair<std::pair<int, int>, int>(tempNodeIDs, ID));
}

std::vector<long> Edge::getUsedNodes(Mesh &mesh) {
    std::vector<long> usedNodeIDs;
    for (long nodeID : nodeIDs) {
        if (mesh.nodes[nodeID].used) {
            usedNodeIDs.push_back(nodeID);
        }
    }
    assert(usedNodeIDs.size() > 0);
    return usedNodeIDs;
}

std::vector<long> Edge::getInnerNodes() {
    std::vector<long>::const_iterator firstInnerNode = nodeIDs.begin() + 1;
    std::vector<long>::const_iterator lastInnerNode = nodeIDs.end() - 1;
    std::vector<long> innerNodeIDs(firstInnerNode, lastInnerNode);
    assert(innerNodeIDs.size() > 0);
    return innerNodeIDs;
}

#include "grid/edge.hpp"

#include <cassert>
#include <cmath>

#include "grid/mesh.hpp"

Edge::Edge(Mesh &mesh, long ID, long startNodeID, long endNodeID, bool boundEdge, int innerNodeNum) {
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
    this->nodeIDs.push_back(startNodeID);
    mesh.nodes[startNodeID].edgeIDs.insert(ID);
    Vector startNode(startNodeX, startNodeY);
    Vector edgeVector(endNodeX - startNodeX, endNodeY - startNodeY);
    for (int i = 1; i < innerNodeNum + 1; i++) {
        Vector nodeCoords = edgeVector / (innerNodeNum + 1) * i + startNode;
        Node *node = new Node(mesh, nodeCoords.x, nodeCoords.y, true,
                              boundEdge, false, false);
        node->edgeIDs.insert(ID);
        mesh.nodes.push_back(*node);
        this->nodeIDs.push_back(node->ID);
        this->innerNodeIDs.push_back(node->ID);
    }
    this->nodeIDs.push_back(endNodeID);
    mesh.nodes[endNodeID].edgeIDs.insert(ID);

    for (long nodeID : this->nodeIDs) {
        if (mesh.nodes[nodeID].used) {
            this->usedNodeIDs.push_back(nodeID);
        }
    }

    std::pair<int, int> tempNodeIDs;
    tempNodeIDs = std::make_pair(startNodeID, endNodeID);
    mesh.map_nodes_with_edge.insert(std::pair<std::pair<int, int>, int>(tempNodeIDs, ID));
    tempNodeIDs = std::make_pair(endNodeID, startNodeID);
    mesh.map_nodes_with_edge.insert(std::pair<std::pair<int, int>, int>(tempNodeIDs, ID));
}

long Edge::getAnotherEndNode(long ID) {
    assert((ID == nodeIDs.front()) || (ID == nodeIDs.back()));
    if (ID == nodeIDs.front()) {
        return nodeIDs.begin()[3];
    } else {
        return nodeIDs.begin()[0];
    }
}

long Edge::getNearInnerNode(long ID) {
    assert((ID == nodeIDs.front()) || (ID == nodeIDs.back()));
    if (ID == nodeIDs.front()) {
        return nodeIDs.begin()[1];
    } else {
        return nodeIDs.begin()[2];
    }
}

long Edge::getFarInnerNode(long ID) {
    assert((ID == nodeIDs.front()) || (ID == nodeIDs.back()));
    if (ID == nodeIDs.front()) {
        return nodeIDs.begin()[2];
    } else {
        return nodeIDs.begin()[1];
    }
}
#include "grid/edge.hpp"

#include <cassert>
#include <cmath>

#include "grid/mesh.hpp"

Edge::Edge(Mesh &mesh, long ID, long startNodeID, long endNodeID, bool boundEdge, int innerNodeNum) {
    this->id = ID;
    this->is_bound = boundEdge;
    end_node_ids.push_back(startNodeID);
    end_node_ids.push_back(endNodeID);

    double startNodeX, startNodeY, endNodeX, endNodeY;
    startNodeX = mesh.nodes[startNodeID].coords.x;
    startNodeY = mesh.nodes[startNodeID].coords.y;
    endNodeX = mesh.nodes[endNodeID].coords.x;
    endNodeY = mesh.nodes[endNodeID].coords.y;

    length = sqrt(pow((endNodeX - startNodeX), 2) + pow((endNodeY - startNodeY), 2));

    // In nodeIDs there will be always nodes order from start node to end node
    this->node_ids.push_back(startNodeID);
    mesh.nodes[startNodeID].edge_ids.insert(ID);
    Vector startNode(startNodeX, startNodeY);
    Vector edgeVector(endNodeX - startNodeX, endNodeY - startNodeY);
    for (int i = 1; i < innerNodeNum + 1; i++) {
        Node *node = new Node(mesh,
                              edgeVector.x / (innerNodeNum + 1) * i + startNode.x,
                              edgeVector.y / (innerNodeNum + 1) * i + startNode.y,
                              true, boundEdge);
        node->edge_ids.insert(ID);
        mesh.nodes.push_back(*node);
        this->node_ids.push_back(node->id);
        this->inner_node_ids.push_back(node->id);
    }
    this->node_ids.push_back(endNodeID);
    mesh.nodes[endNodeID].edge_ids.insert(ID);

    for (long nodeID : this->node_ids) {
        if (mesh.nodes[nodeID].used) {
            this->used_node_ids.push_back(nodeID);
        }
    }

    std::pair<int, int> tempNodeIDs;
    tempNodeIDs = std::make_pair(startNodeID, endNodeID);
    mesh.map_nodes_with_edge.insert(std::pair<std::pair<int, int>, int>(tempNodeIDs, ID));
    tempNodeIDs = std::make_pair(endNodeID, startNodeID);
    mesh.map_nodes_with_edge.insert(std::pair<std::pair<int, int>, int>(tempNodeIDs, ID));
}

long Edge::getAnotherEndNode(long ID) {
    assert((ID == node_ids.front()) || (ID == node_ids.back()));
    if (ID == node_ids.front()) {
        return node_ids.begin()[3];
    } else {
        return node_ids.begin()[0];
    }
}

long Edge::getNearInnerNode(long ID) {
    assert((ID == node_ids.front()) || (ID == node_ids.back()));
    if (ID == node_ids.front()) {
        return node_ids.begin()[1];
    } else {
        return node_ids.begin()[2];
    }
}

long Edge::getFarInnerNode(long ID) {
    assert((ID == node_ids.front()) || (ID == node_ids.back()));
    if (ID == node_ids.front()) {
        return node_ids.begin()[2];
    } else {
        return node_ids.begin()[1];
    }
}
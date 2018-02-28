#include "Cell.hpp"
#include "Node.hpp"
#include "Mesh.hpp"

#include <vector>
#include <algorithm>
#include <cassert>

double determinant(double a1, double a2, double a3, double a4) {
    return a1 * a4 - a2 * a3;
}

Cell::Cell(Mesh &mesh, long id, long node_id_1, long node_id_2, long node_id_3) {
    std::vector<double> intersectionCoords;
    double x1, x2, x3, y1, y2, y3;

    x1 = mesh.nodes[node_id_1].data.coords.x;
    x2 = mesh.nodes[node_id_2].data.coords.x;
    x3 = mesh.nodes[node_id_3].data.coords.x;
    y1 = mesh.nodes[node_id_1].data.coords.y;
    y2 = mesh.nodes[node_id_2].data.coords.y;
    y3 = mesh.nodes[node_id_3].data.coords.y;

    Node *node = new Node(mesh, (x1 + x2 +x3) / 3, (y1 + y2 + y3) / 3);
    mesh.nodes.push_back(*node);

    centerNodeID = node->ID;
    ID = id;
    nodeIDs.push_back(node_id_1);
    nodeIDs.push_back(node_id_2);
    nodeIDs.push_back(node_id_3);
    volume = determinant((x1 - x3), (x2 - x3), (y1 - y3), (y2 - y3))/2;

    std::map<std::pair<int, int>,int>::const_iterator it;
    int tempEdgeID;

    it = mesh.mapNodesWithEdge.find(std::make_pair(node_id_1, node_id_2));
    tempEdgeID = it->second;
    edgeIDs.push_back(tempEdgeID);
    mesh.edges[tempEdgeID].cellIDs.push_back(ID);
    it = mesh.mapNodesWithEdge.find(std::make_pair(node_id_2, node_id_3));
    tempEdgeID = it->second;
    edgeIDs.push_back(tempEdgeID);
    mesh.edges[tempEdgeID].cellIDs.push_back(ID);
    it = mesh.mapNodesWithEdge.find(std::make_pair(node_id_3, node_id_1));
    tempEdgeID = it->second;
    edgeIDs.push_back(tempEdgeID);
    mesh.edges[tempEdgeID].cellIDs.push_back(ID);
}

long Cell::getNextNodeID(unsigned long nodeIDPos) {
    assert(nodeIDPos < nodeIDs.size());

    if (nodeIDPos == nodeIDs.size() - 1) {
        return nodeIDs.front();
    }
    return nodeIDs[nodeIDPos + 1];
}

long Cell::getPrevNodeID(unsigned long nodeIDPos) {
    assert(nodeIDPos < nodeIDs.size());

    if (nodeIDPos == 0) {
        return nodeIDs.back();
    }
    return nodeIDs[nodeIDPos - 1];
}

std::vector<long> Cell::getEdgeOrderedNodeIDs(std::vector<long> edgeUnorderedNodeIDs) {
    assert(edgeUnorderedNodeIDs.size() == 2);

    std::vector<long> edgeOrderedNodeIDs;
    for (unsigned long i = 0; i < nodeIDs.size(); i++) {
        if (nodeIDs[i] == edgeUnorderedNodeIDs[0]) {
            if (getNextNodeID(i) == edgeUnorderedNodeIDs[1]) {
                edgeOrderedNodeIDs = edgeUnorderedNodeIDs;
                break;
            }
            if (getPrevNodeID(i) == edgeUnorderedNodeIDs[1]) {
                std::reverse(edgeUnorderedNodeIDs.begin(), edgeUnorderedNodeIDs.end());
                edgeOrderedNodeIDs = edgeUnorderedNodeIDs;
                break;
            }
        }
    }

    assert(edgeOrderedNodeIDs.size() == 2);
    return edgeOrderedNodeIDs;
}

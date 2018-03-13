#include "Cell.hpp"

#include "Node.hpp"
#include "Mesh.hpp"
#include "Vector.hpp"

#include <vector>
#include <algorithm>
#include <cassert>

double determinant(double a1, double a2, double a3, double a4) {
    return a1 * a4 - a2 * a3;
}

Cell::Cell(Mesh &mesh, long id, long node_id_1, long node_id_2, long node_id_3) {
    std::vector<double> intersectionCoords;
    double x1, x2, x3, y1, y2, y3, xMedian, yMedian;

    x1 = mesh.nodes[node_id_1].data.coords.x;
    x2 = mesh.nodes[node_id_2].data.coords.x;
    x3 = mesh.nodes[node_id_3].data.coords.x;
    y1 = mesh.nodes[node_id_1].data.coords.y;
    y2 = mesh.nodes[node_id_2].data.coords.y;
    y3 = mesh.nodes[node_id_3].data.coords.y;
    xMedian = (x1 + x2 +x3) / 3;
    yMedian = (y1 + y2 + y3) / 3;

    double median1Length = Vector::length(1.5 * (x1 - xMedian), 1.5 * (y1 - yMedian));
    double median2Length = Vector::length(1.5 * (x2 - xMedian), 1.5 * (y2 - yMedian));
    double median3Length = Vector::length(1.5 * (x3 - xMedian), 1.5 * (y3 - yMedian));
    maxH = std::max(std::max(median1Length, median2Length), median3Length);

    Node *node = new Node(mesh, xMedian, yMedian);
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
    edgeToMedianLength[tempEdgeID] = median3Length;
    mesh.edges[tempEdgeID].cellIDs.push_back(ID);
    it = mesh.mapNodesWithEdge.find(std::make_pair(node_id_2, node_id_3));
    tempEdgeID = it->second;
    edgeIDs.push_back(tempEdgeID);
    edgeToMedianLength[tempEdgeID] = median1Length;
    mesh.edges[tempEdgeID].cellIDs.push_back(ID);
    it = mesh.mapNodesWithEdge.find(std::make_pair(node_id_3, node_id_1));
    tempEdgeID = it->second;
    edgeIDs.push_back(tempEdgeID);
    edgeToMedianLength[tempEdgeID] = median2Length;
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

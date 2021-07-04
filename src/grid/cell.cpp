#include "grid/cell.hpp"

#include <algorithm>
#include <cassert>
#include <vector>
#include <limits>

#include "grid/vector.hpp"
#include "grid/mesh.hpp"

double calculateDeterminant(double a11, double a12, double a21, double a22) {
    return std::abs(a11 * a22 - a12 * a21);
}

void Cell::assign_opposite_node_ids(Mesh &mesh) {
    std::map<unsigned long, std::vector<double>> edgeToLineCoef;
    for (unsigned long edgeID : edge_ids) {
        Edge *edge = &mesh.edges[edgeID];
        std::vector<long> orderedEndNodes = get_edge_ordered_node_ids(edge->end_node_ids);
        Node *endNode1 = &mesh.nodes[orderedEndNodes[0]];
        Node *endNode2 = &mesh.nodes[orderedEndNodes[1]];

        // let's build edge line like Ax+By+C=0
        double A = endNode1->coords.y - endNode2->coords.y;
        double B = endNode2->coords.x - endNode1->coords.x;
        double C = endNode1->coords.x * endNode2->coords.y -
            endNode1->coords.y * endNode2->coords.x;

        std::vector<double> coef;
        coef.push_back(A);
        coef.push_back(B);
        coef.push_back(C);
        edgeToLineCoef[edgeID] = coef;
    }

    Node *cellNode = &mesh.nodes[center_node_id];
    for (unsigned long edgeID : edge_ids) {
        Edge *edge = &mesh.edges[edgeID];
        for (unsigned long nodeID : edge->used_node_ids) {
            Node *edgeNode = &mesh.nodes[nodeID];
            // let's build line passing edgeNode and cellNode like Ax+By+C=0
            double lineA = edgeNode->coords.y - cellNode->coords.y;
            double lineB = cellNode->coords.x - edgeNode->coords.x;
            double lineC = edgeNode->coords.x * cellNode->coords.y -
                edgeNode->coords.y * cellNode->coords.x;

            unsigned long oppositeEdgeID = std::numeric_limits<unsigned long>::max();
            double oppositeX = std::numeric_limits<double>::max();
            double oppositeY = std::numeric_limits<double>::max();
            double minCosBetweenLines = 1;
            std::map<unsigned long, std::vector<double>>::iterator it;
            for (it = edgeToLineCoef.begin(); it != edgeToLineCoef.end(); it++) {
                if (it->first == edgeID) {
                    continue;
                }

                double edgeLineA = it->second[0];
                double edgeLineB = it->second[1];
                double edgeLineC = it->second[2];

                Vector nodeLine(lineA, lineB);
                Vector edgeLine(edgeLineA, edgeLineB);
                double dot = nodeLine.x * edgeLine.x + nodeLine.y * edgeLine.y;
                double cosBetweenLines = std::abs((dot) / (nodeLine.length() * edgeLine.length()));
                if (cosBetweenLines > minCosBetweenLines) {
                    continue;
                }

                double determinant = calculateDeterminant(edgeLineA, edgeLineB, lineA, lineB);
                if (determinant == 0) {
                    continue;
                }

                double xCross = (edgeLineC * lineB - lineC * edgeLineB) / determinant;
                double yCross = (lineC * edgeLineA - edgeLineC * lineA) / determinant;

                Edge *oppositeEdge = &mesh.edges[it->first];
                Node *endNode1 = &mesh.nodes[oppositeEdge->end_node_ids[0]];
                Node *endNode2 = &mesh.nodes[oppositeEdge->end_node_ids[1]];
                double maxEndX = std::max(endNode1->coords.x, endNode2->coords.x);
                double minEndX = std::min(endNode1->coords.x, endNode2->coords.x);
                double maxEndY = std::max(endNode1->coords.y, endNode2->coords.y);
                double minEndY = std::min(endNode1->coords.y, endNode2->coords.y);
                bool minAndMaxXEquals = std::abs(minEndX - maxEndX) <= 0.1e-9;
                bool minAndMaxYEquals = std::abs(minEndY - maxEndY) <= 0.1e-9;
                if (minAndMaxXEquals && std::abs(minEndX - xCross) > 0.1e-9) {
                    continue;
                }
                if (minAndMaxYEquals && std::abs(minEndY - yCross) > 0.1e-9) {
                    continue;
                }
                if ((!minAndMaxXEquals && ((minEndX > xCross) || (maxEndX < xCross))) ||
                    (!minAndMaxYEquals && ((minEndY > yCross) || (maxEndY < yCross)))) {
                    continue;
                }

                Vector cross(xCross, yCross);
                double distFromNodeToCellNode = Vector::length(
                    cellNode->coords.x - edgeNode->coords.x,
                    cellNode->coords.y - edgeNode->coords.y
                );
                double distFromCellNodeToCross = Vector::length(
                    cellNode->coords.x - cross.x,
                    cellNode->coords.y - cross.y
                );
                double ratio = distFromCellNodeToCross / distFromNodeToCellNode;
                if ((ratio < 0.95) || (ratio > 1.05)) {
                    continue;
                }

                oppositeEdgeID = it->first;
                oppositeX = xCross;
                oppositeY = yCross;
                minCosBetweenLines = cosBetweenLines;
            }

            if (oppositeEdgeID == std::numeric_limits<unsigned long>::max()) {
                node_id_to_opposite_node_id[nodeID] = -1;
            } else {
                Edge *oppositeEdge = &mesh.edges[oppositeEdgeID];
                unsigned long finalOppositeNodeID = std::numeric_limits<unsigned long>::max();
                double minDist = std::numeric_limits<double>::max();
                Vector opposite(oppositeX, oppositeY);
                for (unsigned long oppositeNodeID : oppositeEdge->node_ids) {
                    Node *oppositeNode = &mesh.nodes[oppositeNodeID];
                    if (!oppositeNode->used) {
                        continue;
                    }

                    double dist = Vector::length(
                        oppositeNode->coords.x - opposite.x,
                        oppositeNode->coords.y - opposite.y
                    );
                    if (dist < minDist) {
                        finalOppositeNodeID = oppositeNodeID;
                        minDist = dist;
                    }
                }
                node_id_to_opposite_node_id[nodeID] = finalOppositeNodeID;
            }
        }

    }
}

Cell::Cell(Mesh &mesh, long id, long node_id_1, long node_id_2, long node_id_3) {
    this->id = id;

    std::vector<double> intersectionCoords;
    double x1, x2, x3, y1, y2, y3, xMedian, yMedian;

    Node *node1 = &mesh.nodes[node_id_1];
    Node *node2 = &mesh.nodes[node_id_2];
    Node *node3 = &mesh.nodes[node_id_3];

    node1->cell_ids.insert(id);
    node2->cell_ids.insert(id);
    node3->cell_ids.insert(id);

    x1 = node1->coords.x;
    x2 = node2->coords.x;
    x3 = node3->coords.x;
    y1 = node1->coords.y;
    y2 = node2->coords.y;
    y3 = node3->coords.y;

    xMedian = (x1 + x2 + x3) / 3;
    yMedian = (y1 + y2 + y3) / 3;

    Node *node = new Node(mesh, xMedian, yMedian, true, false);
    mesh.nodes.push_back(*node);

    this->center_node_id = node->id;
    node_ids.push_back(node_id_1);
    node_ids.push_back(node_id_2);
    node_ids.push_back(node_id_3);
    volume = calculateDeterminant((x1 - x3), (x2 - x3), (y1 - y3), (y2 - y3)) / 2;

    std::map<std::pair<int, int>, int>::const_iterator it;
    int tempEdgeID;

    it = mesh.map_nodes_with_edge.find(std::make_pair(node_id_1, node_id_2));
    tempEdgeID = it->second;
    edge_ids.push_back(tempEdgeID);
    mesh.edges[tempEdgeID].cell_ids.push_back(id);

    it = mesh.map_nodes_with_edge.find(std::make_pair(node_id_2, node_id_3));
    tempEdgeID = it->second;
    edge_ids.push_back(tempEdgeID);
    mesh.edges[tempEdgeID].cell_ids.push_back(id);

    it = mesh.map_nodes_with_edge.find(std::make_pair(node_id_3, node_id_1));
    tempEdgeID = it->second;
    edge_ids.push_back(tempEdgeID);
    mesh.edges[tempEdgeID].cell_ids.push_back(id);

    assign_opposite_node_ids(mesh);
}

long Cell::get_next_node_id(unsigned long nodeIDPos) {
    assert(nodeIDPos < node_ids.size());

    if (nodeIDPos == node_ids.size() - 1) {
        return node_ids.front();
    }
    return node_ids[nodeIDPos + 1];
}

long Cell::get_prev_node_id(unsigned long nodeIDPos) {
    assert(nodeIDPos < node_ids.size());

    if (nodeIDPos == 0) {
        return node_ids.back();
    }
    return node_ids[nodeIDPos - 1];
}

std::vector<long> Cell::get_edge_ordered_node_ids(std::vector<long> edgeUnorderedNodeIDs) {
    assert(edgeUnorderedNodeIDs.size() == 2);

    std::vector<long> edgeOrderedNodeIDs;
    for (unsigned long i = 0; i < node_ids.size(); i++) {
        if (node_ids[i] == edgeUnorderedNodeIDs[0]) {
            if (get_next_node_id(i) == edgeUnorderedNodeIDs[1]) {
                edgeOrderedNodeIDs = edgeUnorderedNodeIDs;
                break;
            }
            if (get_prev_node_id(i) == edgeUnorderedNodeIDs[1]) {
                std::reverse(edgeUnorderedNodeIDs.begin(), edgeUnorderedNodeIDs.end());
                edgeOrderedNodeIDs = edgeUnorderedNodeIDs;
                break;
            }
        }
    }

    assert(edgeOrderedNodeIDs.size() == 2);
    return edgeOrderedNodeIDs;
}

#include "Cell.hpp"

#include "Mesh.hpp"
#include "Node.hpp"
#include "Vector.hpp"

#include <algorithm>
#include <cassert>
#include <vector>
#include <limits>
#include <iostream>

double calculateDeterminant(double a11, double a12, double a21, double a22) {
    return a11 * a22 - a12 * a21;
}

void Cell::assignOppositeNodeIDs(Mesh &mesh) {
    std::map<unsigned long, std::vector<double>> edgeToLineCoef;
    for (unsigned long edgeID : edgeIDs) {
        Edge *edge = &mesh.edges[edgeID];
        std::vector<long> orderedEndNodes = getEdgeOrderedNodeIDs(edge->edgeEndsNodeIDs);
        Node *endNode1 = &mesh.nodes[orderedEndNodes[0]];
        Node *endNode2 = &mesh.nodes[orderedEndNodes[1]];

        // let's build edge line like Ax+By+C=0
        double A = endNode1->data.coords.y - endNode2->data.coords.y;
        double B = endNode2->data.coords.x - endNode1->data.coords.x;
        double C = endNode1->data.coords.x * endNode2->data.coords.y -
                   endNode1->data.coords.y * endNode2->data.coords.x;

        std::vector<double> coef;
        coef.push_back(A);
        coef.push_back(B);
        coef.push_back(C);
        edgeToLineCoef[edgeID] = coef;
    }

    Node *cellNode = &mesh.nodes[centerNodeID];
    for (unsigned long edgeID : edgeIDs) {
        Edge *edge = &mesh.edges[edgeID];
        for (unsigned long nodeID : edge->getUsedNodes(mesh)) {
            // std::cout << "------------------ " << ID << std::endl;
            Node *edgeNode = &mesh.nodes[nodeID];
            // let's build line passing edgeNode and cellNode like Ax+By+C=0
            double lineA = edgeNode->data.coords.y - cellNode->data.coords.y;
            double lineB = cellNode->data.coords.x - edgeNode->data.coords.x;
            double lineC = edgeNode->data.coords.x * cellNode->data.coords.y -
                           edgeNode->data.coords.y * cellNode->data.coords.x;
            
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
                double cosBetweenLines = (nodeLine * edgeLine) / (nodeLine.length() * edgeLine.length());
                // std::cout << "cos " << cosBetweenLines << std::endl;
                if (std::abs(cosBetweenLines) > minCosBetweenLines) {
                    continue;
                }

                double determinant = std::abs(calculateDeterminant(edgeLineA, edgeLineB, lineA, lineB));
                // std::cout << "det " << determinant << std::endl;
                if (determinant == 0) {
                    continue;
                }

                double xCross = (edgeLineC * lineB - lineC * edgeLineB) / determinant;
                double yCross = (lineC * edgeLineA - edgeLineC * lineA) / determinant;
                
                Edge *oppositeEdge = &mesh.edges[it->first];
                Node *endNode1 = &mesh.nodes[oppositeEdge->edgeEndsNodeIDs[0]];
                Node *endNode2 = &mesh.nodes[oppositeEdge->edgeEndsNodeIDs[1]];
                double maxEndX = std::max(endNode1->data.coords.x, endNode2->data.coords.x);
                double minEndX = std::min(endNode1->data.coords.x, endNode2->data.coords.x);
                double maxEndY = std::max(endNode1->data.coords.y, endNode2->data.coords.y);
                double minEndY = std::min(endNode1->data.coords.y, endNode2->data.coords.y);
                if (std::abs(minEndX - maxEndX) <= 0.1e-9 && std::abs(minEndX - xCross) > 0.1e-9) {
                    continue;
                }
                if (std::abs(minEndY - maxEndY) <= 0.1e-9 && std::abs(minEndY - yCross) > 0.1e-9) {
                    continue;
                }
                // std::cout << it->first << " " << minEndX << " " << maxEndX << " " << minEndY << " " << maxEndY << " " << xCross << " " << yCross << std::endl;
                if ((std::abs(minEndX - maxEndX) > 0.1e-9 && ((minEndX > xCross) || (maxEndX < xCross))) || 
                    (std::abs(minEndY - maxEndY) > 0.1e-9 && ((minEndY > yCross) || (maxEndY < yCross)))) {
                    continue;
                }
                
                Vector cross(xCross, yCross);
                double distFromNodeToCellNode = (cellNode->data.coords - edgeNode->data.coords).length();
                double distFromCellNodeToCross = (cellNode->data.coords - cross).length();
                double ratio = distFromCellNodeToCross / distFromNodeToCellNode;
                // std::cout << "ratio = " <<ratio << std::endl;
                if ((ratio < 0.95) || (ratio > 1.05)) {
                    continue;
                }

                oppositeEdgeID = it->first;
                oppositeX = xCross;
                oppositeY = yCross;
                minCosBetweenLines = cosBetweenLines;
            }

            if (oppositeEdgeID == std::numeric_limits<unsigned long>::max()) {
                // std::cout << ID << " " << edgeID << " " <<edgeNode->data.coords.x << " " << edgeNode->data.coords.y << " " << cellNode->data.coords.x << " " << cellNode->data.coords.y << std::endl;
                nodeIDToOppositeNodeID[nodeID] = -1;
            } else {
                Edge *oppositeEdge = &mesh.edges[oppositeEdgeID];
                unsigned long finalOppositeNodeID = std::numeric_limits<unsigned long>::max();
                double minDist = std::numeric_limits<double>::max();
                Vector opposite(oppositeX, oppositeY);
                for (unsigned long oppositeNodeID : oppositeEdge->nodeIDs) {
                    Node *oppositeNode = &mesh.nodes[oppositeNodeID];
                    if (!oppositeNode->used) {
                        continue;
                    }

                    double dist = (oppositeNode->data.coords - opposite).length();
                    if (dist < minDist) {
                        finalOppositeNodeID = oppositeNodeID;
                        minDist = dist;
                    }
                }
                nodeIDToOppositeNodeID[nodeID] = finalOppositeNodeID;
            }
        }
        
    }
}

Cell::Cell(Mesh &mesh, long ID, long nodeID1, long nodeID2, long nodeID3) {
    this->ID = ID;

    std::vector<double> intersectionCoords;
    double x1, x2, x3, y1, y2, y3, xMedian, yMedian;

    x1 = mesh.nodes[nodeID1].data.coords.x;
    x2 = mesh.nodes[nodeID2].data.coords.x;
    x3 = mesh.nodes[nodeID3].data.coords.x;
    y1 = mesh.nodes[nodeID1].data.coords.y;
    y2 = mesh.nodes[nodeID2].data.coords.y;
    y3 = mesh.nodes[nodeID3].data.coords.y;

    xMedian = (x1 + x2 + x3) / 3;
    yMedian = (y1 + y2 + y3) / 3;

    double median1Length = Vector::length(1.5 * (x1 - xMedian), 1.5 * (y1 - yMedian));
    double median2Length = Vector::length(1.5 * (x2 - xMedian), 1.5 * (y2 - yMedian));
    double median3Length = Vector::length(1.5 * (x3 - xMedian), 1.5 * (y3 - yMedian));
    maxH = std::max(std::max(median1Length, median2Length), median3Length);

    Node *node = new Node(mesh, xMedian, yMedian, true);
    mesh.nodes.push_back(*node);

    centerNodeID = node->ID;
    nodeIDs.push_back(nodeID1);
    nodeIDs.push_back(nodeID2);
    nodeIDs.push_back(nodeID3);
    volume = calculateDeterminant((x1 - x3), (x2 - x3), (y1 - y3), (y2 - y3)) / 2;

    std::map<std::pair<int, int>, int>::const_iterator it;
    int tempEdgeID;

    it = mesh.mapNodesWithEdge.find(std::make_pair(nodeID1, nodeID2));
    tempEdgeID = it->second;
    edgeIDs.push_back(tempEdgeID);
    edgeToMedianLength[tempEdgeID] = median3Length;
    mesh.edges[tempEdgeID].cellIDs.push_back(ID);

    it = mesh.mapNodesWithEdge.find(std::make_pair(nodeID2, nodeID3));
    tempEdgeID = it->second;
    edgeIDs.push_back(tempEdgeID);
    edgeToMedianLength[tempEdgeID] = median1Length;
    mesh.edges[tempEdgeID].cellIDs.push_back(ID);

    it = mesh.mapNodesWithEdge.find(std::make_pair(nodeID3, nodeID1));
    tempEdgeID = it->second;
    edgeIDs.push_back(tempEdgeID);
    edgeToMedianLength[tempEdgeID] = median2Length;
    mesh.edges[tempEdgeID].cellIDs.push_back(ID);

    assignOppositeNodeIDs(mesh);
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

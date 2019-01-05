#include "Cell.hpp"

#include "Mesh.hpp"
#include "Node.hpp"
#include "Vector.hpp"
#include "Parameters.hpp"

#include <algorithm>
#include <cassert>
#include <vector>
#include <limits>

double calculateDeterminant(double a11, double a12, double a21, double a22) {
    return std::abs(a11 * a22 - a12 * a21);
}

void Cell::assignOppositeNodeIDs(Mesh &mesh) {
    std::map<unsigned long, std::vector<double>> edgeToLineCoef;
    for (unsigned long edgeID : edgeIDs) {
        Edge *edge = &mesh.edges[edgeID];
        std::vector<long> orderedEndNodes = getEdgeOrderedNodeIDs(edge->endNodeIDs);
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
                double cosBetweenLines = std::abs((nodeLine * edgeLine) / (nodeLine.length() * edgeLine.length()));
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
                Node *endNode1 = &mesh.nodes[oppositeEdge->endNodeIDs[0]];
                Node *endNode2 = &mesh.nodes[oppositeEdge->endNodeIDs[1]];
                double maxEndX = std::max(endNode1->data.coords.x, endNode2->data.coords.x);
                double minEndX = std::min(endNode1->data.coords.x, endNode2->data.coords.x);
                double maxEndY = std::max(endNode1->data.coords.y, endNode2->data.coords.y);
                double minEndY = std::min(endNode1->data.coords.y, endNode2->data.coords.y);
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
                double distFromNodeToCellNode = (cellNode->data.coords - edgeNode->data.coords).length();
                double distFromCellNodeToCross = (cellNode->data.coords - cross).length();
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

void assignRowValues(arma::mat &mat, double x, double y, int row) {
    mat(row, 0) = 1;
    mat(row, 1) = x;
    mat(row, 2) = y;
    mat(row, 3) = x * x;
    mat(row, 4) = y * y;
    mat(row, 5) = x * y;
    mat(row, 6) = mat(row, 3) * x;
    mat(row, 7) = mat(row, 4) * y;
    mat(row, 8) = mat(row, 5) * x;
    mat(row, 9) = mat(row, 5) * y;
}

void Cell::buildInterpolationMat(Mesh &mesh) {
    int matSize = Parameters::EDGE_INNER_NODES_NUMBER * 3 + 4;
    if (matSize != 10) {
        std::cout << matSize << std::endl;
        return;
    }

    arma::mat interpolationMat = arma::mat(matSize, matSize);
    int row = 0;

    for (long nodeID : nodeIDs) {
        Vector coords = mesh.nodes[nodeID].data.coords;
        assignRowValues(interpolationMat, coords.x, coords.y, row);
        row++;
    }
    for (long edgeID : edgeIDs) {
        for (long nodeID : mesh.edges[edgeID].getInnerNodes()) {
            Vector coords = mesh.nodes[nodeID].data.coords;
            assignRowValues(interpolationMat, coords.x, coords.y, row);
            row++;
        }
    }
    Vector coords = mesh.nodes[centerNodeID].data.coords;
    assignRowValues(interpolationMat, coords.x, coords.y, row);
    this->interpolationMat = interpolationMat;
}

Cell::Cell(Mesh &mesh, long ID, long nodeID1, long nodeID2, long nodeID3) {
    this->ID = ID;

    std::vector<double> intersectionCoords;
    double x1, x2, x3, y1, y2, y3, xMedian, yMedian;

    Node *node1 = &mesh.nodes[nodeID1];
    Node *node2 = &mesh.nodes[nodeID2];
    Node *node3 = &mesh.nodes[nodeID3];

    node1->cellIDs.insert(ID);
    node2->cellIDs.insert(ID);
    node3->cellIDs.insert(ID);

    x1 = node1->data.coords.x;
    x2 = node2->data.coords.x;
    x3 = node3->data.coords.x;
    y1 = node1->data.coords.y;
    y2 = node2->data.coords.y;
    y3 = node3->data.coords.y;

    xMedian = (x1 + x2 + x3) / 3;
    yMedian = (y1 + y2 + y3) / 3;

    Node *node = new Node(mesh, xMedian, yMedian, true, false, true);
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
    mesh.edges[tempEdgeID].cellIDs.push_back(ID);

    it = mesh.mapNodesWithEdge.find(std::make_pair(nodeID2, nodeID3));
    tempEdgeID = it->second;
    edgeIDs.push_back(tempEdgeID);
    mesh.edges[tempEdgeID].cellIDs.push_back(ID);

    it = mesh.mapNodesWithEdge.find(std::make_pair(nodeID3, nodeID1));
    tempEdgeID = it->second;
    edgeIDs.push_back(tempEdgeID);
    mesh.edges[tempEdgeID].cellIDs.push_back(ID);

    assignOppositeNodeIDs(mesh);
    buildInterpolationMat(mesh);
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

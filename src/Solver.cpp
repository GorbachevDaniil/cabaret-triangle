#include "Solver.hpp"

#include "Cell.hpp"
#include "Edge.hpp"
#include "Node.hpp"
#include "Vector.hpp"

#include <limits>
#include <cassert>
#include <armadillo>
#include <iomanip>

double Solver::calculateTau() {
    double tau = std::numeric_limits<double>::max();
    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];
        Node *cellNode = &mesh->nodes[cell->centerNodeID];
        for (unsigned long edgeID : cell->edgeIDs) {
            for (unsigned long nodeID : mesh->edges[edgeID].getUsedNodes(*mesh)) {
                Node *node = &mesh->nodes[nodeID];
                double lambda = node->data.u * cell->nodeToTransferVector[nodeID];
                double h = 2 * (node->data.coords - cellNode->data.coords).length();
                double nodeTau = h / std::abs(lambda);
                if (tau > nodeTau) {
                    tau = nodeTau;
                }
            }
        }
    }
    assert(tau != std::numeric_limits<double>::max());
    return cfl * tau;
}

double Solver::calculateDivOnEdge(Edge *edge, int phase) {
    std::vector<long> usedNodeIDs = edge->getUsedNodes(*mesh);
    double div = 0;
    switch (usedNodeIDs.size()) {
        case 1: {
            Node *edgeNode = &mesh->nodes[usedNodeIDs[0]];

            double phi = phase == 1 ? edgeNode->data.phi0 : edgeNode->data.phi2;
            div = phi * (edgeNode->data.u * edge->normal) * edge->length;
            break;
        }
        case 2: {
            Node *edgeNode1 = &mesh->nodes[usedNodeIDs[0]];
            Node *edgeNode2 = &mesh->nodes[usedNodeIDs[1]];

            double phi1 = phase == 1 ? edgeNode1->data.phi0 : edgeNode1->data.phi2;
            double phi2 = phase == 1 ? edgeNode2->data.phi0 : edgeNode2->data.phi2;
            double phi = (phi1 + phi2) / 2;
            div = phi * (edgeNode1->data.u * edge->normal) * edge->length;
            break;
        }
        case 4: {
            Node *edgeNode1 = &mesh->nodes[usedNodeIDs[0]];
            Node *edgeNode2 = &mesh->nodes[usedNodeIDs[1]];
            Node *edgeNode3 = &mesh->nodes[usedNodeIDs[2]];
            Node *edgeNode4 = &mesh->nodes[usedNodeIDs[3]];

            double phi1 = phase == 1 ? edgeNode1->data.phi0 : edgeNode1->data.phi2;
            double phi2 = phase == 1 ? edgeNode2->data.phi0 : edgeNode2->data.phi2;
            double phi3 = phase == 1 ? edgeNode3->data.phi0 : edgeNode3->data.phi2;
            double phi4 = phase == 1 ? edgeNode4->data.phi0 : edgeNode4->data.phi2;

            double div1 = phi1 * (edgeNode1->data.u * edge->normal);
            double div2 = phi2 * (edgeNode2->data.u * edge->normal);
            double div3 = phi3 * (edgeNode3->data.u * edge->normal);
            double div4 = phi4 * (edgeNode4->data.u * edge->normal);

            div = (div1 + 3 * div2 + 3 * div3 + div4) / 8 * edge->length;
            break;
        }
        default:
            assert(false);
            break;
    }

    return div;
}

void Solver::processPhase1(double tau) {
    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];
        Node *cellCenter = &mesh->nodes[cell->centerNodeID];
        double div = 0;

        for (long edgeID : cell->edgeIDs) {
            div += calculateDivOnEdge(&mesh->edges[edgeID], 1) * cell->edgeToNormalDir[edgeID];
        }
        div /= cell->volume;
        cellCenter->data.phi1 = cellCenter->data.phi0 - tau * div / 2;
    }
}

double getNewInvariantValue(Data *data, Data *centerData, Data *oppositeData, 
                            Vector transfer, double tau, bool monotize) {
    double phiOpposite0 = oppositeData->phi0;
    double phiCenter0 = centerData->phi0;
    double phi0 = data->phi0;

    double phiCenter1 = centerData->phi1;

    double phi2 = 2 * phiCenter1 - phiOpposite0;

    if (monotize) {
        double h = (oppositeData->coords - data->coords).length();
        double Q = (phiCenter1 - phiCenter0) / tau / 2 +
                    (centerData->u * transfer) * (phi0 - phiOpposite0) / h;
        double min = std::min(std::min(phi0, phiOpposite0), phiCenter0) + tau * Q;
        double max = std::max(std::max(phi0, phiOpposite0), phiCenter0) + tau * Q;

        phi2 = std::max(phi2, min);
        phi2 = std::min(phi2, max);
    }

    return phi2;
}

void processPhase2BoundEdge(Mesh *mesh, Edge *edge, double tau) {
    Cell *cell = &mesh->cells[edge->cellIDs[0]];
    Data *centerData = &mesh->nodes[cell->centerNodeID].data;

    for (unsigned long nodeID : edge->getInnerNodes()) {
        Node *node = &mesh->nodes[nodeID];
        node->phase2Calculated = true;
        Data *data = &node->data;
        if (centerData->u * cell->nodeToTransferVector[nodeID] < 0) {
            data->phi2 = 0;
            continue;
        }

        assert(cell->nodeIDToOppositeNodeID[nodeID] != -1);
        Data *oppositeData = &mesh->nodes[cell->nodeIDToOppositeNodeID[nodeID]].data;

        data->phi2 = getNewInvariantValue(data, centerData, oppositeData, 
                                          cell->nodeToTransferVector[nodeID], tau, true);
    }
}

void processPhase2InnerEdge(Mesh *mesh, Edge *edge, double tau) {
    Cell *cell1 = &mesh->cells[edge->cellIDs[0]];
    Cell *cell2 = &mesh->cells[edge->cellIDs[1]];

    Data *centerData1 = &mesh->nodes[cell1->centerNodeID].data;
    Data *centerData2 = &mesh->nodes[cell2->centerNodeID].data;

    Vector normal(edge->normal * cell1->edgeToNormalDir[edge->ID]);
    Vector normalizedAvgU((centerData1->u + centerData2->u) / 2);
    normalizedAvgU = normalizedAvgU / normalizedAvgU.length();

    // normalized value of velocity used to not depend on initial value of velocity
    double cos = normal * normalizedAvgU;
    for (unsigned long nodeID : edge->getInnerNodes()) {
        Node *node = &mesh->nodes[nodeID];
        node->phase2Calculated = true;
        Data *data = &node->data;
        if ((cos > -0.1e-9) && (cos < 0.1e-9)) {
            assert(cell1->nodeIDToOppositeNodeID[nodeID] != -1);
            assert(cell2->nodeIDToOppositeNodeID[nodeID] != -1);
            Data *oppositeData1 = &mesh->nodes[cell1->nodeIDToOppositeNodeID[nodeID]].data;
            Data *oppositeData2 = &mesh->nodes[cell2->nodeIDToOppositeNodeID[nodeID]].data;
            double phi21 = getNewInvariantValue(data, centerData1, oppositeData1, 
                                                cell1->nodeToTransferVector[nodeID], tau, false);
            double phi22 = getNewInvariantValue(data, centerData2, oppositeData2, 
                                                cell2->nodeToTransferVector[nodeID], tau, false);
            data->phi2 = (phi21 + phi22) / 2;
            continue;
        }
        if (cos > 0) {
            assert(cell1->nodeIDToOppositeNodeID[nodeID] != -1);
            Data *oppositeData1 = &mesh->nodes[cell1->nodeIDToOppositeNodeID[nodeID]].data;
            data->phi2 = getNewInvariantValue(data, centerData1, oppositeData1, 
                                              cell1->nodeToTransferVector[nodeID], tau, true);
        } else {
            assert(cell2->nodeIDToOppositeNodeID[nodeID] != -1);
            Data *oppositeData2 = &mesh->nodes[cell2->nodeIDToOppositeNodeID[nodeID]].data;
            data->phi2 = getNewInvariantValue(data, centerData2, oppositeData2, 
                                              cell2->nodeToTransferVector[nodeID], tau, true);
        }
    }
}

void Solver::processPhase2(double tau) {
    for (unsigned long i = 0; i < mesh->edges.size(); i++) {
        Edge *edge = &mesh->edges[i];
        if (edge->boundEdge) {
            processPhase2BoundEdge(mesh, edge, tau);
        } else {
            processPhase2InnerEdge(mesh, edge, tau);
        }
    }

    for (unsigned long i = 0; i < mesh->nodes.size(); i++) {
        Node *node = &mesh->nodes[i];
        Data *data = &node->data;
        if (node->phase2Calculated || node->cellCenterNode || !node->used) {
            continue;
        }
        if (node->boundNode) {
            // TODO proper bound condition
            data->phi2 = 0;
            continue;
        }

        double uX = 0;
        double uY = 0;
        for (unsigned long edgeID : node->edgeIDs) {
            Edge *edge = &mesh->edges[edgeID];
            long neighborInnerNodeID = edge->getNearInnerNode(node->ID);
            Vector u = mesh->nodes[neighborInnerNodeID].data.u;
            uX += u.x;
            uY += u.y;
        }
        Vector avgU(uX, uY);
        avgU = avgU / node->edgeIDs.size();
        Vector normalizedAvgU = avgU / avgU.length();

        // --- Calculating new value with transfer vectors from cells and edges ---
        // double maxCosEdge = -1;
        // long edgeIDWithMaxCos = -1;
        // for (unsigned long edgeID : node->edgeIDs) {
        //     Edge *edge = &mesh->edges[edgeID];
        //     Data *oppositeData;
        //     if (edge->nodeIDs.back() == node->ID) {
        //         oppositeData = &mesh->nodes[edge->nodeIDs.front()].data;
        //     } else {
        //         oppositeData = &mesh->nodes[edge->nodeIDs.back()].data;
        //     }
        //     Vector transfer = (data->coords - oppositeData->coords);
        //     transfer = transfer / transfer.length();
        //     double cos = transfer * normalizedAvgU;
        //     assert(cos >= -1 && cos <= 1);
        //     if (cos >= maxCosEdge) {
        //         maxCosEdge = cos;
        //         edgeIDWithMaxCos = edgeID;
        //     }
        // }
        // assert(edgeIDWithMaxCos != -1);

        // double maxCosCell = -1;
        // long cellIDWithMaxCos = -1;
        // for (unsigned long cellID : node->cellIDs) {
        //     Vector transfer = mesh->cells[cellID].nodeToTransferVector[node->ID];
        //     double cos = transfer * normalizedAvgU;
        //     assert(cos >= -1 && cos <= 1);
        //     if (cos >= maxCosCell) {
        //         maxCosCell = cos;
        //         cellIDWithMaxCos = cellID;
        //     }
        // }
        // assert(cellIDWithMaxCos != -1);

        // if (maxCosEdge > maxCosCell) {
        //     Edge *edge = &mesh->edges[edgeIDWithMaxCos];
        //     Data *centerData = &mesh->nodes[edge->getNearInnerNode(node->ID)].data;
        //     Data *oppositeData = &mesh->nodes[edge->getFarInnerNode(node->ID)].data;

        //     double phiOpposite2 = oppositeData->phi2;
        //     double phiCenter2 = centerData->phi2;

        //     double phi2 = 2 * phiCenter2 - phiOpposite2;

        //     double phiOpposite0 = oppositeData->phi0;
        //     double phiCenter0 = centerData->phi0;
        //     double phi0 = data->phi0;
        //     Vector transfer = (data->coords - oppositeData->coords);
        //     transfer = transfer / transfer.length();
        //     double h = (oppositeData->coords - data->coords).length();
        //     double Q = (phiCenter2 - phiCenter0) / tau +
        //                 (centerData->u * transfer) *
        //                 (phi0 - phiOpposite0) / h;
        //     double min = std::min(std::min(phi0, phiOpposite0), phiCenter0) + tau * Q;
        //     double max = std::max(std::max(phi0, phiOpposite0), phiCenter0) + tau * Q;

        //     phi2 = std::max(phi2, min);
        //     phi2 = std::min(phi2, max);

        //     data->phi2 = phi2;
        //     node->phase2Calculated = true;
        // } else {
        //     Cell *cell = &mesh->cells[cellIDWithMaxCos];
        //     Data *centerData = &mesh->nodes[cell->centerNodeID].data;

        //     arma::vec phi0(10);
        //     int pos = 0;
        //     for (long nodeID : cell->nodeIDs) {
        //         phi0[pos] = mesh->nodes[nodeID].data.phi0;
        //         pos++;
        //     }
        //     for (long edgeID : cell->edgeIDs) {
        //         for (long nodeID : mesh->edges[edgeID].getInnerNodes()) {
        //             phi0[pos] = mesh->nodes[nodeID].data.phi0;
        //             pos++;
        //         }
        //     }
        //     phi0[pos] = centerData->phi0;
        //     arma::vec a = arma::solve(cell->interpolationMat, phi0);
            
        //     double coef = -(avgU * cell->nodeToTransferVector[node->ID]) * tau;
        //     Vector intersectCoords = cell->nodeToTransferVector[node->ID] * coef;
        //     intersectCoords = intersectCoords + data->coords;
        //     double intersectX = intersectCoords.x;
        //     double intersectY = intersectCoords.y;

        //     double phi2 = 0;
        //     phi2 += a[0];
        //     phi2 += a[1] * intersectX;
        //     phi2 += a[2] * intersectY;
        //     phi2 += a[3] * intersectX * intersectX;
        //     phi2 += a[4] * intersectY * intersectY;
        //     phi2 += a[5] * intersectX * intersectY;
        //     phi2 += a[6] * intersectX * intersectX * intersectX;
        //     phi2 += a[7] * intersectY * intersectY * intersectY;
        //     phi2 += a[8] * intersectX * intersectY * intersectX;
        //     phi2 += a[9] * intersectX * intersectY * intersectY;

        //     double phiLeft0 = centerData->phi0;
        //     double phiCenter0 = centerData->phi0;
        //     double phiRight0 = data->phi0;

        //     double phiCenter1 = centerData->phi1;

        //     double h = (centerData->coords - data->coords).length();
        //     double Q = (phiCenter1 - phiCenter0) / tau / 2 +
        //                 (((data->u + centerData->u) / 2) * cell->nodeToTransferVector[node->ID]) *
        //                 (phiRight0 - phiLeft0) / h;
        //     double min = std::min(std::min(phiRight0, phiLeft0), phiCenter0) + tau * Q;
        //     double max = std::max(std::max(phiRight0, phiLeft0), phiCenter0) + tau * Q;

        //     phi2 = std::max(phi2, min);
        //     phi2 = std::min(phi2, max);

        //     data->phi2 = phi2;
        //     node->phase2Calculated = true;
        // }
        // --- ---

        // --- Calculating new value with transfer vectors from edges ---
        // double maxCosEdge = -1;
        // long edgeIDWithMaxCos = -1;
        // for (unsigned long edgeID : node->edgeIDs) {
        //     Edge *edge = &mesh->edges[edgeID];
        //     Data *oppositeData;
        //     if (edge->nodeIDs.back() == node->ID) {
        //         oppositeData = &mesh->nodes[edge->nodeIDs.front()].data;
        //     } else {
        //         oppositeData = &mesh->nodes[edge->nodeIDs.back()].data;
        //     }
        //     Vector transfer = (data->coords - oppositeData->coords);
        //     transfer = transfer / transfer.length();
        //     double cos = transfer * normalizedAvgU;
        //     assert(cos >= -(1 + 0.1e-9) && cos <= 1 + 0.1e-9);
        //     if (cos >= maxCosEdge) {
        //         maxCosEdge = cos;
        //         edgeIDWithMaxCos = edgeID;
        //     }
        // }
        // assert(edgeIDWithMaxCos != -1);

        // Edge *edge = &mesh->edges[edgeIDWithMaxCos];
        // Data *centerData = &mesh->nodes[edge->getNearInnerNode(node->ID)].data;
        // Data *oppositeData = &mesh->nodes[edge->getFarInnerNode(node->ID)].data;

        // double phiOpposite2 = oppositeData->phi2;
        // double phiCenter2 = centerData->phi2;

        // double phi2 = 2 * phiCenter2 - phiOpposite2;

        // double phiOpposite0 = oppositeData->phi0;
        // double phiCenter0 = centerData->phi0;
        // double phi0 = data->phi0;
        // Vector transfer = (data->coords - oppositeData->coords);
        // transfer = transfer / transfer.length();
        // double h = (oppositeData->coords - data->coords).length();
        // double Q = (phiCenter2 - phiCenter0) / tau +
        //             (centerData->u * transfer) *
        //             (phi0 - phiOpposite0) / h;
        // double min = std::min(std::min(phi0, phiOpposite0), phiCenter0) + tau * Q;
        // double max = std::max(std::max(phi0, phiOpposite0), phiCenter0) + tau * Q;

        // phi2 = std::max(phi2, min);
        // phi2 = std::min(phi2, max);

        // data->phi2 = phi2;
        // node->phase2Calculated = true;
        // --- ---

        // --- Calculating new value with transfer vectors from cells ---
        double maxCos = -1;
        long cellIDWithMaxCos = -1;
        for (unsigned long cellID : node->cellIDs) {
            Vector transfer = mesh->cells[cellID].nodeToTransferVector[node->ID];
            double cos = transfer * normalizedAvgU;
            assert(cos >= -1 && cos <= 1);
            if (cos >= maxCos) {
                maxCos = cos;
                cellIDWithMaxCos = cellID;
            }
        }
        assert(cellIDWithMaxCos != -1);

        Cell *cell = &mesh->cells[cellIDWithMaxCos];
        Data *centerData = &mesh->nodes[cell->centerNodeID].data;

        arma::vec phi0(10);
        int pos = 0;
        for (long nodeID : cell->nodeIDs) {
            phi0[pos] = mesh->nodes[nodeID].data.phi0;
            pos++;
        }
        for (long edgeID : cell->edgeIDs) {
            for (long nodeID : mesh->edges[edgeID].getInnerNodes()) {
                phi0[pos] = mesh->nodes[nodeID].data.phi0;
                pos++;
            }
        }
        phi0[pos] = centerData->phi0;
        arma::vec a = arma::solve(cell->interpolationMat, phi0);
        
        double coef = -(avgU * cell->nodeToTransferVector[node->ID]) * tau;
        Vector intersectCoords = cell->nodeToTransferVector[node->ID] * coef;
        intersectCoords = intersectCoords + data->coords;
        double intersectX = intersectCoords.x;
        double intersectY = intersectCoords.y;

        double phi2 = 0;
        phi2 += a[0];
        phi2 += a[1] * intersectX;
        phi2 += a[2] * intersectY;
        phi2 += a[3] * intersectX * intersectX;
        phi2 += a[4] * intersectY * intersectY;
        phi2 += a[5] * intersectX * intersectY;
        phi2 += a[6] * intersectX * intersectX * intersectX;
        phi2 += a[7] * intersectY * intersectY * intersectY;
        phi2 += a[8] * intersectX * intersectY * intersectX;
        phi2 += a[9] * intersectX * intersectY * intersectY;

        double phiLeft0 = centerData->phi0;
        double phiCenter0 = centerData->phi0;
        double phiRight0 = data->phi0;

        double phiCenter1 = centerData->phi1;

        double h = (centerData->coords - data->coords).length();
        double Q = (phiCenter1 - phiCenter0) / tau / 2 +
                    (((data->u + centerData->u) / 2) * cell->nodeToTransferVector[node->ID]) *
                    (phiRight0 - phiLeft0) / h;
        double min = std::min(std::min(phiRight0, phiLeft0), phiCenter0) + tau * Q;
        double max = std::max(std::max(phiRight0, phiLeft0), phiCenter0) + tau * Q;

        phi2 = std::max(phi2, min);
        phi2 = std::min(phi2, max);

        data->phi2 = phi2;
        node->phase2Calculated = true;
        // --- ---

        // --- Calculating new value with divergence -- 
        // double div = 0;
        // double volume = 0;
        // for (unsigned long cellID : node->cellIDs) {
        //     Cell *cell = &mesh->cells[cellID];
        //     volume += cell->volume;
        //     for (long edgeID : cell->edgeIDs) {
        //         div += calculateDivOnEdge(&mesh->edges[edgeID], 1) * cell->edgeToNormalDir[edgeID];
        //     }
        // }
        // div /= volume;
        // data->phi2 = data->phi0 - tau * div;
        // node->phase2Calculated = true;
        // --- ---
    }
}

void Solver::processPhase3(double tau) {
    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];
        Node *cellCenter = &mesh->nodes[cell->centerNodeID];

        double div = 0;
        for (long edgeID : cell->edgeIDs) {
            div += calculateDivOnEdge(&mesh->edges[edgeID], 3) * cell->edgeToNormalDir[edgeID];
        }
        div /= cell->volume;

        cellCenter->data.phi2 = cellCenter->data.phi1 - tau * div / 2;
    }
}

void Solver::prepareNextStep() {
    for (unsigned long i = 0; i < mesh->nodes.size(); i++) {
        Node *node = &mesh->nodes[i];
        node->phase2Calculated = false;

        Data *data = &node->data;
        data->phi0 = data->phi2;
    }
}

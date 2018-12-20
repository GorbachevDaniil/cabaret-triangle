#include "Solver.hpp"

#include "Cell.hpp"
#include "Edge.hpp"
#include "Node.hpp"
#include "Parameters.hpp"
#include "Vector.hpp"

#include <limits>
#include <cassert>
#include <armadillo>

double Solver::calculateTau() {
    double tau = std::numeric_limits<double>::max();
    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];
        Node *cellCenter = &mesh->nodes[cell->centerNodeID];
        for (unsigned long edgeID : cell->edgeIDs) {
            for (unsigned long nodeID : mesh->edges[edgeID].getUsedNodes(*mesh)) {
                double lambda = cellCenter->data.u * cell->nodeToTransferVector[nodeID];
                double edgeTau = cell->maxH / std::abs(lambda);
                if (tau > edgeTau) {
                    tau = edgeTau;
                }
            }
        }
    }
    assert(tau != std::numeric_limits<double>::max());
    return Parameters::CFL * tau;
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
            // TODO
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

Cell *getCellToProcessInnerNodes(Mesh *mesh, Edge *edge) {
    if (edge->boundEdge) {
        return &mesh->cells[edge->cellIDs[0]];
    }

    Cell *cell1 = &mesh->cells[edge->cellIDs[0]];
    Cell *cell2 = &mesh->cells[edge->cellIDs[1]];
    Node *cellCenter1 = &mesh->nodes[cell1->centerNodeID];
    Node *cellCenter2 = &mesh->nodes[cell2->centerNodeID];

    Vector uAverage((cellCenter1->data.u + cellCenter2->data.u) / 2);

    Vector normal(edge->normal * cell1->edgeToNormalDir[edge->ID]);
    if (normal * uAverage >= 0) {
        return cell1;
    } else {
        return cell2;
    }
}

void Solver::processPhase2(double tau) {
    for (unsigned long i = 0; i < mesh->edges.size(); i++) {
        Edge *edge = &mesh->edges[i];
        Cell *cell = getCellToProcessInnerNodes(mesh, edge);

        Data *cellCenterData = &mesh->nodes[cell->centerNodeID].data;
        for (unsigned long nodeID : edge->getInnerNodes()) {
            Node *node = &mesh->nodes[nodeID];
            Data *nodeData = &node->data;
            assert(cell->nodeIDToOppositeNodeID[nodeID] != -1);
            Data *oppositeNodeDate = &mesh->nodes[cell->nodeIDToOppositeNodeID[nodeID]].data;

            double phiLeft0 = oppositeNodeDate->phi0;
            double phiCenter0 = cellCenterData->phi0;
            double phiRight0 = nodeData->phi0;

            double phiCenter1 = cellCenterData->phi1;

            double phiRight2 = 2 * phiCenter1 - phiLeft0;

            double h = 2 * (cellCenterData->coords - nodeData->coords).length();
            double Q = (phiCenter1 - phiCenter0) / tau / 2 +
                       (cellCenterData->u * cell->nodeToTransferVector[nodeID]) *
                        (phiRight0 - phiLeft0) / h;
            double min = std::min(std::min(phiRight0, phiLeft0), phiCenter0) + tau * Q;
            double max = std::max(std::max(phiRight0, phiLeft0), phiCenter0) + tau * Q;

            phiRight2 = std::max(phiRight2, min);
            phiRight2 = std::min(phiRight2, max);

            nodeData->phi2 = phiRight2;
            node->phase2Calculated = true;
        }
    }

    for (unsigned long i = 0; i < mesh->nodes.size(); i++) {
        Node *node = &mesh->nodes[i];
        if (node->phase2Calculated || node->cellCenterNode) {
            continue;
        }

        double uX = 0;
        double uY = 0;
        for (unsigned long cellID : node->cellIDs) {
            Vector u = mesh->nodes[mesh->cells[cellID].centerNodeID].data.u;
            uX += u.x;
            uY += u.y;
        }
        Vector avgU(uX, uY);
        avgU = avgU / node->cellIDs.size();
        avgU = avgU / avgU.length();

        double maxCos = -1;
        long cellIDWithMaxCos = -1;
        for (unsigned long cellID : node->cellIDs) {
            Vector transfer = mesh->cells[cellID].nodeToTransferVector[node->ID];
            double cos = transfer * avgU;
            assert(cos >= -1 && cos <= 1);
            if (cos >= maxCos) {
                maxCos = cos;
                cellIDWithMaxCos = cellID;
            }
        }
        assert(cellIDWithMaxCos != -1);

        Cell *cell = &mesh->cells[cellIDWithMaxCos];
        Node *centerNode = &mesh->nodes[cell->centerNodeID];

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
        phi0[pos] = centerNode->data.phi0;
        arma::vec a = arma::solve(cell->interpolationMat, phi0);
        
        double coef = -(centerNode->data.u * cell->nodeToTransferVector[node->ID]) * tau;
        Vector intersectCoords = cell->nodeToTransferVector[node->ID] * coef;
        intersectCoords = intersectCoords + node->data.coords;
        double intersectX = intersectCoords.x;
        double intersectY = intersectCoords.y;

        double newPhi = 0;
        newPhi += a[0];
        newPhi += a[1] * intersectX;
        newPhi += a[2] * intersectY;
        newPhi += a[3] * intersectX * intersectX;
        newPhi += a[4] * intersectY * intersectY;
        newPhi += a[5] * intersectX * intersectY;
        newPhi += a[6] * intersectX * intersectX * intersectX;
        newPhi += a[7] * intersectY * intersectY * intersectY;
        newPhi += a[8] * intersectX * intersectY * intersectX;
        newPhi += a[9] * intersectX * intersectY * intersectY;

        Data *data = &node->data;

        double h = (centerNode->data.coords - data->coords).length();
        double phiLeft0 = centerNode->data.phi0;
        double phiCenter0 = centerNode->data.phi0;
        double phiRight0 = data->phi0;

        double phiCenter1 = centerNode->data.phi1;

        double Q = (phiCenter1 - phiCenter0) / tau / 2 +
                    (centerNode->data.u * cell->nodeToTransferVector[node->ID]) *
                    (phiRight0 - phiLeft0) / h;
        double min = std::min(std::min(phiRight0, phiLeft0), phiCenter0) + tau * Q;
        double max = std::max(std::max(phiRight0, phiLeft0), phiCenter0) + tau * Q;

        newPhi = std::max(newPhi, min);
        newPhi = std::min(newPhi, max);

        data->phi2 = newPhi;
        node->phase2Calculated = true;
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

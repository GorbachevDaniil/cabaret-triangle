#include "Solver.hpp"

#include "Cell.hpp"
#include "Edge.hpp"
#include "Node.hpp"
#include "Parameters.hpp"
#include "Vector.hpp"

#include <limits>
#include <cassert>

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
    Node *edgeNode1 = &mesh->nodes[edge->getUsedNodes(*mesh)[0]];
    Node *edgeNode2 = &mesh->nodes[edge->getUsedNodes(*mesh)[1]];
    double phi1 = phase == 1 ? edgeNode1->data.phi0 : edgeNode1->data.phi2;
    double phi2 = phase == 1 ? edgeNode2->data.phi0 : edgeNode2->data.phi2;
    double phi = (phi1 + phi2) / 2;
    return phi * (edgeNode1->data.u * edge->normal) * edge->length;
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

Cell *getCellToProcessPhase2(Mesh *mesh, Edge *edge) {
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
        Cell *cell = getCellToProcessPhase2(mesh, edge);

        Data *cellCenterData = &mesh->nodes[cell->centerNodeID].data;
        for (unsigned long nodeID : edge->getUsedNodes(*mesh)) {
            Data *nodeData = &mesh->nodes[nodeID].data;
            assert(cell->nodeIDToOppositeNodeID[nodeID] != -1);
            Data *oppositeNodeDate = &mesh->nodes[cell->nodeIDToOppositeNodeID[nodeID]].data;

            double phiLeft0 = oppositeNodeDate->phi0;
            double phiCenter0 = cellCenterData->phi0;
            double phiRight0 = nodeData->phi0;

            double phiCenter1 = cellCenterData->phi1;

            double phiRight2 = 2 * phiCenter1 - phiLeft0;

            double Q = (phiCenter1 - phiCenter0) / tau / 2 +
                       (cellCenterData->u * cell->nodeToTransferVector[nodeID]) *
                           (phiRight0 - phiLeft0) / (2. / 3.) / cell->edgeToMedianLength[edge->ID];
            double min = std::min(std::min(phiRight0, phiLeft0), phiCenter0) + tau * Q;
            double max = std::max(std::max(phiRight0, phiLeft0), phiCenter0) + tau * Q;

            phiRight2 = std::max(phiRight2, min);
            phiRight2 = std::min(phiRight2, max);

            nodeData->phi2 = phiRight2;
            
            // Data *nodeData = &mesh->nodes[nodeID].data;

            // double phiAveraged0 = 0;
            // for (unsigned long edgeID : cell->edgeIDs) {
            //     Data *data = &mesh->nodes[mesh->edges[edgeID].getUsedNodes(*mesh)[0]].data;
            //     phiAveraged0 += data->phi0;
            // }
            // double phiLeft0 = (phiAveraged0 / 3 + cellCenterData->phi0) - nodeData->phi0;
            // double phiCenter0 = cellCenterData->phi0;
            // double phiRight0 = nodeData->phi0;

            // double phiCenter1 = cellCenterData->phi1;

            // double phiRight2 = 2 * phiCenter1 - phiLeft0;

            // double Q = (phiCenter1 - phiCenter0) / tau / 2 +
            //            (cellCenterData->u * cell->nodeToTransferVector[nodeID]) *
            //                (phiRight0 - phiLeft0) / (2. / 3.) / cell->edgeToMedianLength[edge->ID];
            // double min = std::min(std::min(phiRight0, phiLeft0), phiCenter0) + tau * Q;
            // double max = std::max(std::max(phiRight0, phiLeft0), phiCenter0) + tau * Q;

            // phiRight2 = std::max(phiRight2, min);
            // phiRight2 = std::min(phiRight2, max);

            // nodeData->phi2 = phiRight2;
        }
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

// void Solver::harmonise() {
//     for (unsigned long i = 0; i < mesh->edges.size(); i++) {
//         Edge *edge = &mesh->edges[i];
//         Data *data = &mesh->nodes[edge->getUsedNodes(*mesh)].data;
//         Data *dataCell1 = &mesh->nodes[mesh->cells[edge->cellIDs[0]].centerNodeID].data;
//         Data *dataCell2 = &mesh->nodes[mesh->cells[edge->cellIDs[1]].centerNodeID].data;
//         data->phi0 = (dataCell1->phi0 + dataCell2->phi0) / 2;
//     }
// }

void Solver::prepareNextStep() {
    for (unsigned long i = 0; i < mesh->nodes.size(); i++) {
        Data *data = &mesh->nodes[i].data;
        data->phi0 = data->phi2;
    }
}

#include "Solver.hpp"

#include "Cell.hpp"
#include "Edge.hpp"
#include "Node.hpp"
#include "Vector.hpp"
#include "Parameters.hpp"

#include <limits>

double Solver::calculateTau() {
    double tau = std::numeric_limits<double>::max();
    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];
        Node *cellCenter = &mesh->nodes[cell->centerNodeID];
        for (unsigned long edgeID : cell->edgeIDs) {
            double lambda = cellCenter->data.u * cell->edgeToTransportDir[edgeID];
            if (tau > cell->maxH / std::abs(lambda)) {
                tau = cell->maxH / std::abs(lambda);
            }
        }
    }
    return Parameters::CFL * tau;
}

void Solver::processPhase1(double tau) {
    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];
        Node *cellCenter = &mesh->nodes[cell->centerNodeID];
        double div = 0;

        for (long edgeID : cell->edgeIDs) {
            Edge *edge = &mesh->edges[edgeID];
            Node *edgeCenter = &mesh->nodes[edge->centerNodeID];

            double uProj = edgeCenter->data.u * edgeCenter->normal;
            uProj *= cell->edgeToNormalDir[edgeID];
            div += edgeCenter->data.phi0 * uProj * edge->length;
        }
        div /= cell->volume;
        cellCenter->data.phi1 = cellCenter->data.phi0 - tau * div / 2;
    }
}

Cell* getCellToProcessPhase2(Mesh *mesh, Edge *edge) {
    if (edge->boundEdge) {
        return &mesh->cells[edge->cellIDs[0]];
    }

    Cell *cell1 = &mesh->cells[edge->cellIDs[0]];
    Cell *cell2 = &mesh->cells[edge->cellIDs[1]];
    Node *cellCenter1 = &mesh->nodes[cell1->centerNodeID];
    Node *cellCenter2 = &mesh->nodes[cell2->centerNodeID];

    Vector uAverage((cellCenter1->data.u + cellCenter2->data.u) / 2);

    Node *edgeCenter = &mesh->nodes[edge->centerNodeID];
    Vector normal(edgeCenter->normal *  cell1->edgeToNormalDir[edge->ID]);
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

        Data *edgeCenterData = &mesh->nodes[edge->centerNodeID].data;
        Data *cellCenterData = &mesh->nodes[cell->centerNodeID].data;

        double phiAveraged0 = 0;
        for (unsigned long edgeID : cell->edgeIDs) {
            Data *data = &mesh->nodes[mesh->edges[edgeID].centerNodeID].data;
            phiAveraged0 += data->phi0;
        }
        double phiLeft0 = (phiAveraged0 / 3 + cellCenterData->phi0) - edgeCenterData->phi0;
        double phiCenter0 = cellCenterData->phi0;
        double phiRight0 = edgeCenterData->phi0;

        double phiCenter1 = cellCenterData->phi1;

        double phiRight2 = 2 * phiCenter1 - phiLeft0;

        double Q = (phiCenter1 - phiCenter0) / tau / 2 +
                (cellCenterData->u * cell->edgeToTransportDir[edge->ID]) *
                (phiRight0 - phiLeft0) / (2. / 3.) / cell->edgeToMedianLength[edge->ID];
        double min = std::min(std::min(phiRight0, phiLeft0), phiCenter0) + tau * Q;
        double max = std::max(std::max(phiRight0, phiLeft0), phiCenter0) + tau * Q;

        phiRight2 = std::max(phiRight2, min);
        phiRight2 = std::min(phiRight2, max);

        edgeCenterData->phi2 = phiRight2;
    }
}

void Solver::processPhase3(double tau) {
    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];
        Node *cellCenter = &mesh->nodes[cell->centerNodeID];
        double div = 0;

        for (long edgeID : cell->edgeIDs) {
            Edge *edge = &mesh->edges[edgeID];
            Node *edgeCenter = &mesh->nodes[edge->centerNodeID];

            double uProj = edgeCenter->data.u * edgeCenter->normal;
            uProj *= cell->edgeToNormalDir[edgeID];
            div += edgeCenter->data.phi2 * uProj * edge->length;
        }
        div /= cell->volume;
        cellCenter->data.phi2 = cellCenter->data.phi1 - tau * div / 2;
    }
}

void Solver::harmonise() {
    for (unsigned long i = 0; i < mesh->edges.size(); i++) {
        Edge *edge = &mesh->edges[i];
        Data *data = &mesh->nodes[edge->centerNodeID].data;
        Data *dataCell1 = &mesh->nodes[mesh->cells[edge->cellIDs[0]].centerNodeID].data;
        Data *dataCell2 = &mesh->nodes[mesh->cells[edge->cellIDs[1]].centerNodeID].data;
        data->phi0 = (dataCell1->phi0 + dataCell2->phi0) / 2;
    }
}


void Solver::prepareNextStep() {
    for (unsigned long i = 0; i < mesh->nodes.size(); i++) {
        Data *data = &mesh->nodes[i].data;
        data->phi0 = data->phi2;
    }
}

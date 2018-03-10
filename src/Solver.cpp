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

        Node *edgeCenter = &mesh->nodes[edge->centerNodeID];
        Node *cellCenter = &mesh->nodes[cell->centerNodeID];

        double phiRight0 = edgeCenter->data.phi0;
        double phiCenter0 = cellCenter->data.phi0;
        double phiLeft0 = 2 * phiCenter0 - phiRight0;

        double phiCenter1 = cellCenter->data.phi1;

        double phiRight2 = 2 * phiCenter1 - phiLeft0;

        // TODO add monotonization

        edgeCenter->data.phi2 = phiRight2;
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

void Solver::prepareNextStep() {
    for (unsigned long i = 0; i < mesh->nodes.size(); i++) {
        Data *data = &mesh->nodes[i].data;
        data->phi0 = data->phi2;
    }
}

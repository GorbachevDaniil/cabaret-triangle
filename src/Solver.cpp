#include "Solver.hpp"
#include "Cell.hpp"
#include "Edge.hpp"
#include "Node.hpp"
#include "Vector.hpp"
#include "Parameters.hpp"

void Solver::processPhase1() {
    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];
        Node *cellCenter = &mesh->nodes[cell->centerNodeID];
        double div = 0;

        for (long edgeID : cell->edgeIDs) {
            Edge edge = mesh->edges[edgeID];
            Node edgeCenter = mesh->nodes[edge.centerNodeID];

            div += edgeCenter.data.phi0 * edge.length * cell->edgeToNormalDirection[edgeID]
                    * Vector::Scalar(edgeCenter.data.u, edgeCenter.normal);
        }
        cellCenter->data.phi1 = cellCenter->data.phi0 + Parameters::CFL * div / cell->volume;
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

    Vector u1Projection;
    u1Projection.set(cellCenter1->data.u);
    u1Projection.mult(Vector::Scalar(cellCenter1->data.u, cell1->edgeToVectorFromCenter[edge->ID]));
    Vector u2Projection;
    u2Projection.set(cellCenter2->data.u);
    u2Projection.mult(Vector::Scalar(cellCenter2->data.u, cell2->edgeToVectorFromCenter[edge->ID]));

    Vector uAverageProjection;
    uAverageProjection.plus(u1Projection);
    uAverageProjection.plus(u2Projection);
    uAverageProjection.mult(0.5);

    if (Vector::Scalar(u1Projection, uAverageProjection) >= 0) {
        return cell1;
    } else {
        return cell2;
    }
}

void Solver::processPhase2() {
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

void Solver::processPhase3() {
    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];
        Node *cellCenter = &mesh->nodes[cell->centerNodeID];
        double div = 0;

        for (long edgeID : cell->edgeIDs) {
            Edge edge = mesh->edges[edgeID];
            Node edgeCenter = mesh->nodes[edge.centerNodeID];

            div += edgeCenter.data.phi2 * edge.length * cell->edgeToNormalDirection[edgeID]
                    * Vector::Scalar(edgeCenter.data.u, edgeCenter.normal);
        }
        cellCenter->data.phi2 = cellCenter->data.phi1 + Parameters::CFL * div / cell->volume;
    }
}

void Solver::prepareNextStep() {
    for (unsigned long i = 0; i < mesh->nodes.size(); i++) {
        Data *data = &mesh->nodes[i].data;
        data->phi0 = data->phi2;
    }
}

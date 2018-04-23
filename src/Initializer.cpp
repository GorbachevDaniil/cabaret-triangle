#include "Initializer.hpp"

void Initializer::initialize(Mesh &mesh) {
    for (unsigned long i = 0; i < mesh.cells.size(); i++) {
        Data *data = &mesh.nodes[mesh.cells[i].centerNodeID].data;
        if (data->coords.x < -0.5 && data->coords.x > -0.9 &&
                data->coords.y < -0.5 && data->coords.y > -0.9) {
            data->phi0 = 1;
        } else {
            data->phi0 = 0;
        }
        data->u = Vector(1, 1);
    }
    for (unsigned long i = 0; i < mesh.edges.size(); i++) {
        Edge *edge = &mesh.edges[i];
        Data *data = &mesh.nodes[edge->centerNodeID].data;
        if (edge->boundEdge) {
            Data *dataCell = &mesh.nodes[mesh.cells[edge->cellIDs[0]].centerNodeID].data;
            data->phi0 = dataCell->phi0;
            data->u = dataCell->u;
        } else {
            Data *dataCell1 = &mesh.nodes[mesh.cells[edge->cellIDs[0]].centerNodeID].data;
            Data *dataCell2 = &mesh.nodes[mesh.cells[edge->cellIDs[1]].centerNodeID].data;
            data->phi0 = (dataCell1->phi0 + dataCell2->phi0) / 2;
            data->u = (dataCell1->u + dataCell2->u) / 2;
        }
    }
}

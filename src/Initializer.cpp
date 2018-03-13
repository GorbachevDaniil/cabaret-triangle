#include "Initializer.hpp"

void Initializer::initialize(Mesh &mesh) {
    for (unsigned long i = 0; i < mesh.cells.size(); i++) {
        Data *data = &mesh.nodes[mesh.cells[i].centerNodeID].data;
        if (data->coords.x < -0.5 && data->coords.x > -0.9 &&
                data->coords.y < 0.2 && data->coords.y > -0.2) {
            data->phi0 = 1;
        } else {
            data->phi0 = 0;
        }
        data->u = Vector(1, 0);
    }
    for (unsigned long i = 0; i < mesh.edges.size(); i++) {
        Edge *edge = &mesh.edges[i];
        Data *data = &mesh.nodes[edge->centerNodeID].data;
        Data *dataFromCell1 = &mesh.nodes[mesh.cells[edge->cellIDs[0]].centerNodeID].data;
        Data *dataFromCell2 = &mesh.nodes[mesh.cells[edge->cellIDs[0]].centerNodeID].data;
        data->phi0 = (dataFromCell1->phi0 + dataFromCell2->phi0) / 2;
        data->u = Vector(1, 0);
    }
}

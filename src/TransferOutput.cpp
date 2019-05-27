#include "TransferOutput.hpp"

#include <iostream>
#include <string>

void TransferOutput::writeParaview(Mesh *mesh, int step) {
    if (step % writePeriod != 0) {
        return;
    }

    FILE *output_f;

    std::string filename = "./bin/output/result." + std::to_string(step) + ".txt";
    output_f = std::fopen(filename.c_str(), "w");

    std::fprintf(output_f, "x,y,z,phi\n");

    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];
        Data *data = &mesh->nodes[cell->centerNodeID].data;
        std::fprintf(output_f, "%f,%f,%f,%f\n", data->coords.x, data->coords.y, 0.0,
                     data->s0[0]);
        for (unsigned long edgeID : cell->edgeIDs) {
            Edge *edge = &mesh->edges[edgeID];
            for (unsigned long usedNodeID : edge->getUsedNodes(*mesh)) {
                data = &mesh->nodes[usedNodeID].data;
                std::fprintf(output_f, "%f,%f,%f,%f\n", data->coords.x, data->coords.y, 0.0,
                             data->s0[0]);
            }
        }
    }
    std::fclose(output_f);
}

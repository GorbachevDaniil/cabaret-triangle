#include "equations/shallow_water/output.hpp"

#include <iostream>
#include <string>

void ShallowWaterOutput::write_paraview(Mesh *mesh, double time, int step) {
    if (step % write_period_ != 0) {
        return;
    }

    FILE *output_f;

    std::string filename = "./bin/output/result." + std::to_string(step) + ".txt";
    output_f = std::fopen(filename.c_str(), "w");

    std::fprintf(output_f, "x,y,z,t,h,ux,uy,bound\n");

    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];
        Data *data = &mesh->nodes[cell->centerNodeID].data;
        if (write_conservative_) {
            long centerNodeID = cell->centerNodeID;
            std::fprintf(output_f, "%f,%f,%f,%f,%f,%f,%f,%d\n",
                         data->coords.x, data->coords.y,
                         0.0, time, mesh->s0[centerNodeID][0],
                         mesh->v0[centerNodeID][0].x, mesh->v0[centerNodeID][0].y,
                         mesh->nodes[centerNodeID].boundNode);
        }
        if (write_flux_) {
            for (unsigned long edgeID : cell->edgeIDs) {
                Edge *edge = &mesh->edges[edgeID];
                for (unsigned long usedNodeID : edge->usedNodeIDs) {
                    Data *data = &mesh->nodes[usedNodeID].data;
                    std::fprintf(output_f, "%f,%f,%f,%f,%f,%f,%f,%d\n",
                                 data->coords.x, data->coords.y,
                                 0.0, time, mesh->s0[usedNodeID][0],
                                 mesh->v0[usedNodeID][0].x, mesh->v0[usedNodeID][0].y,
                                 mesh->nodes[usedNodeID].boundNode);
                }
            }
        }
    }
    std::fclose(output_f);
}

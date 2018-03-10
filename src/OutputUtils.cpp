#include "OutputUtils.hpp"
#include "Mesh.hpp"

#include <iostream>
#include <string>

void OutputUtils::OutputParaview(Mesh *mesh, int step) {

    FILE *output_f;

    std::string filename = "./bin/output/result." + std::to_string(step) + ".txt";
    output_f = std::fopen(filename.c_str(), "w");

    std::fprintf(output_f, "x,y,phi\n");

    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Data *data = &mesh->nodes[mesh->cells[i].centerNodeID].data;
        std::fprintf(output_f, "%f,%f,%f\n", data->coords.x, data->coords.y, data->phi0);
    }
    std::fclose(output_f);

}

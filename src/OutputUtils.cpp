#include "OutputUtils.hpp"
#include "Mesh.hpp"
#include <iostream>
#include <string>

void OutputUtils::OutputParaview(Mesh *mesh, int step) {

    FILE *output_f;

    std::string filename = "./bin/output/result." + std::to_string(step) + ".txt";
    output_f = std::fopen(filename.c_str(), "w");

    std::fprintf(output_f, "x,y,phi\n");

    for (unsigned long i = 0; i < mesh->nodes.size(); i++) {
        std::fprintf(output_f, "%f,%f,%f\n",
                        mesh->nodes[i].data.coords.x,
                        mesh->nodes[i].data.coords.y,
                        mesh->nodes[i].data.phi0);
    }
    std::fclose(output_f);

}

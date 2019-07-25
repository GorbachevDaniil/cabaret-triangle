#include "ShallowWaterOutput.hpp"

#include <iostream>
#include <string>

void ShallowWaterOutput::writeParaview(Mesh *mesh, int step) {
    if (step % writePeriod != 0) {
        return;
    }

    FILE *output_f;

    std::string filename = "./bin/output/result." + std::to_string(step) + ".txt";
    output_f = std::fopen(filename.c_str(), "w");

    std::fprintf(output_f, "x,y,z,h,ux,uy\n");

    for (unsigned long i = 0; i < mesh->cells.size(); i++) {
        Cell *cell = &mesh->cells[i];
        Data *data = &mesh->nodes[cell->centerNodeID].data;
        long centerNodeID = cell->centerNodeID;
        std::fprintf(output_f, "%f,%f,%f,%f,%f,%f\n", data->coords.x, data->coords.y, 0.0,
                     mesh->s0[centerNodeID][0], mesh->v0[centerNodeID][0].x, mesh->v0[centerNodeID][0].y);
        // for (unsigned long edgeID : cell->edgeIDs) {
        //     Edge *edge = &mesh->edges[edgeID];
        //     for (unsigned long usedNodeID : edge->getUsedNodes(*mesh)) {
        //         data = &mesh->nodes[usedNodeID].data;
        //         std::fprintf(output_f, "%f,%f,%f,%f,%f,%f\n", data->coords.x, data->coords.y, 0.0,
        //                      data->s0[0], data->v0[0].x, data->v0[0].y);
        //     }
        // }


        // if (step % 3 == 1) {
        //     std::fprintf(output_f, "%f,%f,%f,%f,%f,%f\n", data->coords.x, data->coords.y, 0.0,
        //                  data->s1[0], data->v1[0].x, data->v1[0].y);
        //     for (unsigned long edgeID : cell->edgeIDs) {
        //         // if (i == 1) {
        //         //     std::cout << edgeID << " " << cell->edgeToNormalDir[edgeID] << std::endl;
        //         // }
        //         Edge *edge = &mesh->edges[edgeID];
        //         for (unsigned long usedNodeID : edge->getUsedNodes(*mesh)) {
        //             // if (i == 1) {
        //             //     Vector n = cell->nodeToTransferVector[usedNodeID];
        //             //     std::cout << usedNodeID << " " << n.x << " " << n.y << std::endl;
        //             // }
        //             data = &mesh->nodes[usedNodeID].data;
        //             std::fprintf(output_f, "%f,%f,%f,%f,%f,%f\n", data->coords.x, data->coords.y, 0.0,
        //                          data->s0[0], data->v0[0].x, data->v0[0].y);
        //         }
        //     }
        // }
        // if (step % 3 == 2) {
        //     std::fprintf(output_f, "%f,%f,%f,%f,%f,%f\n", data->coords.x, data->coords.y, 0.0,
        //                  data->s1[0], data->v1[0].x, data->v1[0].y);
        //     for (unsigned long edgeID : cell->edgeIDs) {
        //         // if (i == 1) {
        //         //     std::cout << edgeID << " " << cell->edgeToNormalDir[edgeID] << std::endl;
        //         // }
        //         Edge *edge = &mesh->edges[edgeID];
        //         for (unsigned long usedNodeID : edge->getUsedNodes(*mesh)) {
        //             // if (i == 1) {
        //             //     Vector n = cell->nodeToTransferVector[usedNodeID];
        //             //     std::cout << usedNodeID << " " << n.x << " " << n.y << std::endl;
        //             // }
        //             data = &mesh->nodes[usedNodeID].data;
        //             std::fprintf(output_f, "%f,%f,%f,%f,%f,%f\n", data->coords.x, data->coords.y, 0.0,
        //                          data->s2[0], data->v2[0].x, data->v2[0].y);
        //         }
        //     }
        // }
        // if (step % 3 == 0) {
        //     std::fprintf(output_f, "%f,%f,%f,%f,%f,%f\n", data->coords.x, data->coords.y, 0.0,
        //                  data->s2[0], data->v2[0].x, data->v2[0].y);
        //     for (unsigned long edgeID : cell->edgeIDs) {
        //         // if (i == 1) {
        //         //     std::cout << edgeID << " " << cell->edgeToNormalDir[edgeID] << std::endl;
        //         // }
        //         Edge *edge = &mesh->edges[edgeID];
        //         for (unsigned long usedNodeID : edge->getUsedNodes(*mesh)) {
        //             // if (i == 1) {
        //             //     Vector n = cell->nodeToTransferVector[usedNodeID];
        //             //     std::cout << usedNodeID << " " << n.x << " " << n.y << std::endl;
        //             // }
        //             data = &mesh->nodes[usedNodeID].data;
        //             std::fprintf(output_f, "%f,%f,%f,%f,%f,%f\n", data->coords.x, data->coords.y, 0.0,
        //                          data->s2[0], data->v2[0].x, data->v2[0].y);
        //         }
        //     }
        // }
    }
    std::fclose(output_f);
}

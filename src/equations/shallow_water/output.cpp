#include "equations/shallow_water/output.hpp"

#include <string>

void ShallowWaterOutput::write_paraview(Mesh& mesh, double time, int step) {
    FILE *output_f;

    std::string filename = "./bin/output/result." + std::to_string(step) + ".txt";
    output_f = std::fopen(filename.c_str(), "w");

    std::fprintf(output_f, "x,y,z,t,h,ux,uy,bound\n");

    for (unsigned long i = 0; i < mesh.cells.size(); i++) {
        Cell *cell = &mesh.cells[i];
        Node *node = &mesh.nodes[cell->center_node_id];
        if (write_conservative_) {
            long centerNodeID = cell->center_node_id;
            std::fprintf(output_f, "%f,%f,%f,%f,%f,%f,%f,%d\n",
                         node->coords.x, node->coords.y,
                         0.0, time, mesh.s0[centerNodeID][0],
                         mesh.v0[centerNodeID][0].x, mesh.v0[centerNodeID][0].y,
                         mesh.nodes[centerNodeID].is_bound);
        }
        if (write_flux_) {
            for (unsigned long edgeID : cell->edge_ids) {
                Edge *edge = &mesh.edges[edgeID];
                for (unsigned long usedNodeID : edge->used_node_ids) {
                    Node *edge_node = &mesh.nodes[usedNodeID];
                    std::fprintf(output_f, "%f,%f,%f,%f,%f,%f,%f,%d\n",
                                 edge_node->coords.x, edge_node->coords.y,
                                 0.0, time, mesh.s0[usedNodeID][0],
                                 mesh.v0[usedNodeID][0].x, mesh.v0[usedNodeID][0].y,
                                 mesh.nodes[usedNodeID].is_bound);
                }
            }
        }
    }
    std::fclose(output_f);
}

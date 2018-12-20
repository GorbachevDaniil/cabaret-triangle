#include "Initializer.hpp"

void Initializer::initialize(Mesh &mesh) {
    for (unsigned long i = 0; i < mesh.edges.size(); i++) {
        Edge *edge = &mesh.edges[i];
        for (unsigned long usedNodeID : edge->getUsedNodes(mesh)) {
            Data *data = &mesh.nodes[usedNodeID].data;
            if (data->coords.x < -0.2 && data->coords.x > -0.5 &&
                    data->coords.y < -0.2 && data->coords.y > -0.5) {
                data->phi0 = 1;
            } else {
                data->phi0 = 0;
            }
            data->u = Vector(1, 1);
        }
    }
    for (unsigned long i = 0; i < mesh.cells.size(); i++) {
        Data *data = &mesh.nodes[mesh.cells[i].centerNodeID].data;
        double phi = 0;
        double cellNodesNum = 0;
        for (unsigned long nodeID : mesh.cells[i].nodeIDs) {
            phi += mesh.nodes[nodeID].data.phi0;
            cellNodesNum += 1;
        }
        for (unsigned long edgeID : mesh.cells[i].edgeIDs) {
            for (unsigned long innerNodeID : mesh.edges[edgeID].getInnerNodes()) {
                phi += mesh.nodes[innerNodeID].data.phi0;
                cellNodesNum += 1;
            }
        }
        data->phi0 = phi / cellNodesNum;
        data->u = Vector(1, 1);
    }
}

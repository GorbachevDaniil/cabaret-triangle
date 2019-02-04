#include "Initializer.hpp"

#include <cmath>

void Initializer::initialize(Mesh &mesh) {
    for (unsigned long i = 0; i < mesh.edges.size(); i++) {
        Edge *edge = &mesh.edges[i];
        for (unsigned long usedNodeID : edge->getUsedNodes(mesh)) {
            Data *data = &mesh.nodes[usedNodeID].data;
            double x = data->coords.x;
            double y = data->coords.y;

            // data->u = Vector(1, 1);
            // if ((x < -0.4) && (x > -0.6) && (y < -0.4) && (y > -0.6)) {
            //     data->phi0 = 1;
            // } else {
            //     data->phi0 = 0;
            // }

            // data->u = Vector(0.5, sqrt(3) / 2);
            // if ((x < -0.4) && (x > -0.6) && (y < -0.4) && (y > -0.6)) {
            //     data->phi0 = 1;
            // } else {
            //     data->phi0 = 0;
            // }

            // data->u = Vector(0, 1);
            // if ((x < 0.1) && (x > -0.1) && (y < -0.5) && (y > -0.7)) {
            //     data->phi0 = 1;
            // } else {
            //     data->phi0 = 0;
            // }

            // data->u = Vector(1, 0);
            // double gamma = 0.1;
            // data->phi0 = exp(-(pow(x + 0.45, 2) + pow(y, 2)) / (2 * pow(gamma, 2))) / 1;

            data->u = Vector(-y, x);
            double rad = 0.345;
            double sqrNodeRad = (pow(x - 0.35, 2) + pow(y, 2));
            double sqrRad = pow(rad, 2);
            if (sqrNodeRad <= sqrRad) {
                Vector nodeVec(x - 0.35, y);
                data->phi0 = 1 - nodeVec.length() / rad;
            } else {
                data->phi0 = 0;
            }

            // data->u = Vector(-y, x);
            // double gamma = 0.15;
            // data->phi0 = exp(-(pow(x - 0.35, 2) + pow(y, 2)) / (2 * pow(gamma, 2))) / 1;
        }
    }
    for (unsigned long i = 0; i < mesh.cells.size(); i++) {
        Data *data = &mesh.nodes[mesh.cells[i].centerNodeID].data;
        double phi = 0;
        double cellNodesNum = 0;
        Vector u(0,0);
        for (unsigned long nodeID : mesh.cells[i].nodeIDs) {
            phi += mesh.nodes[nodeID].data.phi0;
            u = u + mesh.nodes[nodeID].data.u;
            cellNodesNum += 1;
        }
        for (unsigned long edgeID : mesh.cells[i].edgeIDs) {
            for (unsigned long innerNodeID : mesh.edges[edgeID].getInnerNodes()) {
                phi += mesh.nodes[innerNodeID].data.phi0;
                u = u + mesh.nodes[innerNodeID].data.u;
                cellNodesNum += 1;
            }
        }
        data->phi0 = phi / cellNodesNum;
        data->u = u / cellNodesNum;
    }
}

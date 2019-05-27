#include "TransferInitializer.hpp"

#include <cmath>

void TransferInitializer::initialize(Mesh &mesh) {
    for (unsigned long i = 0; i < mesh.cells.size(); i++) {
        Data *data = &mesh.nodes[mesh.cells[i].centerNodeID].data;
        double x = data->coords.x;
        double y = data->coords.y;

        data->vector[0] = Vector(1, 1);
        if ((x < -0.4) && (x > -0.6) && (y < -0.4) && (y > -0.6)) {
            data->s0[0] = 1;
        } else {
            data->s0[0] = 0;
        }

        // data->vector[0] = Vector(0.5, sqrt(3) / 2);
        // if ((x < -0.4) && (x > -0.6) && (y < -0.4) && (y > -0.6)) {
        //     data->s0[0] = 1;
        // } else {
        //     data->s0[0] = 0;
        // }

        // data->vector[0] = Vector(0, 1);
        // if ((x < 0.1) && (x > -0.1) && (y < -0.5) && (y > -0.7)) {
        //     data->s0[0] = 1;
        // } else {
        //     data->s0[0] = 0;
        // }

        // data->vector[0] = Vector(1, 0);
        // double gamma = 0.1;
        // data->s0[0] = exp(-(pow(x + 0.45, 2) + pow(y, 2)) / (2 * pow(gamma, 2))) / 1;

        // data->vector[0] = Vector(-y, x);
        // double rad = 0.345;
        // double sqrNodeRad = (pow(x - 0.346, 2) + pow(y, 2));
        // double sqrRad = pow(rad, 2);
        // if (sqrNodeRad <= sqrRad) {
        //     Vector nodeVec(x - 0.346, y);
        //     data->s0[0] = 1 - nodeVec.length() / rad;
        // } else {
        //     data->s0[0] = 0;
        // }

        // data->vector[0] = Vector(-y, x);
        // double gamma = 0.15;
        // data->s0[0] = exp(-(pow(x - 0.35, 2) + pow(y, 2)) / (2 * pow(gamma, 2))) / 1;
    }

    for (unsigned long i = 0; i < mesh.nodes.size(); i++) {
        Node *node = &mesh.nodes[i];
        if (!node->isApex) {
            continue;
        }
        double phi = 0;
        Vector u = Vector(0, 0);
        int cellsNum = 0;
        for (unsigned long cellID : node->cellIDs) {
            Data *data = &mesh.nodes[mesh.cells[cellID].centerNodeID].data;
            phi = phi + data->s0[0];
            u = u + data->vector[0];
            cellsNum++;
        }
        Data *data = &node->data;
        data->s0[0] = phi / cellsNum;
        data->vector[0] = u / cellsNum;
    }
    for (unsigned long i = 0; i < mesh.edges.size(); i++) {
        Edge *edge = &mesh.edges[i];
        Data *endData1 = &mesh.nodes[edge->nodeIDs.front()].data;
        Data *endData2 = &mesh.nodes[edge->nodeIDs.back()].data;
        int intervalNum = mesh.edgeInnerNodesNumber + 1;
        for (unsigned long j = 0; j < edge->nodeIDs.size(); j++) {
            Data *data = &mesh.nodes[edge->nodeIDs[j]].data;
            data->s0[0] = endData1->s0[0] + (endData2->s0[0] - endData1->s0[0]) * j / intervalNum;
            data->vector[0] = endData1->vector[0] + (endData2->vector[0] - endData1->vector[0]) * j / intervalNum;
        }
    }
}

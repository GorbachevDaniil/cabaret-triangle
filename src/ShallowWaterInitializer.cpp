#include "ShallowWaterInitializer.hpp"

#include <cmath>

void ShallowWaterInitializer::initialize(Mesh &mesh) {
    for (unsigned long i = 0; i < mesh.edges.size(); i++) {
        long edgeID = mesh.edges[i].ID;
        mesh.edgeIDToDivs[edgeID] = std::vector<double>(3);
    }

    for (unsigned long i = 0; i < mesh.cells.size(); i++) {
        long cellID = mesh.cells[i].ID;
        mesh.cellIDToGrads[cellID] = std::vector<Vector>(3);
    }

    for (unsigned long i = 0; i < mesh.nodes.size(); i++) {
        Data *data = &mesh.nodes[i].data;
        data->s1[0] = 0;
        data->s2[0] = 0;
        data->v1[0] = Vector(0, 0);
        data->v2[0] = Vector(0, 0);
    }

    for (unsigned long i = 0; i < mesh.cells.size(); i++) {
        Data *data = &mesh.nodes[mesh.cells[i].centerNodeID].data;
        double x = data->coords.x;
        double y = data->coords.y;

        // data->v0[0] = Vector(0, 0);
        // data->s0[0] = 1;

        // data->v0[0] = Vector(0, 0);
        // if (x < 0) {
        //     data->s0[0] = 1;
        // } else {
        //     data->s0[0] = 2;
        // }

        // data->v0[0] = Vector(0, 0);
        // double gamma = 0.07;
        // data->s0[0] = 1 - exp(-(pow(x, 2) + pow(y, 2)) / (2 * pow(gamma, 2))) / 10;

        // double alpha = 0.404;
        // double beta = 0.3;
        // double r_0 = 0.12;
        // double r = sqrt(pow(x, 2) + pow(y, 2));
        // data->s0[0] = 1 - pow(alpha, 2) * exp(2 * beta * (1 - pow(r / r_0, 2))) / (4 * beta);
        // data->v0[0] = Vector(y, -x) * (alpha * exp(beta * (1 - pow(r / r_0, 2))) / r_0);

        // double alpha = 0.2;
        // double beta = 0.6;
        // double r_0 = 0.09;
        // double x_step = 0.175;
        // double r_1 = sqrt(pow(x - x_step, 2) + pow(y, 2));
        // double r_2 = sqrt(pow(x + x_step, 2) + pow(y, 2));
        // data->s0[0] = 1 - (pow(alpha, 2) * exp(2 * beta * (1 - pow(r_1 / r_0, 2))) / (4 * beta) + 
        //                    pow(alpha, 2) * exp(2 * beta * (1 - pow(r_2 / r_0, 2))) / (4 * beta));
        // data->v0[0] = Vector(y, -(x - x_step)) * (alpha * exp(beta * (1 - pow(r_1 / r_0, 2))) / r_0) + 
        //               Vector(-y, (x + x_step)) * (alpha * exp(beta * (1 - pow(r_2 / r_0, 2))) / r_0);


        // Positivity-Preserving Well-Balanced Discontinuous Galerkin Methods
        // for the Shallow Water Equations on Unstructured Triangular Meshes.
        // Yulong Xing, Xiangxiong Zhang 
        // data->v0[0] = Vector(0, 0);
        // data->s0[0] = 1 + exp(-100 * (pow(x - 0.5, 2) + pow(y - 0.5, 2))) / 10;

        // Solution of the 2D shallow water equations using the 
        // finite volume method on unstructured triangular meshes. 
        // Anastasiou K., Chan C. T.
        data->v0[0] = Vector(0, 0);
        if ((pow(x - 25, 2) + pow(y - 25, 2)) <= 11 * 11) {
            data->s0[0] = 10;
        } else {
            data->s0[0] = 1;
        }
    }

    for (unsigned long i = 0; i < mesh.nodes.size(); i++) {
        Node *node = &mesh.nodes[i];
        if (!node->isApex) {
            continue;
        }
        double h = 0;
        Vector u = Vector(0, 0);
        int cellsNum = 0;
        for (unsigned long cellID : node->cellIDs) {
            Data *data = &mesh.nodes[mesh.cells[cellID].centerNodeID].data;
            h = h + data->s0[0];
            u = u + data->v0[0];
            cellsNum++;
        }
        Data *data = &node->data;
        data->s0[0] = h / cellsNum;
        data->v0[0] = u / cellsNum;
    }
    for (unsigned long i = 0; i < mesh.edges.size(); i++) {
        Edge *edge = &mesh.edges[i];
        Data *endData1 = &mesh.nodes[edge->nodeIDs.front()].data;
        Data *endData2 = &mesh.nodes[edge->nodeIDs.back()].data;
        int intervalNum = mesh.edgeInnerNodesNumber + 1;
        for (unsigned long j = 0; j < edge->nodeIDs.size(); j++) {
            Data *data = &mesh.nodes[edge->nodeIDs[j]].data;
            data->s0[0] = endData1->s0[0] + (endData2->s0[0] - endData1->s0[0]) * j / intervalNum;
            data->v0[0] = endData1->v0[0] + (endData2->v0[0] - endData1->v0[0]) * j / intervalNum;
        }
    }
}

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
        mesh.s0[i] = std::vector<double>(1);
        mesh.s1[i] = std::vector<double>(1);
        mesh.s2[i] = std::vector<double>(1);

        mesh.v0[i] = std::vector<Vector>(1);
        mesh.v1[i] = std::vector<Vector>(1);
        mesh.v2[i] = std::vector<Vector>(1);

        mesh.s1[i][0] = 0;
        mesh.s2[i][0] = 0;
        mesh.v1[i][0] = Vector(0, 0);
        mesh.v2[i][0] = Vector(0, 0);
    }

    // for (unsigned long i = 0; i < mesh.cells.size(); i++) {
    //     Data *data = &mesh.nodes[mesh.cells[i].centerNodeID].data;
    for (unsigned long i = 0; i < mesh.nodes.size(); i++) {
        Data *data = &mesh.nodes[i].data;
        double x = data->coords.x;
        double y = data->coords.y;

        // long centerNodeID = mesh.cells[i].centerNodeID;
        long centerNodeID = i;

        mesh.v0[centerNodeID][0] = Vector(0, 0);
        mesh.s0[centerNodeID][0] = 1;

        // mesh.v0[centerNodeID][0] = Vector(0, 0);
        // if (x < 25) {
        //     mesh.s0[centerNodeID][0] = 1;
        // } else {
        //     mesh.s0[centerNodeID][0] = 2;
        // }

        // mesh.v0[centerNodeID][0] = Vector(0, 0);
        // double gamma = 5.0;
        // mesh.s0[centerNodeID][0] = 1 - exp(-(pow(x - 25, 2) + pow(y - 25, 2)) / (2 * pow(gamma, 2))) / 1000;

        // double alpha = 0.7;
        // double beta = 0.3;
        // double r_0 = 0.15;
        // double r = sqrt(pow(x, 2) + pow(y, 2));
        // double teta = alpha * exp(beta * (1 - pow(r / r_0, 2)));
        // mesh.s0[centerNodeID][0] = 1 - pow(teta, 2) / (4 * beta);
        // mesh.v0[centerNodeID][0].x = teta * y / r_0;
        // mesh.v0[centerNodeID][0].y = -teta * x / r_0;

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
        // mesh.v0[centerNodeID][0] = Vector(0, 0);
        // mesh.s0[centerNodeID][0] = 1 + exp(-100 * (pow(x - 0.5, 2) + pow(y - 0.5, 2))) / 10;

        // Solution of the 2D shallow water equations using the 
        // finite volume method on unstructured triangular meshes. 
        // Anastasiou K., Chan C. T.
        // Circular dam break
        // mesh.v0[centerNodeID][0] = Vector(0, 0);
        // if ((pow(x - 25, 2) + pow(y - 25, 2)) <= 11 * 11) {
        //     mesh.s0[centerNodeID][0] = 10;
        // } else {
        //     mesh.s0[centerNodeID][0] = 1;
        // }

        // Solution of the 2D shallow water equations using the 
        // finite volume method on unstructured triangular meshes. 
        // Anastasiou K., Chan C. T.
        // Dam break
        // mesh.v0[centerNodeID][0] = Vector(0, 0);
        // if (x <= 95) {
        //     mesh.s0[centerNodeID][0] = 10;
        // } else {
        //     mesh.s0[centerNodeID][0] = 5;
        // }
    }

    // for (unsigned long i = 0; i < mesh.nodes.size(); i++) {
    //     Node *node = &mesh.nodes[i];
    //     if (!node->isApex) {
    //         continue;
    //     }
    //     double h = 0;
    //     Vector u = Vector(0, 0);
    //     int cellsNum = 0;
    //     for (unsigned long cellID : node->cellIDs) {
    //         long centerNodeID = mesh.cells[cellID].centerNodeID;
    //         h = h + mesh.s0[centerNodeID][0];
    //         u = u + mesh.v0[centerNodeID][0];
    //         cellsNum++;
    //     }
    //     mesh.s0[i][0] = h / cellsNum;
    //     mesh.v0[i][0] = u / cellsNum;
    // }
    // for (unsigned long i = 0; i < mesh.edges.size(); i++) {
    //     Edge *edge = &mesh.edges[i];
    //     long endNodeID1 = edge->nodeIDs.front();
    //     long endNodeID2 = edge->nodeIDs.back();
    //     int intervalNum = mesh.edgeInnerNodesNumber + 1;
    //     for (unsigned long j = 0; j < edge->nodeIDs.size(); j++) {
    //         long nodeID = edge->nodeIDs[j];
    //         mesh.s0[nodeID][0] = mesh.s0[endNodeID1][0] + (mesh.s0[endNodeID2][0] - mesh.s0[endNodeID1][0]) * j / intervalNum;
    //         mesh.v0[nodeID][0].x = mesh.v0[endNodeID1][0].x + (mesh.v0[endNodeID2][0].x - mesh.v0[endNodeID1][0].x) * j / intervalNum;
    //         mesh.v0[nodeID][0].y = mesh.v0[endNodeID1][0].y + (mesh.v0[endNodeID2][0].y - mesh.v0[endNodeID1][0].y) * j / intervalNum;
    //     }
    // }
}

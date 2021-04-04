#include "equations/shallow_water/initializer.hpp"

#include <cmath>

void ShallowWaterInitializer::initialize(Mesh &mesh) {
    for (unsigned long i = 0; i < mesh.edges.size(); i++) {
        long edgeID = mesh.edges[i].ID;
        mesh.edgeIDToDivs0[edgeID] = std::vector<double>(3);
        mesh.edgeIDToDivs2[edgeID] = std::vector<double>(3);
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

    for (unsigned long j = 0; j < mesh.edges.size(); j++) {
        Edge *edge = &mesh.edges[j];
        for (unsigned long i : edge->usedNodeIDs) {
            Data *data = &mesh.nodes[i].data;
            double x = data->coords.x;
            double y = data->coords.y;

            // mesh.v0[i][0] = Vector(0, 0);
            // mesh.s0[i][0] = 1;

            // mesh.v0[i][0] = Vector(0, 0);
            // if (x < 25) {
            //     mesh.s0[i][0] = 1;
            // } else {
            //     mesh.s0[i][0] = 2;
            // }

            // mesh.v0[i][0] = Vector(0, 0);
            // double beta = 0.05;
            // mesh.s0[i][0] = 1 + 1e-2 * exp(-beta * (pow(x - 25, 2) + pow(y - 25, 2)));

            double alpha = 0.3;
            double beta = 0.3;
            double r_0 = 0.15;
            double r = sqrt(pow(x, 2) + pow(y, 2));
            double teta = alpha * exp(beta * (1 - pow(r / r_0, 2)));
            mesh.s0[i][0] = 1 - pow(teta, 2) / (4 * beta);
            mesh.v0[i][0].x = teta * y / r_0;
            mesh.v0[i][0].y = -teta * x / r_0;

            // double alpha = 0.2;
            // double beta = 0.6;
            // double r_0 = 0.15;
            // double x_step = 0.28;
            // double r_1 = sqrt(pow(x - x_step, 2) + pow(y, 2));
            // double r_2 = sqrt(pow(x + x_step, 2) + pow(y, 2));
            // mesh.s0[i][0] = 1 - (pow(alpha, 2) * exp(2 * beta * (1 - pow(r_1 / r_0, 2))) / (4 * beta) + 
            //                      pow(alpha, 2) * exp(2 * beta * (1 - pow(r_2 / r_0, 2))) / (4 * beta));
            // mesh.v0[i][0] = Vector(-y, (x - x_step)) * (alpha * exp(beta * (1 - pow(r_1 / r_0, 2))) / r_0) + 
            //                 Vector(y, -(x + x_step)) * (alpha * exp(beta * (1 - pow(r_2 / r_0, 2))) / r_0);


            // Positivity-Preserving Well-Balanced Discontinuous Galerkin Methods
            // for the Shallow Water Equations on Unstructured Triangular Meshes.
            // Yulong Xing, Xiangxiong Zhang 
            // mesh.v0[i][0] = Vector(0, 0);
            // mesh.s0[i][0] = 1 + exp(-100 * (pow(x - 0.5, 2) + pow(y - 0.5, 2))) / 10;

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
            // mesh.v0[i][0] = Vector(0, 0);
            // if (x < 100) {
            //     mesh.s0[i][0] = 10;
            // } else {
            //     mesh.s0[i][0] = 5;
            // }
        }
    }

    for (unsigned long i = 0; i < mesh.cells.size(); i++) {
        Cell *cell = &mesh.cells[i];

        double avgH = 0;
        Vector avgU = Vector(0, 0);
        int num = 0;
        for (unsigned long edgeID : cell->edgeIDs) {
            Edge *edge = &mesh.edges[edgeID];
            for (unsigned long nodeID : edge->usedNodeIDs) {
                avgH = avgH + mesh.s0[nodeID][0];
                avgU = avgU + mesh.v0[nodeID][0];
                num++;
            }
        }

        avgH = avgH / num;
        avgU = avgU / num;
        mesh.s0[cell->centerNodeID][0] = avgH;
        mesh.v0[cell->centerNodeID][0].x = avgU.x;
        mesh.v0[cell->centerNodeID][0].y = avgU.y;
    }
}

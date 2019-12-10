#include <cassert>
#include <cmath>
#include <iostream>

double LEFT_BOUNDARY_X = 0.0;
double RIGHT_BOUNDARY_X = 50.0;
double LEFT_BOUNDARY_Y = 0.0;
double RIGHT_BOUNDARY_Y = 50.0;

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "Number of cells in x direction must be specified" << std::endl;
        return 1;
    }
    // number of cells in x direction
    int nx = std::atoi(argv[1]);
    int nIntervals = nx + 1;
    double hx = 2 * (abs(RIGHT_BOUNDARY_X) + abs(LEFT_BOUNDARY_X)) / nIntervals;
    double hy = sqrt(3) * hx / 2;
    // double hy = (abs(RIGHT_BOUNDARY_Y) + abs(LEFT_BOUNDARY_Y)) / nIntervals;
    // number of cells in y direction
    int ny = (abs(RIGHT_BOUNDARY_Y) + abs(LEFT_BOUNDARY_Y)) / hy;
    // int ny = nx + 1;
    
    int minNodesPerLayer = nx / 2 + 1;
    int maxNodesPerLayer = nx % 2 == 0 ? minNodesPerLayer : minNodesPerLayer + 1;
    unsigned long nNodes = (minNodesPerLayer + maxNodesPerLayer) + 2;
    nNodes = ny % 2 == 1 ? nNodes * ((ny + 1) / 2): nNodes * (ny / 2 + 1);
    nNodes = ny % 2 == 1 ? nNodes : nNodes - maxNodesPerLayer;
    nNodes = nx % 2 == 0 && ny % 2 == 0 ? nNodes - 1 : nNodes;

    std::cout << nx << ' ' << ny << std::endl;
    std::cout << minNodesPerLayer << ' ' << maxNodesPerLayer << std::endl;

    int minNEdgesPerLayer = minNodesPerLayer - 1;
    minNEdgesPerLayer = nx % 2 == 0 ? minNEdgesPerLayer + 1 : minNEdgesPerLayer + 2;
    int maxNEdgesPerLayer = maxNodesPerLayer - 1;
    maxNEdgesPerLayer = nx % 2 == 0 ? maxNEdgesPerLayer + 1 : maxNEdgesPerLayer;
    unsigned long nEdges = minNEdgesPerLayer + maxNEdgesPerLayer;
    nEdges = ny % 2 == 1 ? nEdges * ((ny + 1) / 2): nEdges * (ny / 2 + 1);
    nEdges += ny * (nx + 3);
    nEdges = ny % 2 == 1 ? nEdges : nEdges - maxNEdgesPerLayer;

    unsigned long nCells = (nx + 2) * ny;

    FILE *outNodes = std::fopen("Mesh.node", "w");
    FILE *outEdges = std::fopen("Mesh.edge", "w");
    FILE *outCells = std::fopen("Mesh.ele", "w");

    std::fprintf(outNodes, "%lu  %d  %d  %d\n", nNodes, 2, 0, 1);
    std::fprintf(outEdges, "%lu  %d\n", nEdges, 1);
    std::fprintf(outCells, "%lu  %d\n", nCells, 3);

    // creating node file
    unsigned long nodeID = 1;
    int nNodesPerLayer = minNodesPerLayer;
    for (int i = 0; i <= ny; i++) {
        double offset = i % 2 == 0 ? hx / 2 : 0;
        double yCoord = i * hy + LEFT_BOUNDARY_Y;
        for (int j = 0; j < nNodesPerLayer; j++) {
            double xCoord = offset + j * hx + LEFT_BOUNDARY_X;
            int bound = 0;
            if ((i == 0) || (i == ny)) {
                bound = 1;
            }
            if (j == 0) {
                if (offset > 0) {
                    std::fprintf(outNodes, "%4lu    %.17g  %.17g    %d\n", nodeID, LEFT_BOUNDARY_X, yCoord, 1);
                    nodeID++;
                }
            }
            std::fprintf(outNodes, "%4lu    %.17g  %.17g    %d\n", nodeID, xCoord, yCoord, bound);
            nodeID++;
            if (j == nNodesPerLayer - 1) {
                if (xCoord < RIGHT_BOUNDARY_X) {
                    std::fprintf(outNodes, "%4lu    %.17g  %.17g    %d\n", nodeID, RIGHT_BOUNDARY_X, yCoord, 1);
                    nodeID++;
                }
            }
        }
        nNodesPerLayer = nNodesPerLayer == minNodesPerLayer ? maxNodesPerLayer : minNodesPerLayer;
    }
    assert(nNodes == nodeID - 1);

    // creating edge file
    unsigned long edgeID = 1;
    int nEdgesPerLayer = minNEdgesPerLayer;
    unsigned long firstNodeIDForLayer = 1;
    for (int i = 0; i <= ny; i++) {
        for (int j = 0; j < nEdgesPerLayer; j++) {
            unsigned long nodeID1 = firstNodeIDForLayer + j;
            unsigned long nodeID2 = firstNodeIDForLayer + j + 1;
            int bound = 0;
            if ((i == 0) || (i == ny)) {
                bound = 1;
            }
            std::fprintf(outEdges, "%4lu   %lu  %lu  %d\n", edgeID, nodeID1, nodeID2, bound);
            edgeID++;
        }
        firstNodeIDForLayer += nEdgesPerLayer + 1;
        nEdgesPerLayer = nEdgesPerLayer == minNEdgesPerLayer ? maxNEdgesPerLayer : minNEdgesPerLayer;
    }

    bool layerHasMinNodes = true;
    nNodesPerLayer = minNEdgesPerLayer + 1;
    firstNodeIDForLayer = 1;
    for (int i = 0; i < ny; i++) {
        unsigned long botNodeID = firstNodeIDForLayer;
        unsigned long topNodeID = firstNodeIDForLayer + nNodesPerLayer;
        bool needToIncreaseBot = layerHasMinNodes ? true : false;
        for (int j = 0; j < nx + 3; j++) {
            int bound = 0;
            if ((j == 0) || (j == nx + 2)) {
                bound = 1;
            }
            std::fprintf(outEdges, "%4lu   %lu  %lu  %d\n", edgeID, botNodeID, topNodeID, bound);
            edgeID++;
            if (needToIncreaseBot) {
                botNodeID++;
                needToIncreaseBot = false;
            } else {
                topNodeID++;
                needToIncreaseBot = true;
            }
        }
        firstNodeIDForLayer += nNodesPerLayer;
        nNodesPerLayer = layerHasMinNodes ? maxNEdgesPerLayer + 1 : minNEdgesPerLayer + 1;
        layerHasMinNodes = layerHasMinNodes ? false : true;
    }
    assert(nEdges == edgeID - 1);

    // creating cell file
    unsigned long cellID = 1;
    layerHasMinNodes = true;
    nNodesPerLayer = minNEdgesPerLayer + 1;
    firstNodeIDForLayer = 1;
    for (int i = 0; i < ny; i++) {
        bool isBotTriangle = layerHasMinNodes ? true : false;
        unsigned long nodeID1 = firstNodeIDForLayer;
        unsigned long nodeID3 = firstNodeIDForLayer + nNodesPerLayer;
        unsigned long nodeID2 = isBotTriangle ? nodeID1 + 1 : nodeID3 + 1;
        for (int j = 0; j < nx + 2; j++) {
            std::fprintf(outCells, "%4lu    %4lu  %4lu  %4lu\n", cellID, nodeID1, nodeID2, nodeID3);
            if (isBotTriangle) {
                unsigned long nodeID1Buff = nodeID2;
                unsigned long nodeID2Buff = nodeID3 + 1;
                unsigned long nodeID3Buff = nodeID3;
                nodeID1 = nodeID1Buff;
                nodeID2 = nodeID2Buff;
                nodeID3 = nodeID3Buff;
                isBotTriangle = false;
            } else {
                unsigned long nodeID1Buff = nodeID1;
                unsigned long nodeID2Buff = nodeID1 + 1;
                unsigned long nodeID3Buff = nodeID2;
                nodeID1 = nodeID1Buff;
                nodeID2 = nodeID2Buff;
                nodeID3 = nodeID3Buff;
                isBotTriangle = true;
            }
            cellID++;
        }
        firstNodeIDForLayer += nNodesPerLayer;
        nNodesPerLayer = layerHasMinNodes ? maxNEdgesPerLayer + 1 : minNEdgesPerLayer + 1;
        layerHasMinNodes = layerHasMinNodes ? false : true;
    }
    assert(nCells == cellID - 1);

    std::fclose(outNodes);
    std::fclose(outEdges);
    std::fclose(outCells);

    std::cout << "number of nodes = " << nNodes << std::endl;
    std::cout << "number of edges = " << nEdges << std::endl;
    std::cout << "number of cells = " << nCells << std::endl;

    return 0;
}

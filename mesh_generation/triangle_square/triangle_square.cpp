#include <cassert>
#include <cmath>
#include <iostream>

double LEFT_BOUNDARY_X = -0.7;
double RIGHT_BOUNDARY_X = 0.7;
double LEFT_BOUNDARY_Y = -0.7;
double RIGHT_BOUNDARY_Y = 0.7;

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cout << "Number of cells by dimensions must be specified" << std::endl;
        return 1;
    }
    int nx = std::atoi(argv[1]);
    int ny = std::atoi(argv[2]);

    double hx = abs(RIGHT_BOUNDARY_X - LEFT_BOUNDARY_X) / nx;
    double hy = abs(RIGHT_BOUNDARY_Y - LEFT_BOUNDARY_Y) / ny;

    unsigned long nNodes = (nx + 1) * (ny + 1);
    unsigned long nEdges = 3 * nx * ny + nx + ny;
    unsigned long nCells = 2 * nx * ny;

    FILE *outNodes = std::fopen("Mesh.node", "w");
    FILE *outEdges = std::fopen("Mesh.edge", "w");
    FILE *outCells = std::fopen("Mesh.ele", "w");

    std::fprintf(outNodes, "%ld  %d  %d  %d\n", nNodes, 2, 0, 1);
    std::fprintf(outEdges, "%ld  %d\n", nEdges, 1);
    std::fprintf(outCells, "%ld  %d\n", nCells, 3);

    unsigned long edgeID = 1;
    unsigned long cellID = 1;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            int bound = 0;
            if (j == 0 || i == 0) {
                bound = 1;
            }
            std::fprintf(outNodes, "%4lu    %.17g  %.17g    %d\n", i * (nx + 1) + j + 1,
                         LEFT_BOUNDARY_X + j * hx, LEFT_BOUNDARY_Y + i * hy, bound);

            if (j == nx - 1) {
                bound = 1;
                std::fprintf(outNodes, "%4lu    %.17g  %.17g    %d\n", i * (nx + 1) + j + 2,
                             LEFT_BOUNDARY_X + nx * hx, LEFT_BOUNDARY_Y + i * hy, bound);
            }

            if (j == 0) {
                bound = 1;
            }
            std::fprintf(outEdges, "%4lu   %d  %d  %d\n", edgeID, i * (nx + 1) + j + 1,
                         i * (nx + 1) + j + nx + 2, bound);
            edgeID++;

            bound = 0;
            if (i == 0) {
                bound = 1;
            }
            std::fprintf(outEdges, "%4lu   %d  %d  %d\n", edgeID, i * (nx + 1) + j + 1,
                         i * (nx + 1) + j + 2, bound);
            edgeID++;

            bound = 0;
            std::fprintf(outEdges, "%4lu   %d  %d  %d\n", edgeID, i * (nx + 1) + j + 1,
                         i * (nx + 1) + j + nx + 3, bound);
            edgeID++;

            if (i == ny - 1) {
                bound = 1;
                std::fprintf(outEdges, "%4lu   %d  %d  %d\n", edgeID, i * (nx + 1) + j + nx + 2,
                             i * (nx + 1) + j + nx + 3, bound);
                edgeID++;
            }

            if (j == nx - 1) {
                bound = 1;
                std::fprintf(outEdges, "%4lu   %d  %d  %d\n", edgeID, i * (nx + 1) + j + 2,
                             i * (nx + 1) + j + nx + 3, bound);
                edgeID++;
            }

            std::fprintf(outCells, "%4lu    %4d  %4d  %4d\n", cellID, i * (nx + 1) + j + 1,
                         i * (nx + 1) + j + 2, i * (nx + 1) + j + nx + 3);
            cellID++;

            std::fprintf(outCells, "%4lu    %4d  %4d  %4d\n", cellID, i * (nx + 1) + j + 1,
                         i * (nx + 1) + j + nx + 3, i * (nx + 1) + j + nx + 2);
            cellID++;
        }
    }

    for (int j = 0; j < nx; j++) {
        int i = ny - 1;
        std::fprintf(outNodes, "%4lu    %.17g  %.17g    %d\n", i * (nx + 1) + j + nx + 2,
                     LEFT_BOUNDARY_X + j * hx, LEFT_BOUNDARY_Y + ny * hy, 1);
        if (j == nx - 1) {
            std::fprintf(outNodes, "%4lu    %.17g  %.17g    %d\n", i * (nx + 1) + j + nx + 3,
                         LEFT_BOUNDARY_X + nx * hx, LEFT_BOUNDARY_Y + ny * hy, 1);
        }
    }

    assert(nEdges == edgeID - 1);
    assert(nCells == cellID - 1);

    std::fclose(outNodes);
    std::fclose(outEdges);
    std::fclose(outCells);

    return 0;
}

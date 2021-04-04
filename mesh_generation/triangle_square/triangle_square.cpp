#include <cassert>
#include <cmath>
#include <iostream>

double LEFT_BOUNDARY_X = 0;
double RIGHT_BOUNDARY_X = 50;
double LEFT_BOUNDARY_Y = 0;
double RIGHT_BOUNDARY_Y = 50;

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cout << "Number of cells by dimensions must be specified" << std::endl;
        return 1;
    }
    int nx = std::atoi(argv[1]);
    int ny = std::atoi(argv[2]);
    int mode = std::atoi(argv[3]);
    assert(nx > 0);
    assert(ny > 0);
    assert((mode == 1) || (mode == 0));

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
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < nx; j++) {
            int boundNode = 0;
            if (i == 0 || j == 0 || j == ny - 1 || i == nx - 1) {
                boundNode = 1;
            }
            std::fprintf(outNodes, "%4d    %.17g  %.17g    %d\n", i * (nx + 1) + j + 1,
                         LEFT_BOUNDARY_X + j * hx, LEFT_BOUNDARY_Y + i * hy, boundNode);

            if (j == ny - 1) {
                std::fprintf(outNodes, "%4d    %.17g  %.17g    %d\n", i * (nx + 1) + j + 2,
                             LEFT_BOUNDARY_X + nx * hx, LEFT_BOUNDARY_Y + i * hy, boundNode);
            }

            int boundEdge = 0;
            if (j == 0) {
                boundEdge = 1;
            }
            std::fprintf(outEdges, "%4lu   %d  %d  %d\n", edgeID, i * (nx + 1) + j + 1,
                         i * (nx + 1) + j + nx + 2, boundEdge);
            edgeID++;

            boundEdge = 0;
            if (i == 0) {
                boundEdge = 1;
            }
            std::fprintf(outEdges, "%4lu   %d  %d  %d\n", edgeID, i * (nx + 1) + j + 1,
                            i * (nx + 1) + j + 2, boundEdge);
            edgeID++;

            boundEdge = 0;
            if (mode == 1) {
                if (i % 2 == 0) {
                    if (j % 2 == 0) {
                        std::fprintf(outEdges, "%4lu   %d  %d  %d\n", edgeID, i * (nx + 1) + j + 1,
                                     i * (nx + 1) + j + nx + 3, boundEdge);
                    } else {
                        std::fprintf(outEdges, "%4lu   %d  %d  %d\n", edgeID, i * (nx + 1) + j + 2,
                                     i * (nx + 1) + j + nx + 2, boundEdge);
                    }
                } else {
                    if (j % 2 == 0) {
                        std::fprintf(outEdges, "%4lu   %d  %d  %d\n", edgeID, i * (nx + 1) + j + 2,
                                     i * (nx + 1) + j + nx + 2, boundEdge);
                    } else {
                        std::fprintf(outEdges, "%4lu   %d  %d  %d\n", edgeID, i * (nx + 1) + j + 1,
                                     i * (nx + 1) + j + nx + 3, boundEdge);
                    }
                }
                edgeID++;
            } else {
                std::fprintf(outEdges, "%4lu   %d  %d  %d\n", edgeID, i * (nx + 1) + j + 1,
                         i * (nx + 1) + j + nx + 3, boundEdge);
                edgeID++;
            }

            if (i == ny - 1) {
                boundEdge = 1;
                std::fprintf(outEdges, "%4lu   %d  %d  %d\n", edgeID, i * (nx + 1) + j + nx + 2,
                             i * (nx + 1) + j + nx + 3, boundEdge);
                edgeID++;
            }

            if (j == nx - 1) {
                boundEdge = 1;
                std::fprintf(outEdges, "%4lu   %d  %d  %d\n", edgeID, i * (nx + 1) + j + 2,
                             i * (nx + 1) + j + nx + 3, boundEdge);
                edgeID++;
            }

            if (mode == 1) {
                if (i % 2 == 0) {
                    if (j % 2 == 0) {
                        std::fprintf(outCells, "%4lu    %4d  %4d  %4d\n", cellID, i * (nx + 1) + j + 1,
                                     i * (nx + 1) + j + 2, i * (nx + 1) + j + nx + 3);
                        cellID++;

                        std::fprintf(outCells, "%4lu    %4d  %4d  %4d\n", cellID, i * (nx + 1) + j + 1,
                                     i * (nx + 1) + j + nx + 3, i * (nx + 1) + j + nx + 2);
                        cellID++;
                    } else {
                        std::fprintf(outCells, "%4lu    %4d  %4d  %4d\n", cellID, i * (nx + 1) + j + 1,
                                     i * (nx + 1) + j + 2, i * (nx + 1) + j + nx + 2);
                        cellID++;

                        std::fprintf(outCells, "%4lu    %4d  %4d  %4d\n", cellID, i * (nx + 1) + j + 2,
                                     i * (nx + 1) + j + nx + 3, i * (nx + 1) + j + nx + 2);
                        cellID++;
                    }
                } else {
                    if (j % 2 == 0) {
                        std::fprintf(outCells, "%4lu    %4d  %4d  %4d\n", cellID, i * (nx + 1) + j + 1,
                                     i * (nx + 1) + j + 2, i * (nx + 1) + j + nx + 2);
                        cellID++;

                        std::fprintf(outCells, "%4lu    %4d  %4d  %4d\n", cellID, i * (nx + 1) + j + 2,
                                     i * (nx + 1) + j + nx + 3, i * (nx + 1) + j + nx + 2);
                        cellID++;
                    } else {
                        std::fprintf(outCells, "%4lu    %4d  %4d  %4d\n", cellID, i * (nx + 1) + j + 1,
                                     i * (nx + 1) + j + 2, i * (nx + 1) + j + nx + 3);
                        cellID++;

                        std::fprintf(outCells, "%4lu    %4d  %4d  %4d\n", cellID, i * (nx + 1) + j + 1,
                                     i * (nx + 1) + j + nx + 3, i * (nx + 1) + j + nx + 2);
                        cellID++;
                    }
                }
            } else {
                std::fprintf(outCells, "%4lu    %4d  %4d  %4d\n", cellID, i * (nx + 1) + j + 1,
                             i * (nx + 1) + j + 2, i * (nx + 1) + j + nx + 3);
                cellID++;

                std::fprintf(outCells, "%4lu    %4d  %4d  %4d\n", cellID, i * (nx + 1) + j + 1,
                             i * (nx + 1) + j + nx + 3, i * (nx + 1) + j + nx + 2);
                cellID++;
            }
        }
    }

    for (int j = 0; j < nx; j++) {
        int i = ny - 1;
        std::fprintf(outNodes, "%4d    %.17g  %.17g    %d\n", i * (nx + 1) + j + nx + 2,
                     LEFT_BOUNDARY_X + j * hx, LEFT_BOUNDARY_Y + ny * hy, 1);
        if (j == nx - 1) {
            std::fprintf(outNodes, "%4d    %.17g  %.17g    %d\n", i * (nx + 1) + j + nx + 3,
                         LEFT_BOUNDARY_X + nx * hx, LEFT_BOUNDARY_Y + ny * hy, 1);
        }
    }

    assert(nEdges == edgeID - 1);
    assert(nCells == cellID - 1);

    std::fclose(outNodes);
    std::fclose(outEdges);
    std::fclose(outCells);

    std::cout << "number of nodes = " << nNodes << std::endl;
    std::cout << "number of edges = " << nEdges << std::endl;
    std::cout << "number of cells = " << nCells << std::endl;

    return 0;
}

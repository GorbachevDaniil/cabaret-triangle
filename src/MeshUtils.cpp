#include "MeshUtils.hpp"

#include "Mesh.hpp"
#include "Cell.hpp"
#include "Edge.hpp"
#include "Node.hpp"
#include "Data.hpp"
#include "Vector.hpp"

#include <vector>

void calculateNodeNormal(Mesh &mesh, Cell &cell, long edgeID) {
    Edge edge = mesh.edges[edgeID];

    std::vector<long> nodeIDs = cell.getEdgeOrderedNodeIDs(edge.nodeIDs);

    Vector tangential;
    tangential.set(mesh.nodes[nodeIDs[1]].data.coords);
    tangential.minus(mesh.nodes[nodeIDs[0]].data.coords);

    Vector normal;
    normal.set(tangential.y, -tangential.x);
    normal.mult(1. / edge.length);

    mesh.nodes[edge.centerNodeID].normal.set(normal);
}

void MeshUtils::calculateNodeNormals(Mesh &mesh) {
    std::vector<int> isNormalCalculated(mesh.edges.size(), 0);
    for (unsigned long i = 0; i < mesh.cells.size(); i++) {
        Cell *cell = &mesh.cells[i];

        for (long edgeID : cell->edgeIDs) {
            if (isNormalCalculated[edgeID] == 1) {
                cell->edgeToNormalDirection[edgeID] = -1;
            } else {
                calculateNodeNormal(mesh, *cell, edgeID);
                isNormalCalculated[edgeID] = 1;
                cell->edgeToNormalDirection[edgeID] = 1;
            }
        }
    }
}

void MeshUtils::calculateVectorsFromCenterToEdges(Mesh &mesh) {
    for (unsigned long i = 0; i < mesh.cells.size(); i++) {
        Cell *cell = &mesh.cells[i];
        Node cellCenter = mesh.nodes[cell->centerNodeID];

        for (long edgeID : cell->edgeIDs) {
            Edge edge = mesh.edges[edgeID];
            Node edgeCenter = mesh.nodes[edge.centerNodeID];

            Vector vectorFromCenter;
            vectorFromCenter.set(edgeCenter.data.coords);
            vectorFromCenter.minus(cellCenter.data.coords);
            vectorFromCenter.mult(1. / vectorFromCenter.length());

            cell->edgeToVectorFromCenter[edgeID].set(vectorFromCenter);
        }
    }
}

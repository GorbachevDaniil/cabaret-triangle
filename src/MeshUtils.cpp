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

    Vector tangential(mesh.nodes[nodeIDs[1]].data.coords -
                      mesh.nodes[nodeIDs[0]].data.coords);

    Vector normal(tangential.y, -tangential.x);
    normal /= edge.length;

    mesh.nodes[edge.centerNodeID].normal = normal;
}

void MeshUtils::calculateNodeNormals(Mesh &mesh) {
    std::vector<int> isNormalCalculated(mesh.edges.size(), 0);
    for (unsigned long i = 0; i < mesh.cells.size(); i++) {
        Cell *cell = &mesh.cells[i];

        for (long edgeID : cell->edgeIDs) {
            if (isNormalCalculated[edgeID] == 1) {
                cell->edgeToNormalDir[edgeID] = -1;
            } else {
                calculateNodeNormal(mesh, *cell, edgeID);
                isNormalCalculated[edgeID] = 1;
                cell->edgeToNormalDir[edgeID] = 1;
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

            Vector vectorFromCenter(edgeCenter.data.coords - cellCenter.data.coords);
            vectorFromCenter /= vectorFromCenter.length();

            cell->edgeToTransportDir[edgeID] = vectorFromCenter;
        }
    }
}

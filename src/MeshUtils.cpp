#include "MeshUtils.hpp"

#include "Mesh.hpp"
#include "Cell.hpp"
#include "Edge.hpp"
#include "Node.hpp"
#include "Data.hpp"
#include "Vector.hpp"

#include <vector>

void calculateNormal(Mesh &mesh, Cell &cell, long edgeID) {
    Edge edge = mesh.edges[edgeID];

    std::vector<long> nodeIDs = cell.getEdgeOrderedNodeIDs(edge.nodeIDs);

    Vector tangential;
    tangential.set(mesh.nodes[nodeIDs[1]].data.coords);
    tangential.minus(mesh.nodes[nodeIDs[0]].data.coords);

    Vector normal;
    normal.set(tangential.y, -tangential.x);
    normal.mult(1 / tangential.length());

    mesh.nodes[edge.centerNodeID].normal.set(normal);
}

void MeshUtils::calculateNormals(Mesh &mesh) {
    std::vector<int> isNormalCalculated(mesh.edges.size(), 0);
    for (Cell cell : mesh.cells) {
        for (long edgeID : cell.edgeIDs) {
            if (isNormalCalculated[edgeID] == 1) {
                cell.edgeIDNormalDirection[edgeID] = -1;
            } else {
                calculateNormal(mesh, cell, edgeID);
                isNormalCalculated[edgeID] = 1;
                cell.edgeIDNormalDirection[edgeID] = 1;
            }
        }
    }
}
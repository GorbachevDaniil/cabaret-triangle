#include "MeshUtils.hpp"

#include "Mesh.hpp"
#include "Cell.hpp"
#include "Edge.hpp"
#include "Node.hpp"
#include "Data.hpp"
#include "Vector.hpp"

#include <vector>

void calculateEdgeNormal(Mesh &mesh, Cell &cell, long edgeID) {
    Edge *edge = &mesh.edges[edgeID];

    std::vector<long> nodeIDs = cell.getEdgeOrderedNodeIDs(edge->edgeEndsNodeIDs);

    Vector tangential(mesh.nodes[nodeIDs[1]].data.coords -
                      mesh.nodes[nodeIDs[0]].data.coords);

    Vector normal(tangential.y, -tangential.x);
    normal = normal / edge->length;

    edge->normal = normal;
}

void MeshUtils::calculateEdgesNormals(Mesh &mesh) {
    std::vector<int> isNormalCalculated(mesh.edges.size(), 0);
    for (unsigned long i = 0; i < mesh.cells.size(); i++) {
        Cell *cell = &mesh.cells[i];

        for (long edgeID : cell->edgeIDs) {
            if (isNormalCalculated[edgeID] == 1) {
                cell->edgeToNormalDir[edgeID] = -1;
            } else {
                calculateEdgeNormal(mesh, *cell, edgeID);
                isNormalCalculated[edgeID] = 1;
                cell->edgeToNormalDir[edgeID] = 1;
            }
        }
    }
}

void MeshUtils::calculateTransferVectors(Mesh &mesh) {
    for (unsigned long i = 0; i < mesh.cells.size(); i++) {
        Cell *cell = &mesh.cells[i];
        Node *cellNode = &mesh.nodes[cell->centerNodeID];

        for (unsigned long edgeID : cell->edgeIDs) {
            for (unsigned long nodeID : mesh.edges[edgeID].getUsedNodes(mesh)) {
                Node *node = &mesh.nodes[nodeID];

                Vector transferVector(node->data.coords - cellNode->data.coords);
                transferVector = transferVector / transferVector.length();

                cell->nodeToTransferVector[nodeID] = transferVector;
            }
        }
    }
}

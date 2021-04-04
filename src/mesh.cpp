#include "mesh.hpp"

#include <string>

#include "parser.hpp"

void Mesh::calculate_edge_normal(Cell* cell, long edgeID) {
    Edge *edge = &this->edges[edgeID];

    std::vector<long> nodeIDs = cell->getEdgeOrderedNodeIDs(edge->endNodeIDs);

    Vector tangential(this->nodes[nodeIDs[1]].data.coords -
        this->nodes[nodeIDs[0]].data.coords);

    Vector normal(tangential.y, -tangential.x);
    normal = normal / edge->length;

    edge->normal = normal;
}

void Mesh::calculate_edges_normals() {
    std::vector<int> isNormalCalculated(this->edges.size(), 0);
    for (unsigned long i = 0; i < this->cells.size(); i++) {
        Cell *cell = &this->cells[i];

        for (long edgeID : cell->edgeIDs) {
            if (isNormalCalculated[edgeID] == 1) {
                cell->edgeToNormalDir[edgeID] = -1;
            } else {
                this->calculate_edge_normal(cell, edgeID);
                isNormalCalculated[edgeID] = 1;
                cell->edgeToNormalDir[edgeID] = 1;
            }
        }
    }
}

void Mesh::calculate_transfer_vectors() {
    for (unsigned long i = 0; i < this->cells.size(); i++) {
        Cell *cell = &this->cells[i];
        Node *cellNode = &this->nodes[cell->centerNodeID];

        for (unsigned long edgeID : cell->edgeIDs) {
            for (unsigned long nodeID : this->edges[edgeID].usedNodeIDs) {
                Node *node = &this->nodes[nodeID];

                Vector transferVector(node->data.coords - cellNode->data.coords);
                transferVector = transferVector / transferVector.length();

                cell->nodeToTransferVector[nodeID] = transferVector;
            }
        }
    }
}

void Mesh::init_mesh(const std::string& path_to_noad_file,
                     const std::string& path_to_edge_file,
                     const std::string& path_to_ele_file) {
    Parser parser;
    parser.LoadNodes(this, path_to_noad_file);
    parser.LoadEdges(this, path_to_edge_file);
    parser.LoadCells(this, path_to_ele_file);

    long nodesSize = this->nodes.size();

    this->s0 = std::vector<std::vector<double>>(nodesSize);
    this->s1 = std::vector<std::vector<double>>(nodesSize);
    this->s2 = std::vector<std::vector<double>>(nodesSize);

    this->v0 = std::vector<std::vector<Vector>>(nodesSize);
    this->v1 = std::vector<std::vector<Vector>>(nodesSize);
    this->v2 = std::vector<std::vector<Vector>>(nodesSize);
}
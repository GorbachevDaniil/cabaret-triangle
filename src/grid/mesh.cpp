#include "grid/mesh.hpp"

#include <iostream>
#include <string>

#include "parser.hpp"

void Mesh::calculate_edge_normal(Cell* cell, long edge_id) {
    Edge *edge = &this->edges[edge_id];

    std::vector<long> nodeIDs = cell->get_edge_ordered_node_ids(edge->end_node_ids);

    double tangential_x = this->nodes[nodeIDs[1]].coords.x - this->nodes[nodeIDs[0]].coords.x;
    double tangential_y = this->nodes[nodeIDs[1]].coords.y - this->nodes[nodeIDs[0]].coords.y;

    cell->edge_to_div_coef[edge_id].x = tangential_y;
    cell->edge_to_div_coef[edge_id].y = -tangential_x;
    cell->edge_to_normal[edge_id].x = tangential_y / edge->length;
    cell->edge_to_normal[edge_id].y = -tangential_x / edge->length;
}

void Mesh::calculate_edges_normals() {
    for (unsigned long i = 0; i < this->cells.size(); i++) {
        Cell *cell = &this->cells[i];

        long edge_id_1 = cell->edge_ids[0];
        long edge_id_2 = cell->edge_ids[1];
        long edge_id_3 = cell->edge_ids[2];

        calculate_edge_normal(cell, edge_id_1);
        calculate_edge_normal(cell, edge_id_2);

        cell->edge_to_div_coef[edge_id_3].x = -(cell->edge_to_div_coef[edge_id_1].x +
                                                cell->edge_to_div_coef[edge_id_2].x);
        cell->edge_to_div_coef[edge_id_3].y = -(cell->edge_to_div_coef[edge_id_1].y +
                                                cell->edge_to_div_coef[edge_id_2].y);

        double edge_len_3 = edges[edge_id_3].length;
        cell->edge_to_normal[edge_id_3].x = cell->edge_to_div_coef[edge_id_3].x / edge_len_3;
        cell->edge_to_normal[edge_id_3].y = cell->edge_to_div_coef[edge_id_3].y / edge_len_3;
    }

    bool exit_flag = false;
    for (unsigned long i = 0; i < this->cells.size(); i++) {
        Cell *cell = &this->cells[i];

        long edge_id_1 = cell->edge_ids[0];
        long edge_id_2 = cell->edge_ids[1];
        long edge_id_3 = cell->edge_ids[2];

        double summ_x = cell->edge_to_div_coef[edge_id_1].x +
                        cell->edge_to_div_coef[edge_id_2].x +
                        cell->edge_to_div_coef[edge_id_3].x;
        double summ_y = cell->edge_to_div_coef[edge_id_1].y +
                        cell->edge_to_div_coef[edge_id_2].y +
                        cell->edge_to_div_coef[edge_id_3].y;
        if (summ_x > 0 || summ_y > 0) {
            std::cout << i << " " << summ_x << " " << summ_y << std::endl;
            exit_flag = true;
        }
    }
    if (exit_flag) {
        exit(0);
    }
}

void Mesh::calculate_transfer_vectors() {
    for (unsigned long i = 0; i < this->cells.size(); i++) {
        Cell *cell = &this->cells[i];
        Node *cell_node = &this->nodes[cell->center_node_id];

        for (long edge_id : cell->edge_ids) {
            for (long node_id : this->edges[edge_id].used_node_ids) {
                Node *node = &this->nodes[node_id];

                double transfer_x = cell_node->coords.x - node->coords.x;
                double transfer_y = cell_node->coords.y - node->coords.y;
                double transfer_len = Vector::length(transfer_x, transfer_y);

                cell->node_to_transfer_vector[node_id].x = transfer_x / transfer_len;
                cell->node_to_transfer_vector[node_id].y = transfer_y / transfer_len;
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
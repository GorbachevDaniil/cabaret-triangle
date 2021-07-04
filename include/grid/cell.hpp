#ifndef Cell_hpp
#define Cell_hpp

#include <unordered_map>
#include <vector>

#include "vector.hpp"

class Mesh;

class Cell {
public:
    Cell(Mesh &mesh, long id, long node_id_1, long node_id_2, long node_id_3);

    long id;
    long center_node_id;
    double volume;

    std::vector<long> node_ids;
    std::vector<long> edge_ids;
    std::unordered_map<long, Vector> edge_to_normal;
    std::unordered_map<long, Vector> edge_to_div_coef;
    std::unordered_map<long, Vector> node_to_transfer_vector;
    std::unordered_map<long, long> node_id_to_opposite_node_id;

    long get_next_node_id(unsigned long nodeIDPos);
    long get_prev_node_id(unsigned long nodeIDPos);
    std::vector<long> get_edge_ordered_node_ids(std::vector<long> edgeUnorderedNodeIDs);

private:
    void assign_opposite_node_ids(Mesh &mesh);
};

#endif

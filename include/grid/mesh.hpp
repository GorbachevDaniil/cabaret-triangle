#ifndef Mesh_hpp
#define Mesh_hpp

#include <map>
#include <vector>

#include "node.hpp"
#include "edge.hpp"
#include "cell.hpp"

class Mesh {
public:
    Mesh(int edge_inner_nodes, bool apex_nodes_used) :
        edge_inner_nodes(edge_inner_nodes),
        apex_nodes_used(apex_nodes_used)
    {}

    Mesh(Mesh&&) = delete;
    Mesh(const Mesh&) = delete;
    Mesh& operator=(Mesh&&) = delete;
    Mesh& operator=(const Mesh&) = delete;

    int edge_inner_nodes;
    bool apex_nodes_used;

    // list of scalar variables on n layer of time
    std::vector<std::vector<double>> s0;
    // list of scalar variables on n+1/2 layer of time
    std::vector<std::vector<double>> s1;
    // list of scalar variables on n+1 layer of time
    std::vector<std::vector<double>> s2;

    // list of vector variables on n layer of time
    std::vector<std::vector<Vector>> v0;
    // list of vector variables on n+1/2 layer of time
    std::vector<std::vector<Vector>> v1;
    // list of vector variables on n+1 layer of time
    std::vector<std::vector<Vector>> v2;

    std::vector<Node> nodes;
    std::vector<Edge> edges;
    std::vector<Cell> cells;

    std::map<std::pair<int, int>, int> map_nodes_with_edge;

    inline long get_new_node_id() { return nodes.size(); };

    void init_mesh(const std::string& path_to_noad_file,
                   const std::string& path_to_edge_file,
                   const std::string& path_to_ele_file);

    void calculate_edges_normals();
    void calculate_transfer_vectors();
private:
    void calculate_edge_normal(Cell* cell, long edge_id);
};

#endif
